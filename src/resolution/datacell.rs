use std::sync::{OnceLock, atomic::AtomicPtr};

use ahash::HashMap;
use algebra::{abelian::Abelian, matrix::Matrix, ring::CRing};
use deepsize::DeepSizeOf;
use itertools::Itertools;

use crate::{
    grading::grading::{GradedIndexing, Grading},
    resolution::cokercell::CokerCell,
    traits::Coalgebra,
    types::{AGeneratorIndex, CoalgebraIndexType, CofreeBasis, ComoduleIndexType},
};

#[allow(type_alias_bounds)]
pub type AGen<G, C: Coalgebra<G>> = (
    <C::RingMorph as Abelian<C::BaseRing>>::Generator,
    AGeneratorIndex,
    ComoduleIndexType,
);
#[allow(type_alias_bounds)]
pub type RGen<G, C: Coalgebra<G>> = (
    CofreeBasis<G>,
    <C::RingMorph as Abelian<C::BaseRing>>::Generator,
);
pub type LUT<G> = HashMap<CofreeBasis<G>, u32>;
#[allow(type_alias_bounds)]
pub type Module<G, C: Coalgebra<G>> = Vec<<C::RingMorph as Abelian<C::BaseRing>>::Generator>;
#[allow(type_alias_bounds)]
pub type ToViLUT<G, C: Coalgebra<G>> = HashMap<CofreeBasis<G>, Vec<(C::BaseRing, AGeneratorIndex)>>;

/// Everything below is wrt. some degree g
///
/// to_cofree: morphism from Q_i to A \otimes V_i
///
/// a_gens: V_i := (module_generator, generated_index, coker_id)
///         generated_index should be a unique integer in the whole of s
///         note that V_i can be seen as a submodule of the cokernel using coker_id
///
/// r_gens: A \otimes V_i  
///
/// lut: A \otimes V_i -> to its index in r_gens  
///
/// lut2: For each basis element in  A \otimes V_i-1 wrt (alg_gr, alg_id), v_i index (so not r_gens basis!)
///         gives a list of elements in V_i \subset Q_i to which it maps
///       One could also see this as the R linear map from A \otimes V_i-1 \to V_i
///         given by the adjunction
///
#[derive(Debug, DeepSizeOf)]
pub struct DataCell<G: Grading, C: Coalgebra<G>> {
    pub transposed_to_cofree: C::RingMorph,
    pub a_gens: Vec<AGen<G, C>>,
    pub r_gens: Vec<RGen<G, C>>, // TODO : Seperate RGens with module structure
    pub lut: LUT<G>,
    pub lut2: ToViLUT<G, C>,
}

impl<G: Grading, C: Coalgebra<G>> DataCell<G, C> {
    fn luts<A: Send + Sync>(
        gs: &G::ContiguousMemory<(A, OnceLock<CokerCell<G, C>>, OnceLock<DataCell<G, C>>)>,
        degree: G,
        coalgebra: &C,
    ) -> (Vec<(AGen<G, C>, G)>, Vec<RGen<G, C>>, LUT<G>, ToViLUT<G, C>) {
        // a_gens (without the a_gen found this cycle)
        let a_gens: Vec<(AGen<G, C>, G)> = degree
            .iterator_from_zero(false)
            .iter()
            .flat_map(|g| {
                (gs.get(g.to_index()).2.get_or_init(|| unreachable!()).a_gens)
                    .iter()
                    .map(|x| (x.clone(), *g))
            })
            .collect();

        // r_gens (without 1 \otimes v_i found this compute cycle!)
        let mut r_gens = vec![];
        let mut lut = HashMap::default();
        // sort by module generator (important for k[t] modules)
        for ((a_module_gen, a_generator, _), g) in a_gens.iter().sorted_by_key(|a| a.0.0) {
            let a_g = degree - *g;
            for i in 0..coalgebra.size_in_degree(a_g) {
                let coalg_basis = ((a_g, i as CoalgebraIndexType), *a_generator);
                lut.insert(coalg_basis, r_gens.len() as ComoduleIndexType);
                r_gens.push((coalg_basis, a_module_gen.clone()));
            }
        }

        let mut reduced_to_coker = HashMap::default();

        for g in degree.iterator_from_zero(false) {
            let prev_s_g = gs.get(g.to_index()).2.get_or_init(|| unreachable!());
            reduced_to_coker.extend(prev_s_g.lut2.clone().into_iter());
        }

        debug_assert!(r_gens.is_sorted_by_key(|x| x.1));

        (a_gens, r_gens, lut, reduced_to_coker)
    }

    pub fn map_to_cofree() {}

    pub fn boundary_resolve<A: Send + Sync>(
        gs: &G::ContiguousMemory<(A, OnceLock<CokerCell<G, C>>, OnceLock<DataCell<G, C>>)>,
        degree: G,
        coalgebra: &C,
    ) -> Self {
        
        let (_, r_gens, lut, _) =
            DataCell::luts(gs, degree, coalgebra);
        let transposed_to_cofree = C::RingMorph::zero(r_gens.len(), 0);

        

        DataCell {
            transposed_to_cofree,
            a_gens: vec![],
            r_gens,
            lut,
            lut2: HashMap::default(),
        }
    }

    pub fn coker_to_cofree<A: Send + Sync>(
        prev_s: &DataCell<G, C>,
        coker: &CokerCell<G, C>,
        coalgebra: &C,
        reduced_to_coker: &HashMap<((G, u16), u16), Vec<(<C as Coalgebra<G>>::BaseRing, u16)>>,
        lut: &HashMap<((G, u16), u16), u32>,
        r_gens: &Vec<RGen<G,C>>
    ) -> C::RingMorph {
        let mut small_to_cofree =
            <C::RingMorph as Matrix<C::BaseRing>>::zero(coker.cokernel.len(), r_gens.len());
        // let coact_ref = &mut small_to_cofree as *mut C::RingMorph;
        // let map_ref = AtomicPtr::new(coact_ref);

        // TODO : PAR Iter < might not want this ! as the function is "to light"
        // (0..cokernel.len()).into_iter().collect::<Vec<_>>().into_par_iter().with_min_len(32).for_each(|coker_id| {
        (0..coker.cokernel.len()).into_iter().for_each(|coker_id| {
            for codom_id in 0..coker.repr_vecs.codomain() {
                let value = coker.repr_vecs.get(coker_id, codom_id);
                if value.is_zero() {
                    continue;
                }

                // Now we have an element in A \otimes V_i-1 which maps with value to coker_id
                // can be seen as a_k \otimes v_k (in A \otimes V_i-1)
                // where a_k = alg_gr, alg_id and v_k = gen_id
                let (((alg_gr, alg_id), gen_id), _) = prev_s.r_gens[codom_id];

                // Get all coactions
                // these are \sum b_j \otimes c_j \otimes v_k
                for (alg_l, alg_r, coact_val) in coalgebra.coaction((alg_gr, alg_id)) {
                    // I want to know if there is a c_j \otimes v_k mapping to some_i q_i which generates an A
                    match reduced_to_coker.get(&(*alg_r, gen_id)) {
                        Some(targets) => {
                            debug_assert!(targets.len() > 0);
                            for (v_i_value, a_generator) in targets {
                                // This is the tagret index in r_gens
                                let target_id = lut[&(*alg_l, *a_generator)];

                                // As coker_id is seperate across parallel instances
                                // This unsafe code is fine, AS LONG as the matrix is FlatMatrix :)
                                // TODO : This is not reallly generic, and depends on the underlying implementation
                                // This probably breaks for a F2 matrix implementation.
                                let final_val = value * *coact_val * *v_i_value;
                                small_to_cofree.add_at(coker_id, target_id as usize, final_val);
                                // unsafe {
                                //     (**(map_ref.as_ptr())).add_at(
                                //         coker_id,
                                //         target_id as usize,

                                //     );
                                // }
                            }
                        }
                        None => {}
                    }
                }
            }
        });

        small_to_cofree
    }

    pub fn resolve<A: Send + Sync>(
        gs: &G::ContiguousMemory<(A, OnceLock<CokerCell<G, C>>, OnceLock<DataCell<G, C>>)>,
        prev_s: &DataCell<G, C>,
        coker: &CokerCell<G, C>,
        degree: G,
        coalgebra: &C,
    ) -> Self {
        // Compute Cokernel of prev_s
        // Get all previous A generators and generate corresponding r_gens and lut for this degree
        // The only elements missing from r_gens and lut are the new a generators found later

        let (gs_a_gens, mut r_gens, mut lut, reduced_to_coker) =
            DataCell::luts(gs, degree, coalgebra);

        let small_to_cofree = DataCell::coker_to_cofree::<A>(prev_s, coker, coalgebra, &reduced_to_coker, &lut, &r_gens);

        // Now we need to find the kernel of the map Q_i -> A \otimes V_i
        // And see if we need to add any extra A generators
        let mut a_gens: Vec<AGen<G, C>> = vec![];
        let mut lut2 = HashMap::default();
        for (id, q_index) in small_to_cofree
            .kernel_destroyers(&coker.cokernel, &r_gens.iter().map(|x| x.1).collect())
            .iter()
            .enumerate()
        {
            let a_gen_id = id + gs_a_gens.len();
            let module_gen = coker.cokernel[*q_index];
            a_gens.push((module_gen, a_gen_id as u16, *q_index as u32));

            // Now we should figure out which elements map to q_index
            // Note that repr_vectors is not this!
            // But we can just look in the row of to_cokernel
            for prev_s_v_i in 0..coker.to_cokernel.domain() {
                let val = coker.to_cokernel.get(prev_s_v_i, *q_index);
                let ((r_alg, prev_a_gen), _) = prev_s.r_gens[prev_s_v_i];

                if !val.is_zero() {
                    let v = lut2.entry((r_alg, prev_a_gen)).or_insert(vec![]);
                    v.push((val, a_gen_id as u16));
                }
            }
        }

        let (r_gens_map, a_gen_idxes) = insert_many_by_key(
            &mut r_gens,
            a_gens
                .iter()
                .map(|x| (((G::zero(), 0), x.1), x.0))
                .collect(),
            |x| x.0,
        );

        // to_cofree and lut are not correct right now !
        // we need to include the information from a_gens in the correct place right now
        lut.iter_mut()
            .for_each(|(_, x)| *x = r_gens_map[*x as usize] as u32);
        a_gens.iter().enumerate().for_each(|(a_gen_id, a_gen)| {
            lut.insert(((G::zero(), 0), a_gen.1), a_gen_idxes[a_gen_id] as u32);
        });

        let mut transposed_to_cofree =
            <C::RingMorph as Matrix<C::BaseRing>>::zero(r_gens.len(), coker.cokernel.len());
        for (a_gen_id, target_id) in a_gen_idxes.iter().enumerate() {
            let q_index = a_gens[a_gen_id].2;
            transposed_to_cofree.set(*target_id, q_index as usize, C::BaseRing::one());
        }

        let coact_ref = &mut transposed_to_cofree as *mut C::RingMorph;
        let map_ref = AtomicPtr::new(coact_ref);

        r_gens_map
            .iter()
            .enumerate()
            .for_each(|(old_r_gen, new_r_gen)| {
                for c in 0..coker.cokernel.len() {
                    let row = small_to_cofree.get(c, old_r_gen);
                    unsafe { (**(map_ref.as_ptr())).set(*new_r_gen, c, row) }
                }
            });

        Self {
            transposed_to_cofree,
            a_gens,
            r_gens,
            lut,
            lut2,
        }
    }
    
    pub fn partial_resolve<A: Send + Sync>(
        gs: &G::ContiguousMemory<(A, OnceLock<CokerCell<G, C>>, OnceLock<DataCell<G, C>>)>,
        prev_s: &DataCell<G, C>,
        coker: &CokerCell<G, C>,
        degree: G,
        coalgebra: &C,
    ) -> Self {
        let (_, r_gens, lut, reduced_to_coker) =
            DataCell::luts(gs, degree, coalgebra);

        let small_to_cofree = DataCell::coker_to_cofree::<A>(prev_s, coker, coalgebra, &reduced_to_coker, &lut, &r_gens);

        Self {
            transposed_to_cofree: small_to_cofree.transpose(),
            a_gens: vec![],
            r_gens,
            lut,
            lut2: HashMap::default(),
        }
    }
}

// Helper function

pub fn insert_many_by_key<T, K: Ord + Clone, F>(
    sorted: &mut Vec<T>,
    new_elements: Vec<T>,
    key_fn: F,
) -> (Vec<usize>, Vec<usize>)
where
    F: Fn(&T) -> K,
    T: Clone,
{
    // Attach original indices to new elements
    let mut new_with_idx: Vec<(usize, T)> = new_elements.into_iter().enumerate().collect();

    // Sort new elements by key
    new_with_idx.sort_by_key(|(_, v)| key_fn(v));

    let mut result = Vec::with_capacity(sorted.len() + new_with_idx.len());
    let mut mapping_new = vec![0; new_with_idx.len()];
    let mut mapping_old = vec![0; sorted.len()];

    let mut i = 0; // index in sorted
    let mut j = 0; // index in new_with_idx

    while i < sorted.len() && j < new_with_idx.len() {
        if key_fn(&sorted[i]) <= key_fn(&new_with_idx[j].1) {
            mapping_old[i] = result.len();
            result.push(sorted[i].clone());
            i += 1;
        } else {
            mapping_new[new_with_idx[j].0] = result.len();
            result.push(new_with_idx[j].1.clone());
            j += 1;
        }
    }

    while i < sorted.len() {
        mapping_old[i] = result.len();
        result.push(sorted[i].clone());
        i += 1;
    }

    while j < new_with_idx.len() {
        mapping_new[new_with_idx[j].0] = result.len();
        result.push(new_with_idx[j].1.clone());
        j += 1;
    }

    *sorted = result;
    (mapping_old, mapping_new)
}
