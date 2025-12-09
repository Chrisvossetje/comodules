use std::sync::Arc;

use ahash::HashMap;
use algebra::{abelian::Abelian, field::Field, matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::CRing, rings::univariate_polynomial_ring::UniPolRing};
use itertools::Itertools;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{
    basiselement::kBasisElement, comodule::{rcomodule::RComodule, traits::{Comodule, ComoduleMorphism}}, graded_module_morphism::GradedModuleMap, graded_space::BasisIndex, grading::{Grading, OrderedGrading}, tensor::TensorMap
};



#[derive(Debug, Clone)]
pub struct RComoduleMorphism<G: Grading, F: Field> {
    pub domain: Arc<RComodule<G, F>>,
    pub codomain: Arc<RComodule<G, F>>,
    pub map: GradedModuleMap<G, F>, 
}




impl<G: Grading, F: Field> ComoduleMorphism<G, RComodule<G,F>>
    for RComoduleMorphism<G, F>
{
    fn cokernel(&self) -> Self {
        if cfg!(debug_assertions) {
            self.map.verify(&self.domain.space, &self.codomain.space).unwrap();
            self.codomain.verify().unwrap()
        }
        
        let (coker_to, coker_inv, coker) = self.map.cokernel::<kBasisElement>(&self.codomain.space);

        

        let coalg = self.codomain.coalgebra.as_ref();
        let tensor = TensorMap::generate(&coalg.space, &coker);


        let codom_lut: HashMap<BasisIndex<G>, Vec<(usize, UniPolRing<F>)>> = self
            .codomain
            .tensor
            // the choice for the tensor here is not neccessary
            // Could also be self.codomain.space and iterate over len of the module 
            .construct
            .par_iter()
            .map(|((f_gr, f_id), _)| {
                // Transfer a specific codomain grade and id (f_gr, f_id) to a list of elements which it maps to in the cokernel
                let v = (0..coker_to
                    .maps
                    .get(f_gr)
                    .map(|map| map.codomain())
                    .unwrap_or(0))
                    .filter_map(|q_index| {
                        let val = coker_to.maps.get(f_gr).unwrap().get(*f_id, q_index);
                        match val.is_zero() {
                            true => None,
                            false => Some((q_index, val)),
                        }
                    })
                    .collect();
                ((*f_gr, *f_id), v)
            })
            .collect();

        let coker_tensor_module = coker.generate_tensor_as_module(&self.codomain.coalgebra.space, &tensor);
        
        let coaction = coker
            .0
            .par_iter()
            .map(|(g, v)| {
                let g_tensor_dimen = tensor.get_dimension(g);
                let mut g_coaction = FlatMatrix::zero(v.len(), g_tensor_dimen);
                let coact_size = self.codomain.tensor.dimensions[g];
                let inv_map = coker_inv.maps.get(&g).unwrap();

                for coker_id in 0..v.len() { 
                    for codom_id in 0..inv_map.codomain() {

                        let inv_val = inv_map.get(coker_id, codom_id);
                        if inv_val.is_zero() {
                            continue;
                        }

                        for codom_coact_id in 0..coact_size {
                            let coact_val = self.codomain.coaction.maps[g].get(codom_id, codom_coact_id);

                            if !coact_val.is_zero() {
                                let ((alg_gr, alg_id), (mod_gr, mod_id)) =
                                    self.codomain.tensor.deconstruct[&(*g, codom_coact_id)];

                                for (target_id, val) in codom_lut.get(&(mod_gr, mod_id)).unwrap() {
                                    let (final_gr, final_id) =
                                        tensor.construct[&(mod_gr, *target_id)][&(alg_gr, alg_id)];
                                    
                                    let tensor_el = &coker_tensor_module.0.get(&final_gr).unwrap()[final_id];
                                    
                                    let final_val = inv_val * coact_val * *val;
                                    
                                    if let Some(tens_el_power) = tensor_el.2 {
                                        if final_val.1 >= tens_el_power {
                                            continue;
                                        }
                                    }
                                    g_coaction.add_at(coker_id, final_id, final_val);
                                }
                            }
                        }
                    }
                }
                (*g, g_coaction)
            }).collect();

        let coaction = GradedModuleMap {
            maps: coaction,
        };

        let coker_comod = RComodule {
            coalgebra: self.codomain.coalgebra.clone(),
            space: coker,
            coaction,
            tensor,
        };

        Self {
            domain: self.codomain.clone(),
            codomain: Arc::new(coker_comod),
            map: coker_to,
        }
    }

    fn inject_codomain_to_cofree(&self, limit: G) -> Self 
    where
        G: OrderedGrading {
        if cfg!(debug_assertions) {
            self.map.verify(&self.domain.space, &self.codomain.space).unwrap();
            self.codomain.verify().unwrap();
        }

        let mut growing_map: GradedModuleMap<G, F> =
            GradedModuleMap::zero_codomain(&self.codomain.space);
        let mut growing_comodule = RComodule::zero_comodule(self.codomain.coalgebra.clone());
        let mut iteration = 0;

        if cfg!(debug_assertions) {
            growing_map.verify(&self.codomain.space, &growing_comodule.space).unwrap()       
        }
        
        let grades: Vec<G> = growing_map
            .maps
            .iter()
            .map(|(g, _)| *g)
            .sorted_by(|a,b| {
                a.compare(b)
            })  
            .filter(|&g| g.compare(&limit).is_le())
            .collect();

        let fixed_limit = limit.incr().incr();

        for pivot_grade in grades {
            // Get lowest graded pivot element
            let map = growing_map.maps.get(&pivot_grade).unwrap();
            let empty = vec![];
            let domain: Vec<_> = self.codomain.space.0.get(&pivot_grade).unwrap_or(&empty).iter().map(|x| x.2).collect();
            let codomain: Vec<_> = growing_comodule.space.0.get(&pivot_grade).unwrap_or(&empty).iter().map(|x| x.2).collect(); 
            
            let l = map.kernel_generators(&domain, &codomain);
            for pivot in l {
                let (_, tau_shift, generator_power) = self.codomain.space.0.get(&pivot_grade).unwrap().get(pivot).unwrap();

                let alg_to_tens = self
                    .codomain
                    .tensor
                    .construct
                    .get(&(pivot_grade, pivot))
                    .expect("The tensor should exist on the codomain in this grade");

                let coalg_space = &self.codomain.coalgebra.space;

                let cofree_map: HashMap<G,FlatMatrix<_>> = coalg_space.0.iter().filter_map(|(alg_gr, alg_gr_space)| {
                    let t_gr = *alg_gr + pivot_grade;

                    if t_gr > fixed_limit {
                        return None;
                    }

                    let codomain_len = self.codomain.space.dimension_in_grade(&t_gr);
                    let coalg_len = alg_gr_space.len();

                    if !alg_to_tens.contains_key(&(*alg_gr, 0)) {
                        let zero_map = FlatMatrix::zero(codomain_len, coalg_len);
                        return Some((t_gr, zero_map));
                    };

                    let mut map = FlatMatrix::zero(codomain_len, coalg_len);

                    for a_id in 0..coalg_len {
                        let (t_gr, t_id) = alg_to_tens.get(&(*alg_gr,a_id)).expect("This BasisIndex should exist on the tensor object in the to inject comodule");
                        let slice = self
                            .codomain
                            .coaction
                            .maps
                            .get(t_gr)
                            .expect("This grade should exist on the coaction of the injecting comodule")
                            .get_row(*t_id);
                        map.set_row(a_id, slice);
                    }

                    Some((t_gr, map))
                }).collect();

                let mut f = RComodule::cofree_comodule(
                    self.codomain.coalgebra.clone(),
                    iteration,
                    pivot_grade,
                    fixed_limit,
                    (*tau_shift, *generator_power)
                );

                let mut extend_map = GradedModuleMap {
                    maps: cofree_map,
                };

                growing_comodule.direct_sum(&mut f);
                growing_map.vstack(&mut extend_map);

                if cfg!(debug_assertions) {
                    growing_map.verify(&self.codomain.space, &growing_comodule.space).unwrap();
                }

                iteration += 1;                
            } 
        }

        Self {
            domain: self.codomain.clone(),
            codomain: Arc::new(growing_comodule),
            map: growing_map,
        }
    }
    
    /// Zero should be the domain, the comodule is the codomain
    fn zero_morphism(comodule: Arc<RComodule<G,F>>) -> Self {
        let zero = Arc::new(RComodule::zero_comodule(comodule.coalgebra.clone()));

        let zero_map = GradedModuleMap::zero(&zero.space, &comodule.space);

        Self {
            domain: zero,
            codomain: comodule,
            map: zero_map,
        }
    }

    fn compose(l: &Self, r: &Self) -> Self {
        debug_assert!(
            l.domain == r.codomain,
            "l-domain and r-codomain should be equal when composing"
        );

        let codomain = l.codomain.clone();
        let domain = r.domain.clone();

        let map = l.map.compose(&r.map);

        Self {
            domain,
            codomain,
            map,
        }
    }

    fn get_codomain(&self) -> Arc<RComodule<G,F>> {
        self.codomain.clone()
    }

    // (generated index , Grade , real index )
    fn get_structure_lines(&self) -> Vec<((usize, G, usize), (usize, G, usize), UniPolRing<F>, String)> {
        let mut lines = vec![];

        for (gr, gr_map) in self.map.maps.iter() {
            for el_id in 0..self.domain.space.dimension_in_grade(gr) {
                let els = self.domain.space.0.get(gr).unwrap();
                match els[el_id].0.primitive {
                    Some(prim_id) => {
                        for t_id in 0..gr_map.codomain() {
                            let t_el = &self.codomain.space.0.get(gr).expect("As codomain of the map is non-zero this vector space should contain an element in this grade.")[t_id];
                            if t_el.0.generator {
                                if !gr_map.get(el_id, t_id).is_zero() {
                                    let s_gen_id = els[el_id].0.generated_index;
                                    // find source generator
                                    let (s_gr, s_id) = {
                                        let mut a = None;
                                        'outer: for (gr, module) in &self.domain.space.0 {
                                            for (id,el) in module.iter().enumerate() {
                                                if el.0.generator && el.0.generated_index == s_gen_id {
                                                    a = Some((*gr, id));
                                                    break 'outer;
                                                }
                                            }
                                        }
                                        match a {
                                            Some(b) => b,
                                            None => {panic!()},
                                        }
                                    };

                                    lines.push((
                                        (s_gen_id, s_gr, s_id),
                                        (t_el.0.generated_index, *gr, t_id),
                                        gr_map.get(el_id, t_id),
                                        "h_".to_string() + &prim_id.to_string(),
                                    ));
                                }
                            }
                        }
                    }
                    None => {}
                }
            }
        }

        lines
    }
}