use std::{collections::HashMap, hash::Hash, ops::Add, sync::Arc};

use crate::linalg::{field::Field, graded::{BasisElement, BasisIndex, GradedLinearMap, GradedVectorSpace, Grading}, matrix::{FieldMatrix, Matrix}};

use super::comodule::{Comodule, ComoduleMorphism, Tensor};


#[derive(Debug, Clone)]
pub struct kBasisElement {
    pub name: String,
    pub generator: bool,
    pub primitive: Option<usize>,
    pub generated_index: usize,
}

impl BasisElement for kBasisElement {

}


#[derive(Debug, Clone)]
pub struct kCoalgebra<G: Grading, F: Field> {
    pub space: GradedVectorSpace<G, kBasisElement>,
    pub coaction: GradedLinearMap<G, F, FieldMatrix<F>>,
    pub tensor: kTensor<G>,
} 


#[derive(Debug, Clone)]
pub struct kComodule<G: Grading, F: Field> {
    pub coalgebra: Arc<kCoalgebra<G, F>>,
    pub space: GradedVectorSpace<G, kBasisElement>,
    pub coaction: GradedLinearMap<G, F, FieldMatrix<F>>,
    pub tensor: kTensor<G>,
} 


#[derive(Debug, Clone)]
pub struct kTensor<G: Grading> {
    // # Module Grade + Index -> Algebra Grading + index -> Tensor Grading + index
    pub construct: HashMap<BasisIndex<G>, HashMap<BasisIndex<G>, BasisIndex<G>>>,


    // # Module Grade + Index -> Algebra Grading + index -> Tensor Grading + index
    pub deconstruct: HashMap<BasisIndex<G>,(BasisIndex<G>, BasisIndex<G>)>, 

    pub dimensions: HashMap<G, usize>,
}

// #[derive(Debug, Clone)]
// pub struct kTensorFor<G: Grading> {
//     construct: HashMap<BasisIndex<G>, HashMap<>>,
//     deconstruct: HashMap<BasisIndex<G>,(BasisIndex<G>, BasisIndex<G>)>, 
// }

#[derive(Debug, Clone)]
pub struct kComoduleMorphism<G: Grading, F: Field> {
    domain: Arc<kComodule<G, F>>,
    codomain: Arc<kComodule<G, F>>,

    map: GradedLinearMap<G,F, FieldMatrix<F>> // Question: Shouldn't this be a module morphism?
}


impl<G: Grading> Tensor<G> for kTensor<G> {
    fn tensor_to_base() {
        todo!()
    }

    fn base_to_tensor() {
        todo!()
    }
    
    fn get_dimension(&self, grading: &G) -> usize {
        self.dimensions[grading]
    }
}

impl<G: Grading> Default for kTensor<G> {
    fn default() -> Self {
        Self { construct: Default::default(), deconstruct: Default::default(), dimensions: Default::default() }
    }
}

impl<G: Grading> kTensor<G> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn generate<B: BasisElement>(left: &GradedVectorSpace<G, B>, right: &GradedVectorSpace<G, B>) -> Self {
        let mut construct = HashMap::new();
        let mut deconstruct = HashMap::new();
        let mut dimensions = HashMap::new();

        for (r_grade, r_elements) in right.0.iter() {
            for r_id in 0..r_elements.len() {

                let mut construct_map = HashMap::new();
                
                for (l_grade, l_elements) in left.0.iter() {
                    let t_grade = *l_grade + *r_grade;

                    if !right.0.contains_key(&t_grade) {
                        continue; 
                    }

                    for l_id in 0..l_elements.len() {
                        
                        
                        let t_id = dimensions.entry(t_grade).or_insert(0);
                        construct_map.insert((*l_grade, l_id), (t_grade, *t_id as usize));

                        deconstruct.insert((t_grade, *t_id as usize), ((*l_grade, l_id),(*r_grade, r_id)));
                        *t_id = *t_id + 1;
                    }
                }
                construct.insert((*r_grade, r_id), construct_map);
            }
        }
        

        Self {
            construct,
            deconstruct,
            dimensions,
        }
    }
}



impl<G: Grading, F: Field> Comodule<G> for kComodule<G, F> {
    type Element = kBasisElement;
    type Coalgebra = kCoalgebra<G,F>;

    fn get_generators(&self) -> Vec<(usize, G, Option<String>)> {
        self.space.0.iter().flat_map(|(k,v)| {
            v.iter().filter_map(|b| {
                if b.generator {
                    Some((b.generated_index, *k, Some(b.name.clone())))
                } else {
                    None
                }
            })
        }).collect()
    }
    
    
    fn zero_comodule(coalgebra: Arc<Self::Coalgebra>) -> Self {
        Self {
            coalgebra: coalgebra,
            space: GradedVectorSpace::new(),
            coaction: GradedLinearMap::empty(),
            tensor: kTensor::new()
        }
    }
    
    fn fp_comodule(coalgebra: Arc<Self::Coalgebra>) -> Self {
        let zero = G::zero();

        // (zero, 
        let el = kBasisElement { name: "fp".to_string(), generator: false, primitive: None, generated_index: 0};


        let space_map: HashMap<G, Vec<kBasisElement>> = [(zero,vec![el])].into_iter().collect();
        let space = GradedVectorSpace::from(space_map);


        let coact_map: HashMap<G, FieldMatrix<F>> = [(zero, FieldMatrix::identity(1))].into_iter().collect();
        let coaction = GradedLinearMap::from(coact_map);

        assert_eq!(coalgebra.space.0.get(&zero).expect("Coalgebra has no element in grade zero").len(), 1, "Coalgebra is not a connected coalgebra");

        let mut dimensions = HashMap::new();
        dimensions.insert(zero, 1);
        
        let mut construct = HashMap::new();
        let mut first_entry = HashMap::new();
        first_entry.insert((zero,0), (zero, 0));
        construct.insert((zero, 0), first_entry);

        let mut deconstruct = HashMap::new();
        deconstruct.insert((zero,0), ((zero,0), (zero,0)));

        let tensor = kTensor { construct, deconstruct, dimensions };

        Self { coalgebra, space, coaction, tensor }   
    }
}


impl<G: Grading, F: Field> ComoduleMorphism<G, kComodule<G, F>> for kComoduleMorphism<G, F> {
    fn cokernel(&self) -> Self {
        let cokernel_map = self.map.get_cokernel();
    
        let coker_space = cokernel_map.codomain_space(kBasisElement { 
            name: "".to_string(), generator: false, primitive: None, generated_index: 0 
        });

        
        let coalg = self.codomain.coalgebra.as_ref();
        let tensor = kTensor::generate(&coalg.space, &coker_space);
        
        dbg!(&tensor);

        let pivots = cokernel_map.pivots();
        let coaction: HashMap<G, FieldMatrix<F>> = coker_space.0.iter().map(|(g,v)| {
            let g_tensor_dimen = tensor.get_dimension(g);
            let mut g_coaction = FieldMatrix::<F>::zero(v.len(), g_tensor_dimen);

            // let codom = self.codomain.clone();

            // TODO:
            // THERE MIGHT BE A FASTER VERSION OF THIS, BUT THIS IS SIMPLER
            // # Slightly slower for large coalgebras

            // (Row, Column)
            let g_map = cokernel_map.maps.get(g).expect("cokernel map should exist in this grade");
            for (coker_id, codom_id) in &pivots[g] {
                // F_coact_vec
                let coact_size = self.codomain.tensor.dimensions[g];
                for codom_coact_id in 0..coact_size {
                    let coact_val = self.codomain.coaction.maps[g].data[codom_coact_id][*codom_id];
                    if !coact_val.is_zero() {
                        let ((alg_gr, alg_id), (mod_gr, mod_id)) = self.codomain.tensor.deconstruct[&(*g, codom_coact_id)];
                        
                        for target_id in 0..g_map.domain {
                            let m = g_map.data[target_id][mod_id];
                            let (_, final_id) = tensor.construct[&(mod_gr, target_id)][&(alg_gr, alg_id)];
                            g_coaction.data[final_id][*coker_id] += coact_val * m;    
                        }
                            
                    }
                }
            }
            

            // #     for Q_index, F_index in enumerate(pivots[Q_grade]):
            // #         F_coact_vec = F.codomain.coaction[Q_grade][:,F_index]
            // #         for coact_id in range(F_coact_vec.dimensions()[0]):
            // #             if F_coact_vec[coact_id]:
            // #                 (alg_gr, alg_id), (mod_gr, mod_id) = F.codomain.tensored[Q_grade][coact_id]
            // #                 target = M[mod_gr][:,mod_id]
        
            // #                 for (target_id, value) in lut[mod_gr][mod_id]:
            // #                         tensor_gr, tensor_id = Q_moduled[mod_gr][target_id][alg_gr][alg_id]
            // #                         if globals.TEST:
            // #                             assert tensor_gr == Q_grade, "Resulting tensor grade not equal to original grade"
            // #                         coaction[tensor_gr][tensor_id,Q_index] += target[target_id] * F_coact_vec[coact_id]
            

            (*g, g_coaction)
        }).collect();


        let comodule = kComodule {
            coalgebra: self.codomain.coalgebra.clone(),
            space: coker_space,
            coaction: GradedLinearMap::from(coaction),
            tensor,
        };
        
        Self {
            codomain: Arc::new(comodule),
            domain: self.codomain.clone(),
            map: cokernel_map, 
        }
    }

    fn inject_codomain_to_cofree(&self) -> Self {
        todo!()
    }

    fn zero_morphism(comodule: Arc<kComodule<G, F>>) -> Self {
        let codomain = comodule.clone();
        let zero = Arc::new(kComodule::zero_comodule(comodule.coalgebra.clone()));

        // Verify how we want to handle this zero map
        let mut zero_map = GradedLinearMap::empty();

        for (gr, elements) in codomain.space.0.iter() {
            zero_map.maps.insert(*gr, FieldMatrix::zero(0, elements.len()));
        }

        Self {
            domain: zero,
            codomain: codomain,
            map: zero_map,
        }
    }

    fn compose(l: Self, r: Self) -> Self {
        let codomain = l.codomain;
        let domain = r.domain;

        let map = l.map.compose(r.map);

        Self {
            domain,
            codomain,
            map,
        }
    }
    
    fn get_structure_lines(&self) -> Vec<(BasisIndex<G>, BasisIndex<G>, usize)> {
        todo!()
    }
    
    fn get_codomain(&self) -> Arc<kComodule<G, F>> {
        self.codomain.clone()
    }
}


