// struct PartialModule()

use std::sync::{Mutex, OnceLock};

use ahash::HashMap;
use algebra::{abelian::Abelian, field::Field, rings::univariate_polynomial_ring::UniPolRing};
use rayon::Scope;

use crate::{
    grading::grading::{GradedIndexing, Grading},
    k_comodule::{kcoalgebra::kCoalgebra, kmorphism::kComoduleMorphism},
    k_t_comodule::{k_t_coalgebra::ktCoalgebra, k_t_morphism::ktComoduleMorphism},
    resolution::{cokercell::CokerCell, datacell::DataCell},
    traits::{Coalgebra, ComoduleMorphism},
    types::ComoduleIndexType,
};

#[derive(Default, PartialEq, Debug)]
pub enum ResolveStateEnum {
    #[default]
    Nothing,
    Ready,
    Done,
}

#[derive(Default, PartialEq, Debug)]
// Cannot be that GS is Done and cokernel is NOT
pub struct ResolveState {
    cokernel: ResolveStateEnum,
    gs: ResolveStateEnum,
    degree_count: usize,
}

#[derive(Debug)]
pub struct ParallelResolution<G: Grading, C: Coalgebra<G>> {
    // First vec is s
    // G::ContiguousMemory is the grading, t in the unigrading case or u,v in the bigraded case
    // grid: Mutex<Vec<G::ContiguousMemory<State>>>,
    pub coalgebra: C,
    pub comodule: C::Comod,
    pub max_s: usize,
    pub max_degree: G,
    // All happing at a certrain grade g,
    // A \otimes V_i-1 -> Q, Q -> A\otimes V_i, V_i, A\otimes V_i, Hashmap to find the index in the list for A\otimes V_i
    pub data: Vec<
        G::ContiguousMemory<(
            Mutex<ResolveState>,
            OnceLock<CokerCell<G, C>>,
            OnceLock<DataCell<G, C>>,
        )>,
    >,
}

impl<G: Grading, C: Coalgebra<G>> ParallelResolution<G, C> {
    pub fn init(coalgebra: C, comodule: C::Comod, s: usize, degree: G) -> Self {
        let mut data = vec![];
        for _ in 0..=s {
            data.push(degree.init_memory(|| {
                (
                    Mutex::new(ResolveState::default()),
                    OnceLock::new(),
                    OnceLock::new(),
                )
            }));
        }

        for g in degree.iterator_from_zero(true) {
            for t in 0..=s {
                let inc = g.incomings();
                data[t].get(g.to_index()).0.lock().unwrap().degree_count = g.incomings();
                if inc == 0 {
                    data[t].get(g.to_index()).0.lock().unwrap().gs = ResolveStateEnum::Ready;
                }
            }
        }

        Self {
            coalgebra,
            comodule,
            data,
            max_s: s,
            max_degree: degree,
        }
    }

    pub fn get_data_cell(&self, s: usize, g: G) -> &DataCell<G, C> {
        self.data[s].get(g.to_index()).2.get().unwrap()
    }

    pub fn get_coker_cell(&self, s: usize, g: G) -> &CokerCell<G, C> {
        self.data[s].get(g.to_index()).1.get().unwrap()
    }

    pub fn resolve_coker_at_s_g(&self, s: usize, degree: G) {
        if s == 0 {
            return;
        }
        println!("Calculating Cokernel for {s} | {degree}");

        let prev_s = self.data[s - 1].get(degree.to_index()).2.get_or_init(|| {
            unreachable!();
        });
        let a = CokerCell::cokernel(prev_s);

        // TODO : make some sort "save" function
        self.data[s].get(degree.to_index()).1.set(a).unwrap();
        let mut m = self.data[s].get(degree.to_index()).0.lock().unwrap();
        m.cokernel = ResolveStateEnum::Done;

        println!("Finshed Cokernel for {s} | {degree}");
    }

    pub fn resolve_data_at_s_g(&self, s: usize, degree: G) {
        if s == 0 {
            return;
        }
        println!("Calculating Data for {s} | {degree}");

        let prev_s = self.data[s - 1].get(degree.to_index()).2.get_or_init(|| {
            unreachable!();
        });
        let coker = self.data[s].get(degree.to_index()).1.get_or_init(|| {
            unreachable!();
        });
        let a = DataCell::resolve(&self.data[s], prev_s, coker, degree, &self.coalgebra);

        // TODO : make some sort "save" function
        self.data[s].get(degree.to_index()).2.set(a).unwrap();
        let mut m = self.data[s].get(degree.to_index()).0.lock().unwrap();
        m.gs = ResolveStateEnum::Done;

        println!("Finshed Data for {s} | {degree}");
    }

    // This should initially be called on all cokernel ready elements
    pub fn recursion_solve<'a>(&'a self, i: &Scope<'a>) {
        for s in 0..=self.max_s {
            for degree in self.max_degree.iterator_from_zero(true) {
                let l = self.data[s].get(degree.to_index()).0.lock().unwrap();
                if l.cokernel == ResolveStateEnum::Ready {
                    i.spawn(move |i| {
                        self.resolve_coker_at_s_g(s, degree);
                        self.coker_spawn_task(i, s, degree);
                    });
                }
            }
        }
    }

    pub fn coker_spawn_task<'a>(&'a self, i: &Scope<'a>, s: usize, degree: G) {
        let m = self.data[s].get(degree.to_index()).0.lock().unwrap();
        match m.gs {
            ResolveStateEnum::Nothing => {}
            ResolveStateEnum::Ready => {
                i.spawn(move |i| {
                    self.resolve_data_at_s_g(s, degree);
                    self.data_spawn_task(i, s, degree);
                });
            }
            ResolveStateEnum::Done => {
                panic!("This could not have been done yet ?")
            }
        }
    }

    pub fn data_spawn_task<'a>(&'a self, i: &Scope<'a>, s: usize, degree: G) {
        if s + 1 <= self.max_s {
            let mut l = self.data[s + 1].get(degree.to_index()).0.lock().unwrap();
            match l.cokernel {
                ResolveStateEnum::Nothing => {
                    l.cokernel = ResolveStateEnum::Ready;
                    i.spawn(move |i| {
                        self.resolve_coker_at_s_g(s + 1, degree);
                        self.coker_spawn_task(i, s + 1, degree);
                    });
                }
                ResolveStateEnum::Ready => {
                    panic!("This could not be ready yet")
                }
                ResolveStateEnum::Done => {
                    panic!("This could not be finished yet")
                }
            }
        }

        for n in degree.nexts() {
            if n <= self.max_degree {
                let mut m = self.data[s].get(n.to_index()).0.lock().unwrap();
                m.degree_count -= 1;
                if m.degree_count == 0 {
                    m.gs = ResolveStateEnum::Ready;
                    match m.cokernel {
                        ResolveStateEnum::Nothing => {}
                        ResolveStateEnum::Ready => {}
                        ResolveStateEnum::Done => {
                            i.spawn(move |i| {
                                self.resolve_data_at_s_g(s, n);
                                self.data_spawn_task(i, s, n);
                            });
                        }
                    }
                }
            }
        }
    }
}

impl<G: Grading, F: Field, M: Abelian<F>> ParallelResolution<G, kCoalgebra<G, F, M>> {
    pub fn populate(&self) {
        let (map, cofree_comod) = kComoduleMorphism::inject_codomain_to_cofree(
            &self.coalgebra,
            &self.comodule,
            self.max_degree,
        );

        for g in self.max_degree.iterator_from_zero(true) {
            let r_gens = cofree_comod
                .space
                .0
                .get(&g)
                .map(|x| x.clone())
                .unwrap_or(vec![]);
            let transposed_to_cofree = map
                .map
                .maps
                .get(&g)
                .map(|x| x.clone())
                .unwrap_or(M::zero(0, 0))
                .transpose();

            let mut lut = HashMap::default();
            for (id, r) in r_gens.iter().enumerate() {
                lut.insert(r.0, id as ComoduleIndexType);
            }

            let cokernel = self
                .comodule
                .space
                .0
                .get(&g)
                .map(|x| x.clone())
                .unwrap_or(vec![]);
            let to_cokernel = M::zero(0, cokernel.len());

            let a_gens = r_gens
                .iter()
                .enumerate()
                .filter_map(|(id, x)| {
                    if x.0.0.0 == G::zero() {
                        Some((x.1, x.0.1, id as ComoduleIndexType))
                    } else {
                        None
                    }
                })
                .collect();

            let lut2 = HashMap::default();

            let a = DataCell {
                transposed_to_cofree,
                a_gens,
                r_gens,
                lut,
                lut2,
            };
            let b = CokerCell {
                to_cokernel,
                cokernel,
                repr_vecs: M::zero(0, 0), // We shouldn't need this anymore ?
            };
            self.data[0].get(g.to_index()).1.set(b).unwrap();
            self.data[0].get(g.to_index()).2.set(a).unwrap();

            let mut l = self.data[0].get(g.to_index()).0.lock().unwrap();
            l.cokernel = ResolveStateEnum::Done;
            l.gs = ResolveStateEnum::Done;
            l.degree_count = 0;
        }

        for h in self.max_degree.iterator_from_zero(true) {
            let mut data = self.data[1].get(h.to_index()).0.lock().unwrap();
            data.cokernel = ResolveStateEnum::Ready;
        }
    }

    // pub fn populate_with_basering(&self) {
    //     // s = 0 should be the coalgebra itself
    //     // Non zero grade is just the
    //     for g in self.max_degree.iterator_from_zero(true) {
    //         if g == G::zero() {
    //             continue;
    //         }
    //         let r_size = kCoalgebra::size_in_degree(&self.coalgebra, g);

    //         let r_gens: Vec<_> = (0..r_size)
    //             .map(|x| (((g, x as CoalgebraIndexType), 0), M::Generator::default()))
    //             .collect();

    //         let to_cokernel = M::zero(0, 0);
    //         let to_cofree = M::zero(0, r_gens.len());

    //         let mut lut = HashMap::default();
    //         for (r_id, r) in r_gens.iter().enumerate() {
    //             lut.insert(r.0, r_id as ComoduleIndexType);
    //         }

    //         let a = DataCell {
    //             to_cofree,
    //             a_gens: vec![],
    //             r_gens,
    //             lut,
    //             lut2: HashMap::default(),
    //         };
    //         self.data[0].get(g.to_index()).1.set(a).unwrap();
    //         let mut data = self.data[0].get(g.to_index()).0.lock().unwrap();
    //         *data = ResolveStateEnum::Ready;
    //     }

    //     // Grade zero
    //     let r_gens: Vec<_> = vec![(((G::zero(), 0), 0), M::Generator::default())];
    //     let to_cokernel = M::zero(0, 1);
    //     let to_cofree = M::identity(1);
    //     let cokernel = vec![M::Generator::default()];
    //     let a_gens = vec![(M::Generator::default(), 0, 0)];

    //     let mut lut = HashMap::default();
    //     for (r_id, r) in r_gens.iter().enumerate() {
    //         lut.insert(r.0, r_id as u32);
    //     }

    //     let a = DataCell {
    //         to_cokernel,
    //         to_cofree,
    //         cokernel,
    //         a_gens,
    //         r_gens,
    //         lut,
    //         lut2: HashMap::default(),
    //     };
    //     self.data[0].get(G::zero().to_index()).1.set(a).unwrap();
    //     let mut data = self.data[0].get(G::zero().to_index()).0.lock().unwrap();
    //     *data = ResolveStateEnum::Ready;

    //     for i in 1..=self.max_s {
    //         let mut data = self.data[i].get(G::zero().to_index()).0.lock().unwrap();
    //         *data = ResolveStateEnum::GsAvailable;
    //     }

    //     for h in self.max_degree.iterator_from_zero(true) {
    //         let mut data = self.data[1].get(h.to_index()).0.lock().unwrap();
    //         *data = ResolveStateEnum::CokernelAvailable;
    //     }

    //     let mut data = self.data[1].get(G::zero().to_index()).0.lock().unwrap();
    //     *data = ResolveStateEnum::Ready;
    // }
}

impl<G: Grading, F: Field, M: Abelian<UniPolRing<F>>> ParallelResolution<G, ktCoalgebra<G, F, M>> {
    pub fn populate(&self) {
        let (map, cofree_comod) = ktComoduleMorphism::inject_codomain_to_cofree(
            &self.coalgebra,
            &self.comodule,
            self.max_degree,
        );

        for g in self.max_degree.iterator_from_zero(true) {
            let r_gens = cofree_comod
                .space
                .0
                .get(&g)
                .map(|x| x.clone())
                .unwrap_or(vec![]);
            let transposed_to_cofree = map
                .map
                .maps
                .get(&g)
                .map(|x| x.clone())
                .unwrap_or(M::zero(0, 0))
                .transpose();

            let mut lut = HashMap::default();
            for (id, r) in r_gens.iter().enumerate() {
                lut.insert(r.0, id as ComoduleIndexType);
            }

            let cokernel = self
                .comodule
                .space
                .0
                .get(&g)
                .map(|x| x.clone())
                .unwrap_or(vec![]);
            let to_cokernel = M::zero(0, cokernel.len());

            let a_gens = r_gens
                .iter()
                .enumerate()
                .filter_map(|(id, x)| {
                    if x.0.0.0 == G::zero() {
                        Some((x.1, x.0.1, id as ComoduleIndexType))
                    } else {
                        None
                    }
                })
                .collect();

            let lut2 = HashMap::default();

            let a = DataCell {
                transposed_to_cofree,
                a_gens,
                r_gens,
                lut,
                lut2,
            };
            let b = CokerCell {
                to_cokernel,
                cokernel,
                repr_vecs: M::zero(0, 0), // We shouldn't need this anymore ?
            };
            self.data[0].get(g.to_index()).1.set(b).unwrap();
            self.data[0].get(g.to_index()).2.set(a).unwrap();

            let mut l = self.data[0].get(g.to_index()).0.lock().unwrap();
            l.cokernel = ResolveStateEnum::Done;
            l.gs = ResolveStateEnum::Done;
            l.degree_count = 0;
        }

        for h in self.max_degree.iterator_from_zero(true) {
            let mut data = self.data[1].get(h.to_index()).0.lock().unwrap();
            data.cokernel = ResolveStateEnum::Ready;
        }
    }

    // pub fn populate_with_basering(&self) {
    //     // s = 0 should be the coalgebra itself
    //     // Non zero grade is just the
    //     for g in self.max_degree.iterator_from_zero(true) {
    //         if g == G::zero() {
    //             continue;
    //         }
    //         let r_size = ktCoalgebra::size_in_degree(&self.coalgebra, g);

    //         let r_gens: Vec<_> = (0..r_size)
    //             .map(|x| (((g, x as CoalgebraIndexType), 0), M::Generator::default()))
    //             .collect();

    //         let to_cokernel = M::zero(0, 0);
    //         let to_cofree = M::zero(0, r_gens.len());

    //         let mut lut = HashMap::default();
    //         for (r_id, r) in r_gens.iter().enumerate() {
    //             lut.insert(r.0, r_id as ComoduleIndexType);
    //         }

    //         let a = DataCell {
    //             to_cokernel,
    //             to_cofree,
    //             cokernel: vec![],
    //             a_gens: vec![],
    //             r_gens,
    //             lut,
    //             lut2: HashMap::default(),
    //         };
    //         self.data[0].get(g.to_index()).1.set(a).unwrap();
    //         let mut data = self.data[0].get(g.to_index()).0.lock().unwrap();
    //         *data = ResolveStateEnum::Ready;
    //     }

    //     // Grade zero
    //     let r_gens: Vec<_> = vec![(((G::zero(), 0), 0), M::Generator::default())];
    //     let to_cokernel = M::zero(0, 1);
    //     let to_cofree = M::identity(1);
    //     let cokernel = vec![M::Generator::default()];
    //     let a_gens = vec![(M::Generator::default(), 0, 0)];

    //     let mut lut = HashMap::default();
    //     for (r_id, r) in r_gens.iter().enumerate() {
    //         lut.insert(r.0, r_id as u32);
    //     }

    //     let a = DataCell {
    //         to_cokernel,
    //         to_cofree,
    //         cokernel,
    //         a_gens,
    //         r_gens,
    //         lut,
    //         lut2: HashMap::default(),
    //     };
    //     self.data[0].get(G::zero().to_index()).1.set(a).unwrap();
    //     let mut data = self.data[0].get(G::zero().to_index()).0.lock().unwrap();
    //     *data = ResolveStateEnum::Ready;

    //     for i in 1..=self.max_s {
    //         let mut data = self.data[i].get(G::zero().to_index()).0.lock().unwrap();
    //         *data = ResolveStateEnum::GsAvailable;
    //     }

    //     for h in self.max_degree.iterator_from_zero(true) {
    //         let mut data = self.data[1].get(h.to_index()).0.lock().unwrap();
    //         *data = ResolveStateEnum::CokernelAvailable;
    //     }

    //     let mut data = self.data[1].get(G::zero().to_index()).0.lock().unwrap();
    //     *data = ResolveStateEnum::Ready;
    // }
}
