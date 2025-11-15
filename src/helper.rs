use std::convert::identity;

use ahash::HashMap;

use crate::{grading::Grading, linalg::{flat_matrix::FlatMatrix, graded::BasisIndex, matrix::RModMorphism, ring::CRing}};


pub fn hashmap_add_restrict<G: Grading, V: Clone>(map: &HashMap<G, V>, grade: G, limit: G) -> HashMap<G, V> {
    hashmap_add_restrict_transform(map, grade, limit, identity)
}

pub fn hashmap_add_restrict_transform<G: Grading, V: Clone, W, T: Fn(V) -> W>(map: &HashMap<G, V>, grade: G, limit: G, f: T) -> HashMap<G, W> {
    map.iter().filter_map(|(&g, v)| {
        let sum = g + grade;
        if sum <= limit {
            Some((g + grade, f(v.clone())))
        } else {
            None
        }
    })
    .collect()
}


// TODO ! REMOVE
pub fn add_or_find_row<G: Grading, R: CRing>(mat: &mut FlatMatrix<R>, exp: &mut Vec<((G, usize), usize)>, domain: usize, element: R, grade: BasisIndex<G>, power: usize) {
    let (final_gr, final_id) = grade;
    
    // Maybe we have to check the module structure here as well :(
    // The module structure of the tensor, as that could dissappear
    // AHHHHHHHHHHH
    let mut found = false;
    for (id, j) in exp.iter().enumerate() {
        if (final_gr, final_id) == j.0 {
            found = true;
            mat.add_at(domain, id, element);
        }
    }
    // And we should compensate for some target_power / inv_power stuff
    // TFFOOOOOEEEE
    if !found {
        mat.extend_one_row();
        mat.add_at(domain, exp.len(), element);
        exp.push(((final_gr, final_id), power));
    }
}