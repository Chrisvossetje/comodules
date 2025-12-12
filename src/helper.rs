use std::convert::identity;

use ahash::HashMap;

use crate::grading::Grading;

pub fn hashmap_add_restrict<G: Grading, V: Clone>(
    map: &HashMap<G, V>,
    grade: G,
    limit: G,
) -> HashMap<G, V> {
    hashmap_add_restrict_transform(map, grade, limit, identity)
}

pub fn hashmap_add_restrict_transform<G: Grading, V: Clone, W, T: Fn(V) -> W>(
    map: &HashMap<G, V>,
    grade: G,
    limit: G,
    f: T,
) -> HashMap<G, W> {
    map.iter()
        .filter_map(|(&g, v)| {
            let sum = g + grade;
            if sum <= limit {
                Some((g + grade, f(v.clone())))
            } else {
                None
            }
        })
        .collect()
}
