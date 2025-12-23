use crate::{field::Field, matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::{CRing, ValuationRing}, rings::univariate_polynomial_ring::UniPolRing, snf::SmithNormalForm, unipol::{UniPolMap, UniPolModule, morphism::order_maps}};

/// M -f> N -g> Q
pub(super) fn internal_cohomology<F: Field>(f: &UniPolMap<F>, g: &UniPolMap<F>, n: &UniPolModule, q: &UniPolModule) -> (UniPolMap<F>, UniPolModule) {
    debug_assert_eq!(f.codomain, g.domain);

    if cfg!(debug_assertions) {
        let mut comp = g.compose(&f);
        comp.reduce(&q);
        for x in 0..comp.domain {
            for y in 0..comp.codomain {
                assert!(comp.get(x, y).is_zero());
            }
        }
    }

    // println!("f:\n{:?}g:\n{:?}n: {:?}\nq: {:?}\n", f, g, n,q);

    let (new_f,new_g, new_n, new_q, trans_map_inv) = order_maps(f, g, n, q);
    let (f,g,n,q) = (&new_f, &new_g, &new_n, &new_q);
    
    // println!("f:\n{:?}g:\n{:?}n: {:?}\nq: {:?}\n", f, g, n,q);


    let new_domain = g.domain + g.codomain;
    let mut g_aug = FlatMatrix::zero(new_domain, g.codomain);

    for i in 0..g_aug.codomain {
        g_aug.set_row(i, g.get_row(i));
        
        
        if let Some(power) = q[i] {
            g_aug.set(g.domain + i, i, UniPolRing(F::one().neg(),power));
        }
    }


    let (_,s_aug,v_aug) = g_aug.snf();

    let mut start_zeros = s_aug.codomain;
    // Assume that s is ordered in some way (nonzero first)
    for r in 0..s_aug.codomain {
        let el = s_aug.get(r, r);
        if el.is_zero() {
            start_zeros = r;
            break;
        }
    }

    // println!("g_aug:\n{:?}s_aug:\n{:?}v_aug:\n{:?}", g_aug, s_aug, v_aug);

    let g_ker_size = s_aug.domain - start_zeros;
    let mut g_ker = FlatMatrix::zero(g_ker_size, g.domain);

    
    // g_ker has as columns the vectors which generate the kernel
    for r in 0..g.domain {
        let row = &v_aug.get_row(r)[start_zeros..s_aug.domain];
        g_ker.set_row(r, row);
    }

    
    // TODO : THIS IS OPTIoNAL ??
    // reduce(&mut g_ker, &n);
    // println!("g_ker:\n{:?}", g_ker);
    
    
    let (_,s_ker,_, uinv_ker, _) = g_ker.full_snf();
    
    let mut non_zero_els = 0;
    for r in 0..s_ker.codomain {
        let el = s_ker.get(r, r);
        if el.is_zero() {
            break;
        }
        non_zero_els += 1;
    }

    // println!("s_ker:\n{:?}uinv_ker:\n{:?}v_ker:\n{:?}", s_ker, uinv_ker, v_ker);
    

    let mut vecs = vec![];
    for r in 0..non_zero_els {
        let vector = uinv_ker.get_column(r);
        let mult_vector: Vec<_> = vector.into_iter().map(|x| s_ker.get(r, r) * x).collect();
        vecs.push(mult_vector);
    }



    // TODO : this is ugly
    let mut real_ker_module_structure = vec![];
    vecs = vecs.into_iter().filter(|v| {
        let mut mod_str = Some(0);
        for y in 0..v.len() {
            let el = v[y];
            if el.is_zero() {
                continue;
            }

            let n_x_structure = n[y];
            match n_x_structure {
                Some(n_x_power) => {
                    let new_power = n_x_power - el.1; // If non-negative then the map was unreduced, probably ?

                    if let Some(p) = mod_str {
                        if p < new_power {
                            mod_str = Some(new_power);
                        }
                    }
                },
                None => {
                    mod_str = None;
                    break;
                },
            }
        }
        
        if mod_str == Some(0) {
            false 
        } else {
            real_ker_module_structure.push(mod_str);
            true
        }
    }).collect();

    let mut real_g_ker = FlatMatrix::zero(vecs.len(), g.domain);
    
    
    for (id, column) in vecs.into_iter().enumerate() {
        real_g_ker.set_column(id, &column[..]);
    }

    // TODO : OPTIONAL ?
    // reduce(&mut real_g_ker, n);
    // println!("real_g_ker:\n{:?}", real_g_ker);


    // Check if vectors are REALLY zero
    if cfg!(debug_assertions) {
        for a in 0..real_g_ker.domain {
            let column = real_g_ker.get_column(a);
            let eval = g.eval_vector(&column);

            println!("{:?}{:?}", eval, column);
            
            for (id,r) in eval.iter().enumerate() {
                let power = match q[id] {
                    Some(p) => p,
                    None => u16::MAX,
                };
                if !(r.is_zero() || r.1 >= power) {
                    panic!("OH OH, kernel is NOT zero")
                }
            }

        }
    }

    // println!("real_ker_module: {:?}\n\n", real_ker_module_structure);
    
    let (u_real_ker,s_real_ker,_) = real_g_ker.snf();
    let f_in_kernel = u_real_ker.compose(&f);
    
    // reduce(Y, codomain);
    // println!("u_real_ker:\n{:?}s_real_ker:\n{:?}Y:\n{:?}", u_real_ker,s_real_ker, y);
    

    let mut sol = FlatMatrix::zero(f.domain + real_ker_module_structure.len(), real_g_ker.domain);
    for x in 0..f.domain {
        for y in 0..sol.codomain {
            let y_i = f_in_kernel.get(x, y);
            if y_i.is_zero() {
                continue;
            }
            let d_i = s_real_ker.get(y, y);
            let z_i = y_i.unsafe_divide(d_i);
            sol.set(x, y, z_i);
        }
    }

    for (r,el) in real_ker_module_structure.iter().enumerate() {
        match el {
            Some(power) => {
                sol.set(f.domain + r, r, UniPolRing(F::one().neg(), *power));
                
            },
            None => {},
        }
    }

    
    let (_,sol_s,_, sol_uinv, _) = sol.full_snf();

    // println!("sol_s:\n{:?}sol_uinv:\n{:?}sol:\n{:?}", sol_s, sol_uinv, sol);

    let mut module = Vec::default();
    let mut columns = vec![];
    for r in 0..sol_s.codomain {
        let el = sol_s.get(r, r);
        if !el.is_unit() {
            let structure = {
                if el.is_zero() { None } else { Some(el.1) }
            };
            let col = sol_uinv.get_column(r);

            module.push(structure);
            columns.push(col);
        }
    }
    
    let mut cohom_to_ker = FlatMatrix::zero(module.len(), real_g_ker.domain);
    for (id, column) in columns.into_iter().enumerate() {
        cohom_to_ker.set_column(id, &column[..]);
    }

    let ker_to_n = trans_map_inv.compose(&real_g_ker);

    // println!("ker_to_n:\n{:?}\n\nker_to_cohom:\n{:?}", ker_to_n, cohom_to_ker);

    (ker_to_n.compose(&cohom_to_ker), module)
}
