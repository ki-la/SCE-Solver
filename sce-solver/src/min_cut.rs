use std::cmp::{max, min};
use ahash::AHashSet;
use indexmap::IndexSet;
use autocxx::prelude::*;
use crate::forb_ind_subgraph::{compute_packing, splittance_from_buckets};
use crate::graph_trait::Graph;

include_cpp! {
    #include "find_min_cut.h"
    safety!(unsafe_ffi)
    generate!("min_cut")
}

/// calculates and returns MinCut of graph 'g'
/// -> calls min_cut() in find_min_cut.h, which uses VieCut library
/// (https://github.com/VieCut/VieCut) to calculate the MinCut
pub fn min_cut(g: &impl Graph) -> usize {
    let mut starts:Vec<i64> = Vec::with_capacity(g.edge_count());
    let mut ends:Vec<i64> = Vec::with_capacity(g.edge_count());
    let mut weights: Vec<i64> = vec![1; g.edge_count()];
    for i in 0..g.size() {
        for j in i+1..g.size() {
            if g.has_edge(i, j) {
                starts.push(i.try_into().unwrap());
                ends.push(j.try_into().unwrap());
            }
        }
    }
    assert_eq!(starts.len(), g.edge_count());

    unsafe{ffi::min_cut(starts.as_mut_ptr(),
                        ends.as_mut_ptr(),
                        weights.as_mut_ptr(),
                        g.size().try_into().unwrap(),
                        g.edge_count().try_into().unwrap()) as usize}
}

/// creates induced subgraph of 'g' ignoring the vertices in 'ignore' and returns its MinCut
/// panics if 'm_new' != #edges in the resulting subgraph
fn min_cut_subgraph(g: &impl Graph, ignore: &AHashSet<usize>, m_new: usize) -> usize {
    let mut starts: Vec<i64> = Vec::with_capacity(m_new);
    let mut ends: Vec<i64> = Vec::with_capacity(m_new);
    // decrease vertex IDs s.t. vertices in resulting graph are labeled 0..n'-1
    // if vertex 'v' is ignored, then all higher vertex IDs are decreased by 1 to avoid a gap
    let mut dec_i_by = 0;
    for i in 0..g.size() {
        if ignore.contains(&i) {
            dec_i_by += 1;
            continue;
        }
        let mut dec_j_by = dec_i_by;
        for j in i+1..g.size() {
            if ignore.contains(&j) {
                dec_j_by += 1;
                continue;
            }
            if g.has_edge(i, j) {
                starts.push((i - dec_i_by).try_into().unwrap());
                ends.push((j - dec_j_by).try_into().unwrap());
            }
        }
    }
    let mut weights: Vec<i64> = vec![1; starts.len()];
    debug_assert_eq!(starts.len(), m_new);

    unsafe { ffi::min_cut(starts.as_mut_ptr(),
                          ends.as_mut_ptr(),
                          weights.as_mut_ptr(),
                          (g.size() - ignore.len()).try_into().unwrap(),
                          m_new.try_into().unwrap()).try_into().unwrap() }
}

/// calculates and returns MinCut Lower Bound for graph 'g'
/// repeatedly removes deg-1 vertices until none are left -> calls min_cut() on resulting graph
/// CAUTION: assumes that 'g' is connected (so, MinCut >= 1)!
pub fn min_cut_lower_bound(g: &impl Graph, s: Option<usize>) -> usize {
    let (mut deg_buckets, mut degrees, min_deg) = g.degree_sequence_buckets(None);
    assert!(min_deg > 0);  // no isolated nodes!

    let split = match s {
        Some(spl) => spl,
        None => splittance_from_buckets(&deg_buckets),
    };

    if split < 2 {  // MinCut >= 1 >= splittance
        return min(split, 1);
    }

    if min_deg > 1 {  // no deg-1 vertices -> return minimum of MinCut and splittance
        return min(min_cut(g), split);
    }

    let mut ignore: AHashSet<usize> = AHashSet::new();  // vertices ignored in MinCut calculation!
    while !deg_buckets[1].is_empty() {
        let mut v_curr = deg_buckets[1].pop().unwrap();
        ignore.insert(v_curr);
        // ignoring 'v' might create new deg-1 vertices (later) -> remove them as well:
        let mut new_deg1_vtx = true;
        while new_deg1_vtx && ignore.len() < g.size() {
            for nb in g.neighbors(v_curr, None) {
                if ignore.contains(&nb) {
                    continue;
                }
                // decrease 'nb's degree -> move to preceding bucket
                deg_buckets[degrees[nb]].remove(&nb);  // O(1) on avg.
                degrees[nb] -= 1;
                deg_buckets[degrees[nb]].insert(nb);  // O(1) on avg.
                if degrees[nb] <= 1 {
                    ignore.insert(nb);
                    deg_buckets[degrees[nb]].remove(&nb);
                    v_curr = nb;
                } else {
                    new_deg1_vtx = false;
                }
                break;
            }
        }
    }
    // CAUTION: 'degrees' was not updated for removed vertices -> don't use it!

    // now, either no vertices or at least three are left!
    if g.size() == ignore.len() {
        return 1;  // the input graph is a tree (-> can be solved in polynomial time)
    }
    assert!(g.size() - ignore.len() > 2);
    if g.size() - ignore.len() == 3 {
        return 1;  // g is a triangle + forest -> at least one cut necessary (as original splittance > 0)
    }

    // compute MinCut on subgraph induced by V \ ignore (-> has ignore.len() edges less, as only deg-1 vertices are ignored)
    let mc = min_cut_subgraph(g, &ignore, g.edge_count() - ignore.len());
    // compute splittance on subgraph induced by V \ ignore (deg_buckets has been updated!)
    let s = splittance_from_buckets(&deg_buckets);
    if s == 0 {
        // splittance(original g) >= 2 & at least one deg-1 vertex -> at least one edit necessary
        // => ignore.len() gives an upper bound as graph is split after removing these vertices
        return 1;
    }
    min(mc, s)
}

/// calculates and returns MinCut Lower Bound for graph 'g'
/// repeatedly removes low degree vertices until mc > split. -> calls min_cut() on resulting graph
/// CAUTION: assumes that 'g' is connected (so, MinCut >= 1)!
pub fn min_cut_lower_bound2(g: &impl Graph, s: Option<usize>, comp: Option<&AHashSet<usize>>, ub: Option<usize>) -> usize {
    let (mut deg_buckets, mut degrees, min_deg) = g.degree_sequence_buckets(comp);
    assert!(min_deg > 0);  // no isolated nodes!

    let mut split = match s {
        Some(spl) => spl,
        None => splittance_from_buckets(&deg_buckets),
    };

    if split < 2 {  // MinCut >= 1 >= splittance
        return min(split, 1);
    }

    let mut rem_edges = 0;
    let mut ignore: AHashSet<usize> = AHashSet::with_capacity(g.size());  // vertices ignored in MinCut calculation!
    if let Some(c) = comp {
        // if comp.is_some(), then ignore all vertices outside of the component
        ignore = (0..g.size()).collect::<AHashSet<usize>>().difference(c).map(|x| *x).collect::<AHashSet<usize>>();
        for u in &ignore {
            for v in &ignore {
                if *u < *v && g.has_edge(*u, *v) {
                    rem_edges += 1;
                }
            }
        }
    }

    let mut mc: usize = if min_deg == 1 {
        1
    } else {
        min_cut_subgraph(g, &ignore, g.edge_count() - rem_edges)
    };
    if mc >= split {  // MinCut gives no LB -> splittance gives opt. solution
        return split;
    }
    let mut lb = mc;
    let upper_bound = match ub {
        Some(u) => u,  // ub given => could allow to break look early
        None => split,       // ub not given => splittance is an upper bound
    };

    for deg in min_deg..deg_buckets.len()-1 {  // loop over degrees from smallest to (second) largest
        // remove all vertices w/ degree <= 'deg' (ignore.len() vertices have been removed already)
        if deg_buckets[deg].is_empty() {
            continue;  // no vertex with degree 'deg' -> skip
        }
        while !deg_buckets[deg].is_empty() {
            let mut v_curr = deg_buckets[deg].pop().unwrap();
            ignore.insert(v_curr);
            // "removing" v_curr might create new low degree vertices (later) -> remove them as well:
            let mut to_check = IndexSet::from([v_curr]);
            while !to_check.is_empty() && ignore.len() < g.size() {
                v_curr = to_check.pop().unwrap();
                for nb in g.neighbors_iter(v_curr) {
                    if ignore.contains(&nb) {
                        if to_check.contains(&nb) {
                            rem_edges += 1;
                        }
                        continue;
                    }
                    // decrease 'nb's degree -> move to preceding bucket
                    deg_buckets[degrees[nb]].remove(&nb);
                    degrees[nb] -= 1;
                    deg_buckets[degrees[nb]].insert(nb);
                    rem_edges += 1;
                    if degrees[nb] <= deg {
                        ignore.insert(nb);
                        to_check.insert(nb);
                        deg_buckets[degrees[nb]].remove(&nb);
                    }
                }
            }
        }
        // CAUTION: 'degrees' was not updated for removed vertices -> don't use it!

        // now, either no vertices or at least three are left!
        if g.size() == ignore.len() {  // no vertices left!
            return lb;
        } else {
            mc = min_cut_subgraph(g, &ignore, g.edge_count() - rem_edges);
            split = splittance_from_buckets(&deg_buckets);
        }
        lb = max(lb, min(split, mc));  // update the lower bound if necessary
        if mc >= split {  // as splittance will only decrease, no better LB can be found from now on
            break;
        }
        if lb >= upper_bound {
            break;  // upper bound reached already -> skip further improvement of lb!
        }
        if mc < deg {
            if mc == 0 {  // Graph is disconnected!
                let new_vert = (0..g.size()).collect::<AHashSet<usize>>().difference(&*ignore).map(|x| *x).collect::<AHashSet<usize>>();
                let (_, sg) = g.induced_subgraph_set(&new_vert);
                for comp in sg.connected_components_sets(None) {
                    if comp.len() > 3 {
                        mc += min_cut_lower_bound2(&sg, None, Some(&comp), Some(upper_bound-lb));
                    }
                }
                lb = max(lb, min(split, mc));  // update the lower bound if necessary
            }
            //println!("mc < deg -> break loop?"); // there could be a bridge in the graph. in this case the above procedure does not work
            break;
        }
    }

    // if lb improvement could still be useful (-> upper bound is not yet reached) and at least four
    // (but not all) vertices were removed, then add lb of subgraph induced by ignored vertices
    if lb < upper_bound && ignore.len() < g.size() && ignore.len() > 3 {
        let (_, mut ignored_sg) = g.induced_subgraph_set(&ignore);
        if g.edge_count() >= 4 {
            let packing_lb = compute_packing(&mut ignored_sg, Some(upper_bound-lb));
            if packing_lb > 0 {  // ignored_sg contains a FIS, i.e., it is no split cluster graph
                let mut mc_lb_sgs = 0;
                if packing_lb < upper_bound - lb {  // total lb does not yet exceed upper bound
                    for comp in ignored_sg.connected_components_sets(None) {
                        if comp.len() > 3 {
                            mc_lb_sgs += min_cut_lower_bound2(&ignored_sg, None, Some(&comp), Some(upper_bound-lb));
                        }
                    }
                }
                lb += max(packing_lb, mc_lb_sgs);
            }
        }
    }

    lb
}

/// globally computes MinCut Lower Bound -> loops over all connected components and returns sum of
/// MinCut Lower Bound of every component
/// if splittance 's' is given, it can be used for connected graphs
/// if upper bound 'ub' is given, the algorithm can possibly break early when lb >= ub already
pub fn min_cut_lower_bound2_global(g: &impl Graph, s: Option<usize>, ub: Option<usize>) -> usize {
    let comps = g.connected_components_sets(None);
    if comps.len() == 1 {   // 'g' is connected
        return min_cut_lower_bound2(g, s, None, ub);
    }
    let mut upper_bound = ub;
    let mut sum_lbs = 0;
    for comp in comps {
        if comp.len() <= 3 {  // split graph -> splittance (UB) = 0 -> LB = 0
            continue;
        }
        sum_lbs += min_cut_lower_bound2(g, None, Some(&comp), upper_bound);
        if let Some(u) = ub {
            if sum_lbs >= u {
                break;  // upper bound reached already -> break early!
            }
            upper_bound = Some(u - sum_lbs);
        }
    }
    sum_lbs
}

#[cfg(test)]
mod tests {
    use crate::graph_parser::{read_pgraph_map_from_file};
    use crate::graph_trait::Graph;
    use super::*;

    #[allow(unused)]
    fn setup() -> Vec<impl Graph> {
        let g03 = read_pgraph_map_from_file(&"../test_instances/g03.dimacs").unwrap();
        let g04 = read_pgraph_map_from_file(&"../test_instances/g04.dimacs").unwrap();
        let g05 = read_pgraph_map_from_file(&"../test_instances/g05.dimacs").unwrap();
        let g01 = read_pgraph_map_from_file(&"../test_instances/g01.dimacs").unwrap();
        let g06 = read_pgraph_map_from_file(&"../test_instances/g06.dimacs").unwrap();
        let g02 = read_pgraph_map_from_file(&"../test_instances/g02.dimacs").unwrap();
        let g07 = read_pgraph_map_from_file(&"../test_instances/g07.dimacs").unwrap();
        vec![g03, g04, g05, g01, g06, g02, g07]
    }

    #[test]
    fn min_cut_test() {
        let mut graphs = setup();
        assert_eq!(min_cut(&mut graphs[0]), 1);
        assert_eq!(min_cut(&mut graphs[3]), 0);  // disconnected graph
        assert_eq!(min_cut_lower_bound(&mut graphs[4], Some(2)), 2);
        assert_eq!(min_cut_lower_bound2(&mut graphs[4], Some(2), None, None), 2);
        assert_eq!(min_cut_lower_bound(&mut graphs[5], None), 1);  // P_10 -> Tree!
        assert_eq!(min_cut_lower_bound2(&mut graphs[5], None, None, None), 1);
        assert_eq!(min_cut_lower_bound(&mut graphs[6], None), 3);
        assert_eq!(min_cut_lower_bound2(&mut graphs[6], None, None, None), 3);
    }

    #[test]
    fn min_cut_global() {  // -> test disconnected graph
        let graphs = setup();
        assert_eq!(min_cut_lower_bound2_global(&graphs[3], None, None), 1);
    }
}