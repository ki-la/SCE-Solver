use std::collections::{HashMap, HashSet};
use crate::forb_ind_subgraph::splittance_from_buckets;
use crate::graph_trait::{Graph, VertexType};

/// computes and returns vertices that have to belong to the core or periphery
/// Periphery Rule: Label degree-one vertices as periphery and set their non-edges to forbidden
/// Core Rule: Label vertex with at least as many (or more than k) deg-1 nbs as other nbs as core
/// and set its edges to deg-1 nbs to permanent
/// CAUTION: assumes that 'g' is no split graph (might label both vertices in a K_2 as periphery)
pub fn deg_1_core_rules(g: &mut impl Graph, k: Option<usize>) -> (HashSet<usize>, HashSet<usize>) {
    debug_assert!(g.size() > 3);  // don't apply DR to trivial instances!
    let (buckets, _, _) = g.degree_sequence_buckets(None);
    debug_assert_ne!(splittance_from_buckets(&buckets), 0);
    let k_curr = match k {
        Some(x) => x,
        None => buckets.len(),  // no. of deg-1 neighbors is always smaller than this
    };
    let mut core: HashSet<usize> = HashSet::new();
    let mut periphery: HashSet<usize> = HashSet::new();

    if buckets[1].is_empty() {  // no deg-1 vertex!
        return (core, periphery);
    }

    // maps vertices with deg-1 neighbors to a set containing all their deg-1 neighbors
    let mut deg_1_nbs_map: HashMap<usize, HashSet<usize>> = HashMap::with_capacity(g.size());

    for v in &buckets[1] {  // buckets[1] contains all deg-1 (-> periphery) vertices
        periphery.insert(*v);
        g.label_vertex(*v, VertexType::Periphery);
        // v is a deg-1 vertex --> 'label' as periphery, set its non-edges to forbidden
        for nb_v in 0..g.size() {
            if nb_v == *v {
                continue;
            }
            if !g.has_edge(*v, nb_v) {    // (cheap even for AL based graphs as v has deg. 1)
                if !g.edited(*v, nb_v) {  // don't edit an edge between periphery vertices twice!
                    g.edit_edge(*v, nb_v); // non-edge of degree-one vertex -> forbidden
                }
            } else {
                let nbs = deg_1_nbs_map.entry(nb_v).or_insert(HashSet::from([*v]));
                nbs.insert(*v);
            }
        }
    }

    for (v, deg_1_nbs) in deg_1_nbs_map {
        if (deg_1_nbs.len() as f64 >= g.degree(v) as f64 / 2.0) ||
            (deg_1_nbs.len() > k_curr) {
            // v has to be a core vertex (improvement idea: if v has >k+1 independent (not
            // necessarily deg-1) neighbors, then it can also be labeled as core)
            core.insert(v);
            g.label_vertex(v, VertexType::Core);
            for deg_1_nb in deg_1_nbs {
                g.edit_edge(v, deg_1_nb);  // edge from degree-one to core vertex -> permanent
            }
        }
    }

    if !core.is_empty() {
        // check if there is a vertex that only has core neighbors -> label as periphery vertex
        for v in 0..g.size() {
            if only_core_neighbors(g, v) {
                g.label_vertex(v, VertexType::Periphery);
                periphery.insert(v);
            }
        }
    }


    (core, periphery)
}



/// returns whether vertex 'v' only has neighbors in graph 'g' that are labeled as Core
/// -> 'v' can be labeled as Periphery then
fn only_core_neighbors(g: &impl Graph, v: usize) -> bool {
    if !g.labeled(v) {  // only check yet unlabeled vertices
        if g.degree(v) <= g.labeled_core() {
            let mut only_core_nbs = true;
            for nb_v in g.neighbors_iter(v) {
                if *g.get_label(nb_v) != VertexType::Core {
                    only_core_nbs = false;
                    break;
                }
            }
            if only_core_nbs {  // vertex v only has core neighbors -> label as periphery!
                return true;
            }
        }
    }
    false
}

/// returns whether 'mod_edges' can be flipped (simultaneously) or if data reduction allows us to
/// derive that they are uneditable (see comments)
pub fn uneditable(g: &impl Graph, mod_edges: &[(usize, usize)], k: usize, dr: bool) -> bool {
    for (u, v) in mod_edges {
        if g.edited(*u, *v) {
            return true; // don't edit an edge twice!
        }
        if (*g.get_label(*u) == VertexType::Periphery ||
            *g.get_label(*v) == VertexType::Periphery) && !g.has_edge(*u, *v) {
            return true; // don't insert edges incident to periphery vertices
        }
        /*if g.has_edge(*u, *v) && ((*g.get_label(*u) == VertexType::Core &&
            *g.get_label(*v) == VertexType::Periphery && g.degree(*v) == 1) ||
            (*g.get_label(*v) == VertexType::Core && *g.get_label(*u) == VertexType::Periphery
            && g.degree(*u) == 1)) {
            return true; // edge between core vertex and degree-one periphery vertex is permanent -> don't remove!
        }*/
        if dr {
            if *g.get_label(*u) == VertexType::Core && *g.get_label(*v) == VertexType::Core &&
                g.has_edge(*u, *v) { // del edge between core vertices -> they cannot have common nbs
                let mut u_v= [*u, *v];
                u_v.sort_unstable_by_key(|x| g.degree(*x));
                for nb in g.neighbors_iter(u_v[0]) {
                    let mut cmn_nbs = 0;
                    if g.has_edge(nb, u_v[1]) {
                        if g.edited(*u, nb) && g.edited(*v, nb) {
                            // cannot remove edge uv -> they have a common permanent neighbor
                            return true;
                        }
                        cmn_nbs += 1;
                    }
                    if cmn_nbs > k {  // u and v have more than k common nbs -> cannot remove uv
                        return true;
                    }
                }
            }
        }
    }
    false
}

/// applies data reduction after mod_edge has been edited
/// returns vector of labeled vertices (for undoing the DR again)
pub fn dr_during_branching(g: &mut impl Graph, mod_edge: &(usize, usize), dr: bool) -> Vec<usize> {
    let mut labeled_vertices = Vec::new();
    let vertices = [mod_edge.0, mod_edge.1];
    for v in vertices {
        if g.labeled(v) {
            continue;  // vertex is labeled already
        }
        if g.degree(v) == 0 {
            g.label_vertex(v, VertexType::Periphery);
            labeled_vertices.push(v);
        } else if g.degree(v) == 1 {
            g.label_vertex(v, VertexType::Periphery);  // new periphery vertex!
            labeled_vertices.push(v);
            let mut new_core_vert = false;
            for nb_v in g.neighbors_iter(v) {  // O(1) as deg(v) == 1
                if g.edited(v, nb_v) {  // deg-1 vertex v has permanent edge -> nb has to be core
                    if !g.labeled(nb_v) {
                        labeled_vertices.push(nb_v);
                        new_core_vert = true;
                    } else if *g.get_label(nb_v) == VertexType::Periphery {
                        if g.degree(nb_v) == 1 {
                            new_core_vert = true;  // v has a deg-1 periphery nbs -> relabel v to Core
                        }
                    }
                }
            }
            if new_core_vert {
                g.label_vertex(labeled_vertices[labeled_vertices.len()-1], VertexType::Core);
            }
        }
        if g.has_edge(mod_edge.0, mod_edge.1) {  // edge insertion -> endpoints are core!
            g.label_vertex(v, VertexType::Core);
            labeled_vertices.push(v);
            if dr {
                // new core vertex -> check if any of its neighbors now only has core neighbors
                let labeled_vertices_len = labeled_vertices.len();
                for nb_v in g.neighbors_iter(v) {
                    if only_core_neighbors(g, nb_v) {
                        labeled_vertices.push(nb_v);  // nb_v only has core nbs -> periphery
                    } else if !g.labeled(nb_v) && nb_v != mod_edge.0 && nb_v != mod_edge.1 &&
                        g.edited(mod_edge.0, nb_v) && g.edited(mod_edge.1, nb_v) {
                        if (mod_edge.0 != v && !g.has_edge(mod_edge.0, nb_v)) ||
                            (mod_edge.1 != v && !g.has_edge(mod_edge.1, nb_v)) {
                            // P3 rule: u - v - nb_v is permanent P3 -> nb_v has to be periphery
                            labeled_vertices.push(nb_v);
                        }
                    }
                }
                // label periphery vertices
                if labeled_vertices.len() > labeled_vertices_len {
                    for u in labeled_vertices.iter().skip(labeled_vertices_len) {
                        g.label_vertex(*u, VertexType::Periphery);
                    }
                }
            }
        } else {  // edge deletion (uv removed -> now forbidden non-edge)
            if dr {
                let labeled_vertices_len = labeled_vertices.len();
                for nb_v in g.neighbors_iter(v) {
                    if !g.labeled(nb_v) &&
                        g.edited(nb_v, mod_edge.0) && g.edited(nb_v, mod_edge.1) &&
                        g.has_edge(nb_v, mod_edge.0) && g.has_edge(nb_v, mod_edge.1) {
                        // P3 rule: u -- nb_v -- v is permanent P3 -> nb_v has to be a core vertex
                        labeled_vertices.push(nb_v);
                    }
                }
                if labeled_vertices.len() > labeled_vertices_len {
                    // u and v have a permanent common neighbor -> have to end up in same cluster
                    // both have to be periphery vertices, as deleting an edge between core and
                    // periphery within one cluster is never optimal
                    if !g.labeled(mod_edge.0) {
                        labeled_vertices.push(mod_edge.0)
                    }
                    if !g.labeled(mod_edge.1) {
                        labeled_vertices.push(mod_edge.1)
                    }
                    for u in labeled_vertices.iter().skip(labeled_vertices_len) {
                        if *u == mod_edge.0 || *u == mod_edge.1 {  // u and v -> periphery
                            g.label_vertex(*u, VertexType::Periphery);
                        } else {
                            g.label_vertex(*u, VertexType::Core);  // nb_v -> core
                        }
                    }
                }
            }
        }
    }

    labeled_vertices
}

/// undo data reduction during branching by simply 'unlabeling' the vertices again
pub fn undo_dr_during_branching(g: &mut impl Graph, stack: Vec<usize>) {
    for v in stack {
        g.label_vertex(v, VertexType::Unknown);
    }
}

#[cfg(test)]
mod tests {
    use crate::data_reduction::deg_1_core_rules;
    use crate::graph_parser::{read_pgraph_map_from_file};
    use crate::graph_trait::Graph;
    use crate::graph_trait::VertexType;

    fn setup() -> Vec<impl Graph> {
        let g03 = read_pgraph_map_from_file(&"../test_instances/g03.dimacs").unwrap();
        let g02 = read_pgraph_map_from_file(& "../test_instances/g02.dimacs").unwrap();
        let g06 = read_pgraph_map_from_file(& "../test_instances/g06.dimacs").unwrap();
        let g07 = read_pgraph_map_from_file(& "../test_instances/g07.dimacs").unwrap();
        let g05 = read_pgraph_map_from_file(& "../test_instances/g05.dimacs").unwrap();
        vec![g03, g02, g06, g07, g05]
    }

    #[test]
    fn test_deg_1_core_rule() {
        let mut graphs = setup();
        let (c, p) = deg_1_core_rules(&mut graphs[0], Some(0));
        assert!(c.len() == 1 && p.len() == 1);
        assert_eq!(*graphs[0].get_label(1), VertexType::Periphery);
        assert_eq!(*graphs[0].get_label(0), VertexType::Core);
        assert!(graphs[0].edited(0, 1));
        let (c, p) = deg_1_core_rules(&mut graphs[1], None);
        assert!(c.len() == 2 && p.len() == 2);
        let (c, p) = deg_1_core_rules(&mut graphs[2], None);
        assert!(c.len() == 0 && p.len() == 0);
        let (c, p) = deg_1_core_rules(&mut graphs[3], None);
        assert!(c.len() == 0 && p.len() == 0);
        let (_, p) = deg_1_core_rules(&mut graphs[4], None);
        for p_v in p {
            assert_eq!(graphs[4].degree(p_v), 1);
        }
    }
}