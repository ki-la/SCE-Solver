use std::collections::{BTreeMap, HashMap};
use ahash::AHashSet;
use indexmap::{IndexMap, IndexSet};  // -> used for determinism
use crate::forb_ind_subgraph::FisType::{NONE, _2K2, C4, C5, P5, BOWTIE, NECKTIE};
use crate::graph_trait::Graph;

#[derive(PartialEq, Debug)]
pub enum FisType { NONE, _2K2, C4, C5, P5, BOWTIE, NECKTIE }

#[derive(Debug)]
pub struct FIS {
    pub sg_type: FisType,
    pub v: Vec<usize>,
}

/// calculates splittances of the graph's connected components using (HS81) and returns their sum
pub fn splittance_global(g: &impl Graph) -> usize {
    let comps = g.connected_components(None);

    if comps.len() == 1 {   // graph is connected
        return splittance_sg(g, None);
    }

    let mut s_glob = 0;
    for comp in comps {
        s_glob += splittance_sg(g, Some(&comp));
    }

    s_glob
}

/// calculates and returns splittance of the (sub-)graph 'g' as described in (HS81)
/// only takes vertices in 'vert' into consideration (all vertices if vert.is_none())
/// CAUTION: assumes that 'g' is connected!
/// (if it's not, then the cost to make it ONE split graph is returned)
pub fn splittance_sg(g: &impl Graph, vert: Option<&[usize]>) -> usize {

    if g.size() < 4 || (vert.is_some() && vert.unwrap().len() < 4) {
        return 0;
    }

    let (deg_seq , _) = g.degree_sequence(vert, false, false); // non-increasing

    splittance(&deg_seq)
}

/// calculates and returns splittance of a graph w/ degree sequence 'deg_seq' as described in (HS81)
/// in O(n) time. CAUTION: assumes that 'deg_seq' is non-increasing!
// Source (HS81): Peter L. Hammer and Bruno Simeone. The splittance of a graph. [Combinatorica '81]
pub fn splittance(deg_seq: &[usize]) -> usize {
    let mut m:usize = 0;                                    // (2) in (HS81)
    let mut sum_c:usize = 0;                                // sum(core vertices' degrees)
    let mut sum_p:usize = 0;                                // sum(periphery vertices' degrees)
    for (i, deg_i) in deg_seq.iter().enumerate() {
        if *deg_i < i {
            break;
        }
        m = i;
        sum_c += *deg_i;
    }
    m += 1; // undo index shift ((HS81) starts with 1, we with 0)

    for deg_j in deg_seq.iter().skip(m) {
        sum_p += *deg_j;
    }

    let s = m * (m - 1) + sum_p - sum_c ; // formula (3) in (HS81) for k=m
    assert_eq!(s % 2, 0);
    s / 2
}

/// calculates and returns splittance of a graph using degree_buckets and formula from (HS81) in
/// O(n) time (iterates buckets once)
pub fn splittance_from_buckets(degree_buckets: &[IndexSet<usize>]) -> usize {
    let mut m:usize = 0;                                    // (2) in (HS81)
    let mut sum_c:usize = 0;                                // sum(core vertices' degrees)
    let mut sum_p:usize = 0;                                // sum(periphery vertices' degrees)
    for (deg, bucket) in degree_buckets.iter().enumerate().rev() {
        if bucket.is_empty() { continue; }
        if deg < m {  // all vertices in this bucket belong to the periphery
            sum_p += deg * bucket.len();
        } else if deg < m + (bucket.len() - 1) {  // deg-m+1 vertices in this bucket belong to the core, the rest to the periphery
            sum_c += deg * (deg - m + 1);
            sum_p += deg * (bucket.len() - (deg - m + 1));
            m += deg - m + 1;
            continue;
        } else {  // all vertices in this bucket belong to the core
            m += bucket.len();
            sum_c += deg * bucket.len();
        }
    }

    let s = m * (m - 1) + sum_p - sum_c ; // formula (3) in (HS81) for k=m
    assert_eq!(s % 2, 0);
    s / 2
}

/// finds and returns one forbidden induced subgraph
pub fn find_fis_global(g: &impl Graph) -> FIS {
    let comps = g.connected_components(None);

    if comps.len() == 1 {   // graph is connected
        return find_fis(g, None);
    }

    for comp in comps {
        if comp.len() < 4 {
            continue;  // every FIS has at least four vertices
        }
        let fis = find_fis(g, Some(&comp));
        if fis.sg_type != NONE {
            return fis;
        }
    }

    FIS { sg_type: NONE, v: Vec::new() }  // no FIS in any connected component
}

/// finds and returns a forbidden induced subgraph for a split cluster graph (C4, C5, P5,
/// BOWTIE, NECKTIE) or NONE if the graph is split following the algorithm described in (BHK15)
/// CAUTION: assumes that the graph is connected, only considers 'vert'
// Source (BHK15): Sharon Bruckner, Falk Hüffner, and Christian Komusiewicz. A graph
// modification approach for finding core-periphery structures in protein interaction networks.
// [Algorithms Mol Biol '15]
fn find_fis(g: &impl Graph, vert: Option<&[usize]>) -> FIS {
    let fis: FIS = certifying_split(g, vert);
    if fis.sg_type == NONE {
        return fis;
    }

    if fis.sg_type == C4 || fis.sg_type == C5 {
        return fis;
    }
    debug_assert_eq!(fis.sg_type, _2K2);
    debug_assert_eq!(fis.v.len(), 4);

    // extend 2K2 to P5, bowtie, or necktie
    // find shortest 'path' between any vertices in the two K2s
    let paths = vec![
        g.shortest_path(fis.v[0], fis.v[2], vert).unwrap(),
        g.shortest_path(fis.v[0], fis.v[3], vert).unwrap(),
        g.shortest_path(fis.v[1], fis.v[2], vert).unwrap(),
        g.shortest_path(fis.v[1], fis.v[3], vert).unwrap()
    ];
    let path = paths.iter().min_by_key(|x| x.len()).unwrap();
    if path.len() > 4 {
        // there cannot be edges between the vertices of 'path' ('path' is shortest) --> P5!
        return FIS { sg_type: P5, v: vec![path[0], path[1], path[2], path[3], path[4]] };
    }
    let mut y1 = fis.v[2];
    let mut y2 = fis.v[3];
    if path.len() == 4 {
        // find new 2K2 s.t. there is a vertex that is adjacent to both
        y1 = path[2];
        y2 = path[3];
    }
    assert_eq!(path.len(), 3);
    let x = if fis.v[0] != path[0] { fis.v[0] } else { fis.v[1] };
    let y = if y1 != path[2] { y1 } else { y2 };
    // x -- path[0] -- path[1] -- path[2] -- y
    if g.has_edge(x, path[1]) {
        if g.has_edge(y, path[1]) {
            FIS { sg_type: BOWTIE, v: vec![x, path[0], path[1], path[2], y] }
        } else {
            FIS { sg_type: NECKTIE, v: vec![x, path[0], path[1], path[2], y] }
        }
    } else if g.has_edge(y, path[1]) {
            FIS { sg_type: NECKTIE, v: vec![y, path[2], path[1], path[0], x] }
    } else {
        FIS { sg_type: P5, v: vec![x, path[0], path[1], path[2], y] }
    }
}

/// finds and returns a forbidden induced subgraph for a split graph (_2K2, C4, or C5) -> cf. (HK07)
/// CAUTION: assumes that 'g' is connected, only considers 'vert' (all vertices, if vert.is_none())
// Source (HK07): Pinar Heggernes and Dieter Kratsch. Linear-time certifying recognition
// algorithms and forbidden induced subgraphs. [Nord. J. Comput '07]
// -> also compared with/ inspired by code published along with (BHK15)
fn certifying_split(g: &impl Graph, vert: Option<&[usize]>) -> FIS {
    let (deg_seq, mut deg_ord) = g.degree_sequence(vert, true, false);
    if splittance(&deg_seq) == 0 {  // 'g' is a split graph already -> no FIS can be found
        return FIS { sg_type: NONE, v: vec![] };
    }
    deg_ord.reverse();
    let (is_peo, conflict_triple, max_clique) = is_peo(g, &deg_ord, vert);
    if !is_peo {
        // deg_ord is not a peo -> g cannot be a split graph!
        if let Some((i, j, k)) = conflict_triple {
            assert!(g.has_edge(i,j) && g.has_edge(i, k) && !g.has_edge(j, k));
            let mut x: Option<usize> = None;
            let mut y: Option<usize> = None;
            for z in g.neighbors(j, vert) {
                if z != i && !g.has_edge(i, z) {
                    if g.has_edge(k, z) {
                        return FIS { sg_type: C4, v: vec![i, j, z, k] } ;           // C4 found
                    }
                    x = Some(z);  // x has edge to j but not to i or k
                }
            }
            assert!(x.is_some());
            for z in g.neighbors(k, vert) {
                if z != i && !g.has_edge(z, i) && !g.has_edge(z, j) {
                    y = Some(z);  // y has edge to k but not to i or j
                }
            }
            assert!(y.is_some());
            return if g.has_edge(x.unwrap(), y.unwrap()) {
                FIS { sg_type: C5, v: vec![i, j, x.unwrap(), y.unwrap(), k] }       // C5 found
            } else {
                FIS { sg_type: _2K2, v: vec![j, x.unwrap(), k, y.unwrap()] }        // 2K2 found
            }
        }
        panic!("if deg_ord is no peo, then conflict_triple cannot be None!");
    } else {
        // deg_ord is a peo -> g is chordal (no C_4, C_5, ...)
        let mut idx: isize = (deg_ord.len() - 1) as isize;
        // calculate core 'c' or return a 2K2
        let mut c: IndexSet<usize> = IndexSet::new();
        while c.len() < max_clique {
            let i = idx as usize;
            for y in &c {
                let y = *y;
                if !g.has_edge(y, deg_ord[i]) {
                    // deg_ord[i] has a non-neighbor in c
                    for x in g.neighbors(deg_ord[i], vert) {
                        if !c.contains(&x) {
                            for z in g.neighbors(y, vert) {
                                if !g.has_edge(z, deg_ord[i]) && !g.has_edge(z, x) {
                                    return FIS { sg_type: _2K2, v: vec![x, deg_ord[i], y, z] };
                                }
                            }
                        }
                    }
                    // (HK07) have proven that deg_ord[i] has to have a neighbor outside c
                    // if it has a non-neighbor in c
                    panic!("this code should be unreachable!");
                }
            }
            c.insert(deg_ord[i]);
            idx -= 1;
        }
        // calculate periphery 'p' or return 2K2
        let mut p: IndexSet<usize> = IndexSet::new();
        while idx >= 0 {
            let i = idx as usize;
            for x in &p {
                let x = *x;
                if g.has_edge(x, deg_ord[i]) {
                    // deg_ord[i] has a neighbor in p
                    for y in c {
                        if !g.has_edge(x, y) && !g.has_edge(deg_ord[i], y) {
                            for z in g.neighbors(y, vert) {
                                if !g.has_edge(x, z) {
                                    return FIS { sg_type: _2K2, v: vec![x, deg_ord[i], y, z] };
                                }
                            }
                        }
                    }
                    // (HK07) have proven that if deg_ord[i] has a neighbor in p, then there are
                    // two (adjacent) vertices in c that are non-neighbors to both of these
                    // vertices
                    panic!("this code should be unreachable!");
                }
            }
            p.insert(deg_ord[i]);
            idx -= 1;
        }
    }
    panic!("this code should be unreachable!");
}

/// tests if 'peo' is a perfect elimination scheme for 'g', returns a conflict triple otherwise
/// follows Algorithm 4.2 in (Gol04)
/// CAUTION: assumes that 'g' is connected, only considers 'vert' (all vertices, if vert.is_none())
/// might panic if vertices in 'peo' are not identical with those in 'vert'
// Source (Gol04): Martin Charles Golumbic. Chapter 4 - Triangulated graphs. [Algorithmic Graph
// Theory and Perfect Graphs '04]
// -> also compared with/ inspired by code published along with (BHK15)
fn is_peo(g: &impl Graph, peo: &[usize], vert: Option<&[usize]>) -> (bool, Option<(usize, usize, usize)>, usize) {
    let mut max_clique: usize = 1;

    // calculate the inverse of 'peo'
    let mut peo_inv: HashMap<usize,usize> = HashMap::with_capacity(peo.len());
    for (i, p_i) in peo.iter().enumerate() {
        peo_inv.insert(*p_i, i);
    }
    for v in peo.iter().take(peo.len()-1) {
        let v: usize = *v;
        let mut u: Option<usize> = None;
        let mut peo_inv_u = peo.len();
        let mut nbs: usize = 0;
        for x in g.neighbors(v, vert) {
            if peo_inv[&v] < peo_inv[&x] {
                // x appears later in the peo than v
                nbs += 1;
                if peo_inv[&x] < peo_inv_u {
                    u = Some(x);
                    peo_inv_u = peo_inv[&x];
                }
            }
        }
        if nbs + 1 > max_clique {
            max_clique = nbs + 1;
        }
        if let Some(u_new) = u {
            for w in g.neighbors(v, vert) {
                if peo_inv[&v] < peo_inv[&w] && w != u_new {
                    // u_new and w are neighbors of v and appear later than v in 'peo'
                    if !g.has_edge(u_new, w) {
                        // v's neighbors in peo[i..n] do not form a clique -> no peo!
                        // println!("found conflict triple: {v}, {u_new}, {w}");
                        return (false, Some((v, u_new, w)), max_clique);
                    }
                }
            }
        } else {
            continue;  // v has no neighbors in peo[i+1..n]
        }
    }
    (true, None, max_clique)  // 'peo' is a partial elimination ordering
}

/// checks and returns if 'fis' is a C4, C5, P5, BOWTIE, or NECKTIE in graph 'g'
pub fn is_valid_fis(g: &impl Graph, fis: &FIS) -> bool {
    if fis.sg_type == NONE {
        return false;
    }

    let (edges, non_edges) = match fis.sg_type {
        C4      =>  (vec![(0, 1), (1, 2), (2, 3), (3, 0)],
                     vec![(0, 2), (1, 3)]),
        C5      =>  (vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
                     vec![(0, 2), (0, 3), (1, 3), (1, 4), (2, 4)]),
        P5      =>  (vec![(0, 1), (1, 2), (2, 3), (3, 4)],
                     vec![(0, 2), (0, 3), (1, 3), (1, 4), (2, 4), (4, 0)]),
        BOWTIE  =>  (vec![(0, 1), (1, 2), (2, 3), (3, 4), (0, 2), (2, 4)],
                     vec![(0, 3), (1, 3), (1, 4), (4, 0)]),
        NECKTIE =>  (vec![(0, 1), (1, 2), (2, 3), (3, 4), (0, 2)],
                     vec![(0, 3), (1, 3), (1, 4), (2, 4), (4, 0)]),
        _       =>  return false,
    };
    for (u, v) in edges {
        if !g.has_edge(fis.v[u], fis.v[v]) {
            return false;
        }
    }
    for (u, v) in non_edges {
        if g.has_edge(fis.v[u], fis.v[v]) {
            return false;
        }
    }
    true
}

/// if a packing for 'g' was already computed, it is updated (else, compute_packing is returned)
/// stops early if packing lower bound already exceeds the given upper bound
pub fn update_packing(g: &mut impl Graph, ub: usize) -> usize {
    let mut fis = g.remove_fis(0);
    if fis.is_none() {
        return compute_packing(g, Some(ub));  // no FIS in g's packing -> compute from scratch!
    }
    // at least one FIS in g's packing -> update!
    let mut packing: Vec<FIS> = Vec::new();
    let mut removed_vert: AHashSet<usize> = AHashSet::with_capacity(g.size());
    let mut lb: usize = 0;
    while fis.is_some() {
        let fis_tmp = fis.unwrap();
        let mut only_act_vert = true;
        for v in &fis_tmp.v {
            if removed_vert.contains(v) {
                only_act_vert = false;
                break;
            }
        }
        if only_act_vert && is_valid_fis(g, &fis_tmp) {
            lb += 1;
            if fis_tmp.sg_type == C5 {
                lb += 1;
            }
            if lb > ub {
                break;
            }
            if !remove_fis_vert_from_set(g, &fis_tmp, &mut removed_vert) {
                // 'fis_tmp' cannot be destroyed due to edited edges -> no sol. in current branch
                return ub + 1;
            }
            packing.push(fis_tmp);
        }
        fis = g.remove_fis(0);
    }
    if lb > ub {
        return lb;  // lb already exceeds ub -> stop early and return lb
    }
    if g.size() - removed_vert.len() >= 4 {
        // at least four vertices left -> try to extend packing...
        let active_vert: Vec<usize> = (0..g.size()).filter(|x| !removed_vert.contains(x)).collect();
        let mut comps = g.connected_components(Some(&active_vert));
        while lb <= ub && !comps.is_empty() {
            let comp = comps.pop().unwrap();
            if comp.len() < 4 {
                continue;                   // every FIS has at least four vertices
            }
            // find and 'remove' FISs in 'comp'
            let fis = find_fis(g, Some(&comp));
            if fis.sg_type == NONE {        // no FIS
                continue;
            }                               // FIS found
            lb += 1;
            if fis.sg_type == C5 {
                lb += 1;                    // count C5s twice
            }
            comps.append(&mut remove_fis_vertices(g, &comp, &fis));
            packing.push(fis);
        }
    }

    g.store_packing(packing);
    lb
}

/// computes and returns a FIS packing (and its size w/ C5s counted twice) of the input graph
/// if an upper bound is given, it might stop early once lb > ub
pub fn compute_packing(g: &mut impl Graph, ub: Option<usize>) -> usize {
    let mut packing: Vec<FIS> = Vec::new();
    let mut lb: usize = 0;
    let upper_b = match ub {
        Some(u) => u,
        None    => g.size(),  // as every FIS has 4-5 vertices, the packing will never have |V| FISs
    };
    // compute connected components and find FISs in them independently
    let mut comps = g.connected_components(None);
    while !comps.is_empty() {
        let comp = comps.pop().unwrap();
        if comp.len() < 4 {
            continue;                   // every FIS has at least four vertices
        }
        // find and 'remove' FISs in 'comp'
        let fis = find_fis(g, Some(&comp));
        if fis.sg_type == NONE {        // no FIS
            continue;
        }                               // FIS found
        lb += 1;
        if fis.sg_type == C5 {
            lb += 1;                    // count C5s twice
        }
        if lb > upper_b {               // break early as lb already exceeds ub
            packing.push(fis);
            break;
        }
        comps.append(&mut remove_fis_vertices(g, &comp, &fis));
        packing.push(fis);
    }
    g.store_packing(packing);
    lb
}

/// removes as few vertices as possible from 'fis' and adds them to the set 'removed_vert'
/// returns false, if this is not possible because it is a 'permanent' FIS
fn remove_fis_vert_from_set(g: &impl Graph, fis: &FIS, removed_vert: &mut AHashSet<usize>) -> bool {
    let evp = edited_vertex_pairs(g, fis);
    let poss_edits = (fis.v.len() * (fis.v.len() - 1)) / 2;  // 6 or 10
    debug_assert!(evp.len() <= poss_edits);
    if evp.len() == poss_edits || (fis.sg_type == C5 && evp.len() == poss_edits - 1) {
        return false;  // no edits that could destroy the 'fis' possible anymore!
    }
    let mut max_deg_vertex = (fis.v[0], 0);
    let mut deg = IndexMap::new();
    for u in &fis.v {
        for v in g.neighbors_iter(*u) {
            if !fis.v.contains(&v) && !removed_vert.contains(&v) {
                *deg.entry(*u).or_insert(0) += 1;
            }
        }
        let d = deg.entry(*u).or_insert(0);
        if *d > max_deg_vertex.1 {
            max_deg_vertex = (*u, *d);
        }
    }
    let mut max_deg_pair: Option<((usize,usize), usize)> = None;
    if evp.len() > 0 {  // at least one edited edge pair!
        for (i, j) in evp.iter() {
            if max_deg_pair.is_none() || deg.get(i).unwrap() + deg.get(i).unwrap() > max_deg_pair.unwrap().1 {
                max_deg_pair = Some(((*i, *j), deg.get(i).unwrap() + deg.get(j).unwrap()));
            }
        }
    }
    if fis.sg_type == BOWTIE || fis.sg_type == NECKTIE || fis.sg_type == P5 {
        let active_vert: Vec<usize> = (0..g.size()).filter(|x| !removed_vert.contains(x)).collect();
        let keep = vertices_to_keep(g, fis, &active_vert);
        if max_deg_pair.is_none() || keep.1  > max_deg_pair.unwrap().1 {
            for fis_vertex in &fis.v {
                if !keep.0.contains(fis_vertex) {
                    removed_vert.insert(*fis_vertex);
                }
            }
            return true;
        }
    }
    if let Some(((u, v), _)) = max_deg_pair {
        for fis_vertex in &fis.v {
            if *fis_vertex != u && *fis_vertex != v {
                removed_vert.insert(*fis_vertex);
            }
        }
    } else {
        for fis_vertex in &fis.v {
            if *fis_vertex != max_deg_vertex.0 {
                removed_vert.insert(*fis_vertex);
            }
        }
    }

    true
}

/// removes most vertices in 'fis', so they cannot be included in another FIS
/// returns the new connected components of the graph
fn remove_fis_vertices(g: &mut impl Graph, comp: &[usize], fis: &FIS) -> Vec<Vec<usize>> {
    let mut new_active = Vec::with_capacity(comp.len());
    for c in comp {
        if !fis.v.contains(c) {
            new_active.push(*c);
        }
    }
    if fis.sg_type == BOWTIE || fis.sg_type == NECKTIE || fis.sg_type == P5 {
        new_active.append(&mut vertices_to_keep(g, fis, comp).0);
    } else {
        let degrees = degree_outside_fis(g, fis, comp);
        new_active.push(*degrees.keys().max_by_key(|x| degrees[x]).unwrap());
    }
    g.connected_components(Some(&new_active))
}

/// returns vertex pairs in 'fis' that are marked as edited
fn edited_vertex_pairs(g: &impl Graph, fis: &FIS) -> Vec<(usize, usize)> {
    let mut evp = Vec::new();
    for u in &fis.v {
        for v in &fis.v {
            if *u < *v && g.edited(*u, *v) {
                evp.push((*u, *v));
            }
        }
    }
    evp
}

/// returns which vertices in 'fis' should be kept (based on their degree outside the FIS)
/// -> also returns sum of degrees outside the FIS
fn vertices_to_keep(g: &impl Graph, fis: &FIS, vert: &[usize]) -> (Vec<usize>, usize) {
    // when calculating a packing by repeatedly searching for a FIS and deleting it, up to two
    // vertices of that FIS can be kept as editing their non-edge would create a new FIS
        // P5: 0-1-2-3-4             --> edges 02, 03, 04, 14, 24 create a new FIS when edited
        // nt: tri(0-1-2)-3-4        --> edges 02, 04, 12, 14, 24 create a new FIS when edited
        // bt: tri(0-1-2)tri(2-3-4)  --> edges 02, 12, 23, 24 create a new FIS when edited
    let keep = match fis.sg_type {
        P5      => vec![(0, 3), (0, 4), (1, 4), (2, 4)],
        NECKTIE => vec![(0, 4), (1, 2), (1, 4), (2, 4)],
        BOWTIE  => vec![(1, 2), (2, 3), (2, 4)],
        _       => return (Vec::new(), 0),
    };
    let deg = degree_outside_fis(g, fis, vert);
    let mut max_deg_sum = (0, 2, deg[&fis.v[0]] + deg[&fis.v[2]]);
    for (i, j) in keep {
        if deg[&fis.v[i]] + deg[&fis.v[j]] > max_deg_sum.2 {
            max_deg_sum = (i, j, deg[&fis.v[i]] + deg[&fis.v[j]]);
        }
    }
    // in case of a necktie, it can still be better to keep the single vertex 3
    if fis.sg_type == NECKTIE && deg[&fis.v[3]] > max_deg_sum.2 {
        return (vec![fis.v[3]], deg[&fis.v[3]])
    }
    (vec![fis.v[max_deg_sum.0], fis.v[max_deg_sum.1]], max_deg_sum.2)
}

/// returns a HashMap mapping the nodes in 'fis' to the number of their neighbors outside of 'fis'
pub fn degree_outside_fis(g: &impl Graph, fis: &FIS, vert: &[usize]) -> BTreeMap<usize, usize> {
    let mut degrees = BTreeMap::new();
    for u in &fis.v {
        for v in g.neighbors(*u, Some(vert)) {
            if !fis.v.contains(&v) {
                *degrees.entry(*u).or_insert(0) += 1;
            }
        }
        degrees.entry(*u).or_insert(0);
    }
    degrees
}

/// calculates and returns the edits necessary to turn 'g's connected components into split graphs
pub fn compute_solution_from_splittance_global(g: &impl Graph) -> Vec<(usize, usize)> {
    let comps = g.connected_components(None);
    if comps.len() == 1 {   // graph is connected
        return compute_solution_from_splittance(g, None);
    }
    let mut s_glob:Vec<(usize, usize)> = Vec::new();
    for comp in comps {
        s_glob.append(&mut compute_solution_from_splittance(g, Some(&comp)));
    }
    s_glob
}

/// calculates and returns the edits necessary to turn 'g' into a split graph (cf. (HS81))
/// if !vert.is_none(), then only the vertices in 'vert' are taken into account
/// CAUTION: assumes that 'g' is connected
pub fn compute_solution_from_splittance(g: &impl Graph, vert: Option<&[usize]>) -> Vec<(usize, usize)> {
    let (deg_seq, deg_ord) = g.degree_sequence(vert, true, false);
    let mut solution: Vec<(usize, usize)> = Vec::new();
    let mut c: AHashSet<usize> = AHashSet::new();  // core
    let mut p: AHashSet<usize> = AHashSet::new();  // periphery

    for (i, deg_i) in deg_seq.iter().enumerate() {
        if *deg_i >= i {
            c.insert(deg_ord[i]);
        } else {
            p.insert(deg_ord[i]);
        }
    }
    if let Some(vt) = vert {
        let mut no_pairs = 0;
        for (i, u) in vt.iter().enumerate() {
            for v in vt.iter().skip(i+1) {
                no_pairs += 1;
                if (c.contains(u) && c.contains(v) && !g.has_edge(*u, *v)) ||   // core vertices have to be adjacent
                    (p.contains(u) && p.contains(v) && g.has_edge(*u, *v)) {    // periphery vertices mustn't be adjacent
                    solution.push((*u, *v));
                }
            }
        }
        debug_assert_eq!(2 * no_pairs, vt.len() * (vt.len() - 1));
    } else {
        for u in 0..g.size() - 1 {
            for v in u+1..g.size() {
                if (c.contains(&u) && c.contains(&v) && !g.has_edge(u, v)) ||   // core vertices have to be adjacent
                    (p.contains(&u) && p.contains(&v) && g.has_edge(u, v)) {    // periphery vertices mustn't be adjacent
                    solution.push((u, v));
                }
            }
        }
    }
    solution
}


#[cfg(test)]
mod tests {
    use crate::forb_ind_subgraph::*;
    use crate::graph_parser::*;
    use crate::graph_trait::Graph;

    fn setup() -> Vec<impl Graph> {
        let files = ["../test_instances/g03.dimacs",
                            "../test_instances/g02.dimacs",
                            "../test_instances/g06.dimacs",
                            "../test_instances/g07.dimacs",
                            "../test_instances/g01.dimacs"];
        let mut graphs = Vec::new();
        for f in files {
            graphs.push(read_graph_from_file(f).unwrap());
        }
        graphs
    }

    #[test]
    fn test_splittance() {
        /* SPLITTANCE of special graphs:
            - P_n -> n-4 (for n >= 5)
            - C_n -> n-3
            - K_p,q -> v(v-1) / 2 (v = min{p,q}) */
        let graphs = setup();
        assert_eq!(splittance_global(&graphs[0]), 4);
        assert_eq!(splittance_global(&graphs[1]), 10-4);
        assert_eq!(splittance_global(&graphs[2]), 5-3);
        assert_eq!(splittance_global(&graphs[3]), 3);
        assert_eq!(splittance_global(&graphs[4]), 1);
        assert_eq!(compute_solution_from_splittance_global(&graphs[4]).len(), 1);
        assert_eq!(splittance_sg(&graphs[4], Some(&vec![4, 5, 6, 7, 8, 9])), 1);
        assert_eq!(splittance_sg(&graphs[4], None), 2);
        assert_eq!(compute_solution_from_splittance(&graphs[4], None).len(), 2);
    }

    #[test]
    fn test_splittance_from_buckets() {
        let graphs = setup();
        assert_eq!(splittance_from_buckets(&graphs[0].degree_sequence_buckets(None).0), 4);
        assert_eq!(splittance_from_buckets(&graphs[1].degree_sequence_buckets(None).0), 10-4);
        assert_eq!(splittance_from_buckets(&graphs[2].degree_sequence_buckets(None).0), 5-3);
        assert_eq!(splittance_from_buckets(&graphs[3].degree_sequence_buckets(None).0), 3);
        assert_eq!(splittance_from_buckets(&graphs[4].degree_sequence_buckets(None).0), 2);
    }

    #[test]
    fn test_is_peo() {
        let graphs = setup();
        assert!(!is_peo(&graphs[0], &vec![1, 9, 0, 8, 2, 3, 4, 5, 6, 7], None).0);
        let mut peo = is_peo(&graphs[1], &vec![0,1,2,3,4,5,6,7,8,9], None); // P_10
        assert!(peo.0);  // all vertices in 0..10 are simplical
        assert!(peo.1.is_none());
        assert_eq!(peo.2, 2);
        peo = is_peo(&graphs[1], &vec![0,9,2,1,3,6,7,8], Some(&vec![0,1,2,3,6,7,8,9]));
        assert!(!peo.0);  // vertex 2 is not simplical (nbs 1 and 3 don't form a clique)
        assert_eq!(peo.1.unwrap(), (2, 1, 3));
    }

    #[test]
    fn test_certifying_split() {
        let graphs = setup();
        assert_ne!(certifying_split(&graphs[0], None).sg_type, FisType::NONE);
        assert_eq!(certifying_split(&graphs[1], None).sg_type, FisType::_2K2);  // P_10
        assert_eq!(certifying_split(&graphs[2], None).sg_type, FisType::C5);    // C_5
        assert_eq!(certifying_split(&graphs[3], None).sg_type, FisType::C4);    // K_3_3
        assert_eq!(certifying_split(&graphs[4], None).sg_type, FisType::_2K2);
        assert_eq!(certifying_split(&graphs[4], Some(&vec![1,4,5,6,7,9])).sg_type, FisType::NONE);
    }

    #[test]
    fn test_find_fis() {
        let graphs = setup();
        assert_ne!(find_fis_global(&graphs[0]).sg_type, FisType::NONE);
        assert!(is_valid_fis(&graphs[0], &FIS{sg_type: FisType::BOWTIE, v: vec![0,9,7,4,5]}));
        assert!(!is_valid_fis(&graphs[0], &FIS{sg_type: FisType::BOWTIE, v: vec![0,9,4,5,7]}));
        assert!(!is_valid_fis(&graphs[0], &FIS{sg_type: FisType::C4, v: vec![0,9,7,5,6]}));
        assert!(!is_valid_fis(&graphs[0], &FIS{sg_type: FisType::NECKTIE, v: vec![0,5,6]}));
        assert_eq!(find_fis_global(&graphs[1]).sg_type, FisType::P5);    // P_10
        assert_eq!(find_fis_global(&graphs[2]).sg_type, FisType::C5);    // C_5
        let mut fis = find_fis(&graphs[2], None);
        assert!(is_valid_fis(&graphs[2], &fis));
        fis.v.sort();
        assert_eq!(fis.v, vec![0,1,2,3,4]);
        fis = find_fis(&graphs[3], Some(&vec![0,1,3,4]));            // K_3_3
        assert!(is_valid_fis(&graphs[3], &fis));
        fis.v.sort(); // caution: 'invalidates' FIS as order is changed
        assert_eq!(fis.v, vec![0,1,3,4]);
        assert!(!is_valid_fis(&graphs[3], &fis));
    }

    #[test]
    fn test_degree_outside_fis() {
        let graphs = setup();
        let fis = FIS { sg_type: FisType::BOWTIE , v: vec![0,9,7,4,5] };
        let dof = degree_outside_fis(&graphs[0], &fis, &vec![0,1,2,3,4,5,6,7,8,9]);
        assert!(dof[&0] == 2 && dof[&7] == 2 && dof[&9] == 0 && dof[&4] == 3 && dof[&5] == 3);
        let fis = FIS { sg_type: FisType::C5, v: vec![0,1,2,3,4] };
        let dof = degree_outside_fis(&graphs[2], &fis, &vec![0,1,2,3,4]);
        assert_eq!(dof.values().collect::<Vec<&usize>>(), vec![&0,&0,&0,&0,&0]);
    }

    #[test]
    fn test_remove_fis_vertices() {
        let mut graphs = setup();
        let nt = FIS {sg_type: NECKTIE, v: vec![3,5,6,0,1]};
        let c4 = FIS {sg_type: C4, v: vec![8,2,7,4]};
        let new_active = remove_fis_vertices(&mut graphs[0], &[0,1,2,3,4,5,6,7,8,9], &nt);
        assert!(new_active[0].contains(&6) && (new_active[0].contains(&3) || new_active[0].contains(&5)));
        assert!(vertices_to_keep(&mut graphs[0], &nt, &[0,1,2,3,4,5,6,7,8,9]).0.contains(&6));
        let new_active = remove_fis_vertices(&mut graphs[0], &[0,1,2,3,4,5,6,7,8,9], &c4);
        assert!(!new_active[0].contains(&8)); // 8 only has fewest nbs outside the C4 -> remove!
        assert_eq!(vertices_to_keep(&mut graphs[0], &c4, &[0,1,2,3,4,5,6,7,8,9]).0, vec![]);
        let p5 = FIS {sg_type: P5, v: vec![2,3,4,5,6]};
        let mut new_active = remove_fis_vertices(&mut graphs[1], &[0,1,2,3,4,5,6,7,8,9], &p5);
        assert_eq!(new_active.len(), 2);
        new_active.sort_by_key(|x| x.len());
        new_active[0].sort();
        new_active[1].sort();
        assert!(new_active[0] == vec![0,1,2] && new_active[1] == vec![6,7,8,9]);
        assert_eq!(vertices_to_keep(&mut graphs[1], &p5, &[0,1,2,3,4,5,6,7,8,9]).0, vec![2,6]);
    }

    #[test]
    fn test_compute_packing() {
        let mut graphs = setup();
        graphs[1].remove_edge(7, 8);
        assert_eq!(compute_packing(&mut graphs[1], None), 1);     // P_10
        assert_eq!(graphs[1].get_fis().unwrap().sg_type, FisType::P5);
        assert!(graphs[1].get_fis().is_none());
        assert_eq!(compute_packing(&mut graphs[2], None), 2);     // C_5 (needs two edits)
        assert_eq!(graphs[2].get_fis().unwrap().sg_type, FisType::C5);
        assert!(graphs[2].get_fis().is_none());
        assert_eq!(compute_packing(&mut graphs[3], None), 1);     // K_3_3
        assert_eq!(graphs[3].get_fis().unwrap().sg_type, FisType::C4);
        assert!(graphs[3].get_fis().is_none());
        assert_eq!(compute_packing(&mut graphs[4], None), 1);
        let fis_type = graphs[4].get_fis().unwrap().sg_type;
        assert!(fis_type == FisType::NECKTIE || fis_type == FisType::BOWTIE);
        assert!(graphs[4].get_fis().is_none());
    }

}