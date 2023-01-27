use std::collections::HashMap;
use ahash::{AHashMap, AHashSet};
use indexmap::IndexSet;
use petgraph::algo::astar;
use petgraph::visit::NodeIndexable;
use crate::forb_ind_subgraph::FIS;
use crate::graph_trait::{EdgeType, Graph, VertexType};

type PGraphType = petgraph::graphmap::UnGraphMap<usize, EdgeType>;

pub struct PGraphMap {
    g: PGraphType,
    edited: Vec<Vec<bool>>,
    fis_packing: Vec<FIS>,
    vertices: Vec<VertexType>,
    labeled_core: usize,
    labeled_periphery: usize,
    degrees: Vec<usize>,
}

impl Graph for PGraphMap {
    fn new(n: usize) -> Self {
        assert!(n > 0);
        let mut g_struct = PGraphMap {
            g: PGraphType::new(),
            edited: vec![vec![false; n]; n],
            fis_packing: Vec::new(),
            vertices: vec![VertexType::Unknown; n],
            labeled_core: 0,
            labeled_periphery: 0,
            degrees: vec![0; n],
        };
        for v in 0..n {
            g_struct.g.add_node(v);
        }
        g_struct
    }

    fn size(&self) -> usize {
        self.g.node_count()
    }

    fn edge_count(&self) -> usize {
        self.g.edge_count()
    }

    fn label_vertex(&mut self, v: usize, vt: VertexType) {
        assert_ne!(self.vertices[v], vt);
        match self.vertices[v] {
            VertexType::Core => self.labeled_core -= 1,
            VertexType::Periphery => self.labeled_periphery -= 1,
            _ => {},
        }

        self.vertices[v] = vt;
        match self.vertices[v] {
            VertexType::Core => self.labeled_core += 1,
            VertexType::Periphery => self.labeled_periphery += 1,
            _ => {},
        }
    }

    fn get_label(&self, v: usize) -> &VertexType {
        &self.vertices[v]
    }

    fn labeled_vertices(&self) -> usize {
        self.labeled_core + self.labeled_periphery
    }

    fn labeled_core(&self) -> usize {
        self.labeled_core
    }

    fn labeled_periphery(&self) -> usize {
        self.labeled_periphery
    }

    fn degree(&self, v: usize) -> usize {
        self.degrees[v]
        //self.g.neighbors(self.g.from_index(v)).count() // (slightly) slower
    }

    fn degree_in_set(&self, v: usize, set: &IndexSet<usize>) -> usize {
        let mut deg_in_set = 0;
        for u in self.g.neighbors(v) {
            if set.contains(&u) {
                deg_in_set += 1;
            }
        }
        deg_in_set
    }

    fn add_edge(&mut self, u: usize, v: usize) {
        let e = self.g.add_edge(u, v, EdgeType::Unedited);
        assert!(e.is_none());  // ensure edge didn't exist yet
        self.degrees[u] += 1;
        self.degrees[v] += 1;
    }

    fn has_edge(&self, u: usize, v: usize) -> bool {
        self.g.contains_edge(u, v)
    }

    fn edit_edge(&mut self, u: usize, v: usize) {
        debug_assert!(!self.edited[u][v]);
        debug_assert_eq!(self.edited[u][v], self.edited[u][v]);
        self.edited[u][v] = true;
        self.edited[v][u] = true;
    }

    fn edited(&self, u: usize, v: usize) -> bool {
        debug_assert_eq!(self.edited[u][v], self.edited[v][u]);
        self.edited[u][v]
    }

    fn unedit_edge(&mut self, u: usize, v: usize) {
        debug_assert!(self.edited[u][v]);
        debug_assert_eq!(self.edited[u][v], self.edited[u][v]);
        self.edited[u][v] = false;
        self.edited[v][u] = false;
    }

    fn remove_edge(&mut self, u: usize, v: usize) {
        let e = self.g.remove_edge(u, v);
        assert!(e.is_some());  // ensure edge did exist
        self.degrees[u] -= 1;
        self.degrees[v] -= 1;
    }

    fn neighbors(&self, v: usize, vert: Option<&[usize]>) -> Vec<usize> {
        if let Some(vertices) = vert {
            assert!(vertices.contains(&v));
            return self.g.neighbors(v)
                .filter(|v| vertices.contains(v))
                .collect::<Vec<usize>>();
        }
        self.g.neighbors(v)
            .collect()
    }

    fn neighbors_iter(&self, v: usize) -> Box<dyn Iterator<Item=usize> + '_> {
        Box::new(self.g.neighbors(v))
    }

    fn degree_sequence(&self, vert: Option<&[usize]>, ord: bool, incr: bool) -> (Vec<usize>, Vec<usize>) {
        let mut deg_seq: Vec<usize>;
        let mut deg_ord: Vec<usize> = Vec::new();
        if let Some(vt) = vert {
            let mut vt_set = AHashSet::with_capacity(vt.len());  // HashSet allows for faster 'contains()' !!
            for elem in vt {
                vt_set.insert(*elem);
            }
            deg_seq = vec![0; vt.len()];
            if ord { deg_ord = Vec::with_capacity(vt.len()); }
            for u in 0..vt.len() {
                if ord { deg_ord.push(u); }  // ord contains index in 'vert'
                for v in self.g.neighbors(vt[u]) {                      // O(n²) (faster for sparse graphs)
                    if vt_set.contains(&v) {
                        deg_seq[u] += 1;
                    }
                }
            }
        } else {
            deg_seq = Vec::with_capacity(self.size());
            if ord { deg_ord = Vec::with_capacity(self.size()); }
            for u in 0..self.size() {									// O(n)
                deg_seq.push(self.degree(u));
                if ord { deg_ord.push(u); }
            }
        }

        // sort degree sequence and ordering -> O(n log n)
        if incr {																// non-decreasing
            if ord { deg_ord.sort_by_key(|v| deg_seq[*v]); }				    // sort by degree  (unstable would be deterministic as well)
            deg_seq.sort_unstable();
        } else {																// non-increasing
            if ord { deg_ord.sort_by(|a, b| deg_seq[*b].cmp(&deg_seq[*a])); }   // (unstable would be deterministic as well)
            deg_seq.sort_unstable_by(|a, b| b.cmp(a));
        }
        if let Some(vt) = vert {
            deg_ord = deg_ord.iter().map(|v| vt[*v]).collect();  // map index to actual vertex
        }

        (deg_seq, deg_ord)
    }

    fn connected_components(&self, vert: Option<&[usize]>) -> Vec<Vec<usize>> {
        let all_vert = (0..self.size()).collect::<Vec<usize>>();
        let act_vert = match vert {
            Some(v) => v,
            None => &all_vert,
        };
        let mut components:Vec<Vec<usize>> = Vec::new();
        let mut seen: AHashSet<usize> = AHashSet::with_capacity(act_vert.len());
        for v in act_vert {
            if seen.contains(v) { continue; }

            seen.insert(*v);
            let mut comp = vec![*v];

            // find all other vertices in component c (using bfs)
            let mut queue = vec![*v];
            while !queue.is_empty() {
                let u = queue.pop().unwrap();
                for nb_u in self.neighbors_iter(u) {
                    if !seen.contains(&nb_u) {
                        if vert.is_some() && !act_vert.contains(&nb_u) {
                            continue;
                        }
                        seen.insert(nb_u);
                        queue.push(nb_u);
                        comp.push(nb_u);
                    }
                }
            }
            components.push(comp);
        }
        components
    }

    fn shortest_path(&self, s: usize, t: usize, vert: Option<&[usize]>) -> Result<Vec<usize>, &'static str> {
        if vert.is_none() {
            if let Some(path) = astar(&self.g,
                                      s,
                                      |finish| finish == t,
                                      |_| 1,
                                      |_| 0) {
                return Ok(path.1.iter().map(|v| self.g.to_index(*v)).collect());
            }
            Err("no shortest path from {s} to {t} found!")
        } else {
            if s == t {
                return Ok(vec![s]);
            }
            let mut succ: HashMap<usize, usize> = HashMap::new();
            succ.insert(t, t);                                // t is its own successor
            let mut q: Vec<usize> = vec![t];
            let mut q_idx:usize = 0;
            while q_idx < q.len() {
                let curr_v = q[q_idx];
                for v in self.neighbors(curr_v, vert) {
                    if !succ.contains_key(&v) {                 // v was not visited yet
                        succ.insert(v, curr_v);                 // curr_v is v's successor
                        if v == s {                                // source reached -> create path
                            let mut p: usize = succ[&v];
                            let mut path: Vec<usize> = vec![s, p];
                            while p != t {
                                p = succ[&p];
                                path.push(p);
                            }
                            return Ok(path);
                        }
                        q.push(v);
                    }
                }
                q_idx += 1;
            }
            Err("no shortest path from {s} to {t} found!")
        }
    }

    fn induced_subgraph(&self, comp: &[usize]) -> (HashMap<usize, usize>, Self) {
        let mut g_out: PGraphMap = Graph::new(comp.len());
        let mut node_map: HashMap<usize,usize> = HashMap::with_capacity(comp.len());
        let mut next_id: usize = 0;
        let mut seen: HashMap<usize,usize> = HashMap::with_capacity(comp.len());
        for u in comp.iter() {
            if !seen.contains_key(u) {
                node_map.insert(next_id, *u);
                seen.insert(*u, next_id);
                next_id += 1;
            }
            for v in self.neighbors(*u, Some(comp)).iter() {
                if *v > *u {
                    if !seen.contains_key(v) {
                        node_map.insert(next_id, *v);
                        seen.insert(*v, next_id);
                        next_id += 1;
                    }
                    g_out.add_edge(seen[u], seen[v]);
                }
            }
        }
        (node_map, g_out)
    }

    // CAUTION -> iteration order over HashSet varies in every execution -> node_map looks different
    // every time -> to achieve determinism, 'comp' is transformed into a sorted vector -> this
    // increases runtime and should be undone if determinism is not required!
    fn induced_subgraph_set(&self, comp: &AHashSet<usize>) -> (AHashMap<usize, usize>, Self) {
        let mut g_out: PGraphMap = Graph::new(comp.len());
        let mut node_map: AHashMap<usize,usize> = AHashMap::with_capacity(comp.len());
        let mut next_id: usize = 0;
        let mut seen: AHashMap<usize,usize> = AHashMap::with_capacity(comp.len());
        let mut comp_arr = comp.into_iter().collect::<Vec<&usize>>();
        comp_arr.sort_unstable();  // ensure that g_out always looks the same!
        for u in comp_arr {
            if !seen.contains_key(u) {
                node_map.insert(next_id, *u);
                seen.insert(*u, next_id);
                next_id += 1;
            }
            for v in self.g.neighbors(*u) {
                if v > *u && comp.contains(&v) {
                    if !seen.contains_key(&v) {
                        node_map.insert(next_id, v);
                        seen.insert(v, next_id);
                        next_id += 1;
                    }
                    g_out.add_edge(seen[u], seen[&v]);
                }
            }
        }
        (node_map, g_out)
    }

    fn store_packing(&mut self, p: Vec<FIS>) {
        self.fis_packing = p;
    }

    fn packing_ref(&self) -> &Vec<FIS> {
        &self.fis_packing
    }

    fn get_fis(&mut self) -> Option<FIS> {
        self.fis_packing.pop()  // O(1)
    }

    fn remove_fis(&mut self, i: usize) -> Option<FIS> {
        if self.fis_packing.is_empty() {
            return None;
        }
        Some(self.fis_packing.swap_remove(i))  // removes i'th fis and swaps last one to index i -> O(1) time!
    }

    fn print_graph(&self) {
        println!("{:?}", petgraph::dot::Dot::new(&self.g));
    }
}