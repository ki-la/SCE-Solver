use std::collections::HashMap;
use ahash::{AHashMap, AHashSet};
use indexmap::IndexSet;
use crate::forb_ind_subgraph::FIS;
use crate::graph_trait::{EdgeType, Graph, VertexType};

#[derive(Clone)]
struct Edge {
    is_edge: bool,
    e_type: EdgeType,
}

pub struct GraphAM {
    n: usize,
    m: usize,
    adj_matrix: Vec<Vec<Edge>>,
    vertices: Vec<VertexType>,
    labeled_core: usize,
    labeled_periphery: usize,
    degrees: Vec<usize>,
    fis_packing: Vec<FIS>,
}

impl Graph for GraphAM {
    fn new(n: usize) -> Self {
        assert!(n > 0);
        GraphAM {
            n,
            m: 0,
            adj_matrix: vec![vec![Edge {is_edge: false, e_type: EdgeType::Unedited}; n]; n],
            vertices: vec![VertexType::Unknown; n],
            labeled_core: 0,
            labeled_periphery: 0,
            degrees: vec![0; n],
            fis_packing: Vec::new(),
        }
    }

    fn size(&self) -> usize {
        self.n
    }

    fn edge_count(&self) -> usize {
        self.m
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
    }

    fn degree_in_set(&self, v: usize, set: &IndexSet<usize>) -> usize {
        let mut deg_in_set = 0;
        for u in set {
            if self.has_edge(*u, v) {
                deg_in_set += 1;
            }
        }
        deg_in_set
    }

    fn add_edge(&mut self, u: usize, v: usize) {
        debug_assert!(!self.has_edge(u, v));
        self.adj_matrix[u][v].is_edge = true;
        self.adj_matrix[v][u].is_edge = true;
        self.degrees[u] += 1;
        self.degrees[v] += 1;
        self.m += 1;
    }

    fn has_edge(&self, u:usize, v:usize) -> bool {
        debug_assert_eq!(self.adj_matrix[u][v].is_edge, self.adj_matrix[v][u].is_edge);
        self.adj_matrix[u][v].is_edge
    }

    fn edit_edge(&mut self, u:usize, v:usize) {
        debug_assert_eq!(self.adj_matrix[u][v].e_type, self.adj_matrix[v][u].e_type);
        self.adj_matrix[u][v].e_type = EdgeType::Edited;
        self.adj_matrix[v][u].e_type = EdgeType::Edited;
    }

    fn edited(&self, u:usize, v:usize) -> bool {
        debug_assert_eq!(self.adj_matrix[u][v].e_type, self.adj_matrix[v][u].e_type);
        self.adj_matrix[u][v].e_type == EdgeType::Edited
    }

    fn unedit_edge(&mut self, u:usize, v:usize) {
        debug_assert_eq!(self.adj_matrix[u][v].e_type, self.adj_matrix[v][u].e_type);
        self.adj_matrix[u][v].e_type = EdgeType::Unedited;
        self.adj_matrix[v][u].e_type = EdgeType::Unedited;
    }

    fn remove_edge(&mut self, u: usize, v: usize) {
        debug_assert!(self.has_edge(u, v));
        self.adj_matrix[u][v].is_edge = false;
        self.adj_matrix[v][u].is_edge = false;
        self.degrees[u] -= 1;
        self.degrees[v] -= 1;
        self.m -= 1;
    }

    fn neighbors(&self, v: usize, vert: Option<&[usize]>) -> Vec<usize> {
        let mut nbs: Vec<usize> = Vec::with_capacity(self.degree(v));
        if let Some(vt) = vert {
            debug_assert!(vt.contains(&v));
            for u in vt {
                if self.has_edge(*u, v) {
                    nbs.push(*u);
                }
            }
        } else {
            for u in 0..self.size() {
                if self.has_edge(u, v) {
                    nbs.push(u);
                }
            }
            debug_assert_eq!(nbs.len(), self.degree(v));
        }
        nbs
    }

    fn neighbors_iter(&self, v: usize) -> Box<dyn Iterator<Item=usize> + '_> {
        Box::new(self.neighbors(v, None).into_iter())
    }

    fn shortest_path(&self, s: usize, t: usize, vert: Option<&[usize]>) -> Result<Vec<usize>, &'static str> {
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

    fn induced_subgraph(&self, comp: &[usize]) -> (HashMap<usize, usize>, Self) {
        let mut g_out: GraphAM = Graph::new(comp.len());
        let mut node_map: HashMap<usize,usize> = HashMap::with_capacity(comp.len());
        for (i, u) in comp.iter().enumerate() {
            node_map.insert(i, *u);
            for (j, v) in comp.iter().enumerate().skip(i+1) {
                if self.has_edge(*u, *v) {
                    g_out.add_edge(i, j);
                }
            }
        }
        (node_map, g_out)
    }

    // CAUTION -> iteration order over HashSet varies in every execution -> node_map looks different
    // every time -> to achieve determinism, 'comp' is transformed into a sorted vector -> this
    // increases runtime and should be undone if determinism is not required!
    fn induced_subgraph_set(&self, comp: &AHashSet<usize>) -> (AHashMap<usize, usize>, Self) {
        let mut g_out: GraphAM = Graph::new(comp.len());
        let mut node_map: AHashMap<usize,usize> = AHashMap::with_capacity(comp.len());
        let mut comp_arr = comp.into_iter().collect::<Vec<&usize>>();
        comp_arr.sort_unstable();  // ensure that g_out always looks the same!
        for (i, u) in comp_arr.iter().enumerate() {
            node_map.insert(i, **u);
            for (j, v) in comp_arr.iter().enumerate().skip(i+1) {
                if self.has_edge(**u, **v) {
                    g_out.add_edge(i, j);
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
        self.fis_packing.pop()
    }

    fn remove_fis(&mut self, i: usize) -> Option<FIS> {
        if self.fis_packing.is_empty() {
            return None;
        }
        Some(self.fis_packing.swap_remove(i))  // removes i'th fis, swaps last fis to idx i -> O(1) time!
    }

    fn print_graph(&self) {  // prints graph in DOT format
        println!("graph {{");
        for vertex in 0..self.size() {
            print!("{vertex} [ label = {vertex}");
            match self.vertices[vertex] {
                VertexType::Core => print!(", color = red ]"),
                VertexType::Periphery => print!(", color = blue ]"),
                VertexType::Unknown => print!(" ]"),
            };
            println!();
        }
        for u in 0..self.size() {
            for v in u+1..self.size() {
                if self.adj_matrix[u][v].is_edge {
                    print!("{u} -- {v}");
                    if self.adj_matrix[u][v].e_type == EdgeType::Edited {
                        print!(" [ label = Edited ]");
                    }
                    println!();
                }
            }
        }
        println!("}}");
    }
}