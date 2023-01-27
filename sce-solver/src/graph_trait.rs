use std::collections::HashMap;
use ahash::{AHashSet, AHashMap};
use indexmap::IndexSet;
use crate::forb_ind_subgraph::FIS;

#[derive(Clone, Debug, PartialEq)]
pub enum VertexType {
	Unknown,
	Core,
	Periphery,
}

#[derive(Clone, Debug, PartialEq)]
pub enum EdgeType {
	Unedited,
	Edited,
}

pub trait Graph {
	/// creates a new graph with 'n' vertices and no edges
	fn new(n: usize) -> Self;
	/// returns the size (i.e. number of vertices) of the graph
	fn size(&self) -> usize;
	/// returns number of edges in the graph (in O(1) time)
	fn edge_count(&self) -> usize;
	/// labels vertex 'v' as 'vt' (Unknown, Core, or Periphery)
	fn label_vertex(&mut self, v: usize, vt: VertexType);
	/// returns (reference to) label (Unknown, Core, or Periphery) of vertex 'v'
	fn get_label(&self, v: usize) -> &VertexType;
	/// returns whether vertex 'v' is labeled (as Core or Periphery) already
	fn labeled(&self, v: usize) -> bool {
		*self.get_label(v) != VertexType::Unknown
	}
	// the three functions below can be made faster (O(1)) by storing the number of labeled vertices
	/// returns the number of labeled vertices
	fn labeled_vertices(&self) -> usize {
		let mut labeled = 0;
		for v in 0..self.size() {  					// -> O(n) time
			if *self.get_label(v) != VertexType::Unknown {
				labeled += 1;
			}
		}
		labeled
	}
	/// returns the number of vertices labeled as Core
	fn labeled_core(&self) -> usize {
		let mut core_vertices = 0;
		for v in 0..self.size() {  					// -> O(n) time
			if *self.get_label(v) == VertexType::Core {
				core_vertices += 1;
			}
		}
		core_vertices
	}
	/// returns the number of vertices labeled as Periphery
	fn labeled_periphery(&self) -> usize {
		let mut periphery_vertices = 0;
		for v in 0..self.size() {						// -> O(n) time
			if *self.get_label(v) == VertexType::Periphery {
				periphery_vertices += 1;
			}
		}
		periphery_vertices
	}
	/// returns degree (i.e., number of neighbors) of vertex 'v' (in the original graph)
	fn degree(&self, v: usize) -> usize;
	/// returns number of neighbors that 'v' has in 'set'
	fn degree_in_set(&self, v: usize, set: &IndexSet<usize>) -> usize;
	/// adds edge {u,v}
	fn add_edge(&mut self, u: usize, v: usize);
	/// returns whether edge {u,v} is present in the graph
	fn has_edge(&self, u: usize, v: usize) -> bool;
	/// marks (non-)edge {u,v} as 'edited'
	fn edit_edge(&mut self, u: usize, v: usize);
	/// returns whether edge {u,v} has been edited before
	fn edited(&self, u: usize, v: usize) -> bool;
	/// marks an edge as 'not edited'
	fn unedit_edge(&mut self, u: usize, v: usize);
	/// flips the (non-)edge {u,v} (edge becomes non-edge and the other way around)
	fn flip_edge(&mut self, u: usize, v: usize) {
		if self.has_edge(u, v) {
			self.remove_edge(u, v);
		} else {
			self.add_edge(u, v);
		}
	}
	/// removes the edge {u,v}
	fn remove_edge(&mut self, u: usize, v: usize);
	/// returns a list of vertex v's neighbors
	fn neighbors(&self, v: usize, vert: Option<&[usize]>) -> Vec<usize>;
	/// returns an iterator over vertex v's neighbors
	fn neighbors_iter(&self, v: usize) -> Box<dyn Iterator<Item=usize> + '_>;
	/// returns non-decreasing (non-increasing if !incr) degree sequence and the corresponding
	/// vertex ordering (if ord==true) as (degrees, ordering)
	/// (if vert.is_some(), then only those vertices are regarded)
	fn degree_sequence(&self, vert: Option<&[usize]>, ord: bool, incr: bool) -> (Vec<usize>, Vec<usize>) {
		let mut deg_seq: Vec<usize>;
		let mut deg_ord: Vec<usize> = Vec::new();
		if let Some(vt) = vert {
			deg_seq = vec![0; vt.len()];
			if ord { deg_ord = Vec::with_capacity(vt.len()); }
			for u in 0..vt.len() {
				if ord { deg_ord.push(u); }  // ord contains index in 'vert'
				for v in u+1..vt.len() {								// O(n²)
					if self.has_edge(vt[u], vt[v]) {
						deg_seq[u] += 1;
						deg_seq[v] += 1;
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
			if ord { deg_ord.sort_by_key(|v| deg_seq[*v]); }				// sort by degree  (unstable would be deterministic as well)
			deg_seq.sort_unstable();
		} else {																// non-increasing
			if ord { deg_ord.sort_by(|a, b| deg_seq[*b].cmp(&deg_seq[*a])); }  // (unstable would be deterministic as well)
			deg_seq.sort_unstable_by(|a, b| b.cmp(a));
		}
		if let Some(vt) = vert {
			deg_ord = deg_ord.iter().map(|v| vt[*v]).collect();  // map index to actual vertex
		}

		(deg_seq, deg_ord)
	}
	/// computes and returns (degree_buckets, degrees, min_degree) in O(n) time
	///  - degree_buckets(d) = vertices in 'self' w/ degree d
	///  - degrees(i) = deg(i) for every vertex i
	///  - min_degree = min{deg(v) | v \in V}
	/// if comp.is_some(), then only this component is regarded
	/// CAUTION: 'comp' has to be a single connected component of 'g' w/o edges to outside vertices!
	fn degree_sequence_buckets(&self, comp: Option<&AHashSet<usize>>) -> (Vec<IndexSet<usize>>, Vec<usize>, usize) {
		// create buckets (degree_buckets[i] = vertices w/ degree i):
		let mut degree_buckets: Vec<IndexSet<usize>>;
		// create vector that maps a vertex to its bucket, i.e., degree (degrees[i] = deg(i)):
		let mut degrees: Vec<usize> = vec![0; self.size()];
		let mut min_degree: usize = self.size();
		if let Some(c) = comp {
			degree_buckets = vec![IndexSet::new(); c.len()]; // initialize w/ empty sets
			for v in c {
				if self.degree(*v) < min_degree {
					min_degree = self.degree(*v);
				}
				degrees[*v] = self.degree(*v);
				if let Some(bucket) = degree_buckets.get_mut(self.degree(*v)) {
					bucket.insert(*v);
				}
			}
		} else {
			// create buckets and initialize w/ empty sets -> degree_buckets[i] = vertices in 'self' w/ degree i
			degree_buckets = vec![IndexSet::new(); self.size()];
			for v in 0..self.size() {
				if self.degree(v) < min_degree {
					min_degree = self.degree(v);
				}
				degrees[v] = self.degree(v);
				if let Some(bucket) = degree_buckets.get_mut(self.degree(v)) {
					bucket.insert(v);
				}
			}
		}
		(degree_buckets, degrees, min_degree)
	}

	/// returns vector of connected components
	/// each component is represented as a vector of the vertices it contains
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
				for nb_u in self.neighbors(u, vert) {
					if !seen.contains(&nb_u) {
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
	/// returns vector of connected components
	/// each component is represented as an AHashSet of the vertices it contains
	fn connected_components_sets(&self, vert: Option<&[usize]>) -> Vec<AHashSet<usize>> {
		let all_vert = (0..self.size()).collect::<Vec<usize>>();
		let act_vert = match vert {
			Some(v) => v,
			None => &all_vert,
		};
		let mut components:Vec<AHashSet<usize>> = Vec::new();
		let mut seen: AHashSet<usize> = AHashSet::with_capacity(act_vert.len());
		for v in act_vert {
			if seen.contains(v) { continue; }

			seen.insert(*v);
			let mut comp = AHashSet::from([*v]);

			// find all other vertices in component c (using bfs)
			let mut queue = vec![*v];
			while !queue.is_empty() {
				let u = queue.pop().unwrap();
				for nb_u in self.neighbors(u, vert) {
					if !seen.contains(&nb_u) {
						seen.insert(nb_u);
						queue.push(nb_u);
						comp.insert(nb_u);
					}
				}
			}
			components.push(comp);
		}
		components
	}
	/// finds and returns the shortest s-t path in the graph (or an Error if there is none)
	fn shortest_path(&self, s: usize, t: usize, vert: Option<&[usize]>) -> Result<Vec<usize>, &'static str>;
	/// returns graph that is induced by vertices in array 'comp'
	fn induced_subgraph(&self, comp: &[usize]) -> (HashMap<usize, usize>, Self);
	/// returns graph that is induced by vertices in AHashSet 'comp'
	fn induced_subgraph_set(&self, comp: &AHashSet<usize>) -> (AHashMap<usize, usize>, Self);
	fn store_packing(&mut self, p: Vec<FIS>);
	/// returns a reference to the vector containing a FIS packing
	fn packing_ref(&self) -> &Vec<FIS>;
	/// removes and returns last fis in the packing
	fn get_fis(&mut self) -> Option<FIS>;
	/// removes and returns fis at index 'i' in the packing
	/// (returns None if packing.is_empty() and panics if i >= packing.len())
	fn remove_fis(&mut self, i: usize) -> Option<FIS>;
	fn print_graph(&self);

}

#[cfg(test)]
mod tests {
	use ahash::AHashSet;
	use crate::graph_parser::{read_graph_from_file, read_pgraph_from_file, read_pgraph_map_from_file};
	use crate::graph_trait::{Graph, VertexType};

	fn setup1() -> Vec<impl Graph> {
		let g01 = read_graph_from_file(&"../test_instances/g03.dimacs").unwrap();
		let g02 = read_graph_from_file(& "../test_instances/g02.dimacs").unwrap();
		let g03 = read_graph_from_file(& "../test_instances/g01.dimacs").unwrap();
		let g04 = read_graph_from_file(& "../test_instances/g04.dimacs").unwrap();
		let g05 = read_graph_from_file(& "../test_instances/g05.dimacs").unwrap();
		vec![g01, g02, g03, g04, g05]
	}

	fn setup2() -> Vec<impl Graph> {
		let g01 = read_pgraph_from_file(&"../test_instances/g03.dimacs").unwrap();
		let g02 = read_pgraph_from_file(& "../test_instances/g02.dimacs").unwrap();
		let g03 = read_pgraph_from_file(& "../test_instances/g01.dimacs").unwrap();
		let g04 = read_pgraph_from_file(& "../test_instances/g04.dimacs").unwrap();
		let g05 = read_pgraph_from_file(& "../test_instances/g05.dimacs").unwrap();
		vec![g01, g02, g03, g04, g05]
	}

	fn setup3() -> Vec<impl Graph> {
		let g01 = read_pgraph_map_from_file(&"../test_instances/g03.dimacs").unwrap();
		let g02 = read_pgraph_map_from_file(& "../test_instances/g02.dimacs").unwrap();
		let g03 = read_pgraph_map_from_file(& "../test_instances/g01.dimacs").unwrap();
		let g04 = read_pgraph_map_from_file(& "../test_instances/g04.dimacs").unwrap();
		let g05 = read_pgraph_map_from_file(& "../test_instances/g05.dimacs").unwrap();
		vec![g01, g02, g03, g04, g05]
	}


	#[test]
	fn test_basics_all_graphs() {
		let mut graphs = setup1();
		test_basics(&mut graphs);
		let mut graphs = setup2();
		test_basics(&mut graphs);
		let mut graphs = setup3();
		test_basics(&mut graphs);
	}
	
	fn test_basics(graphs: &mut Vec<impl Graph>) {
		//graphs[3].print_graph();
		assert_eq!(graphs[3].size(), 10);
		assert_eq!(graphs[3].edge_count(), 14);
		assert_eq!(*graphs[3].get_label(7), VertexType::Unknown);
		assert!(!graphs[3].labeled(7));
		assert_eq!(graphs[3].labeled_vertices(), 0);
		graphs[3].label_vertex(7, VertexType::Periphery);
		assert_eq!(*graphs[3].get_label(7), VertexType::Periphery);
		assert!(graphs[3].labeled(7));
		assert_eq!(graphs[3].labeled_vertices(), 1);
		assert_eq!(graphs[3].labeled_periphery(), 1);
		assert_eq!(graphs[3].labeled_core(), 0);
		graphs[3].label_vertex(7, VertexType::Core);
		assert!(graphs[3].labeled_core() == 1 && graphs[3].labeled_periphery() == 0);
		graphs[3].label_vertex(7, VertexType::Unknown);
		assert_eq!(graphs[3].labeled_vertices(), 0);
		assert_eq!(graphs[3].degree(7), 1);
		assert_eq!(graphs[3].degree(0), 5);
		assert!(graphs[3].has_edge(0, 4));
		assert!(!graphs[3].has_edge(0, 1));
		graphs[3].add_edge(0, 1);
		assert!(graphs[3].has_edge(0, 1));
		assert_eq!(graphs[3].degree(0), 6);
		graphs[3].remove_edge(0, 4);
		assert!(!graphs[3].has_edge(0, 4));
		graphs[3].flip_edge(0, 4);
		assert!(graphs[3].has_edge(0, 4));
		assert_eq!(graphs[3].degree(0), 6);
		assert_eq!(graphs[3].edge_count(), 15);
		graphs[3].edit_edge(0, 1);  // edge -> permanent
		graphs[3].edit_edge(0, 2);  // non-edge -> forbidden
		assert!(graphs[3].edited(0, 1));
		assert!(graphs[3].edited(0, 2));
		graphs[3].unedit_edge(0, 1);
		graphs[3].unedit_edge(0, 2);
		assert!(!graphs[3].edited(0, 1));
		assert!(!graphs[3].edited(0, 2));
	}

	#[test]
	fn test_nbs() {
		let graphs = setup1();
		let mut nbs = graphs[3].neighbors(0, None);
		nbs.sort();
		assert_eq!(nbs, vec![4,5,6,7,9]);
		nbs = graphs[3].neighbors(0, Some(&vec![0,1,2,3,4,5,6]));
		nbs.sort();
		assert_eq!(nbs, vec![4,5,6]);
		assert_eq!(graphs[2].neighbors(1, None), vec![2]);
		assert_eq!(graphs[2].neighbors(1, Some(&vec![0,1,3,4,5,6,7])), Vec::new());
		assert_eq!(graphs[2].neighbors(4, None).len(), 3);
		assert_eq!(graphs[4].neighbors(0, None).len(), 131);
	}

	#[test]
	fn test_degree_sequence() {
		let graphs = setup1();
		let (d, o) = graphs[3].degree_sequence(None, false, true);
		assert_eq!(d, vec![1,2,2,2,2,3,3,3,5,5]);
		assert!(o.is_empty());
		let (d, o) = graphs[3].degree_sequence(Some(&vec![0,4,7,1,5]), true, false);
		assert_eq!(d, vec![3,2,1,1,1]);
		assert!(o[0] == 0 && o[1] == 5);
		let (d, _) = graphs[2].degree_sequence(None, true, true);
		assert_eq!(d, vec![0,0,1,1,2,3,3,3,4,5]);
		let (d, o) = graphs[0].degree_sequence(Some(&vec![0,1,6,7,9]), true, false);
		assert_eq!(d, vec![4,3,2,2,1]);
		assert!(o == vec![0,7,6,9,1] || o == vec![0,7,9,6,1]);
		let (d, mut o) = graphs[0].degree_sequence(Some(&vec![9,1,4,2]), true, true);
		assert_eq!(d, vec![0,0,0,0]);
		o.sort_unstable();
		assert_eq!(o, vec![1,2,4,9]);
		let (d, _) = graphs[4].degree_sequence(None, true, true);
		assert_eq!(d.iter().sum::<usize>(), graphs[4].edge_count() * 2); // Handshake Lemma
	}

	#[test]
	fn test_degree_buckets() {
		let graphs = setup3();
		let (buckets, degrees, min_deg) = graphs[3].degree_sequence_buckets(None);
		assert!(buckets[4].is_empty());
		for i in 0..graphs[3].size() {
			assert_eq!(graphs[3].degree(i), degrees[i]);
			assert!(buckets[graphs[3].degree(i)].contains(&i));
		}
		assert_eq!(min_deg, 1);

		let (buckets, degrees, min_deg) = graphs[2].degree_sequence_buckets(None);
		assert_eq!(buckets[0].len(), 2);
		for i in 0..graphs[2].size() {
			assert_eq!(graphs[2].degree(i), degrees[i]);
			assert!(buckets[graphs[2].degree(i)].contains(&i));
		}
		assert_eq!(min_deg, 0);

		// two connected components from g01:
		let comp1 = AHashSet::from([1,2]);
		let comp2 = AHashSet::from([4,5,6,7,8,9]);
		let (buckets, degrees, min_deg) = graphs[2].degree_sequence_buckets(Some(&comp1));
		for i in comp1 {
			assert_eq!(graphs[2].degree(i), degrees[i]);
			assert!(buckets[graphs[2].degree(i)].contains(&i));
		}
		assert_eq!(min_deg, 1);
		let (buckets, degrees, min_deg) = graphs[2].degree_sequence_buckets(Some(&comp2));
		for i in comp2 {
			assert_eq!(graphs[2].degree(i), degrees[i]);
			assert!(buckets[graphs[2].degree(i)].contains(&i));
		}
		assert_eq!(min_deg, 2);
		assert_eq!(degrees, vec![0,0,0,0,3,4,3,5,2,3]);  // CAUTION: v1 and v2 have deg-1, but are not part of the considered component!
	}

	#[test]
	fn test_connected_components() {
		let mut graphs = setup1();
		let mut comps = graphs[3].connected_components(None);
		assert_eq!(comps.len(), 1);
		graphs[3].remove_edge(0, 5);
		comps = graphs[3].connected_components(None);
		assert_eq!(comps.len(), 2);
		let mut comp1 = comps[0].clone();
		let mut comp2 = comps[1].clone();
		comp1.sort();
		comp2.sort();
		assert!((comp1 == vec![0,4,6,7,9] && comp2 == vec![1,2,3,5,8])
			|| (comp2 == vec![0,4,6,7,9] && comp1 == vec![1,2,3,5,8]));
		comps = graphs[3].connected_components(Some(&vec![1,2,3,4,5,6,7,9]));
		assert_eq!(comps.len(), 3);
		assert_eq!(graphs[2].connected_components(None).len(), 4);
	}

	#[test]
	fn test_connected_components_sets() {
		let mut graphs = setup1();
		let mut comps = graphs[3].connected_components_sets(None);
		assert_eq!(comps.len(), 1);
		graphs[3].remove_edge(0, 5);
		comps = graphs[3].connected_components_sets(None);
		assert_eq!(comps.len(), 2);
		assert!((comps[0] == AHashSet::from([0,4,6,7,9]) && comps[1] == AHashSet::from([1,2,3,5,8]))
			|| (comps[1] == AHashSet::from([0,4,6,7,9]) && comps[0] == AHashSet::from([1,2,3,5,8])));
		comps = graphs[3].connected_components_sets(Some(&vec![1,2,3,4,5,6,7,9]));
		assert_eq!(comps.len(), 3);
		assert_eq!(graphs[2].connected_components_sets(None).len(), 4);
	}

	#[test]
	fn test_shortest_path() {
		let graphs = setup1();
		assert_eq!(graphs[1].shortest_path(0, 5, None).unwrap(), vec![0,1,2,3,4,5]);
		assert!(graphs[1].shortest_path(0,9, Some(&vec![0,1,2,3,8,9])).is_err());
		assert!(graphs[2].shortest_path(0, 6, None).is_err());
	}

	#[test]
	fn test_induced_subgraph() {
		let graphs = setup3();
		let (_, sg) = graphs[2].induced_subgraph(&[1,5,8,9]);
		assert_eq!(sg.edge_count(), 2);
		assert_eq!(sg.size(), 4);
	}

	#[test]
	fn test_induced_subgraph_set() {
		let graphs = setup3();
		let (_, sg) = graphs[2].induced_subgraph_set(&AHashSet::from([1,5,8,9]));
		assert_eq!(sg.edge_count(), 2);
		assert_eq!(sg.size(), 4);
	}
}
