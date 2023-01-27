use indexmap::{IndexMap, IndexSet};
use std::{fs, fs::File, io, io::Write, time::Instant};
use std::cmp::max;
use crate::forb_ind_subgraph::{compute_packing, compute_solution_from_splittance, splittance, splittance_global};
use crate::graph_trait::Graph;
use crate::min_cut as mc;
use rand::{Rng, SeedableRng};

struct HeuristicStats {
    file_name: String,              // name of the graph file
    output_file: String,            // name of the output file that these stats are written to
    runs: usize,                    // number of runs (repetitions of the heuristic)
    max_steps: usize,               // max. number of steps executed in every run of the heuristic
    runtime: f32,                   // how long it took to solve SCE heuristically
    size: usize,                    // |V| (no. of vertices in the graph)
    sol_size: usize,                // solution size of the heuristic solution
    k: Option<usize>,               // solution size
    splittance_ub: usize,           // global splittance (splittance calculated for every comp.)
    solved_optimally: bool,         // states whether problem was solved optimally (if OPT known)
}

impl HeuristicStats {
    fn new(in_file: String, out_file: String, g: &impl Graph, r: usize, steps: usize) -> Self {
        HeuristicStats {
            file_name: in_file,
            output_file: out_file,
            runs: r,
            max_steps: steps,
            runtime: 0.0,
            size: g.size(),
            sol_size: g.edge_count(),  // E is always a solution set for SCE
            k: None,
            splittance_ub: splittance_global(g),
            solved_optimally: false,
        }
    }

    fn write_stats(&self) -> Result<(), io::Error> {
        let mut file: File;
        if let Ok(f) = File::options().write(true).append(true).open(&self.output_file) {
            // stats file already exists -> open for writing and appending data
            file = f;
        } else {
            file = File::create(&self.output_file)?;
            file.write_all(b"File\tRuns\tMax. Steps\tTime (s)\tGraph Size\tk\tHeuristic (UB)\tSplittance (UB)\tOptimal?\n")?;
        }
        let k = match self.k {
            Some(sol_size) => sol_size.to_string(),
            None => String::from("unknown"),
        };
        file.write_all(format!("{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\n",
                               self.file_name,
                               self.runs,
                               self.max_steps,
                               self.runtime,
                               self.size,
                               k,
                               self.sol_size,
                               self.splittance_ub,
                               self.solved_optimally).as_bytes())?;
        println!("Statistics written to {:?}", self.output_file);
        Ok(())
    }
}

pub struct Clusters {
    pub clusters: Vec<IndexSet<usize>>,     // list of all clusters
    pub cluster: IndexMap<usize, usize>,    // maps vertex to cluster it is contained in
}

/// calculates the value of 'cluster' as the number of kept edges (from the original graph) minus
/// the necessary insertions, i.e., #edges_within_cluster - splittance(cluster)
/// c.f. 'SCEclusterValue(...)' in 'Heuristic.cc' in (BHK15)
fn value_cluster(g: &impl Graph, cluster: &IndexSet<usize>) -> usize {
    let mut degrees = vec![0;cluster.len()];
    let mut sum_degrees = 0;
    for (i, v) in cluster.iter().enumerate() {
        degrees[i] = g.degree_in_set(*v, cluster);
        sum_degrees += degrees[i];
    }

    let mut s = 0;  // splittance
    if cluster.len() >= 4 {
        degrees.sort_unstable();
        degrees.reverse();
        s = splittance(&degrees);
    }

    (sum_degrees / 2) - s
}

/// calculates the costs (i.e., number of edits) to split graph 'g' into the clusters 'c'
/// calculation: #edges - #kept_edges + splittance(c) f.a. clusters c
/// c.f. 'numEditsSCE()' in 'Heuristic.cc' in (BHK15)
fn clusters_costs(g: &impl Graph, c: &Clusters) -> usize {
    let mut costs = 0;
    for c in &c.clusters {
        costs += value_cluster(g, c);
    }
    g.edge_count() - costs
}

#[allow(unused)]
/// creates an initial solution that places every vertex in 'g' in its own cluster
fn initial_solution(g: &impl Graph) -> Clusters {
    let mut cs = Vec::with_capacity(g.size());
    let mut c_map = IndexMap::new();
    for v in 0..g.size() {
        cs.push(IndexSet::from([v]));
        c_map.insert(v, v);
    }
    Clusters {
        clusters: cs,
        cluster: c_map,
    }
}

#[allow(unused)]
/// creates an initial solution placing a vertex in a cluster w/ all its (yet unassigned) neighbors
fn initial_solution2(g: &impl Graph) -> Clusters {
    let mut cs = Vec::with_capacity(g.size());
    let mut c_map = IndexMap::new();
    let mut c_idx = 0;
    let mut vertex_count = 0;
    for v in 0..g.size() {
        if vertex_count == g.size() {
            break;
        }
        if !c_map.contains_key(&v) {
            cs.push(IndexSet::from([v]));
            c_map.insert(v, c_idx);
            vertex_count += 1;
            for u in g.neighbors_iter(v) {
                if !c_map.contains_key(&u) {
                    cs[c_idx].insert(u);
                    c_map.insert(u, c_idx);
                    vertex_count += 1;
                }
            }
            c_idx += 1;
        }
    }

    // insert enough empty clusters s.t. every vertex could potentially have one
    while cs.len() < g.size() {
        cs.push(IndexSet::new());
    }

        Clusters {
        clusters: cs,
        cluster: c_map,
    }
}

#[allow(unused)]
/// creates an initial solution placing every vertex in 'g' in one cluster (-> value = splittance)
fn initial_solution3(g: &impl Graph) -> Clusters {
    let mut cs = Vec::with_capacity(g.size());
    cs.push((0..g.size()).collect::<IndexSet<usize>>());  // one cluster w/ all vertices
    let mut c_map = IndexMap::new();
    c_map.insert(0, 0);  // (avoid adding too many empty clusters in the loop)
    for v in 1..g.size() {
        c_map.insert(v, 0);  // every vertex is part of the first cluster
        // insert enough empty clusters s.t. every vertex could potentially have one:
        cs.push(IndexSet::new());
    }
    Clusters {
        clusters: cs,
        cluster: c_map,
    }
}

/// moves vertex 'v' from cluster 'from' to cluster 'to' and returns the delta (value(new cluster)
/// - value(old cluster)) -> a positive delta implies an improvement (decrease of costs)
fn move_vertex(g: &impl Graph, c: &mut Clusters, v: usize, from: usize, to: usize) -> f64 {
    let val_old = value_cluster(g, &c.clusters[from]) + value_cluster(g, &c.clusters[to]);
    debug_assert_eq!(*c.cluster.get(&v).unwrap(), from);
    c.clusters[from].remove(&v);
    c.clusters[to].insert(v);
    let f = c.cluster.insert(v, to);
    debug_assert!(f.is_some());
    debug_assert_eq!(f.unwrap(), from);
    let val_new = value_cluster(g, &c.clusters[from]) + value_cluster(g, &c.clusters[to]);
    val_new as f64 - val_old as f64
}

/// moves vertex 'v' back to its old cluster 'old_c'
fn move_vertex_back(c: &mut Clusters, v: usize, old_c: usize, new_c: usize) {
    debug_assert_eq!(*c.cluster.get(&v).unwrap(), new_c);
    c.clusters[new_c].remove(&v);
    c.clusters[old_c].insert(v);
    let f = c.cluster.insert(v, old_c);
    debug_assert!(f.is_some());
    debug_assert_eq!(f.unwrap(), new_c);
}

/// applies 'max_steps' of simulated annealing procedure on graph 'g' starting with clusters 'c'
/// stops early if timeout 't' many seconds have passed
/// cf. 'simulatedAnnealing(...)' in 'Heuristic.cc' in (BHK15)
fn simulated_annealing(g: &impl Graph, t: f32, c: &mut Clusters, max_steps: usize, max_kt: f64, seed: u64) -> usize {
    let start_time = Instant::now();

    let mut poss_moves: Vec<(usize, usize)> = Vec::new();   // all possible moves
    for u in 0..g.size() {
        for v in g.neighbors_iter(u) {
            poss_moves.push((u, v));                  // "move u to v's cluster"
        }
        poss_moves.push((u, g.size()));               // means "move u to empty cluster"
    }

    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let uni_distr = rand::distributions::Uniform::new(0, poss_moves.len());
    let std_distr = rand::distributions::Standard;

    for step in 0..max_steps {
        if start_time.elapsed().as_secs_f32() > t {
            println!("HEURISTIC: TIMEOUT!");
            break;
        }
        let kt: f64 = max_kt * ((max_steps - step) as f64 / max_steps as f64);  // decrease temp.
        let (u, v) = poss_moves[rng.sample(uni_distr)];      // choose random move
        let cluster_old = c.cluster.get(&u).unwrap().clone();
        let mut cluster_new: Option<usize> = None;
        if v == g.size() {                                              // move 'v' to empty cluster
            if c.clusters[cluster_old].len() == 1 {                     // 'v' is alone already
                continue;
            }
            for c_idx in 0..c.clusters.len() {
                if c.clusters[c_idx].is_empty() {                       // this is eventually true!
                    cluster_new = Some(c_idx);
                    break
                }
            }
        } else {
            cluster_new = Some(c.cluster.get(&v).unwrap().clone());
        }
        assert!(cluster_new.is_some());
        let cluster_new = cluster_new.unwrap();
        if cluster_new == cluster_old {                                 // no move
            continue;
        }
        let delta = move_vertex(g, c, u, cluster_old, cluster_new);
        let r: f64 = rng.sample(std_distr);
        if delta < 0.0 && r >= std::f64::consts::E.powf(delta / kt) {
            // a worse solution (delta < 0) is not kept w/ increasing probability
            move_vertex_back(c, u, cluster_old, cluster_new);
        }
    }
    clusters_costs(g, c)
}

/// cf. 'solve(...)' in 'Heuristic.cc' in (BHK15)
/// (not used in search_tree_algo anymore)
pub fn solve_heuristically(g: &impl Graph, timeout: f32, runs: usize, max_steps: usize) -> Vec<(usize, usize)> {
    let start_time = Instant::now();
    assert!(runs > 0 && max_steps > 0);
    let mut min_sol_size: Option<usize> = None;
    let mut best_cluster: Option<Clusters> = None;

    for r in 0..runs {
        if start_time.elapsed().as_secs_f32() > timeout {
            println!("HEURISTIC: TIMEOUT!");
            break;
        }
        let mut clusters =
        if r % 2 == 0 {
            initial_solution(g)
        } else {
            initial_solution2(g)
        };
        let sol_size = simulated_annealing(g,
                                           timeout - start_time.elapsed().as_secs_f32(),
                                           &mut clusters,
                                           max_steps,
                                           0.1,
                                           123456 * (r+1) as u64);  // change seed here
        if min_sol_size.is_none() || sol_size < min_sol_size.unwrap() {
            min_sol_size = Some(sol_size);
            best_cluster = Some(clusters);
        }
    }
    assert!(min_sol_size.is_some() && best_cluster.is_some());
    let best_cluster = best_cluster.unwrap();

    compute_edits(g, &best_cluster, min_sol_size.unwrap())
}

/// runs the heuristic and returns solution size and clusters found
///     -> edits still have to be computed!
/// - stops early if optimal solution is found (sol == lb)
/// - aims to improve upper bound 'ub' -> increases steps if heur_sol > ub (splittance)
pub fn run_heuristic(g: &impl Graph, timeout: f32, runs: usize, steps: usize, lb: usize, ub: usize) -> (usize, Clusters) {
    assert!(timeout < 600.0);  // make sure heuristic does not run 'too long'
    let start_time = Instant::now();
    let mut rem_runs = runs;  // remaining runs (might not be decreased in every iteration)
    let mut max_steps: usize = steps;
    let mut min_sol_size: Option<usize> = None;
    let mut best_cluster: Option<Clusters> = None;
    let mut runs_cntr: u32 = 0;  // counter for the total runs so far

    while rem_runs > 0 {
        if start_time.elapsed().as_secs_f32() > timeout {
            println!("HEURISTIC: TIMEOUT!");
            break;
        }

        let mut clusters =
            if (runs_cntr % 3 == 2) && min_sol_size.is_some() && min_sol_size.unwrap() > ub  {
                initial_solution3(g)
            } else if runs_cntr % 2 == 0 {
                initial_solution(g)
            } else {
                initial_solution2(g)
            };

        let sol_size = simulated_annealing(g,
                                           timeout - start_time.elapsed().as_secs_f32(),
                                           &mut clusters,
                                           max_steps,
                                           0.1,
                                           123456 * (runs_cntr+1) as u64);  // change seed here!
        if min_sol_size.is_none() || sol_size < min_sol_size.unwrap() {
            min_sol_size = Some(sol_size);
            best_cluster = Some(clusters);
        }

        if sol_size == lb {  // OPT found! :)
            break;
        }

        if min_sol_size.unwrap() > ub {
            // run longer and increase max_steps to (hopefully) improve below known upper bound
            max_steps += 200;
        } else {
            rem_runs -= 1;
        }

        runs_cntr += 1;
    }
    debug_assert!(min_sol_size.is_some() && best_cluster.is_some());
    let best_cluster = best_cluster.unwrap();

    (min_sol_size.unwrap(), best_cluster)
}

/// computes necessary edits to transform graph 'g' into a split cluster graph w/ the clusters 'c'
pub fn compute_edits(g: &impl Graph, c: &Clusters, size: usize) ->  Vec<(usize, usize)> {
    let mut edits:  Vec<(usize, usize)> = Vec::with_capacity(size);
    for cluster in &c.clusters {
        if !cluster.is_empty() {
            if cluster.len() > 3 {  // smaller clusters always have 0-splittance!!
                let cluster_vec = cluster.iter().cloned().collect::<Vec<usize>>();
                edits.append(&mut compute_solution_from_splittance(g, Some(&cluster_vec)));
            }
            for u in cluster {
                for v in g.neighbors_iter(*u) {
                    if *u < v && !cluster.contains(&v) {  // u < v -> add every pair only once
                        edits.push((*u, v));
                    }
                }
            }
        }
    }
    assert_eq!(edits.len(), size);
    edits
}

/// starts heuristic solver and writes the statistics to a file
pub fn start_heuristic_solver(in_file: String, out_file: String, timeout: f32, g: &mut impl Graph, runs: usize, max_steps: usize) -> Vec<(usize, usize)>{
    let mut stats = HeuristicStats::new(in_file.clone(), out_file, g, runs, max_steps);
    let now = Instant::now();

    //let edits = solve_heuristically(g, runs, max_steps);
    let ub = splittance_global(g);
    let packing_lb = compute_packing(g, None);
    let mc_lb = mc::min_cut_lower_bound2_global(g, Some(ub), Some(ub));
    let lb = max(packing_lb, mc_lb);
    let (sol_size, c) = run_heuristic(g, timeout, runs, max_steps, lb, ub);
    let edits = compute_edits(g, &c, sol_size);

    stats.sol_size = edits.len();
    stats.runtime = now.elapsed().as_secs_f32();

    // apply edits to the graph
    for (u, v) in &edits {
        g.flip_edge(*u, *v);
    }
    // make sure that the graph is a split cluster graph now
    assert_eq!(splittance_global(g), 0);

    // check for solution file and save OPT
    let sol_file = in_file + ".solution";
    if let Ok(sol) = fs::read_to_string(sol_file) {
        let mut k = 0;
        for line in sol.lines() {
            if !line.is_empty() {
                k += 1;
            }
        }
        stats.k = Some(k);
        assert!(k <= edits.len());
        if k == edits.len() {
            stats.solved_optimally = true;
        }
    } else {
        println!("No .solution file found.");
    }

    stats.write_stats().unwrap_or_else(|err| {
        println!("Problem writing statistics to file: {err}");
    });

    println!("HEURISTIC SOLUTION: k = {}\nEdits:", edits.len());
    for (u, v) in &edits {
        if g.has_edge(*u, *v) {
            println!("Del({}, {})", *u + 1, *v + 1);
        } else {
            println!("Ins({}, {})", *u + 1, *v + 1);
        }
    }

    edits
}

#[cfg(test)]
mod heuristic_tests {
    use crate::graph_parser::*;
    use crate::graph_trait::Graph;
    use crate::heuristic::{compute_edits, run_heuristic, solve_heuristically};
    use crate::forb_ind_subgraph::splittance_global;

    #[allow(unused)]
    fn setup() -> Vec<impl Graph> {
        let files = ["../test_instances/final/random_ER/random_10_3_1.dimacs",
                            "../test_instances/final/random_SBM/sbm_18_4_80_20.dimacs",
                            "../test_instances/final/50_50/cog_31_1602.dimacs",
                            "../test_instances/final/50_50/cog_49_508.dimacs",
                            "../test_instances/final/50_50/cog_50_523.dimacs",
                            "../test_instances/final/50_50/cog_50_981.dimacs",
        ];
        let mut graphs = Vec::new();
        for f in files {
            graphs.push(read_pgraph_map_from_file(f).unwrap()); // change graph type here!
        }
        graphs
    }

    #[test]
    fn test_solve_heuristically() {
        let graphs = setup();
        let edits = solve_heuristically(&graphs[0], 10.0,5, 500);
        assert!(edits.len() >= 4);  // random_10_3_1.dimacs has OPT solution k=4

        let edits = solve_heuristically(&graphs[1], 10.0, 5, 500);
        assert!(edits.len() >= 19);  // sbm_18_4_80_20.dimacs has OPT solution k=19

        for mut graph in graphs {
            let edits = solve_heuristically(&graph, 60.0, 3, 700);
            for (u, v) in edits {
                graph.flip_edge(u, v);
            }
            assert_eq!(splittance_global(&graph), 0);  // check if solution produces valid SCE
        }
    }
    #[test]
    fn test_run_heuristic() {
        let graphs = setup();
        for mut g in graphs {
            let ub = splittance_global(&g);
            let (sol_size, c) = run_heuristic(&g,10.0, 5, 1000, 0, ub);
            let edits = compute_edits(&g, &c, sol_size);
            for (u, v) in edits {
                g.flip_edge(u, v);
            }
            assert_eq!(splittance_global(&g), 0);  // check if solution produces valid SCE
        }
    }
}