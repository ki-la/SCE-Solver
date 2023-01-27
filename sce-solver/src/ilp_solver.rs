extern crate rplex;
use rplex::*;
use std::{fs::File, io, io::Write, time::Instant, cmp::{max, min}, collections::{HashMap, HashSet}};
use crate::{data_reduction as dr, search_tree_algo as sta};
use crate::forb_ind_subgraph::splittance_sg;
use crate::graph_trait::Graph;

struct IlpStats {
    file_name: String,              // name of the graph file
    output_file: String,            // name of the output file that these stats are written to
    runtime: f32,                   // how long it took to solve SCE on the graph
    size: usize,                    // |V| (no. of vertices in the graph)
    kernel_size: usize,             // |V'| (no. of vertices after data reduction)
    k: Option<usize>,               // solution size
    constraints_total: usize,       // total no. of constraints in the ILP
    constraints: usize,             // actual no. of constraints (reduced due to data reduction)
}

impl IlpStats {
    fn new(in_file: String, out_file: String, n: usize) -> Self {
        IlpStats {
            file_name: in_file,
            output_file: out_file,
            runtime: 0.0,
            size: n,
            kernel_size: n,
            k: None,
            constraints_total: ((n * (n-1) * (n-2) ) + (n * (n-1))) / 2,
            constraints: 0,
        }
    }

    fn write_stats(&self) -> Result<(), io::Error>{
        let mut file: File;
        if let Ok(f) = File::options().write(true).append(true).open(&self.output_file) {
            // stats file already exists -> open for writing and appending data
            file = f;
        } else {
            file = File::create(&self.output_file)?;
            file.write_all(b"File\tTime (s)\tGraph Size\tKernel Size\tk\tConstraints (total)\tConstraints (needed)\n")?;
        }
        let k = match self.k {
            Some(sol_size) => sol_size.to_string(),
            None => String::from("unsolved"),
        };
        file.write_all(format!("{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\n",
                           self.file_name,
                           self.runtime,
                           self.size,
                           self.kernel_size,
                           k,
                           self.constraints_total,
                           self.constraints).as_bytes())?;
        println!("Statistics written to {:?}", self.output_file);
        Ok(())
    }
}

/// starts ILP solver (separately for every connected component)
/// applies data reduction rules (to decrease the number of variables and constraints)
/// saves the statistics in a .csv file (see write_stats)
pub fn start_ilp_solver(in_file: String, out_file: String, g: &mut impl Graph, config: sta::SolverConfig) -> Vec<(usize, usize)> {
    assert!(config.apply_dr == 1 || config.apply_dr == u32::MAX);    // either initial DR or no DR
    let now = Instant::now();                                 // start timer
    let mut stats = IlpStats::new(in_file, out_file, g.size());
    let mut solution: Vec<(usize,usize)> = Vec::new();
    let comps = g.connected_components(None);

    if comps.len() == 1 {                                            // 'g' is connected
        if splittance_sg(g, None) == 0 {                        // 'g' is a split graph already
            stats.kernel_size -= g.size();
            solution = Vec::new();
        } else {
            if config.apply_dr == 1 {  // initial data reduction
                let (core, periphery) = dr::deg_1_core_rules(g, None);
                println!("Core: {:?}\nPeriphery: {:?}", core, periphery);
                solution = solve_sce_ilp(g, &mut stats, core, periphery);
            } else {  // no initial data reduction
                solution = solve_sce_ilp(g, &mut stats, HashSet::new(), HashSet::new());
            }
        }
    } else {
        for comp in comps {                                // solve comps independently
            if comp.len() < 4 {                                     // trivial component -> remove!
                stats.kernel_size -= comp.len();
                continue;
            }
            let (node_map, mut g_tmp) = g.induced_subgraph(&comp);
            if splittance_sg(&g_tmp, None) == 0 {            // 'g_tmp' is split -> remove!
                stats.kernel_size -= comp.len();
                continue;
            }
            let solution_subgraph: Vec<(usize,usize)>;
            if config.apply_dr == 1 {
                let (core, periphery) = dr::deg_1_core_rules(&mut g_tmp, None);
                println!("Core: {:?}\nPeriphery: {:?}", core, periphery);
                solution_subgraph = solve_sce_ilp(&g_tmp, &mut stats, core, periphery);
            } else {
                solution_subgraph = solve_sce_ilp(&g_tmp, &mut stats, HashSet::new(), HashSet::new());
            }
            if !solution_subgraph.is_empty() {
                for (u, v) in solution_subgraph {
                    solution.push((node_map[&u], node_map[&v]));
                }
            }
        }
    }
    stats.runtime = now.elapsed().as_secs_f32();
    stats.k = Some(solution.len());
    stats.write_stats().unwrap_or_else(|err| {
        println!("Problem writing statistics to file: {err}");
    });

    solution
}

/// returns (non-)edges that have to be flipped to create a split cluster graph
/// uses extern crate rplex to create an ILP that is solved by CPLEX
/// CAUTION: assumes that vertices are numbered 0..g.size()!
#[allow(unused_parens)]
fn solve_sce_ilp(g: &impl Graph, stats: &mut IlpStats, c: HashSet<usize>, p: HashSet<usize>) -> Vec<(usize, usize)> {

    // create a CPLEX environment
    let mut env = Env::new().unwrap();
    env.set_param(EnvParam::Threads(1)).unwrap();  // -> specify #threads used by CPLEX (remove to apply default -> faster!)

    // populate it with a problem
    let mut prob = Problem::new(&env, "partition_var").unwrap();

    // minimize the objective
    prob.set_objective_type(ObjectiveType::Minimize).unwrap();

    //------------------------- create the variables (incl. constraints) -------------------------//
    // one variable for every vertex
    let mut vertices: Vec<usize> = Vec::with_capacity(g.size());
    for v in 0..g.size() {
        let v_name = String::from("v") + &v.to_string();
        if c.contains(&v) {         // 'v' is a core vertex -> set variable to 1.0
            vertices.push(prob.add_variable(var!(1.0 <= v_name <= 1.0 -> 0.0 as Binary)).unwrap());
        }
        else if p.contains(&v) {    // 'v' is a periphery vertex -> set variable to 0.0
            vertices.push(prob.add_variable(var!(0.0 <= v_name <= 0.0 -> 0.0 as Binary)).unwrap());
        } else {                         // 'v' is still unlabeled -> add binary variable
            vertices.push(prob.add_variable(var!(v_name -> 0.0 as Binary)).unwrap());
        }
    }
    // one variable for every vertex pair (edges and non-edges)
    let number_vertex_pairs = (g.size() * (g.size() - 1)) / 2;
    let mut vertex_pairs: Vec<(usize,usize)> = Vec::with_capacity(number_vertex_pairs);
    let mut edge_dict: HashMap<(usize,usize), usize> = HashMap::with_capacity(number_vertex_pairs);
    for u in 0..g.size()-1 {
        for v in u+1..g.size() {
            let e_name = String::from("e") + &u.to_string() + &v.to_string();
            if g.has_edge(u, v) {
                if g.edited(u, v) {     // uv is a permanent edge -> set variable to 1.0
                    edge_dict.insert((u,v),prob.add_variable(var!(1.0 <= e_name <= 1.0 -> (-1.0) as Binary)).unwrap());
                } else {
                    edge_dict.insert((u,v),prob.add_variable(var!(e_name -> (-1.0) as Binary)).unwrap());
                }
            } else {
                if g.edited(u, v) {     // uv is a forbidden non-edge -> set variable to 0.0
                    edge_dict.insert((u,v),prob.add_variable(var!(0.0 <= e_name <= 0.0 -> (1.0) as Binary)).unwrap());
                } else {
                    edge_dict.insert((u,v),prob.add_variable(var!(e_name -> (1.0) as Binary)).unwrap());
                }
            }
            vertex_pairs.push((u, v));
        }
    }
    let m = g.edge_count();
    let offset = prob.add_variable(var!(1.0 <= "Offset" <= 1.0 -> (m as f64) as Integer)).unwrap();

    //------------------------------------- add constraints --------------------------------------//
    let mut peri_constr = 0;                       // counter to count the periphery constraints
    let mut core_constr = 0;                       // counter to count the (core) constraints
    let mut cnt1 = 0;                              // counter to numerate the constraints
    for u in 0..g.size()-1 {
        for v in u+1..g.size() {
            let con_nm1 = String::from("c") + &cnt1.to_string();       // constraint's name
            let var_u:usize = vertices[u];
            let var_v:usize = vertices[v];
            if let Some(var_uv) = edge_dict.get(&(u, v)) {
                let var_uv:usize = *var_uv;
                if !c.contains(&u) && !c.contains(&v) &&
                    (g.has_edge(u, v) || !g.edited(u, v)) {
                    // periphery vertices are an independent set
                    //prob.add_constraint(con!(con_nm1: (0.0) <= (1.0) var_u + (1.0) var_v + (-1.0) var_uv)).unwrap();
                    prob.add_constraint(con!(con_nm1: (1.0) <= (1.0) var_u + (1.0) var_v + (-1.0) var_uv + (1.0) offset)).unwrap();
                    peri_constr += 1;
                    // (this constraint is skipped if 'u' or 'v' are labeled 'core' already or if
                    // uv is a forbidden non-edge since the constraint is trivially fulfilled then)
                }
                if p.contains(&u) || p.contains(&v) {
                    // 'u' or 'v' are labeled 'periphery' already -> no core constraints needed!
                    continue;
                }
                let mut cnt2 = 0;
                for w in 0..g.size() {
                    if w != u && w != v {
                        let con_nm2 = String::from("c") + &cnt1.to_string() + &cnt2.to_string();
                        cnt2 += 1;
                        let a:(usize, usize) = (min(u, w), max(u, w));
                        let b:(usize, usize) = (min(v, w), max(v, w));
                        if let Some(var_a) = edge_dict.get(&a) {
                            let var_a:usize = *var_a;
                            if let Some(var_b) = edge_dict.get(&b) {
                                let var_b:usize = *var_b;
                                // core vertices in the same connected component form a clique (no induced P_3)
                                //prob.add_constraint(con!(con_nm2: (-3.0) <= (1.0) var_uv + (-1.0) var_a + (-1.0) var_b + (-1.0) var_u + (-1.0) var_v)).unwrap();
                                prob.add_constraint(con!(con_nm2: (1.0) <= (1.0) var_uv + (-1.0) var_a + (-1.0) var_b + (-1.0) var_u + (-1.0) var_v + (4.0) offset)).unwrap();
                                core_constr += 1;
                            }
                        }
                    }
                }
            }
            cnt1 += 1;
        }
    }

    stats.constraints += core_constr + peri_constr;

    // prob.write("sce.lp").unwrap();  // write the ILP to a file

    //------------------------------------ solve the problem -------------------------------------//
    println!("solving the ILP...");
    let sol = prob.solve().unwrap();
    //println!("{:?}", sol);

    //--------------------------------- calculate solution set -----------------------------------//
    let mut v_idx = 0;
    let mut e_idx = 0;
    let mut core_vertices: Vec<usize> = Vec::with_capacity(g.size());
    let mut periphery_vertices: Vec<usize> = Vec::with_capacity(g.size());
    let mut modified_edges: Vec<(usize,usize)> = Vec::new();
    for s in sol.variables {
        if v_idx < vertices.len() {
            if let VariableValue::Binary(value) = s {
                if value {
                    core_vertices.push(v_idx);
                } else {
                    periphery_vertices.push(v_idx);
                }
            }
            v_idx += 1;
        } else {
            if e_idx == vertex_pairs.len() {
                break; // skip variable 'offset'
            }
            if let VariableValue::Binary(value) = s {
                if (value && !(g.has_edge(vertex_pairs[e_idx].0, vertex_pairs[e_idx].1))) ||
                    (!(value) && g.has_edge(vertex_pairs[e_idx].0, vertex_pairs[e_idx].1)) {
                    modified_edges.push(vertex_pairs[e_idx]);
                }
                e_idx += 1;
            }
        }
    }
    assert_eq!(sol.objective.round() as usize, modified_edges.len());
    //println!("core: {:?}, periphery: {:?}", core_vertices, periphery_vertices);

    modified_edges
}