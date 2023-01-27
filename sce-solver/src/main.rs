use std::{env, process, fs};
use sce_solver::{ilp_solver, search_tree_algo, graph_parser, graph_trait::Graph, heuristic};
use sce_solver::forb_ind_subgraph::splittance_global;

fn main() {
    let args: Vec<String> = env::args().collect();
    let (graph_file, out_file, opt, timeout) = graph_parser::parse_input(&args).unwrap_or_else(|err| {
        println!("Problem parsing input: {err}");
        process::exit(1);
    });
    let config = search_tree_algo::SolverConfig::new(&opt);
    if opt[1] == 1 {  // use PGraph data structure
        let mut g = graph_parser::read_pgraph_from_file(&graph_file).unwrap_or_else(|err| {
            println!("Problem reading file {graph_file}: {err}");
            process::exit(1);
        });
        start_solver(&mut g, &opt, graph_file, out_file, timeout, config);
    } else if opt[1] == 2 {  // use PGraphMap data structure
        let mut g = graph_parser::read_pgraph_map_from_file(&graph_file).unwrap_or_else(|err| {
            println!("Problem reading file {graph_file}: {err}");
            process::exit(1);
        });
        start_solver(&mut g, &opt, graph_file, out_file, timeout, config);
    } else {  // default: use own graph data structure
        let mut g = graph_parser::read_graph_from_file(&graph_file).unwrap_or_else(|err| {
            println!("Problem reading file {graph_file}: {err}");
            process::exit(1);
        });
        start_solver(&mut g, &opt, graph_file, out_file, timeout, config);
    }
}

/// starts solver and prints the solution to the console
fn start_solver(g: &mut impl Graph, opt: &Vec<u32>, in_file: String, out_file: String, t: f32, config: search_tree_algo::SolverConfig) {
    let edits = match opt[0] {
        1 => ilp_solver::start_ilp_solver(in_file.clone(), out_file, g, config),
        2 => heuristic::start_heuristic_solver(in_file.clone(), out_file, t, g, 30, 3000),
        _ => search_tree_algo::solve_independently(in_file.clone(), out_file, t, g, config),
    };

    if opt[0] != 2 {
        print_solution(&in_file, &edits);
    }
}

/// prints solution to SCE problem given a vector of vertex pairs
fn print_solution(file: &str, edits: &[(usize, usize)]) {
    // read graph again (makes sure that it's the original graph)
    let mut g = graph_parser::read_graph_from_file(file).unwrap();

    // print solution (+1 -> undoes index shift)
    println!("SOLUTION:\nk = {}\nEdits:", edits.len());
    for (u, v) in edits {
        if g.has_edge(*u, *v) {
            println!("Del({}, {})", *u + 1, *v + 1);
        } else {
            println!("Ins({}, {})", *u + 1, *v + 1);
        }
    }

    // apply edits to the graph
    for (u, v) in edits {
        g.flip_edge(*u, *v);
    }
    // make sure that the graph is a split cluster graph now
    assert_eq!(splittance_global(&g), 0);

    let sol_file = file.to_owned() + ".solution";
    if let Ok(sol) = fs::read_to_string(sol_file) {
        let mut k = 0;
        for line in sol.lines() {
            if !line.is_empty() {
                k += 1;
            }
        }
        assert_eq!(edits.len(), k);
    } else {
        println!("No .solution file found.");
    }
}