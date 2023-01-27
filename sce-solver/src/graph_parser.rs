use std::fs::File;
use std::io::BufRead;
use std::io;
use crate::graph::GraphAM;
use crate::graph_trait::Graph;
use crate::petgraph::PGraph;
use crate::petgraph_hybrid::PGraphMap;

fn error_msg() -> &'static str {
        "Too few/ wrong arguments!\
        \nUsage:\tcargo run -- <path_to_graph_file> <config.txt> -o <output_file.csv>\n\
        OR:\tcargo run -- <path_to_graph_file> <timeout> <OPTIONS_1_to_6> -o <output_file.csv>\n\
        Options: \n 1: solver (0: search-tree, 1: ilp, 2: heuristic)\
                 \n 2: graph (0: graph, 1: pgraph, 2: pgraph_map)\
                 \n 3: branching (0: branch until UB reached, 1: k-loop)\
                 \n 4: apply packing lower bound every x'th recursive step (0: never)\
                 \n 5: apply minimum cut lower bound every x'th recursive step (0: never)\
                 \n 6: apply data reduction every x'th recursive step (0: never)\
                 \n -o <name_output_file> is optional! (default: stats.csv)"
}

pub fn parse_input(args: &[String]) -> Result<(String, String, Vec<u32>, f32), &'static str> {
    if args.len() < 3 {
        return Err(error_msg());
    }

    let mut opt: Vec<u32> = vec![0, 2, 0, 1, 1, 1, 0];  // default config
    let timeout: f32;
    let graph_file = args[1].to_owned();

    // third argument: either config file...
    if args[2].ends_with(".txt") && args[2].contains("config") {
        match read_config(&args[2]) {
            Ok((o, t)) => (opt, timeout) = (o, t),
            Err(_) => return Err("Error reading config file!"),
        }
    } else {
        // ... or timeout (followed by other arguments)
        if args.len() < 9 {
            return Err(error_msg());
        }
        let mut args_iter = args.iter().skip(2);
        timeout = args_iter.next().unwrap().parse::<f32>().unwrap();
        opt[0] = args_iter.next().unwrap().parse::<u32>().unwrap();
        opt[1] = args_iter.next().unwrap().parse::<u32>().unwrap();
        opt[2] = args_iter.next().unwrap().parse::<u32>().unwrap();
        opt[3] = args_iter.next().unwrap().parse::<u32>().unwrap();
        opt[4] = args_iter.next().unwrap().parse::<u32>().unwrap();
        let dr = args_iter.next().unwrap().parse::<i32>().unwrap();
        if dr >= 0 {
            opt[5] = dr as u32;
            opt[6] = 0;
        } else if dr == -1 {
            opt[5] = 0;
            opt[6] = 1;
        } else {
            return Err("Invalid argument!")
        }
    }
    for i in 3..6 {
        if opt[i] == 0 {  // means: no application of LB/DR
            opt[i] = u32::MAX;
        }
    }

    if timeout < 0.0 || timeout > 90000.0 {
        return Err("Timeout should be between 0 and 90000 seconds!");
    }

    let mut out_file = String::from("../results/");
    let default_out_name = String::from("stats.csv");
    if args[args.len()-2] == "-o" && args[args.len()-1].ends_with(".csv") {  // out_file given!
        out_file += args.get(args.len()-1).unwrap_or_else(|| &default_out_name);
    } else {
        out_file += &default_out_name;
    }

    Ok((graph_file, out_file, opt, timeout))
}

/// reads config file and returns solver options and timeout
fn read_config(file: &str) -> Result<(Vec<u32>, f32), io::Error> {
    let mut opt: Vec<u32> = vec![0, 2, 0, 1, 1, 1, 0];                  // default config
    let mut timeout: f32 = 300.0;                                       // default timeout
    let f: File = File::open(file)?;                               // open file or return error
    let mut dr: i32 = 1;
    let lines = io::BufReader::new(f).lines();
    for line in lines {
        let line = line?;
        if !line.is_empty() && !line.starts_with('#') {
            let data: Vec<&str> = line.split(':').collect();
            match data[0] {
                "timeout"               => timeout = data[1].parse().unwrap(),
                "solver"                => opt[0] = data[1].parse().unwrap(),
                "graph_data_structure"  => opt[1] = data[1].parse().unwrap(),
                "branching"             => opt[2] = data[1].parse().unwrap(),
                "apply_packing_lb"      => opt[3] = data[1].parse().unwrap(),
                "apply_mincut_lb"       => opt[4] = data[1].parse().unwrap(),
                "apply_data_reduction"  => dr = data[1].parse().unwrap(),
                _ => return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Wrong format of config file!")),
            }
        }
    }
    if dr >= 0 {
        opt[5] = dr as u32;
        opt[6] = 0;
    } else if dr == -1 {
        opt[5] = 0;
        opt[6] = 1;
    } else {
        return Err(io::Error::new(io::ErrorKind::InvalidData,"Invalid argument"));
    }
    Ok((opt,timeout))
}


fn read_lines(file: &str) -> Result<(io::Lines<io::BufReader<File>>, usize), io::Error> {
    let f = File::open(file)?;                              // open file or return error
    let mut lines = io::BufReader::new(f).lines();
    let mut n: Option<usize> = None;
    //let mut m: Option<usize> = None;
    for line in lines.by_ref() {
        let line = line?;
        match line.bytes().next().unwrap() {
            b'c' => continue,                                        // ignore comments
            b'p' => {
                // descriptor line --> dimacs format from PACE 2021
                let mut split = line.split(' ');
                n = split.nth(2).and_then(|s| s.parse().ok());
                //m = split.nth(3).and_then(|s| s.parse().ok());
                break;
            },
            b'#' => {
                // number of vertices --> format from other instances
                let mut split = line.split(' ');
                n = split.nth(1).and_then(|s| s.parse().ok());
                break;
            },
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Did not find number of vertices"))
            },
        }
    }

    let n = match n {
        Some(no_vert) => no_vert,
        None => return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Could not read number of vertices")),
    };

    Ok((lines, n))
}

fn create_graph(g: &mut impl Graph, lines: io::Lines<io::BufReader<File>>) -> Result<i32, io::Error> {
    let mut max_id = 0;
    let n = g.size();
    for line in lines {
        let line = line?;
        if line.starts_with('c') || line.starts_with('#') || line.is_empty() {
            continue;                                                  // ignore comments
        }
        let edge: Vec<&str> = line.split(' ').collect();
        if edge.len() == 2 {
            let u: usize = (*edge.get(0).unwrap()).parse::<usize>().unwrap();
            let v: usize = (*edge.get(1).unwrap()).parse::<usize>().unwrap();
            if u > n || v > n {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Incorrect number of nodes"));
            }
            g.add_edge(u-1, v-1);  // index shift --> first node has ID 0
            if v > max_id { max_id = v; }
            if u > v { max_id = u; }
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid graph format"));
        }
    }
    if n < max_id {  // n > max_id is no problem -> isolated nodes
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Incorrect number of nodes"));
    }
    Ok(0)
}

pub fn read_graph_from_file(file: &str) -> Result<impl Graph, io::Error> {
    let (lines, n) = read_lines(file).unwrap();
    let mut g: GraphAM = Graph::new(n);
    create_graph(&mut g, lines).unwrap();
    Ok(g)
}

pub fn read_pgraph_from_file(file: &str) -> Result<impl Graph, io::Error> {
    let (lines, n) = read_lines(file).unwrap();
    let mut g: PGraph = Graph::new(n);
    create_graph(&mut g, lines).unwrap();
    Ok(g)
}

pub fn read_pgraph_map_from_file(file: &str) -> Result<impl Graph, io::Error> {
    let (lines, n) = read_lines(file).unwrap();
    let mut g: PGraphMap = Graph::new(n);
    create_graph(&mut g, lines).unwrap();
    Ok(g)
}