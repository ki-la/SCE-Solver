use std::{io, io::Write, process, time::Instant, fs::File, cmp::{max, min}};
use crate::data_reduction as dr;
use crate::forb_ind_subgraph::*;
use crate::min_cut as mc;
use crate::graph_trait::Graph;
use crate::heuristic;

pub struct SolverConfig {
    graph_ds: u32,          // states which graph data structure has been used
    alt_br: bool,           // states whether alternative branching should be used
    apply_packing_lb: u32,  // how often PackingLB is applied during branching (every xth rec. step)
    apply_mc_lb: u32,       // how often MinCutLB -------- " ---------- (every xth rec. step)
    pub apply_dr: u32,      // how often data reduction -------- " ---------- (every xth rec. step)
    no_dr: bool,            // if true, then no dr is applied at all (not even the fastest ones)
}

impl SolverConfig {
    pub fn new(opt: &Vec<u32>) -> Self {
        assert_ne!(opt[3], 0);
        assert_ne!(opt[4], 0);
        assert_ne!(opt[5], 0);
        SolverConfig {
            graph_ds: opt[1],
            alt_br: opt[2] != 1,
            apply_packing_lb: opt[3],
            apply_mc_lb: opt[4],
            apply_dr: opt[5],
            no_dr: opt[6] == 1,
        }
    }
}

struct Solutions {
    curr_sol: Vec<(usize,usize)>,
    ub: Vec<(usize,usize)>,
}

struct SolverStats {
    file_name: String,              // name of the (input) graph file
    output_file: String,            // name of the output file that these stats are written to
    start_time: Instant,            // is set at start of the solver -> allows to stop
    timeout: f32,                   // timeout for the solver
    runtime: f32,                   // how long it took to solve SCE on the graph
    size: usize,                    // |V| (no. of vertices in the graph)
    kernel_size: usize,             // |V'| (no. of vertices after data reduction)
    k: Option<usize>,               // solution size
    clusters: usize,                // number of clusters in the final split cluster graph
    heur_clusters: Option<heuristic::Clusters>,  // clusters found by heuristic
    rec_steps: usize,               // number of recursive steps
    packing_lb: usize,              // size of a FIS-packing lower bound
    min_cut_lb: usize,              // MinCut -> if it's smaller than the splittance it's a LB
    splittance_ub: usize,           // global splittance (splittance calculated for every comp.)
    heuristic_ub: usize,            // upper bound given by heuristic
    labeled_vert: usize,            // #labeled vertices
}


impl SolverStats {
    fn new(in_file: String, out_file: String, n: usize) -> Self {
        SolverStats {
            file_name: in_file,
            output_file: out_file,
            start_time: Instant::now(),
            timeout: 300.0,
            runtime: 0.0,
            size: n,
            kernel_size: n,
            k: None,
            clusters: 0,
            heur_clusters: None,
            rec_steps: 0,
            packing_lb: 0,
            min_cut_lb: 0,
            splittance_ub: 0,
            heuristic_ub: 0,
            labeled_vert: 0,
        }
    }

    fn write_stats(&self, config: &SolverConfig) -> Result<(), io::Error> {
        let mut file: File;
        if let Ok(f) = File::options().write(true).append(true).open(&self.output_file) {
            // stats file already exists -> open for writing and appending data
            file = f;
        } else {
            file = File::create(&self.output_file)?;
            file.write_all(b"File\tTime (s)\tGraph Size\tKernel Size\tk\tClusters\t\
            Recursive Steps\tPacking Size (LB)\tMinCut (LB)\tSplittance (UB)\tHeuristic (UB)\tLabeled Vertices\n")?;
        }
        let k = match self.k {
            Some(sol_size) => sol_size.to_string(),
            None => String::from("unsolved"),
        };
        file.write_all(format!("{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                           self.file_name, self.runtime, self.size, self.kernel_size, k,
                           self.clusters, self.rec_steps, self.packing_lb, self.min_cut_lb,
                           self.splittance_ub, self.heuristic_ub, self.labeled_vert,
                           config.graph_ds, config.alt_br, config.apply_mc_lb,
                           config.apply_packing_lb, config.apply_dr, config.no_dr).as_bytes())?;
        println!("Statistics written to {:?}", self.output_file);
        Ok(())
    }
}

// -------------------------- hardcoded improved branching rules -------------------------------- //
fn branch_on_c4(c4: &FIS) -> Vec<Vec<(usize, usize)>> {  // branching number: 6
    debug_assert_eq!(c4.sg_type, FisType::C4);
    // C4: (0,1,2,3) --> only trivial edits
    vec![vec![(c4.v[0], c4.v[1])], vec![(c4.v[0], c4.v[2])], vec![(c4.v[0], c4.v[3])],
         vec![(c4.v[1], c4.v[2])], vec![(c4.v[1], c4.v[3])], vec![(c4.v[2], c4.v[3])]]
}

fn branch_on_c5(c5: &FIS) -> Vec<Vec<(usize, usize)>> {  // branching number < 5.855
    debug_assert_eq!(c5.sg_type, FisType::C5);
    // C5: (0,1,2,3,4)
    // trivial_edits = [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)]
    // further_edits = [[(0, 2), (0, 3)], [(0, 2), (2, 4)], [(1, 4), (2, 4)], [(0, 3), (1, 3)],
    // [(1, 3), (1, 4)]]
    vec![vec![(c5.v[0], c5.v[1])], vec![(c5.v[1], c5.v[2])], vec![(c5.v[2], c5.v[3])],
         vec![(c5.v[3], c5.v[4])], vec![(c5.v[0], c5.v[4])],
         vec![(c5.v[0], c5.v[2]), (c5.v[0], c5.v[3])], vec![(c5.v[0], c5.v[2]), (c5.v[2], c5.v[4])],
         vec![(c5.v[1], c5.v[4]), (c5.v[2], c5.v[4])], vec![(c5.v[0], c5.v[3]), (c5.v[1], c5.v[3])],
         vec![(c5.v[1], c5.v[3]), (c5.v[1], c5.v[4])]]
}

fn branch_on_p5(p5: &FIS) -> Vec<Vec<(usize, usize)>> {  // branching number < 5.405
    debug_assert_eq!(p5.sg_type, FisType::P5);
    // P5: 0--1--2--3--4
    // trivial_edits = [(0, 1), (1, 2), (2, 3), (3, 4), (1, 3)]
    // further_edits = [[(0, 2), (0, 3)], [(1, 4), (2, 4)], [(0, 2), (0, 4), (2, 4)]]
    vec![vec![(p5.v[0], p5.v[1])], vec![(p5.v[1], p5.v[2])], vec![(p5.v[2], p5.v[3])],
         vec![(p5.v[3], p5.v[4])], vec![(p5.v[1], p5.v[3])],
         vec![(p5.v[0], p5.v[2]), (p5.v[0], p5.v[3])], vec![(p5.v[1], p5.v[4]), (p5.v[2], p5.v[4])],
         vec![(p5.v[0], p5.v[2]), (p5.v[0], p5.v[4]), (p5.v[2], p5.v[4])]]
}

fn branch_on_bt(bt: &FIS) -> Vec<Vec<(usize, usize)>> {  // branching number < 6.317
    debug_assert_eq!(bt.sg_type, FisType::BOWTIE);
    // bowtie: triangle(0--1--2) + triangle(2--3--4)
    // trivial_edits = [(0, 1), (3, 4), (0, 3), (0, 4), (1, 3), (1, 4)]
    // further_edits = [[(0, 2), (1, 2)], [(2, 3), (2, 4)]]
    vec![vec![(bt.v[0], bt.v[1])], vec![(bt.v[3], bt.v[4])], vec![(bt.v[0], bt.v[3])],
         vec![(bt.v[0], bt.v[4])], vec![(bt.v[1], bt.v[3])], vec![(bt.v[1], bt.v[4])],
         vec![(bt.v[0], bt.v[2]), (bt.v[1], bt.v[2])], vec![(bt.v[2], bt.v[3]), (bt.v[2], bt.v[4])]]
}

fn branch_on_nt(nt: &FIS) -> Vec<Vec<(usize, usize)>> {  // branching number < 5.571
    debug_assert_eq!(nt.sg_type, FisType::NECKTIE);
    // necktie: triangle(0--1--2) + P3(2--3--4)
    // trivial_edits = [(0, 1), (0, 3), (1, 3), (2, 3), (3, 4)]
    // further_edits = [[(2, 4), (0, 4)], [(2, 4), (1, 4)], [(2, 4), (0, 2), (1, 2)], [(0, 2), (1, 2)]]
    vec![vec![(nt.v[0], nt.v[1])], vec![(nt.v[0], nt.v[3])], vec![(nt.v[1], nt.v[3])],
         vec![(nt.v[2], nt.v[3])], vec![(nt.v[3], nt.v[4])],
         vec![(nt.v[2], nt.v[4]), (nt.v[0], nt.v[4])], vec![(nt.v[2], nt.v[4]), (nt.v[1], nt.v[4])],
         vec![(nt.v[2], nt.v[4]), (nt.v[0], nt.v[2]), (nt.v[1], nt.v[2])],
         vec![(nt.v[0], nt.v[2]), (nt.v[1], nt.v[2])]]
}
// ---------------------------------------------------------------------------------------------- //

/// sets some (yet unedited) vertex pairs in a P5, NECKTIE, or C5 to edited for the upcoming
/// branches and returns a vector containing these edges s.t. they can be unedited before returning
/// to the parent branch!
/// this is part of the improved branching strategy described in the thesis
fn reedit_edges(g: &mut impl Graph, fis: &FIS, br: usize, edited: &mut Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    let mut edited_edges: Vec<(usize, usize)> = Vec::with_capacity(2);
    if fis.sg_type == FisType::P5 {
        if br == 6 {  // in the last branch (#7) on a P5, (0,3) and (1,4) can be set to forbidden
            if !g.edited(fis.v[0], fis.v[3]) {
                g.edit_edge(fis.v[0], fis.v[3]);
                edited_edges.push((fis.v[0], fis.v[3]));
            }
            if !g.edited(fis.v[1], fis.v[4]) {
                g.edit_edge(fis.v[1], fis.v[4]);
                edited_edges.push((fis.v[1], fis.v[4]));
            }
        }
    } else if fis.sg_type == FisType::C5 {
        if br == 6 && !g.edited(fis.v[0], fis.v[2]) {  // from now on (0,2) -> forbidden
            g.edit_edge(fis.v[0], fis.v[2]);
            edited_edges.push((fis.v[0], fis.v[2]));
        }
        if br == 7 && !g.edited(fis.v[2], fis.v[4]) {  // now, (2,4) -> forbidden
            g.edit_edge(fis.v[2], fis.v[4]);
            edited_edges.push((fis.v[2], fis.v[4]));
        }
        if br == 8 && !g.edited(fis.v[0], fis.v[3]) {  // now, (0,3) -> forbidden
            g.edit_edge(fis.v[0], fis.v[3]);
            edited_edges.push((fis.v[0], fis.v[3]));
        }
    } else if fis.sg_type == FisType::NECKTIE {
        if br == 6 {  // in the next branch (0,4) and (1,4) can be set to forbidden
            if !g.edited(fis.v[0], fis.v[4]) {
                g.edit_edge(fis.v[0], fis.v[4]);
                edited_edges.push((fis.v[0], fis.v[4]));
            }
            if !g.edited(fis.v[1], fis.v[4]) {
                g.edit_edge(fis.v[1], fis.v[4]);
                edited_edges.push((fis.v[1], fis.v[4]));
            }
        } else if br == 7 {  // in the last branch (8): (2,4) -> forbidden + unedit (0,4), (1,4)
            // only if edited contains (0,4) or (1,4) they have to be unedited again
            if edited.len() > 0 {
                if edited[edited.len() - 1] == (fis.v[0], fis.v[4])
                    || edited[edited.len() - 1] == (fis.v[1], fis.v[4]) {
                    let (u, v) = edited.pop().unwrap();
                    g.unedit_edge(u, v);
                    if edited.len() > 0 {
                        if edited[edited.len() - 1] == (fis.v[0], fis.v[4])
                            || edited[edited.len() - 1] == (fis.v[1], fis.v[4]) {
                            let (u, v) = edited.pop().unwrap();
                            g.unedit_edge(u, v);
                        }
                    }
                }
            }
            if !g.edited(fis.v[2], fis.v[4]) {
                g.edit_edge(fis.v[2], fis.v[4]);
                edited_edges.push((fis.v[2], fis.v[4]));
            }
        }
    }
    edited_edges
}

/// approximates branching number of all forbidden induced subgraphs in the graph's current packing
/// by decreasing an initial score depending on the number of uneditable edges in it
/// -> removes and returns fis with best score from the packing
// (this was not very useful in preliminary tests and is therefore currently not used)
#[allow(unused)]
fn get_best_fis(g: &mut impl Graph) -> Option<FIS> {
    let mut min_score = 10.0;  // initial score -> is definitely changed if packing contains a valid fis
    let mut best_fis_idx = 0;
    for (i, fis) in g.packing_ref().iter().enumerate() {
        if is_valid_fis(g, fis) {
            let edits;
            let mut score;
            match fis.sg_type {
                FisType::C4         => { edits = branch_on_c4(fis);
                                         score = 6.0 },
                FisType::C5         => { edits = branch_on_c5(fis);
                                         score = 5.9 },
                FisType::P5         => { edits = branch_on_p5(fis);
                                         score = 5.5 },
                FisType::BOWTIE     => { edits = branch_on_bt(fis);
                                         score = 6.4 },
                FisType::NECKTIE    => { edits = branch_on_nt(fis);
                                         score = 5.6 },
                _                   => { panic!("this should not happen!") },
            }
            for edges in edits {
                if dr::uneditable(g, &edges, 0, false) {  // (k unimportant if dr == false)
                    score -= 1.0 / (edges.len().pow(2) as f64);
                }
            }
            if score <= min_score {
                min_score = score;
                best_fis_idx = i;
            }
        }
    }
    if min_score == 10.0 {  // no valid fis found!
        return None;
    }

    g.remove_fis(best_fis_idx)
}

/// branches on a given forbidden induced subgraph ('fis'), applies data reduction rules,
/// applies edits, decreases 'k' accordingly, and calls 'branch()' on resulting graph
/// -> returns None if none of the possible edge modifications leads to a valid solution!
fn branch_on_fis(g: &mut impl Graph, fis: &FIS, k: isize, stats: &mut SolverStats, config: &SolverConfig, ub: usize) -> Option<Vec<(usize, usize)>> {
    let edits = match fis.sg_type {
        FisType::C4         => branch_on_c4(fis),
        FisType::C5         => branch_on_c5(fis),
        FisType::P5         => branch_on_p5(fis),
        FisType::BOWTIE     => branch_on_bt(fis),
        FisType::NECKTIE    => branch_on_nt(fis),
        _                   => { println!("No FIS to branch on!"); return None; },
    };

    let mut edited_edges: Vec<(usize, usize)> = Vec::with_capacity(10);
    for (br, mut edges) in edits.into_iter().enumerate() {
        let apply_dr = stats.rec_steps % config.apply_dr as usize == 0;  // if true, then apply more expensive data reduction
        if dr::uneditable(g, &edges, k as usize, apply_dr) {
            // although edits were skipped, some edges could still be marked as edited
            let mut edited = reedit_edges(g, &fis, br, &mut edited_edges);
            edited_edges.append(&mut edited);
            continue;  // do not branch on uneditable edges
        }

        let mut dr_stack = Vec::new();
        for edge in &edges {
            let edge = *edge;
            g.flip_edge(edge.0, edge.1);
            g.edit_edge(edge.0, edge.1);
            if !config.no_dr {
                dr_stack.append(&mut dr::dr_during_branching(g, &edge, apply_dr));
            }
        }
        let k_new = k - edges.len() as isize;
        let sol_new = branch(g, k_new, stats, config, ub);
        if let Some(mut solution) = sol_new {           // solution found
            solution.append(&mut edges);
            return Some(solution);
        } else {                                                // no solution found -> undo changes
            let edges_len = edges.len();
            for edge in &edges {
                let edge = *edge;
                if edges_len > 1 {
                    g.unedit_edge(edge.0, edge.1);
                } else {
                    edited_edges.push((edge.0, edge.1));  // unedit only after sibling branches are done
                }
                g.flip_edge(edge.0, edge.1);
            }
            if edges.len() > 1 {
                // due to the improved branching some edges can be set to edited in certain branches
                let mut edited = reedit_edges(g, &fis, br, &mut edited_edges);
                edited_edges.append(&mut edited);
            }
        }

        dr::undo_dr_during_branching(g, dr_stack);
    }

    for edge in edited_edges {
        g.unedit_edge(edge.0, edge.1);
    }

    None  // no solution found
}

/// Search tree algorithm for Split Cluster Editing
/// -> finds a fis and branches on the possibilities to destroy it -> calls branch_on_fis()
/// -> returns None if g has no solution of size k
fn branch(g: &mut impl Graph, k: isize, stats: &mut SolverStats, config: &SolverConfig, ub: usize) -> Option<Vec<(usize, usize)>> {
    stats.rec_steps += 1;

    if stats.start_time.elapsed().as_secs_f32() > stats.timeout {
        println!("TIMEOUT!");
        exit(g, stats, config);
    }

    if k < 0 {
        return None;
    }

    if stats.rec_steps % config.apply_mc_lb as usize == 0 {
        let mc_lb = mc::min_cut_lower_bound2_global(g, None, Some(ub - k as usize));
        if k < mc_lb as isize {  // k is smaller than lower bound -> no solution
            return None;
        }
    }
    let fis: FIS;
    if stats.rec_steps % config.apply_packing_lb as usize == 0 {
        let packing_lb = update_packing(g, k as usize);
        if packing_lb == 0 {  // no FIS found -> g is a split cluster graph!
            return Some(Vec::new());
        }
        if k < packing_lb as isize {  // k is smaller than lower bound -> no solution
            return None;
        }
        fis = g.get_fis().unwrap();
    } else {
        // check if a valid FIS is contained in g.fis_packing...
        let mut fis_opt: Option<FIS> = g.get_fis();
        while fis_opt.is_some() {
            let fis_tmp = fis_opt.as_ref().unwrap();
            if is_valid_fis(g, fis_tmp) {
                break;
            }
            fis_opt = g.get_fis();
        }
        fis = match fis_opt {
            Some(f) => f,
            None => find_fis_global(g)
        };
    }
    if fis.sg_type == FisType::NONE {       // no FIS found --> g is split cluster graph
        return Some(Vec::new());
    } else if k == 0 {
        return None;                        // FIS found --> cannot be destroyed with k=0 edits!
    }

    branch_on_fis(g, &fis, k, stats, config, ub)
}

/// loops over values for k and checks if graph has a SCE with size k
/// if Some solution is found, the minimum k is found and the solution is returned
fn solve(g: &mut impl Graph, stats: &mut SolverStats, config: &SolverConfig, lb: usize, ub: usize) -> Vec<(usize, usize)> {
    let s: Vec<(usize, usize)> = Vec::new();
    let mut k: isize = lb as isize;

    while k <= ub as isize {
        assert!(k as usize <= stats.splittance_ub);
        println!("k:{k}");
        if k == ub as isize {
            println!("Upper Bound is OPT!");
            let sol = if let Some(c) = &stats.heur_clusters { // ub given by heuristic
                heuristic::compute_edits(g, c, ub)
            } else {  // ub given by splittance (as this is the only other ub possible)
                compute_solution_from_splittance(g, None)
            };
            for (u, v) in &sol {  // apply edits to make g a split cluster graph
                g.flip_edge(*u, *v);
            }
            return sol;
        }
        let sol: Option<Vec<(usize, usize)>> = branch(g, k, stats, config, ub);
        if let Some(solution) = sol {
            return solution;
        }
        k += 1;
    }
    println!("Graph has no solution up to k={k}?!");  // should be unreachable
    s
}

/// alternative branching strategy that walks trough the search tree until the optimal solution is
/// found -> cuts search tree when upper bound (s.ub.len()) is reached as OPT cannot be found in
/// that branch and thus avoids that the search tree becomes too big
fn branch2(g: &mut impl Graph, k: usize, stats: &mut SolverStats, config: &SolverConfig, s: &mut Solutions, level: usize) {
    stats.rec_steps += 1;

    if stats.start_time.elapsed().as_secs_f32() > stats.timeout {
        println!("TIMEOUT!");
        exit(g, stats, config);
    }

    if k >= s.ub.len() {
        return;  // cut branch as better (or equally good) solution is known
    }
    if stats.rec_steps % config.apply_mc_lb as usize == 0 {
        let split = splittance_global(g);
        if split == 0 {  // g is a split cluster graph
            if s.curr_sol.len() < s.ub.len() {
                s.ub = s.curr_sol.clone();  // better solution found -> update ub
            }
            return;
        }
        let mc_lb = mc::min_cut_lower_bound2_global(g, Some(split), Some(s.ub.len() - k));
        if k + mc_lb >= s.ub.len() {
            return;  // cut branch due to mc_lb!
        }
    }

    let fis: FIS;
    if stats.rec_steps % config.apply_packing_lb as usize == 0 {
        let packing_lb = update_packing(g, s.ub.len() - k);
        if packing_lb == 0 {  // no FIS found -> g is a split cluster graph
            if s.curr_sol.len() < s.ub.len() {
                s.ub = s.curr_sol.clone();  // better solution found -> update ub
            }
            return;
        }
        if k + packing_lb >= s.ub.len() {
            return;  // cut branch as it cannot lead to opt. solution
        }
        fis = g.get_fis().unwrap();
    } else {
        // check if a valid FIS is contained in g.fis_packing...
        let mut fis_opt: Option<FIS> = g.get_fis();
        while fis_opt.is_some() {
            let fis_tmp = fis_opt.as_ref().unwrap();
            if is_valid_fis(g, fis_tmp) {
                break;
            }
            fis_opt = g.get_fis();
        }
        fis = match fis_opt {
            Some(f) => f,
            None => find_fis_global(g)
        };
        if fis.sg_type == FisType::NONE {       // no FIS found --> g is split cluster graph
            if s.curr_sol.len() < s.ub.len() {
                s.ub = s.curr_sol.clone();      // current solution is better than UB -> update UB
            }
            return;
        }
    }

    debug_assert_ne!(fis.sg_type, FisType::NONE); // split > 0, so there has to be a FIS!!
    let edits = match fis.sg_type {
        FisType::C4         => branch_on_c4(&fis),
        FisType::C5         => branch_on_c5(&fis),
        FisType::P5         => branch_on_p5(&fis),
        FisType::BOWTIE     => branch_on_bt(&fis),
        FisType::NECKTIE    => branch_on_nt(&fis),
        _                   => panic!(),
    };

    let mut edited_edges: Vec<(usize, usize)> = Vec::with_capacity(10);
    for (br, edges) in edits.iter().enumerate() {
        let apply_dr = stats.rec_steps % config.apply_dr as usize == 0;
        if dr::uneditable(g, edges, s.ub.len() - k, apply_dr) {
            if br >= 6 && br < 9 && (fis.sg_type == FisType::NECKTIE || fis.sg_type == FisType::C5
                || fis.sg_type == FisType::P5) {
                // although edits were skipped, some edges could still be marked as edited
                let mut edited = reedit_edges(g, &fis, br, &mut edited_edges);
                edited_edges.append(&mut edited);
            }
            continue;  // do not branch on uneditable edges
        }

        let mut dr_stack = Vec::new();
        for edge in edges {
            let edge = *edge;
            g.flip_edge(edge.0, edge.1);
            g.edit_edge(edge.0, edge.1);
            s.curr_sol.push(edge);
            if !config.no_dr {
                dr_stack.append(&mut dr::dr_during_branching(g, &edge, apply_dr));
            }
        }
        let k_new = k + edges.len();

        branch2(g, k_new, stats, config, s, level + 1);

        let edges_len = edges.len();
        for edge in edges {
            if edges_len > 1 {
                g.unedit_edge(edge.0, edge.1);
            } else {
                edited_edges.push((edge.0, edge.1));  // unedit only after sibling branches are done
            }
            g.flip_edge(edge.0, edge.1);
            s.curr_sol.pop();
        }
        if edges.len() > 1 {
            // due to the improved branching some edges can be set to edited in certain branches
            let mut edited = reedit_edges(g, &fis, br, &mut edited_edges);
            edited_edges.append(&mut edited);
        }

        dr::undo_dr_during_branching(g, dr_stack);
    }

    for edge in edited_edges {
        g.unedit_edge(edge.0, edge.1);
    }
}

/// computes an optimal solution k for SCE by first computing an upper bound and then recursively
/// branching on forbidden induced subgraphs (lower bounds are used to prune the search tree)
/// -> the search tree has size O(6.2^UB), so a good upper bound is essential!
fn solve2(g: &mut impl Graph, stats: &mut SolverStats, config: &SolverConfig, ub: usize) -> Vec<(usize, usize)> {
    println!("start solve2...");
    let ub_edits = if let Some(c) = &stats.heur_clusters { // ub given by heuristic
        heuristic::compute_edits(g, c, ub)
    } else {  // ub given by splittance (as this is the only other ub possible)
        compute_solution_from_splittance(g, None)
    };
    let mut s = Solutions { curr_sol: vec![], ub: ub_edits };
    branch2(g, 0, stats, config, &mut s, 0);
    for (u, v) in &s.ub {  // apply edits to make g a split cluster graph
        g.flip_edge(*u, *v);
    }
    s.ub
}

/// computes and returns lower and upper bounds for connected graph 'g' as (LB, UB, OPT/None)
/// if UB == 0 or LB == UB, then the optimal solution is returned as well (else: None)
/// CAUTION: assumes that 'g' is connected
fn lower_and_upper_bound(g: &mut impl Graph, stats: &mut SolverStats) -> (usize, usize, Option<Vec<(usize, usize)>>) {
    let split = splittance_sg(g, None);
    stats.splittance_ub += split;
    if split == 0 {                                   // 'g' is a split graph already
        println!("G is a split graph already!");
        stats.kernel_size -= g.size();
        return (0, 0, Some(Vec::new()));
    }
    // 'g' is no split graph!
    let packing_lb = compute_packing(g, None);
    stats.packing_lb += packing_lb;

    let mc_lb = mc::min_cut_lower_bound2(g, Some(split), None, Some(split));
    let lb = max(packing_lb, mc_lb);
    let (heur_ub, clusters) = heuristic::run_heuristic(g, 1.0, 7, 500, lb, split);
    stats.heuristic_ub += heur_ub;
    let ub = min(split, heur_ub);
    debug_assert_ne!(heur_ub, 0);
    debug_assert!(packing_lb <= ub);
    if mc_lb >= split {                                             // k = splittance; mc is no LB
        println!("MinCut >= Splittance -> Splittance gives optimal solution!");
        stats.min_cut_lb = split;
        stats.kernel_size -= g.size();
        return (split, split, Some(compute_solution_from_splittance(g, None)));
    }
    stats.min_cut_lb += mc_lb;                                      // mc_lb < splittance
    if lb == split {                                                // -> lb = packing_lb
        println!("LB = Splittance -> Splittance gives optimal solution!");
        return (split, split, Some(compute_solution_from_splittance(g, None)));
    }
    if lb == heur_ub {
        println!("LB = Heuristic -> Heuristic gives optimal solution!");
        return (lb, lb, Some(heuristic::compute_edits(g, &clusters, heur_ub)));
    }
    assert!(lb < ub);

    if stats.heur_clusters.is_some() {  // clusters from different connected component
        stats.heur_clusters = None;
    }
    if heur_ub == ub {
        stats.heur_clusters = Some(clusters);  // save clusters if heuristic gives upper bound
    } else {
        println!("Heuristic is worse than splittance!");
    }

    (lb, ub, None)  // lb < ub; solution cannot be computed from splittance/ heuristic (yet)
}

/// splits disconnected graph into multiple graphs and solves them independently
/// if 'alt_branching' alternative branching (w/o k as input) is executed
pub fn solve_independently(in_file: String, out_file: String, timeout: f32, g: &mut impl Graph, config: SolverConfig) -> Vec<(usize, usize)> {
    let mut stats = SolverStats::new(in_file, out_file, g.size());
    stats.timeout = timeout;
    let now = Instant::now();                                 // start timer
    stats.start_time = now;

    let mut solution: Vec<(usize,usize)> = Vec::new();

    let mut comps = g.connected_components(None);
    if comps.len() == 1 {                                            // 'g' is connected
        let (lb, ub, sol) = lower_and_upper_bound(g, &mut stats);
        if let Some(s) = sol {
            for (u, v) in &s {  // apply edits to make g a split cluster graph
                g.flip_edge(*u, *v);
            }
            solution = s;
        } else {
            if !config.no_dr {
                let (core, periphery) = dr::deg_1_core_rules(g, Some(ub));
                println!("Core: {:?}\nPeriphery: {:?}", core, periphery);
            }
            stats.labeled_vert = g.labeled_vertices();
            solution = if config.alt_br {
                solve2(g, &mut stats, &config, ub)
            } else {
                solve(g, &mut stats, &config, lb, ub)
            };
        }
    } else {
        comps.sort_by_key(|x| x.len());                   // start w/ smallest comps
        for comp in comps {                                // solve comps independently
            if comp.len() < 4 {                                      // trivially split -> remove
                stats.kernel_size -= comp.len();
                continue;
            }
            let (node_map, mut g_tmp) = g.induced_subgraph(&comp);
            let (lb, ub, sol) = lower_and_upper_bound(&mut g_tmp, &mut stats);
            let solution_sg: Vec<(usize,usize)> =
            if let Some(s) = sol {
                s
            } else {
                if !config.no_dr {
                    let (core, periphery) = dr::deg_1_core_rules(&mut g_tmp, Some(ub));
                    println!("Core: {:?}\nPeriphery: {:?}", core, periphery);
                }
                stats.labeled_vert += g_tmp.labeled_vertices();
                if config.alt_br {
                    solve2(&mut g_tmp, &mut stats, &config, ub)
                } else {
                    solve(&mut g_tmp, &mut stats, &config, lb, ub)
                }
            };
            if !solution_sg.is_empty() {
                for (u, v) in solution_sg {
                    solution.push((node_map[&u], node_map[&v]));
                    g.flip_edge(node_map[&u], node_map[&v]);  // apply edits to make g a split cluster graph
                }
            }
        }
    }
    assert_eq!(splittance_global(g), 0);
    stats.clusters = g.connected_components(None).len();
    stats.runtime = now.elapsed().as_secs_f32();
    stats.k = Some(solution.len());
    stats.write_stats(&config).unwrap_or_else(|err| {
        println!("Problem writing statistics to file: {err}");
    });

    solution
}

/// exits the program, but saves statistics before
fn exit(_g: &mut impl Graph, stats: &mut SolverStats, config: &SolverConfig) {
    stats.runtime = stats.start_time.elapsed().as_secs_f32();
    stats.k = None; // no solution found within time limit!
    stats.write_stats(config).unwrap_or_else(|err| {
        println!("Problem writing statistics to file: {err}");
    });
    process::exit(0);
}