import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.ticker import LogLocator, NullFormatter
import tikzplotlib
import numpy as np
import networkx as nx

out_path = ''  # path to the folder where the plots' .pdf files are saved

plt.rcParams["font.family"] = "serif"

label = {'random_ER': 'ER', 'random_SBM': 'SBM', 'bio_33_50': 'Bio33', 'bio_50_50': 'Bio50',
         'bio_66_50': 'Bio66', 'bio33': 'Bio33', 'bio50': 'Bio50', 'bio66': 'Bio66'}

# all possible markers: . , o v ^ < > 1 2 3 4 8 s p P * h H + x X D d | _
marker = {'random_ER': 'p', 'random_SBM': '^', 'bio_33_50': 'o', 'bio_50_50': 'd', 'bio_66_50': 's', 'bio33': 'o',
          'bio50': 'd', 'bio66': 's'}

size = {'random_ER': 15, 'random_SBM': 13, 'bio_33_50': 13, 'bio_50_50': 13, 'bio_66_50': 13, 'bio33': 13, 'bio50': 13,
        'bio66': 13}

# colorblind-safe palette:
# [blue, bright orange, light gray, dark gray, light blue, dark orange, gray, very light blue, light orange, light gray]
c = ['#006ba4', '#ff800e', '#ababab', '#595959', '#5f9ed1', '#c85200', '#898989', '#a2c8ec', '#ffbc79', '#cfcfcf']
color = {'random_ER': c[5], 'random_SBM': c[2], 'bio_33_50': c[4], 'bio_50_50': c[3], 'bio_66_50': c[8],
         'bio33': c[4], 'bio50': c[3], 'bio66': c[8]}

marker_random = {0: 'P', 1: 'o', 2: 'p', 3: 's', 4: 'v', 5: '*', 6: 'X', 7: '^', 8: 'd'}
color_random = {0: c[0], 1: c[1], 2: c[2], 3: c[3], 4: c[4], 5: c[5], 6: c[6], 7: c[7], 8: c[8]}
size_random = {0: 15, 1: 13, 2: 17, 3: 13, 4: 13, 5: 30, 6: 13, 7: 13, 8: 13}

my_cmap = mpl.cm.viridis

# min. and max. graph sizes for the instance types:
min_n_bio, max_n_bio = 6, 50
min_n_random, max_n_random = 10, 30
max_n_default = 50


def jitter_dots(dots):
    # code found at: https://stackoverflow.com/questions/60229586/how-to-make-jitterplot-on-matplolib-python
    offsets = dots.get_offsets()
    jittered_offsets = offsets
    # only jitter in the x-direction
    jittered_offsets[:, 0] += np.random.uniform(-0.2, 0.2, offsets.shape[0])
    jittered_offsets[:, 1] += np.random.uniform(-0.2, 0.2, offsets.shape[0])
    dots.set_offsets(jittered_offsets)


def split_file(file, ilp=False):
    """
    splits 'file' into five files, each containing one group of instances
    :param file: statistics file generated during the experiments (contains statistics of all graphs
    :param ilp: if the file contains the statistics of the ILP solver (-> timed-out instances are not listed)
    :return: None (creates five files)
    """
    file_endings = ['_bio33', '_bio50', '_bio66', '_er', '_sbm']
    no_instances = [90, 90, 90, 198, 198]
    with open(file=file, mode='r', encoding='utf-8') as f:
        all_lines = f.readlines()
        assert '\t' in all_lines[0]
        line0 = all_lines[0]
        lines = [line0]
        last_path = '33_50'
        i = 0
        for j in range(1, len(all_lines)):
            lines.append(all_lines[j])
            if j == len(all_lines) - 1 or all_lines[j + 1].strip().split('\t')[0].split('/')[3] != last_path:
                out_file = file[:-4] + file_endings[i] + '.csv'
                if ilp:
                    if len(lines) - 1 != no_instances[i]:
                        lines = fill_missing_lines_ilp(lines)
                assert len(lines) - 1 == no_instances[i]
                with open(file=out_file, mode='w', encoding='utf-8') as f_out:
                    for line in lines:
                        f_out.write(line)
                print("written to file {}".format(out_file))
                lines = [line0]
                i += 1
                if j != len(all_lines) - 1:
                    last_path = all_lines[j + 1].strip().split('\t')[0].split('/')[3]


def fill_missing_lines_ilp(lines):
    """
    adds lines missing in ILP stats (ILP solver doesn't have timeout option -> doesn't write unsolved instances to file)
    :param lines: statistics of solved instances
    :return: new_lines -> all statistics (solved and unsolved)
    """
    new_lines = [lines[0]]
    file = lines[1].split('\t')[0].split('/')
    dir_path = file[0] + '/' + file[1] + '/' + file[2] + '/' + file[3]
    i = 1
    for filename in sorted(os.listdir(dir_path)):
        if 'solution' not in filename:  # skip solution files
            if filename in lines[i]:  # instance was solved
                new_lines.append(lines[i])
                i += 1
            else:  # instances timed-out -> not listed -> add manually
                print("{} missing -> timeout".format(filename))
                n = int(filename.split('_')[1])  # graph size
                constr = int(((n * (n - 1) * (n - 2)) + (n * (n - 1))) / 2)  # total number of constraints
                line = dir_path + '/' + filename + '\t3600.0\t{}\t{}\tunsolved\t{}\t{}\n'.format(n, n, constr, constr)
                new_lines.append(line)

    return new_lines


def extract_values(file, values):
    """
    extracts all 'values' in given 'file' and returns them in an array of arrays (in the same order)
    :param file: input file
    :param values: array of values to extract from it (e.g. 'Graph Size', 'k', etc.)
    :return: array of arrays of requested values found in 'file'
    """
    assert file[-4:] == '.csv'
    output_arr = []
    with open(file=file, mode='r', encoding='utf-8') as f:
        lines = f.readlines()
        assert '\t' in lines[0]
        line0 = lines[0].strip().split('\t')
        for value in values:
            assert value in ['Time (s)', 'Graph Size', 'Kernel Size', 'k', 'Clusters', 'Recursive Steps', 'MinCut (LB)',
                             'Packing Size (LB)', 'File', 'Splittance (UB)', 'Heuristic (UB)', 'Labeled Vertices',
                             'Constraints (needed)']
            for i in range(0, len(line0)):  # find index of 'value' in first line
                if line0[i] == value:
                    break
            assert line0[i] == value  # make sure index was actually found
            output = []
            for j in range(1, len(lines)):  # collect all values in this column in an array
                line = lines[j].strip().split('\t')
                if value in ['Graph Size', 'Kernel Size', 'Recursive Steps', 'MinCut (LB)', 'Constraints (needed)',
                             'Clusters', 'Splittance (UB)', 'Heuristic (UB)', 'Labeled Vertices', 'Packing Size (LB)']:
                    output.append(int(line[i]))  # convert to int
                elif value == 'Time (s)':
                    output.append(float(line[i]))  # convert time to float
                elif value == 'k':
                    if line[i] == 'unsolved':
                        output.append(-1)  # store unsolved instances as k = -1
                    else:
                        output.append(int(line[i]))
                else:
                    output.append(line[i])  # don't convert file name
            output_arr.append(output)

    # make sure all arrays have the same size
    for i in range(len(output_arr) - 1):
        assert len(output_arr[i]) == len(output_arr[i + 1])

    return output_arr


def split_solved_unsolved(k, values):
    """
    splits 'values' into solved and unsolved instances based on the values in 'k'
    :param k: array of solution sizes (-1, if the instance was not solved)
    :param values: array of arrays containing other values collected for the instances
    :return:
    """
    solved_k, solved_arr = [], [[] for _ in range(len(values))]
    unsolved_k, unsolved_arr = [], [[] for _ in range(len(values))]
    for i in range(len(k)):
        if k[i] != -1:
            solved_k.append(k[i])
            for j in range(len(values)):
                solved_arr[j].append(values[j][i])
        else:
            for j in range(len(values)):
                unsolved_arr[j].append(values[j][i])

    return solved_k, solved_arr, unsolved_arr


def delta_ub_lb(spl, heur, pack, mc):
    """
    computes difference between upper and lower bound
    :param spl: Splittance UB
    :param heur: Heuristic UB
    :param pack: Packing LB
    :param mc: MinCut LB
    :return: delta = min(s,h) - max(p,m)
    """
    # compute UB and LB:
    assert len(spl) == len(heur)
    ub = [min(int(spl[i]), int(heur[i])) for i in range(len(spl))]
    assert len(pack) == len(mc)
    lb = [max(int(pack[i]), int(mc[i])) for i in range(len(pack))]
    assert len(ub) == len(lb)
    delta = [ub[i] - lb[i] for i in range(len(ub))]
    assert min(delta) >= 0  # LB cannot be larger than UB!
    if ub != list(map(int, heur)):
        print("Heuristic is not always better than Splittance!!")

    return delta


def optimality_of_upper_bounds(files, inst):
    """
    prints statistics about the upper bounds (heuristic and splittance): how often they have been optimal, etc.
    :param files: input files from which the statistics shall be read
    :param inst: names of the instances ([bio33, bio50, ...])
    :return: None
    """
    cntr = 0
    for i in range(len(files)):
        print("--------------- {} instances ---------------".format(inst[i]))
        f, k, s, h = extract_values(files[i], ['File', 'k', 'Splittance (UB)', 'Heuristic (UB)'])
        heur_opt, split_opt, solv = 0, 0, 0         # stats for solved instances
        heur_opt_us, split_opt_us, us = 0, 0, 0     # stats for unsolved instances
        heur_off_by, heur_off_by_us = [], []        # delta(heur,OPT)
        for j in range(len(f)):
            if k[j] == -1:                          # unsolved instance
                unsolved_inst = f[j] + '.solution'
                try:
                    with open(file=unsolved_inst, mode='r', encoding='utf-8') as f_us:
                        opt = len(f_us.readlines())
                        us += 1
                        if opt == h[j]:                 # heuristic is optimal
                            heur_opt_us += 1
                        else:
                            heur_off_by_us.append(h[j]-opt)
                        if opt == s[j]:                 # splittance is optimal
                            split_opt_us += 1
                except FileNotFoundError:
                    print("no .solution file for", f[j])
            else:                                   # solved instance
                solv += 1
                if k[j] == h[j]:                        # heuristic is optimal
                    heur_opt += 1
                else:
                    heur_off_by.append(h[j]-k[j])
                if k[j] == s[j]:                        # splittance is optimal
                    split_opt += 1
            assert h[j] <= s[j]
            cntr += 1
        print("STATISTICS FOR ALL SOLVED INSTANCES:")
        print("\tThe Heuristic is optimal for {} instances ({:.2f}%)".format(heur_opt, heur_opt / solv * 100))
        print("\tWhen the Heuristic is not optimal, it is off by: ", heur_off_by)
        print("\tThe Splittance is optimal for {} instances ({:.2f}%)".format(split_opt, split_opt / solv * 100))
        print("STATISTICS FOR UNSOLVED INSTANCES FOR WHICH OPT IS KNOWN:")
        print("\tThe Heuristic is optimal for {} further instances ({:.2f}%)".format(heur_opt_us, heur_opt_us / us * 100))
        print("\tWhen the Heuristic is not optimal, it is off by: ", heur_off_by_us)
        print("\tThe Splittance is optimal for {} further instances ({:.2f}%)".format(split_opt_us, split_opt_us / us * 100))

    assert cntr == 666


def efficiency_of_initial_data_reduction(files_bb, files_ilp, inst):
    """
    prints statistics about the efficiency of the initial data reduction
    :param files_bb: input files  from which the statistics for the BB solver shall be read
    :param files_ilp: input files  from which the statistics for the ILP solver shall be read
    :param inst: names of the instances ([bio33, bio50, ...])
    :return: None
    """
    print("\nEFFICIENCY OF INITIAL DATA REDUCTION:")
    for i in range(len(files_bb)):
        print("--------------- {} instances ---------------".format(inst[i]))
        f, n, ks_bb, lv = extract_values(files_bb[i], ['File', 'Graph Size', 'Kernel Size', 'Labeled Vertices'])
        ks_ilp = extract_values(files_ilp[i], ['Kernel Size'])[0]
        assert min([ks_ilp[x] - ks_bb[x] for x in range(len(ks_bb))]) >= 0
        no_dr = 0
        labeled_vertices = []
        kernel_size = []
        for j in range(len(f)):
            if lv[j] == 0:
                no_dr += 1
            labeled_vertices.append(lv[j] / n[j])
            kernel_size.append(ks_ilp[j] / n[j])
        print("Initial Data Reduction cannot label any vertex for {} of {} instances ({:.2f}%)"
              .format(no_dr, len(f), no_dr/len(f)*100))
        print("Initial Data Reduction labels {:.2f}% of a graph's vertices on average (median={:.2f}, max={:.2f})"
              .format(np.average(labeled_vertices)*100, np.median(labeled_vertices)*100, max(labeled_vertices)*100))
        print("Kernel size is on average {:.2f}% of the initial size (median={:.2f}, min={:.2f})"
              .format(np.average(kernel_size)*100, np.median(kernel_size)*100, min(kernel_size)*100))


def number_of_clusters_2plots(files1, files2, inst1, inst2):
    """
    plots difference between number of connected components in the input graph and clusters in the final split cluster graph
    :param files1: input data files for first (left) plot
    :param files2: input data files for second (right) plot
    :param inst1: list of instance names for first plot (to retrieve markers and labels)
    :param inst2: list of instance names for second plot (to retrieve markers and labels)
    :return: None (shows plot and saves it in a file)
    """
    print("plotting number of clusters (before and after)...")
    plt.figure(figsize=(10, 5))
    plt.suptitle('Number of Connected Components and Clusters')
    my_norm = colors.Normalize(vmin=0, vmax=14)

    files = [files1, files2]
    inst = [inst1, inst2]

    for j in range(2):
        plt.subplot(1, 2, j + 1)
        if j == 0:
            plt.ylabel("#Clusters - #Components")

        # add light gray grid in the background:
        plt.grid(color='0.95', which='minor', zorder=1)
        plt.grid(color='0.95', which='major', zorder=1)

        for i in range(len(files[j])):
            print("--------------- {} instances ---------------".format(inst[j][i]))
            # extract number of clusters from the file
            clusters, instance = extract_values(files[j][i], ['Clusters', 'File'])
            clusters_solved, conn_comps = [], []
            for idx in range(len(instance)):
                if clusters[idx] > 0:  # only consider solved instances
                    clusters_solved.append(clusters[idx])
                    # generate graph and count its connected components:
                    g = nx.Graph()
                    with open(file=instance[idx], mode='r', encoding='utf-8') as f:
                        lines = f.readlines()
                        assert lines[0][0] == '#'
                        n = (int(lines[0].strip().split(' ')[1]))
                        g.add_nodes_from(range(1, n+1))
                        for edges in range(1, len(lines)):
                            if lines[edges][0] != '#' and len(lines[edges]) > 0:
                                edge = lines[edges].strip().split(' ')
                                g.add_edge(int(edge[0]), int(edge[1]))
                    assert g.number_of_nodes() == n
                    conn_comps.append(len(list(nx.connected_components(g))))
            print("max. num of conn comps: {}".format(max(conn_comps)))
            delta = [clusters_solved[i] - conn_comps[i] for i in range(len(clusters_solved))]
            assert min(delta) >= 0  # an optimal solution never connects different connected components

            # sort by (increasing) difference to #clusters:
            conn_comps = [val for (_, val) in sorted(zip(delta, conn_comps), key=lambda z: z[0])]
            delta.sort()
            # create colors from number of connected components
            colors_solved = my_cmap(my_norm(conn_comps))

            # add small offset to reduce overlapping:
            delta = [delta[x]+0.1*i for x in range(len(delta))]
            # plot:
            plt.plot(range(1, len(delta) + 1), delta, c='black', linewidth=0.3, zorder=2)
            plt.scatter([x for x in range(1, len(delta) + 1)], delta, edgecolors='None', c=colors_solved,
                        marker=marker[inst[j][i]], s=size[inst[j][i]], linewidth=.5, label=label[inst[j][i]], zorder=2)
        plt.xlabel("Instance Count")
        plt.yticks([1, 2, 3, 4, 5, 6])
        plt.minorticks_on()
        plt.tick_params(which='both', direction='in')
        plt.legend(loc='upper left')
        plt.tight_layout()
    c_bar = plt.colorbar(mpl.cm.ScalarMappable(norm=my_norm, cmap=my_cmap), ticks=range(1, 14, 2))
    c_bar.set_label('#comp.', labelpad=-26, y=-0.025, rotation=0)
    plt.tight_layout()
    plt.savefig(out_path + 'clusters.pdf')
    # tikzplotlib.save(out_path + 'clusters.tikz')  # -> needs some modifications before it compiles!
    plt.show()


def packing_against_min_cut_2plots(files1, files2, inst1, inst2):
    """
       plots MinCut Lower Bound (x-axis) against Packing Lower Bound (y-axis) and indicates the graph size by the color
       -> plots data from 'files1' and 'files2' side by side
       :param files1: input data files for first (left) plot
       :param files2: input data files for second (right) plot
       :param inst1: list of instance names for first plot (to retrieve markers and labels)
       :param inst2: list of instance names for second plot (to retrieve markers and labels)
       :return: None (shows plot and saves it in a file)
    """
    print("plotting lower bound comparison...")
    plt.figure(figsize=(10, 5))
    plt.suptitle('Size of Packing and MinCut Lower Bound')

    max_val_p, max_val_m = 0, 0
    min_n, max_n = 6, max_n_bio

    my_norm = colors.Normalize(vmin=min_n, vmax=max_n)
    files = [files1, files2]
    inst = [inst1, inst2]

    for j in range(2):
        plt.subplot(1, 2, j+1)
        if j == 0:
            plt.ylabel("Packing Lower Bound")

        # add light gray grid in the background:
        plt.grid(color='0.95', which='minor', zorder=1)
        plt.grid(color='0.95', which='major', zorder=1)

        for i in range(len(files[j])):
            print("--------------- {} instances ---------------".format(inst[j][i]))
            # extract lower bounds from the file
            packing_lb, mincut_lb, n, k = extract_values(files[j][i],
                                                         ['Packing Size (LB)', 'MinCut (LB)', 'Graph Size', 'k'])

            # split into solved and unsolved instances:
            _, [m_solved, p_solved, n_solved], [m_unsolved, p_unsolved, n_unsolved] = \
                split_solved_unsolved(k, [mincut_lb, packing_lb, n])

            min_n = min(min_n, min(n))
            max_val_p = max(max_val_p, max(packing_lb))
            max_val_m = max(max_val_m, max(mincut_lb))
            diff_lbs = [mincut_lb[x] - packing_lb[x] for x in range(len(mincut_lb))]
            print("MinCut LB performed better for {} instance(s)".format(len([x for x in diff_lbs if x > 0])))
            print("LBs performed equally good for {} instance(s)".format(len([x for x in diff_lbs if x == 0])))
            print("Packing LB performed better for {} instance(s)".format(len([x for x in diff_lbs if x < 0])))

            # create colors from graph sizes:
            colors_solved = my_cmap(my_norm(n_solved))
            colors_unsolved = my_cmap(my_norm(n_unsolved))

            # plot solved (filled markers) and unsolved (unfilled markers) instances
            # and jitter slightly to reduce overlapping:
            d1 = plt.scatter(m_unsolved, p_unsolved, edgecolors=colors_unsolved, facecolors='None',
                             marker=marker[inst[j][i]], s=size[inst[j][i]], linewidth=.7, zorder=2, alpha=.85)
            jitter_dots(d1)
            d2 = plt.scatter(m_solved, p_solved, edgecolors=colors_solved, facecolors=colors_solved,
                             marker=marker[inst[j][i]], s=size[inst[j][i]], linewidth=.5, label=label[inst[j][i]],
                             zorder=2, alpha=.85)
            jitter_dots(d2)
        min_max_val = min(max_val_m, max_val_p)
        plt.plot([-2, min_max_val + 5], [-2, min_max_val + 5], ls='dashed', c='gray', alpha=0.5)
        plt.xlim(-1, max_val_m + 1)
        plt.ylim(-1, max_val_p + 1)
        plt.xlabel("MinCut Lower Bound")
        plt.tick_params(which='both', direction='in')
        plt.minorticks_on()
        plt.legend()
        plt.tight_layout()

    c_bar = plt.colorbar(mpl.cm.ScalarMappable(norm=my_norm, cmap=my_cmap), ticks=range(10, max_n + 1, 5))
    c_bar.set_label('n', labelpad=-26, y=-0.025, rotation=0)

    plt.tight_layout()
    plt.savefig(out_path + 'packing_mc_lb.pdf')
    # tikzplotlib.save(out_path + 'packing_mc_lb.tikz')  # -> needs some modifications before it compiles!
    plt.show()


def plot_runtimes_2plots(files_dr, files_no_dr, inst):
    """
    plots running times of two solvers w/ and w/o data reduction applied
    :param files_dr: 2x [bio33,bio50,bio66,ER,SMB] data w/ DR applied (files_dr[0]: left plot, files_dr[1]: right plot)
    :param files_no_dr: same as files_dr but w/o DR
    :param inst: list of instance names (to retrieve markers and labels)
    :return: None (shows plot and saves it in a file)
    """
    plt.figure(figsize=(10, 5))
    plt.suptitle('Running Times with and without Data Reduction (DR)')
    all_ilp, all_bb = [], []
    for j in range(2):
        plt.subplot(1, 2, j+1)
        if j == 0:
            print("plotting running times for ILP solver w/ and w/o data reduction...")
            plt.title('ILP Solver')
            plt.ylabel('Ratio (Without DR / With DR)')
        else:
            print("plotting running times for BB solver w/ and w/o data reduction...")
            plt.title('BB Solver')
        plt.xlabel('With DR [s]')

        # plot horizontal lines marking running time factors:
        plt.plot([0, 3800], [1, 1], ls='solid', linewidth=1, c=c[9], zorder=1)
        plt.plot([0, 3800], [2, 2], ls=(0, (5, 5)), c=c[9], linewidth=1, zorder=1)      # dashed
        plt.plot([0, 3800], [5, 5], ls=(0, (1, 5)), c=c[9], linewidth=1, zorder=1)      # dotted
        plt.plot([0, 3800], [.5, .5], ls=(0, (5, 5)), c=c[9], linewidth=1, zorder=1)
        plt.plot([0, 3800], [.2, .2], ls=(0, (1, 5)), c=c[9], linewidth=1, zorder=1)

        # files_{no_}dr[0]: [bio33,bio50,bio66,ER,SMB] for ILP solver
        # files_{no_}dr[1]: ----------- ''------------ for BB solver
        max_ratio = 0
        all_r_dr = []
        for i in range(len(files_dr[j])):
            print("---------------- {} instances ----------------".format(inst[i]))
            r_dr, r_no_dr, solved_w_dr, solved_wo_dr = [], [], [], []
            if j == 0:  # ILP solver:
                instance, r_dr, c_dr = extract_values(files_dr[j][i], ['File', 'Time (s)', 'Constraints (needed)'])
                r_no_dr, c_no_dr = extract_values(files_no_dr[j][i], ['Time (s)', 'Constraints (needed)'])
                all_r_dr += r_dr
                faster_wo = []
                for idx in range(len(r_dr)):
                    if r_dr[idx] >= 3600.0 and r_no_dr[idx] < 3600.0:
                        solved_wo_dr.append(r_no_dr[idx])
                        print("Instance {} solved w/o DR, but unsolved w/ DR!".format(instance[idx]))
                    if r_dr[idx] < 3600.0 and r_no_dr[idx] >= 3600.0:
                        solved_w_dr.append(r_dr[idx])
                        print("Instance {} solved w/ DR, but unsolved w/o DR!".format(instance[idx]))
                    if r_no_dr[idx] < r_dr[idx]:
                        faster_wo.append(r_dr[idx] - r_no_dr[idx])
                if len(faster_wo) == 0:
                    print("All instances solved faster w/ DR")
                else:
                    print("{} instances solved faster w/o DR (max. time difference = {}s, avg. = {}s"
                          .format(len(faster_wo), max(faster_wo), np.average(faster_wo)))
                no_reduction = [abs(r_no_dr[x]-r_dr[x]) for x in range(len(r_dr)) if c_dr[x] == c_no_dr[x]]
                print("Max. running time difference for an instance to which DR was not applicable: {}s"
                      .format(max(min(no_reduction), max(no_reduction))))
                assert len([c_dr[x] for x in range(len(c_dr)) if c_dr[x] > c_no_dr[x]]) == 0  # DR only decreases #con.
                print("{} instances ({:.2f}%) solved".format(
                    len([r_dr[x] for x in range(len(r_dr)) if r_dr[x] < 3600.0]),
                    len([r_dr[x] for x in range(len(r_dr)) if r_dr[x] < 3600.0]) / len(instance) * 100),
                    )
                # exclude data points for instances to which initial DR was not applicable:
                # -> still plot all data points to illustrate running time of the solver!!
                r_dr2 = [r_dr[x] for x in range(len(r_dr)) if c_dr[x] < c_no_dr[x]]
                #r_no_dr2 = [r_no_dr[x] for x in range(len(r_no_dr)) if c_dr[x] < c_no_dr[x]]
                print("Initial DR was possible for {:.2f}% of the instances".format((len(r_dr2) / len(c_dr)) * 100))
                c_dr = [c_dr[x]/c_no_dr[x] for x in range(len(c_dr)) if c_dr[x] < c_no_dr[x]]
                print("After data reduction only {:.2f}% of constraints were needed on average (compared to no "
                      "data reduction)".format(np.average(c_dr) * 100))
            else:
                instance, r_dr, k = extract_values(files_dr[j][i], ['File', 'Time (s)', 'k'])
                r_no_dr, k = extract_values(files_no_dr[j][i], ['Time (s)', 'k'])
                all_r_dr += r_dr
                faster_wo = []
                for idx in range(len(r_dr)):
                    if r_dr[idx] >= 3600.0 and r_no_dr[idx] < 3600.0:
                        solved_wo_dr.append(r_no_dr[idx])
                        print("Instance {} solved w/o DR, but unsolved w/ DR!".format(instance[idx]))
                    if r_dr[idx] < 3600.0 and r_no_dr[idx] >= 3600.0:
                        solved_w_dr.append(r_dr[idx])
                        print("Instance {} solved w/ DR, but unsolved w/o DR!".format(instance[idx]))
                    if r_no_dr[idx] < r_dr[idx]:
                        faster_wo.append(r_dr[idx] - r_no_dr[idx])
                if len(faster_wo) == 0:
                    print("All instances solved faster w/ DR")
                print("{} instances solved faster w/o DR (max. time difference = {}s, avg. = {}s"
                      .format(len(faster_wo), max(faster_wo), np.average(faster_wo)))
            # add small offset to avoid runtime of 0 and exclude unsolved instances:
            r_dr_solved = [r_dr[x]+0.01 for x in range(len(r_dr)) if r_dr[x] < 3600.0 and r_no_dr[x] < 3600.0]
            r_no_dr_solved = [r_no_dr[x]+0.01 for x in range(len(r_no_dr)) if r_dr[x] < 3600.0 and r_no_dr[x] < 3600.0]
            assert len(r_dr_solved) == len(r_no_dr_solved)
            if j == 1:
                print("{} instances ({:.2f}%) solved".format(len(r_dr_solved)+len(solved_w_dr),
                                                             len(r_dr_solved) / len(r_dr)*100))
            ratio = [r_no_dr_solved[x] / r_dr_solved[x] for x in range(len(r_dr_solved))]
            max_ratio = max(max_ratio, max(ratio))

            plt.scatter(r_dr_solved, ratio, edgecolors='None', marker=marker[inst[i]], label=label[inst[i]], c=color[inst[i]],
                        s=size[inst[i]] + 3, zorder=2, alpha=.75)
            if len(solved_w_dr) != 0:
                plt.scatter(solved_w_dr, [6]*len(solved_w_dr), edgecolors='None', marker=marker[inst[i]], c='black',
                            s=size[inst[i]] + 3, zorder=2)
            if len(solved_wo_dr) != 0:
                plt.scatter(solved_wo_dr, [0.5]*len(solved_wo_dr), edgecolors='black', facecolors='None',
                            marker=marker[inst[i]], c='black', s=size[inst[i]] + 3, zorder=2)

        assert len(all_r_dr) == 666
        if len(all_ilp) == 0:
            all_ilp = all_r_dr
        else:
            all_bb = all_r_dr
        p50 = sorted(all_r_dr)[round(len(all_r_dr) * .5) - 1]
        print("50% of all instances were solved within {} seconds".format(p50))
        plt.plot([p50, p50], [0, 6.3], ls='dashdot', linewidth=1, c=c[9], zorder=1)
        all_r_dr = [all_r_dr[x] for x in range(666) if all_r_dr[x] < 3600]
        print("\n{} ({:.2f}%) of all instances solved within 3600s.".format(len(all_r_dr), len(all_r_dr) / 666 * 100))
        all_r_dr = [all_r_dr[x] for x in range(len(all_r_dr)) if all_r_dr[x] < 1000]
        print("{:.2f}% of all instances solved within 1000s.".format(len(all_r_dr) / 666 * 100))
        all_r_dr = [all_r_dr[x] for x in range(len(all_r_dr)) if all_r_dr[x] < 100]
        print("{:.2f}% of all instances solved within 100s.".format(len(all_r_dr) / 666 * 100))
        all_r_dr = [all_r_dr[x] for x in range(len(all_r_dr)) if all_r_dr[x] < 10]
        print("{:.2f}% of all instances solved within 10s.".format(len(all_r_dr) / 666 * 100))
        all_r_dr = [all_r_dr[x] for x in range(len(all_r_dr)) if all_r_dr[x] < 1]
        print("{:.2f}% of all instances could be solved within 1s.\n".format(len(all_r_dr) / 666 * 100))
        plt.xscale('log')
        plt.xlim(0.005, 3800)
        plt.ylim(0, 6.3)
        plt.legend(loc='upper left')
        plt.tight_layout()

    ilp_faster = [all_bb[x] - all_ilp[x] for x in range(666) if all_bb[x] - all_ilp[x] > 0.0]
    bb_faster = [all_ilp[x] - all_bb[x] for x in range(666) if all_ilp[x] - all_bb[x] > 0.0]
    print("BB solver is faster for {} instances (avg: {}s, max: {}s)"
          .format(len(bb_faster), np.average(bb_faster), max(bb_faster)))
    print("ILP solver is faster for {} instances (avg: {}s, median: {}s, max: {}s)"
          .format(len(ilp_faster), np.average(ilp_faster), np.median(ilp_faster), max(ilp_faster)))

    plt.tight_layout()
    plt.savefig(out_path + 'runtimes.pdf')
    # tikzplotlib.save(out_path + 'runtimes.tikz')
    plt.show()


def ilp_runtime_against_constraints_2plots(files1, files2, inst1, inst2):
    """
        plots the running times of the ILP solver against the number of constraints (x-axis) and optimal solution size
        :param files1: input data files for first (left) plot
        :param files2: input data files for second (right) plot
        :param inst1: list of instance names for first plot (to retrieve markers and labels)
        :param inst2: list of instance names for second plot (to retrieve markers and labels)
        :return: None (shows plot and saves it in a file)
    """
    print("plotting running times and number of constraints for ILP solver...")
    rt, con, k = [], [], []
    files = [files1, files2]
    inst = [inst1, inst2]
    for j in range(2):
        r1, c1, sol_size = [], [], []
        for i in range(len(files[j])):
            runtimes1, constr1, sol = extract_values(files[j][i], ['Time (s)', 'Constraints (needed)', 'k'])
            r1.append(runtimes1)
            c1.append(constr1)
            sol_size.append(sol)
        rt.append(r1)
        con.append(c1)
        k.append(sol_size)

    plt.figure(figsize=(10, 5))
    plt.suptitle('Running Times of the ILP Solver')
    max_sol_size = 0
    for j in range(2):
        for i in range(len(rt[j])):
            max_sol_size = max(max_sol_size, max(k[j][i]))
    my_norm = colors.Normalize(vmin=0, vmax=max_sol_size)

    for j in range(2):
        plt.subplot(1, 2, j+1)
        plt.tick_params(which='both', direction='in')
        if j == 0:
            plt.ylabel("Running Time [s]")

        plt.grid(color='0.95', zorder=1)  # add light gray grid in the background
        for i in range(len(rt[j])):
            rt[j][i] = [rt[j][i][x] + 0.01 for x in range(len(rt[j][i]))]  # avoid runtime of 0 (invisible on log scale)
            k_solv, [rt_solv, con_solv], [rt_unsolv, con_unsolv] = split_solved_unsolved(k[j][i], [rt[j][i], con[j][i]])
            print("{} unsolved instance(s) for {}".format(len(rt_unsolv), label[inst[j][i]]))
            colors_k = my_cmap(my_norm(k_solv))
            plt.scatter(con_solv, rt_solv, edgecolors='None', marker=marker[inst[j][i]], label=label[inst[j][i]], c=colors_k,
                        s=size[inst[j][i]] + 3, zorder=2, alpha=.75)
            plt.scatter(con_unsolv, rt_unsolv, edgecolors=c[3], facecolors='None', marker=marker[inst[j][i]],
                        s=size[inst[j][i]] + 3,
                        zorder=2, alpha=.75)

        plt.xlabel("Number of Constraints")
        plt.yscale('log')
        plt.legend()
        plt.tight_layout()

    c_bar = plt.colorbar(mpl.cm.ScalarMappable(norm=my_norm, cmap=my_cmap), ticks=range(0, max_sol_size + 1, 10))
    c_bar.set_label('sol. size', labelpad=-26, y=-0.025, rotation=0)

    plt.tight_layout()
    plt.savefig(out_path + 'ilp_runtimes_constr.pdf')
    # tikzplotlib.save(out_path + 'ilp_runtimes_constr.tikz')
    plt.show()


def sol_sizes_bio(files, m):
    """
    plots solution sizes in dependence of the graph size for the biological instances
    :param files: [bio33, bio50, bio66]
    :param m: list of instance names (to retrieve markers and labels)
    """
    plt.grid(color='0.95', zorder=1)  # add light gray grid in the background
    min_graph_size, max_graph_size, max_sol_size = 5, 0, 0

    for i in range(len(files)):
        k, n = extract_values(files[i], ['k', 'Graph Size'])

        solved_k, [solved_n], _ = split_solved_unsolved(k, [n])  # remove unsolved instances

        dots = plt.scatter(solved_n, solved_k, c=color[m[i]], marker=marker[m[i]], label=label[m[i]], s=size[m[i]] + 8,
                           zorder=2, alpha=.75, edgecolors='None')
        jitter_dots(dots)

        max_graph_size = max(max_graph_size, max(n))
        max_sol_size = max(max_sol_size, max(k))

    plt.xlabel("Graph Size (n)")
    plt.legend()
    plt.xticks(range(min_graph_size, max_graph_size + 5, 5))
    plt.yticks(range(0, max_sol_size + 1, 5))
    plt.tick_params(which='both', direction='in')


def solution_sizes_random(file, er):
    """
    plots solution sizes in dependence of the graph size for the random instances
    :param file: graph file
    :param er: if the random graphs are ER; else: SBM
    """
    instance, k, n = extract_values(file, ['File', 'k', 'Graph Size'])

    instance_name = [instance[i].split('/')[-1] for i in range(len(instance))]  # only file name
    if er:
        f = list(map(lambda prob: int(prob[10]), instance_name))  # retrieve p value (in {1..9})
        f_map = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8}
        probs = 9
    else:
        f = list(map(lambda prob: int(prob[9:11]), instance_name))  # retrieve p value (in {80,88,95}) -> q = 1 - p
        f_map = {80: 0, 88: 1, 95: 2}
        probs = 3

    # split according to edge insertion probability used:
    n_arr, k_arr = [[] for _ in range(probs)], [[] for _ in range(probs)]
    for i in range(len(n)):
        n_arr[f_map[f[i]]].append(n[i])
        k_arr[f_map[f[i]]].append(k[i])

    plt.grid(color='0.95', zorder=1)  # add light gray grid in the background

    # retrieve solved instances
    solved_k, solved_n = [[] for _ in range(probs)], [[] for _ in range(probs)]
    for i in range(len(k)):
        if k[i] != -1:
            solved_k[f_map[f[i]]].append(int(k[i]))
            solved_n[f_map[f[i]]].append(n[i])

    # plot
    for i in range(len(solved_n)):
        if er:
            prob_label = 'p=0.{}'.format(i + 1)
        else:
            prob_label = 'p=0.{}, q=0.5, r={:.2f}'.format(list(f_map.keys())[i], (100 - list(f_map.keys())[i]) / 100)

        dots = plt.scatter(solved_n[i], solved_k[i], edgecolors='None', c=color_random[i], marker=marker_random[i],
                           label=prob_label, s=size_random[i] + 8, zorder=2, alpha=.75)
        jitter_dots(dots)

    max_k_solved = max([max(solved_k[i]) for i in range(probs)])
    plt.xlabel("Graph Size (n)")
    plt.legend(fontsize=9)
    plt.xticks(range(10, 32, 2))
    plt.yticks(range(0, max_k_solved, 5))
    plt.tick_params(which='both', direction='in')


def rec_steps_ub_lb_bio(files, ublb, m):
    """
    plots number of recursive steps or ub-lb for bio files
    :param files: [bio33, bio50, bio66]
    :param ublb: if ub-lb should be plotted (else rec_steps)
    :param m: list of instance names (to retrieve markers and labels)
    """
    plt.grid(color='0.95', zorder=1)  # add light gray grid in the background
    plt.xlabel("Instance Count")
    my_norm = colors.Normalize(vmin=min_n_bio, vmax=max_n_bio)

    for i in range(len(files)):
        if ublb:
            s, h, p, mc, k, n = extract_values(files[i], ['Splittance (UB)', 'Heuristic (UB)', 'Packing Size (LB)',
                                                          'MinCut (LB)', 'k', 'Graph Size'])
            y_data = delta_ub_lb(s, h, p, mc)  # compute UB - LB
        else:
            y_data, k, n = extract_values(files[i], ['Recursive Steps', 'k', 'Graph Size'])

        # sort by (increasing) number of y_data (rec. steps or UB-LB):
        n = [val for (_, val) in sorted(zip(y_data, n), key=lambda z: z[0])]
        k = [val for (_, val) in sorted(zip(y_data, k), key=lambda z: z[0])]
        x = [i for i in range(1, len(y_data) + 1)]
        y_data.sort()

        # split into solved and unsolved instances:
        _, [x_solved, y_solved, n_solved], [x_unsolved, y_unsolved, n_unsolved] = split_solved_unsolved(k,
                                                                                                        [x, y_data, n])
        if not ublb:
            # add +1 to all recursive steps (in order to display 0 rec_steps on a log scale)
            y_solved = [y_solved[j] + 1 for j in range(len(y_solved))]
            y_unsolved = [y_unsolved[j] + 1 for j in range(len(y_unsolved))]
            y_data = [y_data[j] + 1 for j in range(len(y_data))]

        # create colors from graph sizes
        colors_solved = my_cmap(my_norm(n_solved))
        colors_unsolved = my_cmap(my_norm(n_unsolved))

        # plot:
        plt.plot(range(1, len(y_data) + 1), y_data, c='black', linewidth=0.3, zorder=2)
        plt.scatter(x_unsolved, y_unsolved, edgecolors=colors_unsolved, facecolors='None', marker=marker[m[i]],
                    s=size[m[i]], linewidth=.5, zorder=2)
        plt.scatter(x_solved, y_solved, edgecolors=colors_solved, facecolors=colors_solved, marker=marker[m[i]],
                    s=size[m[i]], linewidth=.5, label=label[m[i]], zorder=2)

    if not ublb:
        # make y-axis log-scale and add major and minor ticks
        ax = plt.gca()
        ax.set_yscale('log')
        y_major = LogLocator(base=10.0, numticks=5)
        ax.yaxis.set_major_locator(y_major)
        y_minor = LogLocator(base=10.0, subs=np.arange(1.0, 10.0) * 0.1, numticks=10)
        ax.yaxis.set_minor_locator(y_minor)
        ax.yaxis.set_minor_formatter(NullFormatter())

    plt.tick_params(which='both', direction='in')
    plt.legend()


def rec_steps_ub_lb_random(file, ublb: bool, er: bool):
    """
    plots number of recursive steps or ub-lb for random files
    :param file: graph file
    :param ublb: if ub-lb should be plotted (else rec_steps)
    :param er: if the random graphs are ER; else: SBM
    """
    plt.xlabel("Instance Count")
    my_norm = colors.Normalize(vmin=min_n_bio, vmax=max_n_bio)
    instance, k, n = extract_values(file, ['File', 'k', 'Graph Size'])

    if ublb:
        s, h, p, m = extract_values(file, ['Splittance (UB)', 'Heuristic (UB)', 'Packing Size (LB)', 'MinCut (LB)'])
        y_data = delta_ub_lb(s, h, p, m)  # compute UB - LB
    else:
        y_data = extract_values(file, ['Recursive Steps'])[0]

    instance_name = [instance[i].split('/')[-1] for i in range(len(instance))]  # only file name
    if er:
        f = list(map(lambda prob: int(prob[10]), instance_name))  # retrieve p value (in {1..9})
        f_map = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8}
        probs = 9
    else:
        f = list(map(lambda prob: int(prob[9:11]), instance_name))  # retrieve p value (in {80,88,95}) -> q = 1 - p
        f_map = {80: 0, 88: 1, 95: 2}
        probs = 3

    # split according to edge insertion probability used:
    n_arr, k_arr, y_arr = [[] for _ in range(probs)], [[] for _ in range(probs)], [[] for _ in range(probs)]
    for i in range(len(n)):
        n_arr[f_map[f[i]]].append(n[i])
        k_arr[f_map[f[i]]].append(k[i])
        y_arr[f_map[f[i]]].append(y_data[i])

    plt.grid(color='0.95', zorder=1)  # add light gray grid in the background

    for i in range(probs):
        # sort by (increasing) number of y_data (rec. steps or UB-LB):
        n_i = [val for (_, val) in sorted(zip(y_arr[i], n_arr[i]), key=lambda z: z[0])]
        k_i = [val for (_, val) in sorted(zip(y_arr[i], k_arr[i]), key=lambda z: z[0])]
        x = [j for j in range(1, len(y_arr[i]) + 1)]
        y = y_arr[i]
        y.sort()

        # split into solved and unsolved instances:
        _, [x_solved, y_solved, n_solved], [x_unsolved, y_unsolved, n_unsolved] = split_solved_unsolved(k_i, [x, y, n_i])

        if not ublb:
            # add +1 to all recursive steps (in order to display 0 rec_steps on a log scale)
            y_solved = [y_solved[j] + 1 for j in range(len(y_solved))]
            y_unsolved = [y_unsolved[j] + 1 for j in range(len(y_unsolved))]
            y = [y[j] + 1 for j in range(len(y))]

        # create colors from graph sizes
        colors_solved = my_cmap(my_norm(n_solved))
        colors_unsolved = my_cmap(my_norm(n_unsolved))

        # plot:
        plt.plot(range(1, len(y) + 1), y, c='black', linewidth=0.3, zorder=2)
        if er:
            prob_label = 'p=0.{}'.format(i + 1)
        else:
            prob_label = 'p=0.{}, q=0.5, r={:.2f}'.format(list(f_map.keys())[i], (100 - list(f_map.keys())[i]) / 100)

        plt.scatter(x_unsolved, y_unsolved, edgecolors=colors_unsolved, facecolors='None', marker=marker_random[i],
                    s=size_random[i], linewidth=.5, zorder=2)
        plt.scatter(x_solved, y_solved, edgecolors=colors_solved, facecolors=colors_solved, marker=marker_random[i],
                    s=size_random[i], linewidth=.5, zorder=2, label=prob_label)

    if not ublb:
        ax = plt.gca()
        ax.set_yscale('log')
        y_major = LogLocator(base=10.0, numticks=5)
        ax.yaxis.set_major_locator(y_major)
        y_minor = LogLocator(base=10.0, subs=np.arange(1.0, 10.0) * 0.1, numticks=10)
        ax.yaxis.set_minor_locator(y_minor)
        ax.yaxis.set_minor_formatter(NullFormatter())

    plt.tick_params(which='both', direction='in')
    plt.legend(fontsize=9)


def nine_plots(bio_files, er, sbm, m_bio, m_er, m_sbm):
    """
    creates 3x3 plot of solutions sizes, #rec.steps, and ub-lb for bio, ER, and SBM instances
    :param bio_files: files w/ bio data
    :param er: files w/ ER data
    :param sbm: files w/ SBM
    :param m_bio, m_er, m_sbm: instance names
    :return: None (shows plot and saves it in a file)
    """
    files = [bio_files, er, sbm]
    m = [m_bio, m_er, m_sbm]
    plt.figure(figsize=(15, 15))
    for row in range(3):            # row 0: solution sizes, row 1: recursive steps, row 2: UB-Lb
        for column in range(3):     # col 0: bio, col 1: ER, col 2: SBM
            plt.subplot(3, 3, (row*3) + column + 1)
            if row == 0:
                if column == 0:
                    print("plotting solution sizes for bio instances...")
                    plt.ylabel("Solution Size")
                    sol_sizes_bio(bio_files, m_bio)
                else:
                    print("plotting solution sizes for {} instances...".format(m[column][0]))
                    solution_sizes_random(files[column][0], column == 1)
            if row == 1:
                if column == 0:
                    print("plotting recursive steps for bio instances...")
                    plt.ylabel("Recursive Steps")
                    rec_steps_ub_lb_bio(bio_files, False, m_bio)
                else:
                    print("plotting recursive steps for {} instances...".format(m[column][0]))
                    rec_steps_ub_lb_random(files[column][0], False, column == 1)
            if row == 2:
                if column == 0:
                    print("plotting ub-lb for bio instances...")
                    plt.ylabel("Upper Bound - Lower Bound")
                    rec_steps_ub_lb_bio(bio_files, True, m_bio)
                else:
                    print("plotting ub-lb for {} instances...".format(m[column][0]))
                    rec_steps_ub_lb_random(files[column][0], True, column == 1)

        if row != 0:
            c_bar = plt.colorbar(mpl.cm.ScalarMappable(norm=colors.Normalize(vmin=min_n_bio, vmax=max_n_bio),
                                                       cmap=my_cmap), ticks=range(10, max_n_bio + 1, 5))
            c_bar.set_label('n', labelpad=-26, y=-0.025, rotation=0)

    plt.tight_layout()
    plt.savefig(out_path + 'graph_type_comparison.pdf')
    # tikzplotlib.save(out_path + 'graph_type_comparison.tikz')
    plt.show()


if __name__ == '__main__':
    # Split output files into five separate files (by instance type):
    split_file('final/final.csv')
    split_file('final/final_no_dr.csv')
    split_file('final/final_ilp.csv', ilp=True)
    split_file('final/final_ilp_no_dr.csv', ilp=True)

    # File Names:
    ilp_dr = ['final/final_ilp_bio33.csv', 'final/final_ilp_bio50.csv', 'final/final_ilp_bio66.csv', 'final/final_ilp_er.csv', 'final/final_ilp_sbm.csv']
    ilp_no_dr = ['final/final_ilp_no_dr_bio33.csv', 'final/final_ilp_no_dr_bio50.csv', 'final/final_ilp_no_dr_bio66.csv', 'final/final_ilp_no_dr_er.csv', 'final/final_ilp_no_dr_sbm.csv']
    bio = ['final/final_bio33.csv', 'final/final_bio50.csv', 'final/final_bio66.csv']
    bio_no_dr = ['final/final_no_dr_bio33.csv', 'final/final_no_dr_bio50.csv', 'final/final_no_dr_bio66.csv']
    all_types = ['bio33', 'bio50', 'bio66', 'random_ER', 'random_SBM']
    rand = ['final/final_er.csv', 'final/final_sbm.csv']
    rand_no_dr = ['final/final_no_dr_er.csv', 'final/final_no_dr_sbm.csv']
    folder = '../test_instances/final/'

    # plot runtimes of BB and ILP solver w/ and w/o data reduction applied:
    plot_runtimes_2plots([ilp_dr, bio+rand], [ilp_no_dr, bio_no_dr+rand_no_dr], all_types)

    # plot running times of ILP solver depending on number of constraints:
    ilp_runtime_against_constraints_2plots(ilp_dr[:3], ilp_dr[3:5], all_types[:3], all_types[3:5])

    # plot comparison of lower bounds:
    packing_against_min_cut_2plots(bio, rand, all_types[:3], all_types[3:5])

    # plot comparison of initial number of connected components and final number of clusters:
    number_of_clusters_2plots(bio, rand, all_types[:3], all_types[3:5])

    # print statistics about upper bounds and initial DR (creates no plots):
    optimality_of_upper_bounds(bio+rand, all_types)
    efficiency_of_initial_data_reduction(bio+rand, ilp_dr, all_types)

    nine_plots(bio, [rand[0]], [rand[1]], all_types[:3], [all_types[3]], [all_types[4]])
