// the code basis for this file comes from the file "library.cpp" in repo https://github.com/venondev/AlmostCliquePoly
#pragma once
#include <stdint.h>
#include <fstream>
#include <ext/alloc_traits.h>
#include <omp.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "common/configuration.h"
#include "common/definitions.h"
#include "tools/random_functions.h"
#include "io/graph_io.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "algorithms/global_mincut/algorithms.h"
#include "algorithms/global_mincut/minimum_cut.h"
#include "algorithms/global_mincut/cactus/cactus_mincut.h"

uint64_t min_cut(int64_t* starts, int64_t* ends, int64_t* weights, uint32_t n, uint32_t m) {
    auto cfg = configuration::getConfig();
    cfg->save_cut = true;
    random_functions::setSeed(cfg->seed);
    std::shared_ptr<mutable_graph> G = graph_io::createGraphFromEdgeList<mutable_graph>(starts, ends, weights, n, m);

    auto mc = selectMincutAlgorithm<std::shared_ptr<mutable_graph>>("noi");
    EdgeWeight cut = mc->perform_minimum_cut(G);

    return cut;
}