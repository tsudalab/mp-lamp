// Copyright (c) 2016, Kazuki Yoshizoe
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
// may be used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// AREDISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/foreach.hpp>

#include "gflags/gflags.h"

#include "gtest/gtest.h"

#include "utils.h"
#include "lcm_graph_vba.h"
#include "lcm_dfs_vba.h"
#include "functions/Functions4fisher.h"

DECLARE_bool(show_progress); // "show progress at each iteration"

using namespace lamp_search;

TEST (LCM_DFS_LAMP_VBA_Test, SmallSearchTest) {
  std::ifstream ifs1, ifs2;
  ifs1.open("../samples/sample_data/sample_item.csv");
  ifs2.open("../samples/sample_data/sample_expression_over1.csv");
  Functions4fisher functions(1);
  TableVBA<uint64> t(ifs1, ifs2, functions);
  ifs1.close();
  ifs2.close();

  LCM_Graph_VBA<uint64> g(t);
  LCM_DFS_VBA search(g);

  search.LAMP();

  search.PrintLAMPResults(std::cout);
}

TEST (LCM_DFS_LAMP_VBA_Test, SmallSingleSearchTest) {
  std::ifstream ifs1, ifs2;
  ifs1.open("../samples/sample_data/sample_item.csv");
  ifs2.open("../samples/sample_data/sample_expression_over1.csv");
  Functions4fisher functions(1);
  TableVBA<uint64> t(ifs1, ifs2, functions);
  ifs1.close();
  ifs2.close();

  LCM_Graph_VBA<uint64> g(t);
  LCM_DFS_VBA search(g);

  search.LAMPSingle();

  search.PrintLAMPResults(std::cout);
  //search.PrintLAMPSingleResults(std::cout);
}

// TEST (LCM_DFS_LAMP_VBA_Test, BugTest20150905) {
//   // $ ./lcm_dfs_lamp_vba --item ../samples/sample_data/yeast_col26/tfsite_both_col26.csv --pos ../samples/sample_data/yeast_col26/yeast_expression_col26_over15.csv --a 0.05 --threshold_init 12 --show_progress > t2h.txt

//   std::ifstream ifs1, ifs2;
//   ifs1.open("../samples/sample_data/yeast_col26/tfsite_both_col26.csv");
//   ifs2.open("../samples/sample_data/yeast_col26/yeast_expression_col26_over15.csv");
//   TableVBA<uint64> t(ifs1, ifs2);
//   ifs1.close();
//   ifs2.close();

//   LCM_Graph_VBA<uint64> g(t);
//   LCM_DFS_VBA search(g);

//   // ----------

//   const int default_value = std::numeric_limits<int>::max();
//   const TableVBA<uint64> & tbl = GetGraph().GetTable();

//   // define init threshold
//   int sup_threshold = 12;
//   int num_iter = 0;
//   total_expand_num_ = 0ll;

//   LAMPInit();
//   LCM_DFS_VBA_State state;
//   state.Init();

//   {
//     expand_num_ = 0ll;
//     closed_set_num_ = 0ll;

//     double pmin_threshold = tbl.PMin( sup_threshold );
//     long long int closed_set_num_threshold =
//         (long long int)(std::min(std::floor(FLAGS_a/pmin_threshold),
//                                  (double)(default_value)));

//     //hoge;
//     sup_stack_.Init(); // make stack for root position
//     LAMPIter(&state, sup_threshold, closed_set_num_threshold, hoge, true);



//   }


// }

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);
  return RUN_ALL_TESTS();
}

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
