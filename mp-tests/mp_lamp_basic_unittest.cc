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

#include "mp_dfs.h"

DEFINE_bool(lcm, false, "item file is lcm style");

DEFINE_string(item, "", "filename of item set");
DEFINE_string(pos, "", "filename of positive / negative file");

DEFINE_int32(n, 1000, "granularity of one Node process");
DEFINE_bool(n_is_ms, true, "true: n is milli sec, false: n is num task");
DEFINE_int32(w, 1, "number of random steal attempts");
DEFINE_int32(m, 3, "number of maximum random steal candidates");
DEFINE_int32(l, 2, "power of lifeline graph");

DECLARE_bool(log); // false, "show log", mp-lamp.cc

using namespace lamp_search;

TEST (MPI_BASIC_TEST, InitTest) {
  int rank, nu_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nu_proc);

  // (int rank, int nu_proc, int n, bool n_is_ms, int w, int l, int m)
  MP_LAMP search(rank, nu_proc, FLAGS_n, FLAGS_n_is_ms, FLAGS_w, FLAGS_l, FLAGS_m);

  if (rank==0) {
    std::ifstream ifs1, ifs2;
    ifs1.open("../../../samples/sample_data/sample_item.csv", std::ios::in);
    ifs2.open("../../../samples/sample_data/sample_expression_over1.csv", std::ios::in);
    search.InitDatabaseRoot(ifs1, ifs2);
    ifs1.close();
    ifs2.close();
  } else {
    search.InitDatabaseSub(true);
  }

  for(int i = 0; i<nu_proc; i++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (i == rank) {
      search.GetDatabase().DumpItems(std::cout);
      search.GetDatabase().ShowInfo(std::cout);
    }
  }
}

TEST (MPI_BASIC_TEST, InitTest2) {
  int rank, nu_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nu_proc);

  // (int rank, int nu_proc, int n, bool n_is_ms, int w, int l, int m)
  MP_LAMP search(rank, nu_proc, FLAGS_n, FLAGS_n_is_ms, FLAGS_w, FLAGS_l, FLAGS_m);

  if (rank==0) {
    std::ifstream ifs1;
    ifs1.open("../../../samples/sample_data/sample_item.csv", std::ios::in);
    search.InitDatabaseRoot(ifs1, 7);
    ifs1.close();
  } else {
    search.InitDatabaseSub(false);
  }

  for(int i = 0; i<nu_proc; i++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (i == rank) {
      search.GetDatabase().DumpItems(std::cout);
      search.GetDatabase().ShowInfo(std::cout);
    }
  }
}

TEST (MPI_BASIC_TEST, InitTestLarge) {
  int rank, nu_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nu_proc);

  // (int rank, int nu_proc, int n, bool n_is_ms, int w, int l, int m)
  MP_LAMP search(rank, nu_proc, FLAGS_n, FLAGS_n_is_ms, FLAGS_w, FLAGS_l, FLAGS_m);

  if (rank==0) {
    std::ifstream ifs1, ifs2;
    ifs1.open("../../../samples/sample_data/yeast_col26/tfsite_both_col26.csv");
    ifs2.open("../../../samples/sample_data/yeast_col26/yeast_expression_col26_over15.csv");
    search.InitDatabaseRoot(ifs1, ifs2);
    ifs1.close();
    ifs2.close();
  } else {
    search.InitDatabaseSub(true);
  }

  for(int i = 0; i<nu_proc; i++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (i == rank) {
      //search.GetDatabase().DumpItems(std::cout);
      search.GetDatabase().ShowInfo(std::cout);
    }
  }  
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);

  int res = RUN_ALL_TESTS();
  MPI_Finalize();
  return res;
}

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
