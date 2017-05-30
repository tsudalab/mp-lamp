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

#include "gtest/gtest.h"

#include "utils.h"
#include "lcm_graph.h"
#include "lcm_dfs.h"
#include "functions/Functions4fisher.h"

using namespace lamp_search;

TEST (LCM_DFSTest, TopKInitTest) {
  std::ifstream ifs1, ifs2;
  ifs1.open("../../../samples/sample_data/sample_item.csv");
  ifs2.open("../../../samples/sample_data/sample_expression_over1.csv");
  Functions4fisher functions(1);
  Table t(ifs1, ifs2, functions);
  ifs1.close();
  ifs2.close();
  LCM_Graph g(t);
  LCM_DFS search(g);

  std::cout << search.TopK();
  //search.DumpTopK(std::cout);

  search.SetTopK(4);
  std::cout << search.TopK();
  //search.DumpTopK(std::cout);

  search.SetTopK(5);
  std::cout << search.TopK();
  //search.DumpTopK(std::cout);
}

TEST (LCM_DFSTest, TopKInitInsert) {
  std::ifstream ifs1, ifs2;
  ifs1.open("../../../samples/sample_data/sample_item.csv");
  ifs2.open("../../../samples/sample_data/sample_expression_over1.csv");
  Functions4fisher functions(1);
  Table t(ifs1, ifs2, functions);
  ifs1.close();
  ifs2.close();
  LCM_Graph g(t);
  LCM_DFS search(g);

  search.SetTopK(10);
  TopKData & topk = search.TopK();
  
  typedef boost::dynamic_bitset<> dbs;
  SortedItemSet ss1;
  ss1.Push(0); ss1.Push(1); ss1.Push(2); ss1.Push(3);
  SortedItemSet ss2;
  ss2.Push(0); ss2.Push(1); ss2.Push(2);
  SortedItemSet ss3;
  ss3.Push(0); ss3.Push(2); ss3.Push(3);
  SortedItemSet ss4;
  ss4.Push(0); ss4.Push(1); ss4.Push(3);
  topk.Insert(0.008, ss1);
  topk.Insert(0.006, ss2);
  topk.Insert(0.002, ss3);
  topk.Insert(0.004, ss4);

  std::cout << search.TopK();

  const TopKData::val_index & p0 = topk.Nth(0);
  const TopKData::val_index & p1 = topk.Nth(1);
  const TopKData::val_index & p2 = topk.Nth(2);
  const TopKData::val_index & p3 = topk.Nth(3);
  std::cout << p0.first << ":" << p0.second << std::endl;
  std::cout << p1.first << ":" << p1.second << std::endl;
  std::cout << p2.first << ":" << p2.second << std::endl;
  std::cout << p3.first << ":" << p3.second << std::endl;

  EXPECT_EQ( p0.first, 0.002);
  EXPECT_EQ( topk.NthItemset(0), ss3);

  EXPECT_EQ( p1.first, 0.004);
  EXPECT_EQ( topk.NthItemset(1), ss4);

  EXPECT_EQ( p2.first, 0.006);
  EXPECT_EQ( topk.NthItemset(2), ss2);

  EXPECT_EQ( p3.first, 0.008);
  EXPECT_EQ( topk.NthItemset(3), ss1);
}

TEST (LCM_DFSTest, SmallSearchTest) {
  std::ifstream ifs1, ifs2;
  ifs1.open("../../../samples/sample_data/sample_item.csv");
  ifs2.open("../../../samples/sample_data/sample_expression_over1.csv");
  Functions4fisher functions(1);
  Table t(ifs1, ifs2, functions);
  ifs1.close();
  ifs2.close();
  LCM_Graph g(t);
  LCM_DFS search(g);

  search.SetTopK(10);
  search.Search();

  std::cout << search.TopK();
}
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
