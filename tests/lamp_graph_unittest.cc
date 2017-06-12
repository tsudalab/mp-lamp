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
#include "timer.h"
#include "random.h"

#include "sorted_itemset.h"
#include "variable_bitset_array.h"
#include "variable_length_itemset.h"

#include "database.h"
#include "lcm_graph_vba.h"
#include "lamp_graph.h"
#include "functions/Functions4fisher.h"

using namespace lamp_search;

TEST (LampGraphTest, ToyDataTest) {
  VariableBitsetHelper<uint64> * bsh = NULL;

  uint64 * data = NULL;
  uint64 * positive = NULL;
  double * pos_val = NULL;
  boost::array<int, 3> counters; // nu_bits, nu_items, max_item_in_transaction
  counters.assign(-1);

  int nu_trans;
  int nu_items;
  int nu_pos_total = 0;
  int max_item_in_transaction;

  std::vector< std::string > * item_names = NULL;
  std::vector< std::string > * transaction_names = NULL;
  item_names = new std::vector< std::string >;
  transaction_names = new std::vector< std::string >;

  DatabaseReader<uint64> reader;

  {
    std::ifstream ifs1;
    ifs1.open("../samples/sample_data/sample_item.csv", std::ios::in);
    // read file into uint64 * array and prepare database
    reader.ReadItems(ifs1, &nu_trans, &nu_items, &bsh,
                     &data, item_names, transaction_names, &max_item_in_transaction);
    counters[0] = (int)(bsh->nu_bits);
    counters[1] = nu_items; 
    counters[2] = max_item_in_transaction;
    ifs1.close();
  }

  {
    std::ifstream ifs2;
    ifs2.open("../samples/sample_data/sample_expression_over1.csv", std::ios::in);
    // read positives into uint64 * array and prepare database
    reader.ReadPosNeg(ifs2, nu_trans, transaction_names,
                      &nu_pos_total, bsh, &positive, &pos_val, false);
    // positive_count = bsh->NuBlocks();
    ifs2.close();
  }

  Functions4fisher functions(nu_trans, nu_pos_total, 1);
  Database<uint64> d(bsh, data, nu_trans, nu_items,
                     positive, pos_val, nu_pos_total,
                     max_item_in_transaction,
                     item_names, transaction_names, functions);
  LampGraph<uint64> g(d);

  //--------

  int x = d.MaxX();
  VariableLengthItemsetStack st(1024*1024);
  // VariableLengthItemsetStack st(1024*1024, x);

  st.PushPre();
  int * i0 = st.Top();
  st.SetSup(i0, x);
  st.PushPostNoSort();

  EXPECT_EQ(-1, g.CoreIndex(st, i0));

  st.PushPre(); // int * i0p1 = st.Top();
  st.PushOneItem(0); st.PushPostNoSort();

  st.PushPre(); // int * i0p2 = st.Top();
  st.PushOneItem(1); st.PushPostNoSort();

  st.PushPre(); // int * i0p3 = st.Top();
  st.PushOneItem(2); st.PushPostNoSort();

  st.PushPre(); // int * i0p4 = st.Top();
  st.PushOneItem(3); st.PushPostNoSort();

  //----

  st.PushPre(); int * i1 = st.Top();
  st.PushOneItem(0);
  st.PushPostNoSort();

  st.PushPre(); int * i12 = st.Top();
  st.PushOneItem(0);
  st.PushOneItem(1);
  st.PushPostNoSort();

  st.PushPre(); int * i13 = st.Top();
  st.PushOneItem(0);
  st.PushOneItem(2);
  st.PushPostNoSort();

  EXPECT_EQ(0, g.CoreIndex(st, i1) );
  EXPECT_EQ(1, g.CoreIndex(st, i12));
  EXPECT_EQ(2, g.CoreIndex(st, i13));



}

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
