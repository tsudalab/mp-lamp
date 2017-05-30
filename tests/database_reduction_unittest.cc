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
#include <cmath>

#include "gtest/gtest.h"

#include "variable_length_itemset.h"
#include "database.h"
#include "functions/Functions4fisher.h"

using namespace lamp_search;

TEST(DatabaseReductionTest, SimpleReduceTest) {
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
    ifs1.open("../../../samples/sample_data/sample_item.csv", std::ios::in);
    // read file into uint64 * array and prepare database
    reader.ReadItems(ifs1, &nu_trans, &nu_items, &bsh,
                     &data, item_names, transaction_names, &max_item_in_transaction);
    counters[0] = (int) (bsh->NuBits());
    counters[1] = nu_items;
    counters[2] = max_item_in_transaction;
    ifs1.close();
  }

  {
    std::ifstream ifs2;
    ifs2.open("../../../samples/sample_data/sample_expression_over1.csv", std::ios::in);
    // read positives into uint64 * array and prepare database
    reader.ReadPosNeg(ifs2, nu_trans, transaction_names,
                      &nu_pos_total, bsh, &positive, &pos_val, false);
    // positive_count = bsh->NuBlocks();
    ifs2.close();
  }

  // std::string str1, str2;
  // str1 =
  //     "010000001110111\n" // 7
  //     "010000101010011\n" // 6
  //     "010001001010011\n" // 6
  //     "001110101000100\n" // 6
  //     ;
  // str2 =
  //     "011000011010011\n"
  //     ;

  Functions4fisher functions(nu_trans, nu_pos_total, 1);
  Database<uint64> d(bsh, data, nu_trans, nu_items,
                     positive, pos_val, nu_pos_total,
                     max_item_in_transaction,
                     item_names, transaction_names, functions);

  //d.SetSigLev(0.05);
  //d.PrepareItemVals();

  BOOST_FOREACH(const ItemInfo item, d.GetItemInfo()) {
    std::cout << item << std::endl;
  }

  std::cout << "nu all zero trans=" << d.NuAllZeroTrans() << std::endl;
  EXPECT_EQ(3, d.NuAllZeroTrans());
  EXPECT_EQ(0, d.NuAllZeroItems());

  std::vector<int> & v1 = d.ReducedList(1);
  std::cout << "v1:";

  BOOST_FOREACH(int i, v1) {
    std::cout << " " << i;
  }
  std::cout << std::endl;

  EXPECT_EQ(4u, v1.size());
  EXPECT_EQ(0, v1[0]);
  EXPECT_EQ(1, v1[1]);
  EXPECT_EQ(2, v1[2]);
  EXPECT_EQ(3, v1[3]);

  std::vector<int> & v7 = d.ReducedList(7);
  std::cout << "v7:";

  BOOST_FOREACH(int i, v7) {
    std::cout << " " << i;
  }
  std::cout << std::endl;

  EXPECT_EQ(1u, v7.size());
  EXPECT_EQ(0, v7[0]);

  EXPECT_EQ(4, d.NuReducedItems(1));
  EXPECT_EQ(-1, d.NuReducedItems(6)); // not prepared
  d.PrepareReducedList(6);
  EXPECT_EQ(4, d.NuReducedItems(6)); // prepared
  EXPECT_EQ(1, d.NuReducedItems(7));

  const std::vector<int> & h = d.SupCumHist();
  std::cout << "hist:";

  BOOST_FOREACH(int i, h) {
    std::cout << " " << i;
  }
  std::cout << std::endl;

  std::cout << "num1=" << d.Count1() << std::endl;
  EXPECT_EQ(25, d.Count1());

  std::cout << "density=" << d.Density() << std::endl;
  EXPECT_EQ((25.0 / 48.0), d.Density());

  uint64 * sup = bsh->New();
  bsh->Delete(sup);
}

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
