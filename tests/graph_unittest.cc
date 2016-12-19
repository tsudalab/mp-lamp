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

#include "sorted_itemset.h"
#include "topk.h"

#include "table.h"
#include "graph.h"

using namespace lamp_search;

TEST (GraphTest, ReadTest) {
  std::ifstream ifs1, ifs2;
  ifs1.open("../../../samples/sample_data/sample_item.csv");
  ifs2.open("../../../samples/sample_data/sample_expression_over1.csv");
  Table t(ifs1, ifs2);
  ifs1.close();
  ifs2.close();

  std::stringstream s1, s2;

  t.DumpItems(s1);
  t.DumpPosNeg(s2);

  std::cout << s1.str();
  std::cout << s2.str();

  std::string str1, str2;
  str1 =
      "010000001110111\n"
      "010000101010011\n"
      "010001001010011\n"
      "001110101000100\n"
      ;
  str2 =
      "011000011010011\n"
      ;

  EXPECT_EQ(str1, s1.str());
  EXPECT_EQ(str2, s2.str());

  EXPECT_EQ(t.MaxItemInTransaction(), 4);
  EXPECT_EQ(t.MaxT(), 5);
}

// P-value
// 0.006993	TF1,TF2,TF3	5	5
// 0.034965	TF2	6	5
// 0.034965	TF3	6	5
// 0.10023	TF1	7	5

// {TF1,2}, {TF1,3}, {TF2,3} are same with {TF1,2,3}

TEST (GraphTest, PValueNanTest) {
  std::ifstream ifs1, ifs2;
  ifs1.open("../../../samples/sample_data/sample_item.csv");
  ifs2.open("../../../samples/sample_data/sample_expression_over1.csv");
  Table t(ifs1, ifs2);
  ifs1.close();
  ifs2.close();

  // NuItems=102
  // NuTransactions=6074
  // postotal=530
  // max x=41
  // max item in transaction=10
  // id=56:pmin=6.58056e-263:pval=nan:sup=226:pos_sup=10
  // id=0:pmin=9.85092e-229:pval=nan:sup=199:pos_sup=8

  double p, m;

  // int nu_item, int nu_transaction, int nu_pos_total
  t.SetValuesForTest(102, 6074, 530);
  t.InitPMinTable();

  m = t.PMinCal(226);
  p = t.PValCal(226, 10);
  std::cout << "p=" << p << "\tm=" << m << std::endl;
  EXPECT_TRUE(!std::isnan(p));

  m = t.PMin(199);
  p = t.PVal(199, 8);
  std::cout << "p=" << p << "\tm=" << m << std::endl;
  EXPECT_TRUE(!std::isnan(p));

  // int nu_item, int nu_transaction, int nu_pos_total
  t.SetValuesForTest(397, 12773, 1129);
  t.InitPMinTable();

  m = t.PMin(785);
  p = t.PVal(785, 93);
  std::cout << "p=" << p << "\tm=" << m << std::endl;
  EXPECT_TRUE(!std::isnan(p));
}

TEST (GraphTest, PValueTest) {
  std::ifstream ifs1, ifs2;
  ifs1.open("../../../samples/sample_data/sample_item.csv");
  ifs2.open("../../../samples/sample_data/sample_expression_over1.csv");
  Table t(ifs1, ifs2);
  ifs1.close();
  ifs2.close();

  t.SetSigLev(0.05);
  t.PrepareItemVals();
  BOOST_FOREACH (const ItemInfo item, t.GetItemInfo()) {
    std::cout << item << std::endl;
  }

  EXPECT_EQ(t.NthData(0).count(), 7u);
  EXPECT_EQ(t.NthData(1).count(), 6u);
  EXPECT_EQ(t.NthData(2).count(), 6u);
  EXPECT_EQ(t.NthData(3).count(), 6u);

  boost::dynamic_bitset<> sup = t.NthData(0);
  sup &= t.NthData(1);
  sup &= t.NthData(2);
  EXPECT_EQ(sup.count(), 5u);
  std::cout << sup << std::endl;

  t.DumpPMinTable(std::cout);
  t.DumpPValTable(std::cout);

  boost::dynamic_bitset<> pos_sup = sup & t.PosNeg();

  int sup_num;
  int pos_sup_num;

  sup_num = sup.count();
  pos_sup_num = pos_sup.count();
  double p = t.PVal(sup_num, pos_sup_num);
  std::cout << "sup=" << sup_num << "\tpos_sup=" << pos_sup_num
            << "\tp=" << p << std::endl;

  double epsilon = 0.00001;
  bool flag = ( p-epsilon <= 0.006993 && 0.006993 <= p+epsilon );
  EXPECT_TRUE(flag);

  sup_num = t.NthData(1).count();
  pos_sup_num = (t.NthData(1) & t.PosNeg()).count();
  p = t.PVal(sup_num, pos_sup_num);;
  std::cout << "sup=" << sup_num << "\tpos_sup=" << pos_sup_num
            << "\tp=" << p << std::endl;
  flag = ( p-epsilon <= 0.034965 && 0.034965 <= p+epsilon );
  EXPECT_TRUE(flag);

  sup_num = t.NthData(0).count();
  pos_sup_num = (t.NthData(0) & t.PosNeg()).count();
  p = t.PVal(sup_num, pos_sup_num);;
  std::cout << "sup=" << sup_num << "\tpos_sup=" << pos_sup_num
            << "\tp=" << p << std::endl;
  flag = ( p-epsilon <= 0.10023 && 0.10023 <= p+epsilon );
  EXPECT_TRUE(flag);
}

// commented out because sample files are removed.

// TEST (GraphTest, LargeReadTest2) {
//   std::ifstream ifs1, ifs2;
//   ifs1.open("../../../samples/mcf7_transcriptome/egf_motif_item.csv");
//   ifs2.open("../../../samples/mcf7_transcriptome/egf_motif_value.csv");
//   Table t(ifs1, ifs2);
//   ifs1.close();
//   ifs2.close();

//   std::stringstream s1, s2;

//   t.DumpItems(s1);
//   t.DumpPosNeg(s2);

//   //std::cout << s1.str();
//   //std::cout << s2.str();

//   // EXPECT_EQ(str1, s1.str());
//   // EXPECT_EQ(str2, s2.str());

//   std::cout << "NuItems=" << t.NuItems() << std::endl;
//   std::cout << "NuTransactions=" << t.NuTransaction() << std::endl;
//   std::cout << "postotal=" << t.PosTotal() << std::endl;

//   std::cout << "max x=" << t.MaxX() << std::endl;
//   std::cout << "max item in transaction=" << t.MaxItemInTransaction() << std::endl;

//   BOOST_FOREACH (const ItemInfo item, t.GetItemInfo()) {
//     std::cout << item << std::endl;
//   }
// }

// commented out because sample files are removed.

// TEST (GraphTest, LargeReadTest) {
//   std::ifstream ifs1, ifs2;
//   ifs1.open("../../../samples/sample_data/yeast_col26/tfsite_both_col26.csv");
//   ifs2.open("../../../samples/sample_data/yeast_col26/yeast_expression_col26_over15.csv");
//   Table t(ifs1, ifs2);
//   ifs1.close();
//   ifs2.close();

//   std::stringstream s1, s2;

//   t.DumpItems(s1);
//   t.DumpPosNeg(s2);

//   //std::cout << s1.str();
//   //std::cout << s2.str();

//   // EXPECT_EQ(str1, s1.str());
//   // EXPECT_EQ(str2, s2.str());

//   std::cout << "NuItems=" << t.NuItems() << std::endl;
//   std::cout << "NuTransactions=" << t.NuTransaction() << std::endl;
//   std::cout << "postotal=" << t.PosTotal() << std::endl;

//   std::cout << "max x=" << t.MaxX() << std::endl;
//   std::cout << "max t=" << t.MaxT() << std::endl;
//   std::cout << "max item in transaction=" << t.MaxItemInTransaction() << std::endl;

//   // BOOST_FOREACH (const Table::ItemInfo item, t.GetItemInfo()) {
//   //   std::cout << item << std::endl;
//   // }
// }

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
