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
#include "table_vba.h"

using namespace lamp_search;

TEST (TableVBATest, ReadTest) {
  std::ifstream ifs1, ifs2;
  ifs1.open("../../../samples/sample_data/sample_item.csv");
  ifs2.open("../../../samples/sample_data/sample_expression_over1.csv");
  TableVBA<uint64> t(ifs1, ifs2);
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

TEST (TableVBATest, PValueTest) {
  typedef VariableBitset<uint64> VBS;
  typedef VariableBitsetArray<uint64> VBA;

  std::ifstream ifs1, ifs2;
  ifs1.open("../../../samples/sample_data/sample_item.csv");
  ifs2.open("../../../samples/sample_data/sample_expression_over1.csv");
  TableVBA<uint64> t(ifs1, ifs2);
  ifs1.close();
  ifs2.close();

  t.SetSigLev(0.05);
  t.PrepareItemVals();
  BOOST_FOREACH (const ItemInfo item, t.GetItemInfo()) {
    std::cout << item << std::endl;
  }

  EXPECT_EQ(t.Data()->Count(t.NthData(0)), 7u);
  EXPECT_EQ(t.Data()->Count(t.NthData(1)), 6u);
  EXPECT_EQ(t.Data()->Count(t.NthData(2)), 6u);
  EXPECT_EQ(t.Data()->Count(t.NthData(3)), 6u);

  VBS sup(t.NuTransaction());
  sup.CopyFrom(t.NthData(0));
  t.Data()->And(t.NthData(1), &sup);
  t.Data()->And(t.NthData(2), &sup);
  EXPECT_EQ(sup.Count(), 5u);
  std::cout << "sup=";
  std::cout << sup << std::endl;

  std::cout << "DumpPMinTable:\n";
  t.DumpPMinTable(std::cout);
  std::cout << "DumpPValTable:\n";
  t.DumpPValTable(std::cout);

  VBS * pos_sup = t.PosNeg();
  pos_sup->And(sup);
  std::cout << "pos_sup=";
  std::cout << *pos_sup << std::endl;

  int sup_num;
  int pos_sup_num;

  sup_num = sup.Count();
  pos_sup_num = pos_sup->Count();
  EXPECT_EQ(pos_sup_num, 5);

  double p = t.PVal(sup_num, pos_sup_num);
  std::cout << "sup=" << sup_num << "\tpos_sup=" << pos_sup_num
            << "\tp=" << p << std::endl;

  double epsilon = 0.00001;
  bool flag = ( p-epsilon <= 0.006993 && 0.006993 <= p+epsilon );
  EXPECT_TRUE(flag);

  sup_num = t.Data()->Count(t.NthData(1));
  pos_sup_num = t.Data()->AndCount(*(t.PosNeg()), t.NthData(1));
  p = t.PVal(sup_num, pos_sup_num);;

  std::cout << "sup=" << sup_num << "\tpos_sup=" << pos_sup_num
            << "\tp=" << p << std::endl;
  flag = ( p-epsilon <= 0.034965 && 0.034965 <= p+epsilon );
  EXPECT_TRUE(flag);

  sup_num = t.Data()->Count(t.NthData(0));
  pos_sup_num = t.Data()->AndCount(*(t.PosNeg()), t.NthData(0));
  p = t.PVal(sup_num, pos_sup_num);;
  std::cout << "sup=" << sup_num << "\tpos_sup=" << pos_sup_num
            << "\tp=" << p << std::endl;
  flag = ( p-epsilon <= 0.10023 && 0.10023 <= p+epsilon );
  EXPECT_TRUE(flag);
}

// TEST (TableVBATest, LargeReadTest) {
//   std::ifstream ifs1, ifs2;
//   ifs1.open("../../../samples/sample_data/yeast_col26/tfsite_both_col26.csv");
//   ifs2.open("../../../samples/sample_data/yeast_col26/yeast_expression_col26_over15.csv");
//   TableVBA<uint64> t(ifs1, ifs2);
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

TEST (TableVBATest, PvalLogTest) {
  typedef VariableBitset<uint64> VBS;
  typedef VariableBitsetArray<uint64> VBA;

  std::ifstream ifs1, ifs2;
  // --item ../local/yeast_lamp/yeast_lamp_data/bind_genes/tfsite_both_col8.csv
  // --pos ../local/yeast_lamp/yeast_lamp_data/exp/yeast_expression_col8_over10.csv

  ifs1.open("../../../local/yeast_lamp/yeast_lamp_data/bind_genes/tfsite_both_col8.csv");
  ifs2.open("../../../local/yeast_lamp/yeast_lamp_data/exp/yeast_expression_col8_over10.csv");
  TableVBA<uint64> t(ifs1, ifs2);
  ifs1.close();
  ifs2.close();

  t.SetSigLev(0.05);
  t.PrepareItemVals();
  
  const double epsilon = 0.0001;

  for (int i=1;i<=t.NuTransaction();i++) {
    double a = t.PMinCal(i);
    double b = exp( t.PMinCalLog(i) );
    bool cond = (1-epsilon <= a/b && a/b <= 1+epsilon);
    //std::cout << "a/b=" << a/b << std::endl;
    if (!cond)
      std::cout << "a=" << a << "\tb=" << b
                << "\ti=" << i
                << std::endl;
    //EXPECT_TRUE(cond);
  }

  for (int i=1;i<=t.MaxX();i++) { // sup == x
    for (int j=1;j<=t.MaxT();j++) { // pos_sup == t
      if (j>i) continue;
      double a = t.PValCal(i, j);
      double b = t.PValCalLog(i, j);
      bool cond = (1-epsilon <= a/b && a/b <= 1+epsilon);
      // std::cout << "a/b=" << a/b << std::endl;
      if (!cond)
        std::cout << "a=" << a << "\tb=" << b
                  << "\ti=" << i << "\tj=" << j
                  << std::endl;
      // EXPECT_TRUE(cond);
    }
  }
}

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
