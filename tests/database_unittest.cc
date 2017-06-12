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

TEST (DatabaseTest, ReadTest) {
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

  std::stringstream s1, s2;

  d.DumpItems(s1);
  d.DumpPosNeg(s2);

  std::cout << s1.str();
  std::cout << s2.str();

  std::string str1, str2;
  // str1 =
  //     "010000001110111\n"
  //     "010000101010011\n"
  //     "010001001010011\n"
  //     "001110101000100\n"
  //      x      x   x  
  //     ;
  // str2 =
  //     "011000011010011\n"
  //     ;
  str1 =
      "100000111111\n"
      "100001101011\n"
      "100010101011\n"
      "011101100100\n"
      ;
  str2 =
      "110000101011\n"
      ;

  EXPECT_EQ(str1, s1.str());
  EXPECT_EQ(str2, s2.str());

  EXPECT_EQ(4, d.MaxItemInTransaction());
  EXPECT_EQ(5, d.MaxT());

  EXPECT_EQ(15, d.NuTransaction());
  EXPECT_EQ(7, d.PosTotal());
}

TEST (DatabaseTest, ReadTest2) {
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

  std::ifstream ifs1;
  ifs1.open("../samples/sample_data/sample_item.csv", std::ios::in);
  std::ifstream ifs2;
  ifs2.open("../samples/sample_data/sample_expression_over1.csv", std::ios::in);

  reader.ReadFiles(&bsh,
                   ifs1, &data, &nu_trans, &nu_items,
                   ifs2, &positive, &pos_val, &nu_pos_total,
                   item_names, transaction_names, &max_item_in_transaction, false);

  ifs1.close();
  ifs2.close();

  counters[0] = (int)(bsh->nu_bits);
  counters[1] = nu_items; 
  counters[2] = max_item_in_transaction;

  Functions4fisher functions(nu_trans, nu_pos_total, 1);
  Database<uint64> d(bsh, data, nu_trans, nu_items,
                     positive, pos_val, nu_pos_total,
                     max_item_in_transaction,
                     item_names, transaction_names, functions);

  std::stringstream s1, s2;

  d.DumpItems(s1);
  d.DumpPosNeg(s2);

  std::cout << s1.str();
  std::cout << s2.str();

  std::string str1, str2;
  // str1 =
  //     "010000001110111\n"
  //     "010000101010011\n"
  //     "010001001010011\n"
  //     "001110101000100\n"
  //     ;
  str1 =
      "100000111111\n"
      "100001101011\n"
      "100010101011\n"
      "011101100100\n"
      ;
  // str1 =
  //     "010000001110111\n"
  //     "010000101010011\n"
  //     "010001001010011\n"
  //     "001110101000100\n"
  //      x      x   x  
  //     ;
  str2 =
      "110000101011\n"
      ;

  EXPECT_EQ(str1, s1.str());
  EXPECT_EQ(str2, s2.str());

  EXPECT_EQ(4, d.MaxItemInTransaction());
  EXPECT_EQ(5, d.MaxT());

  EXPECT_EQ(15, d.NuTransaction());
  EXPECT_EQ(7, d.PosTotal());
}

// commented out because sample files
// "tfsite_both_col26.csv" and "yeast_expression_col26_over15.csv" are removed.

TEST (DatabaseTest, LargeReadTest) {
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
    ifs1.open("../samples/sample_data/yeast_col26/tfsite_both_col26.csv", std::ios::in);
    if (!ifs1.good()) {
      std::cout << "[  SKIPPED ] File not found : sample_data/yeast_col26/tfsite_both_col26.csv" << std::endl;
      return;
    }
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
    ifs2.open("../samples/sample_data/yeast_col26/yeast_expression_col26_over15.csv",
              std::ios::in);
    if (!ifs2.good()) {
      std::cout << "[  SKIPPED ] File not found : sample_data/yeast_col26/yeast_expression_col26_over15.csv" << std::endl;
      return;
    }
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

  std::stringstream s1, s2;

  d.DumpItems(s1);
  d.DumpPosNeg(s2);

  std::cout << "NuItems=" << d.NuItems() << std::endl;
  std::cout << "NuTransactions=" << d.NuTransaction() << std::endl;
  std::cout << "postotal=" << d.PosTotal() << std::endl;

  std::cout << "max x=" << d.MaxX() << std::endl;
  std::cout << "max t=" << d.MaxT() << std::endl;
  std::cout << "max item in transaction=" << d.MaxItemInTransaction() << std::endl;

  // BOOST_FOREACH (const Table::ItemInfo item, d.GetItemInfo()) {
  //   std::cout << item << std::endl;
  // }
}

TEST (DatabaseTest, PValueTest) {
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

  //d.SetSigLev(0.05);
  //d.PrepareItemVals();
  BOOST_FOREACH (const ItemInfo item, d.GetItemInfo()) {
    std::cout << item << std::endl;
  }

  std::cout << "nu all zero trans=" << d.NuAllZeroTrans() << std::endl;
  EXPECT_EQ(3, d.NuAllZeroTrans());

  EXPECT_EQ(bsh->Count(d.NthData(0)), 7u);
  EXPECT_EQ(bsh->Count(d.NthData(1)), 6u);
  EXPECT_EQ(bsh->Count(d.NthData(2)), 6u);
  EXPECT_EQ(bsh->Count(d.NthData(3)), 6u);

  uint64 * sup = bsh->New();
  bsh->Copy(d.NthData(0), sup);
  bsh->And(d.NthData(1), sup);
  bsh->And(d.NthData(2), sup);
  EXPECT_EQ(bsh->Count(sup), 5u);

  std::cout << "sup=";
  bsh->Print(std::cout, sup);
  std::cout << std::endl;

  std::cout << "DumpPMinTable:\n";
  d.DumpPMinTable(std::cout);
  std::cout << "DumpPValTable:\n";
  d.DumpPValTable(std::cout);

  uint64 * pos_sup = d.PosNeg();
  bsh->And(sup, pos_sup);
  std::cout << "pos_sup=";
  bsh->Print(std::cout, pos_sup);
  std::cout << std::endl;

  int sup_num;
  int pos_sup_num;

  sup_num = bsh->Count(sup);
  pos_sup_num = bsh->Count(pos_sup);
  EXPECT_EQ(pos_sup_num, 5);

  double p = d.PVal(sup_num, pos_sup_num, NULL, NULL);
  std::cout << "sup=" << sup_num << "\tpos_sup=" << pos_sup_num
            << "\tp=" << p << std::endl;

  double epsilon = 0.00001;
  bool flag = ( p-epsilon <= 0.006993 && 0.006993 <= p+epsilon );
  EXPECT_TRUE(flag);

  sup_num = bsh->Count(d.NthData(1));
  pos_sup_num = bsh->AndCount(d.PosNeg(), d.NthData(1));
  p = d.PVal(sup_num, pos_sup_num, NULL, NULL);

  std::cout << "sup=" << sup_num << "\tpos_sup=" << pos_sup_num
            << "\tp=" << p << std::endl;
  flag = ( p-epsilon <= 0.034965 && 0.034965 <= p+epsilon );
  EXPECT_TRUE(flag);

  sup_num = bsh->Count(d.NthData(0));
  pos_sup_num = bsh->AndCount(d.PosNeg(), d.NthData(0));
  p = d.PVal(sup_num, pos_sup_num, NULL, NULL);
  std::cout << "sup=" << sup_num << "\tpos_sup=" << pos_sup_num
            << "\tp=" << p << std::endl;
  flag = ( p-epsilon <= 0.10023 && 0.10023 <= p+epsilon );
  EXPECT_TRUE(flag);

  bsh->Delete(sup);
}

TEST (DatabaseTest, PValueTest2) {
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

  // $ ./lamp3 --item ../../../samples/alzheimer_Webster2009/alzheimer_Webster2009_SNPrec_upper05_item.csv --pos ../../../samples/alzheimer_Webster2009/alzheimer_Webster2009_SNPrec_upper05_value.csv --a 0.05 --show_progress --nosecond_phase --nothird_phase

  std::ifstream ifs1;
  ifs1.open("../samples/alzheimer_Webster2009/alzheimer_Webster2009_SNPrec_upper05_item.csv", std::ios::in);
  if (!ifs1.good()) {
    std::cout << "[  SKIPPED ] File not found : alzheimer_Webster2009/alzheimer_Webster2009_SNPrec_upper05_item.csv" << std::endl;
    return;
  }
  std::ifstream ifs2;
  ifs2.open("../samples/alzheimer_Webster2009/alzheimer_Webster2009_SNPrec_upper05_value.csv", std::ios::in);
  if (!ifs2.good()) {
    std::cout << "[  SKIPPED ] File not found : alzheimer_Webster2009/alzheimer_Webster2009_SNPrec_upper05_value.csv" << std::endl;
    return;
  }

  reader.ReadFiles(&bsh,
                   ifs1, &data, &nu_trans, &nu_items,
                   ifs2, &positive, &pos_val, &nu_pos_total,
                   item_names, transaction_names, &max_item_in_transaction, false);

  ifs1.close();
  ifs2.close();

  counters[0] = (int)(bsh->nu_bits);
  counters[1] = nu_items; 
  counters[2] = max_item_in_transaction;

  Functions4fisher functions(nu_trans, nu_pos_total, 1);
  Database<uint64> d(bsh, data, nu_trans, nu_items,
                     positive, pos_val, nu_pos_total,
                     max_item_in_transaction,
                     item_names, transaction_names, functions);

  std::stringstream s1, s2;

  EXPECT_EQ(194, d.MaxItemInTransaction());
  EXPECT_EQ(15, d.MaxX());
  EXPECT_EQ(10, d.MaxT());
  EXPECT_EQ(176, d.PosTotal());

  d.DumpPMinTable(std::cout);
}

TEST (DatabaseTest, PMinTestNoPos) {
  VariableBitsetHelper<uint64> * bsh = NULL;
  uint64 * data = NULL;
  uint64 * positive = NULL;
  double * pos_val = NULL;
  boost::array<int, 3> counters; // nu_bits, nu_items, max_item_in_transaction
  counters.assign(-1);

  int nu_trans;
  int nu_items;
  int nu_pos_total = 7;
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

  Functions4fisher functions(nu_trans, nu_pos_total, 1);
  Database<uint64> d(bsh, data, nu_trans, nu_items,
                     positive, pos_val, nu_pos_total,
                     max_item_in_transaction,
                     item_names, transaction_names, functions);

  // d.SetSigLev(0.05);
  d.PrepareItemVals();
  BOOST_FOREACH (const ItemInfo item, d.GetItemInfo()) {
    std::cout << item << std::endl;
  }

  EXPECT_EQ(bsh->Count(d.NthData(0)), 7u);
  EXPECT_EQ(bsh->Count(d.NthData(1)), 6u);
  EXPECT_EQ(bsh->Count(d.NthData(2)), 6u);
  EXPECT_EQ(bsh->Count(d.NthData(3)), 6u);

  uint64 * sup = bsh->New();
  bsh->Copy(d.NthData(0), sup);
  bsh->And(d.NthData(1), sup);
  bsh->And(d.NthData(2), sup);
  EXPECT_EQ(bsh->Count(sup), 5u);

  std::cout << "sup=";
  bsh->Print(std::cout, sup);
  std::cout << std::endl;

  std::cout << "DumpPMinTable:\n";
  d.DumpPMinTable(std::cout);
  std::cout << "DumpPValTable:\n";
  d.DumpPValTable(std::cout);

  int sup_num;
  sup_num = bsh->Count(sup);
  EXPECT_EQ(5, sup_num);

  double pmin = d.PMin(sup_num);
  std::cout << "sup=" << sup_num
            << "\tpmin=" << pmin << std::endl;

  double epsilon = 0.00001;
  // N=15, n=7, x=5
  // pmin = (7*6) / (2*1) / (15*14*13*12*11) * (5*4*3*2*1) = 0.0069930
  bool flag = ( pmin-epsilon <= 0.006993 && 0.006993 <= pmin+epsilon );
  EXPECT_TRUE(flag);

  sup_num = bsh->Count(d.NthData(1));
  EXPECT_EQ(6, sup_num);
  pmin = d.PMin(sup_num);
  std::cout << "sup=" << sup_num
            << "\tpmin=" << pmin << std::endl;

  // N=15, n=7, x=6
  // pmin = (7) / (1) / (15*14*13*12*11*10) * (6*5*4*3*2*1) = 0.0013986
  flag = ( pmin-epsilon <= 0.0013986 && 0.0013986 <= pmin+epsilon );
  EXPECT_TRUE(flag);

  sup_num = bsh->Count(d.NthData(0));
  EXPECT_EQ(7, sup_num);
  pmin = d.PMin(sup_num);
  std::cout << "sup=" << sup_num
            << "\tpmin=" << pmin << std::endl;

  // N=15, n=7, x=7
  // pmin = 1 / (15*14*13*12*11*10*9) * (7*6*5*4*3*2*1) = 0.00015540
  flag = ( pmin-epsilon <= 0.000155400 && 0.000155400 <= pmin+epsilon );
  EXPECT_TRUE(flag);

  bsh->Delete(sup);
}

TEST (DatabaseTest, ItemReverseLoopTest) {
  int N, p;
  std::vector<int> vec;
  bool flag;

  // ----
  N = 100;
  p = 10;
  vec.clear();
  vec.resize(N, 0);
  for (int r=0;r<p;r++) {
    for (int i = Database<uint64>::NextItemInReverseLoop(true, r, p, N);
         i >= 0 ; i = Database<uint64>::NextItemInReverseLoop(true, r, p, i) ) {
      vec[i] ++;
    }
  }

  std::cout << "vec:";
  for (int i=0;i<N;i++) std::cout << " " << vec[i];
  std::cout << std::endl;

  flag=true;
  for (int i=0;i<N;i++) if (vec[i] != 1) flag = false;
  EXPECT_TRUE(flag);

  // ----
  N = 91;
  p = 10;
  vec.clear();
  vec.resize(N, 0);
  for (int r=0;r<p;r++) {
    for (int i = Database<uint64>::NextItemInReverseLoop(true, r, p, N);
         i >= 0 ; i = Database<uint64>::NextItemInReverseLoop(true, r, p, i) ) {
      vec[i] ++;
    }
  }

  std::cout << "vec:";
  for (int i=0;i<N;i++) std::cout << " " << vec[i];
  std::cout << std::endl;

  flag=true;
  for (int i=0;i<N;i++) if (vec[i] != 1) flag = false;
  EXPECT_TRUE(flag);

  // ----
  N = 9;
  p = 10;
  vec.clear();
  vec.resize(N, 0);
  for (int r=0;r<p;r++) {
    for (int i = Database<uint64>::NextItemInReverseLoop(true, r, p, N);
         i >= 0 ; i = Database<uint64>::NextItemInReverseLoop(true, r, p, i) ) {
      vec[i] ++;
    }
  }

  std::cout << "vec:";
  for (int i=0;i<N;i++) std::cout << " " << vec[i];
  std::cout << std::endl;

  flag=true;
  for (int i=0;i<N;i++) if (vec[i] != 1) flag = false;
  EXPECT_TRUE(flag);

  // ----
  N = 1;
  p = 10;
  vec.clear();
  vec.resize(N, 0);
  for (int r=0;r<p;r++) {
    for (int i = Database<uint64>::NextItemInReverseLoop(true, r, p, N);
         i >= 0 ; i = Database<uint64>::NextItemInReverseLoop(true, r, p, i) ) {
      vec[i] ++;
    }
  }

  std::cout << "vec:";
  for (int i=0;i<N;i++) std::cout << " " << vec[i];
  std::cout << std::endl;

  flag=true;
  for (int i=0;i<N;i++) if (vec[i] != 1) flag = false;
  EXPECT_TRUE(flag);

  // ----
  N = 9;
  p = 1;
  vec.clear();
  vec.resize(N, 0);
  for (int r=0;r<p;r++) {
    for (int i = Database<uint64>::NextItemInReverseLoop(true, r, p, N);
         i >= 0 ; i = Database<uint64>::NextItemInReverseLoop(true, r, p, i) ) {
      vec[i] ++;
    }
  }

  std::cout << "vec:";
  for (int i=0;i<N;i++) std::cout << " " << vec[i];
  std::cout << std::endl;

  flag=true;
  for (int i=0;i<N;i++) if (vec[i] != 1) flag = false;
  EXPECT_TRUE(flag);

}

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
