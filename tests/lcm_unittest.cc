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

#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/foreach.hpp>

#include "gtest/gtest.h"

#include "utils.h"
#include "table.h"
#include "lcm_graph.h"

#include "variable_bitset_array.h"
#include "table_vba.h"
#include "lcm_graph_vba.h"
#include "functions/FunctionsSuper.h"
#include "functions/Functions4fisher.h"

using namespace lamp_search;

namespace {

static const int END = std::numeric_limits<int>::max();
typedef SortedItemSet::ItemSetType ba;

static void PrintItemSet(std::ostream & out, const SortedItemSet & items) {
  std::stringstream s;
  s << "items:";
  for (SortedItemSet::const_iterator it=items.begin();it!=items.end();++it) {
    s << " " << *it;
  }
  // for (int i=0 ; i<kMaxSearchDepth ; i++) {
  //   if (items[i] >= E) break;
  //   s << " " << items[i];
  // }
  out << s.str() << std::endl;
}

static void FillItemSet(SortedItemSet & items) {
  int i=0;
  for (; i<kMaxSearchDepth ; i++) {
    if (items[i] >= END) break;
  }
  for (; i<kMaxSearchDepth ; i++) {
    items[i] = END;
  }
}

static void ReadSampleTestData(Table ** t, FunctionsSuper **functions) {
  std::stringstream item_ss, pos_ss;

  item_ss << "#test,i0,i1,i2,i3,i4,i5,i6\n"
      //         0 1 2 3 4 5 6 (0 is dummy)
          << "t1,0,1,1,1,1,1,1\n" // 1 2 3 4 5 6
          << "t2,0,0,1,1,0,1,0\n" // 2 3 5
          << "t3,0,0,1,0,0,1,0\n" // 2 5
          << "t4,0,1,1,0,1,1,1\n" // 1 2 4 5 6
          << "t5,0,0,1,0,1,0,0\n" // 2 4
          << "t6,0,1,0,0,1,0,1\n" // 1 4 6
          << "t7,0,0,0,1,1,0,1\n";// 3 4 6
  
  // random dummy data
  pos_ss  << "#test,expression\n"
          << "t1,1\n"
          << "t2,1\n"
          << "t3,0\n"
          << "t4,0\n"
          << "t5,0\n"
          << "t6,0\n"
          << "t7,0\n";

  *functions = new Functions4fisher(1);
  *t = new Table(item_ss, pos_ss, **functions);
}

} // anonymous namespace

TEST (LCMTest, SampleDataTest) {
  Table * t;
  FunctionsSuper * functions;
  ReadSampleTestData(&t, &functions);

  std::stringstream ss;
  t->DumpItems(ss);

  std::cout << ss.str();

  std::string ss_comp_str =
      "0000000\n"
      "0101001\n"
      "0011111\n"
      "1000011\n"
      "1111001\n"
      "0001111\n"
      "1101001\n";

  EXPECT_EQ(ss.str(), ss_comp_str);

  delete t;
  delete functions;
}

TEST (LCMTest, ToyDataTest) {
  std::ifstream ifs1, ifs2;
  ifs1.open("../samples/sample_data/sample_item.csv");
  ifs2.open("../samples/sample_data/sample_expression_over1.csv");
  Functions4fisher functions(1);
  Table t(ifs1, ifs2, functions);
  ifs1.close();
  ifs2.close();
  LCM_Graph g(t);
  
  SortedItemSet i0; // -1 (by definition, -1 added for lamp_search)
  ba a0p1  = {{0,END}}; SortedItemSet i0p1(a0p1); // ppc-ext of (i0 \cup 1)
  ba a0p2  = {{1,END}}; SortedItemSet i0p2(a0p2); // ppc-ext of (i0 \cup 2)
  ba a0p3  = {{2,END}}; SortedItemSet i0p3(a0p3); // ppc-ext of (i0 \cup 3)
  ba a0p4  = {{3,END}}; SortedItemSet i0p4(a0p4); // ppc-ext of (i0 \cup 4)
  EXPECT_EQ(g.CoreIndex(i0),   -1);
  
  ba a1   = {{0,END}};           SortedItemSet i1(a1);   // 0
  ba a12  = {{0,1,END}};         SortedItemSet i12(a12); // 1
  ba a13  = {{0,2,END}};         SortedItemSet i13(a13); // 2
  EXPECT_EQ(g.CoreIndex(i1),   0);
  EXPECT_EQ(g.CoreIndex(i12),  1);
  EXPECT_EQ(g.CoreIndex(i13),  2);

  SortedItemSet buf;
  bool res=false;
  boost::dynamic_bitset<> sup(t.NuTransaction());
  //boost::dynamic_bitset<> sup_buf(t.NuTransaction());
  std::string sup_str_buf;
  int ni = 0;

  g.Support(i0, &sup);
  std::cout << "sup=" << sup << std::endl;

  ni = 0; // 1
  res = g.PPCExtension(i0, (sup & t.NthData(ni)), g.CoreIndex(i0), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i0p1, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("010000001110111", sup_str_buf);

  ni = 1; // 2
  res = g.PPCExtension(i0, (sup & t.NthData(ni)), g.CoreIndex(i0), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i0p2, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("010000101010011", sup_str_buf);

  ni = 2; // 3
  res = g.PPCExtension(i0, (sup & t.NthData(ni)), g.CoreIndex(i0), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i0p3, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("010001001010011", sup_str_buf);

  ni = 3; // 4
  res = g.PPCExtension(i0, (sup & t.NthData(ni)), g.CoreIndex(i0), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i0p4, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("001110101000100", sup_str_buf);

  ba a1p2  = {{0,1,2,END}}; SortedItemSet i1p2(a1p2); // ppc-ext of (i1 \cup 1)
  ba a1p3  = {{0,2,END}};   SortedItemSet i1p3(a1p3); // ppc-ext of (i1 \cup 2)
  ba a1p4  = {{0,3,END}};   SortedItemSet i1p4(a1p4); // ppc-ext of (i1 \cup 3)

  g.Support(i1, &sup);
  std::cout << "sup=" << sup << std::endl;

  ni = 1; // 2
  res = g.PPCExtension(i1, (sup & t.NthData(ni)), g.CoreIndex(i1), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i1p2, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("010000001010011", sup_str_buf);

  ni = 2; // 3
  res = g.PPCExtension(i1, (sup & t.NthData(ni)), g.CoreIndex(i1), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_FALSE(res); // core_i condition not satisfied
  // EXPECT_TRUE(res);
  // EXPECT_EQ(i1p3, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  // boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  // EXPECT_EQ("010000001010011", sup_str_buf);

  ni = 3; // 4
  res = g.PPCExtension(i1, (sup & t.NthData(ni)), g.CoreIndex(i1), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i1p4, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("000000001000100", sup_str_buf);
}

TEST (LCMTest, CoreIndexTest) {
  Table * t;
  FunctionsSuper * functions;
  ReadSampleTestData(&t, &functions);
  LCM_Graph g(*t);

  // from fig.2 in [uno2004a]    1 2 3 4 5 6 7 8 9 10      core_i
  SortedItemSet i0; // -1 (originally 0 but changed to -1, by definition)
  EXPECT_EQ(g.CoreIndex(i0),   -1);

  // from fig.2 in [uno2004a]    1 2 3 4 5 6 7 8 9 10      core_i
  // vi d1 = {1,4,6};
  // SortedItemSet i1(d1);
  // int a1[] = {1,4,6};
  // SortedItemSet i1();
  ba a1   = {{1,4,6,END}};       SortedItemSet i1(a1);     // 1 (typon of fig.2)
  ba a12  = {{1,2,4,5,6,END}};   SortedItemSet i12(a12);   // 2 
  ba a123 = {{1,2,3,4,5,6,END}}; SortedItemSet i123(a123); // 3
  EXPECT_EQ(g.CoreIndex(i1),   1);
  EXPECT_EQ(g.CoreIndex(i12),  2);
  EXPECT_EQ(g.CoreIndex(i123), 3);

  ba a2   = {{2,END}};           SortedItemSet i2(a2);     // 2
  ba a23  = {{2,3,5,END}};       SortedItemSet i23(a23);   // 3 
  ba a24  = {{2,4,END}};         SortedItemSet i24(a24);   // 4
  ba a25  = {{2,5,END}};         SortedItemSet i25(a25);   // 5
  EXPECT_EQ(g.CoreIndex(i2),   2);
  EXPECT_EQ(g.CoreIndex(i23),  3);
  EXPECT_EQ(g.CoreIndex(i24),  4);
  EXPECT_EQ(g.CoreIndex(i25),  5);

  delete t;
  delete functions;
}

TEST (LCMTest, VBACoreIndexTest) {
  std::stringstream item_ss, pos_ss;

  item_ss << "#test,i0,i1,i2,i3,i4,i5,i6\n"
      //         0 1 2 3 4 5 6 (0 is dummy)
          << "t1,0,1,1,1,1,1,1\n" // 1 2 3 4 5 6
          << "t2,0,0,1,1,0,1,0\n" // 2 3 5
          << "t3,0,0,1,0,0,1,0\n" // 2 5
          << "t4,0,1,1,0,1,1,1\n" // 1 2 4 5 6
          << "t5,0,0,1,0,1,0,0\n" // 2 4
          << "t6,0,1,0,0,1,0,1\n" // 1 4 6
          << "t7,0,0,0,1,1,0,1\n";// 3 4 6

  // random dummy data
  pos_ss  << "#test,expression\n"
          << "t1,1\n"
          << "t2,1\n"
          << "t3,0\n"
          << "t4,0\n"
          << "t5,0\n"
          << "t6,0\n"
          << "t7,0\n";

  Functions4fisher functions(1);
  TableVBA<uint64> t(item_ss, pos_ss, functions);
  LCM_Graph_VBA<uint64> g(t);

  // from fig.2 in [uno2004a]    1 2 3 4 5 6 7 8 9 10      core_i
  SortedItemSet i0; // -1 (originally 0 but changed to -1, by definition)
  EXPECT_EQ(g.CoreIndex(i0),   -1);

  // from fig.2 in [uno2004a]    1 2 3 4 5 6 7 8 9 10      core_i
  // vi d1 = {1,4,6};
  // SortedItemSet i1(d1);
  // int a1[] = {1,4,6};
  // SortedItemSet i1();
  ba a1   = {{1,4,6,END}};       SortedItemSet i1(a1);     // 1 (typon of fig.2)
  ba a12  = {{1,2,4,5,6,END}};   SortedItemSet i12(a12);   // 2
  ba a123 = {{1,2,3,4,5,6,END}}; SortedItemSet i123(a123); // 3
  EXPECT_EQ(g.CoreIndex(i1),   1);
  EXPECT_EQ(g.CoreIndex(i12),  2);
  EXPECT_EQ(g.CoreIndex(i123), 3);

  ba a2   = {{2,END}};           SortedItemSet i2(a2);     // 2
  ba a23  = {{2,3,5,END}};       SortedItemSet i23(a23);   // 3
  ba a24  = {{2,4,END}};         SortedItemSet i24(a24);   // 4
  ba a25  = {{2,5,END}};         SortedItemSet i25(a25);   // 5
  EXPECT_EQ(g.CoreIndex(i2),   2);
  EXPECT_EQ(g.CoreIndex(i23),  3);
  EXPECT_EQ(g.CoreIndex(i24),  4);
  EXPECT_EQ(g.CoreIndex(i25),  5);
}

#if 0
  // from fig.2 in [uno2004a]    1 2 3 4 5 6 7 8 9 10      core_i
  SortedItemSet i3   =         {{3,D,D,D,D,D,D,D,D,D}}; // 3
  SortedItemSet i34  =         {{3,4,6,D,D,D,D,D,D,D}}; // 4
  FillItemSet(i3); FillItemSet(i34);
  EXPECT_EQ(g.CoreIndex(i3),   3);
  EXPECT_EQ(g.CoreIndex(i34),  4);

  // from fig.2 in [uno2004a]    1 2 3 4 5 6 7 8 9 10      core_i
  SortedItemSet i4   =         {{4,D,D,D,D,D,D,D,D,D}}; // 4
  SortedItemSet i46  =         {{4,6,D,D,D,D,D,D,D,D}}; // 6
  SortedItemSet i6   =         {{6,D,D,D,D,D,D,D,D,D}}; // 6
  FillItemSet(i4); FillItemSet(i46); FillItemSet(i6);
  EXPECT_EQ(g.CoreIndex(i4),   4);
  EXPECT_EQ(g.CoreIndex(i46),  6);
  EXPECT_EQ(g.CoreIndex(i6),   6);

  delete t;
}

TEST (LCMTest, ClosureTest) {
  Table * t;
  FunctionsSuper * functions;
  ReadSampleTestData(&t, &functions);
  LCM_Graph g(*t);

  SortedItemSet buf   = {{D,D,D,D,D,D,D,D,D,D}}; FillItemSet(buf);

  SortedItemSet i1    = {{1,D,D,D,D,D,D,D,D,D}}; FillItemSet(i1);
  SortedItemSet i1c   = {{1,4,6,D,D,D,D,D,D,D}}; FillItemSet(i1c);// closure
  g.Closure(i1, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i1c);

  SortedItemSet i12   = {{1,2,4,6,D,D,D,D,D,D}}; FillItemSet(i12);
  SortedItemSet i12c  = {{1,2,4,5,6,D,D,D,D,D}}; FillItemSet(i12c); // closure
  g.Closure(i12, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i12c);

  SortedItemSet i123  = {{1,2,3,4,6,D,D,D,D,D}}; FillItemSet(i123);
  SortedItemSet i123c = {{1,2,3,4,5,6,D,D,D,D}}; FillItemSet(i123c); // closure
  g.Closure(i123, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i123c);

  SortedItemSet i2    = {{2,D,D,D,D,D,D,D,D,D}}; FillItemSet(i2);
  SortedItemSet i2c   = {{2,D,D,D,D,D,D,D,D,D}}; FillItemSet(i2c); // closure
  g.Closure(i2, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i2c);

  SortedItemSet i23   = {{2,3,D,D,D,D,D,D,D,D}}; FillItemSet(i23);
  SortedItemSet i23c  = {{2,3,5,D,D,D,D,D,D,D}}; FillItemSet(i23c);// closure
  g.Closure(i23, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i23c);

  SortedItemSet i24   = {{2,4,D,D,D,D,D,D,D,D}}; FillItemSet(i24);
  SortedItemSet i24c  = {{2,4,D,D,D,D,D,D,D,D}}; FillItemSet(i24c); // closure
  g.Closure(i24, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i24c);

  SortedItemSet i25   = {{2,5,D,D,D,D,D,D,D,D}}; FillItemSet(i25);
  SortedItemSet i25c  = {{2,5,D,D,D,D,D,D,D,D}}; FillItemSet(i25c);// closure
  g.Closure(i25, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i25c);

  SortedItemSet i3    = {{3,D,D,D,D,D,D,D,D,D}}; FillItemSet(i3);
  SortedItemSet i3c   = {{3,D,D,D,D,D,D,D,D,D}}; FillItemSet(i3c);// closure
  g.Closure(i3, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i3c);

  SortedItemSet i34   = {{3,4,D,D,D,D,D,D,D,D}}; FillItemSet(i34);
  SortedItemSet i34c  = {{3,4,6,D,D,D,D,D,D,D}}; FillItemSet(i34c); // closure
  g.Closure(i34, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i34c);

  SortedItemSet i4    = {{4,D,D,D,D,D,D,D,D,D}}; FillItemSet(i4);
  SortedItemSet i4c   = {{4,D,D,D,D,D,D,D,D,D}}; FillItemSet(i4c); // closure
  g.Closure(i4, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i4c);

  SortedItemSet i46   = {{4,6,D,D,D,D,D,D,D,D}}; FillItemSet(i46);
  SortedItemSet i46c  = {{4,6,D,D,D,D,D,D,D,D}}; FillItemSet(i46c); // closure
  g.Closure(i46, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i46c);

  SortedItemSet i5    = {{5,D,D,D,D,D,D,D,D,D}}; FillItemSet(i5);
  SortedItemSet i5c   = {{2,5,D,D,D,D,D,D,D,D}}; FillItemSet(i5c); // closure
  g.Closure(i5, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i5c);

  SortedItemSet i6    = {{6,D,D,D,D,D,D,D,D,D}}; FillItemSet(i6);
  SortedItemSet i6c   = {{4,6,D,D,D,D,D,D,D,D}}; FillItemSet(i6c); // closure
  g.Closure(i6, &buf);
  PrintItemSet(std::cout, buf);
  EXPECT_EQ(buf, i6c);

  delete t;
  delete functions;
}

TEST (LCMTest, PPCExtentionTest) {
  Table * t;
  FunctionsSuper * functions;
  ReadSampleTestData(&t, &functions);
  LCM_Graph g(*t);

  SortedItemSet buf     = {{D,D,D,D,D,D,D,D,D,D}}; FillItemSet(buf);
  bool res=false;
  boost::dynamic_bitset<> sup(t->NuTransaction());
  //boost::dynamic_bitset<> sup_buf(t->NuTransaction());
  std::string sup_str_buf;
  int ni = 0;

  SortedItemSet i0    = {{D,D,D,D,D,D,D,D,D,D}};
  SortedItemSet i0p1  = {{1,4,6,D,D,D,D,D,D,D}}; // ppc-ext of (i0 \cup 1)
  SortedItemSet i0p2  = {{2,D,D,D,D,D,D,D,D,D}}; // ppc-ext of (i0 \cup 2)
  SortedItemSet i0p3  = {{3,D,D,D,D,D,D,D,D,D}}; // ppc-ext of (i0 \cup 3)
  SortedItemSet i0p4  = {{4,D,D,D,D,D,D,D,D,D}}; // ppc-ext of (i0 \cup 4)
  FillItemSet(i0); FillItemSet(i0p1); FillItemSet(i0p2); FillItemSet(i0p3); FillItemSet(i0p4);
  // the following is a type of fig.2 [uno2004a]
  // SortedItemSet i0p6  = {{6,D,D,D,D,D,D,D,D,D}}; // ppc-ext of (i0 \cup 6)
  g.Support(i0, &sup);
  std::cout << "sup=" << sup << std::endl;

  ni = 1;
  res = g.PPCExtension(i0, (sup & t->NthData(ni)), g.CoreIndex(i0), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i0p1, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("0101001", sup_str_buf);

  ni = 2;
  res = g.PPCExtension(i0, (sup & t->NthData(ni)), g.CoreIndex(i0), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i0p2, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("0011111", sup_str_buf);

  ni = 3;
  res = g.PPCExtension(i0, (sup & t->NthData(ni)), g.CoreIndex(i0), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i0p3, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("1000011", sup_str_buf);

  ni = 4;
  res = g.PPCExtension(i0, (sup & t->NthData(ni)), g.CoreIndex(i0), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i0p4, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("1111001", sup_str_buf);

  ni = 6;
  res = g.PPCExtension(i0, (sup & t->NthData(ni)), g.CoreIndex(i0), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_FALSE(res);

  SortedItemSet i1    = {{1,4,6,D,D,D,D,D,D,D}};
  SortedItemSet i1p2  = {{1,2,4,5,6,D,D,D,D,D}}; // ppc-ext of (i1 \cup 2)
  FillItemSet(i1); FillItemSet(i1p2);
  g.Support(i1, &sup);
  std::cout << "sup=" << sup << std::endl;

  ni = 2;
  res = g.PPCExtension(i1, (sup & t->NthData(ni)), g.CoreIndex(i1), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i1p2, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("0001001", sup_str_buf);

  ni = 3;
  res = g.PPCExtension(i1, (sup & t->NthData(ni)), g.CoreIndex(i1), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_FALSE(res);

  ni = 5;
  res = g.PPCExtension(i1, (sup & t->NthData(ni)), g.CoreIndex(i1), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_FALSE(res);

  SortedItemSet i12   = {{1,2,4,5,6,D,D,D,D,D}};
  SortedItemSet i12p3 = {{1,2,3,4,5,6,D,D,D,D}}; // ppc-ext of (i12 \cup 3)
  FillItemSet(i12); FillItemSet(i12p3);
  std::cout << "sup=" << sup << std::endl;

  ni = 3;
  res = g.PPCExtension(i12, (sup & t->NthData(ni)), g.CoreIndex(i12), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i12p3, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("0000001", sup_str_buf);

  SortedItemSet i2    = {{2,D,D,D,D,D,D,D,D,D}};
  SortedItemSet i2p3  = {{2,3,5,D,D,D,D,D,D,D}}; // ppc-ext of (i2 \cup 3)
  SortedItemSet i2p4  = {{2,4,D,D,D,D,D,D,D,D}}; // ppc-ext of (i2 \cup 4)
  SortedItemSet i2p5  = {{2,5,D,D,D,D,D,D,D,D}}; // ppc-ext of (i2 \cup 5)
  FillItemSet(i2); FillItemSet(i2p3); FillItemSet(i2p4); FillItemSet(i2p5);
  g.Support(i2, &sup);
  std::cout << "sup=" << sup << std::endl;

  ni = 3;
  res = g.PPCExtension(i2, (sup & t->NthData(ni)), g.CoreIndex(i2), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i2p3, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("0000011", sup_str_buf);

  ni = 4;
  res = g.PPCExtension(i2, (sup & t->NthData(ni)), g.CoreIndex(i2), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i2p4, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("0011001", sup_str_buf);

  ni = 5;
  res = g.PPCExtension(i2, (sup & t->NthData(ni)), g.CoreIndex(i2), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i2p5, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("0001111", sup_str_buf);

  ni = 6;
  res = g.PPCExtension(i2, (sup & t->NthData(ni)), g.CoreIndex(i2), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_FALSE(res);

  SortedItemSet i3    = {{3,D,D,D,D,D,D,D,D,D}};
  SortedItemSet i3p4  = {{3,4,6,D,D,D,D,D,D,D}}; // ppc-ext of (i3 \cup 4)
  FillItemSet(i3); FillItemSet(i3p4);
  g.Support(i3, &sup);
  std::cout << "sup=" << sup << std::endl;

  ni = 4;
  res = g.PPCExtension(i3, (sup & t->NthData(ni)), g.CoreIndex(i3), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i3p4, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("1000001", sup_str_buf);

  ni = 5;
  res = g.PPCExtension(i3, (sup & t->NthData(ni)), g.CoreIndex(i3), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_FALSE(res);

  ni = 6;
  res = g.PPCExtension(i3, (sup & t->NthData(ni)), g.CoreIndex(i3), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_FALSE(res);

  SortedItemSet i4    = {{4,D,D,D,D,D,D,D,D,D}};
  SortedItemSet i4p6  = {{4,6,D,D,D,D,D,D,D,D}}; // ppc-ext of (i4 \cup 6)
  FillItemSet(i4); FillItemSet(i4p6);
  g.Support(i4, &sup);
  std::cout << "sup=" << sup << std::endl;

  ni = 5;
  res = g.PPCExtension(i4, (sup & t->NthData(ni)), g.CoreIndex(i4), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_FALSE(res);

  ni = 6;
  res = g.PPCExtension(i4, (sup & t->NthData(ni)), g.CoreIndex(i4), ni, &buf);
  if (res) PrintItemSet(std::cout, buf);
  EXPECT_TRUE(res);
  EXPECT_EQ(i4p6, buf);
  // boost::to_string(sup_buf, sup_str_buf);
  boost::to_string((sup & t.NthData(ni)), sup_str_buf);
  EXPECT_EQ("1101001", sup_str_buf);

  delete t;
  delete functions;
}
#endif

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
