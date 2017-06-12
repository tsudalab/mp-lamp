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

using namespace lamp_search;

TEST (DatabaseLCMTest, LCMOutTest1) {
  VariableBitsetHelper<uint64> * bsh = NULL;

  uint64 * data = NULL;
  boost::array<int, 3> counters; // nu_bits, nu_items, max_item_in_transaction
  counters.assign(-1);

  int nu_trans;
  int nu_items;
  int max_item_in_transaction;

  std::vector< std::string > * item_names = NULL;
  std::vector< std::string > * transaction_names = NULL;
  item_names = new std::vector< std::string >;
  transaction_names = new std::vector< std::string >;

  DatabaseReader<uint64> reader;

  std::ifstream ifs1;
  ifs1.open("../samples/sample_data/sample_item.csv", std::ios::in);
  reader.ReadFiles(&bsh,
                   ifs1, &data, &nu_trans, &nu_items,
                   item_names, transaction_names, &max_item_in_transaction);
  ifs1.close();

  std::stringstream s;
  reader.PrintLCM(s, nu_trans, nu_items, bsh, data);
  std::cout << s.str() << std::flush;

  std::string lcm_str =
      "1 2 3 \n"
      "1 2 3 \n"
      "1 4 \n"
      "\n"
      "1 2 3 \n"
      "1 \n"
      "1 2 3 4 \n"
      "\n"
      "2 4 \n"
      "3 \n"
      "4 \n"
      "4 \n"
      "4 \n"
      "1 2 3 \n"
      "\n"
      ;
  EXPECT_EQ(lcm_str, s.str());

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

  delete item_names;
  delete transaction_names;
}

TEST (DatabaseLCMTest, LCMReadTest1) {
  VariableBitsetHelper<uint64> * bsh = NULL;

  uint64 * data = NULL;
  boost::array<int, 3> counters; // nu_bits, nu_items, max_item_in_transaction
  counters.assign(-1);

  int nu_trans;
  int nu_items;
  int max_item_in_transaction;

  std::vector< std::string > * item_names = NULL;
  std::vector< std::string > * transaction_names = NULL;
  item_names = new std::vector< std::string >;
  transaction_names = new std::vector< std::string >;

  DatabaseReader<uint64> reader;

  std::string lcm_str =
      "1 2 3 \n"
      "1 2 3 \n"
      "1 4 \n"
      "\n"
      "1 2 3 \n"
      "1 \n"
      "1 2 3 4 \n"
      "\n"
      "2 4 \n"
      "3 \n"
      "4 \n"
      "4 \n"
      "4 \n"
      "1 2 3 \n"
      "\n"
      ;

  std::stringstream lcm_ss;
  lcm_ss << lcm_str;

  reader.ReadItemsLCM(lcm_ss, &nu_trans, &nu_items,
                      &bsh, &data, item_names, &max_item_in_transaction);

  EXPECT_EQ(4, nu_items);
  EXPECT_EQ(4, max_item_in_transaction);
  EXPECT_EQ(15, nu_trans);

  std::stringstream s1;
  reader.PrintLCM(s1, nu_trans, nu_items, bsh, data);
  std::cout << s1.str() << std::flush;

  EXPECT_EQ(lcm_str, s1.str());

  std::stringstream s2;
  reader.PrintCSV(s2, nu_trans, nu_items, bsh, item_names, data);
  std::cout << s2.str() << std::flush;

  std::string csv_str =
      "#i,1,2,3,4\n"
      "0,1,1,1,0\n"
      "1,1,1,1,0\n"
      "2,1,0,0,1\n"
      "3,0,0,0,0\n"
      "4,1,1,1,0\n"
      "5,1,0,0,0\n"
      "6,1,1,1,1\n"
      "7,0,0,0,0\n"
      "8,0,1,0,1\n"
      "9,0,0,1,0\n"
      "10,0,0,0,1\n"
      "11,0,0,0,1\n"
      "12,0,0,0,1\n"
      "13,1,1,1,0\n"
      "14,0,0,0,0\n"
      ;

  EXPECT_EQ(csv_str, s2.str());

  std::string item_str =
      "010000001110111\n"
      "010000101010011\n"
      "010001001010011\n"
      "001110101000100\n"
      ;

  // str1 =
  //     "100000111111\n"
  //     "100001101011\n"
  //     "100010101011\n"
  //     "011101100100\n"
  //     ;
  // str1 =
  //     "010000001110111\n"
  //     "010000101010011\n"
  //     "010001001010011\n"
  //     "001110101000100\n"
  //      x      x   x
  //     ;
  // std::str2 =
  //     "110000101011\n"
      // ;

  delete item_names;
  delete transaction_names;
}

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
