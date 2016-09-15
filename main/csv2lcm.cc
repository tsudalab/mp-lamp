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

// analyze problem by trying reduction
#include <iostream>
#include <iomanip>
#include <fstream>
#include "gflags/gflags.h"

#include <sys/time.h>
#include <sys/resource.h>

#include "timer.h"
#include "lcm_graph_vba.h"
#include "lcm_dfs_vba.h"

#include "database.h"
#include "lamp_graph.h"
#include "lamp.h"

DEFINE_string(item, "", "filename of item set");

using namespace lamp_search;

int main(int argc, char ** argv)
{
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_item == "") {
    std::cout << "specify item file with --item" << std::endl;
    return 1;
  }

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
  ifs1.open(FLAGS_item.c_str(), std::ios::in);

  reader.ReadFiles(&bsh,
                   ifs1, &data, &nu_trans, &nu_items,
                   item_names, transaction_names, &max_item_in_transaction);

  ifs1.close();

  reader.PrintLCM(std::cout, nu_trans, nu_items, bsh, data);

  delete item_names;
  delete transaction_names;

  return 0;
}

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
