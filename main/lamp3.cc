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

// LAMP ver.3
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

DEFINE_bool(lcm, false, "item file is lcm style");

DEFINE_string(item, "", "filename of item set");
DEFINE_string(pos, "", "filename of positive / negative file");
DEFINE_int32(posnum, 0, "positive total (used if for 1st phase only)");

DECLARE_bool(second_phase); // true, "do second phase"
DECLARE_bool(third_phase);  // true, "do third phase"

DEFINE_bool(loop, false, "use loop version (similar behavior to MP version)");

DEFINE_bool(log, false, "show log");

using namespace lamp_search;

int main(int argc, char ** argv)
{
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_item == "") {
    std::cout << "specify item file with --item" << std::endl;
    return 1;
  }

  if (FLAGS_pos == "" && FLAGS_posnum == 0) {
    std::cout << "specify --pos or --posnum" << std::endl;
    return 1;
  }

  if (FLAGS_third_phase && FLAGS_pos == "") {
    std::cout << "specify positive file by --pos for third phase" << std::endl;
    return 1;
  }

  // error if at_least_n > at_most_n

  long long int search_start_time, search_end_time;
  Timer::GetInstance()->Start();

  VariableBitsetHelper<uint64> * bsh = NULL;

  uint64 * data = NULL;
  uint64 * positive = NULL;

  int nu_trans;
  int nu_items;
  int nu_pos_total = FLAGS_posnum;
  int max_item_in_transaction;

  std::vector< std::string > * item_names = NULL;
  std::vector< std::string > * transaction_names = NULL;
  item_names = new std::vector< std::string >;
  transaction_names = new std::vector< std::string >;

  DatabaseReader<uint64> reader;

  try
  {
    if (FLAGS_pos != "") {
      std::ifstream ifs1;
      ifs1.open(FLAGS_item.c_str(), std::ios::in);
      if (ifs1.fail())
        throw std::runtime_error(std::string("file not found: ") + FLAGS_item);

      std::ifstream ifs2;
      ifs2.open(FLAGS_pos.c_str(), std::ios::in);
      if (ifs2.fail())
        throw std::runtime_error(std::string("file not found: ") + FLAGS_item);

      if (FLAGS_lcm) {
        reader.ReadFilesLCM(&bsh,
                            ifs1, &data, &nu_trans, &nu_items,
                            ifs2, &positive, &nu_pos_total,
                            item_names, &max_item_in_transaction);
      } else {
        reader.ReadFiles(&bsh,
                         ifs1, &data, &nu_trans, &nu_items,
                         ifs2, &positive, &nu_pos_total,
                         item_names, transaction_names, &max_item_in_transaction);
      }
      ifs2.close();
      ifs1.close();
    }
    else {
      std::ifstream ifs1;
      ifs1.open(FLAGS_item.c_str(), std::ios::in);
      if (ifs1.fail())
        throw std::runtime_error(std::string("file not found: ") + FLAGS_item);

      if (FLAGS_lcm) {
        reader.ReadFilesLCM(&bsh,
                            ifs1, &data, &nu_trans, &nu_items,
                            item_names, &max_item_in_transaction);
      } else {
        reader.ReadFiles(&bsh,
                         ifs1, &data, &nu_trans, &nu_items,
                         item_names, transaction_names, &max_item_in_transaction);
      }
      ifs1.close();
    }
  }
  catch(std::runtime_error & err)
  {
    std::cout << err.what() << std::endl;
    return 1;
  }

  Database<uint64> d(bsh, data, nu_trans, nu_items,
                     positive, nu_pos_total,
                     max_item_in_transaction,
                     item_names, transaction_names);
  LampGraph<uint64> g(d);
  Lamp search(g);

  std::cout << "# ";
  d.ShowInfo(std::cout);

  search_start_time = Timer::GetInstance()->Elapsed();
  try
  {
    if (FLAGS_loop)
      search.SearchLoop();
    else
      search.Search();
    search_end_time = Timer::GetInstance()->Elapsed();
  }
  catch(std::runtime_error & err)
  {
    std::cout << err.what() << std::endl;
    return 1;
  }

  std::cout << "# time all="
            << std::setw(12) << search_end_time / GIGA
            << "\ttime search=" << std::setw(12)
            << (search_end_time - search_start_time) / GIGA
            << std::endl;
  search.PrintResults(std::cout);
  if (FLAGS_log) search.PrintLog(std::cout);

  return 0;
}

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
