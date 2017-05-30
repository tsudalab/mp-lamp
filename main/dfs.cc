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

// depth first search
#include <iostream>
#include <iomanip>
#include <fstream>
#include "gflags/gflags.h"

#include <sys/time.h>
#include <sys/resource.h>

#include "timer.h"
#include "graph.h"
#include "dfs.h"
#include "functions/Functions4fisher.h"

DEFINE_string(item, "", "filename of item set");
DEFINE_string(pos, "", "filename of positive / negative file");

using namespace lamp_search;

namespace {

// long long int gettime_temp() {
//   long long int musec;
//   struct timeval t;
//   struct timezone tz;
//   gettimeofday(&t, &tz);
//   musec = t.tv_sec * 1000000 + t.tv_usec;
//   return musec;
// }

} // namespace

int main(int argc, char ** argv)
{
  google::ParseCommandLineFlags(&argc, &argv, true);

  //long long int start_time, search_start_time, search_end_time;
  //start_time = gettime_temp();
  long long int search_start_time, search_end_time;
  Timer::GetInstance()->Start();

  std::ifstream ifs1, ifs2;
  ifs1.open(FLAGS_item.c_str(), std::ios::in);
  ifs2.open(FLAGS_pos.c_str(), std::ios::in);
  Functions4fisher functions(1);
  Table t(ifs1, ifs2, functions);
  ifs1.close();
  ifs2.close();
  Graph g(t);
  DFS search(g);

  //search_start_time = gettime_temp();
  search_start_time = Timer::GetInstance()->Elapsed();
  search.Search();

  //search_end_time = gettime_temp();
  search_end_time = Timer::GetInstance()->Elapsed();

  std::cout << "time all="
            << std::setw(12) << search_end_time / GIGA
      //<< std::setw(12) << (search_end_time - start_time) / 1000000.0
            << "\ttime search=" << std::setw(12)
            << (search_end_time - search_start_time) / GIGA
            << std::endl;
  search.PrintResults(std::cout);

  return 0;
}
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
