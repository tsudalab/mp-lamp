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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <cmath>

#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
//#include <boost/range/algorithm/for_each.hpp>

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#include "gflags/gflags.h"

#include "table.h"
#include "graph.h"

namespace lamp_search {

Graph::Node::Node() :
    depth_ (0),
    visits_ (0ll),
    pmin_ (0.0), // is this OK?
    pval_ (1.0) // is this OK?
{
  // init as the root node
  for (std::size_t i=0 ; i<kMaxItemSet ; i++) {
    sorted_items_[i] = std::numeric_limits<int>::max();
  }
}

void Graph::Node::Init() {
  // init as the root node
  for (std::size_t i=0 ; i<kMaxItemSet ; i++) {
    sorted_items_[i] = std::numeric_limits<int>::max();
  }
  depth_ = 0;
  visits_ = 0;
  pmin_ = 0.0; // is this OK?
  pval_ = 1.0; // is this OK?
}

//==============================================================================

Graph::Graph(const Table & table) :
    t_ (table)
{
}

Graph::~Graph() {
}

} // namespace lamp_search
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
