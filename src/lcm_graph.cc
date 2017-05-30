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
#include <algorithm>

#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
//#include <boost/range/algorithm/for_each.hpp>

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#include "gflags/gflags.h"

#include "sorted_itemset.h"
#include "table.h"
#include "lcm_graph.h"

DECLARE_double(a); // significance level alpha

namespace {

boost::dynamic_bitset<> support_pre; // buffer for prev support
boost::dynamic_bitset<> support_buf; // buffer for current support

boost::dynamic_bitset<> closure_buf; // buffer for closure calculation

} // anonymous namespace

namespace lamp_search {

LCM_Graph::LCM_Graph(const Table & table) :
    t_ (table) {
  support_pre.resize(t_.NuTransaction(), true); // initialize buffer
  support_buf.resize(t_.NuTransaction(), true); // initialize buffer

  closure_buf.resize(t_.NuTransaction(), true); // initialize buffer
}

LCM_Graph::~LCM_Graph() {
}

int LCM_Graph::CoreIndex(const SortedItemSet & items) const {
  if (items.Empty()) return -1;
  //if (items[0] == std::numeric_limits<int>::max()) return 0;

  support_buf.set(); // set all to true

  // starting from item[0] instead of 0
  // for exceptional case where all bits in table are 1
  int core_i = items[0];
  support_pre = t_.NthData(items[0]);

  int i=0;
  for(SortedItemSet::const_iterator it=items.begin(); it !=items.end(); ++it) {
    if (i >= t_.NuItems() ) break;
    i++;

    support_buf &= t_.NthData(*it); // support for *it item

    // if support didn't change don't udpate core_i
    if (support_buf == support_pre) continue;
    else {
      support_pre = support_buf;
      core_i = *it;
      continue;
    }
  }

  return core_i;
}

void LCM_Graph::Closure(const SortedItemSet & items, SortedItemSet * clo) {
  for (int i=0 ; i<kMaxItemSet ; i++) (*clo)[i] = items[i]; // copy

  // calculate support of items
  int clo_last=0;
  support_buf.set(); // set all to true
  for (int i=0 ; i<kMaxItemSet ; i++) {
    if (i >= t_.NuItems() ) break;
    if (items[i] > t_.NuTransaction()) break; // std::numeric_limits<int>::max()
    support_buf &= t_.NthData(items[i]); // support for item[i]
    clo_last++;
  }  

  int ii = 0;
  for (int i=0 ; i<t_.NuItems() ; i++) {
    if (items[ii] == i) { ii++; continue; }

    closure_buf = t_.NthData(i);
    closure_buf &= support_buf;
    closure_buf ^= support_buf;
    if ( closure_buf.none() ) (*clo)[clo_last++] = i;
  }

  std::sort(clo->begin(), clo->end());
  return;
}

/** return if not ppc extention
 *  asserts new_item > core_i(items) and new_item not in items
 *  and sup already &= ed with support of new_item */
bool LCM_Graph::PPCExtension(const SortedItemSet & items,
                             const boost::dynamic_bitset<> & sup,
                             int core_i, // not needed?
                             int new_item,
                             SortedItemSet * ext) const {
  // todo: reuse this buffer to avoid redundant copy
  if ( items.Full() ) return false;
  (*ext) = items; // copy

  // checking items(i-1) and ext(i-1)
  // e.g. checking (i-1) prefix of items is preserved
  int ii = 0;
  for (int i=0 ; i<new_item ; i++) {
    // skip existing item
    if (items.Size() > ii && items[ii] == i) { ii++; continue; }

    closure_buf = t_.NthData(i);
    closure_buf &= sup;
    closure_buf ^= sup;
    // if closure_buf.none(), sup is subset of t_.NthData(i)
    if ( closure_buf.none() ) return false;
    // means cond (iii) not satisfied [uno et al. 2004a] sec. 4.3
  }

  (*ext).Push(new_item);

  for (int i=new_item+1 ; i<t_.NuItems() ; i++) { // calculating closure
    // skip existing item
    if (items.Size() > ii && items[ii] == i) { ii++; continue; }

    closure_buf = t_.NthData(i);
    closure_buf &= sup;
    closure_buf ^= sup;
    // if closure_buf.none(), sup is subset of t_.NthData(i)
    if ( closure_buf.none() ) {
      if ( !ext->Full() ) (*ext).Push(i); // just skip pushing
    }
  }

  return true;
}

void LCM_Graph::Support(const SortedItemSet & items,
                        boost::dynamic_bitset<> * bs) const {
  bs->set(); // set all to true
  int i=0;
  for(SortedItemSet::const_iterator it=items.begin(); it !=items.end(); ++it) {
    if (i >= t_.NuItems() ) break;
    if (*it > t_.NuTransaction()) break; // std::numeric_limits<int>::max()
    *bs &= t_.NthData(*it); // support for *it
    i++;
  }
}

} // namespace lamp_search

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
