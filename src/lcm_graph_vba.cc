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

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
//#include <boost/range/algorithm/for_each.hpp>

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#include "gflags/gflags.h"

#include "sorted_itemset.h"
#include "variable_bitset_array.h"
#include "table.h"
#include "lcm_graph_vba.h"

DECLARE_double(a); // significance level alpha

namespace lamp_search {

template<typename Block>
LCM_Graph_VBA<Block>::LCM_Graph_VBA(const TableVBA<Block> & table) :
    t_ (table) {
  bsh = new VariableBitsetHelper<Block>(t_.NuTransaction());
  support_pre = bsh->New();
  support_buf = bsh->New();
}

template<typename Block>
LCM_Graph_VBA<Block>::~LCM_Graph_VBA() {
  bsh->Delete(support_pre);
  bsh->Delete(support_buf);
  delete bsh;
}

template<typename Block>
int LCM_Graph_VBA<Block>::CoreIndex(const SortedItemSet & items) const {
  if (items.Empty()) return -1;
  //if (items[0] == std::numeric_limits<int>::max()) return 0;

  bsh->Set(support_buf); // set all to 1

  // starting from item[0] instead of 0
  // for exceptional case where all bits in table are 1
  int core_i = items[0];
  bsh->Copy(t_.Data()->N(core_i), support_pre);
  // support_pre->CopyFrom(t_.Data()->N(core_i));

  int i=0;
  for(SortedItemSet::const_iterator it=items.begin(); it !=items.end(); ++it) {
    if (i >= t_.NuItems() ) break;
    i++;

    bsh->And(t_.Data()->N(*it), support_buf); // support for *it item

    if ( !(bsh->IsEqualTo(support_buf, support_pre)) ) {
      bsh->Copy(support_buf, support_pre);
      core_i = *it;
      continue;
    }

    // // AndEqualToUpdate1(p1, p2, p3) return true if p1 == (p2&p3)
    // // and overwrites p1 with (p2&p3)
    // if ( !(bsh->AndEqualToUpdate1(support_pre, t_.Data()->N(*it), support_buf)) ) {
    //   // update core_i only if support had changed
    //   core_i = *it;
    //   continue;
    // }
  }
  
  return core_i;
}

/** return if not ppc extention
 *  asserts new_item > core_i(items) and new_item not in items
 *  and sup already &= ed with support of new_item */
template<typename Block>
bool LCM_Graph_VBA<Block>::PPCExtension(const SortedItemSet & items,
                                        Block * sup,
                                        int core_i, // not needed?
                                        int new_item,
                                        SortedItemSet * ext_buf) const {
  // todo: reuse this buffer to avoid redundant copy
  if ( items.Full() ) return false;
  (*ext_buf) = items; // copy

  // checking items(i-1) and ext(i-1)
  // e.g. checking (i-1) prefix of items is preserved
  int ii = 0;
  for (int i=0 ; i<new_item ; i++) {
    // skip existing item
    if (items.Size() > ii && items[ii] == i) {
      ii++;
      continue;
    }

    // if sup is subset of t_.Data()->N(i), not PPCExtension
    if ( t_.Data()->IsSubsetOf(sup, t_.Data()->N(i))) return false;
    // means cond (iii) not satisfied [uno et al. 2004a] sec. 4.3
  }

  (*ext_buf).Push(new_item);

  for (int i=new_item+1 ; i<t_.NuItems() ; i++) { // calculating closure
    // skip existing item
    if (items.Size() > ii && items[ii] == i) {
      ii++;
      continue;
    }

    // if closure_buf.none(), sup is subset of t_.Data()->N(i)
    if ( t_.Data()->IsSubsetOf(sup, t_.Data()->N(i))) {
      if ( !ext_buf->Full() ) (*ext_buf).Push(i); // just skip pushing if Full
    }
  }

  return true;
}

template class LCM_Graph_VBA<uint64>;

} // namespace lamp_search

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
