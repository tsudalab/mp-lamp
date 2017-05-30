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
#include "lamp_graph.h"

DECLARE_double(a); // significance level alpha

namespace lamp_search {

template<typename Block>
LampGraph<Block>::LampGraph(const Database<Block> & data) :
    d_ (data), bsh_ ( &(data.VBSHelper()) ) {
  //bsh_ = new VariableBitsetHelper<Block>(d_.NuTransaction());
  support_pre_ = bsh_->New();
  support_buf_ = bsh_->New();
}

template<typename Block>
LampGraph<Block>::~LampGraph() {
  bsh_->Delete(support_pre_);
  bsh_->Delete(support_buf_);
  // delete bsh_;
}

template<typename Block>
int LampGraph<Block>::CoreIndex(const VariableLengthItemsetStack & st,
                                const int * items) const {
  if ( st.GetItemNum(items) == 0 ) return -1;

  bsh_->Set(support_buf_); // set all to 1

  // starting from item[0] instead of 0
  // for exceptional case where all bits in table are 1
  int core_i = st.GetNthItem(items, 0);
  bsh_->Copy(bsh_->N(d_.Data(), core_i), support_pre_);

  for(int i=0 ; i<st.GetItemNum(items) ; i++) {
    int itm = st.GetNthItem(items, i);
    assert(itm < d_.NuItems()); // does this hold?

    bsh_->And(bsh_->N(d_.Data(), itm), support_buf_); // support for i-th item

    if ( !(bsh_->IsEqualTo(support_buf_, support_pre_)) ) {
      bsh_->Copy(support_buf_, support_pre_);
      core_i = itm;
      continue;
    }
  }
  
  return core_i;
}

/** return if not ppc extention
 *  asserts new_item > core_i(items) and new_item not in items
 *  and sup already &= ed with support of new_item */
template<typename Block>
bool LampGraph<Block>::PPCExtension(VariableLengthItemsetStack * st,
                                    const int * items,
                                    Block * sup,
                                    int core_i, // not needed?
                                    int new_item,
                                    int * ext_buf) const {
  st->CopyItem(items, ext_buf);
  // ext_buf should have VariableLengthItemsetStack::kMaxItemSet + 1 size

  // note: do following before call
  // st->PushPre()
  // ext_buf = st->Top()

  // checking items(i-1) and ext(i-1)
  // e.g. checking (i-1) prefix of items is preserved
  int ii = 0;
  for (int i=0 ; i<new_item ; i++) {
    // skip existing item
    if (st->GetItemNum(items) > ii && st->GetNthItem(items, ii) == i) {
      ii++;
      continue;
    }

    // if sup is subset of d_.Data()->N(i), not PPCExtension
    if ( bsh_->IsSubsetOf( sup, bsh_->N(d_.Data(), i) ) ) return false;
    // means cond (iii) not satisfied [uno et al. 2004a] sec. 4.3
  }

  st->PushOneItem(new_item);

  for (int i=new_item+1 ; i<d_.NuItems() ; i++) { // calculating closure
    // skip existing item
    if (st->GetItemNum(items) > ii && st->GetNthItem(items, ii) == i) {
      ii++;
      continue;
    }

    // if closure_buf.none(), sup is subset of t_.Data()->N(i)
    if ( bsh_->IsSubsetOf( sup, bsh_->N(d_.Data(), i) ) ) {
      // todo: maybe success check is needed
      st->PushOneItem(i);
    }
  }

  // note: do following after return
  // st->SetSup(ext, sup_num)
  // st->PushPost(ext)

  return true;
}

template<typename Block>
void LampGraph<Block>::Support(const VariableLengthItemsetStack & st,
                        const int * items, Block * sup) const {
  bsh_->Set(sup); // set all to true
  int n = st.GetItemNum(items);
  for (int i=0;i<n;i++) {
    int item_id = st.GetNthItem(items, i);
    bsh_->And(d_.NthData(item_id), sup);
  }
}

template class LampGraph<uint64>;

} // namespace lamp_search

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
