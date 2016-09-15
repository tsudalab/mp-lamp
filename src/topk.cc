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

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <list>
#include <limits>

#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>

#include "utils.h"
#include "sorted_itemset.h"
#include "topk.h"

namespace lamp_search {

void TopKData::Init() {
  // const std::vector<Table::ItemInfo> & items = t_.GetItemInfo();

  pval_index_pair_.clear();
  itemset_.clear();
  double epsilon = 0.00001;

  SortedItemSet empty_itemset;
  // for (std::size_t d=0;d<kMaxSearchDepth;d++) {
  //   empty_itemset[d] = std::numeric_limits<int>::max();
  // }

  // boost::dynamic_bitset<> tmp_bitset(t_.NuItems());
  // tmp_bitset.reset();

  itemset_.resize(k_, empty_itemset);
  for (std::size_t ki = 0; ki < k_ ; ki++) {
    pval_index_pair_.push_back( std::make_pair(1.0 + epsilon, ki) );
    // itemset_.push_back(empty_itemset);
    // itemset_.push_back(tmp_bitset);
  }

  SortedItemSet tmp_itemset;
  std::size_t i = 0;
  for (std::vector< ItemInfo >::const_iterator it = item_info_.begin();
       it != item_info_.end(); ++it) {
    if (i > k_) break;
    tmp_itemset.Push((*it).id);
    Insert((*it).pval, tmp_itemset);
    tmp_itemset.Clear();
    // tmp_bitset[(*it).id] = true;
    // Insert((*it).pval, tmp_bitset);
    // tmp_bitset.reset();
    i++;
  }
}

std::ostream & operator<<(std::ostream & out, const TopKData & topk) {
  std::stringstream s;

  for(TopKData::const_iterator it = topk.pval_index_pair_.begin();
      it != topk.pval_index_pair_.end(); ++it) {
    s << "pval=" << std::setw(12) << (*it).first
      << ":set=" << topk.itemset_[(*it).second]
        //<< ":set=" << topk.itemset_[(*it).second]
      << ":idx=" << (*it).second
      << std::endl;
  }
  
  out << s.str() << std::flush;
  return out;
}

} // namespace lamp_search
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
