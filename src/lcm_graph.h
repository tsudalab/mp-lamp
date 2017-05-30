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

#ifndef _LAMP_SEARCH_LCM_GRAPH_H_
#define _LAMP_SEARCH_LCM_GRAPH_H_

#include <iostream>
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
#include "table.h"

namespace lamp_search {

class LCM_Graph {
 public:

  LCM_Graph(const Table & t);
  ~LCM_Graph();

  const Table & GetTable() const { return t_; }

  // core_i is same as clo_tail [uno2004a, uno2004b]
  int CoreIndex(const SortedItemSet & items) const;

  void Closure(const SortedItemSet & items, SortedItemSet * clo);

  // return if not ppc extention
  // asserts new_item > core_i(items) and new_item not in items
  // and sup already &= ed with support of new_item
  bool PPCExtension(const SortedItemSet & items, const boost::dynamic_bitset<> & sup,
                    int core_i, int new_item, SortedItemSet * ext) const;

  void Support(const SortedItemSet & items, boost::dynamic_bitset<> * sup) const;

 private:
  const Table & t_;

  // SortedItemSet sorted_items_; // default sort based on index
};

} // namespace lamp_search

#endif // _LAMP_SEARCH_LCM_GRAPH_H_
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
