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

#ifndef _LAMP_SEARCH_TOPK_DATA_H_
#define _LAMP_SEARCH_TOPK_DATA_H_

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

namespace lamp_search {

class TopKData {
 public:
  friend class DFS;
  friend class LCM_DFS;

  typedef std::pair<double, std::size_t> val_index;
  typedef std::vector< val_index >::const_iterator const_iterator;

  TopKData(const std::vector< ItemInfo > & info, int topk_num) :
      item_info_(info), k_(topk_num) {}
  // TopKData(const Table & t, int topk_num) :
  //     t_(t), k_(topk_num) {}
  // TopKData(const Graph & g, const Table & t, int topk_num) :
  //     g_(g), t_(t), k_(topk_num) {}

  void Init();

  void ResetK(int k);

  // void Insert(double pval, const boost::dynamic_bitset<> & item_bitset);
  void Insert(double pval, const SortedItemSet & is);

  double BestPVal() const { return pval_index_pair_[0].first; }
  double KthPVal() const { return pval_index_pair_[k_-1].first; }
  double KthIndex() const {
    return pval_index_pair_[k_-1].second;
  }

  const val_index & Nth(std::size_t n) const { return pval_index_pair_[n]; }
  const SortedItemSet & NthItemset(std::size_t n) const {
    return itemset_[Nth(n).second];
  }
  // const boost::dynamic_bitset<> & NthItemset(std::size_t n) const {
  //   return itemset_[Nth(n).second];
  // }

  //const Graph & GetGraph() const { return g_; }
  // const Table & GetTable() const { return t_; }

  friend std::ostream& operator<<(std::ostream & out, const TopKData & topk);

 private:
  //const Graph & g_;
  // const Table & t_;
  const std::vector< ItemInfo > & item_info_;

  std::size_t k_;

  std::vector< val_index > pval_index_pair_;
  // prune if (pval > worst in sorted_bestK_pval_)
  // tok k itemset sorted with the same order as sorted_bestK_pval_

  std::vector<SortedItemSet> itemset_;
  // todo: replace with SortedItemSet
  //std::vector< boost::dynamic_bitset<> > itemset_; // sorted

  void Sort();
  static bool less(const val_index & lhs,
                   const val_index & rhs) {
    if (lhs.first < rhs.first) return true;
    if (lhs.first == rhs.first) return (lhs.second < rhs.second);
    return false;
  }
};

//==============================================================================

// void TopKData::Insert(double pval, const boost::dynamic_bitset<> & item_bitset) {
//   if (pval > KthPVal()) return;

//   bool found = false;
//   for (std::size_t ki=0; ki < k_ ; ki++) {
//     //if ( pval_index_pair_[ki].second != pval ) continue;
    
//     if ( itemset_[pval_index_pair_[ki].second] == item_bitset) found = true;
//   }
//   if ( found ) return;

//   itemset_[KthIndex()] = item_bitset;
//   pval_index_pair_[k_-1] = std::make_pair(pval, KthIndex());

//   Sort();
// }

inline
void TopKData::Insert(double pval, const SortedItemSet & is) {
  if (pval > KthPVal()) return;

  bool found = false;
  for (std::size_t ki=0; ki < k_ ; ki++) {
    if (itemset_[pval_index_pair_[ki].second] == is) found = true;
  }
  if ( found ) return;

  itemset_[KthIndex()] = is;
  pval_index_pair_[k_-1] = std::make_pair(pval, KthIndex());

  Sort();
}

inline
void TopKData::Sort() {
  std::sort(pval_index_pair_.begin(), pval_index_pair_.end(), less);
}

inline
void TopKData::ResetK(int k) {
  k_ = k;
  Init();
}

//==============================================================================

std::ostream & operator<<(std::ostream & out, const TopKData & topk);

} // namespace lamp_search

#endif // _LAMP_SEARCH_TOPK_DATA_H_
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
