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

#ifndef _LAMP_SEARCH_LCM_DFS_H_
#define _LAMP_SEARCH_LCM_DFS_H_

#include <vector>
#include <algorithm>

#include "utils.h"
#include "sorted_itemset.h"
#include "topk.h"
#include "lcm_graph.h"

#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>

namespace lamp_search {

class LCM_DFS;

class LCMEdge {
 public:
  LCMEdge() {
    depth_=0;
    items_.Clear();
    for (std::size_t i=0;i<kMaxItemSet;i++) {
      items_[i] = std::numeric_limits<uint64>::max();
    }
  }

  LCMEdge(const SortedItemSet & it) {
    items_ = it;
    // depth_ = it.size();
    // for (int i=0;i<depth_;i++) items_[i] = it[i];
  }

  int depth_;
  SortedItemSet items_;
};

/** class for recording search state.
 ** basically, use one State for one thread/process
 ** for parallel version, state will be sent as messages */
class LCM_DFSState {
 public:
  LCM_DFSState();

  void Init(const LCM_DFS & search);

  // push item to pash and update selected_depth_
  //void Push(int edge_index, int item_id, const boost::dynamic_bitset<> & item_sup);

  // new_support should be equal to current support of all items in path including new_items
  void Push(SortedItemSet new_items, const boost::dynamic_bitset<> & new_support);
  void Pop();

  const std::vector< boost::dynamic_bitset<> > & BSDiffStack() { return bitset_diff_stack_; }

  // const SortedItemSet & ItemSet() const { return itemset_; }
  const SortedItemSet & ItemSet() const { return EdgeTop().items_; }
  boost::dynamic_bitset<> * Support() { return &support_; }
  const LCMEdge & EdgeTop() const { return edge_stack_.back(); }

  int Depth() const { return depth_; }

  // todo: if only used as buffers, move to DFS
  //       or, preapre non-bitset current_item_set
  // todo: change variable name
  boost::dynamic_bitset<> sup_bitset_buf_;
  //boost::dynamic_bitset<> item_bitset_buf_;

  // todo: add excluded_item_bitset / transanction_bitset

 private:
  std::vector<LCMEdge> edge_stack_; // sequence of selected edge (stack of SortedItemSet)
  int depth_;

  std::vector< boost::dynamic_bitset<> > bitset_diff_stack_;
  // SortedItemSet itemset_;
  boost::dynamic_bitset<> support_; // current support
  boost::dynamic_bitset<> support_buf_; // buffer for support calculation

  // on iter, record t and support for all items?
};

class LCM_DFS {
 public:
  LCM_DFS(const LCM_Graph & g);
  ~LCM_DFS();

  void Search();

  void Iter(LCM_DFSState * st, int sup_threshold, double minp_threshold);

  void InitDFS();

  std::ostream & PrintResults(std::ostream & out) const;

  void LAMP();

  /** cs_thr = closed_set_num_ threshold */
  void LAMPIter(LCM_DFSState * st, int sup_threshold,
                double sig_level, bool second_phase);

  void InitLAMP();

  std::ostream & PrintLAMPResults(std::ostream & out) const;

  // must call init after this ?
  void SetTopK(int k) { top_k_.ResetK(k); }

  std::ostream & DumpTopK(std::ostream & out) const;

  const LCM_Graph & GetGraph() const { return g_; }

  class Stats {
   public:
    Stats()
        : nu_iteration_ (0ll), nu_node_visits_ (0ll),
          sum_branch_(0ll), sq_sum_branch_(0ll),
          max_branch_(0ll), min_branch_(0ll)
    {}
    
    long long int nu_iteration_;
    long long int nu_node_visits_;

    long long int sum_branch_; // number of branches for which ppcext was true
    long long int sq_sum_branch_; // square sum of branch num (for variance)
    long long int max_branch_; // max of branches
    long long int min_branch_; // min of branches (maybe not interesting)

  };

  Stats stats_;

  // class ItemInfo {
  //  public:
  //   ItemInfo(const Table::ItemInfo & i) :
  //       id (i.id),
  //       pmin (i.pmin),
  //       pval (i.pval)
  //       //selected_depth_ (-1)
  //   {}
    
  //   int id;
  //   double pmin;
  //   double pval;
  // };

  TopKData & TopK() { return top_k_; }
  const TopKData & TopK() const { return top_k_; }

  // std::vector<ItemInfo> & GetItemInfo() { return item_info_; }
  // const std::vector<ItemInfo> & GetItemInfo() const { return item_info_;}

 private:
  const LCM_Graph & g_;

  TopKData top_k_;
  LCM_DFSState root_state_;

  std::vector< std::pair<double, SortedItemSet> > significant_list_;
  static bool less(const std::pair<double, SortedItemSet> & lhs,
                   const std::pair<double, SortedItemSet> & rhs) {
    if (lhs.first < rhs.first) return true;
    return false;
  }

  // expand num of current iteration (== number of Iter calls)
  // without positive check, this should be (the number of closed set) + 1
  long long int expand_num_;
  // total expand num
  long long int total_expand_num_;

  long long int closed_set_num_;

  long long int final_closed_set_num_;
  int final_support_;

  // std::vector<ItemInfo> item_info_; // list of item id sorted by pmin
};

} // namespace lamp_search

#endif // _LAMP_SEARCH_LCM_DFS_H_
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
