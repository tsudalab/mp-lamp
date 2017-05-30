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

#ifndef _LAMP_SEARCH_LCM_DFS_VBA_H_
#define _LAMP_SEARCH_LCM_DFS_VBA_H_

#include <vector>
#include <algorithm>

#include <boost/array.hpp>
#include <map>

#include "utils.h"
#include "sorted_itemset.h"
#include "topk.h"
#include "variable_bitset_array.h"
#include "table_vba.h"
#include "timer.h"

#include "random.h"
#include "lcm_graph_vba.h"
#include "lcm_dfs.h" // needed?

namespace lamp_search {

// used only for sequential version
// support will be reconstructed on each job for parallel version
class SupportStack {
 public:
  SupportStack(const VariableBitsetHelper<uint64> & helper)
      : bsh (helper) {
    stack_ = new VariableBitsetArray<uint64>(bsh.nu_bits, kMaxSearchDepth+1);
  }

  ~SupportStack() {
    delete stack_;
  }

  void Init() {
    num_ = 0; // clear
    // push all 1 support
    stack_->Set(Top());
    // bsh.Print(std::cout, Top()); std::cout << std::endl;
  }

  void Inc() { num_++; } // update Top manually
  void Dec() { num_--; }

  uint64 * Top() { return stack_->N(num_); }

  bool Empty() const { return (num_ == 0); }
  bool Full() const { return (num_ >= kMaxSearchDepth+1); }

 private:
  const VariableBitsetHelper<uint64> & bsh; // bitset accessor

  //boost::array<VariableBitset<uint64>, kMaxSearchDepth+1> * stack_;
  VariableBitsetArray<uint64> * stack_;

  int num_;
};

class LCM_DFS_VBA;

/** class for recording search state.
 ** basically, use one State for one thread/process
 ** for parallel version, state will be sent as messages */
class LCM_DFS_VBA_State {
 public:
  // typedef VariableBitsetTraits<uint64> traits;

  void Init();

  int depth_;

  // note: consider replacing it with path of core_i
  // (which is unordered array of item id's)

  // must have a fixed size (e.g. kMaxSearchDepth+1)
  // empty set is pushed at edge_stack_[0] for root state
  // sequence of selected edge (stack of SortedItemSet)
  SortedItemSet itemset_; // current itemset_

  // nu_bits_ not needed because it should be equal to t_.NuTransaction()
  // note: support should be reconstructed each time for parallel version

  // for hash table version
  // on iter, record support (as short int) for all items
};

class LCM_DFS_VBA {
 public:
  LCM_DFS_VBA(const LCM_Graph_VBA<uint64> & g);
  ~LCM_DFS_VBA();

  // --------
  // DFS

  void Search();
  void DFSInit();

  void Iter(LCM_DFS_VBA_State * st, int sup_threshold, double minp_threshold);

  // must call init after this ?
  void SetTopK(int k) { top_k_.ResetK(k); }

  std::ostream & DumpTopK(std::ostream & out) const;

  std::ostream & PrintResults(std::ostream & out) const;

  // --------
  // LAMP

  void LAMP();
  void LAMPInit();

  void LAMPIter(LCM_DFS_VBA_State * st, int sup_threshold,
                double sig_level, bool second_phase);

  std::ostream & PrintLAMPResults(std::ostream & out) const;

  // --------
  // LAMP Single

  void LAMPSingle();
  void LAMPSingleInit();

  void LAMPSingleIter(LCM_DFS_VBA_State * st);

  bool ExceedCsThr() const;
  void IncLambdaThr();

  void DiscardSignificantList();

  void IncCsAccum(int sup_num);

  std::ostream & PrintLAMPSingleResults(std::ostream & out) const;

  // --------

  std::ostream & PrintSignificantList(std::ostream & out) const;
  std::ostream & PrintSignificantMap(std::ostream & out) const;

  const LCM_Graph_VBA<uint64> & GetGraph() const { return g_; }

  class Stats {
   public:
    Stats() : nu_iteration_ (0ll), nu_node_visits_ (0ll) {}
    
    long long int nu_iteration_;
    long long int nu_node_visits_;
  };

  Timer * timer_;
  Stats stats_;
  long long int search_start_time_;

  TopKData & TopK() { return top_k_; }
  const TopKData & TopK() const { return top_k_; }

  static const int k_int_max;
  // assuming (digits in long long int) > (bits of double mantissa)
  static const long long int k_cs_max;

 private:
  const LCM_Graph_VBA<uint64> & g_;
  const TableVBA<uint64> & tbl_;
  const VariableBitsetHelper<uint64> & bsh; // bitset helper

  SupportStack sup_stack_;
  LCM_DFS_VBA_State current_state_;

  // for topk
  TopKData top_k_;

  // for LAMP
  std::vector< std::pair<double, SortedItemSet> > significant_list_;
  static bool less(const std::pair<double, SortedItemSet> & lhs,
                   const std::pair<double, SortedItemSet> & rhs) {
    if (lhs.first < rhs.first) return true;
    return false;
  }

  // variables for single phase LAMP
  // equals to maximum support of single item
  int lambda_max_;
  // lambda_ initially set to 1. will be incremented to N if cs_thr[N] exceeded
  // in 1st phase, search will be pruned if (sup_num < lambda_)
  int lambda_;
  // todo: initial value can be set after checking all single item itemsets
  double sig_level_; // initially set to 1. set to 0.05 (FLAGS_a) / cs_thr_[global_sup_thr_-1]
  double * pmin_thr_; // pmin_thr[sup] == tbl.PMin(sup), maybe redundant
  long long int * cs_thr_; // cs_thr[sup] shows closed set num threshold
  // cs_num_accum_[sup] shows closed set num with support higher than or equals to sup
  long long int * cs_num_accum_;

  std::multimap< double, SortedItemSet > significant_map_;

  // should be equal to sum of cs_num_[global_sup_thr_..global_sup_max_]
  //long long int global_cs_num_;

  // for LAMP A*
  // replace with variable length itemset array/deque
  std::multimap< int, SortedItemSet > open_list_;

  // expand num of current iteration (== number of Iter calls)
  // without positive check, this should be (the number of closed set) + 1
  long long int expand_num_;
  // total expand num
  long long int total_expand_num_;

  long long int closed_set_num_;

  long long int final_closed_set_num_;
  int final_support_;

  // todo: allocate somewhere outside of state class
  uint64 * sup_bitset_buf_;

  // std::vector<ItemInfo> item_info_; // list of item id sorted by pmin
};

} // namespace lamp_search

#endif // _LAMP_SEARCH_LCM_DFS_VBA_H_

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
