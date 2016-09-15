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

#ifndef _LAMP_SEARCH_MCTS_H_
#define _LAMP_SEARCH_MCTS_H_

#include <vector>
#include <algorithm>

#include "utils.h"
#include "graph.h"
#include "hash.h"

#include <boost/dynamic_bitset.hpp>

namespace lamp_search {

class TDS {
 public:

 private:
  // life line graph for broadcasting information
};

class MCTS;

class MCTSElement {
 public:
  MCTSElement * getNext() { return next_; }
  void setNext(MCTSElement * p) { next_ = p; }

  void Init();

  hash_type getHashcode() { return hash_code_; }

 private:
  MCTSElement * next_;
  hash_type hash_code_;
  Graph::Node node_;

  double bound_;

  double w_; // accumulative reward
  double s_; // number of simulation
};

/** class for recording search state.
 ** basically, use one State for one thread/process
 ** for parallel version, state will be sent as messages */
class MCTSState {
 public:
  MCTSState();

  hash_type getHashcode() const { return hash_code_; }
  void Init(const MCTS & search, hash_type hash_code);

  // push item to pash and update selected_depth_
  void Push(int edge_index, int item_id, boost::dynamic_bitset<> & sup_bs);
  void Pop();
  void PushN(const std::vector<int> & edge_list);
  void PopN(const std::vector<int> & edge_list);

  const std::vector<int> & ItemStack() const { return item_stack_; }
  int last_edge_index_in_item_stack_;
  const std::vector< boost::dynamic_bitset<> > & BSDiffStack() { return bitset_diff_stack_; }

  boost::dynamic_bitset<> * Support() { return &support_; }

  int Depth() const { return depth_; }

  // todo: if only used as buffers, move to DFS
  //       or, preapre non-bitset current_item_set
  boost::dynamic_bitset<> sup_bitset_buf_;
  boost::dynamic_bitset<> item_bitset_buf_;

 private:
  hash_type hash_code_;
  Graph::Node node_;

  std::vector<int> item_stack_; // sequence of selected items
  std::vector<int> edge_stack_; // sequence of selected edge id
  int depth_;

  std::vector< boost::dynamic_bitset<> > bitset_diff_stack_;
  boost::dynamic_bitset<> support_; // current support

  // on iter, record t and support for all items?
};

/** depth first search class */
class MCTS {
 public:
  MCTS(const Graph & g);
  ~MCTS();
  
  void Search();

  void Iter(MCTSState * st, int sup_threshold, double minp_threshold);

  void Init();

  // must call init after this ?
  void SetTopK(int k) { top_k_.ResetK(k); }

  std::ostream & DumpRandomTable(std::ostream & out) const;

  std::ostream & DumpTopK(std::ostream & out) const;

  std::ostream & PrintResults(std::ostream & out) const;

  const Graph & GetGraph() const { return g_; }

  class TopKData {
   public:
    friend class MCTS;

    typedef std::pair<double, std::size_t> val_index;
    typedef std::vector< val_index >::const_iterator const_iterator;

    TopKData(const Graph & g, const Table & t, int topk_num) :
        g_(g), t_(t), k_(topk_num) {}

    void Init();

    void ResetK(int k);

    void Insert(double pval, const boost::dynamic_bitset<> & item_bitset);

    double BestPVal() const { return pval_index_pair_[0].first; }
    double KthPVal() const { return pval_index_pair_[k_-1].first; }
    double KthIndex() const {
      return pval_index_pair_[k_-1].second;
    }

    const val_index & Nth(std::size_t n) const { return pval_index_pair_[n]; }
    const boost::dynamic_bitset<> & NthItemset(std::size_t n) const {
      return itemset_[Nth(n).second];
    }

    const Graph & GetGraph() const { return g_; }
    const Table & GetTable() const { return t_; }

   private:
    const Graph & g_;
    const Table & t_;

    std::size_t k_;

    std::vector< val_index > pval_index_pair_;
    // prune if (pval > worst in sorted_bestK_pval_)
    // tok k itemset sorted with the same order as sorted_bestK_pval_
    std::vector< boost::dynamic_bitset<> > itemset_; // sorted

    void Sort();
    static bool less(const val_index & lhs,
              const val_index & rhs) {
      if (lhs.first < rhs.first) return true;
      if (lhs.first == rhs.first) return (lhs.second < rhs.second);
      return false;
    }
  };

  TopKData & TopK() { return top_k_; }
  const TopKData & TopK() const { return top_k_; }

  class Stats {
   public:
    Stats() : nu_iteration_ (0ll), nu_node_visits_ (0ll) {}

    long long int nu_iteration_;
    long long int nu_node_visits_;
  };

  Stats stats_;

  class EdgeInfo {
   public:
    EdgeInfo(const ItemInfo & i) :
        id (i.id),
        pmin (i.pmin),
        pval (i.pval)
        //selected_depth_ (-1)
    {}
    
    int id;
    double pmin;
    double pval;
    //int selected_depth_; // not used in brute force version

    //bool Selected() const { return (selected_depth_ >= 0); }
  };

  std::vector<EdgeInfo> & GetEdgeInfo() { return edge_info_; }
  const std::vector<EdgeInfo> & GetEdgeInfo() const { return edge_info_;}

 private:
  const Graph & g_;

  //Hash<MCTSElement, hash_type> * tt_;

  TopKData top_k_;

  void MakeRandomTable();
  std::vector<hash_type> random_table_;

  hash_type root_hash_code_;
  MCTSState root_state_;

  //bool expanded_;
  long long int expand_num_;

  std::vector<EdgeInfo> edge_info_; // list of item id sorted by pmin
  //std::vector<int> edge_index_; // needed? item id to edge id conversion

  // should be equal to best_pval_ ?
  //double upper_bound_; // initially 1.0 <- should be replaced with TopKData::worst_pval_ ?
  // pval < lower_bound_ will be ignored. will be used for top K
  //double lower_bound_; // initially zero
};

} // namespace lamp_search

#endif // _LAMP_SEARCH_MCTS_H_
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
