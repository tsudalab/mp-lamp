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

#ifndef _LAMP_SEARCH_DFS_H_
#define _LAMP_SEARCH_DFS_H_

#include <vector>
#include <algorithm>

#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>

#include "utils.h"
#include "variable_bitset_array.h"
#include "table.h"
#include "graph.h"
#include "hash.h"
#include "sorted_itemset.h"
#include "topk.h"

namespace lamp_search {

class DFS;

class DFSElement {
 public:
  DFSElement * getNext() { return next_; }
  void setNext(DFSElement * p) { next_ = p; }

  void Init();

  hash_type getHashcode() { return hash_code_; }

 private:
  DFSElement * next_;
  hash_type hash_code_;
  Graph::Node node_;
};

/** class for recording search state.
 ** basically, use one State for one thread/process
 ** for parallel version, state will be sent as messages */
class DFSState {
 public:
  DFSState();

  hash_type getHashcode() const { return hash_code_; }
  void Init(const DFS & search, hash_type hash_code);

  // push item to path and update selected_depth_
  void Push(int edge_index, int item_id, const boost::dynamic_bitset<> & item_sup);
  void Pop();

  const std::vector<int> & ItemStack() const { return item_stack_; }
  int last_edge_index_in_item_stack_;
  const std::vector< boost::dynamic_bitset<> > & BSDiffStack() { return bitset_diff_stack_; }

  boost::dynamic_bitset<> * Support() { return &support_; }

  int Depth() const { return depth_; }

  // todo: if only used as buffers, move to DFS
  //       or, preapre non-bitset current_item_set
  // todo: change variable name
  boost::dynamic_bitset<> sup_bitset_buf_;
  SortedItemSet itemset_buf_;
  //boost::dynamic_bitset<> item_bitset_buf_;

  // todo: add excluded_item_bitset / transanction_bitset

 private:
  hash_type hash_code_;
  Graph::Node node_;

  std::vector<int> item_stack_; // sequence of selected items
  std::vector<int> edge_stack_; // sequence of selected edge id
  int depth_;

  // hoge;
  // SupportDiffStack support_diff_stack_;
  std::vector< boost::dynamic_bitset<> > bitset_diff_stack_;
  boost::dynamic_bitset<> support_; // current support
  boost::dynamic_bitset<> support_buf_; // buffer for support calculation

  // on iter, record t and support for all items?
};

/** depth first search class */
class DFS {
 public:
  DFS(const Graph & g);
  ~DFS();
  
  void Search();

  void Iter(DFSState * st, int sup_threshold, double minp_threshold);

  void Init();

  // must call init after this ?
  void SetTopK(int k) { top_k_.ResetK(k); }

  std::ostream & DumpRandomTable(std::ostream & out) const;

  //std::ostream & DumpTopK(std::ostream & out) const;

  std::ostream & PrintResults(std::ostream & out) const;

  const Graph & GetGraph() const { return g_; }

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

  TopKData & TopK() { return top_k_; }
  const TopKData & TopK() const { return top_k_; }

  std::vector<EdgeInfo> & GetEdgeInfo() { return edge_info_; }
  const std::vector<EdgeInfo> & GetEdgeInfo() const { return edge_info_;}

 private:
  const Graph & g_;

  Hash<uint8> * tt_;

  TopKData top_k_;

  void MakeRandomTable();
  std::vector<hash_type> random_table_;

  hash_type root_hash_code_;
  DFSState root_state_;

  long long int expand_num_; // of current iteration
  long long int total_expand_num_;

  std::vector<EdgeInfo> edge_info_; // list of item id sorted by pmin
  //std::vector<int> edge_index_; // needed? item id to edge id conversion

  // should be equal to best_pval_ ?
  //double upper_bound_; // initially 1.0 <- should be replaced with TopKData::worst_pval_ ?
  // pval < lower_bound_ will be ignored. will be used for top K
  //double lower_bound_; // initially zero
};

} // namespace lamp_search

#endif // _LAMP_SEARCH_DFS_H_
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
