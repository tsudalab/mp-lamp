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
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>

#include <boost/random.hpp>
#include <boost/foreach.hpp>

#include "gflags/gflags.h"

#include "utils.h"
#include "sorted_itemset.h"
#include "topk.h"
#include "timer.h"
#include "lcm_graph.h"
#include "lcm_dfs.h"

DECLARE_int32(topk); // "find top K"
DECLARE_int32(threshold_delta); // "sup threshold delta for iterative deepening"
DECLARE_int32(threshold_init); // "initial sup threshold. if 0, use max sup(single item)"
DECLARE_bool(iter); // "do iterative deepening"

DEFINE_bool(prof, false, "measure profile (PPCExt branch, Iter time histogram)");

DECLARE_double(a); // significance level alpha

DECLARE_bool(show_progress); // "show progress at each iteration"

namespace lamp_search {

LCM_DFSState::LCM_DFSState()
{}

void LCM_DFSState::Init(const LCM_DFS & search) {
  const Table & tbl = search.GetGraph().GetTable();

  edge_stack_.clear();
  SortedItemSet empty_set;
  edge_stack_.push_back(empty_set);
  depth_ = 0;
  // todo: prepare dynamic bitset array
  bitset_diff_stack_.clear();

  support_.resize(tbl.NuTransaction(), true); // all bit set for root state
  support_buf_.resize(tbl.NuTransaction(), true); // initialize buffer also

  sup_bitset_buf_.resize(tbl.NuTransaction(), true);
  //item_bitset_buf_.resize(tbl.NuItems(), false);
}

void LCM_DFSState::Push(SortedItemSet new_items,
                        const boost::dynamic_bitset<> & new_support) {
  LCMEdge edge(new_items);
  edge_stack_.push_back( edge );
  depth_++;

  // todo: need to avoid malloc here
  boost::dynamic_bitset<> bs_diff = support_ ^ new_support;
  bitset_diff_stack_.push_back(bs_diff);
  support_ = new_support;
}

void LCM_DFSState::Pop() {
  // todo: rewrite to avoid malloc
  support_ ^= bitset_diff_stack_.back();
  bitset_diff_stack_.pop_back();

  depth_--;
  assert(depth_ >= 0);

  // todo:
  // something needed for the following line
  // currently only used for top_k_.Insert()
  //item_bitset_buf_.reset( item_stack_.back() );

  edge_stack_.pop_back();
}

//==============================================================================

LCM_DFS::LCM_DFS(const LCM_Graph & g) :
    g_ (g),
    top_k_ (g.GetTable().GetItemInfo(), FLAGS_topk),
    final_closed_set_num_ (0ll),
    final_support_ (0)
{
  assert(FLAGS_topk > 0);
}

LCM_DFS::~LCM_DFS() {}

void LCM_DFS::InitDFS() {
  //const Table & tbl = GetGraph().GetTable();
  // const std::vector<Table::ItemInfo> & ivec = tbl.GetItemInfo();
  // for (std::vector<Table::ItemInfo>::const_iterator it = ivec.begin();
  //      it != ivec.end(); ++it) {
  //   edge_info_.push_back( ItemInfo(*it) );
  // }

  root_state_.Init(*this);

  expand_num_ = 0ll;
  total_expand_num_ = 0ll;

  top_k_.Init();
}

// simple dfs on LCM space
// no hash table
// can be safely pruned if (pmin > pval_threshold)
void LCM_DFS::Search() {
  InitDFS();
  // prepare root node (done in Init)
  const Table & tbl = GetGraph().GetTable();
  Timer * g_timer = Timer::GetInstance();
  long long int search_start_time = g_timer->Elapsed();

  // define init threshold
  int sup_threshold;
  if (FLAGS_threshold_init == 0) sup_threshold = tbl.MaxT(); // should be MaxX()?
  else sup_threshold = FLAGS_threshold_init;
  double pmin_threshold = tbl.PMin( sup_threshold );

  int num_iter = 0;
  total_expand_num_ = 0ll;

  {// init progress
    if (FLAGS_show_progress) {
      std::cout << "sup_thr=" << sup_threshold
                << "\tpmin_thr=" << pmin_threshold
                << "\tbest_pval=" << top_k_.BestPVal()
                << "\tkth_pval=" << top_k_.KthPVal()
                << "\tnum expand=" << std::setw(12) << expand_num_
                << "\ttotal expand=" << std::setw(12) << total_expand_num_
                << "\telapsed time=" << (g_timer->Elapsed() - search_start_time) / GIGA
                << "\n";
    }
  }

  while (true) {
    expand_num_ = 0ll;

    Iter(&root_state_, sup_threshold, pmin_threshold);
    num_iter++;
    if (FLAGS_show_progress) {
      std::cout << "sup_thr=" << sup_threshold
                << "\tpmin_thr=" << pmin_threshold
                << "\tbest_pval=" << top_k_.BestPVal()
                << "\tkth_pval=" << top_k_.KthPVal()
                << "\tnum expand=" << std::setw(12) << expand_num_
                << "\ttotal expand=" << std::setw(12) << total_expand_num_
                << "\telapsed time=" << (g_timer->Elapsed() - search_start_time) / GIGA
                << "\n";
    }

    total_expand_num_ += expand_num_;

    if (!FLAGS_iter) break;
    if (pmin_threshold >= top_k_.KthPVal()) break;
    if (expand_num_ == 0) break; // useless
    if (sup_threshold == 1) break;
    if (tbl.PMin(sup_threshold - 1) >= top_k_.KthPVal()) break;

    int delta;
    for (delta=1; delta<=FLAGS_threshold_delta; delta++) {
      double next_pmin_threshold = tbl.PMin(sup_threshold - delta);
      if ( next_pmin_threshold >= top_k_.KthPVal() || sup_threshold - delta <= 0)
        break;
    }
    delta--;
    sup_threshold -= delta;
    pmin_threshold = tbl.PMin( sup_threshold );
  }

  // note: for hash table version
  // don't forget to init element after tt_.allocate()
}

//==============================================================================

// dfs over LCM space
void LCM_DFS::Iter(LCM_DFSState * st, int sup_threshold, double minp_threshold) {
  const Table & tbl = GetGraph().GetTable();
  expand_num_++;

  boost::dynamic_bitset<> * sup = st->Support();

  SortedItemSet ppc_ext_buf;
  int core_i = g_.CoreIndex(st->ItemSet());
  for ( int new_item = core_i + 1 ; new_item < tbl.NuItems() ; new_item++ ) {
    // skip existing item
    const SortedItemSet & current_itemset = st->ItemSet();
    if (current_itemset.Exists(new_item)) continue;

    bool res=false;
    // todo:
    // should add pseudo database reduction
    // todo:
    // can skip part of closure/core_i calculation
    // by excluding transactions based on previous iteration results ?

    st->sup_bitset_buf_ = (*sup & tbl.NthData(new_item));
    int sup_num = st->sup_bitset_buf_.count();
    if (sup_num == 0) continue;
    double minp = tbl.PMin(sup_num);
    if (minp > top_k_.KthPVal()) continue;

    // todo: rewrite to avoid malloc, use fixed buffer
    boost::dynamic_bitset<> pos_sup = st->sup_bitset_buf_ & tbl.PosNeg();
    int pos_sup_num = pos_sup.count();
    double pval = tbl.PVal(sup_num, pos_sup_num);
    if (pval < 0.0) continue; // ???

    if (sup_num < sup_threshold) continue; // ???

    res = g_.PPCExtension(current_itemset,
                          st->sup_bitset_buf_, core_i, new_item,
                          &ppc_ext_buf);

    if (res) {
      if ( pval < top_k_.KthPVal() ) top_k_.Insert(pval, ppc_ext_buf);

      assert(sup_num >= sup_threshold);

      if ( !ppc_ext_buf.Full() ) {
        st->Push(ppc_ext_buf, st->sup_bitset_buf_);
        Iter(st, sup_threshold, minp_threshold);
        st->Pop();
      }
    }
  }
}

void LCM_DFS::InitLAMP() {
  root_state_.Init(*this);

  expand_num_ = 0ll;
  total_expand_num_ = 0ll;

  closed_set_num_ = 0ll;

  // init priority queue?
}

void LCM_DFS::LAMP() {
  InitLAMP();
  const int default_value = std::numeric_limits<int>::max();

  const Table & tbl = GetGraph().GetTable();
  Timer * g_timer = Timer::GetInstance();
  long long int search_start_time = g_timer->Elapsed();

  // define init threshold
  int sup_threshold;
  if (FLAGS_threshold_init == 0) sup_threshold = tbl.MaxT(); // should be MaxX()?
  else sup_threshold = FLAGS_threshold_init;
  int num_iter = 0;
  total_expand_num_ = 0ll;

  {// init progress
    if (FLAGS_show_progress) {
      std::cout << "first phase\n";
      std::cout << "sup_thr=" << sup_threshold
                << "\tclosed_set_num=" << std::setw(12) << closed_set_num_
                << "\tcs_thr=" << std::setw(12) << "---";
      if (closed_set_num_)
        std::cout << "\tsig_lev=" << std::setw(12) << (FLAGS_a / closed_set_num_);
      else
        std::cout << "\tsig_lev=" << std::setw(12) << "---";
      std::cout << "\tpmin_thr=" << tbl.PMin( sup_threshold );
      std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                << "\ttotal_expand=" << std::setw(12) << total_expand_num_
                << "\telapsed_time=" << (g_timer->Elapsed() - search_start_time) / GIGA
                << "\n";
    }
  }

  long long int closed_set_num_1st = 0;
  int sup_threshold_1st = default_value;

  // long long int closed_set_num_prev = 0;
  // int sup_threshold_prev = default_value;

  // phase 1
  while (true) {
    expand_num_ = 0ll;
    closed_set_num_ = 0ll;

    // double pmin_threshold = tbl.PMin( sup_threshold );
    // long long int closed_set_num_threshold =
    //     (long long int)(std::min(std::floor(FLAGS_a/pmin_threshold),
    //                              (double)(default_value)));
    double pmin_threshold_dec = tbl.PMin( sup_threshold-1 );
    long long int closed_set_num_threshold =
        (long long int)(std::min(std::floor(FLAGS_a/pmin_threshold_dec),
                                 (double)(default_value)));

    // std::cout << "a=" << FLAGS_a
    //           << "\tspthr=" << sup_threshold
    //           << "\tpmin=" << pmin_threshold
    //           << "\tclthr=" << closed_set_num_threshold
    //           << std::endl;
    assert(closed_set_num_threshold);

    LAMPIter(&root_state_, sup_threshold, 0.0, false);
    num_iter++;
    total_expand_num_ += expand_num_;
    if (FLAGS_show_progress) {
      std::cout << "sup_thr=" << sup_threshold
                << "\tclosed_set_num=" << std::setw(12) << closed_set_num_
                << "\tcs_thr=" << std::setw(12) << closed_set_num_threshold;
      if (closed_set_num_)
        std::cout << "\tsig_lev=" << std::setw(12) << (FLAGS_a / closed_set_num_);
      else
        std::cout << "\tsig_lev=" << std::setw(12) << "---";
      std::cout << "\tpmin_thr=" << pmin_threshold_dec;
      std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                << "\ttotal_expand=" << std::setw(12) << total_expand_num_
                << "\telapsed_time=" << (g_timer->Elapsed() - search_start_time) / GIGA
                << "\n";
    }

    if ( closed_set_num_ > closed_set_num_threshold ) {
      sup_threshold_1st = sup_threshold;
      closed_set_num_1st = closed_set_num_;
      break;
    }

    // todo: adjust delta and do binary search. or use hash table
    int delta = 1;
    sup_threshold -= delta;

    if ( sup_threshold <= 0 ) {
      sup_threshold_1st = 0;
      closed_set_num_1st = 0ll;
      break;
    }
  }

  // phase 2
  long long int closed_set_num_2nd = 0ll;
  if (sup_threshold_1st > 0) {
    if (FLAGS_show_progress) {
      std::cout << "second phase\n";
      std::cout << "sup_thr=" << sup_threshold_1st
                << "\tclosed_set_num=" << std::setw(12) << closed_set_num_1st
                << "\tcs_thr=" << std::setw(12) << "---";
      if (closed_set_num_1st)
        std::cout << "\tsig_lev=" << std::setw(12) << (FLAGS_a / closed_set_num_1st);
      else
        std::cout << "\tsig_lev=" << std::setw(12) << "---";
      std::cout << "\tpmin_thr=" << tbl.PMin( sup_threshold_1st );
      std::cout << "\tnum_expand=" << std::setw(12) << "---"
                << "\ttotal_expand=" << std::setw(12) << total_expand_num_
                << "\telapsed_time=" << (g_timer->Elapsed() - search_start_time) / GIGA
                << "\n";
    }

    closed_set_num_ = 0ll;
    expand_num_ = 0ll;
    if (closed_set_num_1st)
      LAMPIter(&root_state_, sup_threshold_1st,(FLAGS_a / closed_set_num_1st), true);
    total_expand_num_ += expand_num_;
    if (FLAGS_show_progress) {
      std::cout << "sup_thr=" << sup_threshold_1st
                << "\tclosed_set_num=" << std::setw(12) << closed_set_num_1st
                << "\tcs_thr=" << std::setw(12) << "---";
      if (closed_set_num_1st)
        std::cout << "\tsig_lev=" << std::setw(12) << (FLAGS_a / closed_set_num_1st);
      else
        std::cout << "\tsig_lev=" << std::setw(12) << "---";
      std::cout << "\tpmin_thr=" << tbl.PMin( sup_threshold_1st );
      std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                << "\ttotal_expand=" << std::setw(12) << total_expand_num_
                << "\telapsed_time=" << (g_timer->Elapsed() - search_start_time) / GIGA
                << "\n";
    }
    closed_set_num_2nd = closed_set_num_;
  }
  else {
    closed_set_num_2nd = 0ll;
  }

  final_closed_set_num_ = closed_set_num_2nd;
  final_support_ = sup_threshold_1st;

  assert( closed_set_num_1st == closed_set_num_2nd );

  std::sort(significant_list_.begin(), significant_list_.end(), less);

  // note: for hash table version
  // don't forget to init element after tt_.allocate()
}

void LCM_DFS::LAMPIter(LCM_DFSState * st, int sup_threshold,
                       double sig_level, bool second_phase) {
  const Table & tbl = GetGraph().GetTable();
  expand_num_++;

  boost::dynamic_bitset<> * sup = st->Support();

  SortedItemSet ppc_ext_buf;
  int core_i = g_.CoreIndex(st->ItemSet());
  const SortedItemSet & current_itemset = st->ItemSet();

  for ( int new_item = core_i + 1 ; new_item < tbl.NuItems() ; new_item++ ) {
    // skip existing item
    if (current_itemset.Exists(new_item)) continue;

    bool res=false;
    // todo:
    // should add pseudo database reduction
    // todo:
    // can skip part of closure/core_i calculation
    // by excluding transactions based on previous iteration results ?

    st->sup_bitset_buf_ = (*sup & tbl.NthData(new_item));

    int sup_num = st->sup_bitset_buf_.count();

    if (sup_num == 0) continue;
    if (sup_num < sup_threshold) continue;

    res = g_.PPCExtension(current_itemset,
                          st->sup_bitset_buf_, core_i, new_item,
                          &ppc_ext_buf);

    if (res) {
      closed_set_num_++;

      if (second_phase) {
        // todo: rewrite to avoid malloc, use fixed buffer
        boost::dynamic_bitset<> pos_sup = st->sup_bitset_buf_ & tbl.PosNeg();
        int pos_sup_num = pos_sup.count();
        double pval = tbl.PVal(sup_num, pos_sup_num);
        if (pval < 0.0) continue; // ???
        if ( pval < sig_level )
          significant_list_.push_back(std::pair<double, SortedItemSet>(pval, ppc_ext_buf));
      }

      assert(sup_num >= sup_threshold);

      // try skipping if supnum_ == sup_threshold,
      // because if sup_num of a node equals to sup_threshold, children will have smaller sup_num
      // therefore no need to check it's children
      if ( sup_num > sup_threshold && !ppc_ext_buf.Full() ) {
        st->Push(ppc_ext_buf, st->sup_bitset_buf_);
        LAMPIter(st, sup_threshold, sig_level, second_phase);
        st->Pop();
      }
    }
  }
}

std::ostream & LCM_DFS::PrintResults(std::ostream & out) const {
  const Table & tbl = GetGraph().GetTable();
  std::stringstream s;

  for(TopKData::const_iterator it = top_k_.pval_index_pair_.begin();
      it != top_k_.pval_index_pair_.end(); ++it) {
    s << "pval=" << std::setw(16) << std::left << (*it).first << std::right
      << "\titems";

    const SortedItemSet & ss = top_k_.itemset_[(*it).second];
    for(SortedItemSet::const_iterator it=ss.begin(); it !=ss.end(); ++it) {
      s << "\t" << tbl.NthItemName(*it);
    }
    s << std::endl;
  }
  
  out << s.str() << std::flush;
  return out;
}  

std::ostream & LCM_DFS::PrintLAMPResults(std::ostream & out) const {
  const Table & tbl = GetGraph().GetTable();
  std::stringstream s;

  s << "min sup=" << final_support_
    << "\t#freq=" << final_closed_set_num_
    << std::endl;

  for(std::vector< std::pair<double, SortedItemSet> >::const_iterator it = significant_list_.begin();
      it != significant_list_.end(); ++it) {
    s << "pval=" << std::setw(16) << std::left << (*it).first << std::right
      << "\titems";

    const SortedItemSet & ss = (*it).second;
    for(SortedItemSet::const_iterator it=ss.begin(); it !=ss.end(); ++it) {
      s << "\t" << tbl.NthItemName(*it);
    }
    s << std::endl;
  }
  
  out << s.str() << std::flush;
  return out;
}  

} // namespace lamp_search
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
