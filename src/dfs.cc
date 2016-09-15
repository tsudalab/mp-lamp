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

#include <boost/random.hpp>
#include <boost/foreach.hpp>

#include "gflags/gflags.h"

#include "utils.h"
#include "timer.h"
#include "variable_bitset_array.h"
#include "table.h"
#include "graph.h"
#include "sorted_itemset.h"
#include "topk.h"
#include "dfs.h"

DECLARE_int32(ttsize);

DEFINE_int32(topk, 3, "find top K");
DEFINE_int32(threshold_delta, 1, "sup threshold delta for iterative deepening");
DEFINE_int32(threshold_init, 0, "initial sup threshold. if 0, use max sup(single item)");
DEFINE_bool(iter, true, "do iterative deepening");

DEFINE_bool(show_progress, false, "show progress at each iteration");

namespace lamp_search {

DFSState::DFSState()
{}

void DFSState::Init(const DFS & search, hash_type root_hash_code) {
  const Table & tbl = search.GetGraph().GetTable();

  last_edge_index_in_item_stack_ = -1;
  
  hash_code_ = root_hash_code;
  node_.Init();

  // initialize hash_code_, node_ when needed
  item_stack_.clear();
  edge_stack_.clear();
  depth_ = 0;
  // todo: prepare dynamic bitset array
  bitset_diff_stack_.clear();

  support_.resize(tbl.NuTransaction(), true); // all bit set for root state
  support_buf_.resize(tbl.NuTransaction(), true); // initialize buffer also

  sup_bitset_buf_.resize(tbl.NuTransaction(), true);
  //item_bitset_buf_.resize(tbl.NuItems(), false);
}

DFS::DFS(const Graph & g) :
    g_ (g),
    tt_ (0),
    top_k_ (g.GetTable().GetItemInfo(), FLAGS_topk)
{
  assert(FLAGS_topk > 0);
  Init();
}

DFS::~DFS() {
  if (tt_) delete tt_;
}

void DFS::Init() {
  if (tt_) delete tt_;

  std::size_t elm_size = sizeof(Graph::Node);
  std::size_t elm_num = FLAGS_ttsize / elm_size;

  //tt_ = new Hash<uint8>(elm_num, false); // (elm_num, use_gc)

  MakeRandomTable();
  root_hash_code_.set(0ull, 0ull);
  root_state_.Init(*this, root_hash_code_);

  const Table & tbl = GetGraph().GetTable();
  const std::vector< ItemInfo > & ivec = tbl.GetItemInfo();
  for (std::vector< ItemInfo >::const_iterator it = ivec.begin();
       it != ivec.end(); ++it) {
    edge_info_.push_back( EdgeInfo(*it) );
  }

  expand_num_ = 0ll;
  total_expand_num_ = 0ll;

  top_k_.Init();
}

// brute force dfs version using only safe pruning
// no hash table
// can be safely pruned if (pmin > pval_threshold)
// can be unsafely pruned if ( estimated_best_pval > pval_threshold )
void DFS::Search() {
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

  while (true) {
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
    expand_num_ = 0ll;

    Iter(&root_state_, sup_threshold, pmin_threshold);
    num_iter++;

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
    //if (sup_threshold <= 0) sup_threshold = 1;
    pmin_threshold = tbl.PMin( sup_threshold );
  }

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

  // note:
  // don't forget to init element allocate after tt_.allocate()
}

void DFSState::Push(int edge_index, int item_id, const boost::dynamic_bitset<> & item_sup) {
  edge_stack_.push_back( edge_index );
  item_stack_.push_back( item_id );
  assert(last_edge_index_in_item_stack_ < edge_index); // move to unittest?
  last_edge_index_in_item_stack_ = edge_index;
  depth_++;

  // todo: need to avoid malloc here
  boost::dynamic_bitset<> bs_diff = support_ ^ item_sup;
  bitset_diff_stack_.push_back(bs_diff);
  // support_buf_ = support_; // bs_diff
  // support_buf_ ^= item_sup; // bs_diff
  // // it should be like
  // //boost::dynamic_bitset<> & temp = bitset_diff_stack_.hoge();
  // //temp = support_ ^ item_sup;
  // 
  // // stack only need max search depth size (+ 1 ?) (kMaxSearchDepth + 1?)
  // // it is needed to avoid calling constructors when pushing to the stack
  // // replace push_back with bitwise copy
  // // need to implement own stack
  // 
  // bitset_diff_stack_.push_back(support_buf_);
  support_ = item_sup;

  itemset_buf_.Push(item_id);
  //item_bitset_buf_.set( item_id );
}

void DFSState::Pop() {
  // todo: rewrite to avoid malloc
  support_ ^= bitset_diff_stack_.back();
  bitset_diff_stack_.pop_back();

  depth_--;
  assert(depth_ >= 0);

  itemset_buf_.Remove(item_stack_.back());
  //item_bitset_buf_.reset( item_stack_.back() );
  //item_bitset_buf_.reset( GetEdgeInfo()[edge_stack_.back()].id);

  edge_stack_.pop_back();
  item_stack_.pop_back();
  if (edge_stack_.empty())
    last_edge_index_in_item_stack_ = -1;
  else
    last_edge_index_in_item_stack_ = edge_stack_.back();
}

// simple brute force version
void DFS::Iter(DFSState * st, int sup_threshold, double minp_threshold) {
  const Table & tbl = GetGraph().GetTable();
  expand_num_++;

  boost::dynamic_bitset<> * sup = st->Support();
  // for (std::vector<int>::const_iterator it=st->ItemStack().begin();
  //      it != st->ItemStack().end() ; ++it)
  //   st->item_bitset_buf_.set(*it);

  const std::vector<EdgeInfo> & einfo = GetEdgeInfo();
  for (std::size_t ei = st->last_edge_index_in_item_stack_+1; ei < einfo.size(); ei++) {

    // improve this operation
    // continue as soon as minp / 0 condition satisfied
    st->sup_bitset_buf_ = (*sup) & tbl.NthData(einfo[ei].id);

    int sup_num = st->sup_bitset_buf_.count();
    if (sup_num == 0) continue;
    double minp = tbl.PMin(sup_num);
    if (minp > top_k_.KthPVal()) continue;

    // todo: rewrite to avoid malloc, use fixed buffer
    boost::dynamic_bitset<> pos_sup = st->sup_bitset_buf_ & tbl.PosNeg();
    int pos_sup_num = pos_sup.count();
    double pval = tbl.PVal(sup_num, pos_sup_num);
    if (pval < 0.0) continue; // ???

    if ( pval < top_k_.KthPVal() && !st->itemset_buf_.Full() ) {
      st->itemset_buf_.Push(einfo[ei].id);
      top_k_.Insert(pval, st->itemset_buf_);
      st->itemset_buf_.Remove(einfo[ei].id);
      // st->item_bitset_buf_.set(einfo[ei].id);
      // top_k_.Insert(pval, st->item_bitset_buf_);
      // st->item_bitset_buf_.reset(einfo[ei].id);
    }

    if (sup_num >= sup_threshold && !st->itemset_buf_.Full() ) {
      // edge_index, item_index, support_bitset
      st->Push( ei, GetEdgeInfo()[ei].id, st->sup_bitset_buf_ );
      Iter(st, sup_threshold, minp_threshold);
      st->Pop();
    }
  }

  // {
    // if support == parent's support
    // record as candidate for closed itemset (1)

    // record support in (map / hash table) (2)
    // check duplicate support for closed itemset candidate
  // }

  // todo:
  // for (1), create superset and makemoves (and prune rest)
    
  // if duplicate, merge them and do only one makemove

  // check closed pattern
  // if child.support == node.support, record only child

  // topk info should include t (and x?) and num of items (==depth)
  // if t is same, check if these are superset
  // items in these superset should have the highest search priority
  // (maybe explicit depth extension?)
}

void DFS::MakeRandomTable() {
  // todo: prepare method for seed initialization
  uint64 seed = 0x1234ll;

  boost::mt19937 rng(seed);
  boost::uniform_int<uint64> dst(0, std::numeric_limits<uint64>::max());
  boost::variate_generator<boost::mt19937&, boost::uniform_int<uint64> >
      rand_u64(rng, dst);

  std::size_t nu_items = g_.GetTable().NuItems();
  // std::cout << "nu_items=" << nu_items << std::endl;
  random_table_.reserve(nu_items);

  //std::bitset<128> bs;
  //uint128 ui;
  hash_type h;

  for (std::size_t bi = 0 ; bi < nu_items ; bi++ ) {
    // bs.set( rand_u64() );
    // bs <<= 64;
    // bs &= rand_u64();
    // random_table_.push_back(bs);
    uint64 r1, r2;
    r1 = rand_u64(); r2 = rand_u64();
    // ui = 0;
    // ui |= r1;
    // ui <<= 64;
    // ui |= r2;
    // random_table_.push_back(ui);
    h.set(r1, r2);
  }
}

std::ostream & DFS::DumpRandomTable(std::ostream & out) const {
  std::stringstream s;

  std::size_t nu_items = g_.GetTable().NuItems();
  for (std::size_t bi = 0 ; bi < nu_items ; bi++ ) {
    //uint128 ui = random_table_[bi];
    hash_type h = random_table_[bi];
    uint64 ui1, ui2;
    // ui1 = ui & 0xffffffffffffffff;
    // ui2 = (ui >> 64) & 0xffffffffffffffff;
    h.get(&ui2, &ui1);
    s << std::hex << std::setw(16) << ui2
      << std::hex << std::setw(16) << ui1 << std::dec << std::endl;
  }
  out << s.str() << std::flush;
  return out;
}

// std::ostream & DFS::DumpTopK(std::ostream & out) const {
//   std::stringstream s;

//   for(TopKData::const_iterator it = top_k_.pval_index_pair_.begin();
//       it != top_k_.pval_index_pair_.end(); ++it) {
//     s << "pval=" << std::setw(12) << (*it).first
//       << ":set=" << top_k_.itemset_[(*it).second]
//       << ":idx=" << (*it).second
//       << std::endl;
//   }
  
//   out << s.str() << std::flush;
//   return out;
// }  

std::ostream & DFS::PrintResults(std::ostream & out) const {
  const Table & tbl = GetGraph().GetTable();
  std::stringstream s;

  for(TopKData::const_iterator it = top_k_.pval_index_pair_.begin();
      it != top_k_.pval_index_pair_.end(); ++it) {
    s << "pval=" << std::setw(16) << std::left << (*it).first << std::right
      << "\titems";
    const SortedItemSet & ss = top_k_.itemset_[(*it).second];
    for (SortedItemSet::const_iterator it=ss.begin();it!=ss.end();++it) {
      s << "\t" << tbl.NthItemName(*it);
    }
    // const boost::dynamic_bitset<> & bs = top_k_.itemset_[(*it).second];
    // for (std::size_t bi=bs.find_first() ; bi < bs.npos ; bi = bs.find_next(bi) ) {
    //   s << "\t" << tbl.NthItemName(bi);
    // }
    s << std::endl;
  }
  
  out << s.str() << std::flush;
  return out;
}  

} // namespace lamp_search
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
