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

#include "variable_bitset_array.h"
#include "lcm_graph.h"
#include "lcm_dfs.h"
#include "lcm_dfs_vba.h"

DECLARE_int32(topk); // "find top K"
DECLARE_int32(threshold_delta); // "sup threshold delta for iterative deepening"
DECLARE_int32(threshold_init); // "initial sup threshold. if 0, use max sup(single item)"
DECLARE_bool(iter); // "do iterative deepening"

DECLARE_bool(prof); // false, "measure profile (PPCExt branch, Iter time histogram)"

DECLARE_double(a); // significance level alpha

DECLARE_bool(show_progress); // "show progress at each iteration"

DEFINE_bool(second_phase, true, "do second phase");
DEFINE_bool(third_phase, true, "do third phase");

namespace lamp_search {

const int LCM_DFS_VBA::k_int_max = std::numeric_limits<int>::max();
// assuming (digits in long long int) > (bits of double mantissa)
const long long int LCM_DFS_VBA::k_cs_max =
    1ll << (std::numeric_limits<double>::digits - 1);

void LCM_DFS_VBA_State::Init() {
  depth_ = 0;
  itemset_.Clear();
}

// LCM_DFS_VBA_State * LCM_DFS_VBA_State::New(std::size_t numbits) {
//   std::size_t n = traits::calc_nu_blocks(numbits);
//   assert(n>0);
//   std::size_t size = sizeof(LCM_DFS_VBA_State) + sizeof(uint64) * (n-1);
  
//   char * tmp_buf = new char[size];
//   return (LCM_DFS_VBA_State *)tmp_buf;
// }

//==============================================================================

LCM_DFS_VBA::LCM_DFS_VBA(const LCM_Graph_VBA<uint64> & g) :
    g_ (g),
    tbl_ (g_.GetTable()),
    bsh (g.GetTable().VBSHelper()),
    sup_stack_ (bsh),
    top_k_ (g.GetTable().GetItemInfo(), FLAGS_topk),
    lambda_max_ (0),
    lambda_ (1),
    sig_level_ (1.0),
    pmin_thr_ (NULL),
    cs_thr_ (NULL),
    cs_num_accum_ (NULL),
    final_closed_set_num_ (0ll),
    final_support_ (0) {
  assert(FLAGS_topk > 0);
  sup_bitset_buf_ = bsh.New();
}

LCM_DFS_VBA::~LCM_DFS_VBA() {
  bsh.Delete(sup_bitset_buf_);

  if (pmin_thr_ != NULL) delete [] pmin_thr_;
  if (cs_thr_ != NULL) delete [] cs_thr_;
  if (cs_num_accum_ != NULL) delete [] cs_num_accum_;
}

void LCM_DFS_VBA::DFSInit() {
  sup_stack_.Init();

  expand_num_ = 0ll;
  total_expand_num_ = 0ll;

  top_k_.Init();
}

void LCM_DFS_VBA::LAMPInit() {
  expand_num_ = 0ll;
  total_expand_num_ = 0ll;

  closed_set_num_ = 0ll;

  significant_list_.clear();

  // todo: implement
  // init hash table?
}

void LCM_DFS_VBA::LAMP() {
  timer_ = Timer::GetInstance();
  search_start_time_ = timer_->Elapsed();

  // define init threshold
  //int sup_threshold;
  int lambda; // == sup_threshold
  if (FLAGS_threshold_init == 0) lambda = tbl_.MaxT(); // should be MaxX()?
  else lambda = FLAGS_threshold_init;
  int num_iter = 0;
  total_expand_num_ = 0ll;

  LAMPInit();
  LCM_DFS_VBA_State root_state;
  root_state.Init();

  {// init progress
    if (FLAGS_show_progress) {
      std::cout << "first phase\n";
      std::cout << "sup_thr=" << lambda
                << "\tclosed_set_num=" << std::setw(12) << closed_set_num_
                << "\tcs_thr=" << std::setw(16) << "---";
      if (closed_set_num_)
        std::cout << "\tsig_lev=" << std::setw(12) << (FLAGS_a / closed_set_num_);
      else
        std::cout << "\tsig_lev=" << std::setw(12) << "---";
      std::cout << "\tpmin_thr=" << tbl_.PMin( lambda );
      std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                << "\ttotal_expand=" << std::setw(12) << total_expand_num_
                << "\telapsed_time=" << (timer_->Elapsed() - search_start_time_) / GIGA
                << "\n";
    }
  }

  long long int closed_set_num_1st = 0;
  int sup_threshold_1st = k_int_max;

  // phase 1
  while (true) {
    expand_num_ = 0ll;
    closed_set_num_ = 0ll;

    //double pmin_threshold = tbl_.PMin( sup_threshold );
    // long long int closed_set_num_threshold =
    //     (long long int)(std::min(std::floor(FLAGS_a/pmin_threshold),
    //                              (double)(k_cs_max)));
    double pmin_threshold_dec = tbl_.PMin( lambda-1 );
    long long int closed_set_num_threshold =
        (long long int)(std::min(std::floor(FLAGS_a/pmin_threshold_dec),
                                 (double)(k_cs_max)));
    // std::cout << "a=" << FLAGS_a
    //           << "\tspthr=" << sup_threshold
    //           << "\tpmin=" << pmin_threshold
    //           << "\tclthr=" << closed_set_num_threshold
    //           << std::endl;
    // note: closed_set_num_threshold can be zero if lambda==1
    //assert(closed_set_num_threshold);

    sup_stack_.Init(); // make stack for root position
    LAMPIter(&root_state, lambda, 0.0, false);
    num_iter++;
    total_expand_num_ += expand_num_;
    if (FLAGS_show_progress) {
      std::cout << "sup_thr=" << lambda
                << "\tclosed_set_num=" << std::setw(12) << closed_set_num_
                << "\tcs_thr=" << std::setw(16) << closed_set_num_threshold;
      if (closed_set_num_)
        std::cout << "\tsig_lev=" << std::setw(12) << (FLAGS_a / closed_set_num_);
      else
        std::cout << "\tsig_lev=" << std::setw(12) << "---";
      std::cout << "\tpmin_thr=" << pmin_threshold_dec;
      std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                << "\ttotal_expand=" << std::setw(12) << total_expand_num_
                << "\telapsed_time=" << (timer_->Elapsed() - search_start_time_) / GIGA
                << "\n";
    }

    if ( closed_set_num_ > closed_set_num_threshold) {
      sup_threshold_1st = lambda;
      closed_set_num_1st = closed_set_num_;
      break;
    }

    // todo: adjust delta and do binary search. or use hash table
    int delta = 1;
    lambda -= delta;

    if ( lambda <= 0 ) {
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
                << "\tcs_thr=" << std::setw(16) << "---";
      if (closed_set_num_1st)
        std::cout << "\tsig_lev=" << std::setw(12) << (FLAGS_a / closed_set_num_1st);
      else
        std::cout << "\tsig_lev=" << std::setw(12) << "---";
      std::cout << "\tpmin_thr=" << tbl_.PMin( sup_threshold_1st );
      std::cout << "\tnum_expand=" << std::setw(12) << "---"
                << "\ttotal_expand=" << std::setw(12) << total_expand_num_
                << "\telapsed_time=" << (timer_->Elapsed() - search_start_time_) / GIGA
                << "\n";
    }

    closed_set_num_ = 0ll;
    expand_num_ = 0ll;
    if (closed_set_num_1st) {
      sup_stack_.Init(); // make stack for root position
      LAMPIter(&root_state, sup_threshold_1st, (FLAGS_a / closed_set_num_1st), true);
    }
    total_expand_num_ += expand_num_;
    if (FLAGS_show_progress) {
      std::cout << "sup_thr=" << sup_threshold_1st
                << "\tclosed_set_num=" << std::setw(12) << closed_set_num_1st
                << "\tcs_thr=" << std::setw(16) << "---";
      if (closed_set_num_1st)
        std::cout << "\tsig_lev=" << std::setw(12) << (FLAGS_a / closed_set_num_1st);
      else
        std::cout << "\tsig_lev=" << std::setw(12) << "---";
      std::cout << "\tpmin_thr=" << tbl_.PMin( sup_threshold_1st );
      std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                << "\ttotal_expand=" << std::setw(12) << total_expand_num_
                << "\telapsed_time=" << (timer_->Elapsed() - search_start_time_) / GIGA
                << "\n";
    }
    closed_set_num_2nd = closed_set_num_;
  } else {
    closed_set_num_2nd = 0ll;
  }

  final_closed_set_num_ = closed_set_num_2nd;
  final_support_ = sup_threshold_1st;

  assert( closed_set_num_1st == closed_set_num_2nd );

  std::sort(significant_list_.begin(), significant_list_.end(), less);
}

void LCM_DFS_VBA::LAMPIter(LCM_DFS_VBA_State * st, int sup_threshold,
                           double sig_level, bool second_phase) {
  expand_num_++;

  // lookup, allocate
  // note: for hash table version
  // don't forget to init element after tt_.allocate()

  uint64 * sup = sup_stack_.Top();

  // pop stack for sequential version
  // reconstuct support from state for parallel version

  // depth limit check
  if (sup_stack_.Full()) {
    std::cout << "depth limit reached" << std::endl;
    return;
  }

  sup_stack_.Inc();
  uint64 * child_sup = sup_stack_.Top();

  SortedItemSet ppc_ext_buf;
  int core_i = g_.CoreIndex(st->itemset_);
  const SortedItemSet & current_itemset = st->itemset_;

  for ( int new_item = core_i + 1 ; new_item < tbl_.NuItems() ; new_item++ ) {
    // skip existing item
    if (current_itemset.Exists(new_item)) continue;

    // todo:
    // should add pseudo database reduction

    bsh.Copy(sup, child_sup);
    int sup_num = bsh.AndCountUpdate(tbl_.NthData(new_item), child_sup);

    if (sup_num == 0) continue;
    if (sup_num < sup_threshold) continue;

    bool res = g_.PPCExtension(current_itemset,
                               child_sup, core_i, new_item,
                               &ppc_ext_buf);

    if (res) {
      closed_set_num_++;

      if (second_phase) {
        int pos_sup_num = bsh.AndCount(tbl_.PosNeg()->Ptr(), child_sup);
        double pval = tbl_.PVal(sup_num, pos_sup_num);
        assert( pval >= 0.0 );
        if ( pval <= sig_level ) // permits == case?
          significant_list_.push_back(std::pair<double, SortedItemSet>(pval, ppc_ext_buf));
      }

      assert(sup_num >= sup_threshold);
      
      // try skipping if supnum_ == sup_threshold,
      // because if sup_num of a node equals to sup_threshold,
      // children will have smaller sup_num
      // therefore no need to check it's children
      if ( sup_num > sup_threshold && !ppc_ext_buf.Full() ) {
        LCM_DFS_VBA_State child_state; // it should not do any malloc operation
        child_state.itemset_ = ppc_ext_buf;
        child_state.depth_ = st->depth_+1;

        LAMPIter(&child_state, sup_threshold, sig_level, second_phase);
      }
    }
  }

  sup_stack_.Dec();
}

//==============================================================================

void LCM_DFS_VBA::LAMPSingleInit() {
  significant_map_.clear();

  lambda_max_ = tbl_.MaxX(); // the upper bound of sup

  if (pmin_thr_) delete [] pmin_thr_;
  if (cs_thr_) delete [] cs_thr_;
  if (cs_num_accum_) delete [] cs_num_accum_;

  pmin_thr_ = new double[lambda_max_+1];
  cs_thr_ = new long long int[lambda_max_+1];
  cs_num_accum_ = new long long int[lambda_max_+1];

  for (int i=0;i<=lambda_max_;i++)
    pmin_thr_[i] = tbl_.PMin(i);

  cs_thr_[0] = 0ll; // should not be used ???
  for (int i=1;i<=lambda_max_;i++) {
    cs_thr_[i] = (long long int)(std::min(std::floor( FLAGS_a / pmin_thr_[i-1]),
                                          (double)(k_cs_max)));
    // dbg
    // std::cout << "cs_thr[" << i << "]=\t" << cs_thr_[i] << "\tpmin_thr=" << pmin_thr_[i-1] << std::endl;
  }

  for (int i=0;i<=lambda_max_;i++)
    cs_num_accum_[i] = 0ll;
}

void LCM_DFS_VBA::LAMPSingle() {
  timer_ = Timer::GetInstance();
  search_start_time_ = timer_->Elapsed();

  expand_num_ = 0ll;
  closed_set_num_ = 0ll;

  LAMPSingleInit();

  total_expand_num_ = 0ll;
  lambda_ = 1;
  LCM_DFS_VBA_State root_state;
  root_state.Init();

  {// init progress
    if (FLAGS_show_progress) {
      std::cout << "# " << "1st phase start\n";
      std::cout << "# " << "lambda=" << lambda_
                << "\tclosed_set_num[n>=lambda]=" << std::setw(12) << cs_num_accum_[lambda_]
                << "\tcs_thr[lambda]=" << std::setw(16) << cs_thr_[lambda_];
      std::cout << "\tpmin_thr[lambda-1]=" << std::setw(12) << tbl_.PMin( lambda_-1 );
      std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                << "\telapsed_time=" << (timer_->Elapsed() - search_start_time_) / GIGA
                << std::endl;
    }
  }

  {
    sup_stack_.Init(); // make stack for root position
    LAMPSingleIter(&root_state);
    if (FLAGS_show_progress) {
      std::cout << "# " << "1st phase end\n";
      std::cout << "# " << "lambda=" << lambda_;
      std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                << "\telapsed_time=" << (timer_->Elapsed() - search_start_time_) / GIGA
                << std::endl;
    }
  }

  lambda_--;

  // dbg
  // std::cout << "lambda=" << lambda_ << std::endl;
  final_support_ = lambda_;

  if (!FLAGS_second_phase) return;

  if (FLAGS_show_progress) {
    std::cout << "# " << "2nd phase start\n";
    std::cout << "# " << "lambda=" << lambda_
              << "\telapsed_time=" << (timer_->Elapsed() - search_start_time_) / GIGA
              << std::endl;
  }

  LAMPInit();
  sup_stack_.Init(); // make stack for root position
  LAMPIter(&root_state, lambda_, 0.0, false);
  // dbg
  // std::cout << "closed_set_num=" << closed_set_num_ << std::endl;
  final_closed_set_num_ = closed_set_num_;

  if (FLAGS_show_progress) {
    std::cout << "# " << "2nd phase end\n";
    std::cout << "# "
              << "\tclosed_set_num=" << std::setw(12) << closed_set_num_
              << "\tsig_lev=" << (FLAGS_a / closed_set_num_)
              << "\telapsed_time=" << (timer_->Elapsed() - search_start_time_) / GIGA
              << std::endl;
  }

  if (!FLAGS_third_phase) return;

  // todo: replace with version recording closed item sets

  sup_stack_.Init(); // make stack for root position
  LAMPIter(&root_state, lambda_, (FLAGS_a / closed_set_num_), true);

  if (FLAGS_show_progress) {
    std::cout << "# " << "3rd phase end\n";
    std::cout << "# " << "sig_lev=" << (FLAGS_a / final_closed_set_num_)
              << "\telapsed_time=" << (timer_->Elapsed() - search_start_time_) / GIGA
              << std::endl;
  }

  std::sort(significant_list_.begin(), significant_list_.end(), less);
}

void LCM_DFS_VBA::LAMPSingleIter(LCM_DFS_VBA_State * st) {
  expand_num_++;

  uint64 * sup = sup_stack_.Top();

  // pop stack for sequential version
  // can reconstuct support from state for parallel version
  // for saving message size (maybe not needed)

  // depth limit check
  if (sup_stack_.Full()) {
    std::cout << "# " << "depth limit reached" << std::endl;
    return;
  }

  sup_stack_.Inc();
  uint64 * child_sup = sup_stack_.Top();

  // todo: replace with variable length itemset
  SortedItemSet ppc_ext_buf;
  int core_i = g_.CoreIndex(st->itemset_);
  const SortedItemSet & current_itemset = st->itemset_;

  // note: no explicit termination condition is needed, if I understand correctly
  for ( int new_item = core_i + 1 ; new_item < tbl_.NuItems() ; new_item++ ) {
    // skip existing item
    // todo: improve speed here
    if (current_itemset.Exists(new_item)) continue;

    bsh.Copy(sup, child_sup);
    int sup_num = bsh.AndCountUpdate(tbl_.NthData(new_item), child_sup);

    if (sup_num == 0) continue;
    if (sup_num < lambda_) continue;

    bool res = g_.PPCExtension(current_itemset,
                               child_sup, core_i, new_item,
                               &ppc_ext_buf);

    if (res) {
      // increment closed_set_num_array
      IncCsAccum(sup_num);

      // todo: how frequent this check should be done in parallel version?
      // check closed set num threshold
      if ( ExceedCsThr() ) {
        if (FLAGS_show_progress) {
          std::cout << "# " << "lambda=" << lambda_
                    << "\tclosed_set_num[n>=lambda]=" << std::setw(12) << cs_num_accum_[lambda_]
                    << "\tcs_thr[lambda]=" << std::setw(16) << cs_thr_[lambda_];
          // std::cout << "\tsig_lev=" << std::setw(12) << sig_level_;
          std::cout << "\tpmin_thr[lambda-1]=" << std::setw(12) << tbl_.PMin( lambda_-1 );
          std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                    << "\telapsed_time=" << (timer_->Elapsed() - search_start_time_) / GIGA
                    << std::endl;
        }

        IncLambdaThr(); // increment lambda_ (can be more than 1 in parallel)
        // todo: check correctness of the following
        // dbg
        // for (int si=1;si<=lambda_;si++)
        //   std::cout << "i=" << si
        //             << "\tcs=" << cs_num_accum_[si]
        //             << "\tcs_thr=" << cs_thr_[si] << std::endl;

        // if (cs_thr_[lambda_] != 0)
        //   sig_level_ = FLAGS_a / cs_thr_[lambda_];
        // else
        //   sig_level_ = FLAGS_a;
        // discard itemset with pval > updated sig_level from significant_map_
        // if (!FLAGS_count_only) DiscardSignificantList();

        if (sup_num < lambda_) continue;
      }

      // if (!FLAGS_count_only) {
      //   int pos_sup_num = bsh.AndCount(tbl_.PosNeg()->Ptr(), child_sup);
      //   double pval = tbl_.PVal(sup_num, pos_sup_num);
      //   if (pval < 0.0) continue; // ???
      //   if ( pval < sig_level_ ) {
      //     significant_map_.insert(std::pair<double, SortedItemSet>(pval, ppc_ext_buf));
      //   }
      // }

      assert(sup_num >= lambda_);

      // try skipping if supnum_ == sup_threshold,
      // because if sup_num of a node equals to sup_threshold, children will have smaller sup_num
      // therefore no need to check it's children
      if ( sup_num > lambda_ && !ppc_ext_buf.Full() ) {
        LCM_DFS_VBA_State child_state; // it should not do any malloc operation
        child_state.itemset_ = ppc_ext_buf;
        child_state.depth_ = st->depth_+1;

        LAMPSingleIter(&child_state);
      }
    }
  }

  sup_stack_.Dec();
}

void LCM_DFS_VBA::IncCsAccum(int sup_num) {
  for (int i=sup_num;i>=lambda_-1;i--)
    cs_num_accum_[i]++;
}

bool LCM_DFS_VBA::ExceedCsThr() const {
  // note: > is correct. permit ==
  return (cs_num_accum_[lambda_] > cs_thr_[lambda_]);
}

void LCM_DFS_VBA::IncLambdaThr() {
  // lambda_++;

  int si;
  for (si=lambda_max_;si>=lambda_;si--)
    if (cs_num_accum_[si] > cs_thr_[si]) break;
  //assert(lambda_ == si);
  lambda_ = si + 1;
  // it is safe because lambda_ higher than max results in immediate search finish
}

void LCM_DFS_VBA::DiscardSignificantList() {
  std::multimap< double, SortedItemSet>::reverse_iterator rit;
  for (rit=significant_map_.rbegin(); rit != significant_map_.rend();) {
    // permits == case
    if ((*rit).first > sig_level_) {
      significant_map_.erase(--rit.base());
    } else break;
  }
}

//==============================================================================

std::ostream & LCM_DFS_VBA::PrintSignificantList(std::ostream & out) const {
  std::stringstream s;

  for(std::vector< std::pair<double, SortedItemSet> >::const_iterator it = significant_list_.begin();
      it != significant_list_.end(); ++it) {
    s << "pval=" << std::setw(16) << std::left << (*it).first << std::right
      << "\titems";

    const SortedItemSet & ss = (*it).second;
    for(SortedItemSet::const_iterator it=ss.begin(); it !=ss.end(); ++it) {
      s << "\t" << tbl_.NthItemName(*it);
    }
    s << std::endl;
  }

  out << s.str() << std::flush;
  return out;
}

std::ostream & LCM_DFS_VBA::PrintSignificantMap(std::ostream & out) const {
  std::stringstream s;

  for(std::multimap< double, SortedItemSet >::const_iterator it = significant_map_.begin();
      it != significant_map_.end(); ++it) {
    s << "pval=" << std::setw(16) << std::left << (*it).first << std::right
      << "\titems";

    const SortedItemSet & ss = (*it).second;
    for(SortedItemSet::const_iterator it=ss.begin(); it !=ss.end(); ++it) {
      s << "\t" << tbl_.NthItemName(*it);
    }
    s << std::endl;
  }
  
  out << s.str() << std::flush;
  return out;
}

std::ostream & LCM_DFS_VBA::PrintLAMPResults(std::ostream & out) const {
  std::stringstream s;

  s << "min sup=" << final_support_
    << "\tfreq=" << final_closed_set_num_
    << std::endl;

  PrintSignificantList(s);

  out << s.str() << std::flush;
  return out;
}

std::ostream & LCM_DFS_VBA::PrintLAMPSingleResults(std::ostream & out) const {
  std::stringstream s;

  s << "min sup=" << final_support_;
  if (FLAGS_second_phase)
    s << "\tfreq=" << final_closed_set_num_;
  s << std::endl;

  if (FLAGS_third_phase) PrintSignificantList(s);

  out << s.str() << std::flush;
  return out;
}

} // namespace lamp_search

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
