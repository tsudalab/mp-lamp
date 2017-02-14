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
#include "table_vba.h"
#include "timer.h"

#include "sorted_itemset.h"
#include "topk.h"
#include "variable_bitset_array.h"
#include "variable_length_itemset.h"

#include "random.h"
#include "lcm_graph_vba.h"
#include "lamp_graph.h"
#include "lamp.h"

DECLARE_bool(prof); // false, "measure profile (PPCExt branch, Iter time histogram)"

DECLARE_double(a); // significance level alpha

DECLARE_bool(show_progress); // "show progress at each iteration"

DECLARE_bool(second_phase); // true, "do second phase"
DECLARE_bool(third_phase);  // true, "do third phase"

DEFINE_int32(stack_size, 1024*1024*64, "node stack size (integer)");
DEFINE_int32(freq_max, 1024*1024*64, "stack size for holding freq sets");

DEFINE_int32(min_sig_size, 0, "show at least top n combinations");
DEFINE_int32(max_sig_size, 0, "maximum size of significant item set. If 0, show all");
// fixme: add sanity check (or stop using this and implement re-run of 2nd phase)
DEFINE_int32(max_freq_size, 0,
             "maximum size of frequent item set. 0 for all.\n"
             "Used for 2nd phase. Should be significantly greater than max_sig_size");

namespace lamp_search {

const int Lamp::k_int_max = std::numeric_limits<int>::max();
// assuming (digits in long long int) > (bits of double mantissa)
const long long int Lamp::k_cs_max =
    1ll << (std::numeric_limits<double>::digits - 1);

//==============================================================================

Lamp::Lamp(const LampGraph<uint64> & g) :
    timer_ (Timer::GetInstance()),
    g_ (g),
    d_ (g.GetDatabase()),
    bsh_ (g.GetDatabase().VBSHelper()),
    sup_stack_ (bsh_),
    node_stack_ (NULL),
    lambda_max_ (0),
    lambda_ (1),
    sig_level_ (1.0),
    lambda_thr_ (1),
    pmin_thr_ (NULL),
    cs_thr_ (NULL),
    cs_accum_array_ (NULL),
    significant_stack_ (NULL), 
    sigset_comp_ (NULL),
    significant_set_ (NULL),
    sigset_record_mode_ (0),
    insignificant_itemset_num_ (0),
    expand_num_ (0ll),
    total_expand_num_ (0ll),
    closed_set_num_ (0ll),
    final_closed_set_num_ (0ll),
    final_support_ (0),
    sup_buf_ (NULL),
    child_sup_buf_ (NULL)
{
  log_.init_time_ = timer_->Elapsed();
  InitSearch();

  if (FLAGS_min_sig_size < 0)
    throw std::runtime_error(std::string("negative value not allowed for min_sig_size"));

  if (FLAGS_max_sig_size < 0)
    throw std::runtime_error(std::string("negative value not allowed for max_sig_size"));

  if (FLAGS_max_sig_size > 0 && FLAGS_min_sig_size > FLAGS_max_sig_size)
    throw std::runtime_error(std::string("option conflict min_sig_size and max_sig_size"));

  if (FLAGS_min_sig_size == 0) {
    if (FLAGS_max_sig_size == 0) { // default
      sigset_record_mode_ = SigsetRecordMode::NORMAL;
    } else { // show at most n
      sigset_record_mode_ = SigsetRecordMode::AT_MOST_N;
    }
  } else {
    if (FLAGS_max_sig_size == 0) { // show at least n
      sigset_record_mode_ = SigsetRecordMode::AT_LEAST_M;
    } else { // show at least n but at most n
      sigset_record_mode_ = SigsetRecordMode::M_TO_N;
    }
  }

  log_.init_time_ = timer_->Elapsed() - log_.init_time_;
}

Lamp::~Lamp() {
  if (node_stack_) delete node_stack_;
  if (significant_stack_) delete significant_stack_;
  if (significant_set_) delete significant_set_;
  if (sigset_comp_) delete sigset_comp_;
  if (pmin_thr_ != NULL) delete [] pmin_thr_;
  if (cs_thr_ != NULL) delete [] cs_thr_;
  if (cs_accum_array_ != NULL) delete [] cs_accum_array_;
  if (sup_buf_) bsh_.Delete(sup_buf_);
  if (child_sup_buf_) bsh_.Delete(child_sup_buf_);
}

void Lamp::InitSearch() {
  lambda_max_ = d_.MaxX(); // the upper bound of sup

  if (node_stack_) delete node_stack_;
  node_stack_ = new VariableLengthItemsetStack(FLAGS_stack_size);

  if (significant_stack_) delete significant_stack_;
  significant_stack_ = new VariableLengthItemsetStack(FLAGS_freq_max);

  if (sigset_comp_) delete sigset_comp_;
  sigset_comp_ = new sigset_compare(*significant_stack_);

  if (significant_set_)
    significant_set_->clear();
  else
    significant_set_ = new std::set<SignificantSetResult, sigset_compare>(*sigset_comp_);

  if (pmin_thr_) delete [] pmin_thr_;
  if (cs_thr_) delete [] cs_thr_;
  if (cs_accum_array_) delete [] cs_accum_array_;

  pmin_thr_ = new double[lambda_max_+1];
  cs_thr_ = new long long int[lambda_max_+1];
  cs_accum_array_ = new long long int[lambda_max_+1];

  for (int i=0;i<=lambda_max_;i++)
    pmin_thr_[i] = d_.PMin(i);

  cs_thr_[0] = 0ll; // should not be used ???
  for (int i=1;i<=lambda_max_;i++) {
    cs_thr_[i] = (long long int)(std::min(std::floor( FLAGS_a / pmin_thr_[i-1]),
                                          (double)(k_cs_max)));
    // dbg
    // std::cout << "cs_thr[" << i << "]=\t" << cs_thr_[i] << "\tpmin_thr=" << pmin_thr_[i-1] << std::endl;
  }

  for (int i=0;i<=lambda_max_;i++)
    cs_accum_array_[i] = 0ll;

  sup_buf_ = bsh_.New();
  child_sup_buf_ = bsh_.New();
}

void Lamp::Search() {
  log_.Init();
  log_.search_start_time_ = timer_->Elapsed();

  total_expand_num_ = 0ll;
  expand_num_ = 0ll;
  closed_set_num_ = 0ll;

  lambda_ = 1;

  // push root state to stack
  int * root_itemset;
  node_stack_->PushPre();
  root_itemset = node_stack_->Top();
  node_stack_->SetSup(root_itemset, lambda_max_);
  node_stack_->PushPostNoSort();

  {// init progress
    if (FLAGS_show_progress) {
      std::cout << "# " << "1st phase start\n";
      std::cout << "# " << "lambda=" << lambda_
                << "\tclosed_set_num[n>=lambda]=" << std::setw(12) << cs_accum_array_[lambda_]
                << "\tcs_thr[lambda]=" << std::setw(16) << cs_thr_[lambda_];
      std::cout << "\tpmin_thr[lambda-1]=" << std::setw(12) << d_.PMin( lambda_-1 );
      std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                << "\telapsed_time=" << (timer_->Elapsed() - log_.search_start_time_) / GIGA
                << std::endl;
    }
  }

  {
    sup_stack_.Init();

    long long int start_time = timer_->Elapsed();
    FirstPhaseIter(root_itemset);
    log_.process_node_time_ += timer_->Elapsed() - start_time;

    if (FLAGS_show_progress) {
      std::cout << "# " << "1st phase end\n";
      std::cout << "# " << "lambda=" << lambda_;
      std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                << "\telapsed_time=" << (timer_->Elapsed() - log_.search_start_time_) / GIGA
                << std::endl;
    }
  }

  lambda_--;
  final_support_ = lambda_;

  if (!FLAGS_second_phase) return;
  // todo: show progress

  expand_num_ = 0ll;
  closed_set_num_ = 0ll;

  // push root state to stack
  node_stack_->Clear();
  node_stack_->PushPre();
  root_itemset = node_stack_->Top();
  node_stack_->SetSup(root_itemset, lambda_max_);
  node_stack_->PushPostNoSort();

  // store significant closed sets in node_stack_
  // or, prepare std::multimap< double, int * > sorted_closed_set_
  // and sort based on pvalue

  double int_sig_lev = GetInterimSigLevel(lambda_);
  if (FLAGS_show_progress) {
    std::cout << "# " << "2nd phase start\n";
    std::cout << "# " << "lambda=" << lambda_
              << "\tint_tsig_lev=" << int_sig_lev
              << "\telapsed_time=" << (timer_->Elapsed() - log_.search_start_time_) / GIGA
              << std::endl;
  }

  {
    long long int start_time = timer_->Elapsed();
    LCMIter(root_itemset, lambda_, int_sig_lev);
    log_.process_node_time_ += timer_->Elapsed() - start_time;
  }

  final_closed_set_num_ = closed_set_num_;

  if (FLAGS_show_progress) {
    std::cout << "# " << "2nd phase end\n";
    std::cout << "# "
              << "closed_set_num=" << std::setw(12) << closed_set_num_
              << "\tsig_lev=" << (FLAGS_a / closed_set_num_)
              << "\tnum_expand=" << std::setw(12) << expand_num_
              << "\telapsed_time=" << (timer_->Elapsed() - log_.search_start_time_) / GIGA
              << std::endl;
  }

  // todo: show progress

  if (!FLAGS_third_phase) return;

  final_sig_level_ = FLAGS_a / final_closed_set_num_;
  DiscardSignificantList(final_sig_level_);

  if (FLAGS_show_progress) {
    std::cout << "# " << "3rd phase end\n";
    std::cout << "# " << "sig_lev=" << final_sig_level_
              << "\telapsed_time=" << (timer_->Elapsed() - log_.search_start_time_) / GIGA
              << std::endl;
  }
}

void Lamp::FirstPhaseIter(int * itemset) {
  expand_num_++;

  uint64 * sup = sup_stack_.Top();

  // pop stack for sequential version
  // can reconstuct support from state for parallel version
  // for saving message size (maybe not needed)

  // node_stack_ limit check is omitted
  // allocate enough memory!

  sup_stack_.Inc();
  uint64 * child_sup = sup_stack_.Top();

  assert(node_stack_->Top() == itemset);
  // int n = node_stack_->GetItemNum(itemset);
  int core_i = g_.CoreIndex(*node_stack_, itemset);

  int * ppc_ext_buf;
  for ( int new_item = core_i + 1 ; new_item < d_.NuItems() ; new_item++ ) {
    // skip existing item
    // todo: improve speed here
    if (node_stack_->Exist(itemset, new_item)) continue;

    bsh_.Copy(sup, child_sup);
    int sup_num = bsh_.AndCountUpdate(d_.NthData(new_item), child_sup);

    if (sup_num < lambda_) continue;

    node_stack_->PushPre();
    ppc_ext_buf = node_stack_->Top();

    bool res = g_.PPCExtension(node_stack_, itemset,
                               child_sup, core_i, new_item,
                               ppc_ext_buf);

    node_stack_->SetSup(ppc_ext_buf, sup_num);
    // todo: remove redundancy when res==false
    node_stack_->PushPostNoSort();

    if (res) {
      node_stack_->SortTop();
      // increment closed_set_num_array

      IncCsAccum(sup_num);

      // todo: how frequent this check should be done in parallel version?
      // check closed set num threshold
      if ( ExceedCsThr() ) {
        if (FLAGS_show_progress) {
          std::cout << "# " << "lambda=" << lambda_
                    << "\tclosed_set_num[n>=lambda]=" << std::setw(12)
                    << cs_accum_array_[lambda_]
                    << "\tcs_thr[lambda]=" << std::setw(16) << cs_thr_[lambda_];
          // std::cout << "\tsig_lev=" << std::setw(12) << sig_level_;
          std::cout << "\tpmin_thr[lambda-1]=" << std::setw(12) << d_.PMin( lambda_-1 );
          std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                    << "\telapsed_time="
                    << (timer_->Elapsed() - log_.search_start_time_) / GIGA
                    << std::endl;
        }

        IncLambdaThr(); // increment lambda_ (can be more than 1 in parallel)
        if (sup_num < lambda_) continue;
      }

      assert(sup_num >= lambda_);

      // try skipping if supnum_ == sup_threshold,
      // because if sup_num of a node equals to sup_threshold, children will have smaller sup_num
      // therefore no need to check it's children
      // note: skipping node_stack_ full check. allocate enough memory!
      if ( sup_num > lambda_ ) {
        int * child_itemset = node_stack_->Top();
        FirstPhaseIter(child_itemset);
      }
    }

    node_stack_->Pop();
  }
  sup_stack_.Dec();

  log_.process_node_num_ ++;
}

void Lamp::IncCsAccum(int sup_num) {
  for (int i=sup_num;i>=lambda_-1;i--)
    cs_accum_array_[i]++;
}

bool Lamp::ExceedCsThr() const {
  // note: > is correct. permit ==
  return (cs_accum_array_[lambda_] > cs_thr_[lambda_]);
}

void Lamp::IncLambdaThr() {
  // lambda_++;

  int si;
  for (si=lambda_max_;si>=lambda_;si--)
    if (cs_accum_array_[si] > cs_thr_[si]) break;
  //assert(lambda_ == si);
  lambda_ = si + 1;
  // it is safe because lambda_ higher than max results in immediate search finish
}

double Lamp::GetInterimSigLevel(int lambda) const {
  long long int csnum = cs_accum_array_[lambda];
  double lv;
  if (csnum > 0) lv = FLAGS_a / (double)csnum;
  else lv = FLAGS_a;

  return lv;
}

void Lamp::LCMIter(int * itemset, int sup_threshold, double sig_level) {
  expand_num_++;

  uint64 * sup = sup_stack_.Top();

  // pop stack for sequential version
  // can reconstuct support from state for parallel version
  // for saving message size (maybe not needed)

  // node_stack_ limit check is omitted
  // allocate enough memory!

  sup_stack_.Inc();
  uint64 * child_sup = sup_stack_.Top();

  assert(node_stack_->Top() == itemset);
  // int n = node_stack_->GetItemNum(itemset);
  int core_i = g_.CoreIndex(*node_stack_, itemset);

  int * ppc_ext_buf;

  for ( int new_item = core_i + 1 ; new_item < d_.NuItems() ; new_item++ ) {
    // skip existing item
    // todo: improve speed of Exist (but only meaningful for large itemset)
    if (node_stack_->Exist(itemset, new_item)) continue;

    bsh_.Copy(sup, child_sup);
    int sup_num = bsh_.AndCountUpdate(d_.NthData(new_item), child_sup);

    if (sup_num < lambda_) continue;

    node_stack_->PushPre();
    ppc_ext_buf = node_stack_->Top();

    bool res = g_.PPCExtension(node_stack_, itemset,
                               child_sup, core_i, new_item,
                               ppc_ext_buf);

    node_stack_->SetSup(ppc_ext_buf, sup_num);
    // todo: remove redundancy when res==false
    node_stack_->PushPostNoSort();

    if (res) {
      node_stack_->SortTop();
      closed_set_num_++;

      // dbg
      // std::cout << "# found cs ";
      // node_stack_->Print(std::cout, ppc_ext_buf);

      if (FLAGS_third_phase) {
        int pos_sup_num = bsh_.AndCount(d_.PosNeg(), child_sup);
        double pval = d_.PVal(sup_num, pos_sup_num);
        assert( pval >= 0.0 );

        RecordSignificantItemset(pval, sig_level, sup_num, pos_sup_num, ppc_ext_buf);
      }

      assert(sup_num >= sup_threshold);

      // try skipping if supnum_ == sup_threshold,
      // because if sup_num of a node equals to sup_threshold,
      // children will have smaller sup_num
      // therefore no need to check it's children
      if ( sup_num > sup_threshold ) {
        int * child_itemset = node_stack_->Top();
        LCMIter(child_itemset, sup_threshold, sig_level);
      }
    }

    node_stack_->Pop();
  }

  sup_stack_.Dec();

  log_.process_node_num_ ++;
}

void Lamp::DiscardSignificantList(double sig_level) {
  std::set<SignificantSetResult, sigset_compare>::reverse_iterator rit;
  switch(sigset_record_mode_) {
    case SigsetRecordMode::NORMAL:
      {
        for(rit = significant_set_->rbegin(); rit != significant_set_->rend();) {
          // permits == case
          if (rit->pval_ > sig_level) {
            significant_set_->erase(--rit.base());
          }
          else break;
        }
      }
      break;
    case SigsetRecordMode::AT_MOST_N:
      {
        assert( significant_set_->size() <= (std::size_t)FLAGS_max_sig_size); // does this hold?

        for(rit = significant_set_->rbegin(); rit != significant_set_->rend();) {
          // permits == case
          if (rit->pval_ > sig_level || significant_set_->size() > (std::size_t)FLAGS_max_sig_size) {
            significant_set_->erase(--rit.base());
          }
          else break;
        }
      }
      break;
    case SigsetRecordMode::AT_LEAST_M:
      {
        for(rit = significant_set_->rbegin(); rit != significant_set_->rend();) {
          // permits == case
          if (rit->pval_ > sig_level && significant_set_->size() > (std::size_t)FLAGS_min_sig_size) {
            significant_set_->erase(--rit.base());
          }
          else break;
        }
      }
      break;
    case SigsetRecordMode::M_TO_N:
      {
        for(rit = significant_set_->rbegin(); rit != significant_set_->rend();) {
          if ( (rit->pval_ > sig_level && significant_set_->size() > (std::size_t)FLAGS_min_sig_size) ||
               significant_set_->size() > (std::size_t)FLAGS_max_sig_size ) {
            significant_set_->erase(--rit.base());
          }
          else break;
        }
      }
      break;
    default:
      throw std::runtime_error(std::string("unknown mode in RecordFrequentItemset"));
      break;
  }

}

// fixme:
// bug in AT_MOST_N and M_TO_N
// because 2nd phase sig level is relaxed
void Lamp::RecordSignificantItemset(double pval, double sig_level,
                                    int sup_num, int pos_sup_num,
                                    int * itemset) {
  // fixme: change some command line options from int32 to int64

  switch(sigset_record_mode_) {
    case SigsetRecordMode::NORMAL:
      if (pval <= sig_level) { // == means significant
        int * item = PushItemsetNoSort(significant_stack_, itemset);
        significant_set_->insert(
            SignificantSetResult(pval, item, sup_num, pos_sup_num)
                                 );
      }
      break;

    case SigsetRecordMode::AT_MOST_N:
      // todo: add pruning based on top_n

      if (significant_set_->size() < (std::size_t)FLAGS_max_sig_size) {
        if (pval <= sig_level) { // == means significant
          int * item = PushItemsetNoSort(significant_stack_, itemset);
          significant_set_->insert(
              SignificantSetResult(pval, item, sup_num, pos_sup_num)
                                   );
        }
      } else {
        std::set<SignificantSetResult, sigset_compare>::iterator it = significant_set_->end();
        --it;

        // erase, prepare for insert and return true
        if (!sigset_comp_->compare(*it, pval, itemset)) {
          significant_set_->erase(it);
          int * item = PushItemsetNoSort(significant_stack_, itemset);
          significant_set_->insert(
              SignificantSetResult(pval, item, sup_num, pos_sup_num)
                                   );
        }
      }
      break;

    case SigsetRecordMode::AT_LEAST_M:
      if (significant_set_->size() < (std::size_t)FLAGS_min_sig_size) {
        int * item = PushItemsetNoSort(significant_stack_, itemset);
        significant_set_->insert(
            SignificantSetResult(pval, item, sup_num, pos_sup_num)
                                 );

        if (!(pval <= sig_level)) insignificant_itemset_num_++;

      } else {
        std::set<SignificantSetResult, sigset_compare>::iterator it = significant_set_->end();
        --it;

        // insert if new one is significant
        if (pval <= sig_level) { // == means significant
          // discard if worst one is not significant
          if (insignificant_itemset_num_ > 0) {
            std::set<SignificantSetResult, sigset_compare>::iterator it = significant_set_->end();
            --it;
            significant_set_->erase(it);
            assert(it->pval_ > sig_level);
            insignificant_itemset_num_--;
          }

          int * item = PushItemsetNoSort(significant_stack_, itemset);
          significant_set_->insert(
              SignificantSetResult(pval, item, sup_num, pos_sup_num)
                                   );
        }

      }
      break;

    case SigsetRecordMode::M_TO_N:
      if (significant_set_->size() < (std::size_t)FLAGS_min_sig_size) {
        int * item = PushItemsetNoSort(significant_stack_, itemset);
        significant_set_->insert(
            SignificantSetResult(pval, item, sup_num, pos_sup_num)
                                 );

        if (!(pval <= sig_level)) insignificant_itemset_num_++;

      } else if (significant_set_->size() < (std::size_t)FLAGS_max_sig_size) {
        std::set<SignificantSetResult, sigset_compare>::iterator it = significant_set_->end();
        --it;

        // insert if new one is significant
        if (pval <= sig_level) { // == means significant
          // discard if worst one is not significant
          if (insignificant_itemset_num_ > 0) {
            std::set<SignificantSetResult, sigset_compare>::iterator it = significant_set_->end();
            --it;
            significant_set_->erase(it);
            assert(it->pval_ > sig_level);
            insignificant_itemset_num_--;
          }

          int * item = PushItemsetNoSort(significant_stack_, itemset);
          significant_set_->insert(
              SignificantSetResult(pval, item, sup_num, pos_sup_num)
                                   );
        }

      } else { // count >= FLAGS_max_sig_size
        std::set<SignificantSetResult, sigset_compare>::iterator it = significant_set_->end();
        --it;

        // erase, prepare for insert and return true
        if (!sigset_comp_->compare(*it, pval, itemset)) {
          significant_set_->erase(it);
          int * item = PushItemsetNoSort(significant_stack_, itemset);
          significant_set_->insert(
              SignificantSetResult(pval, item, sup_num, pos_sup_num)
                                   );
        }

      }
      // todo: can combine with pruning based on top n?

      break;

    default:
      throw std::runtime_error(std::string("unknown mode in RecordSignificantItemset"));
      break;
  }

  // // dbg
  // std::cout << "insert: "
  //           << "pval=" << pval
  //           << "\titems: ";
  // significant_stack_->Print(std::cout, item);
}

//==============================================================================

void Lamp::SearchLoop() {
  log_.search_start_time_ = timer_->Elapsed();

  total_expand_num_ = 0ll;
  expand_num_ = 0ll;
  closed_set_num_ = 0ll;

  lambda_ = 1;

  // push root state to stack
  int * root_itemset;
  node_stack_->PushPre();
  root_itemset = node_stack_->Top();
  node_stack_->SetSup(root_itemset, lambda_max_);
  node_stack_->PushPostNoSort();

  {// init progress
    if (FLAGS_show_progress) {
      std::cout << "# " << "1st phase start\n";
      std::cout << "# " << "lambda=" << lambda_
                << "\tclosed_set_num[n>=lambda]=" << std::setw(12) << cs_accum_array_[lambda_]
                << "\tcs_thr[lambda]=" << std::setw(16) << cs_thr_[lambda_];
      std::cout << "\tpmin_thr[lambda-1]=" << std::setw(12) << d_.PMin( lambda_-1 );
      std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                << "\telapsed_time=" << (timer_->Elapsed() - log_.search_start_time_) / GIGA
                << std::endl;
    }
  }

  {
    sup_stack_.Init();

    long long int start_time = timer_->Elapsed();
    FirstPhaseLoop();
    log_.process_node_time_ += timer_->Elapsed() - start_time;

    if (FLAGS_show_progress) {
      std::cout << "# " << "1st phase end\n";
      std::cout << "# " << "lambda=" << lambda_;
      std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                << "\telapsed_time=" << (timer_->Elapsed() - log_.search_start_time_) / GIGA
                << std::endl;
    }
  }

  lambda_--;
  final_support_ = lambda_;

  if (!FLAGS_second_phase) return;
  // todo: show progress

  expand_num_ = 0ll;
  closed_set_num_ = 0ll;

  // push root state to stack
  node_stack_->Clear();
  node_stack_->PushPre();
  root_itemset = node_stack_->Top();
  node_stack_->SetSup(root_itemset, lambda_max_);
  node_stack_->PushPostNoSort();

  // store significant closed sets in node_stack_
  // or, prepare std::multimap< double, int * > sorted_closed_set_
  // and sort based on pvalue

  double int_sig_lev = GetInterimSigLevel(lambda_);
  if (FLAGS_show_progress) {
    std::cout << "# " << "2nd phase start\n";
    std::cout << "# " << "lambda=" << lambda_
              << "\tint_sig_lev=" << int_sig_lev
              << "\telapsed_time=" << (timer_->Elapsed() - log_.search_start_time_) / GIGA
              << std::endl;
  }

  {
    sig_level_ = int_sig_lev;
    lambda_thr_ = lambda_;

    long long int start_time = timer_->Elapsed();
    LCMLoop();
    log_.process_node_time_ += timer_->Elapsed() - start_time;
  }

  final_closed_set_num_ = closed_set_num_;

  if (FLAGS_show_progress) {
    std::cout << "# " << "2nd phase end\n";
    std::cout << "# "
              << "\tclosed_set_num=" << std::setw(12) << closed_set_num_
              << "\tsig_lev=" << (FLAGS_a / closed_set_num_)
              << "\telapsed_time=" << (timer_->Elapsed() - log_.search_start_time_) / GIGA
              << std::endl;
  }

  if (!FLAGS_third_phase) return;

  final_sig_level_ = FLAGS_a / final_closed_set_num_;
  DiscardSignificantList(final_sig_level_);

  if (FLAGS_show_progress) {
    std::cout << "# " << "3rd phase end\n";
    std::cout << "# " << "sig_lev=" << final_sig_level_
              << "\telapsed_time=" << (timer_->Elapsed() - log_.search_start_time_) / GIGA
              << std::endl;
  }
}

void Lamp::FirstPhaseLoop() {
  while(!node_stack_->Empty()) {
    expand_num_++;

    node_stack_->CopyItem(node_stack_->Top(), itemset_buf_);
    node_stack_->Pop();

    // calculate support from itemset_buf_
    bsh_.Set(sup_buf_);
    {
      int n = node_stack_->GetItemNum(itemset_buf_);
      for (int i=0;i<n;i++) {
        int item = node_stack_->GetNthItem(itemset_buf_, i);
        bsh_.And(d_.NthData(item), sup_buf_);
      }
    }

    int core_i = g_.CoreIndex(*node_stack_, itemset_buf_);

    int * ppc_ext_buf;
    // todo: use database reduction
    // reverse order
    for ( int new_item = d_.NuItems()-1 ; new_item >= core_i+1 ; new_item-- ) {
      // skip existing item
      // todo: improve speed here
      if (node_stack_->Exist(itemset_buf_, new_item)) continue;

      bsh_.Copy(sup_buf_, child_sup_buf_);
      int sup_num = bsh_.AndCountUpdate(d_.NthData(new_item), child_sup_buf_);
      
      if (sup_num < lambda_) continue;

      node_stack_->PushPre();
      ppc_ext_buf = node_stack_->Top();

      bool res = g_.PPCExtension(node_stack_, itemset_buf_,
                                 child_sup_buf_, core_i, new_item,
                                 ppc_ext_buf);

      node_stack_->SetSup(ppc_ext_buf, sup_num);
      node_stack_->PushPostNoSort();

      if (!res) {// todo: remove this redundancy
        node_stack_->Pop();
      }
      else {
        node_stack_->SortTop();

        // increment closed_set_num_array
        IncCsAccum(sup_num);

        assert(sup_num >= lambda_);

        // try skipping if supnum_ == sup_threshold,
        // because if sup_num of a node equals to sup_threshold, children will have smaller sup_num
        // therefore no need to check it's children
        // note: skipping node_stack_ full check. allocate enough memory!
        if ( sup_num <= lambda_ ) {
          node_stack_->Pop();
        }
      }
    }

    if ( ExceedCsThr() ) {
      if (FLAGS_show_progress) {
        std::cout << "# " << "lambda=" << lambda_
                  << "\tclosed_set_num[n>=lambda]=" << std::setw(12)
                  << cs_accum_array_[lambda_]
                  << "\tcs_thr[lambda]=" << std::setw(16) << cs_thr_[lambda_];
        // std::cout << "\tsig_lev=" << std::setw(12) << sig_level_;
        std::cout << "\tpmin_thr[lambda-1]=" << std::setw(12) << d_.PMin( lambda_-1 );
        std::cout << "\tnum_expand=" << std::setw(12) << expand_num_
                  << "\telapsed_time="
                  << (timer_->Elapsed() - log_.search_start_time_) / GIGA
                  << std::endl;
      }
      IncLambdaThr(); // increment lambda_ (can be more than 1 in parallel)
    }
    
    log_.process_node_num_ ++;
  }
}

void Lamp::LCMLoop() {
  while(!node_stack_->Empty()) {
    expand_num_++;

    node_stack_->CopyItem(node_stack_->Top(), itemset_buf_);
    node_stack_->Pop();

    bsh_.Set(sup_buf_);
    {
      int n = node_stack_->GetItemNum(itemset_buf_);
      for (int i=0;i<n;i++) {
        int item = node_stack_->GetNthItem(itemset_buf_, i);
        bsh_.And(d_.NthData(item), sup_buf_);
      }
    }

    int core_i = g_.CoreIndex(*node_stack_, itemset_buf_);

    int * ppc_ext_buf;
    // reverse order
    for ( int new_item = d_.NuItems()-1 ; new_item >= core_i+1 ; new_item-- ) {
      // skip existing item
      // todo: improve speed here
      if (node_stack_->Exist(itemset_buf_, new_item)) continue;

      bsh_.Copy(sup_buf_, child_sup_buf_);
      int sup_num = bsh_.AndCountUpdate(d_.NthData(new_item), child_sup_buf_);
      
      if (sup_num < lambda_thr_) continue;
      // zoe 2017/01/25
      // to add pruning, add something here
      // like, if (depth_of_child >= threshold) continue;

      node_stack_->PushPre();
      ppc_ext_buf = node_stack_->Top();

      bool res = g_.PPCExtension(node_stack_, itemset_buf_,
                                 child_sup_buf_, core_i, new_item,
                                 ppc_ext_buf);

      node_stack_->SetSup(ppc_ext_buf, sup_num);
      node_stack_->PushPostNoSort();

      if (!res) {// todo: remove this redundancy
        node_stack_->Pop();
      }
      else {
        node_stack_->SortTop();
        closed_set_num_++;

        if (FLAGS_third_phase) {
          int pos_sup_num = bsh_.AndCount(d_.PosNeg(), child_sup_buf_);
          double pval = d_.PVal(sup_num, pos_sup_num);
          assert( pval >= 0.0 );

          RecordSignificantItemset(pval, sig_level_, sup_num, pos_sup_num, ppc_ext_buf);
        }

        assert(sup_num >= lambda_thr_);

        // try skipping if supnum_ == sup_threshold,
        // because if sup_num of a node equals to sup_threshold,
        // children will have smaller sup_num
        // therefore no need to check it's children
        // note: skipping node_stack_ full check. allocate enough memory!
        if ( sup_num <= lambda_thr_ ) {
          node_stack_->Pop();
        }
      }
    }

    log_.process_node_num_ ++;
  }
}

//==============================================================================

std::ostream & Lamp::PrintResults(std::ostream & out) const {
  std::stringstream s;

  s << "# min. sup=" << final_support_;
  if (FLAGS_second_phase)
    s << "\tcorrection factor=" << final_closed_set_num_;
  s << std::endl;

  if (FLAGS_third_phase) PrintSignificantSet(s);

  out << s.str() << std::flush;
  return out;
}

std::ostream & Lamp::PrintSignificantSet(std::ostream & out) const {
  std::stringstream s;

  s << "# number of significant patterns=" << significant_set_->size() << std::endl;
  s << "# pval (raw)    pval (corr)         freq     pos        # items items\n";
  for(std::set<SignificantSetResult, sigset_compare>::const_iterator it
          = significant_set_->begin();
      it != significant_set_->end(); ++it) {

    s << ""   << std::setw(16) << std::left << (*it).pval_ << std::right
      << ""  << std::setw(16) << std::left << (*it).pval_ * final_closed_set_num_ << std::right
      << "" << std::setw(8)  << (*it).sup_num_
      << ""  << std::setw(8)  << (*it).pos_sup_num_
      << "";
    // s << "pval (raw)=" << std::setw(16) << std::left << (*it).pval_ << std::right
    //   << "pval (corr)="  << std::setw(16) << std::left << (*it).pval_ * final_closed_set_num_ << std::right
    //   << "\tfreq=" << std::setw(8)  << (*it).sup_num_
    //   << "\tpos="  << std::setw(8)  << (*it).pos_sup_num_
    //   << "\titems";

    const int * item = (*it).set_;
    significant_stack_->Print(s, d_.ItemNames(), item);
  }
  
  out << s.str() << std::flush;
  return out;
}

Lamp::Log::Log() {
  Init();
}

void Lamp::Log::Init() {
  init_time_ = 0ll;
  search_start_time_ = 0ll;
  process_node_num_ = 0ll;
  process_node_time_ = 0ll;
}

std::ostream & Lamp::PrintLog(std::ostream & out) const {
  std::stringstream s;

  s << "# process_node_num  ="
    << std::setw(16) << log_.process_node_num_
    << std::endl;
  s << "# process_node_time ="
    << std::setw(16) << log_.process_node_time_ / MEGA
    << "(ms)" << std::endl;
  s << "# node / second     ="
    << std::setw(16)
    << (log_.process_node_num_) / (log_.process_node_time_ / GIGA)
    << std::endl;

  out << s.str() << std::flush;
  return out;
}

} // namespace lamp_search

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
