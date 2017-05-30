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

#ifndef _LAMP_SEARCH_LAMP_H_
#define _LAMP_SEARCH_LAMP_H_

#include <vector>
#include <algorithm>

#include <boost/array.hpp>
#include <map>

#include "gflags/gflags.h"

#include "utils.h"
#include "table_vba.h"
#include "timer.h"

#include "sorted_itemset.h"
#include "topk.h"
#include "variable_bitset_array.h"
#include "variable_length_itemset.h"

#include "random.h"
#include "lcm_dfs_vba.h"
#include "lamp_graph.h"

namespace lamp_search {

class Lamp {
 public:
  Lamp(const LampGraph<uint64> & g);
  ~Lamp();

  // call before search
  void InitSearch();
  void Search();

  // 1st phase
  void FirstPhaseIter(int * itemset);

  bool ExceedCsThr() const;
  void IncLambdaThr();
  void IncCsAccum(int sup_num);

  // 2nd phase
  void LCMIter(int * itemset, int sup_threshold, double sig_level);
  double GetInterimSigLevel(int lambda) const;

  // 3rd phase
  // collect significant sets
  void DiscardSignificantList(double sig_level);

  // no iteration version which resembles behavior of mp-lamp3
  void SearchLoop();
  void FirstPhaseLoop();
  void LCMLoop();

  // implement display method for each mode
  std::ostream & PrintResults(std::ostream & out) const;
  // std::ostream & PrintSignificantMap(std::ostream & out) const;
  std::ostream & PrintSignificantSet(std::ostream & out) const;
  std::ostream & PrintLog(std::ostream & out) const;

  const LampGraph<uint64> & Graph() const { return g_; }

  class Log {
   public:
    Log();
    void Init();

    long long int init_time_;
    long long int search_start_time_;

    long long int process_node_num_;
    long long int process_node_time_;
  };

  Timer * timer_;
  Log log_;

  static const int k_int_max;
  // assuming (digits in long long int) > (bits of double mantissa)
  static const long long int k_cs_max;

 private:
  const LampGraph<uint64> & g_;
  const Database<uint64> & d_;
  const VariableBitsetHelper<uint64> & bsh_; // bitset helper

  SupportStack sup_stack_;
  VariableLengthItemsetStack * node_stack_;

  // variables for single phase LAMP
  // equals to maximum support of single item
  int lambda_max_;
  // lambda_ initially set to 1. will be incremented to N if cs_thr[N] exceeded
  // in 1st phase, search will be pruned if (sup_num < lambda_)
  int lambda_;
  // todo: initial value can be set after checking all single item itemsets

  // initially set to 1. set to 0.05 (FLAGS_a) / cs_thr_[global_sup_thr_-1]
  // used for loop version
  double sig_level_;
  int lambda_thr_; // used for loop version

  double * pmin_thr_; // pmin_thr[sup] == tbl.PMin(sup), maybe redundant
  long long int * cs_thr_; // cs_thr[sup] shows closed set num threshold
  // cs_accum_array_[sup] shows closed set num with support higher than or equals to sup
  long long int * cs_accum_array_;      // int array of 0..lambda_max_ (size lambda_max_+1)
  //int * cs_accum_array_recv_; // int array of 0..lambda_max_ (size lambda_max_+1)

  class SignificantSetResult {
   public:
    SignificantSetResult(double p, int * s, int nu_sup, int nu_pos)
        : pval_(p),
          set_(s),
          sup_num_(nu_sup),
          pos_sup_num_(nu_pos)
    {}

    double pval_;
    int * set_;
    int sup_num_;
    int pos_sup_num_;
  };

  struct sigset_compare {
    sigset_compare(VariableLengthItemsetStack & ss) : ss_(ss) {}

    bool operator()(const SignificantSetResult & lhs,
                    const SignificantSetResult & rhs) {
      if (lhs.pval_ < rhs.pval_) return true;
      else if (lhs.pval_ > rhs.pval_) return false;
      else {
        int l_item_num = ss_.GetItemNum(lhs.set_);
        int r_item_num = ss_.GetItemNum(rhs.set_);

        if (l_item_num > r_item_num) return true;
        else if (l_item_num < r_item_num) return false;
        else {
          // sort based on dictionary order of item
          int n = l_item_num;
          for (int i=0;i<n;i++) {
            int l_item = ss_.GetNthItem(lhs.set_, i);
            int r_item = ss_.GetNthItem(rhs.set_, i);
            if (l_item < r_item) return true;
            else if (l_item > r_item) return false;
          }
          throw std::runtime_error("identical duplicate itemsets found");
          return false;
        }
      }
      return false;
    }

    bool compare(const SignificantSetResult & lhs,
                 double rhs_pval, int * rhs_item) {
      if (lhs.pval_ < rhs_pval) return true;
      else if (lhs.pval_ > rhs_pval) return false;
      else {
        int l_item_num = ss_.GetItemNum(lhs.set_);
        int r_item_num = ss_.GetItemNum(rhs_item);

        if (l_item_num > r_item_num) return true;
        else if (l_item_num < r_item_num) return false;
        else {
          // sort based on dictionary order of item
          int n = l_item_num;
          for (int i=0;i<n;i++) {
            int l_item = ss_.GetNthItem(lhs.set_, i);
            int r_item = ss_.GetNthItem(rhs_item, i);
            if (l_item < r_item) return true;
            else if (l_item > r_item) return false;
          }
          throw std::runtime_error("identical duplicate itemsets found");
          return false;
        }
      }
      return false;
    }

    const VariableLengthItemsetStack & ss_;
  };

  /** pushing item, assuming items in itemsets are sorted */
  int * PushItemsetNoSort(VariableLengthItemsetStack * ss, int * itemset) {
    ss->PushPre();
    int * item = ss->Top();
    ss->CopyItem(itemset, item);
    ss->PushPostNoSort();
    return item;
  }

  /** record itemset if needed, following the SigsetRecordMode */
  void RecordSignificantItemset(double pval, double sig_level,
                                int sup_num, int pos_sup_num, int * itemset);

  // std::multimap< double, int * > significant_map_;
  // static bool less(const std::pair<double, int * > & lhs,
  //                  const std::pair<double, int * > & rhs) {
  //   if (lhs.first < rhs.first) return true;
  //   return false;
  // }
  VariableLengthItemsetStack * significant_stack_;

  sigset_compare * sigset_comp_;

  /** holds significant set result. store itemset pointers and sort only the pointers */
  std::set<SignificantSetResult, sigset_compare> * significant_set_;

  struct SigsetRecordMode {
    enum SigsetRecordModeType {
      NORMAL = 0, // default, record all significant patterns
      AT_MOST_N, // record at most N significant patterns
      AT_LEAST_M, // show at least M patterns regardless of significance
      M_TO_N, // show at least M and at most N
    };
  };

  int sigset_record_mode_;
  int insignificant_itemset_num_;

  // expand num of current iteration (== number of Iter calls)
  // without positive check, this should be (the number of closed set) + 1
  long long int expand_num_;
  long long int total_expand_num_;

  long long int closed_set_num_;

  // for printing results
  long long int final_closed_set_num_;
  int final_support_;
  double final_sig_level_;

  int itemset_buf_[VariableLengthItemsetStack::kMaxItemsPerSet];
  uint64 * sup_buf_, * child_sup_buf_;
};

} // namespace lamp_search

#endif // _LAMP_SEARCH_LAMP_H_
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
