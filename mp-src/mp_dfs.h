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

#ifndef _LAMP_SEARCH_MP_DFS_H_
#define _LAMP_SEARCH_MP_DFS_H_

#include <vector>
#include <algorithm>

#include <boost/array.hpp>
#include <map>

#include "mpi.h"

#include "utils.h"
#include "sorted_itemset.h"
#include "topk.h"
#include "variable_bitset_array.h"
#include "variable_length_itemset.h"
#include "table_vba.h"
#include "timer.h"

#include "random.h"
#include "database.h"
#include "lamp_graph.h"

namespace lamp_search {

// use these as mpi tag
struct Tag {
  enum TaskType {
    // control tasks
    CONTROL_TASK_BEGIN=0,
    DTD_REQUEST,
    DTD_REPLY,

    DTD_ACCUM_REQUEST, // request reporting accum count
    DTD_ACCUM_REPLY, // reduce closed set count

    BCAST_FINISH,
    CONTROL_TASK_END,
  
    // basic tasks
    BASIC_TASK_BEGIN,
    
    LAMBDA, // for lamp

    REQUEST, // request tasks
    REJECT, // reject requests
    GIVE, // send tasks

    BASIC_TASK_END,

    // third phase tasks
    THIRD_PHASE_BEGIN,

    RESULT_REQUEST,
    RESULT_REPLY,

    THIRD_PHASE_END,
  };
};

class FixedSizeStack {
 public:
  FixedSizeStack(int capacity) {
    size_ = 0;
    stack_ = new int[capacity];
  }
  ~FixedSizeStack() { delete [] stack_; }
  void Push(int num) {
    stack_[size_] = num;
    size_++;
  }
  int Pop() {
    size_--;
    return stack_[size_];
  }
  int Size() const { return size_; }

  void Clear() { size_ = 0; }

  friend std::ostream& operator<<(std::ostream & out, const FixedSizeStack & st);

 private:
  int size_;
  int * stack_;
};

class MP_LAMP {
 public:
  /** 
   * n: granularity of tasks
   * m: number of random victim candidates
   * w: number of random victim tries
   * l: power of lifeline graph
   */
  MP_LAMP(int rank, int nu_proc, int n, bool n_is_ms, int w, int l, int m);

  ~MP_LAMP();

  // echo operation uses n-array tree (default: 3-ary)
  static const int k_echo_tree_branch;

  static int ComputeZ(int p, int l);

  // call this before mainloop ?
  // set dtd counters and misc.
  void Init();

  void CheckPoint();

  // clear unreceived MPI messages
  void ClearTasks();

  void InitTreeRequest();
  void SetTreeRequest();

  // read file, prepare database, broadcast to all procs
  // for h_==0
  void InitDatabaseRoot(std::istream & is1, std::istream &is2);
  void InitDatabaseRoot(std::istream & is1, int posnum);
  // other
  void InitDatabaseSub(bool pos);

  void InitDatabaseCommon();

  void Search();
  void MainLoop();

  void SearchStraw1(); // no workload distribution
  void MainLoopStraw1();
  void SearchStraw2(); // centralized queue
  void MainLoopStraw2();

  // Mattern 1987
  // bounded clock-counter version of TIME algorithm
  class DTD {
   public:
    DTD() :
        terminated_(false),
        count_(0),
        time_zone_(0),
        time_warp_(false),
        not_empty_(false)
    {
      accum_flag_ = new bool[k_echo_tree_branch];
    }

    ~DTD()
    {
      delete [] accum_flag_;
    }

    void Init() {
      terminated_ = false;

      count_ = 0;
      time_zone_ = 0;
      time_warp_ = false;
      not_empty_ = false;

      // for reduce
      //requesting_ = false;
      ClearAccumFlags();
      ClearReduceVars();
    }

    void CheckPoint() {
      terminated_ = false;
      //assert(requesting_ == false);
      for (int i=0;i<k_echo_tree_branch;i++) assert(accum_flag_[i] == false);
    }

    bool terminated_;

    int count_;
    int time_zone_;
    bool time_warp_;
    bool not_empty_;
    static const int kNuTimezone = 1024*1024;

    // on receiving basic message
    void UpdateTimeZone(int ts) {
      time_warp_ = time_warp_ || (ts == (time_zone_+1) % kNuTimezone);
    }

    void OnSend() { count_++; }
    void OnRecv() { count_--; }

    // on receiving control messages
    void ResetTimeWarp() { time_warp_ = false; }
    void IncTimeZone() {
      time_zone_++;
      time_zone_ = time_zone_ % kNuTimezone;
    }

    // for reduce
    //bool requesting_;
    bool * accum_flag_;

    int reduce_count_;
    bool reduce_time_warp_;
    bool reduce_not_empty_;

    void Reduce(int count, bool time_warp, bool not_empty) {
      reduce_count_ += count;
      reduce_time_warp_ = reduce_time_warp_ || time_warp;
      reduce_not_empty_ = reduce_not_empty_ || not_empty;
    }

    void ClearAccumFlags() {
      for (int i=0;i<k_echo_tree_branch;i++) accum_flag_[i] = false;
    }

    void ClearReduceVars() {
      reduce_count_ = 0;
      reduce_time_warp_ = false;
      reduce_not_empty_ = false;
    }

  };

  DTD dtd_;

  class Log {
   public:
    Log();
    ~Log();
    void Init();

    long long int idle_start_;

    void GatherLog(int nu_proc);

    //--------
    // periodic log

    void InitPeriodicLog();
    void StartPeriodicLog();
    void FinishPeriodicLog();

    long long int periodic_log_start_;
    long long int next_log_time_in_second_;

    // show this separately for phase_ 1 and 2
    void TakePeriodicLog(long long int capacity, int lambda, int phase);
    struct PeriodicLog_T {
      void Clear() {
        seconds_ = 0ll;
        capacity_ = 0ll;
        lambda_ = 0;
        phase_ = 0;
      }
      long long int seconds_;
      long long int capacity_;
      int lambda_;
      int phase_;
    };
    std::vector< PeriodicLog_T > plog_;

    void AggregatePLog(int nu_proc);
    int sec_max_;
    PeriodicLog_T * plog_buf_; // gather results for copying local plog
    PeriodicLog_T * plog_gather_buf_; // gather plog of all ranks here

    //--------
    // other log

    struct LogData {
      void Init();

      long long int search_start_time_;
      long long int search_finish_time_;

      long long int iprobe_num_;
      long long int iprobe_time_;
      long long int iprobe_time_max_;

      long long int iprobe_succ_num_;
      long long int iprobe_succ_time_;
      long long int iprobe_succ_time_max_;

      long long int iprobe_fail_num_;
      long long int iprobe_fail_time_;
      long long int iprobe_fail_time_max_;

      long long int probe_num_;
      long long int probe_time_;
      long long int probe_time_max_;

      long long int recv_num_;
      long long int recv_time_;
      long long int recv_time_max_;

      long long int bsend_num_;
      long long int bsend_time_;
      long long int bsend_time_max_;

      long long int bcast_num_;
      long long int bcast_time_;
      long long int bcast_time_max_;

      // long long int accum_task_time_;
      // long long int accum_task_num_;
      // long long int basic_task_time_;
      // long long int basic_task_num_;
      // long long int control_task_time_;
      // long long int control_task_num_;

      long long int dtd_phase_num_;
      double dtd_phase_per_sec_;
      long long int dtd_accum_phase_num_;
      double dtd_accum_phase_per_sec_;
      //long long int dtd_reply_num_;

      // long long int accum_phase_num_; // completed accum phase num
      // double accum_phase_per_sec_; // completed accum phase per second

      long long int lifeline_steal_num_;
      long long int lifeline_nodes_received_;
      long long int steal_num_;
      long long int nodes_received_;

      long long int lifeline_given_num_;
      long long int lifeline_nodes_given_;
      long long int given_num_;
      long long int nodes_given_;

      long long int process_node_num_;
      long long int process_node_time_;

      long long int preprocess_time_;

      long long int idle_time_;

      long long int pval_table_time_;

      long long int node_stack_max_itm_;
      long long int give_stack_max_itm_;

      long long int node_stack_max_cap_;
      long long int give_stack_max_cap_;

      long long int cleared_tasks_;

    };
    LogData d_;

    LogData * gather_buf_;
    LogData a_; // aggregated
    void Aggregate(int nu_proc);
  };

  const Database<uint64> & GetDatabase() { return *d_; }

  std::ostream & PrintDBInfo(std::ostream & out) const;

  std::ostream & PrintResults(std::ostream & out) const;
  // std::ostream & PrintSignificantMap(std::ostream & out) const;
  std::ostream & PrintSignificantSet(std::ostream & out) const;

  std::ostream & PrintLog(std::ostream & out) const;
  std::ostream & PrintAggrLog(std::ostream & out);

  std::ostream & PrintPLog(std::ostream & out);
  std::ostream & PrintAggrPLog(std::ostream & out);

 private:
  static const int k_int_max;
  // assuming (digits in long long int) > (bits of double mantissa)
  static const long long int k_cs_max;

  static const int k_probe_period;

  int * bsend_buffer_;

  // GLB variables

  int h_; // MPI Rank
  int p_; // total proc number

  int n_; // granularity of tasks
  bool n_is_ms_; // false: n_ is number of task, true: n_ is milli sec
  int w_; // number of random steal trials
  int m_; // number of random steal candidates

  int l_; // power of lifeline graph (length the hypercube edge)
  int z_; // dimension of lifeline (dimension of the hypercube)

  boost::mt19937 rng_; // use seed as rank
  boost::uniform_smallint<int> dst_p_;
  boost::uniform_smallint<int> dst_m_;
  boost::variate_generator<boost::mt19937&, boost::uniform_smallint<int> > rand_p_;
  boost::variate_generator<boost::mt19937&, boost::uniform_smallint<int> > rand_m_;

  int bcast_source_;
  int * bcast_targets_;

  // bool dtd_request_;
  // bool dtd_flag_[3];

  bool echo_waiting_; // new! [2015-10-05 22:23]
  //  bool bcast_requesting_; // new! [2015-10-01 13:49]

  //bool accum_requesting_;
  bool * accum_flag_;

  int * victims_; // proc id of random victims
  int * lifelines_; // proc id of lifeline buddies

  FixedSizeStack * thieves_; // max size == nu_proc_
  FixedSizeStack * lifeline_thieves_; // size == lifelines_ size + 3
  bool * lifelines_activated_;

  // proc id of lifeline thieves
  // bool * thieves_requests_;
  // bool * lifeline_requests_;

  // LAMP variables

  Database<uint64> * d_;
  LampGraph<uint64> * g_;
  VariableBitsetHelper<uint64> * bsh_; // bitset helper

  Log log_;
  Timer * timer_;

  // variables for LAMP
  int lambda_max_; // equals to maximum support of single item
  // initially set to 1. will be incremented to N if cs_thr[N] exceeded
  // in 1st phase, search will be pruned if (sup_num < lambda_)
  int lambda_;
  // todo: initial value can be set after checking all single item itemsets
  double sig_level_; // initially set to 1. set to 0.05 (FLAGS_a) / cs_thr_[global_sup_thr_-1]
  double * pmin_thr_; // pmin_thr[sup] == tbl.PMin(sup), maybe redundant
  long long int * cs_thr_; // cs_thr[sup] shows closed set num threshold

  void CheckCSThreshold();

  bool ExceedCsThr() const;
  int NextLambdaThr() const;
  void IncCsAccum(int sup_num);
  double GetInterimSigLevel(int lambda) const;

  // cs_accum_array is int array of 0..lambda_max_ (size lambda_max_+1)
  // cs_accum_array_[sup] shows closed set num with support higher than or equals to sup

  // note: allocate cs_accum_array_ as
  //   cs_accum_array_base_ = new int[lambda_max_+2]
  // and then
  //   cs_accum_array_ = cs_accum_array_base_+1
  // so that cs_accum_array_[-1] (== cs_accum_array_base_[0])
  // can be used for sending the timestamp counter
  // long long int * cs_accum_array_base_; // int array of -1..lambda_max_ (size lambda_max_+2)

  // 0: count, 1: time warp, 2: empty flag, 3--: array
  long long int * dtd_accum_array_base_; // int array of [-3..lambda_max_] (size lambda_max_+4)
  // -3: count, -2: time warp, -1: empty flag, 0--: array
  long long int * accum_array_;      // int array of [0...lambda_max_] (size lambda_max_+1)

  // 0: count, 1: time warp, 2: empty flag, 3--: array
  long long int * dtd_accum_recv_base_; // int array of [-3..lambda_max_] (size lambda_max_+4)
  // -3: count, -2: time warp, -1: empty flag, 0--: array
  long long int * accum_recv_;      // int array of [0...lambda_max_] (size lambda_max_+1)

  VariableLengthItemsetStack * node_stack_;
  // todo: prepare stack with no sup hist
  // 0: time zone, 1: is_lifeline
  VariableLengthItemsetStack * give_stack_;

  // periodic closed set count reduce.
  void Probe();
  bool processing_node_; // prevent termination while processing node

  // void ProbeAccumTask();
  // void ProbeBasicTask();
  // void ProbeControlTask();

  // first, send to random thief and then to lifeline theives
  // random thief has higher priority
  void Distribute();

  void Give(VariableLengthItemsetStack * st, int steal_num);

  void Deal();

  // send reject to remaining requests
  void Reject();

  // set this in steal, reset this in RecvGive and RecvReject
  // small difference from x10 implementation
  bool waiting_;

  // send steal requests
  // two phase, 1, random, 2, lifeline
  // if succeeds, break
  // note: original x10 GLB implementation always probes from i=0..z, is this OK?
  // how about prepare int steal_id_; and do steal_id_++/ steal_id_ %= z ?
  // note:
  // don't send multiple requests at once
  void Steal();
  // Steal needs change from x10 because of "bool waiting"
  // Steal sends one request each time it is called
  // there should be steal_state and counters c_r and c_l (random and lifeline)
  // c_r=0, c_l=0
  // 0, if (state == RANDOM && w > 0) goto 1, else goto 4
  // 1, send random request, c_r++
  // 2, wait for reject or give
  // 3, if c_r >= w, c_r = 0, set state LIFELINE and return
  //
  // 4, send lifeline request, c_l++
  // 5, wait for reject or give
  // 6, if c_l >= z, c_l = 0, set state RANDOM and return

  // steal with probe
  void Steal2();

  class StealState {
   public:
    enum {
      RANDOM = 0,
      LIFELINE,
    };

    StealState(int w, int z)
        : w_(w), z_(z)
    {
      Init();
    }

    int w_;
    int z_;

    void Init() {
      steal_phase_started_ = true;

      random_counter_ = 0;
      lifeline_counter_ = 0;
      next_lifeline_victim_ = 0;

      if (w_) state_ = RANDOM;
      else state_ = LIFELINE;
      requesting_ = false;
    }      

    void Finish() {
      steal_phase_started_ = false;

      random_counter_ = 0;
      lifeline_counter_ = 0;

      next_lifeline_victim_ = 0;
      // this is currently needed to resset lifeline victim
      // when lifeline has -1 (if p_ is not powers of 2)

      if (w_) state_ = RANDOM;
      else state_ = LIFELINE;
      requesting_ = false;
    }

    void ResetCounters() {
      random_counter_ = 0;
      lifeline_counter_ = 0;
    }

    void IncRandomCount() {
      random_counter_++;
      if (random_counter_ >= w_) random_counter_ = 0;
    }
    int RandomCount() const { return random_counter_; }

    void IncLifelineCount() { lifeline_counter_++; }
    int LifelineCounter() const { return lifeline_counter_; }

    void IncLifelineVictim() {
      next_lifeline_victim_++;
      if (next_lifeline_victim_ >= z_) next_lifeline_victim_ = 0;
    }
    int LifelineVictim() const { return next_lifeline_victim_; }

    void SetState(int state) { state_ = state; }
    int State() const {
      if (w_==0) return LIFELINE;
      else return state_;
    }

    bool Requesting() const { return requesting_; }
    void SetRequesting(bool flag=true) { requesting_ = flag; }
    void ResetRequesting() { requesting_ = false; }

    bool StealStarted() const { return steal_phase_started_; }
    void SetStealStart() { steal_phase_started_ = true; }
    void ResetStealStart() { steal_phase_started_ = true; }

    bool steal_phase_started_;

    int random_counter_; // random victim counter
    int lifeline_counter_; // lifeline victim counter
    int next_lifeline_victim_; // next lifeline victim range [0..z-1]

    int state_; // 0: random steal, 1: lifeline steal
    bool requesting_; // true if requesting. become false if rejected or given
    // after sending request, wait until give or reject
  };

  // will be false if w random steal and all lifeline steal finished
  // will be true if RecvGive
  StealState stealer_;

  //--------
  // control

  // 0: count, 1: time warp flag, 2: empty flag
  void SendDTDRequest();
  void RecvDTDRequest(int src);

  bool DTDReplyReady() const;
  void DTDCheck();

  // 0: count, 1: time warp flag, 2: empty flag
  void SendDTDReply();
  void RecvDTDReply(int src);

  bool DTDAccumReady() const;

  // 0: count, 1: time warp flag, 2: empty flag, 3--: data
  void SendDTDAccumRequest();
  void RecvDTDAccumRequest(int src);

  // 0: count, 1: time warp flag, 2: empty flag, 3--: data
  void SendDTDAccumReply();
  void RecvDTDAccumReply(int src);

  void SendBcastFinish();
  void RecvBcastFinish(int src);

  //--------

  int phase_; // 1, 2, 3

  //--------
  // basic

  // send recv functions
  void SendRequest(int dst, int is_lifeline); // for random thieves, is_lifeline = -1
  void RecvRequest(int src);

  // 0: time zone, 1: is_lifeline
  void SendReject(int dst);
  void RecvReject(int src);

  // 1: time zone
  void SendGive(VariableLengthItemsetStack * st, int dst, int is_lifeline);

  // sets lifelines_activated_ = false
  // lifelines_activated_ becomes false only in this case (reject does NOT)
  void RecvGive(int src, MPI_Status status);

  void SendLambda(int lambda);
  void RecvLambda(int src);

  // 0: time zone
  // search for depth 1 and get initial lambda
  void PreProcessRootNode();

  // provide int n and bool n_is_ms_ (if n_is_ms_==false, it shows number of nodes)
  bool ProcessNode(int n);
  bool CheckProcessNodeEnd(int n, bool n_is_ms, int processed, long long int start_time);

  bool ProcessNodeStraw1(int n);
  bool ProcessNodeStraw2(int n);

  int itemset_buf_[VariableLengthItemsetStack::kMaxItemsPerSet];
  uint64 * sup_buf_, * child_sup_buf_;

  /* sorting itemset
     1, ascending order of pval
     2, descending order of item numbers
     3, dictionary order of items
   */
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

  // hoge;
  // modify lamp.cc also

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
          throw std::runtime_error("identical duplicate itemsets found at operator()");
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
          throw std::runtime_error("identical duplicate itemsets found at compare");
          return false;
        }
      }
      return false;
    }

    const VariableLengthItemsetStack & ss_;
  };

  /** record frequent enough itemset during 2nd phase (candidates for significant itemsets)
   *  (if third phase is enabled) */
  VariableLengthItemsetStack * freq_stack_; // record freq itemsets
  // record significant set per process. gathered in the final phase
  // checkme: not neede because we can use significant_set_ for both purpose?
  // std::set<SignificantSetResult, sigset_compare> significant_set_per_proc_;
  // std::multimap< double, int * > freq_map_; // record (pval, *itemsets)
  /** record final significant itemsets in 3rd phase */
  VariableLengthItemsetStack * final_significant_stack_;

  sigset_compare * freq_comp_;
  sigset_compare * final_sigset_comp_;

  /** holds final significant set result. store itemset pointers and sort only the pointers */
  std::set<SignificantSetResult, sigset_compare> * freq_set_;
  std::set<SignificantSetResult, sigset_compare> * final_significant_set_;

  // add pos_sup_num info in significant map
  // add stable sort
  // std::multimap< double, int * > significant_map_;

  /** pushing item, assuming items in itemsets are sorted */
  int * PushItemsetNoSort(VariableLengthItemsetStack * ss, int * itemset) {
    ss->PushPre();
    int * item = ss->Top();
    ss->CopyItem(itemset, item);
    ss->PushPostNoSort();
    return item;
  }

  /** record itemset if needed, following the SigsetRecordMode */
  void RecordFrequentItemset(double pval, double sig_level,
                             int sup_num, int pos_sup_num, int * itemset);

  struct SigsetRecordMode {
    enum SigsetRecordModeType {
      NORMAL = 0, // default, record all significant patterns
      AT_MOST_N, // record at most N significant patterns
      AT_LEAST_M, // show at least M patterns regardless of significance
      M_TO_N, // show at least M and at most N
    };
  };

  int sigset_record_mode_;
  long long int infrequent_itemset_num_;

  long long int total_expand_num_;
  long long int expand_num_;
  long long int closed_set_num_;

  //--------
  // third phase

  // void ProbeThirdPhaseTask();

  bool AccumCountReady() const;

  void SendResultRequest();
  void RecvResultRequest(int src);

  void SendResultReply();
  void RecvResultReply(int src, MPI_Status status);

  void ExtractSignificantSet();

  // insert pointer into significant_map_ (do not sort the stack itself)
  void SortSignificantSets();

  /** calculate support from itemset */
  void CalculatePval(const VariableLengthItemsetStack * ss, const int * itemset,
                     int * sup_num, int * pos_sup_num, double * pval) {
    bsh_->Set(sup_buf_);
    {
      int n = ss->GetItemNum(itemset);
      for (int i=0;i<n;i++) {
        int item = ss->GetNthItem(itemset, i);
        bsh_->And(d_->NthData(item), sup_buf_);
      }
    }

    *sup_num = bsh_->Count(sup_buf_);
    *pos_sup_num = bsh_->AndCount(d_->PosNeg(), sup_buf_);
    *pval = d_->PVal(*sup_num, *pos_sup_num);
  }

  //--------
  // for printing results

  long long int final_closed_set_num_;
  int final_support_;
  double final_sig_level_;

  // true if bcast_targets_ are all -1
  bool IsLeaf() const;

  // return flag. if (flag), it is ready to receive
  int CallIprobe(MPI_Status * status, int * count, int * src);

  int CallRecv(void * buffer, int count, MPI_Datatype type,
               int src, int tag, MPI_Status * status);
  int CallBsend(void * buffer, int count_int, MPI_Datatype type,
                int dest, int tag);

  int CallBcast(void * buffer, int data_count, MPI_Datatype type);
  // todo: implement call reduce, call gather

  //--------
  // for debug

  void SetLogFileName();

  std::string log_file_name_;
  std::ofstream lfs_;

  /** Null stream.
      This file stream will never be opened and acts as a null stream
      immitating Fuego's SgDebug(). */
  static std::ofstream null_stream_;
  std::ostream& D(int level, bool show_phase=true);

  // static std::stringstream alert_ss;
  // std::ostream& Alert();

  //--------
  // testing [2015-10-01 13:54]
  bool last_bcast_was_dtd_;
  
};

} // namespace lamp_search

#endif // _LAMP_SEARCH_MP_DFS_H_

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
