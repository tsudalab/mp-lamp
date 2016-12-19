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

#ifndef _LAMP_SEARCH_VARIABLE_LENGTH_ITEMSET_H_
#define _LAMP_SEARCH_VARIABLE_LENGTH_ITEMSET_H_

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>

namespace lamp_search {

/** variable length itemset data
 *  packed to a continuous memory region */
class VariableLengthItemsetStack {
 public:

  static const int TIMESTAMP = 0;
  static const int FLAG = 1; // used for lifeline flag
  static const int SENTINEL = 2;

  // data structure for one itemset
  // 0: [NUM] = (-1) * (number of items + 1)
  //    negated for indicating start of an itemset. see GetItemNum for details
  // 1: [SUP] number of support
  // 2: [ITM] points to the first item in the itemset
  static const int NUM = 0;
  static const int SUP = 1;
  static const int ITM = 2;

  // todo: add check
  static const int kMaxItemsPerSet = 1024 * 1024 - ITM; // max 4MB

  // todo: add check. size must be greater than kMaxItemsPerSet
  // should be like 256 * kMaxItemsPerSet or something
  VariableLengthItemsetStack(std::size_t size);

  // VariableLengthItemsetStack(std::size_t size, int sup_max);

  ~VariableLengthItemsetStack();

  // move top_ and inc nu_itemset_;
  void PushPre();
  // update 
  void PushPost();
  void PushPostNoSort();
  void Pop();

  bool SetItemNum(int * index, int num);
  bool SetSup(int * index, int sup);
  bool SetOneItem(int * index, int item_num, int item);

  bool PushOneItem(int item);
  bool IncItemNum(int * index);
  bool DecItemNum(int * index);
  void Clear();

  static int GetItemNum(const int * index);
  static int GetSup(const int * index);
  static int * GetItemArray(int * index);
  static const int * GetItemArray(const int * index);

  int GetNthItem(const int * index, std::size_t n) const;
  void CopyItem(const int * src, int * dst);
  // return next item, if not next item, return NULL
  int * NextItemset(int * index) const;
  // todo: test
  bool Exist(const int * index, int item) const;

  // return ptr to stack_
  int * Stack() const { return stack_; }

  // return ptr to top item
  int * Top() const;
  int * FirstItemset() const;
  
  void RemoveOneItemset();

  void SortOneSet(int * index);
  void SortTop();

  // void IncSupHistogram(int * index);
  // void DecSupHistogram(int * index);

  // give half of entries to dst (dst should be empty)
  // return given itemset num
  int Split(VariableLengthItemsetStack * dst);
  // merge with entries in another stack
  bool Merge(VariableLengthItemsetStack * src);
  // merge with entries in array
  // stack points to first itemset, size is equal to used capacity - 2
  bool MergeStack(int * stack, int size);

  // used to discard top itemset
  // mainly used when removing top after used up detected ???
  // will be called instead of PushPost
  void DiscardTop();

  // discard elements with less than support threshold
  void GC(int sup_threshold);

  bool Empty() const { return nu_itemset_ == 0; }

  // when to call this? every time adding item?
  bool Full() const;

  int TotalCapacity() const { return total_capacity_; }
  int UsedCapacity() const { return used_capacity_; }

  int NuItemset() const { return nu_itemset_; }

  // functions for mp-lamp
  void SetTimestamp(int ts) { stack_[TIMESTAMP] = ts; }
  int Timestamp() { return stack_[TIMESTAMP]; }

  void SetFlag(int flg) { stack_[FLAG] = flg; }
  int Flag() { return stack_[FLAG]; }

  // idea
  // push procedure
  // PushPre() (just increase nu_itemset and return next pointer)
  // ptr = top()
  // set sup, push items, set itemnum, update others
  // PushPost()

  // pop procedure
  // Pop()

  std::ostream& PrintAll(std::ostream & out) const;

  std::ostream& Print(std::ostream & out, const int * index) const;

  std::ostream& Print(std::ostream & out, const std::vector< std::string > * item_names,
                      const int * index) const;

 private:
  int * stack_;
  // note: top_ points to the head of the last item
  //       stack_[used_capacity_] points to the head of the next item
  int * top_;

  // note: is this enough?
  int total_capacity_;
  int used_capacity_;

  std::size_t nu_itemset_;

  // int sup_max_;
  // int * sup_hist_; // sup histogram
};

} // namespace lamp_search

#endif // _LAMP_SEARCH_VARIABLE_LENGTH_ITEMSET_H_
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
