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
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>

#include "utils.h"
#include "variable_length_itemset.h"

namespace lamp_search {

VariableLengthItemsetStack::VariableLengthItemsetStack(std::size_t size)
    : stack_ (NULL), top_ (NULL)
{
  //assert(size >= kMaxItemsPerSet + SENTINEL);

  total_capacity_ = size;
  stack_ = new int[size];
  stack_[TIMESTAMP] = 0; // timestamp
  stack_[SENTINEL] = -1; // sentinel
  top_ = &(stack_[SENTINEL]); // is this correct?
  nu_itemset_ = 0;
  used_capacity_ = SENTINEL+1; // including sentinel

  // sup_hist_ = new int[sup_max_];
  // for (int i=0;i<sup_max_;i++)
  //   sup_hist_[i] = 0;
}

VariableLengthItemsetStack::~VariableLengthItemsetStack() {
  if (stack_) delete [] stack_;
  // if (sup_hist_) delete [] sup_hist_;
}

void VariableLengthItemsetStack::PushPre() {
  top_ = &(stack_[used_capacity_]);
  SetItemNum(top_, 0); // clear item num
  SetSup(top_, 0); // clear sup
  nu_itemset_++;
}

void VariableLengthItemsetStack::PushPost() {
  int * index = Top();
  used_capacity_ += ITM + GetItemNum(index);
  SortOneSet(index);
  // IncSupHistogram(index);
}

void VariableLengthItemsetStack::PushPostNoSort() {
  int * index = Top();
  used_capacity_ += ITM + GetItemNum(index);
  // IncSupHistogram(index);
}

void VariableLengthItemsetStack::SortTop() {
  int * index = Top();
  SortOneSet(index);
}

void VariableLengthItemsetStack::SortOneSet(int * index) {
  int * array = GetItemArray(index);
  int num = GetItemNum(index);
  std::sort( array, array + num );
}

void VariableLengthItemsetStack::Pop() {
  RemoveOneItemset();
}

bool VariableLengthItemsetStack::SetItemNum(int * index, int num) {
  if (used_capacity_ + NUM >= total_capacity_) return false;
  index[NUM] = (-1) * (num + 1);
  return false;
}

int VariableLengthItemsetStack::GetItemNum(const int * index) {
  return ((-1) * index[NUM] -1 );
}

bool VariableLengthItemsetStack::SetSup(int * index, int sup) {
  if (used_capacity_ + SUP >= total_capacity_) return false;
  index[SUP] = sup;
  return true;
}

int VariableLengthItemsetStack::GetSup(const int * index) {
  return index[SUP];
}

bool VariableLengthItemsetStack::SetOneItem(int * index, int item_num, int item) {
  if (used_capacity_ + item_num + ITM >= total_capacity_) return false;
  index[ITM + item_num] = item;
  return true;
}

bool VariableLengthItemsetStack::PushOneItem(int item) {
  int * index = Top();
  int item_num = GetItemNum(index);
  index[ITM + item_num] = item;
  bool check = IncItemNum(index);
  if (!check) return false;
  return true;
}

bool VariableLengthItemsetStack::IncItemNum(int * index) {
  int num = GetItemNum(index);
  num++;
  if (used_capacity_ + num + ITM > total_capacity_) return false;
  SetItemNum(index, num);
  return true;
}

bool VariableLengthItemsetStack::DecItemNum(int * index) {
  int num = GetItemNum(index);
  num--;
  if (num < 0) return false;
  SetItemNum(index, num);
  return false;
}

int * VariableLengthItemsetStack::GetItemArray(int * index) {
  return &(index[ITM]);
}

const int * VariableLengthItemsetStack::GetItemArray(const int * index) {
  return &(index[ITM]);
}

int VariableLengthItemsetStack::GetNthItem(const int * index, std::size_t n) const {
  return GetItemArray(index)[n];
}

void VariableLengthItemsetStack::CopyItem(const int * src, int * dst) {
  int num = GetItemNum(src);

  dst[NUM] = src[NUM];
  dst[SUP] = src[SUP];
  for (int i=0 ; i<num ; i++)
    dst[i+ITM] = src[i+ITM];
}

int * VariableLengthItemsetStack::NextItemset(int * index) const {
  if (index == top_) return NULL;
  int num = GetItemNum(index);
  return index + num + ITM;
}

bool VariableLengthItemsetStack::Exist(const int * index, int item) const {
  int n = GetItemNum(index);

  for (int i=0;i<n;i++) {
    if (index[ITM+i] == item) return true;
    else if (index[ITM+i] > item) return false; // assert sorted
  }

  return false;
}

int * VariableLengthItemsetStack::Top() const {
  return top_;
}

int * VariableLengthItemsetStack::FirstItemset() const {
  if (NuItemset() == 0) return NULL;
  return &(stack_[SENTINEL+1]);
}

void VariableLengthItemsetStack::RemoveOneItemset() {
  // todo: add exception or exit with error message
  assert(top_ != stack_);
  int * index = Top();
  // DecSupHistogram(index);
  used_capacity_ -= ITM + GetItemNum(index);
  nu_itemset_--;
  top_--;
  while(*top_ >=0) top_--;
  assert( Top()[NUM] < 0);
}

// void VariableLengthItemsetStack::IncSupHistogram(int * index) {
//   sup_hist_[ GetSup(index) ]++;
// }

// void VariableLengthItemsetStack::DecSupHistogram(int * index) {
//   sup_hist_[ GetSup(index) ]--;
// }

int VariableLengthItemsetStack::Split(VariableLengthItemsetStack * dst) {
  if (NuItemset() == 0) return 0;

  int * first = FirstItemset();
  int * next_top = first;
  int * head = NextItemset(first);
  int * give = head;
  if (give == NULL) return 0;
  int * keep = NextItemset(give);
  int * next;
  if (keep == NULL) next = NULL;
  else next = NextItemset(keep);

  int * dst_head = dst->FirstItemset();
  assert(dst->NuItemset() == 0);

  int given_num = 0;

  while (true) {
    dst->PushPre();
    dst_head = dst->Top();
    CopyItem(give, dst_head);
    dst->PushPostNoSort();
    given_num++;

    // DecSupHistogram(give);
    used_capacity_ -= ITM + GetItemNum(give);
    nu_itemset_--;

    if (keep == NULL) break;

    next_top = head;
    CopyItem(keep, head);
    head = NextItemset(head);

    if (next == NULL) break;

    give = next;
    keep = NextItemset(give);
    if (keep == NULL) next = NULL;
    else next = NextItemset(keep);
  }
  top_ = next_top;
  // update
  // top_, used_capacity_, nu_itemset_, sup_hist_
  return given_num;
}

bool VariableLengthItemsetStack::Merge(VariableLengthItemsetStack * src) {
  // check capacity
  if (this->UsedCapacity() + (src->UsedCapacity()-SENTINEL-1) > this->TotalCapacity())
    return false;
  
  int * next_top = this->top_;
  int di=this->used_capacity_;
  // note: si starts from (SENTINEL+1) to skip timestamp and sentinel
  for (int si=SENTINEL+1;si<src->used_capacity_;si++) {
    this->stack_[di] = src->stack_[si];
    if (this->stack_[di] < 0) next_top = &(this->stack_[di]);
    di++;
  }

  this->top_ = next_top;
  this->used_capacity_ += src->UsedCapacity() - SENTINEL-1;
  this->nu_itemset_ += src->nu_itemset_;

  // for (int i=0;i<sup_max_;i++)
  //   this->sup_hist_[i] += src->sup_hist_[i];
  
  return true;

  // todo: limit maximum amount of give to (stack size) / z ???
}

bool VariableLengthItemsetStack::MergeStack(int * src_st, int size) {
  // check capacity
  if ( (this->UsedCapacity() + size) > this->TotalCapacity() )
    return false;
  
  int * next_top = this->top_;
  int di=this->used_capacity_;

  for (int si=0;si<size;si++) {
    this->stack_[di] = src_st[si];

    // update nu_itemset_ and sup_hist_
    if (src_st[si] < 0) {
      next_top = &(this->stack_[di]);
      // IncSupHistogram(&(src_st[si]));
      nu_itemset_++;
    }

    di++;
  }

  this->top_ = next_top;
  this->used_capacity_ += size;

  return true;

  // todo: limit maximum amount of give to (stack size) / z ???
}

void VariableLengthItemsetStack::Clear() {
  nu_itemset_ = 0;
  used_capacity_ = SENTINEL+1; // including sentinel
  top_ = &(stack_[SENTINEL]); // is this correct?

  // for (int i=0;i<sup_max_;i++)
  //   sup_hist_[i] = 0;
}

std::ostream& VariableLengthItemsetStack::PrintAll(std::ostream & out) const {
  std::stringstream s;

  int * p = FirstItemset();
  while (p != NULL) {
    Print(s, p);
    p = NextItemset(p);
  }

  out << s.str() << std::flush;
  return out;
}

std::ostream& VariableLengthItemsetStack::Print(std::ostream & out,
                                                const int * index) const {
  std::stringstream s;

  int n = GetItemNum(index);

  s << "num=" << n;
  s << "\tsup=" << GetSup(index);
  s << "\titems:";
  for (int i=0;i<n;i++)
    s << " " << index[ITM+i];
  s << std::endl;
  
  out << s.str() << std::flush;
  return out;
}

std::ostream& VariableLengthItemsetStack::Print(std::ostream & out,
                                                const std::vector< std::string > * item_names,
                                                const int * index) const {
  std::stringstream s;

  int n = GetItemNum(index);

  // s << "num=" << n;
  // s << "\tsup=" << GetSup(index);
  // s << "\titems:";
  s << "\t" << n;
  if (item_names != NULL) {
    for (int i=0;i<n;i++)
      s << "\t" << (*item_names)[ index[ITM+i] ];
    s << std::endl;
  } else {
    for (int i=0;i<n;i++)
      s << "\t" << index[ITM+i];
    s << std::endl;
  }
  
  out << s.str() << std::flush;
  return out;
}

} // namespace lamp_search
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
