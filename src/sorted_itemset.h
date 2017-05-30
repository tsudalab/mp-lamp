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

#ifndef _LAMP_SEARCH_SORTED_ITEMSET_H_
#define _LAMP_SEARCH_SORTED_ITEMSET_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <list>
#include <limits>

#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>

#include "utils.h"

namespace lamp_search {

class ItemInfo {
 public:
  int id;
  double pmin;
  double pval;

  int sup;
  int pos_sup;

  bool operator<(const ItemInfo& rhs) const {
    return this->pmin < rhs.pmin;
  }
  friend std::ostream& operator<<(std::ostream & out, const ItemInfo & i);
};

class SortedItemSet {
 public:
  SortedItemSet() : num_(0) {}

  typedef boost::array<int, kMaxItemSet> ItemSetType;

  SortedItemSet(const std::vector<int> & v) {
    num_ = v.size();
    int s=0;
    for (std::vector<int>::const_iterator it=v.begin();it!=v.end();++it) {
      set_[s++] = *it;
    }
    assert(s==num_);
  }

  SortedItemSet(const ItemSetType & a) {
    int s=0;
    for (ItemSetType::const_iterator it=a.begin();it!=a.end();++it) {
      if (*it >= std::numeric_limits<int>::max()) break;
      set_[s++] = *it;
    }
    num_=s;
  }

  SortedItemSet(const int * a) {
    int s=0;
    for (int i=0;i<kMaxItemSet;i++) {
      if (a[i] >= std::numeric_limits<int>::max()) break;
      set_[s++] = a[i];
    }
    num_=s;
  }

  bool operator==(const SortedItemSet & rhs) const {
    if (this->num_ != rhs.num_) return false;
    bool ret = true;
    for (int d=0 ; d<num_ ; d++) {
      if (this->set_[d] != rhs.set_[d]) {
        ret = false;
        break;
      }
    }
    return ret;
  }

  SortedItemSet & operator=(const SortedItemSet & rhs) {
    // std::cout << "d: " << rhs << std::endl;
    // std::cout << "d: " << rhs.num_;
    this->num_ = rhs.num_;
    for (int d=0 ; d<this->num_ ; d++) this->set_[d] = rhs.set_[d];
    return *this;
  }

  // todo: test Merge
  void Merge(const SortedItemSet & ss) {
    if (this->num_ + ss.num_ <= kMaxItemSet) {
      int ti=this->num_;
      for (int si=0; si<ss.num_; si++) {
        this->set_[ti] = ss[si];
        ti++;
      }
      this->num_ += ss.num_;
      std::sort(set_.begin(), set_.begin()+num_);
    }      
    else {
      int ti=this->num_;
      for (int si=0; (this->num_ + si)<kMaxItemSet; si++) {
        this->set_[ti] = ss[si];
        ti++;
      }
      this->num_ = kMaxItemSet;
      std::sort(set_.begin(), set_.begin()+num_);
    }
  }

  int & operator[](std::size_t i) { return set_[i]; }
  const int & operator[](std::size_t i) const { return set_[i]; }

  void Push(int id) {
    //std::cout << "d2: num_=" << num_ << std::endl;
    assert(num_<kMaxItemSet);
    set_[num_++]=id;
    std::sort(set_.begin(), set_.begin()+num_);
  }

  bool Exists(int id) const {
    bool ret = false;
    for (int d=0 ; d<num_ ; d++) {
      if (set_[d] == id) {
        ret = true;
        break;
      }
    }
    return ret;
  }

  bool Remove(int id) {
    bool ret = false;
    for (int d=0 ; d<num_ ; d++) {
      if (set_[d] == id) {
        num_--; assert(num_ >= 0);
        for (int dd=d;dd<num_;dd++) set_[dd] = set_[dd+1];
        ret = true;
        break;
      }
    }
    return ret;
  }

  // void Pop() {
  //   num_--;
  //   assert(num_ >= 0);
  // }

  void Clear() { num_=0; }

  bool Empty() const { return num_ <= 0; }
  bool Full()  const { return num_ == kMaxItemSet; }

  int Size() const { return num_; }

  friend std::ostream& operator<<(std::ostream & out, const SortedItemSet & ss);

  typedef ItemSetType::iterator iterator;
  typedef ItemSetType::const_iterator const_iterator;
  iterator begin() { return set_.begin(); }
  iterator end() { return set_.begin() + num_; }
  const_iterator begin() const { return set_.begin(); }
  const_iterator end() const { return set_.begin() + num_; }

 private:
  ItemSetType set_;
  int num_;
};

//==============================================================================
std::ostream & operator<<(std::ostream & out, const ItemInfo & i);

std::ostream & operator<<(std::ostream & out, const SortedItemSet & ss);

} // namespace lamp_search

#endif // _LAMP_SEARCH_SORTED_ITEMSET_H_

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
