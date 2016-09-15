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

#ifndef _LAMP_SEARCH_HASH_H_
#define _LAMP_SEARCH_HASH_H_

#include <cstdio>
#include <cstddef>
#include <cstring>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <vector>

#include "utils.h"

#include "gflags/gflags.h"

namespace lamp_search {

// template parameter ELEMENT, hash value type

template<typename VAL>
struct HashEntry {
 public:
  uint64 code;
  HashEntry<VAL> * next;
  VAL sup[1];
 private:
};

template<typename VAL>
class Hash {
 private:
  /** number of elements checked for replace. must be greater than 2. */
  static const int kReplaceEPS = 5;
  /** element with small number of visits will be replaced.  */
  static const int kReplaceVisitNum = 3; // what number is appropriate?

  bool single_byte_;

  HashEntry<VAL> * element_table_;
  HashEntry<VAL> ** table_;
  HashEntry<VAL> * free_list_;

  int table_size_;
  bool use_gc_;
  int num_used_;
  int smalltree_gc_threshold_;
  bool is_used_up_;

  //CArray<uint64, Game::kMaxBranch*3> random_table_;
  //static uint64 random_table_[Game::kMaxBranch*3];

 public:
  Hash(std::size_t nu_item, int size, bool use_gc, bool size_in_bytes);
  ~Hash();

  int ComputeHashKey(uint64 code) {
    int key = code % table_size_;
    return key;
  }

  int Size() { return table_size_; }
  int NumUsed() { return num_used_; }

  void Init();

  HashEntry<VAL> * Lookup(uint64 code);
  HashEntry<VAL> * Allocate(uint64 code);
  HashEntry<VAL> * Replace(int key);

  /** dump all element (including unused black elements) */
  void DumpAll(std::ostream & out);

  /** dump element with visit >= size */
  void DumpUsed(std::ostream & out, std::size_t size);
};

template<typename VAL>
Hash<VAL>::Hash(std::size_t nu_item, int size, bool use_gc, bool size_in_bytes)
    : table_size_(size), use_gc_(use_gc) {
  smalltree_gc_threshold_ = 0;
  // hoge
  // (HashEntry<VAL>*)hoge
  // 
  // element_table_ = new char[table_size_];
  // table_ = new char*[table_size_];

  // todo: implement GC
  // just discard too deep entries?
  assert(0);
}

template<typename VAL>
Hash<VAL>::~Hash() {
  assert(element_table_);
  assert(table_);
  // delete [] element_table_;
  // delete [] table_;
}

template<typename VAL>
void Hash<VAL>::Init() {
  is_used_up_ = false;
  num_used_ = 0;
  smalltree_gc_threshold_ = 0;

  //std::memset(static_cast<void *>(table_), 0, table_size_ * sizeof(ELEMENT *));
  //std::memset(static_cast<void *>(element_table_), 0, table_size_ * sizeof(ELEMENT));
  free_list_ = (HashEntry<VAL> *)NULL;
}

template<typename VAL>
HashEntry<VAL> *Hash<VAL>::Replace(int key) {
  int index = key;
  // victim1 : smallest visits
  HashEntry<VAL> * victimElement = 0, * victimPrevElement = 0;
  int victimIndex = -1;

  // compare visits and value upto kReplaceEPS entries

  // todo: implement
  assert(0);

  for (int i=0;i<kReplaceEPS;) {
    for (HashEntry<VAL> * element = table_[index], *prevElement = 0; element;
         prevElement = element, element = element->getNext()) {
      if (element->isSearched())
        continue;
      i++;
      if (i >= kReplaceEPS) {
        assert (victimElement);
        assert (victimIndex != -1);
        if (!victimPrevElement)
          table_[victimIndex] = victimElement->getNext();
        else victimPrevElement->setNext(victimElement->getNext());
        return victimElement;
      }
      if (victimElement == 0) {
        victimElement = element;
        victimPrevElement = prevElement;
        victimIndex = index;
        if (element->getImportance() < kReplaceVisitNum) {
          assert(victimElement);
          assert(victimIndex != -1);
          if (!victimPrevElement)
            table_[victimIndex] = victimElement->getNext();
          else
            victimPrevElement->setNext(victimElement->getNext());
          return victimElement;
        }
      }
      else {
        if (element->getImportance() < victimElement->getImportance() ) {
          victimElement = element; victimPrevElement = prevElement;
          victimIndex = index;
          if (element->getImportance() < kReplaceVisitNum) {
            assert(victimElement);
            assert(victimIndex != -1);
            if (!victimPrevElement)
              table_[victimIndex] = victimElement->getNext();
            else
              victimPrevElement->setNext(victimElement->getNext());
            return victimElement;
          }
        }
      }
      assert(victimElement);
    }
    index++;
    index %= table_size_;
  }
  assert(victimElement); assert(victimIndex != -1);
  if (!victimPrevElement) table_[victimIndex] = victimElement->getNext();
  else victimPrevElement->setNext(victimElement->getNext());
  return victimElement;
}

template<typename VAL>
HashEntry<VAL> *Hash<VAL>::Allocate(uint64 code) {
  int key = ComputeHashKey(code);
  HashEntry<VAL> * element;
  for (element = table_[key]; element; element = element->getNext()) {
    if (element->getHashcode() == code) {
      return element;
    }
  }
  if (!is_used_up_) {
    element = &element_table_[num_used_++];
    if (num_used_ == table_size_) is_used_up_ = true;
  }
  else {
    if (!use_gc_) {
      element = Replace(key);
    }
    else {
      assert(0);
      if (free_list_ == 0) {
        // fixme: implement smalltreegc
        // smalltregc();
      }
      element = free_list_;
      free_list_ = free_list_->getNext();
      num_used_++;
    }
  }
  assert (element);

  element->setNext(table_[key]);
  table_[key] = element;
  return element;
}

template<typename VAL>
HashEntry<VAL> * Hash<VAL>::Lookup(uint64 code) {
  // todo:implement
  assert(0);
  int key = ComputeHashKey(code);
  for (HashEntry<VAL> * element = table_[key]; element; element = element->getNext()) {
    if (element->getHashcode() == code) return element;
  }
  return 0;
}

template<typename VAL>
void Hash<VAL>::DumpAll(std::ostream & out) {
  out << "====" << std::endl;
  for(int i=0;i<table_size_;i++) {
    out << "key=" << i << ":";
    out << ":pointer=" << std::hex << table_[i] << ":";
    if (table_[i]) out << table_[i];
    else out << std::endl;
  }
}

template<typename VAL>
void Hash<VAL>::DumpUsed(std::ostream & out, std::size_t size) {
  for(int i=0;i<table_size_;i++) {
    HashEntry<VAL> * elm = table_[i];
    while(elm) {
      if ( (std::size_t)elm->getImportance() >= size) {
        out << "key=" << std::setw(6) << i;
        out << ":pointer=" << std::hex << elm << ":" << std::dec;
        out << *elm;
      }
      elm = elm->getNext();
    }
  }
}

} // namespace lamp_search

#endif // _LAMP_SEARCH_HASH_H_
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
