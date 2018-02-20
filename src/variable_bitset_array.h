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

#ifndef _LAMP_SEARCH_VARIABLE_BITSET_ARRAY_H_
#define _LAMP_SEARCH_VARIABLE_BITSET_ARRAY_H_

#include <iostream>

#include <assert.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <limits>
//#include <climits> // for CHAR_BIT ?

#include <boost/array.hpp>
#include <boost/static_assert.hpp>

#include "utils.h"

namespace lamp_search {

namespace {

// from boost/dynamic_bitset
// for static_asserts
template <typename T>
struct allowed_block_type {
  enum { value = T(-1) > T(0) }; // ensure T has no sign
};

template <>
struct allowed_block_type<bool> {
  enum { value = false };
};

} // namespace anonymous

template<typename Block>
class VariableBitsetTraits {
 public:
  static const std::size_t bits_per_block;
  static const std::size_t bytes_per_block;
  // static const std::size_t block_index_shift;
  // static const std::size_t block_index_mask;

  // static const std::size_t bytes_per_block = sizeof(Block);
  // static const std::size_t bits_per_block = std::numeric_limits<Block>::digits;

  BOOST_STATIC_ASSERT((bool)allowed_block_type<Block>::value); // Block is unsigned

  // enum { val = bits_per_block && !(bits_per_block & (bits_per_block - 1)) };
  // BOOST_STATIC_ASSERT(val); // bits_per_block is powers of 2

  // bits_per_block is powers of 2
  BOOST_STATIC_ASSERT( std::numeric_limits<Block>::digits != 0);
  BOOST_STATIC_ASSERT(
      !(std::numeric_limits<Block>::digits & (std::numeric_limits<Block>::digits - 1))
                      );

  static std::size_t calc_msb_position(Block num) {
    // std::cout << "calc_msb_position() called" << std::endl;
    assert(num != 0);
    std::size_t pos = 0;
    while (true) {
      if (num == 1) break;
      num >>= 1;
      pos++;
    }
    return pos;
  }

  static Block calc_msb_mask(std::size_t pos) {
    // std::cout << "calc_msb_mask() called" << std::endl;
    assert(pos < bits_per_block);
    if (pos == 0) return 0;
    Block mask = 1;
    for (std::size_t i=1 ; i<pos ; i++) {
      mask <<= 1u;
      mask |= 1u;
    }
    return mask;
  }

  // init from string?

  static std::size_t calc_nu_blocks(std::size_t nu_bits) {
    return nu_bits / bits_per_block
        + static_cast<std::size_t>( nu_bits % bits_per_block !=0 );
  }

  static const std::size_t & msb_pos();

  static const Block & msb_mask();

  static std::size_t block_index(std::size_t pos) {
    //return pos >> block_index_shift;
    return pos >> msb_pos();
  }

  static std::size_t bit_index(std::size_t pos) {
    //return pos & block_index_mask;
    return pos & msb_mask();
  }

  static Block bit_mask(std::size_t pos) {
    return Block(1) << bit_index(pos);
  }

  static std::size_t count_extra_bits(std::size_t nu_bits) {
    return bit_index(nu_bits);
  }

  static Block extra_bits_mask(std::size_t nu_bits) {
    std::size_t n = count_extra_bits(nu_bits);
    Block mask;
    if (n==0) mask = (~static_cast<Block>(0));
    else mask = (static_cast<Block>(1u) << count_extra_bits(nu_bits)) - 1;
    return mask;
  }

  static std::size_t pop_count(Block num);

  static std::size_t ctz(Block num); // count trailing zero

  // need extra_bits_mask
  // for all, set(), reset(), (for none() also?)

};

template<typename Block> inline
const std::size_t & VariableBitsetTraits<Block>::msb_pos() {
  static const std::size_t msb_pos_ret = calc_msb_position(bits_per_block);
  return msb_pos_ret;
}

template<typename Block> inline
const Block & VariableBitsetTraits<Block>::msb_mask() {
  static const Block msb_mask_ret = calc_msb_mask(msb_pos());
  return msb_mask_ret;
}

template<> inline
std::size_t VariableBitsetTraits<unsigned int>::pop_count(unsigned int num) {
  return __builtin_popcountl(num);
}

template<> inline
std::size_t VariableBitsetTraits<uint64>::pop_count(uint64 num) {
  return __builtin_popcountll(num);
}

template<> inline
std::size_t VariableBitsetTraits<unsigned int>::ctz(unsigned int num ) {
  return __builtin_ctzl(num);
}

template<> inline
std::size_t VariableBitsetTraits<uint64>::ctz(uint64 num) {
  return __builtin_ctzll(num);
}

//==============================================================================

template<typename Block> class VariableBitset;
template<typename Block> class VariableBitsetArray;

//==============================================================================

template<typename Block>
class VariableBitsetHelper {
 public:
  typedef VariableBitsetTraits<Block> traits;

  VariableBitsetHelper(std::size_t n)
      : nu_bits_(n),
        nu_blocks_(traits::calc_nu_blocks(nu_bits_))
  {}

  // allocate single bitset
  Block * New() const;
  // allocate array of bitset
  Block * NewArray(std::size_t size) const;
  int NewArraySize(std::size_t size) const { return NuBlocks() * size; }

  void Delete(Block * a) const;

  Block * N(Block * a, std::size_t i) const;
  const Block * N(const Block * a, std::size_t i) const;

  void Copy(const VariableBitset<Block> & src, Block * elm) const;
  void Copy(const Block * src, Block * elm) const;

  void And(const Block * src, Block * elm) const;
  void Or(const Block * src, Block * elm) const;
  void Xor(const Block * src, Block * elm) const;
  void AndNot(const Block * src, Block * elm) const;
  void Flip(Block * elm) const;

  // update and return count (popcnt)
  std::size_t AndCountUpdate(const Block * src, Block * elm) const;
  std::size_t OrCountUpdate(const Block * src, Block * elm) const;
  std::size_t XorCountUpdate(const Block * src, Block * elm) const;
  std::size_t AndNotCountUpdate(const Block * src, Block * elm) const;
  std::size_t FlipCountUpdate(Block * elm) const;

  // only does count (popcnt)
  std::size_t AndCount(const Block * src, const Block * elm) const;
  std::size_t OrCount(const Block * src, const Block * elm) const;
  std::size_t XorCount(const Block * src, const Block * elm) const;
  std::size_t AndNotCount(const Block * src, const Block * elm) const;
  std::size_t FlipCount(const Block * elm) const;

  // for core_i calculation with shortcut
  // return true if (op & elm) == comp
  bool AndEqualTo(const Block * comp, const Block * op, const Block * elm) const;
  // for core_i calculation with shortcut
  // return true if (op & elm) == comp
  // and also overwrite elm with (op & elm)
  bool AndEqualToUpdate3(const Block * comp, const Block * op, Block * elm) const;
  // for core_i calculation with shortcut
  // return true if (op & elm) == comp
  // and also overwrite comp with (op & elm)
  bool AndEqualToUpdate1(Block * comp, const Block * op, const Block * elm) const;

  bool IsSubsetOf(const Block * src, const Block * elm) const;
  bool IsSupersetOf(const Block * src, const Block * elm) const;
  bool IsEqualTo(const Block * src, const Block * elm) const;

  Block * Set(std::size_t pos, bool val, Block * elm) const;
  Block * Doset(std::size_t pos, Block * elm) const;
  Block * Reset(std::size_t pos, Block * elm) const;

  void Set(Block * elm) const; // set elm all 1
  void Reset(Block * elm) const; // set elm all 0

  // return nu_bits if not found
  std::size_t FindFirst(const Block * src) const;
  std::size_t FindNext(const Block * src, std::size_t prev_pos) const;

  bool Test(const Block * elm, std::size_t pos) const;

  std::size_t Count(const Block * elm) const;

  std::ostream& Print(std::ostream & out, const Block * elm) const;

  void InitFromString(const std::string & str, Block * elm) const;

  std::size_t NuBits() const { return nu_bits_; }
  std::size_t NuBlocks() const { return nu_blocks_; }

  std::size_t nu_bits_;
  std::size_t nu_blocks_;

 private:
};

template<typename Block>
void VariableBitsetHelper<Block>::InitFromString(const std::string & str, Block * elm) const {
  std::size_t pos = nu_bits_ - 1;
  
  for (std::size_t i=0 ; i<str.length() ; i++) {
    if (str[i] == '1') Doset(pos, elm);
    else Reset(pos, elm);
    pos--;
  }
}

template<typename Block>
std::ostream& VariableBitsetHelper<Block>::Print(std::ostream & out,
                                                 const Block * elm) const {
  std::stringstream s;

  for (std::size_t i=nu_bits_; i > 0; --i) {
    if (Test(elm, i-1)) s << "1";
    else s << "0";
  }
  
  out << s.str() << std::flush;
  return out;
}

template<typename Block> inline
Block * VariableBitsetHelper<Block>::New() const {
  if (nu_bits_ == 0) return NULL;
  Block * elm = new Block[ NuBlocks() ];
  for (std::size_t i=0;i<NuBlocks();i++) elm[i] = static_cast<Block>(0);
  return elm;
}

template<typename Block> inline
Block * VariableBitsetHelper<Block>::NewArray(std::size_t size) const {
  if (nu_bits_ == 0) return NULL;
  if (size == 0) return NULL;
  Block * array = new Block[NuBlocks() * size];
  for (std::size_t i=0; i<NuBlocks()*size ;i++) array[i] = static_cast<Block>(0);
  return array;
}

template<typename Block> inline
void VariableBitsetHelper<Block>::Delete(Block * a) const {
  delete[] a;
}

template<typename Block> inline
Block * VariableBitsetHelper<Block>::N(Block * array, std::size_t i) const {
  return &(array[i * NuBlocks()]);
}

template<typename Block> inline
const Block * VariableBitsetHelper<Block>::N(const Block * array, std::size_t i) const {
  return &(array[i * NuBlocks()]);
}

template<typename Block> inline
Block * VariableBitsetHelper<Block>::Set(std::size_t pos, bool val, Block * elm) const {
  assert(pos < nu_bits_);
  if (val)
    Doset(pos, elm);
  else
    Reset(pos, elm);
  return elm;
}

template<typename Block> inline
Block * VariableBitsetHelper<Block>::Doset(std::size_t pos, Block * elm) const {
  elm[traits::block_index(pos)] |= traits::bit_mask(pos);
  return elm;
}

template<typename Block> inline
Block * VariableBitsetHelper<Block>::Reset(std::size_t pos, Block * elm) const {
  elm[traits::block_index(pos)] &= ~(traits::bit_mask(pos));
  return elm;
}

template<typename Block> inline
void VariableBitsetHelper<Block>::Set(Block * elm) const {
  std::size_t i = NuBlocks();
  {
    if (i==0) return;
    i--;
    elm[i] = ~(static_cast<Block>(0));
    elm[i] &= traits::extra_bits_mask(nu_bits_); // do sanitize
  }
  while(true) {
    if (i==0) return;
    i--;
    elm[i] = ~(static_cast<Block>(0));
  }
}

template<typename Block> inline
void VariableBitsetHelper<Block>::Reset(Block * elm) const {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] = static_cast<Block>(0);
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::FindFirst(const Block * src) const {
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp = src[i];
    if ( tmp != static_cast<Block>(0) )
      return (i * traits::bits_per_block + traits::ctz(tmp));
  }
  return nu_bits_; // not found
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::FindNext(const Block * src,
                                                  std::size_t prev_pos) const {
  prev_pos++; // increment

  if (prev_pos >= nu_bits_) 
	return nu_bits_; // not found

  std::size_t i = traits::block_index(prev_pos);
  { // first word
    Block tmp = src[i];
    tmp &= (~static_cast<Block>(0)) << traits::bit_index(prev_pos);
    if (tmp != static_cast<Block>(0))
      return (i * traits::bits_per_block + traits::ctz(tmp));
  }
  i++;

  // rest words
  for ( ;i<NuBlocks();i++) {
    Block tmp = src[i];
    if ( tmp != static_cast<Block>(0) )
      return (i * traits::bits_per_block + traits::ctz(tmp));
  }

  return nu_bits_; // not found
}

template<typename Block> inline
bool VariableBitsetHelper<Block>::Test(const Block * elm, std::size_t pos) const {
  assert(pos < nu_bits_);
  return (elm[traits::block_index(pos)] & traits::bit_mask(pos)) != 0;
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::Count(const Block * src) const {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++)
    count += traits::pop_count(src[i]);
  return count;
}

template<typename Block> inline
void VariableBitsetHelper<Block>::Copy(const VariableBitset<Block> & src, Block * elm) const {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] = src.elm_[i];
}

template<typename Block> inline
void VariableBitsetHelper<Block>::Copy(const Block * src, Block * elm) const {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] = src[i];
}

template<typename Block> inline
void VariableBitsetHelper<Block>::And(const Block * src, Block * elm) const {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] &= src[i];
}

template<typename Block> inline
void VariableBitsetHelper<Block>::Or(const Block * src, Block * elm) const {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] |= src[i];
}

template<typename Block> inline
void VariableBitsetHelper<Block>::Xor(const Block * src, Block * elm) const {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] ^= src[i];
}

template<typename Block> inline
void VariableBitsetHelper<Block>::AndNot(const Block * src, Block * elm) const {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] &= ~(src[i]);
}

template<typename Block> inline
void VariableBitsetHelper<Block>::Flip(Block * elm) const {
  std::size_t i = NuBlocks();
  {
    if (i==0) return;
    i--;
    elm[i] = ~(elm[i]);
    elm[i] &= traits::extra_bits_mask(nu_bits_); // do sanitize
  }
  while(true) {
    if (i==0) return;
    i--;
    elm[i] = ~(elm[i]);
  }
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::AndCountUpdate(const Block * src, Block * elm) const {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] &= src[i];
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::OrCountUpdate(const Block * src, Block * elm) const {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] |= src[i];
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::XorCountUpdate(const Block * src, Block * elm) const {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] ^= src[i];
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::AndNotCountUpdate(const Block * src, Block * elm) const {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] &= ~(src[i]);
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::FlipCountUpdate(Block * elm) const {
  std::size_t count=0;
  std::size_t i = NuBlocks();
  {
    if (i==0) return count;
    i--;
    elm[i] = ~(elm[i]);
    elm[i] &= traits::extra_bits_mask(nu_bits_); // do sanitize
    count += traits::pop_count(elm[i]);
  }
  while(true) {
    if (i==0) return count;
    i--;
    elm[i] = ~(elm[i]);
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::AndCount(const Block * src, const Block * elm) const {
  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm[i] & src[i];
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::OrCount(const Block * src, const Block * elm) const {
  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm[i] | src[i];
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::XorCount(const Block * src, const Block * elm) const {
  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm[i] ^ src[i];
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::AndNotCount(const Block * src, const Block * elm) const {
  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm[i] & ~(src[i]);
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetHelper<Block>::FlipCount(const Block * elm) const {
  std::size_t count=0;
  std::size_t i = NuBlocks();
  Block tmp;
  {
    if (i==0) return count;
    i--;
    tmp = ~(elm[i]);
    tmp &= traits::extra_bits_mask(nu_bits_); // do sanitize
    count += traits::pop_count(tmp);
  }
  while(true) {
    if (i==0) return count;
    i--;
    tmp = ~(elm[i]);
    count += traits::pop_count(tmp);
  }
  return count;
}

// for core_i calculation with shortcut
// return true if (src1 & src2) == elm
template<typename Block>
bool VariableBitsetHelper<Block>::AndEqualTo(const Block * comp, const Block * op,
                                             const Block * elm) const {
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp;
    tmp = op[i] & elm[i];
    if (tmp != comp[i]) return false;
  }
  return true;
}

// for core_i calculation with shortcut
// return true if (op & elm) == comp
// and also overwrite elm with (op & elm)
template<typename Block>
bool VariableBitsetHelper<Block>::AndEqualToUpdate3(const Block * comp,
                                                    const Block * op,
                                                    Block * elm) const {
  bool res = true;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] &= op[i];
    // note: can skip rest of this comparison if res becomes false
    if (elm[i] != comp[i]) res = false;
  }
  return res;
}

// for core_i calculation with shortcut
// return true if (op & elm) == comp
// and also overwrite comp with (op & elm)
template<typename Block>
bool VariableBitsetHelper<Block>::AndEqualToUpdate1(Block * comp,
                                                    const Block * op,
                                                    const Block * elm) const {
  bool res = true;
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp;
    tmp = elm[i] & op[i];
    // note: can skip rest of this comparison if res becomes false
    if (tmp != comp[i]) res = false;
    comp[i] = tmp;
  }
  return res;
}

// src is subset of elm
template<typename Block>
bool VariableBitsetHelper<Block>::IsSubsetOf(const Block * src, const Block * elm) const {
  // example
  // this = 11000, src = 11111
  // 11000 & 11111 = 11000
  // 11000 ^ 11000 = 00000
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp;
    tmp = elm[i] & src[i];
    tmp ^= src[i];
    if (tmp != 0) return false;
  }
  return true;
}

// src is superset of elm
template<typename Block>
bool VariableBitsetHelper<Block>::IsSupersetOf(const Block * src, const Block * elm) const {
  // example
  // this = 11111, src = 11000
  // 11111 & 11000 = 11000
  // 11000 ^ 11000 = 00000
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp;
    tmp = elm[i] & src[i];
    tmp ^= elm[i];
    if (tmp != 0) return false;
  }
  return true;
}

template<typename Block>
bool VariableBitsetHelper<Block>::IsEqualTo(const Block * src, const Block * elm) const {
  for (std::size_t i=0;i<NuBlocks();i++) {
    if (elm[i] != src[i]) return false;
  }
  return true;
}

//==============================================================================

template<typename Block>
class VariableBitset {
 public:
  typedef VariableBitsetTraits<Block> traits;
  typedef std::vector<Block> buffer_type;

  // note this is not Block
  template <typename B>
  friend std::ostream& operator<<(std::ostream & out, const VariableBitset<B> & vbs);

  friend class VariableBitsetArray<Block>;

  VariableBitset() {
    nu_bits = 0;
    elm_ = NULL;
  }

  /** nu_bits: number of bits in one bitset */
  VariableBitset(int num) {
    nu_bits = num;
    if (nu_bits == 0) return;
    elm_ = new Block[NuBlocks()];
    for (std::size_t i=0;i<NuBlocks();i++) elm_[i] = 0;
  }

  VariableBitset(const std::string & str) {
    nu_bits = str.length();
    if (nu_bits == 0) return;
    elm_ = new Block[NuBlocks()];
    for (std::size_t i=0;i<NuBlocks();i++) elm_[i] = 0;
    InitFromString(str);
  }

  VariableBitset(const VariableBitset<Block> & src) {
    nu_bits = src.nu_bits;
    if (nu_bits == 0) return;
    elm_ = new Block[NuBlocks()];
    for (std::size_t i=0;i<NuBlocks();i++) elm_[i] = src.elm_[i];
  }

  ~VariableBitset() {
    if (nu_bits > 0) delete elm_;
  }
  
  // takes parameter by value (e.g. copy) and then swap
  VariableBitset<Block> & operator=(VariableBitset<Block> src) {
    std::swap(nu_bits, src.nu_bits);
    std::swap(elm_, src.elm_);
    return *this;
  }

  void CopyFrom(Block * src) {
    for (std::size_t i=0;i<NuBlocks();i++)
      elm_[i] = src[i];
  }

  void CopyFrom(const Block * src) {
    for (std::size_t i=0;i<NuBlocks();i++)
      elm_[i] = src[i];
  }

  Block * Ptr() { return elm_; }

  void And(const VariableBitset<Block> & src);
  void Or(const VariableBitset<Block> & src);
  void Xor(const VariableBitset<Block> & src);
  void AndNot(const VariableBitset<Block> & src);
  void Flip();

  // update and return count (popcnt)
  std::size_t AndCountUpdate(const VariableBitset<Block> & src);
  std::size_t OrCountUpdate(const VariableBitset<Block> & src);
  std::size_t XorCountUpdate(const VariableBitset<Block> & src);
  std::size_t AndNotCountUpdate(const VariableBitset<Block> & src);
  std::size_t FlipCountUpdate();

  // only does count (popcnt)
  std::size_t AndCount(const VariableBitset<Block> & src) const;
  std::size_t OrCount(const VariableBitset<Block> & src) const;
  std::size_t XorCount(const VariableBitset<Block> & src) const;
  std::size_t AndNotCount(const VariableBitset<Block> & src) const;
  std::size_t FlipCount() const;

  bool IsSubsetOf(const VariableBitset<Block> & src) const;
  bool IsSupersetOf(const VariableBitset<Block> & src) const;
  bool IsEqualTo(const VariableBitset<Block> & src) const;
  bool operator==(const VariableBitset<Block> & src) const;

  VariableBitset<Block> & Set(std::size_t pos, bool val);
  VariableBitset<Block> & Doset(std::size_t pos);
  VariableBitset<Block> & Reset(std::size_t pos);

  // todo: test
  void Set(); // set all 1
  void Reset(); // set all 0

  std::size_t FindFirst() const; // return nu_bits if not found
  std::size_t FindNext(std::size_t prev_pos) const; // return nu_bits if not found

  std::size_t nu_bits;

  bool Test(std::size_t i) const;
  std::size_t Count() const;

  // todo: prepare bit_mask(pos) as a static function

  // below boost::dynamic_bitset source
  // m_bits[block_index(pos)] |= bit_mask(pos);

  //   private (static) functions
  // block_width_type count_extra_bits() const BOOST_NOEXCEPT { return bit_index(nu_bits); }

 private:
  Block * elm_;

  std::size_t NuBits() const { return nu_bits; }
  std::size_t NuBlocks() const { return traits::calc_nu_blocks(nu_bits); }

  void InitFromString(const std::string & str);

};

template<typename Block>
void VariableBitset<Block>::InitFromString(const std::string & str) {
  std::size_t pos = nu_bits - 1;
  
  for (std::size_t i=0 ; i<str.length() ; i++) {
    if (str[i] == '1') Doset(pos);
    else Reset(pos);
    pos--;
  }
}

template<typename Block> inline
VariableBitset<Block> & VariableBitset<Block>::Set(std::size_t pos, bool val) {
  assert(pos < nu_bits);
  if (val)
    Doset(pos);
  else
    Reset(pos);
  return *this;
}

template<typename Block> inline
VariableBitset<Block> & VariableBitset<Block>::Doset(std::size_t pos) {
  assert(pos < nu_bits);
  elm_[traits::block_index(pos)] |= traits::bit_mask(pos);
  return *this;
}

template<typename Block> inline
VariableBitset<Block> & VariableBitset<Block>::Reset(std::size_t pos) {
  assert(pos < nu_bits);
  elm_[traits::block_index(pos)] &= ~(traits::bit_mask(pos));
  return *this;
}

template<typename Block> inline
std::size_t VariableBitset<Block>::FindFirst() const {
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp = elm_[i];
    if ( tmp != static_cast<Block>(0) )
      return (i * traits::bits_per_block + traits::ctz(tmp));
  }
  return nu_bits; // not found
}

template<typename Block> inline
std::size_t VariableBitset<Block>::FindNext(std::size_t prev_pos) const {
  prev_pos++; // increment

  if (prev_pos >= nu_bits) 
	return nu_bits; // not found

  std::size_t i = traits::block_index(prev_pos);
  { // first word
    Block tmp = elm_[i];
    tmp &= (~static_cast<Block>(0)) << traits::bit_index(prev_pos);
    if (tmp != static_cast<Block>(0))
      return (i * traits::bits_per_block + traits::ctz(tmp));
  }
  i++;

  // rest words
  for ( ;i<NuBlocks();i++) {
    Block tmp = elm_[i];
    if ( tmp != static_cast<Block>(0) )
      return (i * traits::bits_per_block + traits::ctz(tmp));
  }

  return nu_bits; // not found
}

template<typename Block> inline
bool VariableBitset<Block>::Test(std::size_t pos) const {
  assert(pos < nu_bits);
  return (elm_[traits::block_index(pos)] & traits::bit_mask(pos)) != 0;
}

template<typename Block> inline
void VariableBitset<Block>::And(const VariableBitset<Block> & src) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm_[i] &= src.elm_[i];
}

template<typename Block> inline
void VariableBitset<Block>::Or(const VariableBitset<Block> & src) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm_[i] |= src.elm_[i];
}

template<typename Block> inline
void VariableBitset<Block>::Xor(const VariableBitset<Block> & src) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm_[i] ^= src.elm_[i];
}

template<typename Block> inline
void VariableBitset<Block>::AndNot(const VariableBitset<Block> & src) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm_[i] &= ~(src.elm_[i]);
}

template<typename Block> inline
void VariableBitset<Block>::Flip() {
  std::size_t i = NuBlocks();
  {
    if (i==0) return;
    i--;
    elm_[i] = ~(elm_[i]);
    elm_[i] &= traits::extra_bits_mask(nu_bits); // do sanitize
  }
  while(true) {
    if (i==0) return;
    i--;
    elm_[i] = ~(elm_[i]);
  }
}

template<typename Block> inline
void VariableBitset<Block>::Set() {
  std::size_t i = NuBlocks();
  {
    if (i==0) return;
    i--;
    elm_[i] = ~(static_cast<Block>(0));
    elm_[i] &= traits::extra_bits_mask(nu_bits); // do sanitize
  }
  while(true) {
    if (i==0) return;
    i--;
    elm_[i] = ~(static_cast<Block>(0));
  }
}

template<typename Block> inline
void VariableBitset<Block>::Reset() {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm_[i] = static_cast<Block>(0);
}

template<typename Block> inline
std::size_t VariableBitset<Block>::AndCountUpdate(const VariableBitset<Block> & src) {
  assert(nu_bits == src.nu_bits);

  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm_[i] &= src.elm_[i];
    count += traits::pop_count(elm_[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitset<Block>::OrCountUpdate(const VariableBitset<Block> & src) {
  assert(nu_bits == src.nu_bits);

  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm_[i] |= src.elm_[i];
    count += traits::pop_count(elm_[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitset<Block>::XorCountUpdate(const VariableBitset<Block> & src) {
  assert(nu_bits == src.nu_bits);

  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm_[i] ^= src.elm_[i];
    count += traits::pop_count(elm_[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitset<Block>::AndNotCountUpdate(const VariableBitset<Block> & src) {
  assert(nu_bits == src.nu_bits);

  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm_[i] &= ~(src.elm_[i]);
    count += traits::pop_count(elm_[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitset<Block>::FlipCountUpdate() {
  std::size_t count = 0;
  std::size_t i = NuBlocks();
  {
    if (i==0) return 0;
    i--;
    elm_[i] = ~(elm_[i]);
    elm_[i] &= traits::extra_bits_mask(nu_bits); // do sanitize
    count += traits::pop_count(elm_[i]);
  }

  while(true) {
    if (i==0) break;
    i--;
    elm_[i] = ~(elm_[i]);
    count += traits::pop_count(elm_[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitset<Block>::AndCount(const VariableBitset<Block> & src) const {
  assert(nu_bits == src.nu_bits);

  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm_[i] & src.elm_[i];
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitset<Block>::OrCount(const VariableBitset<Block> & src) const {
  assert(nu_bits == src.nu_bits);

  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm_[i] | src.elm_[i];
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitset<Block>::XorCount(const VariableBitset<Block> & src) const {
  assert(nu_bits == src.nu_bits);

  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm_[i] ^ src.elm_[i];
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitset<Block>::AndNotCount(const VariableBitset<Block> & src) const {
  assert(nu_bits == src.nu_bits);

  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm_[i] & ~(src.elm_[i]);
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitset<Block>::FlipCount() const {
  std::size_t count = 0;
  std::size_t i = NuBlocks();
  Block tmp;
  {
    if (i==0) return 0;
    i--;
    tmp = ~(elm_[i]);
    tmp &= traits::extra_bits_mask(nu_bits); // do sanitize
    count += traits::pop_count(tmp);
  }

  while(true) {
    if (i==0) break;
    i--;
    tmp = ~(elm_[i]);
    count += traits::pop_count(tmp);
  }
  return count;
}

// *this is subset of src
template<typename Block> inline
bool VariableBitset<Block>::IsSubsetOf(const VariableBitset<Block> & src) const {
  assert(nu_bits == src.nu_bits);
  // example
  // this = 11000, src = 11111
  // 11000 & 11111 = 11000
  // 11000 ^ 11000 = 00000
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp;
    tmp = elm_[i] & src.elm_[i];
    tmp ^= elm_[i];
    if (tmp != 0) return false;
  }
  return true;
}

template<typename Block> inline
bool VariableBitset<Block>::IsSupersetOf(const VariableBitset<Block> & src) const {
  assert(nu_bits == src.nu_bits);
  // example
  // this = 11111, src = 11000
  // 11111 & 11000 = 11000
  // 11000 ^ 11000 = 00000
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp;
    tmp = elm_[i] & src.elm_[i];
    tmp ^= src.elm_[i];
    if (tmp != 0) return false;
  }
  return true;
}

template<typename Block> inline
bool VariableBitset<Block>::IsEqualTo(const VariableBitset<Block> & src) const {
  if (nu_bits != src.nu_bits) return false;
  for (std::size_t i=0;i<NuBlocks();i++) {
    if (elm_[i] != src.elm_[i]) return false;
  }
  return true;
}

template<typename Block> inline
bool VariableBitset<Block>::operator==(const VariableBitset<Block> & src) const {
  return IsEqualTo(src);
}

template<typename Block> inline
std::size_t VariableBitset<Block>::Count() const {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    // std::cout << "elm=0x" << std::hex << elm_[i]
    //           << "\5pop=" << std::dec << traits::pop_count(elm_[i]) << std::endl;
    count += traits::pop_count(elm_[i]);
  }
  return count;
}

//==============================================================================

template<typename Block>
class VariableBitsetArray {
  typedef VariableBitsetTraits<Block> traits;
  typedef std::vector<Block> buffer_type;

 public:
  /** parameter
   *  nu_bits: number of bits in one bitset
   *  size:    number of array entries */
  VariableBitsetArray(std::size_t nu_bits, int size);
  // todo: implement
  VariableBitsetArray(std::size_t nu_bits, int size, bool cache_line_align);
  
  ~VariableBitsetArray();

  Block * N(std::size_t i);
  const Block * N(std::size_t i) const;

  void Copy(const VariableBitset<Block> & src, Block * elm);
  void Copy(const Block * src, Block * elm);

  // elm is the entry in array_, 2nd operand is modified
  void And(const VariableBitset<Block> & src, Block * elm);
  void And(const Block * src, Block * elm);
  void And(const Block * src, VariableBitset<Block> * dst);

  void Or(const VariableBitset<Block> & src, Block * elm);
  void Or(const Block * src, Block * elm);
  void Or(const Block * src, VariableBitset<Block> * dst);

  void Xor(const VariableBitset<Block> & src, Block * elm);
  void Xor(const Block * src, Block * elm);
  void Xor(const Block * src, VariableBitset<Block> * dst);

  void AndNot(const VariableBitset<Block> & src, Block * elm);
  void AndNot(const Block * src, Block * elm);
  void AndNot(const Block * src, VariableBitset<Block> * dst);

  void Flip(Block * elm);

  // update and return count (popcnt)
  // elm is the entry in array_, 2nd operand is modified
  std::size_t AndCountUpdate(const VariableBitset<Block> & src, Block * elm);
  std::size_t AndCountUpdate(const Block * src, Block * elm);
  std::size_t OrCountUpdate(const VariableBitset<Block> & src, Block * elm);
  std::size_t OrCountUpdate(const Block * src, Block * elm);
  std::size_t XorCountUpdate(const VariableBitset<Block> & src, Block * elm);
  std::size_t XorCountUpdate(const Block * src, Block * elm);
  std::size_t AndNotCountUpdate(const VariableBitset<Block> & src, Block * elm);
  std::size_t AndNotCountUpdate(const Block * src, Block * elm);
  std::size_t FlipCountUpdate(Block * elm);

  // only does count (popcnt)
  std::size_t AndCount(const VariableBitset<Block> & src, const Block * elm) const;
  std::size_t AndCount(const Block * src, const Block * elm) const;
  std::size_t OrCount(const VariableBitset<Block> & src, const Block * elm) const;
  std::size_t OrCount(const Block * src, const Block * elm) const;
  std::size_t XorCount(const VariableBitset<Block> & src, const Block * elm) const;
  std::size_t XorCount(const Block * src, const Block * elm) const;
  std::size_t AndNotCount(const VariableBitset<Block> & src, const Block * elm) const;
  std::size_t AndNotCount(const Block * src, const Block * elm) const;
  std::size_t FlipCount(const Block * elm) const;

  bool IsSubsetOf(const VariableBitset<Block> & src, const Block * elm) const;
  bool IsSubsetOf(const Block * src, const Block * elm) const;
  bool IsSupersetOf(const VariableBitset<Block> & src, const Block * elm) const;
  bool IsSupersetOf(const Block * src, const Block * elm) const;
  bool IsEqualTo(const VariableBitset<Block> & src, const Block * elm) const;
  bool IsEqualTo(const Block * src, const Block * elm) const;

  Block * Set(std::size_t pos, bool val, Block * elm);
  Block * Doset(std::size_t pos, Block * elm);
  Block * Reset(std::size_t pos, Block * elm);

  void Set(Block * elm); // set elm all 1
  void Reset(Block * elm); // set elm all 0

  std::size_t FindFirst(const Block * src) const;
  std::size_t FindNext(const Block * src, std::size_t prev_pos) const;

  bool Test(const Block * elm, std::size_t pos) const;
  std::size_t Count(const Block * elm) const;

  std::size_t nu_bits;
  std::size_t size;

  std::ostream& ShowNth(std::ostream & out, const Block * elm) const;

 private:
  Block * array_;

  std::size_t NuBlocks() const { return traits::calc_nu_blocks(nu_bits); }

};

template<typename Block>
VariableBitsetArray<Block>::VariableBitsetArray(std::size_t p_nu_bits, int p_size) {
  nu_bits = p_nu_bits;
  if (nu_bits == 0) return;
  size = p_size;
  array_ = new Block[NuBlocks() * size];
  for (std::size_t i=0; i<NuBlocks()*size ;i++) array_[i] = 0;
}

template<typename Block>
VariableBitsetArray<Block>::~VariableBitsetArray() {
  if (nu_bits > 0) delete array_;
}

template<typename Block> inline
Block * VariableBitsetArray<Block>::Set(std::size_t pos, bool val, Block * elm) {
  assert(pos < nu_bits);
  if (val)
    Doset(pos, elm);
  else
    Reset(pos, elm);
  return elm;
}

template<typename Block> inline
Block * VariableBitsetArray<Block>::Doset(std::size_t pos, Block * elm) {
  elm[traits::block_index(pos)] |= traits::bit_mask(pos);
  return elm;
}

template<typename Block> inline
Block * VariableBitsetArray<Block>::Reset(std::size_t pos, Block * elm) {
  elm[traits::block_index(pos)] &= ~(traits::bit_mask(pos));
  return elm;
}

template<typename Block> inline
void VariableBitsetArray<Block>::Set(Block * elm) {
  std::size_t i = NuBlocks();
  {
    if (i==0) return;
    i--;
    elm[i] = ~(static_cast<Block>(0));
    elm[i] &= traits::extra_bits_mask(nu_bits); // do sanitize
  }
  while(true) {
    if (i==0) return;
    i--;
    elm[i] = ~(static_cast<Block>(0));
  }
}

template<typename Block> inline
void VariableBitsetArray<Block>::Reset(Block * elm) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] = static_cast<Block>(0);
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::FindFirst(const Block * src) const {
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp = src[i];
    if ( tmp != static_cast<Block>(0) )
      return (i * traits::bits_per_block + traits::ctz(tmp));
  }
  return nu_bits; // not found
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::FindNext(const Block * src,
                                                 std::size_t prev_pos) const {
  prev_pos++; // increment

  if (prev_pos >= nu_bits) 
	return nu_bits; // not found

  std::size_t i = traits::block_index(prev_pos);
  { // first word
    Block tmp = src[i];
    tmp &= (~static_cast<Block>(0)) << traits::bit_index(prev_pos);
    if (tmp != static_cast<Block>(0))
      return (i * traits::bits_per_block + traits::ctz(tmp));
  }
  i++;

  // rest words
  for ( ;i<NuBlocks();i++) {
    Block tmp = src[i];
    if ( tmp != static_cast<Block>(0) )
      return (i * traits::bits_per_block + traits::ctz(tmp));
  }

  return nu_bits; // not found
}

template<typename Block> inline
bool VariableBitsetArray<Block>::Test(const Block * elm, std::size_t pos) const {
  assert(pos < nu_bits);
  return (elm[traits::block_index(pos)] & traits::bit_mask(pos)) != 0;
}

template<typename Block> inline
Block * VariableBitsetArray<Block>::N(std::size_t i) {
  return &(array_[i * NuBlocks()]);
}

template<typename Block> inline
const Block * VariableBitsetArray<Block>::N(std::size_t i) const {
  return &(array_[i * NuBlocks()]);
}

template<typename Block> inline
void VariableBitsetArray<Block>::Copy(const VariableBitset<Block> & src, Block * elm) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] = src.elm_[i];
}

template<typename Block> inline
void VariableBitsetArray<Block>::Copy(const Block * src, Block * elm) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] = src[i];
}

template<typename Block> inline
void VariableBitsetArray<Block>::And(const Block * src, Block * elm) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] &= src[i];
}

template<typename Block> inline
void VariableBitsetArray<Block>::And(const VariableBitset<Block> & src, Block * elm) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] &= src.elm_[i];
}

template<typename Block> inline
void VariableBitsetArray<Block>::And(const Block * src, VariableBitset<Block> * dst) {
  for (std::size_t i=0;i<NuBlocks();i++)
    (dst->elm_)[i] &= src[i];
}

template<typename Block> inline
void VariableBitsetArray<Block>::Or(const Block * src, Block * elm) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] |= src[i];
}

template<typename Block> inline
void VariableBitsetArray<Block>::Or(const VariableBitset<Block> & src, Block * elm) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] |= src.elm_[i];
}

template<typename Block> inline
void VariableBitsetArray<Block>::Or(const Block * src, VariableBitset<Block> * dst) {
  for (std::size_t i=0;i<NuBlocks();i++)
    (dst->elm_)[i] |= src[i];
}

template<typename Block> inline
void VariableBitsetArray<Block>::Xor(const Block * src, Block * elm) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] ^= src[i];
}

template<typename Block> inline
void VariableBitsetArray<Block>::Xor(const VariableBitset<Block> & src, Block * elm) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] ^= src.elm_[i];
}

template<typename Block> inline
void VariableBitsetArray<Block>::Xor(const Block * src, VariableBitset<Block> * dst) {
  for (std::size_t i=0;i<NuBlocks();i++)
    (dst->elm_)[i] ^= src[i];
}

template<typename Block> inline
void VariableBitsetArray<Block>::AndNot(const Block * src, Block * elm) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] &= ~(src[i]);
}

template<typename Block> inline
void VariableBitsetArray<Block>::AndNot(const VariableBitset<Block> & src, Block * elm) {
  for (std::size_t i=0;i<NuBlocks();i++)
    elm[i] &= ~(src.elm_[i]);
}

template<typename Block> inline
void VariableBitsetArray<Block>::AndNot(const Block * src, VariableBitset<Block> * dst) {
  for (std::size_t i=0;i<NuBlocks();i++)
    (dst->elm_)[i] &= ~(src[i]);
}

template<typename Block> inline
void VariableBitsetArray<Block>::Flip(Block * elm) {
  std::size_t i = NuBlocks();
  {
    if (i==0) return;
    i--;
    elm[i] = ~(elm[i]);
    elm[i] &= traits::extra_bits_mask(nu_bits); // do sanitize
  }
  while(true) {
    if (i==0) return;
    i--;
    elm[i] = ~(elm[i]);
  }
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::Count(const Block * src) const {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++)
    count += traits::pop_count(src[i]);
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::AndCountUpdate(const Block * src, Block * elm) {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] &= src[i];
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::AndCountUpdate(const VariableBitset<Block> & src, Block * elm) {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] &= src.elm_[i];
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::OrCountUpdate(const Block * src, Block * elm) {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] |= src[i];
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::OrCountUpdate(const VariableBitset<Block> & src, Block * elm) {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] |= src.elm_[i];
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::XorCountUpdate(const Block * src, Block * elm) {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] ^= src[i];
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::XorCountUpdate(const VariableBitset<Block> & src, Block * elm) {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] ^= src.elm_[i];
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::AndNotCountUpdate(const Block * src, Block * elm) {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] &= ~(src[i]);
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::AndNotCountUpdate(const VariableBitset<Block> & src, Block * elm) {
  std::size_t count = 0;
  for (std::size_t i=0;i<NuBlocks();i++) {
    elm[i] &= ~(src.elm_[i]);
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::FlipCountUpdate(Block * elm) {
  std::size_t count=0;
  std::size_t i = NuBlocks();
  {
    if (i==0) return count;
    i--;
    elm[i] = ~(elm[i]);
    elm[i] &= traits::extra_bits_mask(nu_bits); // do sanitize
    count += traits::pop_count(elm[i]);
  }
  while(true) {
    if (i==0) return count;
    i--;
    elm[i] = ~(elm[i]);
    count += traits::pop_count(elm[i]);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::AndCount(const Block * src, const Block * elm) const {
  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm[i] & src[i];
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::AndCount(const VariableBitset<Block> & src, const Block * elm) const {
  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm[i] & src.elm_[i];
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::OrCount(const Block * src, const Block * elm) const {
  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm[i] | src[i];
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::OrCount(const VariableBitset<Block> & src, const Block * elm) const {
  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm[i] | src.elm_[i];
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::XorCount(const Block * src, const Block * elm) const {
  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm[i] ^ src[i];
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::XorCount(const VariableBitset<Block> & src, const Block * elm) const {
  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm[i] ^ src.elm_[i];
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::AndNotCount(const Block * src, const Block * elm) const {
  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm[i] & ~(src[i]);
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::AndNotCount(const VariableBitset<Block> & src, const Block * elm) const {
  std::size_t count = 0;
  Block tmp;
  for (std::size_t i=0;i<NuBlocks();i++) {
    tmp = elm[i] & ~(src.elm_[i]);
    count += traits::pop_count(tmp);
  }
  return count;
}

template<typename Block> inline
std::size_t VariableBitsetArray<Block>::FlipCount(const Block * elm) const {
  std::size_t count=0;
  std::size_t i = NuBlocks();
  Block tmp;
  {
    if (i==0) return count;
    i--;
    tmp = ~(elm[i]);
    tmp &= traits::extra_bits_mask(nu_bits); // do sanitize
    count += traits::pop_count(tmp);
  }
  while(true) {
    if (i==0) return count;
    i--;
    tmp = ~(elm[i]);
    count += traits::pop_count(tmp);
  }
  return count;
}

// src is subset of elm
template<typename Block>
bool VariableBitsetArray<Block>::IsSubsetOf(const Block * src, const Block * elm) const {
  // example
  // this = 11000, src = 11111
  // 11000 & 11111 = 11000
  // 11000 ^ 11000 = 00000
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp;
    tmp = elm[i] & src[i];
    tmp ^= src[i];
    if (tmp != 0) return false;
  }
  return true;
}

// src is subset of elm
template<typename Block>
bool VariableBitsetArray<Block>::IsSubsetOf(const VariableBitset<Block> & src, const Block * elm) const {
  assert(nu_bits == src.nu_bits);
  // example
  // this = 11000, src = 11111
  // 11000 & 11111 = 11000
  // 11000 ^ 11000 = 00000
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp;
    tmp = elm[i] & src.elm_[i];
    tmp ^= src.elm_[i];
    if (tmp != 0) return false;
  }
  return true;
}

// src is superset of elm
template<typename Block>
bool VariableBitsetArray<Block>::IsSupersetOf(const Block * src, const Block * elm) const {
  // example
  // this = 11111, src = 11000
  // 11111 & 11000 = 11000
  // 11000 ^ 11000 = 00000
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp;
    tmp = elm[i] & src[i];
    tmp ^= elm[i];
    if (tmp != 0) return false;
  }
  return true;
}

// src is superset of elm
template<typename Block>
bool VariableBitsetArray<Block>::IsSupersetOf(const VariableBitset<Block> & src, const Block * elm) const {
  assert(nu_bits == src.nu_bits);
  // example
  // this = 11111, src = 11000
  // 11111 & 11000 = 11000
  // 11000 ^ 11000 = 00000
  for (std::size_t i=0;i<NuBlocks();i++) {
    Block tmp;
    tmp = elm[i] & src.elm_[i];
    tmp ^= elm[i];
    if (tmp != 0) return false;
  }
  return true;
}

template<typename Block>
bool VariableBitsetArray<Block>::IsEqualTo(const Block * src, const Block * elm) const {
  for (std::size_t i=0;i<NuBlocks();i++) {
    if (elm[i] != src[i]) return false;
  }
  return true;
}

template<typename Block>
bool VariableBitsetArray<Block>::IsEqualTo(const VariableBitset<Block> & src, const Block * elm) const {
  assert(nu_bits == src.nu_bits);
  for (std::size_t i=0;i<NuBlocks();i++) {
    if (elm[i] != src.elm_[i]) return false;
  }
  return true;
}

//==============================================================================

template <typename Block>
std::ostream& operator<<(std::ostream & out, const VariableBitset<Block> & vbs) {
  std::stringstream s;

  for (std::size_t i=vbs.nu_bits; i > 0; --i) {
    if (vbs.Test(i-1)) s << "1";
    else s << "0";
  }

  out << s.str() << std::flush;
  return out;
}

template<typename Block>
std::ostream& VariableBitsetArray<Block>::ShowNth(std::ostream & out,
                                                  const Block * elm) const {
  std::stringstream s;

  for (std::size_t i=nu_bits; i > 0; --i) {
    if (Test(elm, i-1)) s << "1";
    else s << "0";
  }
  
  out << s.str() << std::flush;
  return out;
}

} // namespace lamp_search

#endif // _LAMP_SEARCH_VARIABLE_BITSET_ARRAY_H_
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
