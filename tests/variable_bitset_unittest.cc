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

#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/foreach.hpp>

#include "gtest/gtest.h"

#include "utils.h"
#include "variable_bitset_array.h"

using namespace lamp_search;

TEST (VariableBitsetTest, Uint64TraitsTest) {
  typedef VariableBitsetTraits<uint64> traits;
  EXPECT_EQ(64u, traits::bits_per_block);
  EXPECT_EQ( 8u, traits::bytes_per_block);

  EXPECT_EQ( 0u, traits::calc_msb_position(  0x1));
  EXPECT_EQ( 8u, traits::calc_msb_position(0x100));
  EXPECT_EQ(63u, traits::calc_msb_position(0xf000000000000000));

  EXPECT_EQ(               0x0u, traits::calc_msb_mask( 0u));
  EXPECT_EQ(               0x1u, traits::calc_msb_mask( 1u));
  EXPECT_EQ(               0x3u, traits::calc_msb_mask( 2u));
  EXPECT_EQ(0x7fffffffffffffffu, traits::calc_msb_mask(63u));

  std::cout << "explicit calls to calc_ functions end here" << std::endl;

  EXPECT_EQ(6u, traits::msb_pos());
  EXPECT_EQ(0x3Fu, traits::msb_mask());

  std::cout << "1st explicit calls to msb_* functions end" << std::endl;

  EXPECT_EQ(6u, traits::msb_pos());
  EXPECT_EQ(0x3Fu, traits::msb_mask());

  std::cout << "2nd explicit calls to msb_* functions end" << std::endl;

  EXPECT_EQ(1u, traits::calc_nu_blocks( 1u));
  EXPECT_EQ(1u, traits::calc_nu_blocks(64u));
  EXPECT_EQ(2u, traits::calc_nu_blocks(65u));
  EXPECT_EQ(2u, traits::calc_nu_blocks(128u));
  EXPECT_EQ(3u, traits::calc_nu_blocks(129u));

  // EXPECT_EQ(6u, traits::block_index_shift);
  // EXPECT_EQ(0x3Fu, traits::block_index_mask);

  EXPECT_EQ(0u, traits::bit_index(0));
  EXPECT_EQ(1u, traits::bit_index(1));
  EXPECT_EQ(7u, traits::bit_index(0x13154bb47));
  EXPECT_EQ(15u, traits::bit_index(0x13154008F));
  EXPECT_EQ(63u, traits::bit_index(0x15311abFF));

  EXPECT_EQ(0u, traits::count_extra_bits(64));
  //EXPECT_EQ(64u, traits::count_extra_bits(64));

  EXPECT_EQ(traits::extra_bits_mask(64), traits::extra_bits_mask(128));
  EXPECT_EQ(0xffffffffffffffffu, traits::extra_bits_mask(64));

  EXPECT_EQ(0x7fffffffffffffffu, traits::extra_bits_mask(127));
}

TEST (VariableBitsetTest, InitTest) {
  VariableBitset<uint64> bs1(70);
  bs1.Set(0, true);
  bs1.Doset(1);
  bs1.Set(2, false);
  bs1.Doset(63);
  bs1.Set(64, true);
  bs1.Doset(69);
  
  VariableBitset<uint64> bs2(
      "100001"
      "10000000" "00000000"
      "00000000" "00000000"
      "00000000" "00000000"
      "00000000" "00000011"
                             );

  std::stringstream ss1, ss2;
  ss1 << bs1;
  ss2 << bs2;

  EXPECT_EQ(ss1.str(), ss2.str());
}

TEST (VariableBitsetTest, SetResetTest) {
  VariableBitset<uint64> bs1(128);
  bs1.Set();
  std::cout << bs1 << std::endl;
  EXPECT_EQ(128u, bs1.Count());
  bs1.Reset();
  std::cout << bs1 << std::endl;
  EXPECT_EQ(  0u, bs1.Count());

}

TEST (VariableBitsetTest, OperationTest) {
  VariableBitset<uint64> bs1(
      "100001" // 2
      "10000000" "00000000" // 1
      "00001100" "00000110" // 4
      "00110000" "00000001" // 3
      "00000000" "00000011" // 2
                             );

  VariableBitset<uint64> bs2(
      "100001"
      "01111100" "00000000"
      "11111111" "00001110"
      "00110000" "00000011"
      "00000011" "00000000"
                             );

  VariableBitset<uint64> bs3(
      "111001"
      "10011000" "11111111"
      "00001100" "00000110"
      "00110011" "01111001"
      "11111100" "00000011"
                             );

  std::cout << bs1 << std::endl;

  EXPECT_TRUE(bs1.Test(0));
  EXPECT_TRUE(bs1.Test(1));
  EXPECT_TRUE(bs1.Test(64));
  EXPECT_TRUE(bs1.Test(69));

  EXPECT_FALSE(bs1.Test(8));
  EXPECT_FALSE(bs1.Test(68));

  EXPECT_EQ((std::size_t)12, bs1.Count());
  
  VariableBitset<uint64> bs1c(bs1);
  EXPECT_TRUE(bs1.IsEqualTo(bs1c));
  EXPECT_EQ(bs1, bs1c);
  std::cout << bs1 << std::endl;
  std::cout << bs1c << std::endl;

  VariableBitset<uint64> bs2c;
  bs2c = bs2;
  EXPECT_EQ(bs2, bs2c);
  std::cout << bs2 << std::endl;
  std::cout << bs2c << std::endl;

  VariableBitset<uint64> tmp;
  std::size_t c;

  // bs1 and bs2
  VariableBitset<uint64> tmp_and(
      "100001"
      "00000000" "00000000"
      "00001100" "00000110"
      "00110000" "00000001"
      "00000000" "00000000"
                             );
  tmp = bs1;
  tmp.And(bs2);
  EXPECT_EQ(tmp_and, tmp);
  std::cout << tmp << std::endl;

  tmp = bs1;
  c = tmp.AndCount(bs2);
  EXPECT_EQ(bs1, tmp);
  EXPECT_EQ((std::size_t)9, c);
  c = tmp.AndCountUpdate(bs2);
  EXPECT_EQ(tmp_and, tmp);
  EXPECT_EQ((std::size_t)9, c);

  // bs1 or bs2
  VariableBitset<uint64> tmp_or(
      "100001"              // 2
      "11111100" "00000000" // 6
      "11111111" "00001110" // 11
      "00110000" "00000011" // 4
      "00000011" "00000011" // 4
                             );

  tmp = bs1;
  tmp.Or(bs2);
  EXPECT_EQ(tmp_or, tmp);
  std::cout << tmp << std::endl;

  tmp = bs1;
  c = tmp.OrCount(bs2);
  EXPECT_EQ(bs1, tmp);
  EXPECT_EQ((std::size_t)27, c);
  c = tmp.OrCountUpdate(bs2);
  EXPECT_EQ(tmp_or, tmp);
  EXPECT_EQ((std::size_t)27, c);

  // bs1 xor bs2
  VariableBitset<uint64> tmp_xor(
      "000000"              // 0
      "11111100" "00000000" // 6
      "11110011" "00001000" // 7
      "00000000" "00000010" // 1
      "00000011" "00000011" // 4
                             );

  tmp = bs1;
  tmp.Xor(bs2);
  EXPECT_EQ(tmp_xor, tmp);
  std::cout << tmp << std::endl;

  tmp = bs1;
  c = tmp.XorCount(bs2);
  EXPECT_EQ(bs1, tmp);
  EXPECT_EQ((std::size_t)18, c);
  c = tmp.XorCountUpdate(bs2);
  EXPECT_EQ(tmp_xor, tmp);
  EXPECT_EQ((std::size_t)18, c);

  // bs1 and_not bs2
  VariableBitset<uint64> tmp_andnot(
      "000000"              // 0
      "10000000" "00000000" // 1
      "00000000" "00000000" // 0
      "00000000" "00000000" // 0
      "00000000" "00000011" // 2
                             );

  tmp = bs1;
  tmp.AndNot(bs2);
  EXPECT_EQ(tmp_andnot, tmp);
  std::cout << tmp << std::endl;

  tmp = bs1;
  c = tmp.AndNotCount(bs2);
  EXPECT_EQ(bs1, tmp);
  EXPECT_EQ((std::size_t)3, c);
  c = tmp.AndNotCountUpdate(bs2);
  EXPECT_EQ(tmp_andnot, tmp);
  EXPECT_EQ((std::size_t)3, c);

  // bs1 flip
  VariableBitset<uint64> tmp_flip(
      "011110"              // 4
      "01111111" "11111111" // 15
      "11110011" "11111001" // 12
      "11001111" "11111110" // 13
      "11111111" "11111100" // 14
                             );

  tmp = bs1;
  tmp.Flip();
  EXPECT_EQ(tmp_flip, tmp);
  std::cout << tmp << std::endl;

  tmp = bs1;
  c = tmp.FlipCount();
  EXPECT_EQ(bs1, tmp);
  EXPECT_EQ((std::size_t)58, c);
  c = tmp.FlipCountUpdate();
  EXPECT_EQ(tmp_flip, tmp);
  EXPECT_EQ((std::size_t)58, c);

  EXPECT_TRUE (bs1.IsSubsetOf(bs3));
  EXPECT_FALSE(bs1.IsSupersetOf(bs3));

  EXPECT_FALSE(bs3.IsSubsetOf(bs1));
  EXPECT_TRUE (bs3.IsSupersetOf(bs1));

  EXPECT_TRUE (bs1.IsSubsetOf(bs1));
  EXPECT_TRUE (bs1.IsSupersetOf(bs1));

  EXPECT_TRUE (bs3.IsSubsetOf(bs3));
  EXPECT_TRUE (bs3.IsSupersetOf(bs3));

  std::size_t pos=0;
  pos = bs2.FindFirst();   EXPECT_EQ( 8u, pos);
  pos = bs2.FindNext(pos); EXPECT_EQ( 9u, pos);
  pos = bs2.FindNext(pos); EXPECT_EQ(16u, pos);
  pos = bs2.FindNext(pos); EXPECT_EQ(17u, pos);
  pos = bs2.FindNext(pos); EXPECT_EQ(28u, pos);
  pos = bs2.FindNext(pos); EXPECT_EQ(29u, pos);

  pos = bs2.FindNext(69u); EXPECT_EQ(bs2.nu_bits, pos);
}

TEST (VariableBitsetArrayTest, InitTest) {
  VariableBitsetArray<uint64> vba(70, 10);

  std::cout << vba.N(0) << std::endl;

  uint64 *a, *b, *c;
  a = vba.N(0);
  b = vba.N(1);
  c = vba.N(2);

  EXPECT_EQ(a+2, b);
  EXPECT_EQ(a+4, c);

  std::cout << std::hex << a << std::endl;
  std::cout << std::hex << b << std::endl;
  std::cout << std::hex << c << std::endl;

  VariableBitset<uint64> bs0(
      "100001" // 2
      "10000000" "00000000" // 1
      "00001100" "00000110" // 4
      "00110000" "00000001" // 3
      "00000000" "00000011" // 2
                             );

  VariableBitset<uint64> bs1(
      "100001"
      "01111100" "00000000"
      "11111111" "00001110"
      "00110000" "00000011"
      "00000011" "00000000"
                             );

  VariableBitset<uint64> bs2(
      "111001"
      "10011000" "11111111"
      "00001100" "00000110"
      "00110011" "01111001"
      "11111100" "00000011"
                             );

  vba.Copy(bs0, vba.N(0));
  vba.Copy(bs1, vba.N(1));
  vba.Copy(bs2, vba.N(2));
  vba.Copy(vba.N(0), vba.N(3));
  vba.Copy(vba.N(1), vba.N(4));
  vba.Copy(vba.N(2), vba.N(5));
  vba.Copy(vba.N(0), vba.N(6));
  vba.Copy(vba.N(1), vba.N(7));
  vba.Copy(vba.N(2), vba.N(8));
  vba.Copy(bs0, vba.N(9));
  
  std::cout << bs0 << std::endl;
  vba.ShowNth(std::cout, vba.N(0)); std::cout << std::endl;
  vba.ShowNth(std::cout, vba.N(3)); std::cout << std::endl;

  std::cout << bs1 << std::endl;
  vba.ShowNth(std::cout, vba.N(1)); std::cout << std::endl;
  vba.ShowNth(std::cout, vba.N(4)); std::cout << std::endl;

  std::cout << bs2 << std::endl;
  vba.ShowNth(std::cout, vba.N(2)); std::cout << std::endl;
  vba.ShowNth(std::cout, vba.N(5)); std::cout << std::endl;

  EXPECT_TRUE(vba.Test(vba.N(0), 0u ));
  EXPECT_TRUE(vba.Test(vba.N(0), 1u ));
  EXPECT_TRUE(vba.Test(vba.N(0), 64u));
  EXPECT_TRUE(vba.Test(vba.N(0), 69u));

  EXPECT_FALSE(vba.Test(vba.N(0), 3u ));
  EXPECT_FALSE(vba.Test(vba.N(0), 68u));

  EXPECT_TRUE(vba.IsEqualTo(vba.N(0), vba.N(3)));
  EXPECT_FALSE(vba.IsEqualTo(vba.N(0), vba.N(4)));
}

TEST (VariableBitsetArrayTest, SetResetTest) {
  VariableBitsetArray<uint64> vba(128, 10);

  vba.Set(vba.N(0));
  EXPECT_EQ(128u, vba.Count(vba.N(0)));
  vba.Reset(vba.N(0));
  EXPECT_EQ(  0u, vba.Count(vba.N(0)));
}

TEST (VariableBitsetArrayTest, OperationTest) {
  VariableBitsetArray<uint64> vba(70, 10);

  std::cout << vba.N(0) << std::endl;

  uint64 *a0, *a1, *a2;
  a0 = vba.N(0);
  a1 = vba.N(1);
  a2 = vba.N(2);

  EXPECT_EQ(a0+2, a1);
  EXPECT_EQ(a0+4, a2);

  std::cout << std::hex << a0 << std::endl;
  std::cout << std::hex << a1 << std::endl;
  std::cout << std::hex << a2 << std::endl;

  VariableBitset<uint64> bs0(
      "100001" // 2
      "10000000" "00000000" // 1
      "00001100" "00000110" // 4
      "00110000" "00000001" // 3
      "00000000" "00000011" // 2
                             );

  VariableBitset<uint64> bs1(
      "100001"
      "01111100" "00000000"
      "11111111" "00001110"
      "00110000" "00000011"
      "00000011" "00000000"
                             );

  VariableBitset<uint64> bs2(
      "111001"
      "10011000" "11111111"
      "00001100" "00000110"
      "00110011" "01111001"
      "11111100" "00000011"
                             );

  vba.Set(vba.N(0));
  EXPECT_EQ(70u, vba.Count(vba.N(0)));
  vba.Reset(vba.N(0));
  EXPECT_EQ( 0u, vba.Count(vba.N(0)));

  vba.Copy(bs0, vba.N(0));
  vba.Copy(bs1, vba.N(1));
  vba.Copy(bs2, vba.N(2));
  vba.Copy(vba.N(0), vba.N(3));
  vba.Copy(vba.N(1), vba.N(4));
  vba.Copy(vba.N(2), vba.N(5));
  vba.Copy(vba.N(0), vba.N(6));
  vba.Copy(vba.N(1), vba.N(7));
  vba.Copy(vba.N(2), vba.N(8));
  vba.Copy(bs0, vba.N(9));

  VariableBitset<uint64> bs_tmp(70);
  bs_tmp.CopyFrom(vba.N(3));
  EXPECT_EQ(bs0, bs_tmp);
  EXPECT_TRUE(vba.IsEqualTo(bs_tmp, vba.N(6)));

  bs_tmp.Set();
  std::cout << bs_tmp << std::endl;
  EXPECT_EQ(70u, bs_tmp.Count());
  bs_tmp.Reset();
  std::cout << bs_tmp << std::endl;
  EXPECT_EQ( 0u, bs_tmp.Count());

  uint64 * tmp;
  std::size_t c;

  // bs0 and bs1
  VariableBitset<uint64> tmp_and(
      "100001"
      "00000000" "00000000"
      "00001100" "00000110"
      "00110000" "00000001"
      "00000000" "00000000"
                             );

  // bs0 or bs1
  VariableBitset<uint64> tmp_or(
      "100001"              // 2
      "11111100" "00000000" // 6
      "11111111" "00001110" // 11
      "00110000" "00000011" // 4
      "00000011" "00000011" // 4
                             );

  // bs0 xor bs1
  VariableBitset<uint64> tmp_xor(
      "000000"              // 0
      "11111100" "00000000" // 6
      "11110011" "00001000" // 7
      "00000000" "00000010" // 1
      "00000011" "00000011" // 4
                             );

  // bs0 and_not bs0
  VariableBitset<uint64> tmp_andnot(
      "000000"              // 0
      "10000000" "00000000" // 1
      "00000000" "00000000" // 0
      "00000000" "00000000" // 0
      "00000000" "00000011" // 2
                             );

  // bs0 flip
  VariableBitset<uint64> tmp_flip(
      "011110"              // 4
      "01111111" "11111111" // 15
      "11110011" "11111001" // 12
      "11001111" "11111110" // 13
      "11111111" "11111100" // 14
                             );

  vba.Copy(bs0, vba.N(9));
  tmp = vba.N(9);
  vba.And(bs1, tmp);
  EXPECT_TRUE(vba.IsEqualTo(tmp_and, tmp));
  vba.Copy(tmp, vba.N(8));
  EXPECT_TRUE(vba.IsEqualTo(tmp, vba.N(8)));

  vba.Copy(bs0, vba.N(9));
  tmp = vba.N(9);
  c = vba.AndCount(bs1, tmp);
  EXPECT_TRUE(vba.IsEqualTo(bs0, tmp));
  EXPECT_EQ((std::size_t)9, c);
  c = vba.AndCountUpdate(bs1, tmp);
  EXPECT_TRUE(vba.IsEqualTo(tmp_and, tmp));
  EXPECT_EQ((std::size_t)9, c);

  vba.Copy(bs0, vba.N(9));
  tmp = vba.N(9);
  vba.Or(bs1, tmp);
  EXPECT_TRUE(vba.IsEqualTo(tmp_or, tmp));
  vba.Copy(tmp, vba.N(8));
  EXPECT_TRUE(vba.IsEqualTo(tmp, vba.N(8)));

  vba.Copy(bs0, vba.N(9));
  tmp = vba.N(9);
  c = vba.OrCount(bs1, tmp);
  EXPECT_TRUE(vba.IsEqualTo(bs0, tmp));
  EXPECT_EQ((std::size_t)27, c);
  c = vba.OrCountUpdate(bs1, tmp);
  EXPECT_TRUE(vba.IsEqualTo(tmp_or, tmp));
  EXPECT_EQ((std::size_t)27, c);

  vba.Copy(bs0, vba.N(9));
  tmp = vba.N(9);
  vba.Xor(bs1, tmp);
  EXPECT_TRUE(vba.IsEqualTo(tmp_xor, tmp));
  vba.Copy(tmp, vba.N(8));
  EXPECT_TRUE(vba.IsEqualTo(tmp, vba.N(8)));

  vba.Copy(bs0, vba.N(9));
  tmp = vba.N(9);
  c = vba.XorCount(bs1, tmp);
  EXPECT_TRUE(vba.IsEqualTo(bs0, tmp));
  EXPECT_EQ((std::size_t)18, c);
  c = vba.XorCountUpdate(bs1, tmp);
  EXPECT_TRUE(vba.IsEqualTo(tmp_xor, tmp));
  EXPECT_EQ((std::size_t)18, c);

  vba.Copy(bs0, vba.N(9));
  tmp = vba.N(9);
  vba.AndNot(bs1, tmp);
  EXPECT_TRUE(vba.IsEqualTo(tmp_andnot, tmp));
  vba.Copy(tmp, vba.N(8));
  EXPECT_TRUE(vba.IsEqualTo(tmp, vba.N(8)));

  vba.Copy(bs0, vba.N(9));
  tmp = vba.N(9);
  c = vba.AndNotCount(bs1, tmp);
  EXPECT_TRUE(vba.IsEqualTo(bs0, tmp));
  EXPECT_EQ((std::size_t)3, c);
  c = vba.AndNotCountUpdate(bs1, tmp);
  EXPECT_TRUE(vba.IsEqualTo(tmp_andnot, tmp));
  EXPECT_EQ((std::size_t)3, c);

  vba.Copy(bs0, vba.N(9));
  tmp = vba.N(9);
  vba.Flip(tmp);
  EXPECT_TRUE(vba.IsEqualTo(tmp_flip, tmp));
  vba.Copy(tmp, vba.N(8));
  EXPECT_TRUE(vba.IsEqualTo(tmp, vba.N(8)));

  vba.Copy(bs0, vba.N(9));
  tmp = vba.N(9);
  c = vba.FlipCount(tmp);
  EXPECT_TRUE(vba.IsEqualTo(bs0, tmp));
  EXPECT_EQ((std::size_t)58, c);
  c = vba.FlipCountUpdate(tmp);
  EXPECT_TRUE(vba.IsEqualTo(tmp_flip, tmp));
  EXPECT_EQ((std::size_t)58, c);

  // vba.Copy(bs0, vba.N(0));
  // vba.Copy(bs1, vba.N(1));
  // vba.Copy(bs2, vba.N(2));
  // vba.Copy(vba.N(0), vba.N(3));
  // vba.Copy(vba.N(1), vba.N(4));
  // vba.Copy(vba.N(2), vba.N(5));

  std::cout << "bs0 :" << bs0 << std::endl;
  std::cout << "vba3:";
  vba.ShowNth(std::cout, vba.N(3)); std::cout << std::endl;

  std::cout << "bs2 :" << bs2 << std::endl;
  std::cout << "vba5:";
  vba.ShowNth(std::cout, vba.N(5)); std::cout << std::endl;

  // subset, superset

  // bs0 (vba3) \in bs2 (vba5)
  EXPECT_TRUE (vba.IsSubsetOf(bs0, vba.N(5)));
  EXPECT_FALSE(vba.IsSupersetOf(bs0, vba.N(5)));

  EXPECT_TRUE (vba.IsSubsetOf(vba.N(3), vba.N(5)));
  EXPECT_FALSE(vba.IsSupersetOf(vba.N(3), vba.N(5)));

  // bs0 (vba3) \in bs1
  EXPECT_TRUE (vba.IsSubsetOf(bs0, vba.N(3)));
  EXPECT_TRUE (vba.IsSupersetOf(bs0, vba.N(3)));

  EXPECT_TRUE (vba.IsSubsetOf(vba.N(5), vba.N(5)));
  EXPECT_TRUE (vba.IsSupersetOf(vba.N(5), vba.N(5)));

  // find first, find next
  std::size_t pos=0;

  pos = vba.FindFirst(vba.N(1));     EXPECT_EQ( 8u, pos);
  pos = vba.FindNext(vba.N(1), pos); EXPECT_EQ( 9u, pos);
  pos = vba.FindNext(vba.N(1), pos); EXPECT_EQ(16u, pos);
  pos = vba.FindNext(vba.N(1), pos); EXPECT_EQ(17u, pos);
  pos = vba.FindNext(vba.N(1), pos); EXPECT_EQ(28u, pos);
  pos = vba.FindNext(vba.N(1), pos); EXPECT_EQ(29u, pos);

  pos = vba.FindNext(vba.N(1), 69u); EXPECT_EQ(vba.nu_bits, pos);
}

TEST (VariableBitsetHelperTest, HelperInitTest) {
  VariableBitsetHelper<uint64> ac(70);
  uint64 * el = ac.New();
  uint64 * ar = ac.NewArray(10);

  // bs0
  ac.InitFromString(
      "100001" // 2
      "10000000" "00000000" // 1
      "00001100" "00000110" // 4
      "00110000" "00000001" // 3
      "00000000" "00000011" // 2
      , el);

  // bs1
  ac.InitFromString(
      "100001"
      "01111100" "00000000"
      "11111111" "00001110"
      "00110000" "00000011"
      "00000011" "00000000"
      , ac.N(ar, 0));

  // bs2
  ac.InitFromString(
      "111001"
      "10011000" "11111111"
      "00001100" "00000110"
      "00110011" "01111001"
      "11111100" "00000011"
      , ac.N(ar, 1));

  
  // bs0 and bs1
  uint64 * tmp_and = ac.New();
  ac.InitFromString(
      "100001"
      "00000000" "00000000"
      "00001100" "00000110"
      "00110000" "00000001"
      "00000000" "00000000"
      , tmp_and);

  // bs0 or bs1
  uint64 * tmp_or = ac.New();
  ac.InitFromString(
      "100001"              // 2
      "11111100" "00000000" // 6
      "11111111" "00001110" // 11
      "00110000" "00000011" // 4
      "00000011" "00000011" // 4
      , tmp_or);

  // bs0 xor bs1
  uint64 * tmp_xor = ac.New();
  ac.InitFromString(
      "000000"              // 0
      "11111100" "00000000" // 6
      "11110011" "00001000" // 7
      "00000000" "00000010" // 1
      "00000011" "00000011" // 4
      , tmp_xor);

  // bs0 and_not bs0
  uint64 * tmp_andnot = ac.New();
  ac.InitFromString(
      "000000"              // 0
      "10000000" "00000000" // 1
      "00000000" "00000000" // 0
      "00000000" "00000000" // 0
      "00000000" "00000011" // 2
      , tmp_andnot);

  // bs0 flip
  uint64 * tmp_flip = ac.New();
  ac.InitFromString(
      "011110"              // 4
      "01111111" "11111111" // 15
      "11110011" "11111001" // 12
      "11001111" "11111110" // 13
      "11111111" "11111100" // 14
      , tmp_flip);

  ac.Copy(ac.N(ar, 1), ac.N(ar, 2));

  EXPECT_FALSE( ac.IsEqualTo(el, ac.N(ar, 1)) );
  EXPECT_FALSE( ac.IsEqualTo(ac.N(ar, 0), ac.N(ar, 1)) );
  EXPECT_TRUE ( ac.IsEqualTo(ac.N(ar, 1), ac.N(ar, 2)) );

  std::size_t count;

  ac.Print(std::cout, el); std::cout << std::endl;
  ac.Print(std::cout, ac.N(ar, 0)); std::cout << std::endl;
  ac.Print(std::cout, ac.N(ar, 1)); std::cout << std::endl;
  ac.Print(std::cout, ac.N(ar, 2)); std::cout << std::endl;
  ac.Print(std::cout, ac.N(ar, 3)); std::cout << std::endl;

  std::cout << "and\n";

  ac.Copy( ac.N(ar, 0), ac.N(ar, 5) ); // copy bs1 to ar[4]
  count = ac.AndCount( el, ac.N(ar, 5) ); // (ar[5] & bs0).Count()
  EXPECT_TRUE( ac.IsEqualTo(ac.N(ar, 0), ac.N(ar, 5)) ); // unchanged
  EXPECT_EQ((std::size_t)9, count);
  ac.And( el, ac.N(ar, 5) ); // ar[5] &= bs0
  ac.Print(std::cout, ac.N(ar, 5)); std::cout << std::endl;
  EXPECT_TRUE ( ac.IsEqualTo(tmp_and, ac.N(ar, 5)) );
  ac.Copy( ac.N(ar, 0), ac.N(ar, 5) ); // copy bs1 to ar[4]
  count = ac.AndCountUpdate( el, ac.N(ar, 5) ); // ar[5] &= bs0
  EXPECT_EQ((std::size_t)9, count);
  ac.Print(std::cout, ac.N(ar, 5)); std::cout << std::endl;
  EXPECT_TRUE ( ac.IsEqualTo(tmp_and, ac.N(ar, 5)) );

  std::cout << "or\n";
  ac.Copy( ac.N(ar, 0), ac.N(ar, 4) ); // copy bs1 to ar[4]
  count = ac.OrCount( el, ac.N(ar, 4) ); // (ar[4] | bs0).Count()
  EXPECT_TRUE( ac.IsEqualTo(ac.N(ar, 0), ac.N(ar, 4)) ); // unchanged
  EXPECT_EQ((std::size_t)27, count);
  ac.Or( el, ac.N(ar, 4) ); // ar[4] |= bs0
  ac.Print(std::cout, ac.N(ar, 4)); std::cout << std::endl;
  EXPECT_TRUE ( ac.IsEqualTo(tmp_or, ac.N(ar, 4)) );
  ac.Copy( ac.N(ar, 0), ac.N(ar, 4) ); // copy bs1 to ar[4]
  count = ac.OrCountUpdate( el, ac.N(ar, 4) ); // ar[4] |= bs0
  EXPECT_EQ((std::size_t)27, count);
  ac.Print(std::cout, ac.N(ar, 4)); std::cout << std::endl;
  EXPECT_TRUE ( ac.IsEqualTo(tmp_or, ac.N(ar, 4)) );

  std::cout << "xor\n";
  ac.Copy( ac.N(ar, 0), ac.N(ar, 6) ); // copy bs1 to ar[6]
  count = ac.XorCount( el, ac.N(ar, 6) ); // (ar[6] ^ bs0).Count()
  EXPECT_TRUE( ac.IsEqualTo(ac.N(ar, 0), ac.N(ar, 6)) ); // unchanged
  EXPECT_EQ((std::size_t)18, count);
  ac.Xor( el, ac.N(ar, 6) ); // ar[6] ^= bs0
  ac.Print(std::cout, ac.N(ar, 6)); std::cout << std::endl;
  EXPECT_TRUE ( ac.IsEqualTo(tmp_xor, ac.N(ar, 6)) );
  ac.Copy( ac.N(ar, 0), ac.N(ar, 6) ); // copy bs1 to ar[6]
  count = ac.XorCountUpdate( el, ac.N(ar, 6) ); // ar[6] ^= bs0
  EXPECT_EQ((std::size_t)18, count);
  ac.Print(std::cout, ac.N(ar, 6)); std::cout << std::endl;
  EXPECT_TRUE ( ac.IsEqualTo(tmp_xor, ac.N(ar, 6)) );

  std::cout << "flip\n";
  ac.Copy( el, ac.N(ar, 7) ); // copy bs1 to ar[7]
  count = ac.FlipCount( ac.N(ar, 7) ); // (~ar[7]).Count()
  EXPECT_TRUE( ac.IsEqualTo(el, ac.N(ar, 7)) ); // unchanged
  EXPECT_EQ((std::size_t)58, count);
  ac.Flip( ac.N(ar, 7) ); // ar[7] = ~ar[7]
  ac.Print(std::cout, ac.N(ar, 7)); std::cout << std::endl;
  EXPECT_TRUE ( ac.IsEqualTo(tmp_flip, ac.N(ar, 7)) );
  ac.Copy( el, ac.N(ar, 7) ); // copy bs1 to ar[7]
  count = ac.FlipCountUpdate( ac.N(ar, 7) ); // ar[7] = ~ar[7]
  EXPECT_EQ((std::size_t)58, count);
  ac.Print(std::cout, ac.N(ar, 7)); std::cout << std::endl;
  EXPECT_TRUE ( ac.IsEqualTo(tmp_flip, ac.N(ar, 7)) );

  // todo: IsSubsetOf, IsSupersetOf
  // todo: FindFirst, FindNext

  ac.Delete(el);
  ac.Delete(ar);
  ac.Delete(tmp_and);
  ac.Delete(tmp_or);
  ac.Delete(tmp_xor);
  ac.Delete(tmp_andnot);
  ac.Delete(tmp_flip);
}
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
