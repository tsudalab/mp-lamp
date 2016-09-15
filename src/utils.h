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

#ifndef _LAMP_SEARCH_UTILS_H_
#define _LAMP_SEARCH_UTILS_H_

#include <stdint.h>
#include <exception>
#include <stdexcept>
#include <boost/array.hpp>

namespace lamp_search {

static const double KILO = 1000.0;
static const double MEGA = 1000.0*1000;
static const double GIGA = 1000.0*1000*1000;
static const double TERA = 1000.0*1000*1000*1000;

typedef unsigned char uint8;
typedef char int8;

typedef unsigned long long int uint64;
typedef long long int int64;

#ifdef OPTLOG
#  define LOG(x) do {} while (0)
#else
#  define LOG(x) do { x; } while (0)
#endif

#ifdef NDEBUG
#  define DBG(x) do {} while (0)
#else
#  define DBG(x) do { x; } while (0)
#endif

const int kMaxSearchDepth = 30; // is this enough?
const int kMaxItemSet = 100; // is this enough?

//typedef std::bitset<128> hash_type;
// typedef __uint128_t hash_type;
// typedef __uint128_t uint128;

class hash_type {
 public:
  boost::array<uint64, 2> v_;

  void set(uint64 u, uint64 l) { v_[1] = u; v_[0] = l; }
  void get(uint64 *u, uint64 *l) { *u = v_[1]; *l = v_[0]; }
};

inline int operator%(hash_type h, int i) {
  return (h.v_[0] % i);
}

} // namespace lamp_search

#endif // _LAMP_SEARCH_UTILS_H_

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
