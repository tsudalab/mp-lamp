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
#include <sstream>
#include <limits>

#include <boost/integer/static_log2.hpp>

#include "utils.h"
#include "variable_bitset_array.h"

namespace lamp_search {

template<typename Block>
const std::size_t VariableBitsetTraits<Block>::bits_per_block =
    std::numeric_limits<Block>::digits;

template<typename Block>
const std::size_t VariableBitsetTraits<Block>::bytes_per_block = sizeof(Block);

// template<>
// const std::size_t VariableBitsetTraits<uint64>::block_index_shift = 6u;
// template<>
// const std::size_t VariableBitsetTraits<uint64>::block_index_mask = 0x3fu;

// template<>
// const std::size_t VariableBitsetTraits<unsigned int>::block_index_shift = 5u;
// template<>
// const std::size_t VariableBitsetTraits<unsigned int>::block_index_mask = 0x1fu;

template class VariableBitsetTraits<uint64>;
template class VariableBitset<uint64>;
template class VariableBitsetArray<uint64>;

} // namespace lamp_search
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
