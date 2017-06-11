/*
 * Copyright (c) 2013, LAMP development team
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the LAMP development team nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL LAMP DEVELOPMENT TEAM BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
/* 
 * File:   FunctionsSuper.cpp
 * Author: miura
 * 
 * Created on 2016/01/06, 15:28
 */

#include "FunctionsSuper.h"

/**
 * Constructor
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
FunctionsSuper::FunctionsSuper(int alternative) :
    alternative(alternative), calTime(0) {
  __t_size = -1;
  __f_size = -1;
  __minpvalTable[0] = 1.0;
}

/**
 * Constructor
 * @param transaction_size all transaction size
 * @param n1_count transaction size which have flag = 1 (n1)
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
FunctionsSuper::FunctionsSuper(int transaction_size, int n1_count, int alternative) :
    alternative(alternative), calTime(0) {
  __t_size = transaction_size;
  __f_size = n1_count;
  __minpvalTable[0] = 1.0;
}

/**
 * Destructor
 */
FunctionsSuper::~FunctionsSuper() {
}

/**
 * Calculate probability of standard normal distribution.
 * this function returns the probability of one-sided test.
 * @param x
 * @return 
 */
double FunctionsSuper::stdNorDistribution(double x) {
  double is_value = -1;
  double y = std::abs(x);
  double c = y * y;
  double p = 0.0f;
  double z = std::exp(-c * 0.5) * pi2;
  if (y < 2.5) {
    for (double i = 20.0; 0 < i; i -= 1.0f) {
      p = i * c / (i * 2.0 + 1.0 + is_value * p);
      is_value = -is_value;
    }
    p = 0.5 - z * y / (1.0 - p);
  } else {
    for (double i = 20.0; 0 < i; i -= 1.0f) {
      p = i / (y + p);
    }
    p = z / (y + p);
  }
  return p;
}
