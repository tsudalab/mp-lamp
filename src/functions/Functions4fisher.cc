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
 * File:   Functions4fisher.cpp
 * Author: miura
 * 
 * This class calculate function f that define in paper. (C(n1, x)/C(n0+n1, x))
 * 
 * Created on 2016/01/06, 16:23
 */

#include <string>
#include "Functions4fisher.h"

/**
 * Constructor
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
Functions4fisher::Functions4fisher(int alternative) :
    FunctionsSuper(std::abs(alternative)) {
  reverse = (alternative < 0);
}

/**
 * Constructor
 * @param transaction_size all transaction size
 * @param n1_count transaction size which have flag = 1 (n1)
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
Functions4fisher::Functions4fisher(int transaction_size, int n1_count, int alternative) :
    FunctionsSuper(transaction_size, n1_count, std::abs(alternative)) {
  reverse = (alternative < 0);
}

/**
 * Destructor
 */
Functions4fisher::~Functions4fisher() {
}

/**
 * This function calculates the minimum p-value which support size is x.
 * @param x support size
 * @return 
 */
double Functions4fisher::funcF(int x) {
  if (__minpvalTable.find(x) != __minpvalTable.end()) {
    return __minpvalTable[x];
  } else {
    double minp_i_log = 0.0;

    unsigned int uplim = x;
    if (x > __f_size) uplim = __f_size;
    for (unsigned int i = 0; i < uplim; ++i) {
      minp_i_log += log((double) (__f_size - i));
      minp_i_log -= log((double) (__t_size - i));
    }
    minp_i_log = exp(minp_i_log);
    __minpvalTable[x] = minp_i_log;
    return minp_i_log;
  }
}

/**
 * Calculate probability of occurrence probability about table.
 * @param x total of the first row (the number of targetting gene)
 * @param a the top-left of the table
 * @return 
 */
double Functions4fisher::__probability(int x, int a) {
  double p = __occrTable.getValue(x, a);
  if (p < 0) {
    int n = __t_size;
    int n1 = __f_size;
    int n0 = n - n1;
    int b = x - a;
    int i = 0;
    p = 1.0;
    while (i < a) {
      p *= (double) (n1 - i) / (double) (a - i); // c(n1, a)
      p *= (double) (x - i) / (double) (n - i); // c(n1+n0, x)
      i++;
    }
    i = 0;
    while (i < b) {
      p *= (double) (n0 - i) / (double) (b - i); // c(n0, b)
      double minus_denominator = a + i;
      p *= (double) (x - minus_denominator) / (double) (n - minus_denominator); // c(n1+n0, x)
      i++;
    }
    __occrTable.putValue(x, a, p);
  }
  return p;
}

/**
 * Calculate p-value by using fisher's exact test.
 * @param ovalues
 * @param score return value
 * @return 
 */
double Functions4fisher::calPValue(const int sup, const int pos_sup,
                                   const std::vector<double>& tgroup_x, const std::vector<double>& tgroup_y) {
  (void)tgroup_x; (void)tgroup_y;

  int total_col1 = __f_size;
  int total_row1 = sup; //sum( ovalues[0] );
  double p = __pvalTable.getValue(total_row1, pos_sup);
  if (p < 0) { // calculate P-value and save to the table
    double p0 = __probability(total_row1, pos_sup);
    p = p0;
    int pos_max = std::min(total_row1, total_col1);
    // when the alternative hypothesis is "two.sided",
    // the lower case probability is cumulated. 
    if (alternative < 1) {
      int a = 0;
      // cumulate the lower case probability.
      while (a < pos_sup) {
        double pa = __probability(total_row1, a);
        if (1.E-16 < pa - p0) // pa > p0
          break;
        p += pa;
        a++;
      }
      // cumulate the upper case probability.
      a = pos_max;
      while (pos_sup < a) {
        double pa = __probability(total_row1, a);
        if (1.E-16 < pa - p0) // pa > p0
          break;
        p += pa;
        a--;
      }
    }			// when the alternative hypothesis is "greater" or "less",
    // the higher/less case probability is cumulated.  
    else {
      int a = pos_sup + 1;
      while (a <= pos_max) {
        double pa = __probability(total_row1, a);
        p += pa;
        a++;
      }
    }
    __pvalTable.putValue(total_row1, pos_sup, p);
    calTime++;
  }
  return p;
}
