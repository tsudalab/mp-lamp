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
 * File:   Functions4chi.cpp
 * Author: miura
 * 
 * This source includes calculate P-value and MASL of the chi-square test.
 * 
 * Created on 2016/01/06, 17:43
 */

#include <vector>

#include "Functions4chi.h"

/**
 * Constructor
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
Functions4chi::Functions4chi(int alternative) :
    FunctionsSuper(std::abs(alternative)) {
  reverse = (alternative < 0);
}

/**
 * Constructor
 * @param transaction_size all transaction size
 * @param n1_count transaction size which have flag = 1 (n1)
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
Functions4chi::Functions4chi(int transaction_size, int n1_count, int alternative) :
    FunctionsSuper(transaction_size, n1_count, std::abs(alternative)) {
  reverse = (alternative < 0);
}

/**
 * Destructor
 */
Functions4chi::~Functions4chi() {
}

/**
 * This function calculates the minimum p-value which support size is x.
 * @param x
 * @return 
 */
double Functions4chi::funcF(int x) {
  if (__minpvalTable.find(x) != __minpvalTable.end()) {
    return __minpvalTable[x];
  } else {
    double p1 = 1.0f, p2 = 1.0f;
    double chi1 = 0.0f, chi2 = 0.0f;
    double total_row1 = __f_size;
    double total = __t_size;
    // when x < n_u
    if (x < total_row1) {
      double ovalues[2][2] = {
        {(double) x, 0},
        {total_row1 - (double) x, total - total_row1}};
      chi1 = __probabilityTable(ovalues);
      p1 = __chi2pval(chi1);
      ovalues[0][0] = 0;
      ovalues[0][1] = (double) x;
      ovalues[1][0] = total_row1;
      ovalues[1][1] = total - total_row1 - (double) x;
      chi2 = __probabilityTable(ovalues);
      p2 = __chi2pval(chi2);
    }			// when x >= n_u
    else {
      double ovalues[2][2] = {
        {total_row1, (double) x - total_row1},
        {0, total - (double) x}};
      chi1 = __probabilityTable(ovalues);
      p1 = __chi2pval(chi1);
      ovalues[0][0] = 0;
      ovalues[0][1] = (double) x;
      ovalues[1][0] = total_row1;
      ovalues[1][1] = total - total_row1 - (double) x;
      chi2 = __probabilityTable(ovalues);
      p2 = __chi2pval(chi2);
    }
    if (alternative == 0) {
      p1 = std::min(p1 * 2.0, 1.0);
      p2 = std::min(p2 * 2.0, 1.0);
    }
    if (p1 < p2) {
      __minpvalTable[x] = p1;
      return p1;
    } else {
      __minpvalTable[x] = p2;
      return p2;
    }
  }
}

/**
 * Calculate probability of occurrence probability about table.
 * @param ovalues
 * @return 
 */
double Functions4chi::__probabilityTable(const double (&ovalues)[2][2]) {
  double means[2][2] = {
    {0, 0},
    {0, 0}};
  __calMeans(ovalues, means); // calculate the exception value

  // Yate continuity correction
  double yate_corr = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      if (means[i][j] < 5) {
        yate_corr = 0.5;
        break;
      }
    }
  }

  double chi = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      double row = ovalues[i][j];
      double mean = means[i][j];
      double po = (std::abs(row - mean) - yate_corr);
      chi += (po * po) / mean;
    }
  }
  return chi;
}

/**
 * Calculate means
 * @param ovalues
 * @param means
 */
void Functions4chi::__calMeans(const double (&ovalues)[2][2], double (&means)[2][2]) {
  double total = __t_size;
  double total_col1 = __f_size; // the number of all flag 1 transaction (n1)
  double total_col2 = total - total_col1; // the number of all flag 0 transactio (n0)
  double total_row1 = ovalues[0][0] + ovalues[0][1];
  double total_row2 = ovalues[1][0] + ovalues[1][1];
  means[0][0] = (total_row1 * total_col1) / total;
  means[0][1] = (total_row1 * total_col2) / total;
  means[1][0] = (total_row2 * total_col1) / total;
  means[1][1] = (total_row2 * total_col2) / total;
}

/**
 * 
 * @param chi
 * @return 
 */
double Functions4chi::__chi2pval(double chi) {
  if (chi == 0.0) {
    return 1.0;
  } else { // dimension = 1
    return stdNorDistribution(std::pow(chi, 0.5));
  }
}

/**
 * Calculate p-value by using chi-square test.
 * @param ovalues
 * @param score return value
 * @return 
 */
double Functions4chi::calPValue(const int sup, const int pos_sup, const std::vector<double>& tgroup_x, const std::vector<double>& tgroup_y) {
  (void)tgroup_x; (void)tgroup_y;

  double total_row1 = sup; //sum( ovalues[0] );
  double p = __pvalTable.getValue(sup, pos_sup);
  double chi = __chiTable.getValue(sup, pos_sup);
  if (p < 0) { // calculate P-value and save to the table
    double ovalues[2][2];
    ovalues[0][0] = pos_sup;
    ovalues[0][1] = sup - pos_sup;
    ovalues[1][0] = __f_size - pos_sup;
    ovalues[1][1] = __t_size - sup - ovalues[1][0];
    chi = __probabilityTable(ovalues);
    p = __chi2pval(chi);
    if (0 < alternative) {
      if (pos_sup < (std::min((double) __f_size, total_row1) / 2.0))
        p = 1. - p;
    }			// when the alternative hypothesis is "two.sided", 
    // the P-value is doubled. 
    else {
      p = std::min(p * 2., 1.0);
    }
    __pvalTable.putValue(sup, pos_sup, p);
    __chiTable.putValue(sup, pos_sup, chi);
    ++calTime;
  }
  return p;
}
