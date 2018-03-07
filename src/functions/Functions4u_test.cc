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
 * File:   Functions4u_test.cpp
 * Author: miura
 * 
 * This class calculate function f that means minimum p-value (MASL).
 * 
 * Created on 2016/01/06, 17:27
 */

#include <vector>
#include <algorithm>

#include "Functions4u_test.h"

/**
 * Constructor
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
Functions4u_test::Functions4u_test(int alternative) :
    FunctionsSuper(alternative) {
  reverse = (alternative < 0);
}

/**
 * Constructor
 * @param transaction_size all transaction size
 * @param n1_count transaction size which have flag = 1 (n1)
 * @param alternative hypothesis, 1 -> greater, 0 -> two sided, -1 -> less
 */
Functions4u_test::Functions4u_test(int transaction_size, int n1_count, int alternative) :
    FunctionsSuper(transaction_size, n1_count, alternative) {
  reverse = (alternative < 0);
}

/**
 * Destructor
 */
Functions4u_test::~Functions4u_test() {
}

/**
 * This function calculates the minimum p-value which support size is x.
 * That is, calculates MASL.
 * The z-value that minimum p-value is mean/var
 * @param x
 * @return 
 */
double Functions4u_test::funcF(int x) {
  if (__minpvalTable.find(x) != __minpvalTable.end()) {
    return __minpvalTable[x];
  } else {
    // calculate p-value if the group_x is consisted from max x transactions.

    // mean and variance which is used in U test.
    double size_x = x;
    if (__t_size < size_x)
      size_x = __t_size;
    double size_y = __t_size - x;
    if (size_y < 0)
      size_y = 0;
    double mean_u = (size_x * size_y) / 2;
    double var_u = size_x * size_y * (size_x + size_y + 1) / 12;

    double min_z = mean_u / std::sqrt(var_u); // minimum z-value limited x.
    double p = stdNorDistribution(min_z); // p-value if transaction divided into max x and other.
    __minpvalTable[x] = p;
    return p;
  }
}

/**
 * calculate p-value by u test.
 * @param tgroup_x test group 1. This is consisted of transaction values.
 * @param tgroup_y test group 2. This is consisted of transaction values.
 * @param z_value return value
 * @param statistic return value
 * @return
 */
double Functions4u_test::calPValue(const std::vector<double>& tgroup_x,
                                   const std::vector<double>& tgroup_y,
                                   double& score, double& statistic) {
  double z_value;
  double p_value = __uTest(tgroup_x, tgroup_y, z_value, statistic);
  if (alternative == 0) {
    p_value = std::min(p_value * 2., 1.0);
  } else {
    if (z_value < 0)
      p_value = 1. - p_value;
    if (alternative < 0)
      p_value = 1. - p_value;
  }
  calTime = calTime + 1;
  score = z_value;
  return p_value;
}

/**
 * Calculate p-value by using Mann-Whitney U test
 * @param tgroup_x test group 1. This is consisted of transaction values.
 * @param tgroup_y test group 2. This is consisted of transaction values.
 * @param z_value return value
 * @param u_value return value
 * @return
 */
double Functions4u_test::__uTest(const std::vector<double>& tgroup_x,
                                 const std::vector<double>& tgroup_y, double& z_value, double& u_value) {
  u_value = __uValue(tgroup_x, tgroup_y); // u-value of two groups.
  // z value of u-value
  //mean_u, var_u = __calStatValue(tgroup_x, tgroup_y);
  double size_x = tgroup_x.size();
  double size_y = tgroup_y.size();
  double mean_u = (size_x * size_y) / 2.0;
  //	double var_u = size_x * size_y * (size_x + size_y + 1.0) / 12.0;
  double var_u = mean_u * (size_x + size_y + 1.0) / 6.0;
  if (var_u == 0) {
    z_value = 0;
    return 1.0;
  }
  z_value = (u_value - mean_u) / std::sqrt(var_u);

  // calculate p-value from z_value
  // this value approximation of standard normal distribution
  return stdNorDistribution(z_value);
}

/**
 * Calculate u value which measurs difference rank sum of two groups.
 * @param tgroup_x test group 1. This is consisted of transaction values.
 * @param tgroup_y test group 2. This is consisted of transaction values.
 * tgroup_x and t_group_y already sorted by transaction value.
 * @return 
 */
double Functions4u_test::__uValue(const std::vector<double>& tgroup_x, const std::vector<double>& tgroup_y) {
  double u_value = 0.0;
  int previous_u_x_min = 0; // The rank of transaction in previous search
  int previous_u_x_max = 0; // The rank of transaction in previous search
  double previous_value = std::numeric_limits<double>::quiet_NaN(); // The previous expression value
  int left_index = 0; // The start point of searching value.
  int right_index = (int) tgroup_y.size() - 1; // The end point of searching value.
  for (std::size_t i = 0; i < tgroup_x.size(); i++) {
    double t_x = tgroup_x[i];
    // u_x_min: rank sum of transaction which the value < t_x in tgroup_y
    // u_x_max: rank sum of transaction which the value <= t_x in tgroup_y
    int u_x_min = std::numeric_limits<int>::quiet_NaN();
    int u_x_max = std::numeric_limits<int>::quiet_NaN();
    // If t_x.value is equal to previous one, u_x_min and u_x_max are also equals.
    if (t_x == previous_value) {
      u_x_min = previous_u_x_min;
      u_x_max = previous_u_x_max;
    }			// Caluclate u_value because tgroup_x value exists between tgroup_y range
    else {
      __binarySearch(t_x, tgroup_y, left_index, right_index, u_x_min, u_x_max);
      left_index = u_x_max;
    }
    // Add rank of t_x to u_value
    u_value = u_value + ((double) u_x_min + (double) u_x_max) / 2.0;
    previous_u_x_min = u_x_min;
    previous_u_x_max = u_x_max;
    previous_value = t_x;
  }
  return u_value;
}

/**
 * Search the group which t.value less than threshold.
 * @param threshold search value
 * @param tgroup list to search threshold
 * @param left_index start index of tgroup to search
 * @param right_index end index of tgroup to search
 * @param u_x_min return value
 * @param u_x_max return value
 */
void Functions4u_test::__binarySearch(double threshold, const std::vector<double>& tgroup,
                                      int left_index, int right_index, int& u_x_min, int& u_x_max) {
  if ((int) tgroup.size() <= left_index) {
    u_x_min = (int) tgroup.size();
    u_x_max = (int) tgroup.size();
    return;
  }

  // compare threshold to min and max value
  if (threshold < tgroup[left_index]) {
    u_x_min = left_index;
    u_x_max = left_index;
    return;
  }
  if (tgroup[right_index] < threshold) {
    u_x_min = right_index + 1;
    u_x_max = right_index + 1;
    return;
  }

  // serach the index which all larger indexes are more than threshold.
  int mid_index = -1;
  while (left_index <= right_index) {
    mid_index = (left_index + right_index) / 2;
    double mid_transaction = tgroup[mid_index];

    // When the mid.value = threshold, finish the search.
    if (mid_transaction == threshold)
      break;
    // When the check value less than threshod, go to serach right.
    else if (mid_transaction < threshold)
      left_index = mid_index + 1;
    // When the check value >= threshold, go to search left.
    else
      right_index = mid_index - 1;
  }

  // search the same range of the threshold
  double mid_transaction = tgroup[mid_index];
  if (mid_transaction == threshold) {
    int min_index = mid_index;
    int max_index = mid_index;
    double min_transaction = tgroup[min_index];
    double max_transaction = tgroup[max_index];
    while (threshold <= min_transaction) {
      min_index = min_index - 1;
      if (min_index < 0)
        break;
      min_transaction = tgroup[min_index];
    }
    while (max_transaction <= threshold) {
      max_index = max_index + 1;
      if ((int) tgroup.size() <= max_index)
        break;
      max_transaction = tgroup[max_index];
    }
    u_x_min = min_index + 1;
    u_x_max = max_index;
  }		// not found the threshold in tgroup.
  // in this case, min_index > max_index
  else if (mid_transaction < threshold) {
    while (mid_transaction < threshold) {
      mid_index = mid_index + 1;
      mid_transaction = tgroup[mid_index];
    }
    u_x_min = mid_index;
    u_x_max = mid_index;
  }		// not found the threshod in tgroup and the case of mid_transaction value > threshold
  // In this case, min_index > max_index
  else {
    while (threshold < mid_transaction) {
      mid_index = mid_index - 1;
      mid_transaction = tgroup[mid_index];
    }
    u_x_min = mid_index + 1;
    u_x_max = mid_index + 1;
  }
  return;
}

double Functions4u_test::calPValue(const int sup, const int pos_sup,
                                   const std::vector<double>& tgroup_x, const std::vector<double>& tgroup_y) {
  (void)sup; (void)pos_sup;

  double statistic;
  double z_value;
  double p_value = __uTest(tgroup_x, tgroup_y, z_value, statistic);
  if (alternative == 0) {
    p_value = std::min(p_value * 2., 1.0);
  } else {
    if (z_value < 0)
      p_value = 1. - p_value;
    if (alternative < 0)
      p_value = 1. - p_value;
  }
  calTime = calTime + 1;
  return p_value;
}
