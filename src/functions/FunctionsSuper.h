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
 * File:   FunctionsSuper.h
 * Author: miura
 *
 * Created on 2016/01/06, 15:28
 */

#ifndef FUNCTIONSSUPER_H
#define FUNCTIONSSUPER_H

#include <cmath>
#include <string>
#include <vector>
#include <map>

/** @file */

/** this class is base to calculate P-value
 */
class FunctionsSuper {
 public:
  FunctionsSuper(int alternative);
  FunctionsSuper(int transaction_size, int n1_count, int alternative);
  virtual ~FunctionsSuper();

  /**
   * get N1
   * @return N1
   */
  int getN1() {
    return __f_size;
  }

  /**
   * set N1
   * @param n1_count transaction size which have flag = 1 (n1)
   */
  void setN1(int n1_count) {
    __f_size = n1_count;
  }

  /**
   * get size of all.
   * @return size of all
   */
  int getAllSize() {
    return __t_size;
  }

  /**
   * set size of all.
   * @param transaction_size all transaction size
   */
  void setAllSize(int transaction_size) {
    __t_size = transaction_size;
  }

  double stdNorDistribution(double x);

  /**
   * This function calculates the minimum p-value which support size is x.
   * @param x support size
   * @return  
   */
  virtual double funcF(int x) = 0;

  /**
   * calculate p-value
   * @param ovalues
   * @param score
   * @return 
   */
  virtual double calPValue(const int sup, const int pos_sup,
                           const std::vector<double>& tgroup_x, const std::vector<double>& tgroup_y) {
    return -1;
  };

  /**
   * get total number of calculate P-value.
   * @return total number of calculate P-value
   **/
  int getCalTime() {
    return calTime;
  }

  /**
   * set total number of calculate P-value to zero.
   **/
  void resetCalTime() {
    calTime = 0;
  }

  /**
   * get name of statistical hypothesis test
   * @return name of method
   */
  virtual const std::string getName() = 0;

  /**
   * get alternative setting
   * @return 
   */
  int getAlternative() {
    return alternative;
  }

  /**
   * Reverse the observed values for alternative = 'less'.
   * @param val observed value.
   */
  virtual double reverseValue(double val) = 0;

  /**
   * get reverse condition.
   **/
  bool isReverse() {
    return reverse;
  }

  /**
   * get binary method.
   **/
  virtual bool isBinary() = 0;

 protected:
  int __f_size; /**< transaction size which have flag = 1 (n1) */
  int __t_size; /**< all transaction size */
  int alternative; /**< alternative hypothesis. greater or less -> 1, two.sided -> 0. */
  int calTime; /**< total number of calculate P-value */
  bool reverse; /**< alternative hypothesis. less -> true, other -> false. */
  std::map<int, double> __minpvalTable; /**< the minimum p-value which support size is x */

 private:
  const static double pi2; /**< constant value */

};

#endif /* FUNCTIONSSUPER_H */
