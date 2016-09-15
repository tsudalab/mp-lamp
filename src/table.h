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

#ifndef _LAMP_SEARCH_TABLE_H_
#define _LAMP_SEARCH_TABLE_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <list>
#include <limits>

#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>

#include "utils.h"
#include "sorted_itemset.h"

namespace lamp_search {

class Table {
 public:
  Table(std::istream & item_file, std::istream & posneg_file);
  ~Table();

  void ReadItems(std::istream & is);
  void ReadPosNeg(std::istream & is);

  const std::vector< std::string > & ItemNames() const { return item_names_; }
  const std::string & NthItemName(std::size_t i) const { return ItemNames()[i]; }

  int NuItems() const { return nu_items_; }
  int NuTransaction() const { return nu_transactions_; }
  int PosTotal() const { return nu_pos_total_; }

  std::ostream & DumpItems(std::ostream & out) const;
  std::ostream & DumpPosNeg(std::ostream & out) const;
  std::ostream & ShowInfo(std::ostream & out) const;

  const boost::dynamic_bitset<> & NthData(std::size_t i) const {
    return (*data_)[i];
  }
  const boost::dynamic_bitset<> & PosNeg() const {
    return *posneg_;
  }

  //        pos   neg     freq
  //---------------------------
  // item |  t   (x-t)  |   x
  // rest | n-t         |  N-x
  //---------------------------
  //      |  n   (N-n)  |   N

  // lampelr local   my code                   lampelr global
  // size_all     == total        == N      == total_size_lst (const) == NuTransaction()
  // pos_size_all == pos_total    == n      == pos_size_lst   (const) == PosTotal()
  // support_all  == sup          == x      == group_sup
  // obs_t        == pos_sup      == t      == group_pos_sup

  double PMin(int sup) const { return pmin_table_[sup]; }
  double PVal(int sup, int pos_sup) const { return pval_table_[sup * (max_t_+1) + pos_sup];}

  // todo: prepare confound factor version
  // double PVal2(int sup, int pos_sup) const; // temp
  double PMinCal(int sup) const;
  double PMinCalSub(int sup) const;
  double PValCal(int sup, int pos_sup) const;

  void InitPMinTable();
  void InitPValTable();

  std::ostream & DumpPMinTable(std::ostream & out) const;
  std::ostream & DumpPValTable(std::ostream & out) const;

  int MaxX() const { return max_x_; }
  int MaxT() const { return max_t_; }
  int MaxItemInTransaction() const { return max_item_in_transaction_; }

  // for test
  void SetSigLev(double d) { siglev = d; }

  void PrepareItemVals();

  struct PValCache_t {
    long long int n; // number of usage
    double v; // value
  };

  const std::vector<ItemInfo> & GetItemInfo() const { return item_info_; }

  void SetValuesForTest(int nu_item, int nu_transaction, int nu_pos_total);

 private:
  int nu_items_;
  std::vector< std::string > item_names_;

  int nu_transactions_;
  std::vector< std::string > transaction_names_;
  int nu_pos_total_;
  
  //std::vector< boost::dynamic_bitset<> > v(10, boost::dynamic_bitset<>(3));
  std::vector< boost::dynamic_bitset<> > * data_;
  // std::vector< boost::dynamic_bitset<> > ** data_; // for confsize > 1

  boost::dynamic_bitset<> * posneg_;
  int max_x_;
  int max_t_;
  int max_item_in_transaction_;

  // ----
  // these are following lampeler variable naming. no trailing _. be careful
  double siglev;

  // reused during PVal calculation 
  double * pval_cal_buf; // originally pval_table
  // ----

  // stores calculated pmin value
  std::vector<double> pmin_table_;
  std::vector<double> pval_table_;

  // curerntly not used
  std::vector<PValCache_t> pval_cache_;
  std::vector<PValCache_t> pmin_cache_;

  std::vector<ItemInfo> item_info_; // list of item id sorted by pmin (sup)
};

} // namespace lamp_search

#endif // _LAMP_SEARCH_TABLE_H_
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
