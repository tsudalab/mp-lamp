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

#ifndef _LAMP_SEARCH_DATABASE_H_
#define _LAMP_SEARCH_DATABASE_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <list>
#include <limits>
#include <set>
#include <map>

#include <boost/array.hpp>

#include "utils.h"
#include "sorted_itemset.h"
#include "variable_length_itemset.h"
#include "variable_bitset_array.h"
#include "functions/FunctionsSuper.h"

namespace lamp_search {

// database with database reduction
// todo: implement based on table_vba

class DatabaseParameter {
 public:

  DatabaseParameter() :
      nu_trans_(-1),
      nu_items_(-1),
      nu_non_zero_trans_(-1),
      nu_non_zero_items_(-1),
      nu_pos_total_(-1),
      max_x_(-1),
      max_t_(-1),
      max_item_in_transaction_(-1),
      num_1_(-1ll) {
  }

  int nu_trans_;
  int nu_items_;
  int nu_non_zero_trans_;
  int nu_non_zero_items_;
  int nu_pos_total_;
  int max_x_;
  int max_t_;
  int max_item_in_transaction_;
  long long int num_1_;
};

template<typename Block>
class DatabaseReader {
 public:

  int CountLines(std::istream & is) const;

  bool ReadFirstPhase(std::istream & is,
                      int * nu_trans,
                      std::vector< std::string > * transaction_names,
                      int * nu_non_zero_trans);

  std::vector<int> non_zero_trans_list_;
  std::map< std::string, int > item_name_id_map_;

  // given stream, allocate Block array and read values,
  // return nu_items, bitset helper and pointer to allocated Block array
  void ReadItems(std::istream & is,
                 int * nu_trans,
                 int * nu_items,
                 VariableBitsetHelper<Block> ** bsh, Block ** data,
                 std::vector< std::string > * item_names,
                 std::vector< std::string > * transaction_names,
                 int * max_item_in_transaction);
  // given stream, allocate Block array and read values
  void ReadPosNeg(std::istream & is, int nu_trans,
                  std::vector< std::string > * trans_names,
                  int * nu_pos_total,
                  const VariableBitsetHelper<Block> * bsh, Block ** positive, double ** pos_val,
                  bool reverse);

  bool ReadFirstPhaseLCM(std::istream & is,
                         int * nu_trans,
                         int * nu_items,
                         int * nu_non_zero_trans);
  // read item for LCM format
  void ReadItemsLCM(std::istream & is,
                    int * nu_trans,
                    int * nu_items,
                    VariableBitsetHelper<Block> ** bsh, Block ** data,
                    std::vector< std::string > * item_names,
                    int * max_item_in_transaction);

  // if item name is NULL, use numbers
  std::ostream & PrintLCM(std::ostream & out,
                          int nu_trans,
                          int nu_items,
                          const VariableBitsetHelper<Block> * bsh,
                          const Block * data) const;
  std::ostream & PrintCSV(std::ostream & out,
                          int nu_trans,
                          int nu_items,
                          const VariableBitsetHelper<Block> * bsh,
                          std::vector< std::string > * item_names,
                          const Block * data) const;

  void ReadFiles(VariableBitsetHelper<Block> ** bsh,
                 std::istream & item_file,
                 Block ** item,
                 int * nu_trans,
                 int * nu_items,
                 std::istream & pos_file,
                 Block ** positive, double ** pos_val,
                 int * nu_pos_total,
                 std::vector< std::string > * item_names,
                 std::vector< std::string > * transaction_names,
                 int * max_item_in_transaction,
                 bool reverse);

  void ReadFiles(VariableBitsetHelper<Block> ** bsh,
                 std::istream & item_file,
                 Block ** item,
                 int * nu_trans,
                 int * nu_items,
                 std::vector< std::string > * item_names,
                 std::vector< std::string > * transaction_names,
                 int * max_item_in_transaction);

  void ReadFilesLCM(VariableBitsetHelper<Block> ** bsh,
                    std::istream & item_file,
                    Block ** item,
                    int * nu_trans,
                    int * nu_items,
                    std::istream & pos_file,
                    Block ** positive, double ** pos_val,
                    int * nu_pos_total,
                    std::vector< std::string > * item_names,
                    int * max_item_in_transaction,
                    bool reverse);

  void ReadFilesLCM(VariableBitsetHelper<Block> ** bsh,
                    std::istream & item_file,
                    Block ** item,
                    int * nu_trans,
                    int * nu_items,
                    std::vector< std::string > * item_names,
                    int * max_item_in_transaction);

};

template<typename Block>
class Database {
 public:
  typedef VariableBitsetHelper<Block> VBH;

  Database(VariableBitsetHelper<Block> * bsh,
           Block * item_array,
           std::size_t nu_trans,
           std::size_t nu_items,
           Block * pos_array, double * pos_val, std::size_t nu_pos_total,
           int max_item_in_transaction,
           std::vector< std::string > * item_names,
           std::vector< std::string > * trans_names,
           FunctionsSuper & functions);

  ~Database();

  void Init();

  const std::vector< std::string > * ItemNames() const {
    return item_names_;
  }

  const std::string & NthItemName(std::size_t i) const {
    return (*ItemNames())[i];
  }

  const std::vector< std::string > * TransactionNames() const {
    return transaction_names_;
  }

  const std::string & NthTransactionNames(std::size_t i) const {
    return (*TransactionNames())[i];
  }

  int NuItems() const {
    return nu_items_;
  }

  int NuTransaction() const {
    return nu_transactions_;
  }

  int PosTotal() const {
    return nu_pos_total_;
  }

  int NuReducedItems(int lambda) const;
  int FirstReducedItemInReverseLoop(int lambda, bool is_root_node,
                                    int rank, int nu_proc) const;
  int NextReducedItemInReverseLoop(int lambda, bool is_root_node,
                                   int rank, int nu_proc, int prev) const;

  // return largest int j < i such that j%nu_proc == rank
  // for first item in loop, set i == NuItems()

  static int NextItemInReverseLoop(bool is_root_node, int rank, int nu_proc, int prev) {
    // return largest int i < prev such that i%nu_proc == rank
    if (is_root_node) {
      int i = prev - prev % nu_proc + rank;
      if (i >= prev) i -= nu_proc;
      return i;
    } else
      return prev - 1;
  }

  std::ostream & DumpItems(std::ostream & out) const;
  std::ostream & DumpPosNeg(std::ostream & out) const;
  std::ostream & ShowInfo(std::ostream & out) const;

  const VBH & VBSHelper() const {
    return *bsh_;
  }

  Block * NthData(std::size_t i) {
    return bsh_->N(data_, i);
  }

  const Block * NthData(std::size_t i) const {
    return bsh_->N(data_, i);
  }

  Block * PosNeg() const {
    return posneg_;
  }

  double * PosVal() const {
    return pos_val_;
  }

  Block * Data() {
    return data_;
  }

  const Block * Data() const {
    return data_;
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

  double PMin(int sup) const {
    return functions.funcF(sup);
  }

  double PMinLog(int sup) const {
    return log(PMin(sup));
  }

  double PVal(int sup, int pos_sup, Block * dt, double * pos_val) const {
    //	  return pval_table_[sup * (max_t_+1) + pos_sup];
    double p = 1.01;
    try {
      std::vector<double> tgroup_x;
      std::vector<double> tgroup_y;
      if (!functions.isBinary()) {
        std::size_t ptr = 0;
        for (std::size_t i = 0; i < bsh_->NuBlocks(); i++) {
          Block dtVal = dt[i];
          Block mask = 1u;
          for (int j = 0; j < VariableBitsetTraits<Block>::bits_per_block; j++) {
            if ((dtVal & mask) > 0)
              tgroup_x.push_back(pos_val[ptr]);
            else
              tgroup_y.push_back(pos_val[ptr]);
            ptr++;
            if (ptr == bsh_->NuBits())
              break;
            mask <<= 1u;
          }
        }
        for (int i = bsh_->NuBits(); i < nu_transactions_; i++) {
          tgroup_y.push_back(pos_val[i]);
        }
        std::sort(tgroup_x.begin(), tgroup_x.end());
        std::sort(tgroup_y.begin(), tgroup_y.end());
      }
      p = functions.calPValue(sup, pos_sup, tgroup_x, tgroup_y);
    } catch (...) {
    }
    return p;
  }

  std::ostream & DumpPMinTable(std::ostream & out) const;
  std::ostream & DumpPValTable(std::ostream & out) const;

  int MaxX() const {
    return max_x_;
  }

  int MaxT() const {
    return max_t_;
  }

  int MaxItemInTransaction() const {
    return max_item_in_transaction_;
  }

  long long int Count1() const;
  double Density() const;

  // for test
  // void SetSigLev(double d) { siglev = d; }

  void PrepareItemVals();

  const std::vector<ItemInfo> & GetItemInfo() const {
    return item_info_;
  }

  int PrepareReducedList(int lambda);
  std::vector<int> & ReducedList(int lambda);

  const std::vector<int> & SupHist() const {
    return sup_hist_;
  }

  const std::vector<int> & SupCumHist() const {
    return sup_cum_hist_;
  }

  int NuAllZeroTrans() const;
  int NuAllZeroItems() const;

  const bool isReverse() const {
    return functions.isReverse();
  }

  void SetValuesForTest(int nu_item, int nu_transaction, int nu_pos_total);

 private:
  VBH * bsh_; // biset helper
  // VBH ** bsh_; // biset helper for each lambda

  int nu_items_;
  std::vector< std::string > * item_names_;

  int nu_transactions_;
  std::vector< std::string > * transaction_names_;
  int nu_pos_total_;

  // todo: prepare multiple data_ for confounder

  Block * data_;
  int max_x_;

  bool has_positives_;
  Block * posneg_;
  double * pos_val_;
  int max_t_;

  int max_item_in_transaction_;

  // ----
  // these are following lampeler variable naming. no trailing _. be careful
  // double siglev;

  // support histogram
  std::vector<int> sup_hist_;
  std::vector<int> sup_cum_hist_;
  std::vector< std::vector<int> > reduced_item_list_;
  std::vector<bool> reduced_item_list_prepared_;

  std::vector<ItemInfo> item_info_; // list of item id sorted by pmin (sup)

  FunctionsSuper & functions; // the function to calculate P-value
};

} // namespace lamp_search

#endif // _LAMP_SEARCH_DATABASE_H_

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
