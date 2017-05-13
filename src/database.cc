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

#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <cmath>
#include <set>
#include <map>

#include <boost/array.hpp>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
//#include <boost/range/algorithm/for_each.hpp>

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#include "gflags/gflags.h"

#include "utils.h"
#include "sorted_itemset.h"
#include "variable_bitset_array.h"
#include "database.h"

// DECLARE_double(a); // significance level alpha
// DEFINE_int32(posnum, 0, "positive total (used if for 1st phase only)");

//==============================================================================

namespace lamp_search {

template<typename Block>
void DatabaseReader<Block>::ReadFiles(VariableBitsetHelper<Block> ** bsh,
                                      std::istream & item_file,
                                      Block ** item,
                                      int * nu_trans,
                                      int * nu_items,
                                      std::istream & pos_file,
                                      Block ** positive,
                                      int * nu_pos_total,
                                      std::vector< std::string > * item_names,
                                      std::vector< std::string > * transaction_names,
                                      int * max_item_in_transaction)
{
  ReadItems(item_file, nu_trans, nu_items, bsh,
            item, item_names, transaction_names, max_item_in_transaction);
  ReadPosNeg(pos_file, *nu_trans, transaction_names, nu_pos_total, *bsh, positive);
}

template<typename Block>
void DatabaseReader<Block>::ReadFiles(VariableBitsetHelper<Block> ** bsh,
                                      std::istream & item_file,
                                      Block ** item,
                                      int * nu_trans,
                                      int * nu_items,
                                      std::vector< std::string > * item_names,
                                      std::vector< std::string > * transaction_names,
                                      int * max_item_in_transaction)
{
  ReadItems(item_file, nu_trans, nu_items, bsh,
            item, item_names, transaction_names, max_item_in_transaction);
}

template<typename Block>
void DatabaseReader<Block>::ReadFilesLCM(VariableBitsetHelper<Block> ** bsh,
                                         std::istream & item_file,
                                         Block ** item,
                                         int * nu_trans,
                                         int * nu_items,
                                         std::istream & pos_file,
                                         Block ** positive,
                                         int * nu_pos_total,
                                         std::vector< std::string > * item_names,
                                         int * max_item_in_transaction)
{
  ReadItemsLCM(item_file, nu_trans, nu_items, bsh,
               item, item_names, max_item_in_transaction);
  ReadPosNeg(pos_file, *nu_trans, NULL, nu_pos_total, *bsh, positive);
}

template<typename Block>
void DatabaseReader<Block>::ReadFilesLCM(VariableBitsetHelper<Block> ** bsh,
                                         std::istream & item_file,
                                         Block ** item,
                                         int * nu_trans,
                                         int * nu_items,
                                         std::vector< std::string > * item_names,
                                         int * max_item_in_transaction)
{
  ReadItemsLCM(item_file, nu_trans, nu_items, bsh,
               item, item_names, max_item_in_transaction);
}

template<typename Block>
int DatabaseReader<Block>::CountLines(std::istream & is) const {
  int nu_lines = 0;
  std::string line;
  while (true) {
    std::getline(is, line);
    if ( ! is.good() ) {
      is.clear();
      is.seekg(0); // rewind
      break;
    }
    nu_lines++;
  }
  return nu_lines;
}

template<typename Block>
bool DatabaseReader<Block>::ReadFirstPhase(std::istream & is,
                                           int * nu_trans,
                                           std::vector< std::string > * transaction_names,
                                           int * nu_non_zero_trans) {
  std::string line;
  std::string trimmed_line;
  boost::char_separator<char> sep(", ");

  {
    std::getline(is, line);
    if ( ! is.good() ) {
      *nu_trans = 0;
      *nu_non_zero_trans = 0;
      is.clear();
      is.seekg(0, is.beg); // rewind
      return false;
    }
  }

  *nu_trans = 0;
  *nu_non_zero_trans = 0;
  while (1) {
    std::getline(is, line);
    trimmed_line = boost::algorithm::trim_copy(line);
    // eof, fail, bad
    if ( ! is.good() ) break;

    boost::tokenizer< boost::char_separator<char> > tokens(trimmed_line, sep);
    boost::tokenizer< boost::char_separator<char> >::iterator tok_iter=tokens.begin();
    transaction_names->push_back(*tok_iter); // first element is transaction name
    ++ tok_iter;
    int counter = 0; // item counter for this transaction
    int pos_counter = 0; // pos counter for this transaction
    for (; tok_iter != tokens.end(); ++tok_iter) {
      if (*tok_iter != "0") pos_counter++;
      counter++;
    }
    if (pos_counter > 0) {
      non_zero_trans_list_.push_back(*nu_trans);
      (*nu_non_zero_trans)++;
    }
    (*nu_trans)++;
  }

  assert(non_zero_trans_list_.size() == (std::size_t)(*nu_non_zero_trans));

  is.clear();
  is.seekg(0, is.beg); // rewind
  return true;
}

template<typename Block>
void DatabaseReader<Block>::ReadItems(std::istream & is,
                                      int * nu_trans,
                                      int * nu_items,
                                      VariableBitsetHelper<Block> ** bsh,
                                      Block ** data,
                                      std::vector< std::string > * item_names,
                                      std::vector< std::string > * transaction_names,
                                      int * max_item_in_transaction) {
  *nu_items = 0;
  int nu_transactions = 0;
  int nu_non_zero_transactions = 0;
  
  // read 1st line, tokenize, count items
  // read line, tokenize

  std::string line;
  std::string trimmed_line;
  boost::char_separator<char> sep(", ");

  bool ok = ReadFirstPhase(is, &nu_transactions, transaction_names,
                           &nu_non_zero_transactions);
  if (!ok || nu_transactions == 0) {
    bsh = NULL;
    data = NULL;
    *nu_trans = nu_transactions;
    *nu_items = 0;
    *max_item_in_transaction = 0;
  }
  transaction_names->reserve(nu_non_zero_transactions);
  *max_item_in_transaction = -1;
  *nu_trans = nu_transactions;

  {
    std::getline(is, line);
    trimmed_line = boost::algorithm::trim_copy(line);
    if ( ! is.good() ) {
      bsh = NULL;
      data = NULL;
      *nu_trans = nu_transactions;
      *nu_items = 0;
      *max_item_in_transaction = 0;
    }
    boost::tokenizer< boost::char_separator<char> > tokens(trimmed_line, sep);
    
    boost::tokenizer< boost::char_separator<char> >::iterator tok_iter=tokens.begin();
    ++ tok_iter; // skip first element. should be "#gene"
    int counter = 0;
    for (; tok_iter != tokens.end(); ++tok_iter) {
      item_names->push_back(*tok_iter);
      counter++;
    }
    *nu_items = counter;
  }

  *bsh = new VariableBitsetHelper<Block>(nu_non_zero_transactions);
  // bitset array for num bits=nu_transactions and entries=nu_items
  *data = (*bsh)->NewArray( (std::size_t)(*nu_items));

  int trans_counter = 0;
  int non_zero_trans_counter = 0;
  while (1) {
    std::getline(is, line);
    trimmed_line = boost::algorithm::trim_copy(line);
    // eof, fail, bad
    if ( ! is.good() || non_zero_trans_counter == nu_non_zero_transactions ) break;

    boost::tokenizer< boost::char_separator<char> > tokens(trimmed_line, sep);
    boost::tokenizer< boost::char_separator<char> >::iterator tok_iter=tokens.begin();
    // read transaction name in 1st phase, skipping
    ++ tok_iter;
    int counter = 0; // item counter for this transaction
    int pos_counter = 0; // pos counter for this transaction
    for (; tok_iter != tokens.end(); ++tok_iter) {
      if (*tok_iter == "0") {
        (*bsh)->Reset(non_zero_trans_counter, (*bsh)->N((*data), counter)); // (pos, *elm)
      }
      else {
        (*bsh)->Doset(non_zero_trans_counter, (*bsh)->N((*data), counter)); // (pos, *elm)
        pos_counter++;
      }
      counter++;
    }
    *max_item_in_transaction = std::max(*max_item_in_transaction, pos_counter);
    trans_counter++;
    if (pos_counter > 0) non_zero_trans_counter++;

    assert( *nu_items == counter );
  }
  // assert( nu_transactions == trans_counter ); // does not hold if last lines are skipped
  assert( nu_non_zero_transactions == non_zero_trans_counter );
  if (nu_non_zero_transactions != non_zero_trans_counter)
    throw std::runtime_error("read 1st 2nd phase mismatch (shouldn't happen)");
}

template<typename Block>
bool DatabaseReader<Block>::ReadFirstPhaseLCM(std::istream & is,
                                              int * nu_trans,
                                              int * nu_items,
                                              int * nu_non_zero_trans) {
  std::string line;
  std::string trimmed_line;
  boost::char_separator<char> sep(" \t");

  std::set< std::string > item_names;

  *nu_trans = 0;
  *nu_non_zero_trans = 0;
  while (1) {
    std::getline(is, line);
    trimmed_line = boost::algorithm::trim_copy(line);
    // eof, fail, bad
    if ( ! is.good() ) break;

    boost::tokenizer< boost::char_separator<char> > tokens(trimmed_line, sep);
    boost::tokenizer< boost::char_separator<char> >::iterator tok_iter=tokens.begin();

    int pos_counter = 0; // pos counter for this transaction
    for (; tok_iter != tokens.end(); ++tok_iter) {
      std::istringstream in(*tok_iter);
      std::string name;
      in >> name;
      if (!in) std::cout << "file format error? token: " << *tok_iter << std::endl;

      item_names.insert(name);

      pos_counter++;
      if (*tok_iter != "0") pos_counter++;
    }
    if (pos_counter > 0) {
      non_zero_trans_list_.push_back(*nu_trans);
      (*nu_non_zero_trans)++;
    }
    (*nu_trans)++;
  }
  *nu_items = item_names.size();
  assert(non_zero_trans_list_.size() == (std::size_t)(*nu_non_zero_trans));

  is.clear();
  is.seekg(0, is.beg); // rewind

  int item_count = 0;
  for (std::set< std::string >::const_iterator it = item_names.begin();
       it != item_names.end(); ++it) {
    item_name_id_map_.insert( std::pair<std::string, int>(*it, item_count) );
    item_count++;
  }

  return true;
}

template<typename Block>
void DatabaseReader<Block>::ReadItemsLCM(std::istream & is,
                                         int * nu_trans,
                                         int * nu_items,
                                         VariableBitsetHelper<Block> ** bsh,
                                         Block ** data,
                                         std::vector< std::string > * item_names,
                                         int * max_item_in_transaction) {
  *nu_items = 0;
  int nu_non_zero_transactions = 0;

  std::string line;
  std::string trimmed_line;
  boost::char_separator<char> sep(" \t");

  // read line, tokenize
  // (istream, nu_trans, nu_items, nu_non_zero_transactions)
  bool ok = ReadFirstPhaseLCM(is, nu_trans, nu_items, &nu_non_zero_transactions);

  if (!ok || *nu_trans == 0) {
    bsh = NULL;
    data = NULL;
    *max_item_in_transaction = 0;
  }
  *max_item_in_transaction = -1;

  *bsh = new VariableBitsetHelper<Block>(nu_non_zero_transactions);
  // bitset array for num bits=nu_non_zero_transactions and entries=nu_items
  *data = (*bsh)->NewArray( (std::size_t)(*nu_items));

  for(std::map< std::string, int >::const_iterator it = item_name_id_map_.begin();
      it != item_name_id_map_.end(); ++it) {
    item_names->push_back( (*it).first );
  }

  int trans_counter = 0;
  int non_zero_trans_counter = 0;
  while (1) {
    std::getline(is, line);

    trimmed_line = boost::algorithm::trim_copy(line);
    // eof, fail, bad
    if ( ! is.good() || non_zero_trans_counter == nu_non_zero_transactions ) break;

    boost::tokenizer< boost::char_separator<char> > tokens(trimmed_line, sep);
    boost::tokenizer< boost::char_separator<char> >::iterator tok_iter=tokens.begin();
    int pos_counter = 0; // pos counter for this transaction

    for (; tok_iter != tokens.end(); ++tok_iter) {
      std::istringstream in(*tok_iter);
      std::string name;
      in >> name;

      std::map< std::string, int >::const_iterator it = item_name_id_map_.find(name);
      assert(it != item_name_id_map_.end());
      int item_id = (*it).second;

      assert(item_id >= 0); assert(item_id < *nu_items);
      (*bsh)->Doset(non_zero_trans_counter, (*bsh)->N((*data), item_id)); // (pos, *elm)
      pos_counter++;
    }
    *max_item_in_transaction = std::max(*max_item_in_transaction, pos_counter);
    trans_counter++;
    if (pos_counter > 0) non_zero_trans_counter++;
  }

  assert( nu_non_zero_transactions == non_zero_trans_counter );
  if (nu_non_zero_transactions != non_zero_trans_counter)
    throw std::runtime_error("read 1st 2nd phase mismatch (shouldn't happen)");
}

template<typename Block>
void DatabaseReader<Block>::ReadPosNeg(std::istream & is,
                                       int nu_trans,
                                       std::vector< std::string > * trans_names,
                                       int * nu_pos_total,
                                       const VariableBitsetHelper<Block> * bsh,
                                       Block ** positive) {
  std::string line;
  std::string trimmed_line;
  boost::char_separator<char> sep(", ");

  {
    std::getline(is, line); // skip 1st line
  }

  *positive = bsh->New();

  int trans_counter = 0;
  int non_zero_trans_counter = 0;
  (*nu_pos_total) = 0;
  while (1) {
    std::getline(is, line);
    trimmed_line = boost::algorithm::trim_copy(line);
    // eof, fail, bad
    if ( ! is.good()  ) break;

    boost::tokenizer< boost::char_separator<char> > tokens(trimmed_line, sep);
    boost::tokenizer< boost::char_separator<char> >::iterator tok_iter=tokens.begin();
    if ( trans_names != NULL && (*trans_names)[trans_counter] != (*tok_iter) )
      throw std::runtime_error(
          std::string("item file / positive file trans name mismatch ")
          + (*trans_names)[trans_counter]
          + std::string(" : ")
          + (*tok_iter));
    ++ tok_iter;

    if (non_zero_trans_list_[non_zero_trans_counter] == trans_counter) {
      if (*tok_iter == "0") {
        bsh->Reset(non_zero_trans_counter, *positive);
      }
      else {// "1"
        bsh->Doset(non_zero_trans_counter, *positive);
        (*nu_pos_total)++;
      }
      non_zero_trans_counter++;
    } else {
      if (*tok_iter != "0") {// "1"
        (*nu_pos_total)++;
      }
    }
    trans_counter++;
  }

  if (non_zero_trans_list_.size() != (std::size_t)non_zero_trans_counter)
    throw std::runtime_error("item file / positive file non zero trans mismatch");
  if (nu_trans != trans_counter)
    throw std::runtime_error("item file / positive file trans mismatch");

}

template<typename Block>
std::ostream & DatabaseReader<Block>::PrintLCM(std::ostream & out,
                                               int nu_trans,
                                               int nu_items,
                                               const VariableBitsetHelper<Block> * bsh,
                                               const Block * data) const {
  std::stringstream s;

  int non_zero_index=0;
  for (int ti=0;ti<nu_trans;ti++) {
    if (non_zero_trans_list_[non_zero_index] == ti) {
      for (int ii=0;ii<nu_items;ii++) {
        if ( bsh->Test(bsh->N(data, ii), non_zero_index) )
          s << (ii+1) << " ";
      }
      s << std::endl;
      non_zero_index++;
    } else {
      // add blank line for all zero transactions
      s << std::endl;
    }
  }

  assert(non_zero_index = bsh->nu_bits);

  out << s.str() << std::flush;
  return out;
}

template<typename Block>
std::ostream & DatabaseReader<Block>::PrintCSV(std::ostream & out,
                                               int nu_trans,
                                               int nu_items,
                                               const VariableBitsetHelper<Block> * bsh,
                                               std::vector< std::string > * item_names,
                                               const Block * data) const {
  std::stringstream s;

  // header line
  s << "#i";
  if (item_names == NULL) for (int ii=0;ii<nu_items;ii++) s << "," << ii;
  else {
    for (std::map< std::string, int >::const_iterator it = item_name_id_map_.begin();
         it != item_name_id_map_.end(); ++it)
    {
      s << "," << (*it).first;
    }
  }
  s << std::endl;

  int non_zero_index=0;
  for (int ti=0;ti<nu_trans;ti++) {
    s << ti;
    if (non_zero_trans_list_[non_zero_index] == ti) {
      for (int ii=0;ii<nu_items;ii++) {
        if ( bsh->Test(bsh->N(data, ii), non_zero_index) ) s << ",1";
        else s << ",0";
      }
      s << std::endl;
      non_zero_index++;
    } else {
      // add all zero line for all zero transactions
      for (int ii=0;ii<nu_items;ii++) s << ",0";
      s << std::endl;
    }
  }

  assert(non_zero_index = bsh->nu_bits);

  out << s.str() << std::flush;
  return out;
}

//==============================================================================

template<typename Block>
Database<Block>::Database(VariableBitsetHelper<Block> * bsh,
                          Block * item_array,
                          std::size_t nu_trans,
                          std::size_t nu_items,
                          Block * pos_array, std::size_t nu_pos_total,
                          int max_item_in_transaction,
                          std::vector< std::string > * item_names,
                          std::vector< std::string > * trans_names) :
    bsh_ (bsh),
    nu_items_ (nu_items),
    item_names_ (item_names),
    nu_transactions_ (nu_trans),
    transaction_names_ (trans_names),
    nu_pos_total_(nu_pos_total),
    data_ (item_array),
    max_x_ (-1),
    has_positives_ ( !(pos_array == NULL) ),
    posneg_ (pos_array),
    max_t_ (-1),
    max_item_in_transaction_ (max_item_in_transaction)
{
  assert(nu_pos_total > 0);
  Init();
}

template<typename Block>
Database<Block>::~Database() {
  if (data_)   bsh_->Delete(data_);
  if (posneg_) bsh_->Delete(posneg_);

  if (item_names_) delete item_names_;
  if (transaction_names_) delete transaction_names_;
}

template<typename Block>
void Database<Block>::Init() {
  // SetSigLev(FLAGS_a); // double
  pval_cal_buf = new double[NuTransaction()];
  pval_log_cal_buf = new double[NuTransaction()];
  PrepareItemVals();
}

template<typename Block>
void Database<Block>::SetValuesForTest(int nu_item, int nu_transaction, int nu_pos_total) {
  nu_items_ = nu_item;
  nu_transactions_ = nu_transaction;
  nu_pos_total_ = nu_pos_total;
  pval_cal_buf = new double[NuTransaction()];
  pval_log_cal_buf = new double[NuTransaction()];

  // should be given as paramters?
  max_x_ = nu_transaction;
  max_t_ = nu_pos_total;
}

template<typename Block>
void Database<Block>::InitPMinLogTable() {
  pmin_table_.resize(max_x_+1);
  pmin_log_table_.resize(max_x_+1);
  for (int i=0;i<=max_x_;i++) {
    pmin_log_table_[i] = PMinCalLog(i);
    pmin_table_[i] = exp(pmin_log_table_[i]);
  }
}

template<typename Block>
std::ostream & Database<Block>::DumpPMinTable(std::ostream & out) const {
  std::stringstream s;

  for (int i=0;i<=max_x_;i++)
    s << "i=" << i
      << "\tpmin=" << pmin_table_[i]
      << "\tpmin_log=" << pmin_log_table_[i] << std::endl;

  out << s.str() << std::flush;
  return out;
}

template<typename Block>
void Database<Block>::InitPValTableLog() {
  pval_table_.resize( (max_x_+1) * (max_t_+1) );
  for (int i=0;i<=max_x_;i++) { // sup == x
    for (int j=0;j<=max_t_;j++) { // pos_sup == t
      if (j>i) continue;
      pval_table_[i*(max_t_+1) + j] = PValCalLog(i, j);
    }
  }
}

template<typename Block>
std::ostream & Database<Block>::DumpPValTable(std::ostream & out) const {
  std::stringstream s;

  for (int i=0;i<=max_x_;i++) { // sup == x
    for (int j=0;j<=max_t_;j++) { // pos_sup == t
      if (j>i) continue;
      s << "i=" << i
        << "\tj=" << j
        << "\tij" << i * j
        << "\tpval=" << pval_table_[i*(max_t_+1) + j]
        << "\tpval=" << PVal(i, j)
        << std::endl;
    }
  }

  out << s.str() << std::flush;
  return out;
}

template<typename Block>
std::ostream & Database<Block>::DumpItems(std::ostream & out) const {
  std::stringstream s;

  for (int i=0 ; i<NuItems() ; i++) {
    bsh_->Print(s, bsh_->N(data_, i));
    s << "\n";
  }

  out << s.str() << std::flush;
  return out;
}

template<typename Block>
std::ostream & Database<Block>::DumpPosNeg(std::ostream & out) const {
  std::stringstream s;
  bsh_->Print(s, posneg_);
  s << std::endl;

  out << s.str() << std::flush;
  return out;
}

template<typename Block>
std::ostream & Database<Block>::ShowInfo(std::ostream & out) const {
  std::stringstream s;

  s << "# of transactions=" << std::setw(12) << nu_transactions_
    << "\t# of items="      << std::setw(12) << nu_items_
    << "\t# of total positives="  << std::setw(12) << nu_pos_total_
    << "\tmax freq="         << std::setw(12) << max_x_
    << "\tmax positive="         << std::setw(12) << max_t_
    << "\tmax items in trans.=" << std::setw(12) << max_item_in_transaction_
    << std::endl;
  
  out << s.str() << std::flush;
  return out;
}

//==============================================================================

template<typename Block>
double Database<Block>::PMinCalLog(int sup) const {
  //if (sup == 0) return -(std::numeric_limits<double>::infinity());
  double minp_log = 0.0;

  { // for (int i = 0; i < confsize; ++i
    minp_log += PMinCalLogSub(sup);
  }
  return minp_log;
}

template<typename Block>
double Database<Block>::PMinCalLogSub(int sup) const {
  double minp_i_log = 0.0;

  unsigned int uplim = sup;
  if (sup > PosTotal()) uplim = PosTotal();
  for (unsigned int i = 0; i < uplim; ++i) {
    minp_i_log += log((double)(PosTotal() - i));
    minp_i_log -= log((double)(NuTransaction() - i));
  }
  assert(!std::isnan(minp_i_log));
  return minp_i_log;
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

template<typename Block>
double Database<Block>::PValCalLog(int sup, int pos_sup) const {
  int uplim = sup;
  if (PosTotal() < uplim) uplim = PosTotal();
  int lowlim = 0;
  if ((NuTransaction() - PosTotal() - sup) < 0) lowlim = sup - NuTransaction() + PosTotal();

  // calculate probability of each table

  // buffer is prepared as a member variable in class Graph
  // maybe doing zero clear is safer
  for (int ti=0;ti<NuTransaction();ti++)
    pval_log_cal_buf[ti] = -(std::numeric_limits<double>::infinity());

  double p1_log = PMinLog(sup);

  if (sup > PosTotal()){
    for (int j = 0.0; j < sup - PosTotal(); ++j){
      p1_log += log((double)(sup-j));
      p1_log -= log((double)(j + 1));
    }
  }
  double p2_log = 0.0;
  pval_log_cal_buf[uplim - lowlim] = p1_log;
  int neg_size = NuTransaction() - PosTotal();

  double p1p2_log = p1_log + p2_log;

  for (int j = uplim - 1; j >= lowlim; --j) {
    p1p2_log += log((double)(j + 1));
    p1p2_log -= log((double)(PosTotal() - j));
    p1p2_log += log((double)(neg_size - sup + j + 1));
    p1p2_log -= log((double)(sup - j));

    assert(!std::isnan(p1p2_log));
    pval_log_cal_buf[j-lowlim] = p1p2_log;
  }
  
  double p = 0.0;
  int t = lowlim;
  while (t <= uplim){
    if (t >= pos_sup){
      double p_i_log = 0.0;
      p_i_log += pval_log_cal_buf[t - lowlim];
      assert(!std::isnan(p_i_log));
      p += exp(p_i_log);
    }
    // update t_list
    t += 1;
  }
  assert(!std::isnan(p));
  return p;
}

template<typename Block>
void Database<Block>::PrepareItemVals() {
  item_info_.clear();

  // todo: do these in read item and store in struct DatabaseParameter
  int max_sup_count = -1;
  int max_pos_count = -1;

  for (std::size_t i=0 ; i < (std::size_t)NuItems() ; i++) {
    int sup = bsh_->Count(bsh_->N(data_, i));
    max_sup_count = std::max(max_sup_count, sup);
    if (has_positives_) {
      int pos_sup = bsh_->AndCount(posneg_, bsh_->N(data_, i));
      max_pos_count = std::max(max_pos_count, pos_sup);
    }
  }

  max_x_ = max_sup_count;
  if (has_positives_) max_t_ = max_pos_count;

  sup_hist_.resize(max_x_+1);
  sup_cum_hist_.resize(max_x_+1);
  reduced_item_list_.resize(max_x_+1);
  reduced_item_list_prepared_.resize(max_x_+1, false);

  InitPMinLogTable();
  if (has_positives_) InitPValTableLog();

  for (std::size_t i=0 ; i < (std::size_t)NuItems() ; i++) {
    int sup = bsh_->Count(bsh_->N(data_, i));

    sup_hist_[sup]++;

    ItemInfo new_item;
    new_item.id = i;
    new_item.sup = sup;
    new_item.pmin = PMin(sup);
    if (has_positives_) {
      int pos_sup = bsh_->AndCount(posneg_, bsh_->N(data_, i));
      new_item.pos_sup = pos_sup;
      new_item.pval = PVal(sup, pos_sup);
    }
    else new_item.pval = -1.0; // maybe wrong

    item_info_.push_back(new_item);
  }

  int total = 0;
  for (int l=max_x_;l>=0;l--) {
    total += sup_hist_[l];
    sup_cum_hist_[l] = total;
  }
}

template<typename Block>
long long int Database<Block>::Count1() const {
  long long int total=0ll;
  for (int l=0;l<=max_x_;l++) {
    total += sup_hist_[l] * l;
  }
  return total;
}

template<typename Block>
double Database<Block>::Density() const {
  long long int num1 = Count1();
  long long int all = bsh_->nu_bits * NuItems();

  return (double)num1 / (double)all;
}

template<typename Block>
int Database<Block>::PrepareReducedList(int lambda) {
  std::vector<int> & vec = reduced_item_list_[lambda];
  if (reduced_item_list_prepared_[lambda]) return vec.size();

  for (int i=0 ; i < NuItems() ; i++) {
    assert( item_info_[i].id == i);
    if (item_info_[i].sup >= lambda) {
      vec.push_back(i);
    }
  }

  reduced_item_list_prepared_[lambda] = true;
  return vec.size();
}

template<typename Block>
std::vector<int> & Database<Block>::ReducedList(int lambda) {
  if (!reduced_item_list_prepared_[lambda]) PrepareReducedList(lambda);
  std::vector<int> & vec = reduced_item_list_[lambda];
  return vec;
}

template<typename Block>
int Database<Block>::NuAllZeroTrans() const {
  Block * sup_buf = bsh_->New();

  bsh_->Reset(sup_buf);
  for (std::size_t i=0 ; i < (std::size_t)NuItems() ; i++) {
    bsh_->Or(bsh_->N(data_, i), sup_buf);
    // std::cout << "i=" << i << " ";
    // bsh_->Print(std::cout, bsh_->N(data_, i));
    // std::cout << " ";
    // bsh_->Print(std::cout, sup_buf);
    // std::cout << std::endl;
  }
  int count = bsh_->Count(sup_buf);

  // assert all one (e.g. no all-zero transaction)
  assert(bsh_->nu_bits == (std::size_t)count);

  bsh_->Delete(sup_buf);
  return (nu_transactions_ - count);
}

template<typename Block>
int Database<Block>::NuAllZeroItems() const {
  return sup_hist_[0];
}

template<typename Block>
int Database<Block>::NuReducedItems(int lambda) const {
  if (!reduced_item_list_prepared_[lambda]) return -1;
  const std::vector<int> & vec = reduced_item_list_[lambda];
  return vec.size();
}

template class DatabaseReader<uint64>;
template class Database<uint64>;

} // namespace lamp_search

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
