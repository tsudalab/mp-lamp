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

#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
//#include <boost/range/algorithm/for_each.hpp>

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#include "gflags/gflags.h"

#include "table.h"

DEFINE_double(a, 0.05, "significance level alpha");

namespace {

int CountLines(std::istream & is) {
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

} // namespace

//==============================================================================

namespace lamp_search {

Table::Table(std::istream & item_file, std::istream & posneg_file, FunctionsSuper & functions) :
    nu_items_ (0),
    nu_transactions_ (0),
    data_ (0), // null
    posneg_ (0), // null
    max_x_ (-1),
    max_t_ (-1),
	functions(functions) {
  ReadItems(item_file);
  ReadPosNeg(posneg_file);

  SetSigLev(FLAGS_a); // double

  functions.setAllSize(NuTransaction());
  functions.setN1(PosTotal());
  PrepareItemVals();
}

Table::~Table() {
  if (data_) delete data_;
  if (posneg_) delete posneg_;
}

void Table::SetValuesForTest(int nu_item, int nu_transaction, int nu_pos_total) {
  nu_items_ = nu_item;
  nu_transactions_ = nu_transaction;
  nu_pos_total_ = nu_pos_total;

  // should be given as paramters?
  max_x_ = nu_transaction;
  max_t_ = nu_pos_total;

  functions.setAllSize(NuTransaction());
  functions.setN1(PosTotal());
  // InitPMinTable();
  // InitPValTable();
}

void Table::InitPMinTable() {
  pmin_table_.resize(max_x_ + 1);
  pmin_table_[0] = 1.0;
  for (int i = 1; i <= max_x_; i++)
    //    pmin_table_[i] = PMinCal(i);
    pmin_table_[i] = functions.funcF(i);
}

std::ostream & Table::DumpPMinTable(std::ostream & out) const {
  std::stringstream s;

  for (int i=0;i<=max_x_;i++)
    s << "i=" << i << "\tpmin=" << pmin_table_[i] << std::endl;

  out << s.str() << std::flush;
  return out;
}

void Table::InitPValTable() {
  pval_table_.resize( (max_x_+1) * (max_t_+1) );
  for (int i=0;i<=max_x_;i++) { // sup == x
    for (int j=0;j<=max_t_;j++) { // pos_sup == t
      if (j>i) continue;
      pval_table_[i*(max_t_+1) + j] = PValCal(i, j);
    }
  }
}

std::ostream & Table::DumpPValTable(std::ostream & out) const {
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

void Table::ReadItems(std::istream & is) {
  nu_items_ = 0;
  nu_transactions_ = 0;
  max_x_ = -1;
  max_t_ = -1;
  // read 1st line, tokenize, count items
  // read line, tokenize

  std::string line;
  std::string trimmed_line;
  boost::char_separator<char> sep(", ");

  int nu_lines = CountLines(is);
  if (nu_lines <= 0) {
    std::cout << "file line <= 0" << std::endl;
    exit(1);
  }

  nu_transactions_ = nu_lines - 1;
  transaction_names_.reserve(nu_transactions_);
  max_item_in_transaction_ = -1;

  {
    std::getline(is, line);
    trimmed_line = boost::algorithm::trim_copy(line);
    if ( ! is.good() ) {
      std::cout << "error in the 1st line: " << trimmed_line << std::endl;
      exit(1);
    }
    boost::tokenizer< boost::char_separator<char> > tokens(trimmed_line, sep);
    
    boost::tokenizer< boost::char_separator<char> >::iterator tok_iter=tokens.begin();
    ++ tok_iter; // skip first element. should be "#gene"
    int counter = 0;
    for (; tok_iter != tokens.end(); ++tok_iter) {
      item_names_.push_back (*tok_iter);
      counter++;
    }
    nu_items_ = counter;
  }

  const boost::dynamic_bitset<> empty_bitset(nu_transactions_);
  data_ = new std::vector< boost::dynamic_bitset<> >(nu_items_, empty_bitset);

  int transaction_counter = 0;
  while (1) {
    std::getline(is, line);
    trimmed_line = boost::algorithm::trim_copy(line);
    // eof, fail, bad
    if ( ! is.good() ) break;

    boost::tokenizer< boost::char_separator<char> > tokens(trimmed_line, sep);
    boost::tokenizer< boost::char_separator<char> >::iterator tok_iter=tokens.begin();
    transaction_names_.push_back(*tok_iter); // first element is transaction name
    ++ tok_iter;
    int counter = 0;
    int pos_counter = 0;
    for (; tok_iter != tokens.end(); ++tok_iter) {
      if (*tok_iter == "0")
        (*data_)[counter][transaction_counter] = false;
      else {
        (*data_)[counter][transaction_counter] = true;
        pos_counter++;
      }
      counter++;
    }
    max_item_in_transaction_ = std::max(max_item_in_transaction_, pos_counter);
    transaction_counter++;
    assert( nu_items_ == counter );
  }
  assert( nu_transactions_ == transaction_counter );
}

void Table::ReadPosNeg(std::istream & is) {
  std::string line;
  std::string trimmed_line;
  boost::char_separator<char> sep(", ");

#ifndef NDEBUG
  int nu_lines = CountLines(is);
  assert( nu_transactions_ == nu_lines - 1 );
#endif
  
  {
    std::getline(is, line); // skip 1st line
  }

  posneg_ = new boost::dynamic_bitset<>(nu_transactions_);

  int transaction_counter = 0;
  while (1) {
    std::getline(is, line);
    trimmed_line = boost::algorithm::trim_copy(line);
    // eof, fail, bad
    if ( ! is.good() ) break;

    boost::tokenizer< boost::char_separator<char> > tokens(trimmed_line, sep);
    boost::tokenizer< boost::char_separator<char> >::iterator tok_iter=tokens.begin();
    transaction_names_.push_back(*tok_iter); // first element is transaction name
    ++ tok_iter;

    {
      if (*tok_iter == "0")
        (*posneg_)[transaction_counter] = false;
      else
        (*posneg_)[transaction_counter] = true;
    }
    transaction_counter++;
  }

  nu_pos_total_ = posneg_->count();
  assert( nu_transactions_ == transaction_counter );
}

std::ostream & Table::DumpItems(std::ostream & out) const {
  std::stringstream s;

  BOOST_FOREACH (const boost::dynamic_bitset<> bs, *data_) {
    s << bs << std::endl;
  }  

  out << s.str() << std::flush;
  return out;
}

std::ostream & Table::DumpPosNeg(std::ostream & out) const {
  std::stringstream s;
  s << (*posneg_) << std::endl;

  out << s.str() << std::flush;
  return out;
}

std::ostream & Table::ShowInfo(std::ostream & out) const {
  std::stringstream s;

  s << "nu_transactions=" << std::setw(12) << nu_transactions_
    << "\tnu_items="      << std::setw(12) << nu_items_
    << "\tnu_pos_total="  << std::setw(12) << nu_pos_total_
    << "\tmax_x="         << std::setw(12) << max_x_
    << "\tmax_t="         << std::setw(12) << max_t_
    << "\tmax_1in_trans=" << std::setw(12) << max_item_in_transaction_
    << std::endl;
  
  out << s.str() << std::flush;
  return out;
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

double Table::PValCal(int sup, int pos_sup) const {
  double p = functions.calPValue(sup, pos_sup, empty, empty);
  return p;
}

void Table::PrepareItemVals() {
  item_info_.clear();

  int max_sup_count = -1;
  int max_pos_count = -1;

  for (std::size_t i=0 ; i < (std::size_t)NuItems() ; i++) {
    int sup = NthData(i).count();
    max_sup_count = std::max(max_sup_count, sup);
    int pos_sup = (NthData(i) & PosNeg()).count();
    max_pos_count = std::max(max_pos_count, pos_sup);
  }

  max_x_ = max_sup_count;
  max_t_ = max_pos_count;

  InitPMinTable();
  InitPValTable();

  for (std::size_t i=0 ; i < (std::size_t)NuItems() ; i++) {
    int sup = NthData(i).count();
    int pos_sup = (NthData(i) & PosNeg()).count();

    ItemInfo new_item;
    new_item.id = i;
    new_item.sup = sup;
    new_item.pos_sup = pos_sup;
    new_item.pmin = PMin(sup);
    new_item.pval = PVal(sup, pos_sup);

    assert(!std::isnan(new_item.pval));

    item_info_.push_back(new_item);
  }

  std::sort(item_info_.begin(), item_info_.end() );  
  double best_pval = 1.0;
  // int best_pval_item = -1;
  for (std::size_t i=0 ; i < (std::size_t)NuItems() ; i++) {
    if (item_info_[i].pval < best_pval) {
      best_pval = item_info_[i].pval;
      // best_pval_item = i;
    }
  }  
}

} // namespace lamp_search
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
