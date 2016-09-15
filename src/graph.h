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

#ifndef _LAMP_SEARCH_GRAPH_H_
#define _LAMP_SEARCH_GRAPH_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <list>
#include <limits>

#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>

#include "table.h"
#include "utils.h"

namespace lamp_search {

class Graph {
 public:
  Graph(const Table & t);
  ~Graph();

  const Table & GetTable() const {return t_; }

  class Node {
   public:
    Node();
 
    /** used for replacement scheme in Hash */
    int getImportance() const { return visits_; }

    // iterator of edges

    void Init();

    friend class ChildIterator;

    class ChildIterator {
     public:
      ChildIterator(const Graph::Node & nd, int nu_items) :
          node_ (nd),
          item_id_ (0),
          sorted_item_id_ (0),
          nu_items_ (nu_items)
      {
        SkipDuplicates();
      }

      void operator++() {
        item_id_ ++;
        SkipDuplicates();
      }

      int operator*() const {
        return item_id_;
      }

      //bool operator!=();

      // return true while valid. return false if end
      operator bool() const {
        assert(item_id_ >= 0);
        return (item_id_ < nu_items_);
      }

     private:
      const Node & node_;

      void SkipDuplicates() {
        while(true) {
          if (node_.sorted_items_[sorted_item_id_] > item_id_) break;
          if (item_id_ >= nu_items_) break;
          sorted_item_id_ ++;
          if (node_.sorted_items_[sorted_item_id_] == item_id_)
            item_id_ ++;
        }
      }

      int item_id_;
      int sorted_item_id_;
      int nu_items_;
    }; // class ChildIterator

   private:
    boost::array<int, kMaxItemSet> sorted_items_; // default sort based on index
    // todo:
    // implement heuristically_sorted_items_
    int depth_;
    long long int visits_;

    double pmin_;
    double pval_;
  }; // class Node


 private:
  // hash table of Node

  double Evaluate(const Node & n);

  const Table & t_;
};

} // namespace lamp_search

#endif // _LAMP_SEARCH_GRAPH_H_
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
