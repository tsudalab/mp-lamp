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

#include <iostream>
#include <iomanip>

#include "gtest/gtest.h"

#include "utils.h"
#include "variable_length_itemset.h"

using namespace lamp_search;

TEST (VariableLengthItemsetTest, InitTest) {
  VariableLengthItemsetStack * vli1;
  vli1 = new VariableLengthItemsetStack(10);
  // vli1 = new VariableLengthItemsetStack(10, 5);
  delete vli1;

  VariableLengthItemsetStack vli2(7);
  // VariableLengthItemsetStack vli2(7, 4);
}

TEST (VariableLengthItemsetTest, OperationTest) {
  VariableLengthItemsetStack * s;
  s = new VariableLengthItemsetStack(16); // capacity
  //s = new VariableLengthItemsetStack(16, 5); // capacity, sup_max

  int * p;

  s->PushPre();
  p = s->Top();
  s->SetOneItem(p, 0, 1); // index, item_num, item_id
  s->SetOneItem(p, 1, 5); // index, item_num, item_id
  s->SetItemNum(p, 2);
  s->SetSup(p, 3);
  s->PushPost();

  s->Print(std::cout, p);
  EXPECT_EQ(2, s->GetItemNum(p));
  EXPECT_EQ(3, s->GetSup(p));
  EXPECT_EQ(1, s->NuItemset());
  EXPECT_EQ(7, s->UsedCapacity()); // 4 + 3 for sentinel and timestamp

  s->PushPre();
  p = s->Top();
  s->SetOneItem(p, 0, 2); // index, item_num, item_id
  s->SetOneItem(p, 1, 3); // index, item_num, item_id
  s->SetOneItem(p, 2, 4); // index, item_num, item_id
  s->SetItemNum(p, 3);
  s->SetSup(p, 2);
  s->PushPost();

  s->Print(std::cout, p);
  EXPECT_EQ(3, s->GetItemNum(p));
  EXPECT_EQ(2, s->GetSup(p));
  EXPECT_EQ(2, s->NuItemset());
  EXPECT_EQ(12, s->UsedCapacity()); // 9 + 3 for sentinel and timestamp

  s->Pop();
  EXPECT_EQ(1, s->NuItemset());

  p = s->Top();
  s->Print(std::cout, p);
  EXPECT_EQ(2, s->GetItemNum(p));
  EXPECT_EQ(3, s->GetSup(p));
  EXPECT_EQ(1, s->NuItemset());
  EXPECT_EQ(7, s->UsedCapacity()); // 4 + 3 for sentinel and timestamp

  s->Pop();
  EXPECT_EQ(0, s->NuItemset());

  delete s;
}

TEST (VariableLengthItemsetTest, PushTest) {
  VariableLengthItemsetStack * s;
  s = new VariableLengthItemsetStack(12); // capacity
  // s = new VariableLengthItemsetStack(12, 5); // capacity, sup_max

  int * p;

  std::cout << "top=" << s->Top() << std::endl;

  s->PushPre();
  p = s->Top();
  std::cout << "0:\t";
  std::cout << "top=" << s->Top() << "\t";
  s->Print(std::cout, p);
  s->PushOneItem(5);
  std::cout << "1:\t";
  std::cout << "top=" << s->Top() << "\t";
  s->Print(std::cout, p);
  s->PushOneItem(1);
  std::cout << "2:\t";
  std::cout << "top=" << s->Top() << "\t";
  s->Print(std::cout, p);
  s->SetSup(p, 3);
  s->PushPost();

  s->Print(std::cout, p);
  EXPECT_EQ(2, s->GetItemNum(p));
  EXPECT_EQ(3, s->GetSup(p));
  EXPECT_EQ(1, s->NuItemset());
  EXPECT_EQ(7, s->UsedCapacity()); // 4 + 3 for sentinel and timestamp

  s->PushPre();
  p = s->Top();
  s->PushOneItem(2);
  s->PushOneItem(3);
  s->PushOneItem(4);
  s->SetSup(p, 2);
  s->PushPost();

  s->Print(std::cout, p);
  EXPECT_EQ(3, s->GetItemNum(p));
  EXPECT_EQ(2, s->GetSup(p));
  EXPECT_EQ(2, s->NuItemset());
  EXPECT_EQ(12, s->UsedCapacity()); // 9 + 3 for sentinel and timestamp

  s->Pop();
  EXPECT_EQ(1, s->NuItemset());

  p = s->Top();
  s->Print(std::cout, p);
  EXPECT_EQ(2, s->GetItemNum(p));
  EXPECT_EQ(3, s->GetSup(p));
  EXPECT_EQ(1, s->NuItemset());
  EXPECT_EQ(7, s->UsedCapacity()); // 4 + 3 for sentinel and timestamp

  s->Pop();
  EXPECT_EQ(0, s->NuItemset());

  delete s;
}

TEST (VariableLengthItemsetTest, SplitTestEven) {
  VariableLengthItemsetStack * src;
  src = new VariableLengthItemsetStack(16); // capacity
  // src = new VariableLengthItemsetStack(16, 5); // capacity, sup_max
  VariableLengthItemsetStack * dst;
  dst = new VariableLengthItemsetStack(16); // capacity
  // dst = new VariableLengthItemsetStack(16, 5); // capacity, sup_max
  int * p;

  // even number items
  src->PushPre();
  p = src->Top();
  src->PushOneItem(5);
  src->PushOneItem(1);
  src->SetSup(p, 3);
  src->PushPost();

  src->PushPre();
  p = src->Top();
  src->PushOneItem(2);
  src->PushOneItem(3);
  src->PushOneItem(4);
  src->SetSup(p, 2);
  src->PushPost();

  EXPECT_EQ(3, src->GetItemNum(p));
  EXPECT_EQ(2, src->GetSup(p));
  EXPECT_EQ(2, src->NuItemset());
  EXPECT_EQ(12, src->UsedCapacity()); // 9 + 3 for sentinel and timestamp

  src->Split(dst);

  std::cout << "src\n";

  p = src->FirstItemset();
  src->Print(std::cout, p);
  EXPECT_EQ(2, src->GetItemNum(p));
  EXPECT_EQ(3, src->GetSup(p));

  p = src->NextItemset(p);
  EXPECT_EQ(NULL, p);

  EXPECT_EQ(1, src->NuItemset());
  EXPECT_EQ(7, src->UsedCapacity()); // 4 + 3 for sentinel and timestamp

  std::cout << "dst\n";

  p = dst->FirstItemset();
  dst->Print(std::cout, p);
  EXPECT_EQ(3, dst->GetItemNum(p));
  EXPECT_EQ(2, dst->GetSup(p));

  p = dst->NextItemset(p);
  EXPECT_EQ(NULL, p);

  EXPECT_EQ(1, dst->NuItemset());
  EXPECT_EQ(8, dst->UsedCapacity()); // 5 + 3 for sentinel and timestamp

  delete src;
  delete dst;
}

TEST (VariableLengthItemsetTest, SplitTestOdd) {
  VariableLengthItemsetStack * src;
  src = new VariableLengthItemsetStack(20); // capacity
  // src = new VariableLengthItemsetStack(20, 5); // capacity, sup_max
  VariableLengthItemsetStack * dst;
  dst = new VariableLengthItemsetStack(20); // capacity
  // dst = new VariableLengthItemsetStack(20, 5); // capacity, sup_max
  int * p;

  // odd number items
  src->PushPre();
  p = src->Top();
  src->PushOneItem(5);
  src->PushOneItem(1);
  src->SetSup(p, 3);
  src->PushPost();

  src->PushPre();
  p = src->Top();
  src->PushOneItem(2);
  src->PushOneItem(3);
  src->PushOneItem(4);
  src->SetSup(p, 2);
  src->PushPost();

  src->PushPre();
  p = src->Top();
  src->PushOneItem(6);
  src->PushOneItem(2);
  src->PushOneItem(4);
  src->PushOneItem(1);
  src->SetSup(p, 1);
  src->PushPost();

  EXPECT_EQ(4, src->GetItemNum(p));
  EXPECT_EQ(1, src->GetSup(p));
  EXPECT_EQ(3, src->NuItemset());
  EXPECT_EQ(18, src->UsedCapacity());

  src->Split(dst);

  std::cout << "src\n";

  p = src->FirstItemset();
  src->Print(std::cout, p);
  EXPECT_EQ(2, src->GetItemNum(p));
  EXPECT_EQ(3, src->GetSup(p));

  p = src->NextItemset(p);
  src->Print(std::cout, p);
  EXPECT_EQ(4, src->GetItemNum(p));
  EXPECT_EQ(1, src->GetSup(p));

  p = src->NextItemset(p);
  EXPECT_EQ(NULL, p);

  EXPECT_EQ(2, src->NuItemset());
  EXPECT_EQ(13, src->UsedCapacity());

  std::cout << "dst\n";

  p = dst->FirstItemset();
  dst->Print(std::cout, p);
  EXPECT_EQ(3, dst->GetItemNum(p));
  EXPECT_EQ(2, dst->GetSup(p));

  p = dst->NextItemset(p);
  EXPECT_EQ(NULL, p);

  EXPECT_EQ(1, dst->NuItemset());
  EXPECT_EQ(8, dst->UsedCapacity());

  delete src;
  delete dst;
}

TEST (VariableLengthItemsetTest, SplitTestOne) {
  VariableLengthItemsetStack * src;
  src = new VariableLengthItemsetStack(20); // capacity
  // src = new VariableLengthItemsetStack(20, 5); // capacity, sup_max
  VariableLengthItemsetStack * dst;
  dst = new VariableLengthItemsetStack(20); // capacity
  // dst = new VariableLengthItemsetStack(20, 5); // capacity, sup_max
  int * p;

  // odd number items
  src->PushPre();
  p = src->Top();
  src->PushOneItem(5);
  src->PushOneItem(1);
  src->SetSup(p, 3);
  src->PushPost();

  EXPECT_EQ(2, src->GetItemNum(p));
  EXPECT_EQ(3, src->GetSup(p));
  EXPECT_EQ(1, src->NuItemset());
  EXPECT_EQ(6, src->UsedCapacity());

  src->Split(dst);

  std::cout << "src\n";

  p = src->FirstItemset();
  src->Print(std::cout, p);
  EXPECT_EQ(2, src->GetItemNum(p));
  EXPECT_EQ(3, src->GetSup(p));

  p = src->NextItemset(p);
  EXPECT_EQ(NULL, p);

  EXPECT_EQ(1, src->NuItemset());
  EXPECT_EQ(7, src->UsedCapacity());

  std::cout << "dst\n";

  p = dst->FirstItemset();
  EXPECT_EQ(NULL, p);

  EXPECT_EQ(0, dst->NuItemset());
  EXPECT_EQ(3, dst->UsedCapacity());

  delete src;
  delete dst;
}

TEST (VariableLengthItemsetTest, MergeTest) {
  VariableLengthItemsetStack * src;
  src = new VariableLengthItemsetStack(20); // capacity
  // src = new VariableLengthItemsetStack(20, 5); // capacity, sup_max
  VariableLengthItemsetStack * dst;
  dst = new VariableLengthItemsetStack(20); // capacity
  // dst = new VariableLengthItemsetStack(20, 5); // capacity, sup_max
  int * p;

  // src

  src->PushPre();
  p = src->Top();
  src->PushOneItem(5);
  src->PushOneItem(1);
  src->SetSup(p, 3);
  src->PushPost();

  src->PushPre();
  p = src->Top();
  src->PushOneItem(2);
  src->PushOneItem(3);
  src->PushOneItem(4);
  src->SetSup(p, 2);
  src->PushPost();

  EXPECT_EQ(3, src->GetItemNum(p));
  EXPECT_EQ(2, src->GetSup(p));
  EXPECT_EQ(2, src->NuItemset());
  EXPECT_EQ(12, src->UsedCapacity());

  // dst

  dst->PushPre();
  p = dst->Top();
  dst->PushOneItem(6);
  dst->PushOneItem(2);
  dst->PushOneItem(4);
  dst->PushOneItem(1);
  dst->SetSup(p, 1);
  dst->PushPost();

  EXPECT_EQ(4, dst->GetItemNum(p));
  EXPECT_EQ(1, dst->GetSup(p));
  EXPECT_EQ(1, dst->NuItemset());
  EXPECT_EQ(9, dst->UsedCapacity());

  // merge

  src->Merge(dst);

  EXPECT_EQ(3, src->NuItemset());
  EXPECT_EQ(18, src->UsedCapacity());

  p = src->FirstItemset();
  EXPECT_NE((int *)NULL, p);
  src->Print(std::cout, p);
  EXPECT_EQ(2, src->GetItemNum(p));
  EXPECT_EQ(3, src->GetSup(p));

  p = src->NextItemset(p);
  EXPECT_NE((int *)NULL, p);
  src->Print(std::cout, p);
  EXPECT_EQ(3, src->GetItemNum(p));
  EXPECT_EQ(2, src->GetSup(p));

  p = src->NextItemset(p);
  EXPECT_NE((int *)NULL, p);
  src->Print(std::cout, p);
  EXPECT_EQ(4, src->GetItemNum(p));
  EXPECT_EQ(1, src->GetSup(p));

  p = src->NextItemset(p);
  EXPECT_EQ(NULL, p);

  // clear

  dst->Clear();

  p = dst->FirstItemset();
  EXPECT_EQ(NULL, p);

  EXPECT_EQ(0, dst->NuItemset());
  EXPECT_EQ(3, dst->UsedCapacity());

  delete src;
  delete dst;
}

TEST (VariableLengthItemsetTest, MergeStackTest) {
  VariableLengthItemsetStack * src;
  src = new VariableLengthItemsetStack(20); // capacity
  // src = new VariableLengthItemsetStack(20, 5); // capacity, sup_max
  VariableLengthItemsetStack * dst;
  dst = new VariableLengthItemsetStack(20); // capacity
  // dst = new VariableLengthItemsetStack(20, 5); // capacity, sup_max
  int * p;

  // src

  src->PushPre();
  p = src->Top();
  src->PushOneItem(5);
  src->PushOneItem(1);
  src->SetSup(p, 3);
  src->PushPost();

  src->PushPre();
  p = src->Top();
  src->PushOneItem(2);
  src->PushOneItem(3);
  src->PushOneItem(4);
  src->SetSup(p, 2);
  src->PushPost();

  EXPECT_EQ(3, src->GetItemNum(p));
  EXPECT_EQ(2, src->GetSup(p));
  EXPECT_EQ(2, src->NuItemset());
  EXPECT_EQ(12, src->UsedCapacity());

  // dst

  dst->PushPre();
  p = dst->Top();
  dst->PushOneItem(6);
  dst->PushOneItem(2);
  dst->PushOneItem(4);
  dst->PushOneItem(1);
  dst->SetSup(p, 1);
  dst->PushPost();

  EXPECT_EQ(4, dst->GetItemNum(p));
  EXPECT_EQ(1, dst->GetSup(p));
  EXPECT_EQ(1, dst->NuItemset());
  EXPECT_EQ(9, dst->UsedCapacity());

  // merge

  src->MergeStack((dst->Stack())+VariableLengthItemsetStack::SENTINEL+1,
                  dst->UsedCapacity()-VariableLengthItemsetStack::SENTINEL-1);

  EXPECT_EQ(3, src->NuItemset());
  EXPECT_EQ(18, src->UsedCapacity());

  p = src->FirstItemset();
  EXPECT_NE((int *)NULL, p);
  src->Print(std::cout, p);
  EXPECT_EQ(2, src->GetItemNum(p));
  EXPECT_EQ(3, src->GetSup(p));

  p = src->NextItemset(p);
  EXPECT_NE((int *)NULL, p);
  src->Print(std::cout, p);
  EXPECT_EQ(3, src->GetItemNum(p));
  EXPECT_EQ(2, src->GetSup(p));

  p = src->NextItemset(p);
  EXPECT_NE((int *)NULL, p);
  src->Print(std::cout, p);
  EXPECT_EQ(4, src->GetItemNum(p));
  EXPECT_EQ(1, src->GetSup(p));

  p = src->NextItemset(p);
  EXPECT_EQ(NULL, p);

  // clear

  dst->Clear();

  p = dst->FirstItemset();
  EXPECT_EQ(NULL, p);

  EXPECT_EQ(0, dst->NuItemset());
  EXPECT_EQ(3, dst->UsedCapacity());

  delete src;
  delete dst;
}

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
