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

#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include "gflags/gflags.h"

#include "timer.h"

namespace lamp_search {

//Timer g_timer; // global timer
static Timer &timer = *(Timer::GetInstance());

long long int Timer::Start() {
  start_time_ = GetTimeInternal();
  lap_time_ = start_time_;
  return start_time_;
}

long long int Timer::Lap() {
  long long int ct = GetTimeInternal();
  long long int lt = ct - lap_time_;
  lap_time_ = ct;
  return lt;
}

long long int Timer::Elapsed() const {
  long long int ct = GetTimeInternal();
  return ct - start_time_;
}

long long int Timer::GetTimeInternal() {
  long long int nsec;
#ifdef HAVE_CLOCK_GETTIME
  struct timespec ct;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ct);
  // clock_gettime(CLOCK_MONOTONIC, &ct);
  nsec = ct.tv_sec * 1000000000 + (ct.tv_nsec);
#else
  struct timeval t;
  struct timezone tz;
  gettimeofday(&t, &tz);
  nsec = t.tv_sec * 1000000000 + t.tv_usec * 1000;
#endif
  return nsec;
}



} // namespace lamp_search
/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */
