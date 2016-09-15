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
#include <fstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>

#include "gflags/gflags.h"

#include <sys/time.h>
#include <sys/resource.h>
#include <cstring>

#include "mpi.h"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank, nu_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nu_proc);

  struct addrinfo hints, *info, *p;
  int gai_result;

  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);

  memset(&hints, 0, sizeof hints);
  hints.ai_family = AF_UNSPEC; /*either IPV4 or IPV6*/
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_flags = AI_CANONNAME;

  if ((gai_result = getaddrinfo(hostname, "http", &hints, &info)) != 0) {
    fprintf(stderr, "getaddrinfo: %s\n", gai_strerror(gai_result));
    return 1;
  }

  for(p = info; p != NULL; p = p->ai_next) {
    printf("rank: %4d\thostname: %s\n", rank, p->ai_canonname);
    break; // hack
  }

  freeaddrinfo(info);

  MPI_Finalize();
  return 0;
}

/* Local Variables:  */
/* compile-command: "scons -u" */
/* End:              */

