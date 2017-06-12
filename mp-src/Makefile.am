CXX=$(MPICXX)
INCLUDES = -I$(top_srcdir)/src/

noinst_LIBRARIES = libmplampsearch.a
libmplampsearch_a_SOURCES = mp_dfs.cc mp_dfs.h
libmplampsearch_a_CPPFLAGS = $(MPI_CXXFLAGS) $(BOOST_CPPFLAGS) $(GFLAGS_CPPFLAGS)

