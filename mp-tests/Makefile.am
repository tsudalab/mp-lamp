noinst_LIBRARIES = libgtest.a
libgtest_a_SOURCES = gtest/gtest-all.cc gtest/gtest_main.cc gtest/gtest.h

INCLUDES = -I$(top_srcdir)/src/ -I$(top_srcdir)/mp-src/ $(MPI_CXXFLAGS) $(BOOST_CPPFLAGS) $(GFLAGS_CPPFLAGS)
LDADD = ../mp-src/libmplampsearch.a ../src/liblampsearch.a libgtest.a $(GFLAGS_LIBS)
CXX = $(MPICXX)

TESTS = mpi_basic_unittest mp_lamp_basic_unittest

check_PROGRAMS = mpi_basic_unittest mp_lamp_basic_unittest

mpi_basic_unittest_SOURCES = mpi_basic_unittest.cc

mp_lamp_basic_unittest_SOURCES = mp_lamp_basic_unittest.cc
