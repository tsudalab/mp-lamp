INCLUDES = -I$(top_srcdir)/src/ -I$(top_srcdir)/mp-src/ $(MPI_CXXFLAGS) $(BOOST_CPPFLAGS) $(GFLAGS_CPPFLAGS)
LDADD = ../mp-src/libmplampsearch.a ../src/liblampsearch.a $(GFLAGS_LIBS)
CXX = $(MPICXX)

bin_PROGRAMS = mp-lamp show_hostname

mp_lamp_SOURCES = mp-lamp.cc
mp_lamp_LDFLAGS = $(MPI_CXXLDFLAGS)

show_hostname_SOURCES = show_hostname.cc
show_hostname_LDFLAGS = $(MPI_CXXLDFLAGS)
