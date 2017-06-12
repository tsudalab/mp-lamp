#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([MP-LAMP], [3.2.4], [lamp_staff@googlegroups.com])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
AC_PROG_RANLIB
AC_LANG([C++])
ACX_MPI(,[AC_MSG_FAILURE([MPI compiler requested, but could not use MPI.])])

# Checks for libraries.
BOOST_REQUIRE([1.41.0])

# Check if there is google-gflags library installed.
#AC_CHECK_LIB([gflags], [main],[], [AC_MSG_FAILURE([gflags test failed])])
SAVE_CXXFLAGS="$CXXFLAGS"
SAVE_LIBS="$LIBS"
AC_ARG_WITH([gflags], [AS_HELP_STRING([--with-gflags=DIR],[prefix of Google gflags library])],
  [GFLAGS_CPPFLAGS="-I${with_gflags}/include"
   GFLAGS_LIBS="-L${with_gflags}/lib -lgflags -lpthread"
   CXXFLAGS="$CXXFLAGS $GFLAGS_CPPFLAGS"
   LIBS="$LIBS $GFLAGS_LIBS"]
)
AC_CHECK_LIB([gflags], [main],
  [if test x"$GFLAGS_LIBS" = x""; then
     GFLAGS_LIBS="-lgflags -lpthread"
   fi],
  AC_MSG_FAILURE([gflags required !!]))
CXXFLAGS="$SAVE_CXXFLAGS"
LIBS="$SAVE_LIBS"
AC_SUBST(GFLAGS_CPPFLAGS)
AC_SUBST(GFLAGS_LIBS)

# Checks for header files.
AC_CHECK_HEADERS([arpa/inet.h fcntl.h float.h limits.h netdb.h stddef.h stdint.h stdlib.h string.h strings.h sys/socket.h sys/time.h sys/timeb.h unistd.h wchar.h wctype.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_FUNC_STRERROR_R
AC_CHECK_FUNCS([clock_gettime dup2 fchdir floor getcwd gethostname getpagesize gettimeofday memset mkdir munmap regcomp rmdir socket sqrt strcasecmp strchr strdup strerror strrchr strstr strtol strtoull])

AC_CONFIG_FILES([Makefile
                 main/Makefile
                 mp-main/Makefile
                 mp-src/Makefile
                 mp-tests/Makefile
                 src/Makefile
                 tests/Makefile])
AC_OUTPUT
