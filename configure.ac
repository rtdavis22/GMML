#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([gmml], [1.0], [rtdavis22@gmail.com])
AM_INIT_AUTOMAKE([foreign silent-rules])
AC_CONFIG_SRCDIR([src/parameter_file.cc])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SUBDIRS([gtest])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
LT_INIT([dlopen])
AC_CHECK_PROG([have_sander], [sander], [yes])

if test "x${have_sander}" = xyes; then
  AC_DEFINE([HAVE_SANDER], 1, [Define to 1 if sander is in the search path])
fi

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FORK
AC_CHECK_FUNCS([gettimeofday localtime_r])

AC_CONFIG_FILES([Makefile
                 test/Makefile])

AC_OUTPUT

if test "x${have_sander}" != xyes; then
echo \
"--------------------------------------------------
Warning: sander was not found in your environment.
You will not be able to use structure minimization.
---------------------------------------------------"
fi
