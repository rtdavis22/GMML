check_PROGRAMS = gmml-test

gmml_test_LDADD = $(top_builddir)/libgmml.a \
                  $(top_builddir)/gtest/lib/libgtest.la \
                  $(top_builddir)/gtest/lib/libgtest_main.la
gmml_test_CPPFLAGS = -I$(top_builddir)/include \
                     -I$(top_builddir)/gtest/include
gmml_test_SOURCES = \
  atom_test.cc \
  geometry_test.cc \
  prep_file_test.cc

TESTS = gmml-test

MAINTAINERCLEANFILES = Makefile.in
