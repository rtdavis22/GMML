check_PROGRAMS = gmml-test

gmml_test_LDADD = $(top_builddir)/libgmml.a \
                  $(top_builddir)/gtest/lib/libgtest.la \
                  $(top_builddir)/gtest/lib/libgtest_main.la
gmml_test_CPPFLAGS = -I$(top_builddir)/include \
                     -I$(top_builddir)/gtest/include
gmml_test_SOURCES = \
  atom_test.cc \
  coordinate_file_test.cc \
  fasta_sequence_test.cc \
  geometry_test.cc \
  glycam_parser_test.cc \
  graph_test.cc \
  netoglyc_test.cc \
  pdb_test.cc \
  prep_file_test.cc \
  sander_minimize_test.cc

TESTS = gmml-test

EXTRA_DIST = \
  dat/1RVZ_New.pdb \
  dat/bad_mdout \
  dat/netoglyc.out \
  dat/mdout \
  dat/test.prep

MAINTAINERCLEANFILES = Makefile.in
