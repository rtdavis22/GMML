SUBDIRS = src# gtest

gmml_baseincludedir = $(includedir)/gmml
gmml_baseinclude_HEADERS = include/gmml/gmml.h

gmmlincludedir = $(includedir)/gmml/internal
gmmlinclude_HEADERS = \
  include/gmml/internal/amber_top_builder.h \
  include/gmml/internal/amber_top_file.h \
  include/gmml/internal/atom.h \
  include/gmml/internal/boxed_structure.h \
  include/gmml/internal/coordinate_file.h \
  include/gmml/internal/coordinate_grid.h \
  include/gmml/internal/environment.h \
  include/gmml/internal/generic_type.h \
  include/gmml/internal/geometry.h \
  include/gmml/internal/geometry-inl.h \
  include/gmml/internal/glycam_code_set.h \
  include/gmml/internal/glycam_parser.h \
  include/gmml/internal/graph.h \
  include/gmml/internal/library_file.h \
  include/gmml/internal/parameter_file.h \
  include/gmml/internal/parameter_file-inl.h \
  include/gmml/internal/pdb_file_builder.h \
  include/gmml/internal/pdb_file.h \
  include/gmml/internal/prep_file.h \
  include/gmml/internal/residue_classification.h \
  include/gmml/internal/residue.h \
  include/gmml/internal/sander_minimize.h \
  include/gmml/internal/sequence_parser.h \
  include/gmml/internal/solvated_structure.h \
  include/gmml/internal/structure.h \
  include/gmml/internal/torsion_combination_builder.h \
  include/gmml/internal/tree.hh \
  include/gmml/internal/tree_residue.h \
  include/gmml/internal/utilities.h

#Look at protobuf/src/Makefile.am for how I want to do this.
#check_PROGRAMS = gmml-test

#gmml_test_LDADD = $(top_builddir)/gtest/lib/libgtest.la \
#                  $(top_builddir)/gtest/lib/libgtest_main.la
#gmml_test_CPPFLAGS = -I$(top_builddir)/gtest/include
#gmml_test_SOURCES = \
#  test/geometry_test.cc

uninstall-hook:
	rmdir $(includedir)/gmml/internal
	rmdir $(includedir)/gmml

MAINTAINERCLEANFILES = \
  aclocal.m4 \
  autom4te.cache \
  config.h \
  config.h.in \
  config.log \
  config.status \
  configure depcomp \
  INSTALL \
  install-sh \
  Makefile \
  Makefile.in \
  missing stamp-h1
