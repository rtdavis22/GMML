SUBDIRS =  . test

DIST_SUBDIRS = ${subdirs} test

ACLOCAL_AMFLAGS = -I m4

# Build gtest before we build gmml tests.  We don't add gtest to SUBDIRS
# because then "make check" would also build and run all of gtest's own tests,
# which takes a lot of time and is generally not useful to us.  Also, we don't
# want "make install" to recurse into gtest since we don't want to overwrite
# the installed version of gtest if there is one.
check-local:
	@echo "Making lib/libgtest.a lib/libgtest_main.a in gtest"
	@cd gtest && $(MAKE) $(AM_MAKEFLAGS) lib/libgtest.la lib/libgtest_main.la

# We would like to clean gtest when "make clean" is invoked.  But we have to
# be careful because clean-local is also invoked during "make distclean", but
# "make distclean" already recurses into gtest because it's listed among the
# DIST_SUBDIRS.  distclean will delete gtest/Makefile, so if we then try to
# cd to the directory again and "make clean" it will fail.  So, check that the
# Makefile exists before recursing.
clean-local:
	@if test -e gtest/Makefile; then \
	  echo "Making clean in gtest"; \
	  cd gtest && $(MAKE) $(AM_MAKEFLAGS) clean; \
	fi

gmml_baseincludedir = $(includedir)/gmml
gmml_baseinclude_HEADERS = include/gmml/gmml.h

gmml_includedir = $(includedir)/gmml/internal
gmml_include_HEADERS = \
  include/gmml/internal/amber_top_builder.h \
  include/gmml/internal/amber_top_file.h \
  include/gmml/internal/array_tree.h \
  include/gmml/internal/atom.h \
  include/gmml/internal/bimap.h \
  include/gmml/internal/boxed_structure.h \
  include/gmml/internal/carbohydrates.h \
  include/gmml/internal/complete_residue.h \
  include/gmml/internal/coordinate_file.h \
  include/gmml/internal/coordinate_grid.h \
  include/gmml/internal/element.h \
  include/gmml/internal/fasta_sequence.h \
  include/gmml/internal/environment.h \
  include/gmml/internal/generic_type.h \
  include/gmml/internal/geometry.h \
  include/gmml/internal/glycam_code_set.h \
  include/gmml/internal/glycam_parser.h \
  include/gmml/internal/glycan_conformation_builder.h \
  include/gmml/internal/glycan_drawer.h \
  include/gmml/internal/graph.h \
  include/gmml/internal/library_file.h \
  include/gmml/internal/naccess.h \
  include/gmml/internal/netoglyc.h \
  include/gmml/internal/parameter_file.h \
  include/gmml/internal/parameter_file-inl.h \
  include/gmml/internal/pdb_file.h \
  include/gmml/internal/pdb_file_builder.h \
  include/gmml/internal/pdb_file_structure.h \
  include/gmml/internal/prep_file.h \
  include/gmml/internal/proteins.h \
  include/gmml/internal/residue_classification.h \
  include/gmml/internal/residue.h \
  include/gmml/internal/sander_minimize.h \
  include/gmml/internal/scoped_ptr.h \
  include/gmml/internal/sequence_parser.h \
  include/gmml/internal/solvated_structure.h \
  include/gmml/internal/structure.h \
  include/gmml/internal/structure-inl.h \
  include/gmml/internal/tree.h \
  include/gmml/internal/tree_residue.h

gmml_stubs_includedir = $(includedir)/gmml/internal/stubs
gmml_stubs_include_HEADERS = \
  include/gmml/internal/stubs/common.h \
  include/gmml/internal/stubs/file.h \
  include/gmml/internal/stubs/logging.h \
  include/gmml/internal/stubs/utils.h


lib_LIBRARIES = libgmml.a
libgmml_a_SOURCES = \
  src/amber_top_builder.cc \
  src/amber_top_file.cc \
  src/boxed_structure.cc \
  src/complete_residue.cc \
  src/coordinate_file.cc \
  src/element.cc \
  src/fasta_sequence.cc \
  src/environment.cc \
  src/generic_type.cc \
  src/geometry.cc \
  src/glycam_code_set.cc \
  src/glycam_parser.cc \
  src/glycan_conformation_builder.cc \
  src/glycan_drawer.cc \
  src/graph.cc \
  src/library_file.cc \
  src/logging.cc \
  src/naccess.cc \
  src/netoglyc.cc \
  src/parameter_file.cc \
  src/pdb_file.cc \
  src/pdb_file_builder.cc \
  src/pdb_file_structure.cc \
  src/prep_file.cc \
  src/proteins.cc \
  src/residue.cc \
  src/sander_minimize.cc \
  src/solvated_structure.cc \
  src/structure.cc \
  src/stubs/file.cc \
  src/utilities.h

libgmml_a_CPPFLAGS = -I$(top_srcdir)/include

uninstall-hook:
	rm -rf $(includedir)/gmml


EXTRA_DIST = \
  README \
  autogen.sh \
  examples/Makefile \
  examples/amber_topology_file_ex.cc \
  examples/carbohydrates_ex.cc \
  examples/dat/1RVZ_New.pdb \
  examples/dat/all_amino94.lib \
  examples/dat/all_aminoct94.lib \
  examples/dat/all_aminont94.lib \
  examples/dat/Glycam_06h.dat \
  examples/dat/Glycam_06.prep \
  examples/dat/HOH.prep \
  examples/dat/min.in \
  examples/dat/parm99.dat \
  examples/dat/tip3pbox.off \
  examples/minimization_ex.cc \
  examples/parameter_set_ex.cc \
  examples/pdb_ex.cc \
  examples/solvation_ex.cc \
  Makefile

MAINTAINERCLEANFILES = \
  aclocal.m4 \
  config.h \
  config.h.in \
  config.h.in~ \
  config.log \
  config.status \
  configure depcomp \
  INSTALL \
  install-sh \
  Makefile \
  Makefile.in \
  missing stamp-h1 \
  config.guess \
  config.sub \
  ltmain.sh \
  config.g~ \
  m4/ltsugar.m4 \
  m4/libtool.m4 \
  m4/ltversion.m4 \
  m4/lt~obsolete.m4 \
  m4/ltoptions.m4
