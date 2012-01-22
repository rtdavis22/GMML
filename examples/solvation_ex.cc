// This program shows how to solvate structures.
// See include/gmml/internal/solvated_structure.h for more information.

#include "gmml/gmml.h"

using gmml::build;
using gmml::load_library_file;
using gmml::load_parameter_file;
using gmml::load_prep_file;
using gmml::Structure;

// This builds the structure that corresponds to the following sequence in
// GLYCAM condensed nomenclature:
// DGalpNAcb1-4DGlcpNAcb1-2DManpa1-6[DGlcpNAcb1-4[DGalpNAcb1-4DGlcpNAcb1-2]DManpa1-3]DManpb1-4DGlcpNAcb1-4DGlcpNAcb1-OH
Structure *build_glycan() {
    Structure *substructure = build("4YB");
    substructure->attach("0VB");

    Structure *structure = build("ROH");

    structure->attach("4YB");

    structure->attach("4YB");

    structure->set_tail(2, "O4");
    structure->attach("VMB");

    structure->set_tail(3, "O3");
    structure->attach("YMA");

    structure->set_tail(4, "O2");
    structure->attach(substructure);

    structure->set_tail(4, "O4");
    structure->attach("0YB");

    structure->set_tail(3, "O6");
    structure->attach("2MA");

    structure->set_tail(8, "O2");
    structure->attach(substructure);

    delete substructure;

    return structure;
}

int main() {
    load_library_file("dat/tip3pbox.off");
    load_parameter_file("dat/parm99.dat");
    load_parameter_file("dat/Glycam_06h.dat");
    load_prep_file("dat/Glycam_06.prep");

    Structure *glycan = build_glycan();

    Structure *solvent = build("TIP3PBOX");

    Structure *solvated = solvate(*glycan, *solvent, 10.0, 1.5);

    solvated->print_amber_top_file("solvated.top");
    solvated->print_coordinate_file("solvated.rst");

    delete solvated;
    delete glycan;

    return 0;
}
