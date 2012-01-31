// This program shows how to use the interface with SANDER.

#include <iostream>

#include "gmml/gmml.h"

using std::cout;
using std::endl;

using gmml::build;
using gmml::MinimizationResults;
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
    gmml::load_parameter_file("dat/parm99.dat");
    gmml::load_parameter_file("dat/Glycam_06h.dat");
    gmml::load_prep_file("dat/Glycam_06.prep");

    Structure *glycan = build_glycan();

    // The argument is a SANDER input file.
    MinimizationResults *results = glycan->minimize("dat/min.in");

    if (results == NULL) {
        cout << "Minimization failed" << endl;
    } else {
        cout << "Minimized energy: " << results->energy() << endl;
        cout << "Bond energy: " << results->bond_energy() << endl;
        cout << "VDW energy: " << results->vdw_energy() << endl;
        delete results;
    }

    delete glycan;

    return 0;
}
