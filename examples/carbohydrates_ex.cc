// This file shows carbohydrate-specific functionality.
// See include/gmml/internal/glycam_code_set.h and
// include/gmml/internal/glycan_conformation_builder.h for more information.

#include <list>
#include <string>

#include "gmml/gmml.h"

using std::list;
using std::string;

using gmml::glycam_build;
using gmml::load_parameter_file;
using gmml::load_prep_file;
using gmml::Structure;

using gmml::carbohydrate::GCBStructure;
using gmml::carbohydrate::GlycanConformationBuilder;
using gmml::carbohydrate::set_phi;

// This function builds all conformations of a structure in which all
// 2-6 linkages with an Neu5Ac at the non-reducing end can take a phi value
// of -60 or 180 and all omega values that are deemed "likely" are possible.
// See the header file referenced above for more ways to specify which
// glycosidic torsions you want to build.
list<GCBStructure*> *build_conformations(const Structure *structure) {
    GlycanConformationBuilder builder(*structure);
    builder.add_phi_value("Neu5Ac", 2, "*", 6, -60.0);
    builder.add_phi_value("Neu5Ac", 2, "*", 6, 180.0);
    builder.add_likely_omega_values();
    return builder.build();
}

int main() {
    load_parameter_file("dat/parm99.dat");
    load_parameter_file("dat/Glycam_06h.dat");
    load_prep_file("dat/Glycam_06.prep");

    // This is a sequence in GLYCAM condensed nomenclature.
    std::string sequence = "DNeup5Aca2-6DGalpb1-4DGlcpNAcb1-2DManpa1-6[DNeup5Aca2-3DGalpb1-4DGlcpNAcb1-2DManpa1-3]DManpb1-4DGlcpNAcb1-4DGlcpNAcb1-OME";

    Structure *glycan = glycam_build(sequence);

    // This sets the phi torsion between the residue with index 3 and the
    // residue at it's reducing end. When building structures from a sequence,
    // the indices of the residues are in the reverse order as their
    // appearance in the sequence. Thus, OME has index 0 and the residue with
    // index 3 is a mannose.
    set_phi(glycan, 3, -60.0);

    list<GCBStructure*> *structures = build_conformations(glycan);
    for (list<GCBStructure*>::iterator it = structures->begin();
            it != structures->end(); ++it) {
        string name = (*it)->name();
        (*it)->print_coordinate_file(name + ".rst");
        delete *it;
    }
    delete structures;

    glycan->print_amber_top_file("glycan.top");

    delete glycan;

    return 0;
}
