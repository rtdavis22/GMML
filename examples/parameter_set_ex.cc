// This program shows how to query a parameter set for force field information.
// See include/gmml/internal/parameter_file.h for documentation.

#include <cmath>

#include <iostream>

#include "gmml/gmml.h"

using std::cout;
using std::endl;

using gmml::kNotSet;
using gmml::ParameterFileAtom;
using gmml::ParameterFileBond;
using gmml::ParameterFileDihedral;
using gmml::ParameterFileDihedralTerm;
using gmml::ParameterFileImproperDihedral;
using gmml::ParameterFileSet;

void print_atom_type_info(const ParameterFileAtom *atom) {
    cout << "Printing info for atom type " << atom->type << ":" << endl;

    if (atom->mass == kNotSet) {
        cout << " Mass not set." << endl;
    } else {
        cout << " Mass: " << atom->mass << endl;
    }

    if (atom->radius == kNotSet) {
        cout << " Radius not set." << endl;
    } else {
        cout << " Radius: " << atom->radius << endl;
    }

    if (atom->well_depth == kNotSet) {
        cout << " Well depth not set." << endl;
    } else {
        cout << " Well depth: " << atom->well_depth << endl;
    }
}

void print_bond_type_info(const ParameterFileBond *bond) {
    cout << "Printing info for bond type " << bond->types[0] << "-" <<
           bond->types[1] << ":" << endl;

    cout << " force constant: " << bond->force_constant << endl;

    cout << " length: " << bond->length << endl;
}

void print_dihedral_term(const ParameterFileDihedralTerm *term) {
    cout << "factor: " << term->factor <<
            ", force constant: " << term->force_constant <<
            ", phase: " << term->phase <<
            ", periodicity: " << std::abs(term->periodicity) << endl;
}

void print_dihedral_type_info(const ParameterFileDihedral *dihedral) {
    cout << "Printing info for dihedral type " << dihedral->types[0] << "-" <<
            dihedral->types[1] << "-" << dihedral->types[2] << "-" <<
            dihedral->types[3] << ":" << endl;

    for (int i = 0; i < dihedral->terms.size(); i++) {
        cout << " ";
        print_dihedral_term(&dihedral->terms[i]);
    }

    if (dihedral->scee != kNotSet) {
        cout << " scee: " << dihedral->scee << endl;
    }

    if (dihedral->scnb != kNotSet) {
        cout << " scnb: " << dihedral->scnb << endl;
    }
}

int main() {
    ParameterFileSet *parmset = new ParameterFileSet;

    parmset->load("dat/parm99.dat");
    parmset->load("dat/Glycam_06h.dat");

    const ParameterFileAtom *atom = parmset->lookup("CG");
    if (atom == NULL) {
        cout << "Atom type not found." << endl;
    } else {
        print_atom_type_info(atom);
    }
    cout << endl;

    const ParameterFileBond *bond = parmset->lookup("H2", "CG");
    if (bond == NULL) {
        cout << "Bond type not found." << endl;
    } else {
        print_bond_type_info(bond);
    }
    cout << endl;

    const ParameterFileDihedral *dihedral =
            parmset->lookup("SM", "CG", "CG", "H1");
    if (dihedral == NULL) {
        cout << "Dihedral type not found." << endl;
    } else {
        print_dihedral_type_info(dihedral);
    }
    cout << endl;

    delete parmset;

    return 0;
}
