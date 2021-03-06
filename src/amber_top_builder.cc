// Copyright (c) 2012 The University of Georgia
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// Author: Robert Davis

#include "gmml/internal/amber_top_builder.h"

#include <cmath>

#include <algorithm>
#include <deque>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "gmml/internal/amber_top_file.h"
#include "gmml/internal/atom.h"
#include "gmml/internal/boxed_structure.h"
#include "gmml/internal/element.h"
#include "gmml/internal/environment.h"
#include "gmml/internal/geometry.h"
#include "gmml/internal/graph.h"
#include "gmml/internal/parameter_file.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/solvated_structure.h"
#include "gmml/internal/structure.h"
#include "utilities.h"

using std::deque;
using std::string;
using std::vector;

namespace gmml {
namespace {

// The names of the topology file sections, in the order they should appear in
// the file.
const char *kSectionList[] = {
    "TITLE",
    "POINTERS",
    "ATOM_NAME",
    "CHARGE",
    "MASS",
    "ATOM_TYPE_INDEX",
    "NUMBER_EXCLUDED_ATOMS",
    "NONBONDED_PARM_INDEX",
    "RESIDUE_LABEL",
    "RESIDUE_POINTER",
    "BOND_FORCE_CONSTANT",
    "BOND_EQUIL_VALUE",
    "ANGLE_FORCE_CONSTANT",
    "ANGLE_EQUIL_VALUE",
    "DIHEDRAL_FORCE_CONSTANT",
    "DIHEDRAL_PERIODICITY",
    "DIHEDRAL_PHASE",
    "SCEE_SCALE_FACTOR",
    "SCNB_SCALE_FACTOR",
    "SOLTY",
    "LENNARD_JONES_ACOEF",
    "LENNARD_JONES_BCOEF",
    "BONDS_INC_HYDROGEN",
    "BONDS_WITHOUT_HYDROGEN",
    "ANGLES_INC_HYDROGEN",
    "ANGLES_WITHOUT_HYDROGEN",
    "DIHEDRALS_INC_HYDROGEN",
    "DIHEDRALS_WITHOUT_HYDROGEN",
    "EXCLUDED_ATOMS_LIST",
    "HBOND_ACOEF",
    "HBOND_BCOEF",
    "HBCUT",
    "AMBER_ATOM_TYPE",
    "TREE_CHAIN_CLASSIFICATION",
    "JOIN_ARRAY",
    "IROTAT",
    "RADIUS_SET",
    "RADII",
    "SCREEN"
};

bool compare_sections(const AmberTopSection *s1, const AmberTopSection *s2) {
    static int array_size = ARRAY_SIZE(kSectionList);
    return std::find(kSectionList, kSectionList + array_size, s1->name()) <
           std::find(kSectionList, kSectionList + array_size, s2->name());
}

}  // namespace

typedef Structure::AtomList AtomList;

const char *AmberTopBuilder::kDefaultTitle =
        "GENERATED BY GMML (http://code.google.com/p/gmml)";

AmberTopBuilder::AmberTopBuilder()
        : parameter_file_set_(*kDefaultEnvironment.parm_set()) {}

AmberTopBuilder::AmberTopBuilder(const ParameterSet& parameter_file_set)
        : parameter_file_set_(parameter_file_set) {}

AmberTopFile *AmberTopBuilder::build(const Structure& structure,
                                     const string& title) const {
    AmberTopFile *file = build_common_sections(structure, title);
    file->sort(&compare_sections);
    return file;
}

AmberTopFile *AmberTopBuilder::build(const BoxedStructure& boxed_structure,
                                     const string& title) const {
    AmberTopFile *file = build_common_sections(boxed_structure, title);
    build_box_section(boxed_structure, file);
    file->sort(&compare_sections);
    return file;
}

AmberTopFile *AmberTopBuilder::build(const SolvatedStructure& structure,
                                     const string& title) const {
    AmberTopFile *file = build_common_sections(structure, title);
    build_box_section(structure, file);
    build_solvation_sections(structure, file);

    file->sort(&compare_sections);

    return file;
}

AmberTopFile *AmberTopBuilder::build_common_sections(
        const Structure& structure, const string& title) const {
    AmberTopFile *file = new AmberTopFile;

    AmberTopStringSection *title_section =
            file->create_string_section("TITLE", "20a4");
    title_section->insert(title);

    AmberTopIntSection *pointer_section =
            file->create_int_section("POINTERS", "10I8");

    for (int i = 0; i < kPointerSectionSize; i++) {
        pointer_section->insert(0);
    }

    build_residues(structure, file);
    build_atoms(structure, file);
    build_type_info(structure, file);
    build_excluded_atoms(structure, file);
    build_bonds(structure, file);
    build_angles(structure, file);
    build_dihedrals(structure, file);
    build_radii_and_screen(structure, file);
    build_garbage_sections(structure.size(), file);

    int ep_atom_count = 0;
    for (int i = 0; i < structure.size(); i++) {
        if (structure.atoms(i)->type() == "EP")
            ep_atom_count++;
    }
    pointer_section->set(30, ep_atom_count);

    return file;
}

void AmberTopBuilder::build_residues(const Structure& structure,
                                     AmberTopFile *file) const {
    AmberTopStringSection *labels =
            file->create_string_section("RESIDUE_LABEL", "20a4");
    AmberTopIntSection *pointers =
            file->create_int_section("RESIDUE_POINTER", "10I8");

    int count = 1;
    int largest_residue_size = 0;
    for (size_t i = 0; i < structure.residue_count(); i++) {
        labels->insert(structure.residues(i)->name());
        pointers->insert(count);
        int residue_size = structure.residues(i)->size();
        count += residue_size;
        if (residue_size > largest_residue_size)
            largest_residue_size = residue_size;
    }

    AmberTopIntSection *pointer_section = file->get_int_section("POINTERS");
    pointer_section->set(11, structure.residue_count());
    pointer_section->set(28, largest_residue_size);
}

void AmberTopBuilder::build_atoms(const Structure& structure,
                                  AmberTopFile *file) const {
    AmberTopStringSection *types =
            file->create_string_section("AMBER_ATOM_TYPE", "20a4");
    AmberTopStringSection *names =
            file->create_string_section("ATOM_NAME", "20a4");
    AmberTopDoubleSection *charges =
            file->create_double_section("CHARGE", "5E16.8");
    AmberTopDoubleSection *masses =
            file->create_double_section("MASS", "5E16.8");

    for (int i = 0; i < structure.size(); i++) {
        const Atom *atom = structure.atoms(i);
        std::string type = atom->type();
        const ParameterFileAtom *parameter_atom = 
                parameter_file_set_.lookup(type);
        if (parameter_atom == NULL)
            type_error(type);
        types->insert(type);
        names->insert(atom->name());
        charges->insert(atom->charge()*kChargeFactor);
        masses->insert(parameter_atom->mass);
    }

    file->get_int_section("POINTERS")->set(0, structure.size());
}

void AmberTopBuilder::build_type_info(const Structure& structure,
                                      AmberTopFile *file) const {
    AmberTopIntSection *indices =
            file->create_int_section("ATOM_TYPE_INDEX", "10I8");

    vector<string> types;
    for (Structure::const_iterator it = structure.begin();
            it != structure.end(); ++it) {
        string type = (*it)->type();
        vector<string>::iterator it2 =
            std::find(types.begin(), types.end(), type);
        if (it2 == types.end()) {
            types.push_back(type);
            indices->insert(static_cast<int>(types.size()));
        } else {
            int index = std::distance(types.begin(), it2) + 1;
            indices->insert(index);
        }
    }

    build_nonbonded(file, types);

    file->get_int_section("POINTERS")->set(1, types.size());
}

void AmberTopBuilder::build_nonbonded(AmberTopFile *file,
                                      const vector<string>& type_names) const {
    vector<TypeInfoType*> type_info_types;
    type_info_types.reserve(type_names.size());
    vector<string>::const_iterator it;
    for (it = type_names.begin(); it != type_names.end(); ++it) {
        const ParameterFileAtom *parameter_atom =
            parameter_file_set_.lookup(*it);
        if (parameter_atom == NULL)
            type_error(*it);
        type_info_types.push_back(new TypeInfoType(parameter_atom->radius,
                                                   parameter_atom->well_depth));
    }

    AmberTopDoubleSection *a_section =
            file->create_double_section("LENNARD_JONES_ACOEF", "5E16.8");
    AmberTopDoubleSection *b_section =
            file->create_double_section("LENNARD_JONES_BCOEF", "5E16.8");

    int parm_index = 1;
    vector<vector<int> > parm_matrix(type_info_types.size());
    for (int i = 0; i < parm_matrix.size(); i++) {
        parm_matrix[i].resize(type_info_types.size());
        for (int j = 0; j <= i; j++) {
            parm_matrix[i][j] = parm_index;
            parm_matrix[j][i] = parm_index;
            parm_index++;
            std::pair<double, double> lj_coefficients =
                get_lennard_jones_coefficients(type_info_types[i],
                                               type_info_types[j]);
            a_section->insert(lj_coefficients.first);
            b_section->insert(lj_coefficients.second);
        }
    }

    std::for_each(type_info_types.begin(), type_info_types.end(), DeletePtr());

    AmberTopIntSection *nonbonded_indices =
            file->create_int_section("NONBONDED_PARM_INDEX", "10I8");
    for (int i = 0; i < parm_matrix.size(); i++)
        for (int j = 0; j < parm_matrix[i].size(); j++)
            nonbonded_indices->insert(parm_matrix[i][j]);
}

void AmberTopBuilder::build_excluded_atoms(const Structure& structure,
                                           AmberTopFile *file) const {
    AmberTopIntSection *count_section =
            file->create_int_section("NUMBER_EXCLUDED_ATOMS", "10I8");
    AmberTopIntSection *list_section =
            file->create_int_section("EXCLUDED_ATOMS_LIST", "10I8");

    int num_atoms = structure.size();
    for (int i = 0; i < num_atoms; i++) {
        vector<size_t> *atoms = structure.bonds()->bfs(i, 3);

        atoms->erase(std::remove_if(atoms->begin(), atoms->end(),
                                    std::bind2nd(std::less_equal<int>(), i)),
                     atoms->end());
        std::sort(atoms->begin(), atoms->end());

        if (atoms->empty()) {
            // It appears that there is no good reason for these values.
            // See http://archive.ambermd.org/200611/0140.html.
            count_section->insert(1);
            list_section->insert(0);
        } else {
            count_section->insert(atoms->size());
            for (int j = 0; j < atoms->size(); j++)
                list_section->insert((*atoms)[j] + 1);
        }
        delete atoms;
    }

    file->get_int_section("POINTERS")->set(10, count_section->sum());
}

void AmberTopBuilder::build_bonds(const Structure& structure,
                                  AmberTopFile *file) const {
    AmberTopIntSection *bonds_with_hydrogen =
            file->create_int_section("BONDS_INC_HYDROGEN", "10I8");
    AmberTopIntSection *bonds_without_hydrogen =
            file->create_int_section("BONDS_WITHOUT_HYDROGEN", "10I8");

    vector<BondType*> bond_types;
    for (int i = 0; i < structure.size(); i++) {
        const vector<size_t>& adj_atoms = structure.bonds(i);
        for (int j = 0; j < adj_atoms.size(); j++) {
            if (adj_atoms[j] < i)
                continue;
            const Atom *atom1 = structure.atoms(i);
            const Atom *atom2 = structure.atoms(adj_atoms[j]);
            int type_index = get_bond_type_index(bond_types, atom1->type(),
                                                 atom2->type());
            AmberTopIntSection *section;
            Element hydrogen("H");
            if (atom1->element() == hydrogen || atom2->element() == hydrogen)
                section = bonds_with_hydrogen;
            else
                section = bonds_without_hydrogen;
            section->insert(i*3);
            section->insert(static_cast<int>(adj_atoms[j])*3);
            section->insert(type_index);
        }
    }

    AmberTopDoubleSection *force_constants =
            file->create_double_section("BOND_FORCE_CONSTANT", "5E16.8");
    AmberTopDoubleSection *equil_values =
            file->create_double_section("BOND_EQUIL_VALUE", "5E16.8");
    for (vector<BondType*>::iterator it = bond_types.begin();
            it != bond_types.end(); ++it) {
        force_constants->insert((*it)->force_constant);
        equil_values->insert((*it)->equil_value);
    }

    std::for_each(bond_types.begin(), bond_types.end(), DeletePtr());

    AmberTopIntSection *pointer_section = file->get_int_section("POINTERS");
    pointer_section->set(2, bonds_with_hydrogen->size()/3);
    pointer_section->set(3, bonds_without_hydrogen->size()/3);
    pointer_section->set(12, pointer_section->get(3));
    pointer_section->set(15, bond_types.size());
}

void AmberTopBuilder::build_angles(const Structure& structure,
                                   AmberTopFile *file) const {
    AmberTopIntSection *angles_with_hydrogen =
            file->create_int_section("ANGLES_INC_HYDROGEN", "10I8");
    AmberTopIntSection *angles_without_hydrogen =
            file->create_int_section("ANGLES_WITHOUT_HYDROGEN", "10I8");

    const Graph *bonds = structure.bonds();
    vector<AngleType*> angle_types;
    for (int i = 0; i < structure.size(); i++) {
        const vector<size_t>& adj_atoms = structure.bonds(i);
        for (int j = 0; j < adj_atoms.size(); j++) {
            for (int k = j + 1; k < adj_atoms.size(); k++) {
                int adj_atom1 = adj_atoms[j];
                int adj_atom2 = adj_atoms[k];
                // This ensures we don't process the same angle twice.
                if (adj_atom2 < adj_atom1)
                    continue;
                // This is a check to make sure we're not looking at a triangle.
                if (std::find(bonds->edges(adj_atom1).begin(),
                              bonds->edges(adj_atom1).end(),
                              adj_atom2) !=
                        bonds->edges(adj_atom1).end()) {
                    continue;
                }
                const Atom *atom1 = structure.atoms(adj_atom1);
                const Atom *atom2 = structure.atoms(i);
                const Atom *atom3 = structure.atoms(adj_atom2);
                if (atom1->type() == "EP" || atom2->type() == "EP" ||
                        atom3->type() == "EP")
                    continue;
                int type_index = get_angle_type_index(angle_types,
                                                      atom1->type(),
                                                      atom2->type(),
                                                      atom3->type());
                AmberTopIntSection *section;
                Element hydrogen("H");
                if (atom1->element() == hydrogen ||
                        atom2->element() == hydrogen ||
                        atom3->element() == hydrogen) {
                    section = angles_with_hydrogen;
                } else {
                    section = angles_without_hydrogen;
                }

                section->insert(adj_atom1*3);
                section->insert(i*3);
                section->insert(adj_atom2*3);
                section->insert(type_index);
            }
        }
    }

    AmberTopDoubleSection *force_constants =
        file->create_double_section("ANGLE_FORCE_CONSTANT", "5E16.8");
    AmberTopDoubleSection *equil_values =
        file->create_double_section("ANGLE_EQUIL_VALUE", "5E16.8");
    for (vector<AngleType*>::iterator it = angle_types.begin();
            it != angle_types.end(); ++it) {
        force_constants->insert((*it)->force_constant);
        equil_values->insert(to_radians((*it)->equil_value));
    }

    std::for_each(angle_types.begin(), angle_types.end(), DeletePtr());

    AmberTopIntSection *pointer_section = file->get_int_section("POINTERS");
    pointer_section->set(4, angles_with_hydrogen->size()/4);
    pointer_section->set(5, angles_without_hydrogen->size()/4);
    pointer_section->set(13, pointer_section->get(5));
    pointer_section->set(16, angle_types.size());
}

void AmberTopBuilder::build_dihedrals(const Structure& structure,
                                      AmberTopFile *file) const {
    AmberTopIntSection *dihedrals_with_hydrogen =
            file->create_int_section("DIHEDRALS_INC_HYDROGEN", "10I8");
    AmberTopIntSection *dihedrals_without_hydrogen =
            file->create_int_section("DIHEDRALS_WITHOUT_HYDROGEN", "10I8");

    vector<DihedralType*> dihedral_types;

    for (int i = 0; i < structure.size(); i++) {
        const Structure::AdjList& adj_atoms = structure.bonds(i);
        for (int j = 0; j < adj_atoms.size(); j++) {
            if (i < adj_atoms[j])
                insert_dihedrals(structure, i, adj_atoms[j],
                                 dihedral_types,
                                 dihedrals_with_hydrogen,
                                 dihedrals_without_hydrogen);
        }
    }

    for (int i = 0; i < structure.size(); i++) {
        insert_improper_dihedrals(structure, i, dihedral_types,
                                  dihedrals_with_hydrogen,
                                  dihedrals_without_hydrogen);
    }

    AmberTopDoubleSection *force_constants =
            file->create_double_section("DIHEDRAL_FORCE_CONSTANT", "5E16.8");
    AmberTopDoubleSection *periodicities =
            file->create_double_section("DIHEDRAL_PERIODICITY", "5E16.8");
    AmberTopDoubleSection *phases =
            file->create_double_section("DIHEDRAL_PHASE", "5E16.8");
    AmberTopDoubleSection *scees =
            file->create_double_section("SCEE_SCALE_FACTOR", "5E16.8");
    AmberTopDoubleSection *scnbs =
            file->create_double_section("SCNB_SCALE_FACTOR", "5E16.8");

    for (vector<DihedralType*>::iterator it = dihedral_types.begin();
            it != dihedral_types.end(); ++it) {
        force_constants->insert((*it)->force_constant);
        periodicities->insert(fabs((*it)->periodicity));
        phases->insert(to_radians((*it)->phase));
        scees->insert((*it)->scee);
        scnbs->insert((*it)->scnb);
    }

    AmberTopIntSection *pointer_section = file->get_int_section("POINTERS");
    pointer_section->set(6, dihedrals_with_hydrogen->size()/5);
    pointer_section->set(7, dihedrals_without_hydrogen->size()/5);
    pointer_section->set(14, pointer_section->get(7));
    pointer_section->set(17, dihedral_types.size());

    std::for_each(dihedral_types.begin(), dihedral_types.end(), DeletePtr());
}

void AmberTopBuilder::build_radii_and_screen(const Structure& structure,
                                             AmberTopFile *file) const {
    AmberTopStringSection *radius_set =
            file->create_string_section("RADIUS_SET", "1a80");
    string radius_name = "modified Bondi radii (mbondi)";
    radius_set->insert(radius_name);

    AmberTopDoubleSection *radii_section =
            file->create_double_section("RADII", "5E16.8");
    AmberTopDoubleSection *screen_section =
            file->create_double_section("SCREEN", "5E16.8");
    const Structure::AtomList& atoms = structure.atoms();
    for (int i = 0; i < atoms.size(); i++) {
        std::pair<double, double> radius_and_screen =
            get_radius_and_screen(structure, i);
        radii_section->insert(radius_and_screen.first);
        screen_section->insert(radius_and_screen.second);
    }
}

void AmberTopBuilder::build_garbage_sections(int atom_count,
                                             AmberTopFile *file) const {
    file->create_double_section("SOLTY", "5E16.8");
    file->create_double_section("HBOND_ACOEF", "5E16.8");
    file->create_double_section("HBOND_BCOEF", "5E16.8");
    file->create_double_section("HBCUT", "5E16.8");
    AmberTopStringSection *chain_section =
            file->create_string_section("TREE_CHAIN_CLASSIFICATION", "20a4");
    AmberTopIntSection *join_array =
            file->create_int_section("JOIN_ARRAY", "10I8");
    AmberTopIntSection *i_rotat = file->create_int_section("IROTAT", "10I8");
    for (int i = 0; i < atom_count; i++) {
        chain_section->insert("M");
        join_array->insert(0);
        i_rotat->insert(0);
    }
}

void AmberTopBuilder::build_box_section(const BoxedStructure& boxed_structure,
                                        AmberTopFile *file) const {
    const Box *box = boxed_structure.box();
    if (box == NULL)
        return;
    AmberTopDoubleSection *box_section =
            file->create_double_section("BOX_DIMENSIONS", "5E16.8");
    box_section->insert(box->angle);
    box_section->insert(box->length);
    box_section->insert(box->width);
    box_section->insert(box->height);

    file->get_int_section("POINTERS")->set(27, 1);
}

void AmberTopBuilder::build_solvation_sections(
        const SolvatedStructure& structure, AmberTopFile *file) const {
    AmberTopIntSection *atoms_per_molecule =
            file->create_int_section("ATOMS_PER_MOLECULE", "10I8");
    int last_solute_atom = structure.last_solute_atom();
    int first_solvent_atom = last_solute_atom + 1;

    int first_solvent_molecule = kNotSet;
    int num_molecules = 0;
    deque<bool> marked(structure.size(), false);
    for (int i = 0; i < marked.size(); i++) {
        if (!marked[i]) {
            num_molecules++;
            vector<size_t> *component = structure.bonds()->bfs(i);
            atoms_per_molecule->insert(static_cast<int>(component->size()));
            for (int j = 0; j < component->size(); j++) {
                int atom = (*component)[j];
                marked[atom] = true;
                if (atom == first_solvent_atom)
                    first_solvent_molecule = num_molecules;
            }
            delete component;
        }
    }

    AmberTopIntSection *solvent_pointers =
            file->create_int_section("SOLVENT_POINTERS", "3I8");
    int last_solute_residue = structure.get_residue_index(last_solute_atom);
    solvent_pointers->insert(last_solute_residue + 1);
    solvent_pointers->insert(num_molecules);
    solvent_pointers->insert(first_solvent_molecule);
}

int AmberTopBuilder::get_bond_type_index(vector<BondType*>& bonds,
                                         const string& type1,
                                         const string& type2) const {
    const ParameterFileBond *parameter_bond =
        parameter_file_set_.lookup(type1, type2);
    if (parameter_bond == NULL || !is_set(parameter_bond->force_constant) ||
            !is_set(parameter_bond->length)) {
        type_error(type1, type2);
        return -1;
    }

    BondType *type = new BondType(parameter_bond->force_constant,
                                  parameter_bond->length);
    vector<BondType*>::iterator it = bonds.begin();
    while (it != bonds.end()) {
        if (*type == (**it)) {
            delete type;
            break;
        }
        ++it;
    }
    if (it == bonds.end()) {
        bonds.push_back(type);
        return bonds.size();
    }
    return std::distance(bonds.begin(), it) + 1;
}

int AmberTopBuilder::get_angle_type_index(vector<AngleType*>& angles,
                                          const string& type1,
                                          const string& type2,
                                          const string& type3) const {
    const ParameterFileAngle *parameter_angle =
        parameter_file_set_.lookup(type1, type2, type3);
    if (parameter_angle == NULL || !is_set(parameter_angle->force_constant) ||
            !is_set(parameter_angle->angle)) {
        type_error(type1, type2, type3);
        return -1;
    }

    AngleType *type = new AngleType(parameter_angle->force_constant,
                                    parameter_angle->angle);
    vector<AngleType*>::iterator it = angles.begin();
    while (it != angles.end()) {
        if (*type == (**it)) {
            delete type;
            break;
        }
        ++it;
    }
    if (it == angles.end()) {
        angles.push_back(type);
        return angles.size();
    }
    return std::distance(angles.begin(), it) + 1;
}

void AmberTopBuilder::insert_dihedrals(
        const Structure& structure,
        int atom1_index, int atom2_index,
        vector<DihedralType*>& dihedral_types,
        AmberTopIntSection *with_hydrogen,
        AmberTopIntSection *without_hydrogen) const {
    const Structure::AdjList& adj_atoms1 = structure.bonds(atom1_index);
    const Structure::AdjList& adj_atoms2 = structure.bonds(atom2_index);

    for (int i = 0; i < adj_atoms1.size(); i++) {
        for (int j = 0; j < adj_atoms2.size(); j++) {
            if (adj_atoms1[i] == adj_atoms2[j] ||
                    adj_atoms1[i] == atom2_index ||
                    adj_atoms2[j] == atom1_index) {
                continue;
            }

            const Atom *atom1 = structure.atoms(adj_atoms1[i]);
            const Atom *atom2 = structure.atoms(atom1_index);
            const Atom *atom3 = structure.atoms(atom2_index);
            const Atom *atom4 = structure.atoms(adj_atoms2[j]);

            if (atom1->type() == "EP" || atom2->type() == "EP" ||
                    atom3->type() == "EP" || atom4->type() == "EP")
                continue;

            const ParameterFileDihedral *parameter_dihedral =
                parameter_file_set_.lookup(atom1->type(), atom2->type(),
                                           atom3->type(), atom4->type());
            if (parameter_dihedral == NULL)
                type_error(atom1->type(), atom2->type(), atom3->type(),
                           atom4->type());
            const vector<ParameterFileDihedralTerm>& terms =
                parameter_dihedral->terms;
            for (int k = 0; k < terms.size(); k++) {
                const ParameterFileDihedralTerm& term = terms[k];
                DihedralType *type =
                    new DihedralType(term.force_constant, term.periodicity,
                                     term.phase, parameter_dihedral->scee,
                                     parameter_dihedral->scnb);
                if (!is_set(type->scee))
                    type->scee = kDefaultScee;
                if (!is_set(type->scnb))
                    type->scnb = kDefaultScnb;

                int type_index = get_dihedral_type_index(dihedral_types, type);
                AmberTopIntSection *section;
                Element hydrogen("H");
                if (atom1->element() == hydrogen ||
                        atom2->element() == hydrogen ||
                        atom3->element() == hydrogen ||
                        atom4->element() == hydrogen) {
                    section = with_hydrogen;
                } else {
                    section = without_hydrogen;
                }

                section->insert(static_cast<int>(adj_atoms1[i])*3);
                section->insert(atom1_index*3);
                section->insert((k != 0)?atom2_index*-3:atom2_index*3);
                section->insert(static_cast<int>(adj_atoms2[j])*3);
                section->insert(type_index);
            }
        }
    }
}

void AmberTopBuilder::insert_improper_dihedrals(
        const Structure& structure,
        int atom_index,
        vector<DihedralType*>& dihedral_types,
        AmberTopIntSection *dihedrals_with_hydrogen,
        AmberTopIntSection *dihedrals_without_hydrogen) const {
    const Structure::AdjList& adj_atoms = structure.bonds(atom_index);

    const Atom *center_atom = structure.atoms(atom_index);
    for (int i = 0; i < adj_atoms.size(); i++) {
        const Atom *adj_atom1 = structure.atoms(adj_atoms[i]);
        for (int j = i + 1; j < adj_atoms.size(); j++) {
            const Atom *adj_atom2 = structure.atoms(adj_atoms[j]);
            for (int k = j + 1; k < adj_atoms.size(); k++) {
                const Atom *adj_atom3 = structure.atoms(adj_atoms[k]);
                std::pair<ParameterFileDihedralTerm, bool> parameter_term =
                    parameter_file_set_.lookup_improper_dihedral(
                            center_atom->type(),
                            adj_atom1->type(),
                            adj_atom2->type(),
                            adj_atom3->type());
                if (!parameter_term.second)
                    continue;
                const ParameterFileDihedralTerm& term = parameter_term.first;
                DihedralType *type =
                    new DihedralType(term.force_constant, term.periodicity,
                                     term.phase, 0.0, 0.0);

                int type_index = get_dihedral_type_index(dihedral_types, type);
                AmberTopIntSection *section;
                Element hydrogen("H");
                if (center_atom->element() == hydrogen ||
                        adj_atom1->element() == hydrogen ||
                        adj_atom2->element() == hydrogen ||
                        adj_atom3->element() == hydrogen) {
                    section = dihedrals_with_hydrogen;
                } else {
                    section = dihedrals_without_hydrogen;
                }

                section->insert(static_cast<int>(adj_atoms[i])*3);
                section->insert(static_cast<int>(adj_atoms[j])*3);
                section->insert(atom_index*-3);
                section->insert(static_cast<int>(adj_atoms[k])*-3);
                section->insert(type_index);
            }
        }
    }
}

int AmberTopBuilder::get_dihedral_type_index(vector<DihedralType*>& types,
                                             DihedralType *type) const {
    vector<DihedralType*>::iterator it = types.begin();
    while (it != types.end()) {
        if (*type == (**it)) {
            delete type;
            break;
        }
        ++it;
    }
    if (it == types.end()) {
        types.push_back(type);
        return types.size();
    }
    return std::distance(types.begin(), it) + 1;
}

std::pair<double, double> AmberTopBuilder::get_lennard_jones_coefficients(
        TypeInfoType *type1,
        TypeInfoType *type2) const {
    double w = sqrt(type1->well_depth*type2->well_depth);
    double r6 = pow(type1->radius + type2->radius, 6.0);
    return std::make_pair(w*r6*r6, 2*w*r6);
}

std::pair<double, double> AmberTopBuilder::get_radius_and_screen(
        const Structure& structure,
        int atom_index) const {
    const Atom *atom = structure.atoms(atom_index);
    const Element& element = atom->element();
    double radius = kNotSet;
    double screen = kNotSet;
    if (element == Element("H")) {
        double screen = 0.85;
        const Structure::AdjList& adj_atoms = structure.bonds(atom_index);
        if (adj_atoms.empty()) {
            radius = 1.2;
        } else {
            const Element& adj_el = structure.atoms(adj_atoms[0])->element();
            if (adj_el == Element("C") || adj_el == Element("N")) {
                radius = 1.3;
            } else if (adj_el == Element("H") || adj_el == Element("O") ||
                    adj_el == Element("S")) {
                    radius = 0.8;
            }
        }
        return std::make_pair(radius, screen);
    }
    if (element == Element("C")) {
        string type = atom->type();
        if (type == "C1" || type == "C2" || type == "C3")
            radius = 2.2;
        else
            radius = 1.7;
        screen = 0.72;
    } else if (element == Element("N")) {
        radius = 1.55;
        screen = 0.79;
    } else if (element == Element("O")) {
        radius = 1.5;
        screen = 0.85;
    } else if (element == Element("F")) {
        radius = 1.5;
    } else if (element == Element("Si")) {
        radius = 2.1;
    } else if (element == Element("P")) {
        radius = 1.85;
        screen = 0.86;
    } else if (element == Element("S")) {
        radius = 1.8;
        screen = 0.96;
    } else if (element == Element("Cl")) {
        radius = 1.7;
    }

    if (!is_set(radius))
        radius = 1.5;
    if (!is_set(screen))
        screen = 0.8;

    return std::make_pair(radius, screen);
}

class BuildTopologyFile {
  public:
    explicit BuildTopologyFile(const AmberTopFile& file) : file_(file) {
    }

    Structure *operator()() {
        structure_ = new Structure;

        add_residues();
        add_bonds();

        return structure_;
    }

  private:
    void add_residues() {
        const AmberTopStringSection *residue_labels =
                file_.get_string_section("RESIDUE_LABEL");
        const AmberTopStringSection *atom_names =
                file_.get_string_section("ATOM_NAME");
        const AmberTopStringSection *type_names =
                file_.get_string_section("AMBER_ATOM_TYPE");
        const AmberTopIntSection *type_indices =
                file_.get_int_section("ATOM_TYPE_INDEX");
        const AmberTopDoubleSection *charges =
                file_.get_double_section("CHARGE");

        vector<int> *residue_sizes = get_residue_sizes();
        int cur_atom = 0;
        for (int i = 0; i < residue_sizes->size(); i++) {
            Residue residue(residue_labels->get(i));
            for (int j = 0; j < (*residue_sizes)[i]; j++) {
                Atom atom(guess_element_by_name(atom_names->get(cur_atom)),
                          Coordinate(0.0, 0.0, 0.0),
                          atom_names->get(cur_atom),
                          type_names->get(cur_atom),
                          charges->get(cur_atom));

                cur_atom++;
                residue.append(&atom);
            }
            structure_->append(&residue, true);
        }
        delete residue_sizes;
    }

    vector<int> *get_residue_sizes() const {
        const AmberTopIntSection *residue_pointers =
                file_.get_int_section("RESIDUE_POINTER");

        vector<int> *sizes = new vector<int>;

        for (int i = 1; i < residue_pointers->size(); i++) {
            sizes->push_back(residue_pointers->get(i) -
                             residue_pointers->get(i - 1));
        }

        int last = file_.get_string_section("ATOM_NAME")->size() -
                   residue_pointers->get(residue_pointers->size() - 1) + 1;

        sizes->push_back(last);

        return sizes;
    }

    void add_bonds() {
        add_bond_section(file_.get_int_section("BONDS_INC_HYDROGEN"));
        add_bond_section(file_.get_int_section("BONDS_WITHOUT_HYDROGEN"));
    }

    void add_bond_section(const AmberTopIntSection *bond_section) {
        for (int i = 0; i < bond_section->size(); i += 3) {
            structure_->add_bond(bond_section->get(i)/3,
                                 bond_section->get(i + 1)/3);
        }
    }

    Structure *structure_;
    const AmberTopFile& file_;
};

Structure *build_topology_file(const AmberTopFile& file) {
    return BuildTopologyFile(file)();
}

}  // namespace gmml
