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

#include "gmml/internal/amber_top_file_comparer.h"

#include <iostream>
#include <set>
#include <vector>

#include "gmml/internal/amber_top_builder.h"
#include "gmml/internal/amber_top_file.h"
#include "gmml/internal/get_residue_mapping.h"
#include "gmml/internal/structure.h"

using std::cout;
using std::endl;
using std::set;
using std::vector;

namespace gmml {
namespace {

class CompareBonds {
  public:
    CompareBonds(const AmberTopFile& file1, const AmberTopFile& file2,
                 const vector<int>& mapping) : mapping_(mapping) {
        init(file1, file2);
    }

    void operator()() {
        for (set<const IndexedBondType*>::iterator it = file1_bonds_.begin();
                it != file1_bonds_.end(); ++it) {
            IndexedBondType bond_type(mapping_[(*it)->atom1_index],
                                      mapping_[(*it)->atom2_index],
                                      (*it)->type_index);
            search_and_remove_bond(&bond_type);
        }
        for (set<const IndexedBondType*>::iterator it = file2_bonds_.begin();
                it != file2_bonds_.end(); ++it) {
            cout << "Bond " << (*it)->atom1_index << "-" <<
                    (*it)->atom2_index << " not present in file 1." << endl;
        }
    }

    ~CompareBonds() {
        // Need to implement this.
    }

  private:
    struct IndexedBondType {
        IndexedBondType(int atom1_index, int atom2_index, int type_index)
                : atom1_index(atom1_index), atom2_index(atom2_index),
                  type_index(type_index) {
        }

        struct PtrLess {
            bool operator()(const IndexedBondType *lhs,
                            const IndexedBondType *rhs) {
                vector<int> lhs_types;
                load_bond_type_in_vector(&lhs_types, lhs);

                vector<int> rhs_types;
                load_bond_type_in_vector(&rhs_types, rhs);

                return lhs_types < rhs_types;
            }

          private:
            void load_bond_type_in_vector(vector<int> *vec,
                                          const IndexedBondType *type) {
                vec->push_back(type->atom1_index);
                vec->push_back(type->atom2_index);
                std::sort(vec->begin(), vec->end());
            }
        };

        int atom1_index;
        int atom2_index;
        int type_index;
    };

    void init(const AmberTopFile& file1, const AmberTopFile& file2) {
        init_bond_types(file1, &file1_bonds_, &file1_types_);
        init_bond_types(file2, &file2_bonds_, &file2_types_);
    }

    static void init_bond_types(const AmberTopFile& file,
                                set<const IndexedBondType*,
                                    IndexedBondType::PtrLess> *bonds,
                                vector<amber_top_builder::BondType*> *types) {
        init_bonds(bonds, file.get_int_section("BONDS_INC_HYDROGEN"));
        init_bonds(bonds, file.get_int_section("BONDS_WITHOUT_HYDROGEN"));

        const AmberTopDoubleSection *force_constants =
                file.get_double_section("BOND_FORCE_CONSTANT");
        const AmberTopDoubleSection *equil_values =
                file.get_double_section("BOND_EQUIL_VALUE");
        for (int i = 0; i < force_constants->size(); i++) {
            types->push_back(new amber_top_builder::BondType(
                    force_constants->get(i),
                    equil_values->get(i)));
        }
    }

    static void init_bonds(set<const IndexedBondType*,
                               IndexedBondType::PtrLess> *bonds,
                           const AmberTopIntSection *bond_section) {
        for (int i = 0; i < bond_section->size(); i += 3) {
            bonds->insert(new IndexedBondType(bond_section->get(i)/3,
                                              bond_section->get(i + 1)/3,
                                              bond_section->get(i + 2) - 1));
        }
    }

    void search_and_remove_bond(const IndexedBondType *bond_type) {
        set<const IndexedBondType*>::iterator it =
                file2_bonds_.find(bond_type);
        if (it == file2_bonds_.end()) {
            cout << "Bond" << bond_type->atom1_index << "-" <<
                    bond_type->atom2_index << " not found in file 2." << endl;
        } else {
            compare_bond_types(bond_type->type_index, (*it)->type_index);
            file2_bonds_.erase(it);
        }
    }

    void compare_bond_types(int file1_index, int file2_index) {
        amber_top_builder::BondType *type1 = file1_types_[file1_index];
        amber_top_builder::BondType *type2 = file2_types_[file2_index];
        if (type1->force_constant != type2->force_constant) {
            cout << "Different force constants for bond.";
            cout << " File 1 index: " << file1_index;
            cout << " File 2 index: " << file2_index << endl;
        }
        if (type1->equil_value != type2->equil_value) {
            cout << "Different equil values for bond.";
            cout << " File 1 index: " << file1_index;
            cout << " File 2 index: " << file2_index << endl;
        }
    }

    set<const IndexedBondType*, IndexedBondType::PtrLess> file1_bonds_;
    set<const IndexedBondType*, IndexedBondType::PtrLess> file2_bonds_;
    vector<amber_top_builder::BondType*> file1_types_;
    vector<amber_top_builder::BondType*> file2_types_;
    const vector<int>& mapping_;
};

class CompareAngles {
  public:
    CompareAngles(const AmberTopFile& file1, const AmberTopFile& file2,
                  const vector<int>& mapping) : mapping_(mapping) {
        init(file1, file2);
    }

    void operator()() {
        for (set<const IndexedAngleType*>::iterator it = file1_angles_.begin();
                it != file1_angles_.end(); ++it) {
            IndexedAngleType angle_type(mapping_[(*it)->atom1_index],
                                        mapping_[(*it)->atom2_index],
                                        mapping_[(*it)->atom3_index],
                                        (*it)->type_index);
            search_and_remove_angle(&angle_type);
        }
        for (set<const IndexedAngleType*>::iterator it = file2_angles_.begin();
                it != file2_angles_.end(); ++it) {
            cout << "Angle " << (*it)->atom1_index << "-" <<
                    (*it)->atom2_index << "-" << (*it)->atom3_index <<
                    " not present in file 1." << endl;
        }
    }

    ~CompareAngles() {
        // Need to implement this.
    }

  private:
    struct IndexedAngleType {
        IndexedAngleType(int atom1_index, int atom2_index, int atom3_index,
                         int type_index)
                : atom1_index(atom1_index), atom2_index(atom2_index),
                  atom3_index(atom3_index), type_index(type_index) {
        }

        struct PtrLess {
            bool operator()(const IndexedAngleType *lhs,
                            const IndexedAngleType *rhs) {
                vector<int> lhs_types;
                load_angle_type_in_vector(&lhs_types, lhs);
                vector<int> rhs_types;
                load_angle_type_in_vector(&rhs_types, rhs);
                return lhs_types < rhs_types;
            }

          private:
            void load_angle_type_in_vector(vector<int> *vec,
                                           const IndexedAngleType *type) {
                vec->push_back(type->atom1_index);
                vec->push_back(type->atom2_index);
                vec->push_back(type->atom3_index);
                std::sort(vec->begin(), vec->end());
            }
        };

        int atom1_index;
        int atom2_index;
        int atom3_index;
        int type_index;
    };

    void init(const AmberTopFile& file1, const AmberTopFile& file2) {
        init_angle_types(file1, &file1_angles_, &file1_types_);
        init_angle_types(file2, &file2_angles_, &file2_types_);
    }

    static void init_angle_types(const AmberTopFile& file,
                                 set<const IndexedAngleType*,
                                     IndexedAngleType::PtrLess> *angles,
                                 vector<amber_top_builder::AngleType*> *types) {
        init_angles(angles, file.get_int_section("ANGLES_INC_HYDROGEN"));
        init_angles(angles, file.get_int_section("ANGLES_WITHOUT_HYDROGEN"));

        const AmberTopDoubleSection *force_constants =
                file.get_double_section("ANGLE_FORCE_CONSTANT");
        const AmberTopDoubleSection *equil_values =
                file.get_double_section("ANGLE_EQUIL_VALUE");
        for (int i = 0; i < force_constants->size(); i++) {
            types->push_back(new amber_top_builder::AngleType(
                    force_constants->get(i),
                    equil_values->get(i)));
        }
    }

    static void init_angles(set<const IndexedAngleType*,
                                IndexedAngleType::PtrLess> *angles,
                           const AmberTopIntSection *angle_section) {
        for (int i = 0; i < angle_section->size(); i += 4) {
            angles->insert(new IndexedAngleType(angle_section->get(i)/3,
                                                angle_section->get(i + 1)/3,
                                                angle_section->get(i + 2)/3,
                                                angle_section->get(i + 3) - 1));
        }
    }

    void search_and_remove_angle(const IndexedAngleType *angle_type) {
        set<const IndexedAngleType*>::iterator it =
                file2_angles_.find(angle_type);
        if (it == file2_angles_.end()) {
            cout << "Angle" << angle_type->atom1_index << "-" <<
                    angle_type->atom2_index << "-" << angle_type->atom3_index <<
                    " not found in file 2." << endl;
        } else {
            compare_angle_types(angle_type->type_index, (*it)->type_index);
            file2_angles_.erase(it);
        }
    }

    void compare_angle_types(int file1_index, int file2_index) {
        amber_top_builder::AngleType *type1 = file1_types_[file1_index];
        amber_top_builder::AngleType *type2 = file2_types_[file2_index];
        if (type1->force_constant != type2->force_constant) {
            cout << "Different force constants for angle.";
            cout << " File 1 index: " << file1_index;
            cout << " File 2 index: " << file2_index << endl;
        }
        if (type1->equil_value != type2->equil_value) {
            cout << "Different equil values for angle.";
            cout << " File 1 index: " << file1_index;
            cout << " File 2 index: " << file2_index << endl;
        }
    }

    set<const IndexedAngleType*, IndexedAngleType::PtrLess> file1_angles_;
    set<const IndexedAngleType*, IndexedAngleType::PtrLess> file2_angles_;
    vector<amber_top_builder::AngleType*> file1_types_;
    vector<amber_top_builder::AngleType*> file2_types_;
    const vector<int>& mapping_;
};

}  // namespace

void AmberTopFileComparer::operator()() const {
    Structure *structure1 = build_topology_file(file1_);
    Structure *structure2 = build_topology_file(file2_);

    vector<int> *residue_mapping = get_residue_mapping(structure1, structure2);
    vector<int> *atom_mapping = get_atom_mapping(structure1, structure2,
                                                 *residue_mapping);

    CompareBonds(file1_, file2_, *atom_mapping)();
    CompareAngles(file1_, file2_, *atom_mapping)();

    delete atom_mapping;
    delete residue_mapping;
}

}  // namespace gmml
