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
//
// This file contains classes to build the sections of an AMBER topology file
// for a given structure.

#ifndef GMML_INTERNAL_AMBER_TOP_BUILDER_H_
#define GMML_INTERNAL_AMBER_TOP_BUILDER_H_

#include <exception>
#include <string>
#include <utility>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "gmml/internal/stubs/common.h"

namespace gmml {

class AmberTopFile;
class AmberTopIntSection;
class AmberTopSection;
class BoxedStructure;
class ParameterSet;
class SolvatedStructure;
class Structure;

namespace amber_top_builder {

struct TypeInfoType;
struct BondType;
struct AngleType;
struct DihedralType;

}

// This class creates the sections of an AMBER topology file. The file
// specification can be found here: http://ambermd.org/formats.html#topology.
class AmberTopBuilder {
  public:
    typedef amber_top_builder::TypeInfoType TypeInfoType;
    typedef amber_top_builder::BondType BondType;
    typedef amber_top_builder::AngleType AngleType;
    typedef amber_top_builder::DihedralType DihedralType;

    static const int kPointerSectionSize = 32;

    // Factor by which to adjust the atomic charges.
    // See http://ambermd.org/Questions/units.html for an explanation.
    static const double kChargeFactor = 18.2223;

    // These values are used if they aren't found in the parameter file set.
    static const double kDefaultScee = 1.2;
    static const double kDefaultScnb = 2.0;

    static const char *kDefaultTitle;

    // Almost all information contained in the topology file comes from
    // a set of parameter files.
    AmberTopBuilder();
    explicit AmberTopBuilder(const ParameterSet& parameter_file_set);

    // Builds the topology file for a structure, and optionally includes a title
    // for the TITLE section.
    AmberTopFile *build(const Structure& structure,
                        const std::string& title) const;
    AmberTopFile *build(const Structure& structure) const {
        return build(structure, kDefaultTitle);
    }

    // Boxed structures include a boxed section.
    AmberTopFile *build(const BoxedStructure& boxed_structure,
                        const std::string& title) const;
    AmberTopFile *build(const BoxedStructure& boxed_structure) const {
        return build(boxed_structure, kDefaultTitle);
    }

    // Solvated structures include some special sections. See below.
    AmberTopFile *build(const SolvatedStructure& solvated_structure,
                        const std::string& title) const;
    AmberTopFile *build(const SolvatedStructure& solvated_structure) const {
        return build(solvated_structure, kDefaultTitle);
    }

  private:
    // It may be useful to make everything prefixed with build_ part of the
    // public interface, if only for testing purposes.

    // This builds a file with sections that relate to all structures.
    AmberTopFile *build_common_sections(const Structure& structure,
                                        const std::string& title) const;

    // Builds the RESIDUE_LABEL and RESIDUE_POINTER sections. Also sets the
    // 12th and 29th pointers.
    void build_residues(const Structure&, AmberTopFile *file) const;

    // Builds the AMBER_ATOM_TYPE, ATOM_NAME, CHARGE, and MASS sections.
    // Also sets the 1st pointer.
    void build_atoms(const Structure&, AmberTopFile *file) const;

    // Builds the ATOM_TYPE_INDEX section and sets the 2nd pointers.
    void build_type_info(const Structure&, AmberTopFile *file) const;

    // Builds the LENNARD_JONES_ACOEF, LENNARD_JONES_BCOEF, and
    // NONBONDED_PARM_INDEX sections from a list of unique atom type names.
    void build_nonbonded(AmberTopFile *file,
                         const std::vector<std::string>& types) const;

    // Builds NUMBER_EXCLUDED_ATOMS and EXCLUDED_ATOMS_LIST sections and sets
    // the 11th pointer.
    void build_excluded_atoms(const Structure&, AmberTopFile *file) const;

    // Builds BONDS_INC_HYDROGEN, BONDS_WITHOUT_HYDROGEN, BOND_FORCE_CONSTANT,
    // and BOND_EQUIL_VALUE and sets pointers 3, 4, 13, and 16.
    void build_bonds(const Structure&, AmberTopFile *file) const;

    // Builds ANGLES_INC_HYDROGEN, ANGLES_WITHOUT_HYDROGEN,
    // ANGLE_FORCE_CONSTANT, ANGLE_EQUIL_VALUE and sets pointers 5, 6, 14,
    // and 17.
    void build_angles(const Structure&, AmberTopFile *file) const;

    // Builds DIHEDRALS_INC_HYDROGEN, DIHEDRALS_WITHOUT_HYDROGEN,
    // DIHEDRAL_FORCE_CONSTANT, DIHEDRAL_PERIODICITY, DIHEDRAL_PHASE,
    // SCEE_SCALE_FACTOR, and SCNB_SCALE_FACTOR and set pointers 7, 8, 15,
    // and 18.
    void build_dihedrals(const Structure&, AmberTopFile *file) const;

    // Builds RADIUS_SET, RADII, and SCREEN, where the RADIUS_SET is
    // "modified Bondi radii (mbondi)". The values used come from LEaP. This
    // is the only radius set we've needed, so it's not configurable. Let
    // us know if you would like to use a different set.
    void build_radii_and_screen(const Structure&, AmberTopFile *file) const;

    // Builds sections with dummy values to appease programs that expect to
    // find them. The section are SOLTY, HBOND_ACOEF, HBOND_BCOEF, HBCUT,
    // TREE_CHAIN_CLASSIFICATION, JOIN_ARRAY, and IROTAT.
    void build_garbage_sections(int atom_count, AmberTopFile *file) const;

    // This builds BOX_DIMENSIONS and sets pointer 28.
    void build_box_section(const BoxedStructure& structure,
                           AmberTopFile *file) const;

    // This builds ATOMS_PER_MOLECULE and SOLVENT_POINTERS.
    void build_solvation_sections(const SolvatedStructure& structure,
                                  AmberTopFile *file) const;

    int get_bond_type_index(std::vector<BondType*>&, const std::string& type1,
                            const std::string& type2) const;

    int get_angle_type_index(std::vector<AngleType*>&,
                             const std::string& type1,
                             const std::string& type2,
                             const std::string& type3) const;

    // A helper to insert all the dihedrals that exist where around a specific
    // bond (the bond from atom1_index to atom2_index).
    void insert_dihedrals(const Structure&, int atom1_index, int atom2_index,
                          std::vector<DihedralType*>& dihedral_types,
                          AmberTopIntSection *dihedrals_with_hydrogen,
                          AmberTopIntSection *dihedrals_without_hydrogen) const;

    // A helper to insert all improper dihedrals, where atom_index is the
    // center atom.
    void insert_improper_dihedrals(
            const Structure&, int atom_index,
            std::vector<DihedralType*>& dihedral_types,
            AmberTopIntSection *dihedrals_with_hydrogen,
            AmberTopIntSection *dihedrals_without_hydrogen) const;

    int get_dihedral_type_index(std::vector<DihedralType*>& dihedral_types,
                                DihedralType *type) const;

    // Returns the A and B coefficients for the type info types.
    std::pair<double, double> get_lennard_jones_coefficients(
            TypeInfoType *type1, TypeInfoType *type2) const;

    // The values returned are pulled from the source of LEaP, from a file
    // called unitio.c, I think.
    std::pair<double, double> get_radius_and_screen(const Structure& structure,
                                                    int atom_index) const;

    // These are used to format the error message. They are called when there
    // is insufficient information in the parameter file set for
    // atom types (1 argument), bonds (2 arguments), angles (3 arguments),
    // and dihedrals (4 arguments).
    void type_error(const std::string& type) const;
    void type_error(const std::string& type1, const std::string& type2) const;
    void type_error(const std::string& type1, const std::string& type2,
                    const std::string& type3) const;
    void type_error(const std::string& type1, const std::string& type2,
                    const std::string& type3, const std::string& type4) const;

    const ParameterSet& parameter_file_set_;

    DISALLOW_COPY_AND_ASSIGN(AmberTopBuilder);
};


// These structs are namespaced because of their general names. They are used
// to build sections of the topology file. An equality relation is defined on
// them to determine the correct "type" of bonds, angles, etc.
namespace amber_top_builder {

// These represent information for atom types and are used to to calculate
// Lennard-Jones coefficients.
struct TypeInfoType {
    TypeInfoType(double radius, double well_depth) : radius(radius),
                                                     well_depth(well_depth) {
    }

    double radius;
    double well_depth;
};

inline bool operator==(const TypeInfoType& lhs, const TypeInfoType& rhs) {
    return lhs.radius == rhs.radius &&
           lhs.well_depth == rhs.well_depth;
}

// These correspond to values in BOND_FORCE_CONSTANT and BOND_EQUIL_VALUE.
struct BondType {
    BondType(double force_constant, double equil_value)
            : force_constant(force_constant), equil_value(equil_value) {
    }

    double force_constant;
    double equil_value;
};

inline bool operator==(const BondType& lhs, const BondType& rhs) {
    return lhs.force_constant == rhs.force_constant &&
           lhs.equil_value == rhs.equil_value;
}

// These correspond to values in ANGLE_FORCE_CONSTANT and ANGLE_EQUIL_VALUE.
struct AngleType {
    AngleType(double force_constant, double equil_value)
            : force_constant(force_constant), equil_value(equil_value) {
    }

    double force_constant;
    double equil_value;
};

inline bool operator==(const AngleType& lhs, const AngleType& rhs) {
    return lhs.force_constant == rhs.force_constant &&
           lhs.equil_value == rhs.equil_value;
}

// These correspond to values in DIHEDRAL_FORCE_CONSTANT, DIHEDRAL_PERIODICITY,
// DIHEDRAL_PHASE, SCEE_SCALE_FACTOR, AND SCNB_SCALE_FACTOR.
struct DihedralType {
    DihedralType(double force_constant, double periodicity, double phase,
                 double scee, double scnb)
            : force_constant(force_constant), periodicity(periodicity),
              phase(phase), scee(scee), scnb(scnb) {
    }

    double force_constant;
    double periodicity;
    double phase;
    double scee;
    double scnb;
};

inline bool operator==(const DihedralType& lhs, const DihedralType& rhs) {
    return lhs.force_constant == rhs.force_constant &&
           lhs.periodicity == rhs.periodicity &&
           lhs.phase == rhs.phase &&
           lhs.scee == rhs.scee &&
           lhs.scnb == rhs.scnb;
}

}  // namespace amber_top_builder


class InsufficientParameterException : public std::exception {
  public:
    explicit InsufficientParameterException(const std::string& what) {
        what_ = what;
    }

    virtual const char *what() const throw() { return what_.c_str(); }

    virtual ~InsufficientParameterException() throw() {}

  private:
    std::string what_;
};

inline void AmberTopBuilder::type_error(const std::string& type) const {
    throw InsufficientParameterException(
            "Insufficient parameters for atom type " + type);
}

inline void AmberTopBuilder::type_error(const std::string& type1,
                                        const std::string& type2) const {
    throw InsufficientParameterException(
            "Insufficient parameters for bond " + type1 + "-" + type2);
}

inline void AmberTopBuilder::type_error(const std::string& type1,
                                        const std::string& type2,
                                        const std::string& type3) const {
    throw InsufficientParameterException(
            "Insufficient parameters for angle  " +
             type1 + "-" + type2 + "-" + type3);
}

inline void AmberTopBuilder::type_error(const std::string& type1,
                                        const std::string& type2,
                                        const std::string& type3,
                                        const std::string& type4) const {
    throw InsufficientParameterException(
            "Insufficient parameters for dihedral " +
            type1 + "-" + type2 + "-" + type3 + "-" + type4);
}

Structure *build_topology_file(const AmberTopFile& file);

}  // namespace gmml

#endif  // GMML_INTERNAL_AMBER_TOP_BUILDER_H_
