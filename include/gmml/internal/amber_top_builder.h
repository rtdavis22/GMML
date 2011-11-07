// Author: Robert Davis
//
// This file contains classes to build the sections of an AMBER topology file
// for a given structure.

#ifndef AMBER_TOP_BUILDER_H
#define AMBER_TOP_BUILDER_H

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include <boost/shared_ptr.hpp>

namespace gmml {

class AmberTopFile;
class AmberTopSection;
class ParameterFileSet;
class Structure;

// These structs are namespaced because of their general names. They are used
// to build sections of the topology file. An equality relation is defined on
// them to determine the correct "type" of bonds, angles, etc.
namespace amber_top_builder {

// These represent information for atom types and are used to to calculate
// Lennard-Jones coefficients.
struct TypeInfoType {
    TypeInfoType(double radius, double well_depth) : radius(radius),
                                                     well_depth(well_depth) {}

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
            : force_constant(force_constant), equil_value(equil_value) {}

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
            : force_constant(force_constant), equil_value(equil_value) {}

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
              phase(phase), scee(scee), scnb(scnb) {}

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

// This class creates the sections of an AMBER topology file. The file
// specification can be found here: http://ambermd.org/formats.html#topology.
class AmberTopBuilder {
  public:
    typedef boost::shared_ptr<AmberTopSection> SectionPtr;
    typedef amber_top_builder::TypeInfoType TypeInfoType;
    typedef amber_top_builder::BondType BondType;
    typedef amber_top_builder::AngleType AngleType;
    typedef amber_top_builder::DihedralType DihedralType;

    // Factor by which to adjust the atomic charges.
    // See http://ambermd.org/Questions/units.html for an explanation.
    static const double kChargeFactor = 18.2223;

    // These values are used if they aren't found in the parameter file set.
    static const double kDefaultScee = 1.2;
    static const double kDefaultScnb = 2.0;

    // Almost all information contained in the topology file comes from
    // a set of parameter files.
    AmberTopBuilder();
    explicit AmberTopBuilder(const ParameterFileSet& parameter_file_set);

    // Builds the topology file for a structure, and optionally includes a title
    // for the TITLE section.
    AmberTopFile *build(const Structure& structure,
                        const std::string& title) const;
    AmberTopFile *build(const Structure& structure) const {
        return build(structure, "");
    }

  private:
    // It may be useful to make everything prefixed with build_ part of the
    // public interface.

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
    // "modified Bondi radii (mbondi)"
    void build_radii_and_screen(const Structure&, AmberTopFile *file) const;

    // Builds sections with dummy values to appease programs that expect to
    // find them. The section are SOLTY, HBOND_ACOEF, HBOND_BCOEF, HBCUT,
    // TREE_CHAIN_CLASSIFICATION, JOIN_ARRAY, and IROTAT.
    void build_garbage_sections(int atom_count, AmberTopFile *file) const;

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
                          SectionPtr dihedrals_with_hydrogen,
                          SectionPtr dihedrals_without_hydrogen) const;

    // A help to insert all improper dihedrals, where atom_index is the
    // center atom.
    void insert_improper_dihedrals(const Structure&, int atom_index,
                                   std::vector<DihedralType*>& dihedral_types,
                                   SectionPtr dihedrals_with_hydrogen,
                                   SectionPtr dihedrals_without_hydrogen) const;

    int get_dihedral_type_index(std::vector<DihedralType*>& dihedral_types,
                                DihedralType *type) const;

    // Returns the A and B coefficients for the type info types.
    std::pair<double, double> get_lennard_jones_coefficients(
            TypeInfoType *type1, TypeInfoType *type2) const;

    // The values returned are pulled from the source of LEaP, from a file
    // called unitio.c, I think.
    std::pair<double, double> get_radius_and_screen(const Structure& structure,
                                                    int atom_index) const;

    // Invokes gmml::error
    void error(const std::string& message) const;

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

    const ParameterFileSet& parameter_file_set_;
};

inline void AmberTopBuilder::type_error(const std::string& type) const {
    error("Insufficient parameters for atom type " + type);
}

inline void AmberTopBuilder::type_error(const std::string& type1,
                                        const std::string& type2) const {
    error("Insufficient parameters for bond " + type1 + "-" + type2);
}

inline void AmberTopBuilder::type_error(const std::string& type1,
                                        const std::string& type2,
                                        const std::string& type3) const {
    error("Insufficient parameters for angle  " +
          type1 + "-" + type2 + "-" + type3);
}

inline void AmberTopBuilder::type_error(const std::string& type1,
                                        const std::string& type2,
                                        const std::string& type3,
                                        const std::string& type4) const {
    error("Insufficient parameters for dihedral " +
          type1 + "-" + type2 + "-" + type3 + "-" + type4);
}

// This class is populated with a list of topology file section names in the
// order they should appear in the file.
class SectionComparer {
  public:
    SectionComparer();

    // Each comparison takes O(n) time, so it may be preferable to use a
    // different data structure.
    bool operator()(const std::string& s1, const std::string& s2) const;

  private:
    // A smart pointer is used here because this functor will be passed to a
    // sort function and we want to avoid constructing the section list
    // multiple times.
    const boost::shared_ptr<std::vector<std::string> > section_list_;
};

}  // namespace gmml

#endif  // AMBER_TOP_BUILDER_H
