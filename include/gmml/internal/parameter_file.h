// Author: Robert Davis
//
// TODO: pimplize this file.

#ifndef GMML_INTERNAL_PARAMETER_FILE_H_
#define GMML_INTERNAL_PARAMETER_FILE_H_

#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

struct ParameterFileAtom {
    ParameterFileAtom() : type(""), mass(kNotSet), polarizability(kNotSet),
                          radius(kNotSet), well_depth(well_depth) {}
    ParameterFileAtom(const std::string& type, double mass,
                      double polarizability, double radius, double well_depth)
            : type(type), mass(mass), polarizability(polarizability),
              radius(radius), well_depth(well_depth) {}
    ParameterFileAtom(const std::string& type, double mass,
                      double polarizability)
            : type(type), mass(mass), polarizability(polarizability),
              radius(kNotSet), well_depth(kNotSet) {}

    std::string type;
    double mass;
    double polarizability;
    double radius;
    double well_depth;
};

struct ParameterFileAtomPtrLess
        : public std::binary_function<const ParameterFileAtom*,
                                      const ParameterFileAtom*, bool> {
    bool operator()(const ParameterFileAtom *lhs,
                    const ParameterFileAtom *rhs) const {
        return lhs->type < rhs->type;
    }
};

struct ParameterFileBond {
    ParameterFileBond() : force_constant(kNotSet), length(kNotSet) {}
    ParameterFileBond(const std::string& type1, const std::string& type2)
            : force_constant(kNotSet), length(kNotSet) {
        types.reserve(2);
        types.push_back(type1);
        types.push_back(type2);
    }
    ParameterFileBond(const std::vector<std::string>& types,
                      double force_constant, double length)
            : types(types), force_constant(force_constant), length(length) {}

    std::vector<std::string> types;
    double force_constant;
    double length;
};

struct ParameterFileBondPtrLess
        : public std::binary_function<const ParameterFileBond*,
                                      const ParameterFileBond*, bool> {
    bool operator()(const ParameterFileBond *lhs,
                    const ParameterFileBond *rhs) const {
        return lhs->types < rhs->types;
    }
};

struct ParameterFileAngle {
    ParameterFileAngle() : force_constant(kNotSet), angle(kNotSet) {}
    ParameterFileAngle(const std::string& type1, const std::string& type2,
                       const std::string& type3)
            : force_constant(kNotSet), angle(kNotSet) {
        types.reserve(3);
        types.push_back(type1);
        types.push_back(type2);
        types.push_back(type3);
    }
    ParameterFileAngle(const std::vector<std::string>& types,
                       double force_constant, double angle)
            : types(types), force_constant(force_constant), angle(angle) {}

    std::vector<std::string> types;
    double force_constant;
    double angle;
};

struct ParameterFileAnglePtrLess
        : public std::binary_function<const ParameterFileAngle*,
                                      const ParameterFileAngle*, bool> {
    bool operator()(const ParameterFileAngle *lhs,
                    const ParameterFileAngle *rhs) const {
        return lhs->types < rhs->types;
    }
};

struct ParameterFileDihedralTerm {
    ParameterFileDihedralTerm() : factor(kNotSet), force_constant(kNotSet),
                                  phase(kNotSet), periodicity(kNotSet) {}

    double factor;
    double force_constant;
    double phase;
    double periodicity;
};

struct ParameterFileDihedral {
    ParameterFileDihedral() : scee(kNotSet), scnb(kNotSet) {}
    ParameterFileDihedral(const std::string& type1, const std::string& type2,
                          const std::string& type3, const std::string& type4)
            : scee(kNotSet), scnb(kNotSet) {
        types.reserve(4);
        types.push_back(type1);
        types.push_back(type2);
        types.push_back(type3);
        types.push_back(type4);
    }
    ParameterFileDihedral(std::vector<std::string> types,
                          const ParameterFileDihedralTerm& initial_term,
                          double scee, double scnb)
            : types(types), scee(scee), scnb(scnb) {
        add_term(initial_term);
    }

    void add_term(const ParameterFileDihedralTerm& term) {
        terms.push_back(term);
    }

    std::vector<std::string> types;
    std::vector<ParameterFileDihedralTerm> terms;
    double scee;
    double scnb;
};

struct ParameterFileDihedralPtrLess
        : std::binary_function<const ParameterFileDihedral*,
                               const ParameterFileDihedral*, bool> {
    bool operator()(const ParameterFileDihedral *lhs,
                    const ParameterFileDihedral *rhs) const {
        return lhs->types < rhs->types;
    }
};

struct ParameterFileImproperDihedral {
    ParameterFileImproperDihedral() {}
    ParameterFileImproperDihedral(const std::string& center_type,
                                  const std::vector<std::string>& types,
                                  const ParameterFileDihedralTerm& term)
            : center_type(center_type), types(types), term(term) { init(); }

    void init() {
        types.erase(std::remove(types.begin(), types.end(), "X"),
                    types.end());
        std::sort(types.begin(), types.end());
    }

    bool matches(const std::vector<std::string>& rhs) const;

    std::string center_type;
    std::vector<std::string> types;
    ParameterFileDihedralTerm term;
};

class ImproperDihedralCollection {
  public:
    typedef std::map<std::string,
                     std::vector<ParameterFileImproperDihedral*>
                    > ImproperDihedralMap;
    typedef ImproperDihedralMap::const_iterator const_iterator;

    ImproperDihedralCollection() {}
    ~ImproperDihedralCollection() {
        for (const_iterator it = begin(); it != end(); ++it)
            std::for_each(it->second.begin(), it->second.end(), DeletePtr());
    }

    const_iterator begin() const { return dihedrals_.begin(); }
    const_iterator end() const { return dihedrals_.end(); }

    void insert(ParameterFileImproperDihedral *dihedral) {
        dihedrals_[dihedral->center_type].push_back(dihedral);
    }

    std::pair<ParameterFileDihedralTerm, bool> lookup(
            const std::string& center,
            const std::vector<std::string>& types) const;

    void append(const ImproperDihedralCollection& rhs);

  private:
    ImproperDihedralMap dihedrals_;

    DISALLOW_COPY_AND_ASSIGN(ImproperDihedralCollection);
};

struct HbondParameter {
    std::string type1, type2;
    double A, B;
    std::string description;
};

class ParameterFile {
  public:
    typedef std::map<std::string, ParameterFileAtom*> AtomTypeMap;
    typedef std::set<ParameterFileBond*, ParameterFileBondPtrLess> BondSet;
    typedef std::set<ParameterFileAngle*, ParameterFileAnglePtrLess> AngleSet;
    typedef std::set<ParameterFileDihedral*,
                     ParameterFileDihedralPtrLess> DihedralSet;

    ParameterFile();
    explicit ParameterFile(const std::string& file_name) { read(file_name); }

    ~ParameterFile();

    const AtomTypeMap& atom_types() const { return atom_types_; }
    const BondSet& bonds() const { return bonds_; }
    const AngleSet& angles() const { return angles_; }
    const DihedralSet& dihedrals() const { return dihedrals_; }
    const DihedralSet& generic_dihedrals() const { return generic_dihedrals_; }
    const ImproperDihedralCollection& improper_dihedrals() const {
        return improper_dihedrals_;
    }

  private:
    void read(const std::string& file_name);
    void read(std::istream&);

    void process_atom_type(const std::string& line);
    void process_hydrophilic(const std::string& line);
    void process_bond(const std::string& line);
    void process_angle(const std::string& line);
    void process_dihedral(std::string& line, int& line_index, std::istream&);
    void process_improper_dihedral(const std::string& line);
    void process_hbond_parameters(const std::string& line);
    void process_equivalent_symbol_list(const std::string& line);
    void process_6_12_parameters(const std::string& line);

    double extract_kv_double(const std::string& str, const std::string& key);

    std::string title_;
    AtomTypeMap atom_types_;
    std::vector<std::string> hydrophilic_types_;
    BondSet bonds_;
    AngleSet angles_;
    DihedralSet dihedrals_;
    DihedralSet generic_dihedrals_;
    ImproperDihedralCollection improper_dihedrals_;
    std::vector<HbondParameter> hbond_parameters_;
    std::map<std::string, std::vector<std::string> > equivalent_symbol_lists_;

    DISALLOW_COPY_AND_ASSIGN(ParameterFile);
};

class ParameterFileProcessingException : public std::exception {
  public:
    explicit ParameterFileProcessingException(const std::string& message)
            : line_number_(kNotSet), message_(message) {}
    ParameterFileProcessingException(int line_number,
                                     const std::string& message)
            : line_number_(line_number), message_(message) {}

    virtual const char *what() const throw();

    virtual ~ParameterFileProcessingException() throw() {}

  private:
    int line_number_;
    std::string message_;
};

class ParameterFileSet {
  public:
    ParameterFileSet() {}

    void load(const ParameterFile& parameter_file);
    void load(const std::string& file_name) { load(ParameterFile(file_name)); }

    const ParameterFileAtom *lookup(const std::string& type) const;
    const ParameterFileBond *lookup(const std::string& type1,
                                    const std::string& type2) const;
    const ParameterFileAngle *lookup(const std::string& type1,
                                     const std::string& type2,
                                     const std::string& type3) const;
    const ParameterFileDihedral *lookup(const std::string& type1,
                                        const std::string& type2,
                                        const std::string& type3,
                                        const std::string& type4) const;
    std::pair<ParameterFileDihedralTerm, bool> lookup_improper_dihedral(
            const std::string& center_type,
            const std::string& type1,
            const std::string& type2,
            const std::string& type3) const;

  private:
    void load_atom_types(ParameterFile::AtomTypeMap::const_iterator first,
                         ParameterFile::AtomTypeMap::const_iterator last);
    void load_bonds(ParameterFile::BondSet::const_iterator first,
                    ParameterFile::BondSet::const_iterator last);
    void load_angles(ParameterFile::AngleSet::const_iterator first,
                     ParameterFile::AngleSet::const_iterator last);
    void load_dihedrals(ParameterFile::DihedralSet::const_iterator first,
                        ParameterFile::DihedralSet::const_iterator last);
    void load_improper_dihedrals(const ImproperDihedralCollection& coll) {
        improper_dihedrals_.append(coll);
    }

    ParameterFile::BondSet::const_iterator find(const std::string& type1,
                                                const std::string& type2) const;
    ParameterFile::BondSet::const_iterator find(
            const ParameterFileBond *bond) const {
        return find(bond->types[0], bond->types[1]);
    }
    ParameterFile::AngleSet::const_iterator find(
            const std::string& type1, const std::string& type2,
            const std::string& type3) const;
    ParameterFile::AngleSet::const_iterator find(
            const ParameterFileAngle *angle) const {
        return find(angle->types[0], angle->types[1], angle->types[2]);
    }
    ParameterFile::DihedralSet::const_iterator find(
            const std::string& type1, const std::string& type2,
            const std::string& type3, const std::string& type4) const;
    ParameterFile::DihedralSet::const_iterator find_generic(
            const std::string& type1, const std::string& type2,
            const std::string& type3, const std::string& type4) const;

    void warning(const std::string& message) const;

    ParameterFile::AtomTypeMap atom_types_;
    ParameterFile::BondSet bonds_;
    ParameterFile::AngleSet angles_;
    ParameterFile::DihedralSet dihedrals_;
    ParameterFile::DihedralSet generic_dihedrals_;
    ImproperDihedralCollection improper_dihedrals_;

    DISALLOW_COPY_AND_ASSIGN(ParameterFileSet);
};

inline const ParameterFileAtom *ParameterFileSet::lookup(
        const std::string& type) const {
    ParameterFile::AtomTypeMap::const_iterator it;
    if ((it = atom_types_.find(type)) != atom_types_.end())
        return it->second;
    else
        return NULL;
}

inline const ParameterFileBond *ParameterFileSet::lookup(
        const std::string& type1,
        const std::string& type2) const {
    ParameterFile::BondSet::const_iterator it;
    if ((it = find(type1, type2)) != bonds_.end())
        return *it;
    else
        return NULL;
}

inline const ParameterFileAngle *ParameterFileSet::lookup(
        const std::string& type1,
        const std::string& type2,
        const std::string& type3) const {
    ParameterFile::AngleSet::const_iterator it;
    if ((it = find(type1, type2, type3)) != angles_.end())
        return *it;
    else
        return NULL;
}

inline const ParameterFileDihedral *ParameterFileSet::lookup(
        const std::string& type1,
        const std::string& type2,
        const std::string& type3,
        const std::string& type4) const {
    ParameterFile::DihedralSet::const_iterator it;
    if ((it = find(type1, type2, type3, type4)) == dihedrals_.end()) {
        if ((it = find_generic("X", type2, type3, "X")) ==
                 generic_dihedrals_.end())
             return NULL;
    }
    return *it;
}

inline std::pair<ParameterFileDihedralTerm, bool>
ParameterFileSet::lookup_improper_dihedral(
        const std::string& center_type,
        const std::string& type1,
        const std::string& type2,
        const std::string& type3) const {
    std::vector<std::string> types;
    types.reserve(3);
    types.push_back(type1);
    types.push_back(type2);
    types.push_back(type3);
    return improper_dihedrals_.lookup(center_type, types);
}

}  // namespace gmml

#include "gmml/internal/parameter_file-inl.h"

#endif  // GMML_INTERNAL_PARAMETER_FILE_H_
