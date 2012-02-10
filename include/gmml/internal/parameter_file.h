// Author: Robert Davis
//
// This file includes structures representing data in AMBER force field
// parameter files. See http://ambermd.org/formats.html#parm.dat for the
// file specification.

#ifndef GMML_INTERNAL_PARAMETER_FILE_H_
#define GMML_INTERNAL_PARAMETER_FILE_H_

#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

struct ParameterFileAtom;
struct ParameterFileBond;
struct ParameterFileAngle;
struct ParameterFileDihedral;
struct ParameterFileDihedralTerm;
struct ParameterFileImproperDihedral;
class ImproperDihedralCollection;
class ParameterFileBondPtrLess;
class ParameterFileAnglePtrLess;
class ParameterFileDihedralPtrLess;

class ParameterFile {
  public:
    typedef std::map<std::string, ParameterFileAtom*> AtomTypeMap;
    typedef std::set<ParameterFileBond*, ParameterFileBondPtrLess> BondSet;
    typedef std::set<ParameterFileAngle*, ParameterFileAnglePtrLess> AngleSet;
    typedef std::set<ParameterFileDihedral*,
                     ParameterFileDihedralPtrLess> DihedralSet;

    // Throws ParameterFileProcessingException.
    explicit ParameterFile(const std::string& file_name);

    virtual ~ParameterFile();

    const AtomTypeMap& atom_types() const;
    const BondSet& bonds() const;
    const AngleSet& angles() const;
    const DihedralSet& dihedrals() const;
    const DihedralSet& generic_dihedrals() const;
    const ImproperDihedralCollection& improper_dihedrals() const;

  private:
    class Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(ParameterFile);
};

class ParameterSet {
  public:
    ParameterSet();

    void load(const ParameterFile& parameter_file);

    // Throws ParameterFileProcessingException.
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
    class Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(ParameterSet);
};

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

struct ParameterFileBond {
    ParameterFileBond() : force_constant(kNotSet), length(kNotSet) {}

    ParameterFileBond(const std::string& type1, const std::string& type2);

    ParameterFileBond(const std::vector<std::string>& types,
                      double force_constant, double length)
            : types(types), force_constant(force_constant), length(length) {}

    std::vector<std::string> types;
    double force_constant;
    double length;
};

struct ParameterFileAngle {
    ParameterFileAngle() : force_constant(kNotSet), angle(kNotSet) {}

    ParameterFileAngle(const std::string& type1, const std::string& type2,
                       const std::string& type3);

    ParameterFileAngle(const std::vector<std::string>& types,
                       double force_constant, double angle);

    std::vector<std::string> types;
    double force_constant;
    double angle;
};

// In the file, the periodicity of a term is negative if the next term is to
// be grouped with it. In this case the periodicity in this data structure
// is negative as well.
struct ParameterFileDihedralTerm {
    ParameterFileDihedralTerm() : factor(kNotSet), force_constant(kNotSet),
                                  phase(kNotSet), periodicity(kNotSet) {}

    ParameterFileDihedralTerm(double factor, double force_constant,
                              double phase, double periodicity)
            : factor(factor), force_constant(force_constant), phase(phase),
              periodicity(periodicity) {}

    double factor;
    double force_constant;
    double phase;
    double periodicity;
};

struct ParameterFileDihedral {
    ParameterFileDihedral() : scee(kNotSet), scnb(kNotSet) {}

    ParameterFileDihedral(const std::string& type1, const std::string& type2,
                          const std::string& type3, const std::string& type4);

    ParameterFileDihedral(std::vector<std::string> types,
                          const ParameterFileDihedralTerm& initial_term,
                          double scee, double scnb);

    void add_term(const ParameterFileDihedralTerm& term) {
        terms.push_back(term);
    }

    std::vector<std::string> types;
    std::vector<ParameterFileDihedralTerm> terms;
    double scee;
    double scnb;
};

struct ParameterFileImproperDihedral {
    //ParameterFileImproperDihedral() {}
    ParameterFileImproperDihedral(const std::string& center_type,
                                  const std::vector<std::string>& types,
                                  const ParameterFileDihedralTerm& term)
            : center_type(center_type), types(types), term(term) { init(); }

    bool matches(const std::vector<std::string>& rhs) const;

    std::string center_type;
    std::vector<std::string> types;
    ParameterFileDihedralTerm term;

  private:
    void init();
};

class ImproperDihedralCollection {
  public:
    typedef std::map<std::string,
                     std::vector<ParameterFileImproperDihedral*>
                    > ImproperDihedralMap;

    typedef ImproperDihedralMap::const_iterator const_iterator;

    ImproperDihedralCollection() {}

    ~ImproperDihedralCollection();

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

// This isn't currently part of a public interface because it doesn't seem
// to be relevant (to us, at least). Let us know if you want this information.
struct HbondParameter {
    std::string type1, type2;
    double A, B;
    std::string description;
};

// This exception is throw when a parameter file is formatted incorrectly.
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
    mutable std::string what_;
};

}  // namespace gmml

#include "gmml/internal/parameter_file-inl.h"

#endif  // GMML_INTERNAL_PARAMETER_FILE_H_
