// Author: Robert Davis
//
// This file includes data structures that represent AMBER prep files. It
// includes functionality to parse prep files and build them into Residue
// objects. See http://ambermd.org/doc/prep.html for more information.

#ifndef GMML_INTERNAL_PREP_FILE_H_
#define GMML_INTERNAL_PREP_FILE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "gmml/internal/stubs/common.h"

namespace gmml {

class Residue;

struct PrepFileAtom;
struct PrepFileResidue;

class PrepFile {
  public:
    typedef boost::shared_ptr<PrepFileResidue> ResiduePtr;

    typedef std::map<std::string, ResiduePtr>::iterator iterator;
    typedef std::map<std::string, ResiduePtr>::const_iterator
            const_iterator;

    explicit PrepFile(const std::string& file_name);

    virtual ~PrepFile();

    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;

  private:
    class Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(PrepFile);
};

class PrepFileSet {
  public:
    PrepFileSet();

    virtual ~PrepFileSet();

    PrepFile::ResiduePtr lookup(const std::string& name) const;
    bool exists(const std::string& name) const;

    void load(const PrepFile&);
    void load(const std::string& file) { load(PrepFile(file)); }

    PrepFileResidue& operator[](const std::string& name);
    const PrepFileResidue& operator[](const std::string& name) const;

    PrepFile::iterator begin();
    PrepFile::const_iterator begin() const;
    PrepFile::iterator end();
    PrepFile::const_iterator end() const;

  private:
    class Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(PrepFileSet);
};

struct PrepFileResidue {
    struct ImproperDihedral {
        std::string atom_names[4];
    };

    struct Loop {
        Loop(int from, int to) : from(from), to(to) {}
        int from;
        int to;
    };

    enum CoordinateType { kINT, kXYZ };
    enum GeometryType { kGeometryCorrect, kGeometryChange };
    enum DummyAtomOmission { kOmit, kNomit };
    enum DummyAtomPosition { kPositionAll, kPositionBeg };
    enum OutputFormat { kFormatted, kBinary };

    PrepFileResidue() {}

    ~PrepFileResidue();

    int find(const std::string& name) const;

    std::string header;
    std::string file;
    std::string name;
    CoordinateType coordinate_type;
    OutputFormat output_format;
    GeometryType geometry_type;
    DummyAtomOmission dummy_atom_omission;
    std::string dummy_atom_type;
    DummyAtomPosition dummy_atom_position;
    double cutoff;
    std::vector<PrepFileAtom*> atoms;
    std::vector<ImproperDihedral> improper_dihedrals;
    std::vector<Loop> loops;

  private:
    DISALLOW_COPY_AND_ASSIGN(PrepFileResidue);
};

struct PrepFileAtom {
    enum TopologicalType { kTopTypeM, kTopTypeS, kTopTypeB, kTopTypeE,
                           kTopType3 };

    PrepFileAtom() {}

    int index;
    std::string name;
    std::string type;
    TopologicalType topological_type;
    int bond_index;
    int angle_index;
    int dihedral_index;
    double bond_length;
    double angle;
    double dihedral;
    double charge;

  private:
    DISALLOW_COPY_AND_ASSIGN(PrepFileAtom);
};

class BuildPrepFileResidue {
  public:
    Residue *operator()(const PrepFileResidue& residue) const;
};

inline Residue *build_prep_file(const PrepFileResidue& residue) {
    return BuildPrepFileResidue()(residue);
}

}  // namespace gmml

#endif  // GMML_INTERNAL_PREP_FILE_H_
