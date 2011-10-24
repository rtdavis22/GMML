#ifndef PREP_FILE_H
#define PREP_FILE_H

#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "residue.h"  // try to get rid of this
#include "utilities.h"  // just for DISALLOW_COPY_AND_ASSIGN, I think

namespace gmml {

class Atom;
class Environment;
class Residue;

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

struct PrepFileResidue {
    enum CoordinateType { kINT, kXYZ };
    enum GeometryType { kGeometryCorrect, kGeometryChange };
    enum DummyAtomOmission { kOmit, kNomit };
    enum DummyAtomPosition { kPositionAll, kPositionBeg };
    enum OutputFormat { kFormatted, kBinary };

    PrepFileResidue() {}

    ~PrepFileResidue() {
        std::for_each(atoms.begin(), atoms.end(), DeletePtr());
    }

    struct ImproperDihedral {
        std::string atom_names[4];
    };
    struct Loop {
        Loop(int from, int to) : from(from), to(to) {}
        int from;
        int to;
    };

    int find(const std::string& name) const {
        for (size_t i = 0; i < atoms.size(); i++)
            if (name == atoms[i]->name)
                return i;
        return -1;
    }

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

class PrepFile {
  public:
    typedef boost::shared_ptr<PrepFileResidue> ResiduePtr;

    typedef std::map<std::string, ResiduePtr>::iterator iterator;
    typedef std::map<std::string, ResiduePtr>::const_iterator
            const_iterator;

    enum OtherSection { kSectionLoop, kSectionImproper, kSectionDone,
                        kSectionOther };

    // Read the file from standard input. This will probably need to be changed
    // if we want to output prep files.
    PrepFile();
    explicit PrepFile(const std::string& file_name) { read(file_name); }

    iterator begin() { return residues_.begin(); }
    const_iterator begin() const { return residues_.begin(); }
    iterator end() { return residues_.end(); }
    const_iterator end() const { return residues_.end(); }

  private:
    void read(const std::string& file_name);
    void read(std::istream&);
    bool process_residue(std::istream&);
    PrepFileResidue::CoordinateType extract_coordinate_type(
            std::istream&) const;
    PrepFileResidue::DummyAtomOmission extract_dummy_omission(
            std::istream&) const;
    PrepFileResidue::DummyAtomPosition extract_dummy_position(
            std::istream&) const;
    PrepFileResidue::OutputFormat extract_output_format(
            std::istream&) const;
    PrepFileResidue::GeometryType extract_geometry_type(
            std::istream&) const;
    PrepFileAtom::TopologicalType extract_topological_type(
            std::istream&) const;
    OtherSection get_other_section(const std::string& line) const;

    std::map<std::string, ResiduePtr> residues_;

    DISALLOW_COPY_AND_ASSIGN(PrepFile);
};

class PrepFileSet {
  public:
    PrepFileSet() {}

    PrepFile::ResiduePtr lookup(const std::string& name) const;
    bool exists(const std::string& name) const;

    void load(const PrepFile&);
    void load(const std::string& file_name) { load(PrepFile(file_name)); }

    PrepFileResidue& operator[](const std::string& name) {
        return *residues_[name];
    }
    const PrepFileResidue& operator[](const std::string& name) const {
        return *(residues_.find(name)->second);
    }

    PrepFile::iterator begin() { return residues_.begin(); }
    PrepFile::const_iterator begin() const { return residues_.begin(); }
    PrepFile::iterator end() { return residues_.end(); }
    PrepFile::const_iterator end() const { return residues_.end(); }

  private:
    std::map<std::string, PrepFile::ResiduePtr> residues_;

    DISALLOW_COPY_AND_ASSIGN(PrepFileSet);
};

inline PrepFile::ResiduePtr PrepFileSet::lookup(const std::string& name) const {
    PrepFile::const_iterator it;
    if ((it = residues_.find(name)) == residues_.end())
        return PrepFile::ResiduePtr();
    else
        return it->second;
}

inline bool PrepFileSet::exists(const std::string& name) const {
    return residues_.find(name) != residues_.end();
}

class BuildPrepFileResidue {
  public:
    Residue *operator()(const PrepFileResidue& residue) const;

  private:
    static void set_parent_list(const std::vector<PrepFileAtom*>& atoms,
                                std::vector<int>& parent_list);
    static void set_dummy_coordinates(const std::vector<PrepFileAtom*>& atoms,
                                      std::vector<Coordinate*>& coordinates);
};

inline Residue *build_prep_file(const PrepFileResidue& residue) {
    BuildPrepFileResidue b;
    return b(residue);
}

}  // namespace gmml

#endif  // PREP_FILE_H
