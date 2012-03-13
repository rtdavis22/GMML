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
#include "gmml/internal/stubs/file.h"

namespace gmml {

class Residue;

struct PrepFileAtom;
struct PrepFileResidue;

class PrepFile : public Readable, public Writeable {
  public:
    typedef boost::shared_ptr<PrepFileResidue> ResiduePtr;

    typedef std::map<std::string, ResiduePtr>::iterator iterator;
    typedef std::map<std::string, ResiduePtr>::const_iterator
            const_iterator;

    PrepFile() {}

    explicit PrepFile(const std::string& file_name) {
        Readable::read(file_name);
    }

    virtual ~PrepFile() {}

    virtual void write(std::ostream&) const;

    void add_residue(ResiduePtr residue);

    void remove_residue(const std::string& code) { residues_.erase(code); }

    iterator begin() { return residues_.begin(); }
    const_iterator begin() const { return residues_.begin(); }

    iterator end() { return residues_.end(); }
    const_iterator end() const { return residues_.end(); }

    void set_header1(const std::string& header) { header1_ = header; }
    void set_header2(const std::string& header) { header2_ = header; }

    std::string header1() { return header1_; }
    std::string header2() { return header2_; }

  private:
    virtual void read(std::istream&);

    std::map<std::string, ResiduePtr> residues_;
    std::string header1_;
    std::string header2_;

    DISALLOW_COPY_AND_ASSIGN(PrepFile);
};

class PrepFileSet {
  public:
    PrepFileSet() {}

    virtual ~PrepFileSet() {}

    PrepFile::ResiduePtr lookup(const std::string& name) const;

    bool exists(const std::string& name) const {
        return lookup(name) != PrepFile::ResiduePtr();
    }

    void load(const PrepFile&);
    void load(const std::string& file) { load(PrepFile(file)); }

    // Why are these here? There lookup().
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

class PrepFileResidue {
  public:
    enum CoordinateType { kINT, kXYZ };
    enum GeometryType { kGeometryCorrect, kGeometryChange };
    enum DummyAtomOmission { kOmit, kNomit };
    enum DummyAtomPosition { kPositionAll, kPositionBeg };
    enum OutputFormat { kFormatted, kBinary };
    
    struct ImproperDihedral {
        std::string atom_names[4];
    };

    class Loop;

    explicit PrepFileResidue(const Residue& residue);

    ~PrepFileResidue();

    static PrepFileResidue *read_from_stream(std::istream& in);

    void write(std::ostream&) const;

    // change name, what uses this?, maybe should be moved/removed
    int find(const std::string& name) const;

    // Accessors
    int atom_count() const;
    const PrepFileAtom *atoms(int index) const;

    int loop_count() const;
    const Loop *loops(int index) const;

    void set_header(const std::string& header);

    std::string header() const;
    std::string file() const;
    std::string name() const;
    CoordinateType coordinate_type() const;
    OutputFormat output_format() const;
    GeometryType geometry_type() const;
    DummyAtomOmission dummy_atom_omission() const;
    std::string dummy_atom_type() const;
    DummyAtomPosition dummy_atom_position() const;
    double cutoff() const;

  private:
    struct Impl;
    std::auto_ptr<Impl> impl_;

    class CreatePrepFile;
    // This next line is not necessary in C++0x.
    friend class CreatePrepFile;

    PrepFileResidue();

    DISALLOW_COPY_AND_ASSIGN(PrepFileResidue);
};

class PrepFileResidue::Loop {
  public:
    int from() const { return from_; }
    int to() const { return to_; }

  public:
    // Make PrepFileResidue a friend?
    Loop(int from, int to) : from_(from), to_(to) {}

    int from_;
    int to_;

    //just make constructor a friend
    friend class PrepFileResidue;
};

class PrepFileAtom {
  public:
    enum TopologicalType { kTopTypeE, kTopTypeS, kTopTypeB, kTopType3,
                           kTopType4, kTopTypeM };

    explicit PrepFileAtom(const std::string& line);

    ~PrepFileAtom();

    static TopologicalType get_topological_type(int count);
    static std::string get_top_type_string(TopologicalType type);

    void write(std::ostream&) const;

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

  public:
    struct Impl;
    std::auto_ptr<Impl> impl_;

    PrepFileAtom() {}

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
