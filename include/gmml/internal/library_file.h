// Author: Robert Davis

#ifndef LIBRARY_FILE_H
#define LIBRARY_FILE_H

#include "structure.h"

#include <iosfwd>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "boost/shared_ptr.hpp"

namespace gmml {

class LibraryFileStructure;

// This class represents an AMBER OFF library file. The file specification can
// be found here: library file http://ambermd.org/doc/OFF_file_format.txt.
class LibraryFile {
  public:
    typedef boost::shared_ptr<LibraryFileStructure> StructurePtr;
    typedef std::map<std::string, StructurePtr>::iterator iterator;
    typedef std::map<std::string, StructurePtr>::const_iterator
            const_iterator;

    explicit LibraryFile(const std::string& file_name) { read(file_name); }

    const_iterator begin() const { return structures_.begin(); }
    const_iterator end() const { return structures_.end(); }

    // This returns the structure with given name. The name comes from the
    // listing at the top of the file. If the structure isn't present,
    // StructurePtr() is returned.
    const StructurePtr operator[](const std::string& name) const;

  private:
    void read(std::istream&);
    void read(const std::string& file_name);

    // A mapping from the structure names at the top of the file to the
    // corresponding structure in the file.
    std::map<std::string, StructurePtr> structures_;

    // Evil constructors
    LibraryFile(const LibraryFile&);
    void operator=(const LibraryFile&);
};

inline const LibraryFile::StructurePtr LibraryFile::operator[](
        const std::string& name) const {
    const_iterator it = structures_.find(name);
    if (it != structures_.end())
        return it->second;
    else
        return StructurePtr();
}

class LibraryFileSet {
  public:
    typedef LibraryFile::iterator iterator;
    typedef LibraryFile::const_iterator const_iterator;

    LibraryFileSet() {}

    const_iterator begin() const { return structures_.begin(); }
    const_iterator end() const { return structures_.end(); }

    // These add the structures in a library file to the current set of
    // structure. If a structure name already exists in the set, it is
    // overwritten and a warning is displayed.
    void load(const LibraryFile& file);
    void load(const std::string& file_name) { load(LibraryFile(file_name)); }

    // This returns the structure with the given name. If the structure is
    // not present, StructurePtr() is returned.
    const LibraryFile::StructurePtr operator[](const std::string& name) const;

  private:
    // A mapping from the structure names found in the library files their
    // corresonding structures.
    std::map<std::string, LibraryFile::StructurePtr> structures_;

    // Evil constructors
    LibraryFileSet(const LibraryFileSet&);
    void operator=(const LibraryFileSet&);
};

inline const LibraryFile::StructurePtr LibraryFileSet::operator[](
        const std::string& name) const {
    const_iterator it = structures_.find(name);
    if (it != structures_.end())
        return it->second;
    else
        return LibraryFile::StructurePtr();
}

class LibraryFileStructure : public Structure {
  public:
    // This is a representation of the box information found in the file.
    struct Box {
        double angle;
        double length;
        double width;
        double height;
    };

    // This constructor reads a library file structure from a stream. The first
    // line of the stream must be the first atom of the structure.
    // The stream will read until the first atom of the next structure in the
    // stream.
    explicit LibraryFileStructure(std::istream& in) : Structure(),
                                                      box_(NULL) { read(in); }

    virtual ~LibraryFileStructure();

    virtual LibraryFileStructure *clone() const;

    // If the box is not set, NULL is returned.
    const Box *box() const { return box_; }

  private:
    LibraryFileStructure() : Structure(), box_(NULL) {}

    void clone_from(const LibraryFileStructure& structure);

    void read(std::istream& in);
    std::pair<AtomPtr, int> read_atom(std::istream& in) const;
    void read_box(std::istream& in);
    void read_connectivity_info(std::istream& in,
                                const std::map<int, int>& atom_map);
    void read_positions(std::istream& in, const std::map<int, int>& atom_map);
    void read_residue_info(std::istream& in,
                           const std::map<int, int>& residue_map);

    Box *box_;

    // Evil constructors
    LibraryFileStructure(const LibraryFileStructure&);
    void operator=(const LibraryFileStructure&);
};

inline LibraryFileStructure::~LibraryFileStructure() {
    if (box_ != NULL)
        delete box_;
}

inline LibraryFileStructure *LibraryFileStructure::clone() const {
    LibraryFileStructure *structure = new LibraryFileStructure;
    structure->clone_from(*this);
    return structure;
}

inline LibraryFileStructure *build_library_file_structure(
        const LibraryFileStructure& structure) {
    return structure.clone();
}

}  // namespace gmml

#endif  // LIBRARY_FILE_H
