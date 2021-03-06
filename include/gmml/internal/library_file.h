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

#ifndef GMML_INTERNAL_LIBRARY_FILE_H_
#define GMML_INTERNAL_LIBRARY_FILE_H_

#include <iosfwd>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "gmml/internal/boxed_structure.h"
#include "gmml/internal/stubs/common.h"
#include "gmml/internal/stubs/file.h"

namespace gmml {

class LibraryFileStructure;

// This class represents an AMBER OFF library file. The file specification can
// be found here: library file http://ambermd.org/doc/OFF_file_format.txt.
class LibraryFile : public Readable {
  public:
    typedef boost::shared_ptr<LibraryFileStructure> StructurePtr;
    typedef std::map<std::string, StructurePtr>::iterator iterator;
    typedef std::map<std::string, StructurePtr>::const_iterator
            const_iterator;

    explicit LibraryFile(const std::string& file_name) {
        Readable::read(file_name);
    }

    const_iterator begin() const { return structures_.begin(); }
    const_iterator end() const { return structures_.end(); }

    // This returns the structure with given name. The name comes from the
    // listing at the top of the file. If the structure isn't present,
    // StructurePtr() is returned.
    const StructurePtr operator[](const std::string& name) const;

  private:
    void read(std::istream&);

    // A mapping from the structure names at the top of the file to the
    // corresponding structure in the file.
    std::map<std::string, StructurePtr> structures_;

    DISALLOW_COPY_AND_ASSIGN(LibraryFile);
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

    DISALLOW_COPY_AND_ASSIGN(LibraryFileSet);
};

inline const LibraryFile::StructurePtr LibraryFileSet::operator[](
        const std::string& name) const {
    const_iterator it = structures_.find(name);
    if (it != structures_.end())
        return it->second;
    else
        return LibraryFile::StructurePtr();
}

class LibraryFileStructure : public BoxedStructure {
  public:
    // This constructor reads a library file structure from a stream. The first
    // line of the stream must be the first atom of the structure.
    // The stream will read until the first atom of the next structure in the
    // stream.
    explicit LibraryFileStructure(std::istream& in) : BoxedStructure() {
        read(in);
    }

    virtual ~LibraryFileStructure() {}

    virtual LibraryFileStructure *clone() const;

  private:
    LibraryFileStructure() : BoxedStructure() {}

    void clone_from(const LibraryFileStructure& structure);

    void read(std::istream& in);
    std::pair<Atom*, int> read_atom(std::istream& in) const;
    void read_box(std::istream& in);
    void read_connect_atoms(std::istream& in);
    void read_connectivity_info(std::istream& in,
                                const std::map<int, int>& atom_map);
    void read_positions(std::istream& in, const std::map<int, int>& atom_map);
    void read_residue_info(std::istream& in,
                           const std::map<int, int>& residue_map);

    DISALLOW_COPY_AND_ASSIGN(LibraryFileStructure);
};

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

#endif  // GMML_INTERNAL_LIBRARY_FILE_H_
