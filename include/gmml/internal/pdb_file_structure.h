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

#ifndef GMML_INTERNAL_PDB_FILE_STRUCTURE_H_
#define GMML_INTERNAL_PDB_FILE_STRUCTURE_H_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "gmml/internal/structure.h"
#include "gmml/internal/stubs/common.h"
#include "gmml/internal/stubs/utils.h"

namespace gmml {

class Atom;
class File;
class PdbAtomCard;
class PdbConnectCard;
class PdbFile;

class PdbChain;
class PdbMappingInfo;
class PdbMappingResults;
class PdbStructureBuilder;

struct PdbResidueId {
    struct Less;
    struct PtrLess;

    PdbResidueId() : chain_id(' '), res_num(0), i_code(' ') {}

    PdbResidueId(char chain_id, int res_num, char i_code)
            : chain_id(chain_id), res_num(res_num), i_code(i_code) {}

    PdbResidueId(char chain_id, int res_num)
            : chain_id(chain_id), res_num(res_num), i_code(' ') {}

    explicit PdbResidueId(int res_num)
            : chain_id(' '), res_num(res_num), i_code(' ') {}

    bool equals(const PdbResidueId *rhs) const {
        return chain_id == rhs->chain_id &&
               res_num == rhs->res_num &&
               i_code == rhs->i_code;
    }

    char chain_id;
    int res_num;
    char i_code;
};

struct NamedPdbResidueId {
  public:

    NamedPdbResidueId() {}
    NamedPdbResidueId(std::string res_name, const PdbResidueId& residue) 
            : res_name(res_name), pdb_residue_id(residue) {}

    std::string res_name;
    PdbResidueId pdb_residue_id;
};

struct PdbResidueId::Less {
    bool operator()(const PdbResidueId& lhs, const PdbResidueId& rhs) const {
        if (lhs.chain_id == rhs.chain_id) {
            if (lhs.res_num == rhs.res_num) {
                return lhs.i_code < rhs.i_code;
            } else {
                return lhs.res_num < rhs.res_num;
            }
        } else {
            return lhs.chain_id < rhs.chain_id;
        }
    }
};

struct PdbResidueId::PtrLess {
    bool operator()(const PdbResidueId *lhs, const PdbResidueId *rhs) const {
        return Less()(*lhs, *rhs);
    }
};

class PdbFileStructure : public Structure {
  public:
    explicit PdbFileStructure(const PdbFile& pdb_file,
                              const PdbStructureBuilder& builder);

    virtual ~PdbFileStructure();

    // TODO: document how head/tail residues are inferred and how the mapping
    // is done.
    // build() is the only way to create a PdbFileStructure. The structure is
    // built with the following guarantees:
    // 1. The residues in the structure are ordered by chain_id (chain
    //    identifier), res_num (residue index), and i_code (insertion code) in
    //    that order. chain_id and i_code are interpreted as ASCII chars.
    //    Thus, the residues with chain_id 'Z' will come before the residues
    //    with chain_id 'a', and the residues with no chain_id (' ') will come
    //    before all residues with an alphanumeric chain_id.
    // 2. If a residue is mapped to another residue (via a global mapping or a
    //    mapping in PdbStructureBuilder), the atoms in the residue will be
    //    in the same order as the mapped residue. Otherwise, the atoms are in
    //    the order that they're found in the file.


    // Returns the index of the atom corresponding to the atom with the given
    // pdb sequence number. -1 is returned if the sequence number is not in
    // the pdb file.
    int map_atom(int pdb_index) const;

    const PdbResidueId *map_residue_index(int index) const;

    using Structure::residues;
   
    const Residue *residues(const PdbResidueId& pdb_id) const {
        int index = map_residue(pdb_id);
        return (index == -1)?NULL:residues(index);
    }

    int map_residue(const PdbResidueId& pdb_id) const;

    void append(const IndexedResidue *residue,
                const PdbResidueId *pdb_residue_id);

    const PdbMappingResults *get_mapping_results() const;

    // Returns the number of PDB chains in the file.
    int chain_count() const;

    // Returns the i'th PDB chain in the file. A PDB chain is defined to be a
    // sequence of residues between TER cards.
    const PdbChain *chains(int i) const;

  private:
    struct Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(PdbFileStructure);
};

class PdbMappingInfo {
  public:
    typedef std::map<std::string, std::string> MapType;

    PdbMappingInfo() {}

    std::pair<std::string, bool> lookup_residue_mapping(
            const std::string& residue_name) const;

    std::pair<std::string, bool> lookup_head_mapping(
            const std::string& residue_name) const;

    std::pair<std::string, bool> lookup_tail_mapping(
            const std::string& residue_name) const;

    void add_residue_mapping(const std::string& from, const std::string& to) {
        residue_map_[from] = to;
    }

    void add_head_mapping(const std::string& from, const std::string& to) {
        head_map_[from] = to;
    }

    void add_tail_mapping(const std::string& from, const std::string& to) {
        tail_map_[from] = to;
    }

  private:
    MapType residue_map_;
    MapType head_map_;
    MapType tail_map_;
};

// This class allows the user to configure how a particular PDB file is built.
// It follows the Builder pattern.
class PdbStructureBuilder {
  public:
    // These initialize mapping_info to the mapping info of the default
    // environment, which are typically specified via the global functions
    // add_mapping, add_head_mapping, and add_tail_mapping. Mappings defined
    // in this class override the global mappings.
    PdbStructureBuilder();

    virtual ~PdbStructureBuilder();

    PdbFileStructure *build(const File& file) const;

    PdbFileStructure *build(const PdbFile& pdb_file) const {
        return new PdbFileStructure(pdb_file, *this);
    }

    //
    // Functions to add mappings.
    //
    void add_mapping(const PdbResidueId& pdb_id, const std::string& name);

    // Adds a mapping from all residues with a specified name.
    void add_mapping(const std::string& from, const std::string& to) {
        mapping_info_.add_residue_mapping(from, to);
    }

    // Adds a mapping from all head residues with a specified name. See above
    // for what constitues a head residue.
    void add_head_mapping(const std::string& from, const std::string& to) {
        mapping_info_.add_head_mapping(from, to);
    }

    // Adds a mapping from all tail residues with a specified name. See above
    // for what constitues a tail residues.
    void add_tail_mapping(const std::string& from, const std::string& to) {
        mapping_info_.add_tail_mapping(from, to);
    }


    //
    // Functions for setting other configurable parameters.
    //
    void ignore_residue_map() { is_residue_map_used_ = false; }

    void use_residue_map() { is_residue_map_used_ = true; }

    // If residue mappings are used, this function can be called to remove
    // hydrogens that are present in the PDB file but not the mapped structure.
    void remove_unknown_hydrogens() { are_unknown_hydrogens_removed_ = true; }

    void keep_unknown_hydrogens() { are_unknown_hydrogens_removed_ = false; }

    // The residue with the given identifier will not be included in the
    // structure, along with all associated bonds. It will, however, remain
    // in a PdbChain.
    void add_residue_to_remove(const PdbResidueId *residue) {
        residues_to_remove_.push_back(new PdbResidueId(*residue));
    }


    //
    // Functions for querying the information in the builder.
    //
    // This function determines what a particular residue is mapped to.
    // The function considers mappings in this order:
    // 1. Individual residue mappings.
    // 2. Head mappings (if is_head is true).
    // 3. Tail mappings (if is_tail is true).
    // 4. Generic mappings. 
    // If a mapping is not found, the specified residue name is returned.
    std::string map_pdb_residue(const PdbResidueId *pdb_residue_id,
                                const std::string& residue_name,
                                bool is_head, bool is_tail) const;

    const PdbResidueId *residues_to_remove(int index) const {
        return residues_to_remove_.at(index);
    }

    int residues_to_remove_count() const {
        return residues_to_remove_.size();
    }

    //
    // Accessors
    //
    bool is_residue_map_used() const { return is_residue_map_used_; }

    bool are_unknown_hydrogens_removed() const {
        return are_unknown_hydrogens_removed_;
    }

  private:
    PdbMappingInfo mapping_info_;
    std::map<PdbResidueId*, std::string,
             PdbResidueId::PtrLess> pdb_residue_map_;
    std::vector<PdbResidueId*> residues_to_remove_;
    bool is_residue_map_used_;
    bool are_unknown_hydrogens_removed_;

    DISALLOW_COPY_AND_ASSIGN(PdbStructureBuilder);
};

// This represents an atom that is present in the PDB file but was removed
// when the structure was being built.
class PdbRemovedAtom {
  public:
    PdbRemovedAtom(int serial, const PdbResidueId& residue, const Atom& atom);

    ~PdbRemovedAtom();

    int serial() const { return serial_; }

    const PdbResidueId& residue() const { return residue_; }

    const Atom& atom() const { return *atom_; }

  private:
    const int serial_;
    const PdbResidueId residue_;
    const Atom *atom_;
};

class PdbMappingResults {
  public:
    ~PdbMappingResults() {
        STLDeleteContainerPointers(unknown_residues_.begin(),
                                   unknown_residues_.end());
        STLDeleteContainerPointers(removed_hydrogens_.begin(),
                                   removed_hydrogens_.end());
    }

    void add_unknown_residue(const PdbResidueId *pdb_residue_id) {
        unknown_residues_.push_back(new PdbResidueId(*pdb_residue_id));
    }

    void add_unknown_atom(int serial) {
        unknown_atoms_.push_back(serial);
    }

    void add_removed_hydrogen(int serial, const PdbResidueId& residue,
                              const Atom& atom) {
        removed_hydrogens_.push_back(new PdbRemovedAtom(serial, residue, atom));
    }

    int unknown_residue_count() const { return unknown_residues_.size(); }

    const PdbResidueId *get_unknown_residue(int i) const {
        return unknown_residues_.at(i);
    }

    int unknown_atom_count() const { return unknown_atoms_.size(); }

    int get_unknown_atom(int i) const { return unknown_atoms_.at(i); }

    int removed_hydrogen_count() const { return removed_hydrogens_.size(); }

    const PdbRemovedAtom *get_removed_hydrogen(int i) const {
        return removed_hydrogens_.at(i);
    }

  private:
    std::vector<PdbResidueId*> unknown_residues_;
    std::vector<int> unknown_atoms_;
    std::vector<PdbRemovedAtom*> removed_hydrogens_;
};

// This represents a sequence of residues between TER cards.
class PdbChain {
  public:
    ~PdbChain() {
        STLDeleteContainerPointers(residues_.begin(), residues_.end());
    }

    void append_if_new(const PdbResidueId *pdb_residue_id) {
        if (empty() || !residues_.back()->equals(pdb_residue_id)) {
            append(pdb_residue_id);
        }
    }

    int size() const { return residues_.size(); }

    bool empty() const { return residues_.empty(); }

    const PdbResidueId *at(int i) const { return residues_.at(i); }

    const PdbResidueId *head() const { return residues_.front(); }

    const PdbResidueId *tail() const { return residues_.back(); }

  private:
    void append(const PdbResidueId *pdb_residue_id) {
        residues_.push_back(new PdbResidueId(*pdb_residue_id));
    }

    std::vector<PdbResidueId*> residues_;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_PDB_FILE_STRUCTURE_H_
