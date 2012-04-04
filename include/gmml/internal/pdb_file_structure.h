// Author: Robert Davis

#ifndef GMML_INTERNAL_PDB_FILE_STRUCTURE_H_
#define GMML_INTERNAL_PDB_FILE_STRUCTURE_H_

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "gmml/internal/bimap.h"
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

// TODO: This should be used throughout in place of Triplet<int>.
struct PdbResidueId {
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

// Modify to use TripletLess
struct PdbResidueIdLess {
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

struct PdbResidueIdPtrLess {
    bool operator()(const PdbResidueId *lhs, const PdbResidueId *rhs) const {
        return PdbResidueIdLess()(*lhs, *rhs);
    }
};

class PdbFileStructure : public Structure {
  public:
    explicit PdbFileStructure(const PdbStructureBuilder& builder);

    //explicit PdbFileStructure(const File& file);

    //explicit PdbFileStructure(const PdbFile& pdb_file);

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

    int map_residue(const PdbResidueId *pdb_id) const;

    // Convenience functions, maybe shouldn't include them.
    int map_residue(char chain_id, int residue_number,
                    char insertion_code) const {
        PdbResidueId pdb_id(chain_id, residue_number, insertion_code);
        return map_residue(&pdb_id);
    }

    int map_residue(char chain_id, int residue_number) const {
        return map_residue(chain_id, residue_number, ' ');
    }

    int map_residue(int residue_number) const {
        return map_residue(' ', residue_number);
    }

    void append(const IndexedResidue *residue,
                const PdbResidueId *pdb_residue_id);

    const PdbMappingResults *get_mapping_results() const;

    int chain_count() const;

    const PdbChain *chains(int index) const;

  private:
    struct Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(PdbFileStructure);
};

struct PdbMappingInfo {
    Bimap<std::string, std::string> residue_map;
    Bimap<std::string, std::string> head_map;
    Bimap<std::string, std::string> tail_map;
    Bimap<std::string, std::string> atom_map;
};

// This class allows the user to configure how a particular pdb file is built.
// It follows the Builder pattern.
class PdbStructureBuilder {
  public:
    // These initialize mapping_info to the mapping info of the default
    // environment, which are typically specified via the global functions
    // add_mapping, add_head_mapping, and add_tail_mapping. Mappings defined
    // in this class override the global mappings.
    //explicit PdbStructureBuilder(const std::string& pdb_file);
    explicit PdbStructureBuilder(const PdbFile& pdb_file);

    virtual ~PdbStructureBuilder();

    // The following 3 functions add a mapping from a particular residue in the
    // pdb file.
    void add_mapping(char chain_id, int residue_number,
                     char insertion_code, const std::string& name);

    void add_mapping(char chain_id, int residue_number,
                     const std::string& name) {
        add_mapping(chain_id, residue_number, ' ', name);
    }

    void add_mapping(int residue_number, const std::string& name) {
        add_mapping(' ', residue_number, name);
    }

    // Adds a mapping from all residues with a specified name.
    void add_mapping(const std::string& from, const std::string& to) {
        mapping_info_.residue_map.put(from, to);
    }

    // Adds a mapping from all head residues with a specified name. See above
    // for what constitues a head residue.
    void add_head_mapping(const std::string& from, const std::string& to) {
        mapping_info_.head_map.put(from, to);
    }

    // Adds a mapping from all tail residues with a specified name. See above
    // for what constitues a tail residues.
    void add_tail_mapping(const std::string& from, const std::string& to) {
        mapping_info_.tail_map.put(from, to);
    }

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

    PdbFileStructure *build() { return new PdbFileStructure(*this); }

    const PdbFile& pdb_file() const { return pdb_file_; }

    bool use_residue_map() const { return use_residue_map_; }

    bool unknown_hydrogens_removed() const {
        return unknown_hydrogens_removed_;
    }

    void remove_unknown_hydrogens() { unknown_hydrogens_removed_ = true; }

    void keep_unknown_hydrogens() { unknown_hydrogens_removed_ = false; }

    void add_residue_to_remove(const PdbResidueId *residue) {
        residues_to_remove_.insert(new PdbResidueId(*residue));
    }

    bool is_to_be_removed(const PdbResidueId *residue) const {
        return residues_to_remove_.find(residue) != residues_to_remove_.end();
    }

  private:
    const PdbFile& pdb_file_;
    PdbMappingInfo mapping_info_;
    // Change this to PdbResidueId.
    std::map<Triplet<int>*, std::string,
             TripletPtrLess<int> > pdb_residue_map_;
    std::set<const PdbResidueId*, PdbResidueIdPtrLess> residues_to_remove_;
    bool use_residue_map_;
    bool unknown_hydrogens_removed_;

    DISALLOW_COPY_AND_ASSIGN(PdbStructureBuilder);
};

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

    const PdbResidueId *at(int index) const { return residues_.at(index); }

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
