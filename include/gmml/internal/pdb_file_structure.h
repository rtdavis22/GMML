// Author: Robert Davis

#ifndef GMML_INTERNAL_PDB_FILE_STRUCTURE_H_
#define GMML_INTERNAL_PDB_FILE_STRUCTURE_H_

#include <map>
#include <memory>
#include <string>

#include "gmml/internal/bimap.h"
#include "gmml/internal/structure.h"
#include "gmml/internal/stubs/common.h"
#include "gmml/internal/stubs/utils.h"

namespace gmml {

class File;
class PdbAtomCard;
class PdbConnectCard;
class PdbFile;

class PdbMappingInfo;
class PdbStructureBuilder;

// TODO: This should be used throughout in place of Triplet<int>.
struct PdbResidueId {
    PdbResidueId(char chain_id, int res_num, char i_code)
            : chain_id(chain_id), res_num(res_num), i_code(i_code) {}

    char chain_id;
    int res_num;
    char i_code;
};

class PdbFileStructure : public Structure {
  public:
    typedef std::map<Triplet<int>*, int>::const_iterator pdb_iterator;

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
    static PdbFileStructure *build(const PdbStructureBuilder& builder);

    static PdbFileStructure *build(const File& file);

    static PdbFileStructure *build(const PdbFile& pdb_file);

    pdb_iterator pdb_begin() const;
    pdb_iterator pdb_end() const;

    // Returns the index of the atom corresponding to the atom with the given
    // pdb sequence number. -1 is returned if the sequence number is not in
    // the pdb file.
    int map_atom(int pdb_index) const;

    // Returns the index of the residue corresponding to the residue with the
    // given chain identifier, residue sequence number, and insertion code in
    // the pdb file.
    int map_residue(char chain_id, int residue_number,
                    char insertion_code) const;

    int map_residue(char chain_id, int residue_number) {
        return map_residue(chain_id, residue_number, ' ');
    }

    int map_residue(int residue_number) {
        return map_residue(' ', residue_number);
    }

  private:
    PdbFileStructure();

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
    std::string map_pdb_residue(Triplet<int> *pdb_index,
                                const std::string& residue_name,
                                bool is_head, bool is_tail) const;

    PdbFileStructure *build() { return PdbFileStructure::build(*this); }

    const PdbFile& pdb_file() const { return pdb_file_; }

  private:
    // TODO: make this a defensive copy of the PdbFile.
    const PdbFile& pdb_file_;
    PdbMappingInfo mapping_info_;
    std::map<Triplet<int>*, std::string,
             TripletPtrLess<int> > pdb_residue_map_;
};

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

}  // namespace gmml

#endif  // GMML_INTERNAL_PDB_FILE_STRUCTURE_H_
