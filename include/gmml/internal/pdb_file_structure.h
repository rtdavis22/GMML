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

class PdbAtomCard;
class PdbConnectCard;
class PdbFile;

class PdbMappingInfo;
class PdbStructureBuilder;

// TODO: This should be used throughout in place of Triplet<int>.
struct PdbIndex {
    char chain_id;
    int res_num;
    char i_code;
};

class PdbFileStructure : public Structure {
  public:
    typedef std::map<Triplet<int>*, int>::const_iterator pdb_iterator;

    virtual ~PdbFileStructure();

    static PdbFileStructure *build(const PdbStructureBuilder& builder);

    static PdbFileStructure *build(const std::string& file);

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
    // environment. I think this is always the desirable behaviour, since
    // mappings can be overridden in this class.
    explicit PdbStructureBuilder(const std::string& pdb_file);
    explicit PdbStructureBuilder(const PdbFile& pdb_file);

    virtual ~PdbStructureBuilder();

    void add_mapping(char chain_id, int residue_number,
                     char insertion_code, const std::string& name);

    void add_mapping(char chain_id, int residue_number,
                     const std::string& name) {
        add_mapping(chain_id, residue_number, ' ', name);
    }

    void add_mapping(int residue_number, const std::string& name) {
        add_mapping(' ', residue_number, name);
    }

    void add_mapping(const std::string& from, const std::string& to) {
        mapping_info_.residue_map.put(from, to);
    }

    void add_head_mapping(const std::string& from, const std::string& to) {
        mapping_info_.head_map.put(from, to);
    }

    void add_tail_mapping(const std::string& from, const std::string& to) {
        mapping_info_.tail_map.put(from, to);
    }

    std::string map_pdb_residue(Triplet<int> *pdb_index,
                                const std::string& residue_name,
                                bool is_head, bool is_tail) const;

    PdbFileStructure *build() { return PdbFileStructure::build(*this); }

    const PdbFile& pdb_file() const { return pdb_file_; }

  private:
    const PdbFile& pdb_file_;
    PdbMappingInfo mapping_info_;
    std::map<Triplet<int>*, std::string,
             TripletPtrLess<int> > pdb_residue_map_;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_PDB_FILE_STRUCTURE_H_
