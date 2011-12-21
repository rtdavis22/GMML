// Author: Robert Davis

#ifndef GMML_INTERNAL_PDB_FILE_STRUCTURE_H_
#define GMML_INTERNAL_PDB_FILE_STRUCTURE_H_

#include <memory>
#include <string>

#include "gmml/internal/bimap.h"
#include "gmml/internal/structure.h"
#include "gmml/internal/stubs/common.h"

namespace gmml {

class PdbAtomCard;
class PdbConnectCard;
class PdbFile;

class PdbMappingInfo;

class PdbFileStructure : public Structure {
  public:
    virtual ~PdbFileStructure();

    static PdbFileStructure *build(const PdbFile& pdb_file);

    static PdbFileStructure *build(const PdbFile& pdb_file,
                                   const PdbMappingInfo& mapping_info);

    // Returns the index of the atom corresponding to the atom with the given
    // pdb sequence number. -1 is returned if the sequence number is not in
    // the pdb file.
    int map_atom(int pdb_index) const;

    // Returns the index of the residue corresponding to the residue with the
    // given chain identifier, residue sequence number, and insertion code in
    // the pdb file.
    int map_residue(char chain_id, int residue_number,
                    char insertion_code) const;

  private:
    PdbFileStructure();

    struct Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(PdbFileStructure);
};

struct PdbMappingInfo {
    Bimap<std::string, std::string> residue_map;
    Bimap<std::string, std::string> n_terminus_map;
    Bimap<std::string, std::string> c_terminus_map;
    Bimap<std::string, std::string> atom_map;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_PDB_FILE_STRUCTURE_H_
