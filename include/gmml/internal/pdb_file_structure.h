#ifndef GMML_INTERNAL_PDB_FILE_STRUCTURE_H_
#define GMML_INTERNAL_PDB_FILE_STRUCTURE_H_

#include <map>
#include <string>
#include <vector>

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
    PdbFileStructure() {}
    virtual ~PdbFileStructure() {}

    static PdbFileStructure *build(const PdbFile& pdb_file);

    static PdbFileStructure *build(const PdbFile& pdb_file,
                                   const PdbMappingInfo& mapping_info);

  private:
    struct RelevantPdbInfo {
        std::vector<PdbAtomCard*> atom_cards;
        std::vector<PdbConnectCard*> connect_cards;
    };

    struct PdbIndexedResidue : public IndexedResidue {
        PdbIndexedResidue() : prev_residue(NULL) {}

        PdbIndexedResidue *prev_residue;
    };

    static RelevantPdbInfo *get_relevant_pdb_info(const PdbFile& pdb_file);

    static std::vector<PdbIndexedResidue*> *get_indexed_residues(
            const std::vector<PdbAtomCard*>& atom_cards);

    void add_protein_bonds(const std::map<PdbIndexedResidue*, int>& index_map);

    DISALLOW_COPY_AND_ASSIGN(PdbFileStructure);
};

class PdbMappingInfo {
    Bimap<std::string, std::string> residue_map;
    Bimap<std::string, std::string> n_terminus_map;
    Bimap<std::string, std::string> c_terminus_map;
    Bimap<std::string, std::string> atom_map;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_PDB_FILE_STRUCTURE_H_
