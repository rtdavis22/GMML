#ifndef GMML_INTERNAL_PDB_FILE_STRUCTURE_H_
#define GMML_INTERNAL_PDB_FILE_STRUCTURE_H_

#include <map>
#include <vector>

#include "structure.h"

namespace gmml {

class PdbAtomCard;
class PdbConnectCard;
class PdbFile;

class PdbFileStructure : public Structure {
  public:
    PdbFileStructure() {}
    virtual ~PdbFileStructure() {}

    static PdbFileStructure *build(const PdbFile& pdb_file);

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
};

}  // namespace gmml

#endif  // GMML_INTERNAL_PDB_FILE_STRUCTURE_H_
