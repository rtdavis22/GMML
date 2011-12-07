#include "gmml/internal/pdb_file_structure.h"

#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "gmml/internal/pdb_file.h"
#include "gmml/internal/stubs/utils.h"
#include "utilities.h"

using std::cout;
using std::endl;
using std::map;
using std::vector;

namespace gmml {
namespace {



}  // namespace

PdbFileStructure *PdbFileStructure::build(const PdbFile& pdb_file) {
    RelevantPdbInfo *pdb_info = get_relevant_pdb_info(pdb_file);

    vector<PdbIndexedResidue*> *residues =
        get_indexed_residues(pdb_info->atom_cards);

 
    PdbFileStructure *structure = new PdbFileStructure;
    map<int, int> atom_map;
    map<Triplet<int>*, PdbIndexedResidue*>::const_iterator residue_it;

    // A map from the indexed residues to their indices within the structure.
    map<PdbIndexedResidue*, int> index_map;

    for (int i = 0; i < residues->size(); i++) {
        PdbIndexedResidue *cur_residue = (*residues)[i];
        int index = structure->add_indexed_residue(atom_map, *cur_residue);
        index_map[cur_residue] = index;
    }

    structure->add_protein_bonds(index_map);

    for (vector<PdbConnectCard*>::const_iterator it =
                pdb_info->connect_cards.begin();
            it != pdb_info->connect_cards.end(); ++it) {
        int source = atom_map[(*it)->connect1];
        if ((*it)->connect2 != kNotSet && atom_map[(*it)->connect2] < source)
             structure->add_bond(source, atom_map[(*it)->connect2]);
        if ((*it)->connect3 != kNotSet && atom_map[(*it)->connect3] < source)
             structure->add_bond(source, atom_map[(*it)->connect3]);
        if ((*it)->connect4 != kNotSet && atom_map[(*it)->connect4] < source)
             structure->add_bond(source, atom_map[(*it)->connect4]);
        if ((*it)->connect5 != kNotSet && atom_map[(*it)->connect5] < source)
             structure->add_bond(source, atom_map[(*it)->connect5]);
    }

    return structure;
}

vector<PdbFileStructure::PdbIndexedResidue*>*
PdbFileStructure::get_indexed_residues(
        const vector<PdbAtomCard*>& atom_cards) {
    // Residues are identified uniquely by their chain id, sequence number, and
    // insertion code.
    map<Triplet<int>*, PdbIndexedResidue*, TripletPtrLess<int> > residue_map;

    typedef map<Triplet<int>*, PdbIndexedResidue*>::iterator iterator;

    // We keep track of the previous residue in the file, so that we can infer
    // bonds between the proteins.
    PdbIndexedResidue *prev_residue = NULL;

    for (vector<PdbAtomCard*>::const_iterator it = atom_cards.begin();
            it != atom_cards.end(); ++it) {
        // A NULL pointer indicates a TER card.
        if (*it == NULL) {
            prev_residue = NULL;
            continue;
        }
        AtomPtr new_atom(new Atom(get_element_by_char((*it)->element[0]),
                                  Coordinate((*it)->x, (*it)->y, (*it)->z),
                                  (*it)->name, "CG", (*it)->charge));

        // Get the triple representing this atom's residue.
        Triplet<int> *triple = new Triplet<int>((*it)->chain_id,
                                                (*it)->res_seq,
                                                (*it)->i_code);

        std::pair<iterator, bool> ret = residue_map.insert(
                std::make_pair(triple, static_cast<PdbIndexedResidue*>(NULL)));
        if (!ret.second)
            delete triple;
        else
            ret.first->second = new PdbIndexedResidue;

        PdbIndexedResidue *cur_residue = ret.first->second;
        cur_residue->atoms.push_back(new IndexedAtom(new_atom,
                                                     (*it)->serial));
        if (cur_residue->name == "")
            cur_residue->name = (*it)->res_name;

        if (cur_residue != prev_residue) {
            cur_residue->prev_residue = prev_residue;
            prev_residue = cur_residue;
        }
    }

    // We don't really care about the triples that identify the residues,
    // so we'll just grab the residues themselves.
    vector<PdbIndexedResidue*> *residues =
            new vector<PdbIndexedResidue*>(residue_map.size());
    int cur_index = 0;
    for (iterator it = residue_map.begin(); it != residue_map.end(); ++it) {
        delete it->first;
        (*residues)[cur_index++] = it->second;
    }

    return residues;
}

PdbFileStructure::RelevantPdbInfo *PdbFileStructure::get_relevant_pdb_info(
        const PdbFile& pdb_file) {
    RelevantPdbInfo *pdb_info = new RelevantPdbInfo;

    PdbFile::const_iterator it = pdb_file.begin();
    while (it != pdb_file.end()) {
        switch (PdbFile::get_card_type(*it)) {
            case PdbFile::ATOM:
                pdb_info->atom_cards.push_back(new PdbAtomCard(*it));
                break;
            case PdbFile::CONECT:
                pdb_info->connect_cards.push_back(new PdbConnectCard(*it));
                break;
            case PdbFile::TER:
                pdb_info->atom_cards.push_back(NULL);
                break;
            default:
                // We don't care about any other card types.
                break;
        }
        ++it;
    }

    return pdb_info;
}

void PdbFileStructure::add_protein_bonds(
        const map<PdbIndexedResidue*, int>& index_map) {
/*
    for (map<PdbIndexedResidue*, int>::const_iterator it = index_map.begin();
            it != index_map.end(); ++it) {
        PdbIndexedResidue *cur_residue = it->first;
        PdbIndexedResidue *prev_residue = cur_residue->prev_residue;
        if (prev_residue != NULL) {
            cout << "add bond between " << it->second << " " <<
                  "(" << cur_residue->name << ")  and " <<
                    index_map[prev_residue] << " (" << prev_residue->name <<
                ")" << endl;
        }
    }
*/
}

}  // namespace gmml
