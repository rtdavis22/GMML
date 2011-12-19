#include "gmml/internal/pdb_file_structure.h"

#include <deque>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "gmml/internal/pdb_file.h"
#include "gmml/internal/standard_proteins.h"
#include "gmml/internal/stubs/utils.h"
#include "utilities.h"

using std::cout;
using std::deque;
using std::endl;
using std::map;
using std::pair;
using std::vector;

namespace gmml {
namespace {


}  // namespace

PdbFileStructure *PdbFileStructure::build(const PdbFile& pdb_file,
                                          const PdbMappingInfo& mapping_info) {
    RelevantPdbInfo *pdb_info = get_relevant_pdb_info(pdb_file);

    vector<PdbIndexedResidue*> *residues =
        get_indexed_residues(pdb_info->atom_cards);

    vector<pair<int, int> > bonds_to_add;

    StandardProteins proteins;
    for (vector<PdbIndexedResidue*>::iterator it = residues->begin();
            it != residues->end(); ++it) {
        PdbIndexedResidue *cur_residue = *it;
        PdbIndexedResidue *prev_residue = (*it)->prev_residue;

        const Structure *protein = proteins.get_protein(cur_residue->name);
        if (protein != NULL) {
            vector<int> atom_indices(cur_residue->atoms.size(), -1);
            deque<bool> marked(cur_residue->atoms.size(), false);

            if (protein->size() >= cur_residue->atoms.size()) {

                vector<IndexedAtom*> new_atoms(protein->size(),
                                               static_cast<IndexedAtom*>(NULL));
                bool success = true;
                for (int i = 0; i < cur_residue->atoms.size(); i++) {
                    bool found = false;
                    for (int j = 0; j < protein->size(); j++) {
                        if (protein->atoms(j)->name() ==
                                cur_residue->atoms[i]->atom->name()) {
                            new_atoms[j] = cur_residue->atoms[i];
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        // An atom in the residue did not have a corresponding
                        // atom in the protein, so we quit.
                        success = false;
                        break;
                    }
                    if (!success)
                         break;
                }
                if (success) {
                    // Pick 3 representative points to set the positions of the
                    // hydrogens.
                    int repr_atoms[3];
                    int num_found = 0;
                    for (int i = 0; i < new_atoms.size(); i++) {
                        if (new_atoms[i] != NULL) {
                            repr_atoms[num_found++] = i;
                            if (num_found == 3)
                                break;
                        }
                    }
                    for (int i = 0; i < new_atoms.size(); i++) {
                        if (new_atoms[i] == NULL) {
                            
                            int atom1 = kNotSet;
                            const AdjList& adj_atoms = protein->bonds(i);
                            for (int j = 0; j < adj_atoms.size(); j++) {
                                if (new_atoms[adj_atoms[j]] != NULL) {
                                    atom1 = adj_atoms[j];
                                    break;
                                }
                            }
                            int atom2 = kNotSet;
                            int atom3 = kNotSet;
                            const AdjList& adj_atoms2 = protein->bonds(atom1);
                            for (int j = 0; j < adj_atoms2.size(); j++) {
                                if (new_atoms[adj_atoms2[j]] != NULL) {
                                    bool found = false;
                                    atom2 = adj_atoms2[j];
                                    const AdjList& adj_atoms4 = protein->bonds(
                                        atom2);
                                    for (int k = 0; k < adj_atoms4.size(); k++){
                                        if (adj_atoms4[k] != atom1 &&
                                            new_atoms[adj_atoms4[k]] != NULL &&
                                            adj_atoms4[k] != atom1) {
                                            atom3 = adj_atoms4[k];
                                            found = true;
                                            break;
                                        }
                                    }
                                    if (found)
                                        break;
                                }
                            }
                            if (atom3 == kNotSet) {
                                cout << "AHHH" << endl;
                            }
/*
                            int atom3 = kNotSet;
                            const AdjList& adj_atoms3 = protein->bonds(atom2);
                            for (int j = 0; j < adj_atoms3.size(); j++) {
                                if (adj_atoms3[j] != atom1 &&
                                        new_atoms[adj_atoms3[j]] != NULL &&
                                        new_atoms[adj_atoms3[j]]->atom->element() != kElementH) {
                                    atom3 = adj_atoms3[j];
                                    break;
                                }
                            }
*/

                            AtomPtr new_atom(protein->atoms(i)->clone());




                            double distance = measure(
                                protein->atoms(atom1)->coordinate(),
                                protein->atoms(i)->coordinate());
                            double angle = measure(
                                protein->atoms(atom2)->coordinate(),
                                protein->atoms(atom1)->coordinate(),
                                protein->atoms(i)->coordinate());
                            double dihedral = measure(
                                protein->atoms(atom3)->coordinate(),
                                protein->atoms(atom2)->coordinate(),
                                protein->atoms(atom1)->coordinate(),
                                protein->atoms(i)->coordinate());
                            Coordinate new_crd = calculate_point(
                                new_atoms[atom3]->atom->coordinate(),
                                new_atoms[atom2]->atom->coordinate(),
                                new_atoms[atom1]->atom->coordinate(),
                                angle, dihedral, distance);

                            vector<Coordinate> close_crds;
                            const AdjList& atom1_bonds = protein->bonds(atom1);
                            for (int j = 0; j < atom1_bonds.size(); j++) {
                                if (new_atoms[atom1_bonds[j]] != NULL) {
                                    close_crds.push_back(
                                        new_atoms[atom1_bonds[j]]->atom->coordinate());
                                    
                                }
                            }


                            bool too_close = false;
                            for (int j = 0; j < close_crds.size(); j++) {
                                if (measure(close_crds[j], new_crd) < 1.0) {
                                    too_close = true;
                                    dihedral += 2.0*kPi/3.0;
                                    if (dihedral > kPi)
                                        dihedral -= 2*kPi;
                                    break;
                                }
                            }

                            if (too_close) {
                                new_crd = calculate_point(
                                    new_atoms[atom3]->atom->coordinate(),
                                    new_atoms[atom2]->atom->coordinate(),
                                    new_atoms[atom1]->atom->coordinate(),
                                    angle, dihedral, distance);

                                too_close = false;
                                for (int j = 0; j < close_crds.size(); j++) {
                                    if (measure(close_crds[j], new_crd) < 1.0) {
                                        too_close = true;
                                        dihedral += 2.0*kPi/3.0;
                                        if (dihedral > kPi)
                                            dihedral -= 2*kPi;
                                        break;
                                    }
                                }
                            }

                            if (too_close) {
                                new_crd = calculate_point(
                                    new_atoms[atom3]->atom->coordinate(),
                                    new_atoms[atom2]->atom->coordinate(),
                                    new_atoms[atom1]->atom->coordinate(),
                                    angle, dihedral, distance);

                                too_close = false;
                                for (int j = 0; j < close_crds.size(); j++) {
                                    if (measure(close_crds[j], new_crd) < 1.0) {
                                        too_close = true;
                                        dihedral += 2.0*kPi/3.0;
                                        if (dihedral > kPi)
                                            dihedral -= 2*kPi;
                                        break;
                                    }
                                }
                            }

                            if (too_close)
                                cout << "always too close =(" << endl;
                            new_atom->set_coordinate(new_crd);
                            new_atoms[i] = new IndexedAtom(new_atom, -1);
                        } else {
                            AtomPtr atom = new_atoms[i]->atom;
                            AtomPtr protein_atom = protein->atoms(i);
                            atom->set_type(protein_atom->type());
                            atom->set_charge(protein_atom->charge());
                            atom->set_element(protein_atom->element());
                        }
                        new_atoms[i]->atom->set_type("CG");
                    }
                    cur_residue->atoms = new_atoms;
                    cur_residue->bonds = protein->bonds()->clone();
                }
            }
        }

        // the inter-protein N-C bonding
        if (proteins.is_standard(cur_residue->name) && prev_residue != NULL &&
                proteins.is_standard(prev_residue->name)) {
            int carbon = cur_residue->get_atom_index_by_name("N");
            int nitrogen = prev_residue->get_atom_index_by_name("C");
            if (carbon != -1 && nitrogen != -1)
                bonds_to_add.push_back(std::make_pair(carbon, nitrogen));
        }



    }
 

 
    PdbFileStructure *structure = new PdbFileStructure;
    map<int, int> atom_map;
    map<Triplet<int>*, PdbIndexedResidue*>::const_iterator residue_it;

    // A map from the indexed residues to their indices within the structure.
    map<PdbIndexedResidue*, int> index_map;
    typedef map<PdbIndexedResidue*, int>::iterator residue_iterator;


    for (int i = 0; i < residues->size(); i++) {
        PdbIndexedResidue *cur_residue = (*residues)[i];
        int index = structure->add_indexed_residue(atom_map, *cur_residue);
        index_map[cur_residue] = index;
    }




    structure->add_protein_bonds(index_map);

    for (int i = 0; i < bonds_to_add.size(); i++) {
        int atom1 = atom_map[bonds_to_add[i].first];
        int atom2 = atom_map[bonds_to_add[i].second];
        structure->add_bond(atom1, atom2);
    }

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
