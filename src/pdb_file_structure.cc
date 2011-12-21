// Author: Robert Davis

#include "gmml/internal/pdb_file_structure.h"

#include <deque>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "gmml/internal/complete_residue.h"
#include "gmml/internal/pdb_file.h"
#include "gmml/internal/standard_proteins.h"
#include "gmml/internal/stubs/utils.h"
#include "utilities.h"

using boost::shared_ptr;

using std::deque;
using std::map;
using std::pair;
using std::string;
using std::vector;

namespace gmml {

// Private implementation
struct PdbFileStructure::Impl {
    struct RelevantPdbInfo {
        std::vector<PdbAtomCard*> atom_cards;
        std::vector<PdbConnectCard*> connect_cards;
    };

    struct PdbIndexedResidue : public IndexedResidue {
        PdbIndexedResidue() : prev_residue(NULL), next_residue(NULL) {}

        // These help to determine if the residue is a c-terminus or
        // an n-terminus.
        PdbIndexedResidue *prev_residue;
        PdbIndexedResidue *next_residue;
    };

    ~Impl() {
        map<Triplet<int>*, int>::const_iterator it;
        for (it = residue_map.begin(); it != residue_map.end(); ++it)
            delete it->first;
    }

    static RelevantPdbInfo *get_relevant_pdb_info(const PdbFile& pdb_file);

    static map<Triplet<int>*, PdbIndexedResidue*, TripletPtrLess<int> >*
    get_indexed_residues(const vector<PdbAtomCard*>& atom_cards);

    map<int, int> atom_map;
    map<Triplet<int>*, int, TripletPtrLess<int> > residue_map;
};

map<Triplet<int>*, PdbFileStructure::Impl::PdbIndexedResidue*,
    TripletPtrLess<int> >*
PdbFileStructure::Impl::get_indexed_residues(
        const vector<PdbAtomCard*>& atom_cards) {
    // Residues are identified uniquely by their chain id, sequence number, and
    // insertion code.
    map<Triplet<int>*, PdbIndexedResidue*, TripletPtrLess<int> > *residue_map =
            new map<Triplet<int>*, PdbIndexedResidue*, TripletPtrLess<int> >;

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
                                  (*it)->name, "", (*it)->charge));

        // Get the triple representing this atom's residue.
        Triplet<int> *triple = new Triplet<int>((*it)->chain_id,
                                                (*it)->res_seq,
                                                (*it)->i_code);

        std::pair<iterator, bool> ret = residue_map->insert(
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
            if (prev_residue != NULL)
                prev_residue->next_residue = cur_residue;
            prev_residue = cur_residue;
        }
    }

    return residue_map;
}

PdbFileStructure::Impl::RelevantPdbInfo*
PdbFileStructure::Impl::get_relevant_pdb_info(
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

// Public implementation
PdbFileStructure::PdbFileStructure() : impl_(new Impl) {}

int PdbFileStructure::map_atom(int pdb_index) const {
    map<int, int>::iterator it;
    if ((it = impl_->atom_map.find(pdb_index)) != impl_->atom_map.end())
        return it->second;
    return -1;
}

int PdbFileStructure::map_residue(char chain_id, int residue_number,
                                  char insertion_code) const {
    Triplet<int> *triplet = new Triplet<int>(chain_id, residue_number,
                                             insertion_code);
    map<Triplet<int>*, int>::iterator it;
    int ret = -1;
    if ((it = impl_->residue_map.find(triplet)) != impl_->residue_map.end())
        ret = it->second;
    delete triplet;
    return ret;
}

PdbFileStructure::~PdbFileStructure() { }

PdbFileStructure *PdbFileStructure::build(const PdbFile& pdb_file) {
    build(pdb_file, PdbMappingInfo());
}

namespace {

vector<int> *get_connect_atoms(const PdbConnectCard *card) {
    vector<int> *atoms = new vector<int>;
    if (card->connect2 != kNotSet)
        atoms->push_back(card->connect2);
    if (card->connect3 != kNotSet)
        atoms->push_back(card->connect3);
    if (card->connect4 != kNotSet)
        atoms->push_back(card->connect4);
    if (card->connect5 != kNotSet)
        atoms->push_back(card->connect5);
    return atoms;
}

}  // namespace

// This should probably be trimmed down.
PdbFileStructure *PdbFileStructure::build(const PdbFile& pdb_file,
                                          const PdbMappingInfo& mapping_info) {
    typedef Impl::PdbIndexedResidue PdbIndexedResidue;

    Impl::RelevantPdbInfo *pdb_info = Impl::get_relevant_pdb_info(pdb_file);

    map<Triplet<int>*, PdbIndexedResidue*, TripletPtrLess<int> > *residues =
        Impl::get_indexed_residues(pdb_info->atom_cards);

    vector<pair<int, int> > bonds_to_add;

    StandardProteins proteins;
    
    for (map<Triplet<int>*, PdbIndexedResidue*>::iterator it =
            residues->begin(); it != residues->end(); ++it) {
        PdbIndexedResidue *cur_residue = it->second;
        PdbIndexedResidue *prev_residue = cur_residue->prev_residue;
        PdbIndexedResidue *next_residue = cur_residue->next_residue;

        // We first convert the current PdbIndexedResidue to a Residue so
        // we can call CompleteResidue().
        vector<shared_ptr<Atom> > *atoms = new vector<shared_ptr<Atom> >;
        for (int i = 0; i < cur_residue->atoms.size(); i++) {
            atoms->push_back(cur_residue->atoms[i]->atom);
        }

        Residue *residue = new Residue(cur_residue->bonds, cur_residue->name,
                                       atoms->begin(), atoms->end());
        pair<string, bool> mapped_value;
        if (prev_residue == NULL) {
            mapped_value = mapping_info.n_terminus_map.get(cur_residue->name);
        } else if (next_residue == NULL) {
            mapped_value = mapping_info.c_terminus_map.get(cur_residue->name);
        } else {
            mapped_value = mapping_info.residue_map.get(cur_residue->name);
        }
        string mapped_name;
        if (mapped_value.second)
            mapped_name = mapped_value.first;
        else
            mapped_name = cur_residue->name;

        bool success = CompleteResidue()(residue, mapped_name);

        // Now we apply the information in the Residue to the PdbIndexedResidue.
        // We need to set the atom indices of the original IndexedAtoms.
        if (success) {
            cur_residue->bonds = residue->bonds()->clone();
            vector<IndexedAtom*> new_atoms;
            for (int i = 0; i < residue->size(); i++) {
                new_atoms.push_back(new IndexedAtom(residue->atoms(i), -1));
            }

            for (int i = 0; i < cur_residue->atoms.size(); i++) {
                AtomPtr old_atom = cur_residue->atoms[i]->atom;
                for (int j = 0; j < new_atoms.size(); j++) {
                    if (old_atom->name() == new_atoms[j]->atom->name()) {
                        new_atoms[j]->index = cur_residue->atoms[i]->index;
                        break;
                    }
                }
            }
            cur_residue->atoms = new_atoms;
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
    
    for (map<Triplet<int>*, PdbIndexedResidue*>::iterator it =
            residues->begin(); it != residues->end(); ++it) {
        PdbIndexedResidue *cur_residue = it->second;
        int index = structure->add_indexed_residue(structure->impl_->atom_map,
                                                   *cur_residue);
        structure->impl_->residue_map[it->first] = index;
    }

    for (int i = 0; i < bonds_to_add.size(); i++) {
        int atom1 = structure->impl_->atom_map[bonds_to_add[i].first];
        int atom2 = structure->impl_->atom_map[bonds_to_add[i].second];
        structure->add_bond(atom1, atom2);
    }

    for (vector<PdbConnectCard*>::const_iterator it =
                pdb_info->connect_cards.begin();
            it != pdb_info->connect_cards.end(); ++it) {
        int source = structure->impl_->atom_map[(*it)->connect1];

        vector<int> *atoms = get_connect_atoms((*it));
        for (int i = 0; i < atoms->size(); i++) {
            int mapped_atom = structure->impl_->atom_map[(*atoms)[i]];
            if (mapped_atom < source)
                structure->add_bond(source, mapped_atom);
        }
    }

    return structure;
}

}  // namespace gmml
