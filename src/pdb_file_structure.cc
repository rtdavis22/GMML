// Author: Robert Davis

#include "gmml/internal/pdb_file_structure.h"

#include <cassert>

#include <deque>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "gmml/internal/atom.h"
#include "gmml/internal/complete_residue.h"
#include "gmml/internal/environment.h"
#include "gmml/internal/pdb_file.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/standard_proteins.h"
#include "gmml/internal/stubs/file.h"
#include "gmml/internal/stubs/utils.h"
#include "utilities.h"

using std::deque;
using std::map;
using std::pair;
using std::string;
using std::vector;

namespace gmml {

// Private implementation
struct PdbFileStructure::Impl {
    // The information in the PDB file that is needed for building the
    // structure. A PdbAtomCard that is NULL indicates a TER card.
    struct RelevantPdbInfo : public PdbCardVisitor {
        virtual void visit(const PdbAtomCard *card) {
            atom_cards.push_back(new PdbAtomCard(*card));
        }

        virtual void visit(const PdbConnectCard *card) {
            connect_cards.push_back(new PdbConnectCard(*card));
        }

        virtual void visit(const PdbTerCard *card) {
            atom_cards.push_back(NULL);
        }

        std::vector<PdbAtomCard*> atom_cards;
        std::vector<PdbConnectCard*> connect_cards;
    };

    // These represent the residues in the PDB file. The index of the atoms is
    // the PDB atom sequence number. These are relevant because they are used
    // in the CONECT cards. This also includes the previous residue and next
    // residue to determine if it's a head residue or a tail residue.
    struct PdbIndexedResidue : public IndexedResidue {
        PdbIndexedResidue(const Residue *residue, PdbIndexedResidue *prev,
                          PdbIndexedResidue *next) : IndexedResidue(residue),
                                                     prev_residue(prev),
                                                     next_residue(next) {}

        PdbIndexedResidue(const std::string& name) : IndexedResidue(name),
                                                     prev_residue(NULL),
                                                     next_residue(NULL) {
            set_bonds(NULL);
        }

        // The residue immediately before this residue in the file. If a TER
        // card precedes this residue, the value is NULL.
        PdbIndexedResidue *prev_residue;

        // The residue immediately after this residue in the file. If a TER
        // card follows this residue, the value is NULL.
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

    // This is map from the PDB atom sequence numbers to the indices of the
    // atoms in the structure.
    map<int, int> atom_map;

    // This is a map from PDB residue identifiers (chain_id, res_num, i_code) to
    // the indices of the residues in the file.
    map<Triplet<int>*, int, TripletPtrLess<int> > residue_map;
};

// This function takes a sequence of PdbAtomCards and returns a map from
// PDB residues identfiers to the residues in the PDB. If the PdbAtomCard is
// NULL, this indicates a TER card.
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
    // head-tail bonding.
    PdbIndexedResidue *prev_residue = NULL;

    for (vector<PdbAtomCard*>::const_iterator it = atom_cards.begin();
            it != atom_cards.end(); ++it) {
        // A NULL pointer indicates a TER card.
        if (*it == NULL) {
            prev_residue = NULL;
            continue;
        }

        Atom *new_atom = new Atom((*it)->element(),
                                  (*it)->coordinate(),
                                  (*it)->name(), "", (*it)->charge());

        // Get the triple representing this atom's residue.
        Triplet<int> *triple = new Triplet<int>((*it)->chain_id(),
                                                (*it)->res_seq(),
                                                (*it)->i_code());

        std::pair<iterator, bool> ret = residue_map->insert(
                std::make_pair(triple, static_cast<PdbIndexedResidue*>(NULL)));
        if (!ret.second) {
            delete triple;
        } else {
            ret.first->second = new PdbIndexedResidue((*it)->res_name());
        }

        PdbIndexedResidue *cur_residue = ret.first->second;
        cur_residue->append(new_atom, (*it)->serial());

        // This may be problematic if the residues in the PDB are mixed up.
        if (cur_residue != prev_residue) {
            cur_residue->prev_residue = prev_residue;
            if (prev_residue != NULL) {
                prev_residue->next_residue = cur_residue;
            }
            prev_residue = cur_residue;
        }
    }

    return residue_map;
}

// Returns the information in the PDB file that we use to build the PDB
// structure.
PdbFileStructure::Impl::RelevantPdbInfo*
PdbFileStructure::Impl::get_relevant_pdb_info(
        const PdbFile& pdb_file) {
    RelevantPdbInfo *info = new RelevantPdbInfo;
    pdb_file.accept(info);
    return info;
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

PdbFileStructure::~PdbFileStructure() {}

PdbFileStructure *PdbFileStructure::build(const PdbFile& pdb_file) {
    return PdbStructureBuilder(pdb_file).build();
}

PdbFileStructure::pdb_iterator PdbFileStructure::pdb_begin() const {
    return impl_->residue_map.begin();
}

PdbFileStructure::pdb_iterator PdbFileStructure::pdb_end() const {
    return impl_->residue_map.end();
}

PdbFileStructure *PdbFileStructure::build(const File& file) {
    return build(PdbFile(file));
}

// This should probably be trimmed down.
PdbFileStructure *PdbFileStructure::build(const PdbStructureBuilder& builder) {
    typedef Impl::PdbIndexedResidue PdbIndexedResidue;
    Impl::RelevantPdbInfo *pdb_info =
            Impl::get_relevant_pdb_info(builder.pdb_file());
    map<Triplet<int>*, PdbIndexedResidue*, TripletPtrLess<int> > *residues =
        Impl::get_indexed_residues(pdb_info->atom_cards);

    vector<pair<int, int> > bonds_to_add;

    StandardProteins proteins;
    for (map<Triplet<int>*, PdbIndexedResidue*>::iterator it =
            residues->begin(); it != residues->end(); ++it) {
        PdbIndexedResidue *cur_residue = it->second;
        PdbIndexedResidue *prev_residue = cur_residue->prev_residue;
        PdbIndexedResidue *next_residue = cur_residue->next_residue;

        bool is_head = (prev_residue == NULL);
        bool is_tail = (next_residue == NULL);

        std::string mapped_name =
                builder.map_pdb_residue(it->first, cur_residue->name(),
                                        is_head, is_tail);
        Residue *result = CompleteResidue()(cur_residue, mapped_name);

        // Now we apply the information in the Residue to the PdbIndexedResidue.
        // We need to set the atom indices of the original IndexedAtoms.
        if (result != NULL) {
            map<string, int> name_map;
            for (int i = 0; i < cur_residue->size(); i++) {
                name_map[cur_residue->atoms(i)->name()] =
                        cur_residue->get_atom_index(i);
            }

            PdbIndexedResidue *new_residue =
                    new (cur_residue) PdbIndexedResidue(result, prev_residue,
                                                        next_residue);
            
            cur_residue = new_residue;

            for (map<string, int>::iterator it = name_map.begin();
                    it != name_map.end(); ++it) {
                new_residue->set_atom_index(it->first, it->second);
            }
        }

        // the inter-protein N-C bonding, change to use head/tail
        if (proteins.is_standard(cur_residue->name()) && prev_residue != NULL &&
                proteins.is_standard(prev_residue->name())) {
            int carbon = cur_residue->get_atom_index("N");
            int nitrogen = prev_residue->get_atom_index("C");
            if (carbon != -1 && nitrogen != -1)
                bonds_to_add.push_back(std::make_pair(carbon, nitrogen));
        }

    }
 
    PdbFileStructure *structure = new PdbFileStructure;
    
    for (map<Triplet<int>*, PdbIndexedResidue*>::iterator it =
            residues->begin(); it != residues->end(); ++it) {
        PdbIndexedResidue *cur_residue = it->second;
        int index = structure->append(structure->impl_->atom_map, cur_residue);
        structure->impl_->residue_map[it->first] = index;
    }

    for (int i = 0; i < bonds_to_add.size(); i++) {
        assert(structure->impl_->atom_map.find(bonds_to_add[i].first) !=
                structure->impl_->atom_map.end());
        assert(structure->impl_->atom_map.find(bonds_to_add[i].second) !=
                structure->impl_->atom_map.end());
        int atom1 = structure->impl_->atom_map[bonds_to_add[i].first];
        int atom2 = structure->impl_->atom_map[bonds_to_add[i].second];
        assert(0 <= atom1 && atom1 < structure->size());
        assert(0 <= atom2 && atom2 < structure->size());
        // Only add head/tail bonds if the atoms are close. This protects
        // against missing TER cards.
        // NOTE: figure out better constant.
        if (measure(structure->atoms(atom1)->coordinate(),
                    structure->atoms(atom2)->coordinate()) < 3.5) {
            structure->add_bond(atom1, atom2);
        }
    }

    for (vector<PdbConnectCard*>::const_iterator it =
                pdb_info->connect_cards.begin();
            it != pdb_info->connect_cards.end(); ++it) {
        int source = structure->impl_->atom_map[(*it)->source()];

        for (int i = 0; i < (*it)->bonded_atom_count(); i++) {
            int bonded_atom = (*it)->get_bonded_atom(i);
            int mapped_atom = structure->impl_->atom_map[bonded_atom];
            if (mapped_atom < source)
                structure->add_bond(source, mapped_atom);
        }
    }

    return structure;
}

/*
PdbStructureBuilder::PdbStructureBuilder(const string& pdb_file)
        : pdb_file_(PdbFile(pdb_file)),
          mapping_info_(*kDefaultEnvironment.pdb_mapping_info()) {
    PdbFile pdb(pdb_file);
    
}
*/

PdbStructureBuilder::PdbStructureBuilder(const PdbFile& pdb_file)
        : pdb_file_(pdb_file),
          mapping_info_(*kDefaultEnvironment.pdb_mapping_info()) {}

PdbStructureBuilder::~PdbStructureBuilder() {
    std::map<Triplet<int>*, std::string>::iterator it;
    for (it = pdb_residue_map_.begin(); it != pdb_residue_map_.end(); ++it) {
        delete it->first;
    }
}

void PdbStructureBuilder::add_mapping(char chain_id, int residue_number,
                                      char insertion_code, const string& name) {
    Triplet<int> *pdb_index = new Triplet<int>(chain_id, residue_number,
                                               insertion_code);
    typedef std::map<Triplet<int>*, string>::iterator iterator;
    pair<iterator, bool> ret =
            pdb_residue_map_.insert(std::make_pair(pdb_index, name));
    // The mapping already exists.
    if (!ret.second) {
        ret.first->second = name;
        delete pdb_index;
    }
}

string PdbStructureBuilder::map_pdb_residue(Triplet<int> *pdb_index,
                                            const string& residue_name,
                                            bool is_head, bool is_tail) const {
    // The pdb_residue_map has the highest priority.
    map<Triplet<int>*, string>::const_iterator it =
            pdb_residue_map_.find(pdb_index);
    if (it != pdb_residue_map_.end())
        return it->second;

    // The head and tail map have the second highest priority.
    if (is_head) {
        pair<string, bool> name = mapping_info_.head_map.get(residue_name);
        if (name.second)
            return name.first;
    }
    if (is_tail) {
        pair<string, bool> name = mapping_info_.tail_map.get(residue_name);
        if (name.second)
            return name.first;
    }

    pair<string, bool> name = mapping_info_.residue_map.get(residue_name);
    if (name.second)
        return name.first;

    return residue_name;
}

}  // namespace gmml
