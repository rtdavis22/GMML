// Author: Robert Davis

#include "gmml/internal/pdb_file_structure.h"

#include <cassert>

#include <deque>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "boost/bimap.hpp"
#include "boost/bimap/set_of.hpp"

#include "gmml/internal/atom.h"
#include "gmml/internal/complete_residue.h"
#include "gmml/internal/environment.h"
#include "gmml/internal/pdb_file.h"
#include "gmml/internal/proteins.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/stubs/file.h"
#include "gmml/internal/stubs/utils.h"
#include "utilities.h"

using std::deque;
using std::map;
using std::pair;
using std::string;
using std::vector;

using boost::bimap;

namespace gmml {
namespace {

class PdbData : public PdbCardVisitor {
  public:
    explicit PdbData(const PdbFile& pdb_file) : cur_chain_(new PdbChain),
                                                ignore_remaining_atoms_(false) {
        chains_.push_back(cur_chain_);
        pdb_file.accept(this);
    }

    ~PdbData() {
        ResidueMapType::const_iterator it = pdb_residues_.begin();
        while (it != pdb_residues_.end()) {
            delete it->first;
            delete it->second;
            ++it;
        }
        std::for_each(chains_.begin(), chains_.end(), DeletePtr());
    }

    virtual void visit(const PdbAtomCard *card) {
        if (ignore_remaining_atoms_) {
            return;
        }
        Atom *atom = new Atom(card->element(), card->coordinate(),
                              card->name(), "", card->charge());
        PdbResidueId *residue_id = new PdbResidueId(card->chain_id(),
                                                    card->res_seq(),
                                                    card->i_code());

        typedef map<PdbResidueId*, IndexedResidue*>::iterator iterator;
        std::pair<iterator, bool> ret = pdb_residues_.insert(
                std::make_pair(residue_id, static_cast<IndexedResidue*>(NULL)));
        if (!ret.second) {
            delete residue_id;
        } else {
            ret.first->second = new IndexedResidue(card->res_name());
            ret.first->second->set_bonds(NULL);
        }

        IndexedResidue *cur_residue = ret.first->second;
        cur_residue->append(atom, card->serial());

        cur_chain_->append_if_new(ret.first->first);
    }

    virtual void visit(const PdbConnectCard *card) {
        int source = card->source();
        for (int i = 0; i < card->bonded_atom_count(); i++) {
            pdb_bonds_.push_back(std::make_pair(source,
                                                card->get_bonded_atom(i)));
        }
    }

    virtual void visit(const PdbTerCard* /* card */) {
        cur_chain_ = new PdbChain;
        chains_.push_back(cur_chain_);
    }

    virtual void visit(const PdbEndMdlCard* /* card */) {
        ignore_remaining_atoms_ = true;
    }

    PdbMappingResults *apply_residue_map(const PdbStructureBuilder& builder) {
        PdbMappingResults *results = new PdbMappingResults;
        map<PdbResidueId*, IndexedResidue*>::iterator it;
        for (it = pdb_residues_.begin(); it != pdb_residues_.end(); ++it) {
            string mapped_name = builder.map_pdb_residue(it->first,
                                                         it->second->name(),
                                                         is_head(it->first),
                                                         is_tail(it->first));
            Structure *s = build(mapped_name);
            if (s == NULL) {
                results->add_unknown_residue(it->first);
                continue;
            }

            bool unknown_atom_found = false;
            for (int i = 0; i < it->second->size(); i++) {
                bool found = false;
                for (int j = 0; j < s->size(); j++) {
                    if (s->atoms(j)->name() == it->second->atoms(i)->name()) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    int serial = it->second->get_atom_index(i);
                    if (it->second->atoms(i)->element() == kElementH) {
                        results->add_removed_hydrogen(serial, *it->first,
                                                      *it->second->atoms(i));
                        it->second->remove_atom(i);
                        i--;
                    } else {
                        unknown_atom_found = true;
                        results->add_unknown_atom(serial);
                    }
                }
            }

            if (unknown_atom_found) continue;

            Residue *r = CompleteResidue()(it->second, s);
            delete s;

            if (r == NULL) {
                results->add_unknown_residue(it->first);
                continue;
            }

            IndexedResidue *new_residue = add_indices(r, it->second);
            delete r;
            delete it->second;
            it->second = new_residue;
        }
        return results;
    }

    // Returns a new residue with the structure of the first argument and
    // the atom indices of the second argument.
    IndexedResidue *add_indices(const Residue *residue,
                                const IndexedResidue *indexed_residue) {
        map<string, int> name_map;
        for (int i = 0; i < indexed_residue->size(); i++) {
            name_map[indexed_residue->atoms(i)->name()] =
                    indexed_residue->get_atom_index(i);
        }

        IndexedResidue *new_residue = new IndexedResidue(residue);
            
        for (map<string, int>::iterator it = name_map.begin();
                it != name_map.end(); ++it) {
            new_residue->set_atom_index(it->first, it->second);
        }

        return new_residue;
    }

    void append_to_structure(PdbFileStructure *structure) {
        map<PdbResidueId*, IndexedResidue*, PdbResidueIdPtrLess>::iterator it;
        for (it = pdb_residues_.begin(); it != pdb_residues_.end(); ++it) {
            structure->append(it->second, it->first);
        }
        add_pdb_bonds(structure);
        add_amino_acid_bonds(structure);
    }

    int chain_count() const { return chains_.size(); }

    const PdbChain *get_chain(int index) const { return chains_.at(index); }

  private:
    typedef map<PdbResidueId*, IndexedResidue*, PdbResidueIdPtrLess>
            ResidueMapType;

    void add_pdb_bonds(PdbFileStructure *structure) {
        for (int i = 0; i < pdb_bonds_.size(); i++) {
            int from = structure->map_atom(pdb_bonds_[i].first);
            int to = structure->map_atom(pdb_bonds_[i].second);
            if (from != -1 && to != -1)
                structure->add_bond(from, to);
        }
    }

    void add_amino_acid_bonds(PdbFileStructure *structure) {
        AminoAcidCodeSet amino_acids;
        for (int i = 0; i < chains_.size(); i++) {
            for (int j = 1; j < chains_[i]->size(); j++) {
                int index1 = map_structure_residue(*structure,
                                                   chains_[i]->at(j - 1));
                int index2 = map_structure_residue(*structure,
                                                   chains_[i]->at(j));
                if (index1 == -1 || index2 == -1) {
                    continue;
                }
                if (amino_acids.lookup(structure->residues(index1)->name()) &&
                        amino_acids.lookup(structure->residues(index2)->name())) {
                    int carbon_id = structure->get_atom_index(index1, "C");
                    int nitrogen_id = structure->get_atom_index(index2, "N");
                    if (carbon_id != -1 && nitrogen_id != -1) {
                        const Atom *carbon = structure->atoms(carbon_id);
                        const Atom *nitrogen = structure->atoms(nitrogen_id);
                        double distance = measure(carbon->coordinate(),
                                                  nitrogen->coordinate());
                        if (distance < 4.0)
                            structure->add_bond(carbon_id, nitrogen_id);
                    }
                }
            }
        }
    }

    static int map_structure_residue(const PdbFileStructure& structure,
                                     const PdbResidueId *id) {
        return structure.map_residue(id->chain_id, id->res_num, id->i_code);
    }

    bool is_head(const PdbResidueId *id) const {
        for (int i = 0; i < chains_.size(); i++) {
            if (chains_[i]->size() > 0 && chains_[i]->head()->equals(id))
                return true;
        }
        return false;
    }

    bool is_tail(const PdbResidueId *id) const {
        for (int i = 0; i < chains_.size(); i++) {
            if (chains_[i]->size() > 0 && chains_[i]->tail()->equals(id))
                return true;
        }
        return false;
    }

    bool ignore_remaining_atoms_;
    ResidueMapType pdb_residues_;
    vector<pair<int, int> > pdb_bonds_;
    vector<PdbChain*> chains_;
    PdbChain *cur_chain_;
};

}  // namespace


struct PdbFileStructure::Impl {
    typedef bimap<int, int> AtomMapType;
    typedef bimap<boost::bimaps::set_of<PdbResidueId*, PdbResidueIdPtrLess>,
                  int> ResidueMapType;

    explicit Impl(const PdbStructureBuilder& builder)
            : pdb_data_(builder.pdb_file()),
              mapping_results(NULL) {
        if (builder.use_residue_map()) {
            mapping_results = pdb_data_.apply_residue_map(builder);
        }
    }

    ~Impl() {
        ResidueMapType::left_const_iterator it = residue_map.left.begin();
        while (it != residue_map.left.end()) {
            delete (it++)->first;
        }
        if (mapping_results != NULL) {
            delete mapping_results;
        }
    }

    // This is map from the PDB atom sequence numbers to the indices of the
    // atoms in the structure.
    AtomMapType atom_map;
    ResidueMapType residue_map;
    PdbData pdb_data_;
    PdbMappingResults *mapping_results;
};


PdbFileStructure::PdbFileStructure(const PdbStructureBuilder& builder)
        : impl_(new Impl(builder)) {
    impl_->pdb_data_.append_to_structure(this);
}

PdbFileStructure::~PdbFileStructure() {
}

int PdbFileStructure::map_atom(int pdb_index) const {
    Impl::AtomMapType::left_map::const_iterator it =
            impl_->atom_map.left.find(pdb_index);
    if (it != impl_->atom_map.left.end()) {
        return it->second;
    }
    return -1;
}

int PdbFileStructure::map_residue(const PdbResidueId *pdb_id) const {
    Impl::ResidueMapType::left_const_iterator it =
            impl_->residue_map.left.find(const_cast<PdbResidueId*>(pdb_id));
    if (it != impl_->residue_map.left.end())
        return it->second;
    return -1;
}

const PdbResidueId *PdbFileStructure::map_residue_index(int index) const {
    Impl::ResidueMapType::right_const_iterator it =
            impl_->residue_map.right.find(index);
    if (it != impl_->residue_map.right.end())
        return it->second;
    return NULL;
}

void PdbFileStructure::append(const IndexedResidue *residue,
                              const PdbResidueId *pdb_residue_id) {
    int residue_index = Structure::append(impl_->atom_map.left, residue);
    PdbResidueId *copy = new PdbResidueId(*pdb_residue_id);
    Impl::ResidueMapType::value_type entry(copy, residue_index);
    impl_->residue_map.insert(entry);
}

const PdbMappingResults *PdbFileStructure::get_mapping_results() const {
    return impl_->mapping_results;
}

int PdbFileStructure::chain_count() const {
    return impl_->pdb_data_.chain_count();
}

const PdbChain *PdbFileStructure::chains(int index) const {
    return impl_->pdb_data_.get_chain(index);
}


PdbStructureBuilder::PdbStructureBuilder(const PdbFile& pdb_file)
        : pdb_file_(pdb_file),
          mapping_info_(*kDefaultEnvironment.pdb_mapping_info()),
          use_residue_map_(true),
          unknown_hydrogens_removed_(true) {}

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

string PdbStructureBuilder::map_pdb_residue(const PdbResidueId *pdb_residue_id,
                                            const string& residue_name,
                                            bool is_head, bool is_tail) const {
    Triplet<int> *pdb_index = new Triplet<int>(pdb_residue_id->chain_id,
                                               pdb_residue_id->res_num,
                                               pdb_residue_id->i_code);
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

PdbRemovedAtom::PdbRemovedAtom(int serial, const PdbResidueId& residue,
                               const Atom& atom)
        : serial_(serial), residue_(residue), atom_(atom.clone()) {
}

PdbRemovedAtom::~PdbRemovedAtom() {
    delete atom_;
}

}  // namespace gmml
