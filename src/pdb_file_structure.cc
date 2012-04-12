// Author: Robert Davis

#include "gmml/internal/pdb_file_structure.h"

#include <cassert>

#include <deque>
#include <utility>

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
using std::set;
using std::string;
using std::vector;

using boost::bimap;

namespace gmml {
namespace {

class PdbData;

class AddAminoAcidBonds {
  public:
    static const double kDistanceCutoff = 4.0;

    explicit AddAminoAcidBonds(PdbFileStructure *structure)
            : structure_(structure) {}

    void operator()() {
        for (int i = 0; i < structure_->chain_count(); i++) {
            add_chain(*structure_->chains(i));
        }
    }

  private:
    void add_chain(const PdbChain& chain) {
        deque<bool> amino_acid(chain.size(), false);
        for (int i = 0; i < chain.size(); i++) {
            if (is_amino_acid(*chain.at(i))) {
                amino_acid[i] = true;
            }
        }
        for (int i = 1; i < amino_acid.size(); i++) {
            if (amino_acid[i - 1] && amino_acid[i])
                bond_residues(*chain.at(i - 1), *chain.at(i));
        }
    }

    bool is_amino_acid(const PdbResidueId& residue) {
        int index = structure_->map_residue(residue);
        if (index == -1) {
            return false;
        }
        return amino_acids_.lookup(structure_->residues(index)->name());
    }

    void bond_residues(const PdbResidueId& n_side,
                       const PdbResidueId& c_side) {
        int n_side_index = structure_->map_residue(n_side);
        int c_side_index = structure_->map_residue(c_side);
        if (n_side_index != -1 && c_side_index != -1)
            bond_residues(n_side_index, c_side_index);
    }

    void bond_residues(int n_side_index, int c_side_index) {
        int carbon_index = structure_->get_atom_index(n_side_index, "C");
        int nitrogen_index = structure_->get_atom_index(c_side_index, "N");
        if (nitrogen_index != -1 && carbon_index != -1)
            bond_atoms(carbon_index, nitrogen_index);
    }

    void bond_atoms(int atom1_index, int atom2_index) {
        const Atom *atom1 = structure_->atoms(atom1_index);
        const Atom *atom2 = structure_->atoms(atom2_index);
        double distance = measure(atom1->coordinate(), atom2->coordinate());
        if (distance < kDistanceCutoff)
            structure_->add_bond(atom1_index, atom2_index);
    }

    PdbFileStructure *structure_;
    AminoAcidCodeSet amino_acids_;
};

class ApplyResidueMap {
  public:
    ApplyResidueMap(PdbData *data, const PdbStructureBuilder& builder);

    ~ApplyResidueMap();

    PdbMappingResults *operator()();

  private:
    struct Impl;
    std::auto_ptr<Impl> impl_;
};


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
        // The card should contain the pdb id.
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

    void remove_residues(const PdbStructureBuilder& builder) {
        for (int i = 0; i < builder.residues_to_remove_count(); i++) {
            PdbResidueId id(*builder.residues_to_remove(i));
            ResidueMapType::iterator it = pdb_residues_.find(&id);
            if (it != pdb_residues_.end()) {
                delete it->first;
                delete it->second;
                pdb_residues_.erase(it);
            }
        }
    }

    void append_to_structure(PdbFileStructure *structure) {
        map<PdbResidueId*, IndexedResidue*, PdbResidueIdPtrLess>::iterator it;
        for (it = pdb_residues_.begin(); it != pdb_residues_.end(); ++it) {
            structure->append(it->second, it->first);
        }
        add_pdb_bonds(structure);
        AddAminoAcidBonds x(structure);
        x();
        // Why doesn't AddAminoAcidBonds(structure)(); work!!!
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

    friend class ApplyResidueMap;
};

struct ApplyResidueMap::Impl {
    typedef map<PdbResidueId*, IndexedResidue*>::iterator MapIterator;

    Impl(PdbData *data, const PdbStructureBuilder& builder)
            : data_(data), builder_(builder) {}

    PdbMappingResults *operator()() {
        results_ = new PdbMappingResults;

        for (MapIterator it = data_->pdb_residues_.begin();
                it != data_->pdb_residues_.end(); ++it) {
            apply(it);
        }

        return results_;
    }

    void apply(MapIterator it) {
        string mapped_name = get_mapped_name(it);
        Structure *mapped_structure = build(mapped_name);
        if (mapped_structure == NULL) {
            results_->add_unknown_residue(it->first);
            return;
        }
        if (builder_.are_unknown_hydrogens_removed()) {
            remove_unknown_hydrogens(it, *mapped_structure);
        }

        bool found = check_for_unknown_atoms(it->second, *mapped_structure);
        if (found) {
            return;
        }
        Residue *new_residue = CompleteResidue()(it->second,
                                                 mapped_structure);
        delete mapped_structure;

        if (new_residue == NULL) {
            results_->add_unknown_residue(it->first);
            return;
        }

        IndexedResidue *new_indexed_residue = add_indices(new_residue,
                                                          it->second);
        delete new_residue;
        delete it->second;
        it->second = new_indexed_residue;
    }

    void remove_unknown_hydrogens(MapIterator it,
                                  const Structure& mapped_structure) {
        IndexedResidue *residue = it->second;
        for (int i = 0; i < residue->size(); i++) {
            const Atom *atom = residue->atoms(i);
            if (atom->element() != kElementH)
                continue;
            if (!is_atom_in_structure(mapped_structure, atom->name())) {
                int serial = residue->get_atom_index(i);
                results_->add_removed_hydrogen(serial, *it->first,
                                               *residue->atoms(i));
                residue->remove_atom(i);
                i--;
            }
        }
    }

    bool check_for_unknown_atoms(const IndexedResidue *residue,
                                 const Structure& mapped_structure) {
        bool found = false;
        for (int i = 0; i < residue->size(); i++) {
            const Atom *atom = residue->atoms(i);
            if (!is_atom_in_structure(mapped_structure, atom->name())) {
                found = true;
                results_->add_unknown_atom(residue->get_atom_index(i));
            }
        }
        return found;
    }

    bool is_atom_in_structure(const Structure& structure,
                              const string& name) {
        for (int i = 0; i < structure.size(); i++) {
            if (structure.atoms(i)->name() == name)
                return true;
        }
        return false;
    }

    std::string get_mapped_name(MapIterator it) {
        return builder_.map_pdb_residue(it->first, it->second->name(),
                                        data_->is_head(it->first),
                                        data_->is_tail(it->first));
    }

    // Returns a new residue with the structure of the first argument and
    // the atom indices of the second argument.
    // Note: Maybe this should be included in CompleteResidue?
    static IndexedResidue *add_indices(const Residue *residue,
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

    PdbData *data_;
    const PdbStructureBuilder& builder_;
    PdbMappingResults *results_;
};


ApplyResidueMap::ApplyResidueMap(PdbData *data,
                                 const PdbStructureBuilder& builder)
        : impl_(new Impl(data, builder)) {
}

ApplyResidueMap::~ApplyResidueMap() {
}

PdbMappingResults *ApplyResidueMap::operator()() {
    return impl_->operator()();
}

}  // namespace


struct PdbFileStructure::Impl {
    typedef bimap<int, int> AtomMapType;
    typedef bimap<boost::bimaps::set_of<PdbResidueId*, PdbResidueIdPtrLess>,
                  int> ResidueMapType;

    explicit Impl(const PdbStructureBuilder& builder)
            : pdb_data_(builder.pdb_file()),
              mapping_results(NULL) {
        pdb_data_.remove_residues(builder);
        if (builder.is_residue_map_used()) {
            mapping_results = ApplyResidueMap(&pdb_data_, builder)();
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

int PdbFileStructure::map_residue(const PdbResidueId& pdb_id) const {
    Impl::ResidueMapType::left_const_iterator it =
            impl_->residue_map.left.find(const_cast<PdbResidueId*>(&pdb_id));
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
          is_residue_map_used_(true),
          are_unknown_hydrogens_removed_(true) {}

PdbStructureBuilder::~PdbStructureBuilder() {
    map<Triplet<int>*, std::string>::iterator it;
    for (it = pdb_residue_map_.begin(); it != pdb_residue_map_.end(); ++it) {
        delete it->first;
    }
    std::for_each(residues_to_remove_.begin(), residues_to_remove_.end(),
                  DeletePtr());
}

void PdbStructureBuilder::add_mapping(const PdbResidueId& pdb_id,
                                      const string& name) {
    Triplet<int> *pdb_index = new Triplet<int>(pdb_id.chain_id,
                                               pdb_id.res_num,
                                               pdb_id.i_code);
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
