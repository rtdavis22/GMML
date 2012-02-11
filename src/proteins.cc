// Author: Robert Davis

#include "gmml/internal/proteins.h"

#include "gmml/internal/residue.h"
#include "gmml/internal/structure.h"
#include "gmml/internal/stubs/common.h"
#include "gmml/internal/stubs/logging.h"
#include "utilities.h"

#include <map>

using std::map;
using std::string;
using std::vector;

namespace gmml {
namespace {

const char *kAminoAcidCodes[] = {
  "ALA", "A",
  "ARG", "R",
  "ASN", "N",
  "ASP", "D",
  "CYS", "C",
  "GLN", "Q",
  "GLU", "E",
  "GLY", "G",
  "HIS", "H",
  "ILE", "I",
  "LEU", "L",
  "LYS", "K",
  "MET", "M",
  "PHE", "F",
  "PRO", "P",
  "SER", "S",
  "THR", "T",
  "TRP", "W",
  "TYR", "Y",
  "VAL", "V"
};

}  // namespace gmml

struct AminoAcidCodeSet::Impl {
    Impl();

    // Note: a sorted vector could be used here instead.
    map<string, char> amino_acids;
};

AminoAcidCodeSet::Impl::Impl() {
    for (int i = 0; i < ARRAY_SIZE(kAminoAcidCodes); i += 2) {
        add_or_update_map(amino_acids, kAminoAcidCodes[i],
                          kAminoAcidCodes[i + 1][0]);
    }
}

AminoAcidCodeSet::AminoAcidCodeSet() : impl_(new Impl) {
}

AminoAcidCodeSet::~AminoAcidCodeSet() {
}

bool AminoAcidCodeSet::lookup(const string& code) const {
    return impl_->amino_acids.find(code) != impl_->amino_acids.end();
}

char AminoAcidCodeSet::get_fasta_letter(const std::string& code) const {
    map<string, char>::const_iterator it = impl_->amino_acids.find(code);
    return (it != impl_->amino_acids.end())?it->second:'X';
}

vector<vector<int>*> *find_proteins(const Structure& structure) {
    AminoAcidCodeSet codes;
    vector<vector<int>*> *chains = new vector<vector<int>*>;
    for (int i = 0; i < structure.residue_count(); i++) {
        if (!codes.lookup(structure.residues(i)->name())) {
            continue;
        }

        vector<int> *residues =
                structure.get_adjacent_residues_by_atom(i, "N");
        if (residues->size() > 0) {
            continue;
        }
        
        vector<int> *new_chain = new vector<int>;
        new_chain->push_back(i);
        while (true) {
            // This should never happen.
            if (new_chain->size() > structure.residue_count()) {
                LOG(ERROR) << "Amino acids form a cycle.";
                delete new_chain;
                break;
            }
            vector<int> *next_residue =
                    structure.get_adjacent_residues_by_atom(new_chain->back(),
                                                            "C");
            if (next_residue->size() == 0) {
                break;
            }
            new_chain->push_back((*next_residue)[0]);
            delete next_residue;
        }
        chains->push_back(new_chain);
    }
    return chains;
}

vector<int> *get_asns_with_sequon(const Structure& structure) {
    vector<vector<int>*> *chains = find_proteins(structure);
    vector<int> *asns = new vector<int>;
    for (int i = 0; i < chains->size(); i++) {
        const vector<int>& chain = *(*chains)[i];
        for (int j = 2; j < chain.size(); j++) {
            if (structure.residues(chain[j])->name() != "ASN" ||
                    structure.residues(chain[j - 1])->name() == "PRO") {
                continue;
            }
            string third = structure.residues(chain[j - 2])->name();
            if (third == "SER" || third == "THR") {
                asns->push_back(chain[j]);
            }
        }
    }
    return asns;
}

}  // namespace gmml
