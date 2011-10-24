#include "gmml/internal/glycam_code_set.h"

#include <bitset>
#include <map>
#include <string>
#include <vector>

#include "gmml/internal/residue_classification.h"
#include "gmml/internal/sequence_parser.h"
#include "gmml/internal/utilities.h"

using std::bitset;
using std::map;
using std::string;
using std::vector;

namespace gmml {

const char *kCodeMap[] = {
    "All", "N",
    "Alt", "E",
    "Ara", "A",
    "Fru", "C",
    "Fuc", "F",
    "Gal", "L",
    "GalA", "O",
    "GalNAc", "V",
    "Glc", "G",
    "GlcA", "Z",
    "GlcNAc", "Y",
    "Gul", "K",
    "Ido", "I",
    "IdoA", "U",
    "Lyx", "D",
    "Man", "M",
    "ManNAc", "W",
    "Neu", "",
    "Neu5Ac", "S",
    "Neu5Gc", "GL",
    "NeuNAc", "S",
    "NeuNGc", "GL",
    "Psi", "P",
    "Qui", "Q",
    "Rha", "H",
    "Rib", "R",
    "Sor", "B",
    "Tag", "J",
    "Tal", "T",
    "Xyl", "X"
};

GlycamCodeSet::GlycamCodeSet() {
    vector<string> map(kCodeMap, kCodeMap + GOOGLE_ARRAYSIZE(kCodeMap));
    for (int i = 0; i < map.size(); i += 2) {
        code_to_letter_.insert(std::make_pair(map[i], map[i + 1]));
        letter_to_code_.insert(std::make_pair(map[i + 1], map[i]));
    }
}

string GlycamCodeSet::get_code(const ParsedResidue& parsed_residue,
                               const vector<int>& open_valences) const {
}

string GlycamCodeSet::get_code(const string& residue_name,
                               ResidueClassification::Isomer isomer,
                               ResidueClassification::RingType ring_type,
                               ResidueClassification::Configuration config,
                               const vector<int>& open_valences) const {


}

string GlycamCodeSet::get_code(const string& residue_name,
                               ResidueClassification::Isomer isomer,
                               ResidueClassification::Configuration config,
                               const vector<int>& open_valences) const {

}

string GlycamCodeSet::get_terminal_code(const string& terminal_name) const {

}

string GlycamCodeSet::get_first_letter(const bitset<10>& open_valences) const {

}

string GlycamCodeSet::get_second_letter(
        const string& residue_name,
        ResidueClassification::Isomer) const {

}

string GlycamCodeSet::get_third_letter(
        ResidueClassification::Configuration configuration,
        ResidueClassification::RingType ring_type) const {

}

}  // namespace gmml
