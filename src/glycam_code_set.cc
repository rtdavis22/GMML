#include "gmml/internal/glycam_code_set.h"

#include <algorithm>
#include <bitset>
#include <map>
#include <stack>
#include <string>
#include <sstream>
#include <vector>

#include "gmml/internal/atom.h"
#include "gmml/internal/glycam_parser.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/residue_classification.h"
#include "gmml/internal/sequence_parser.h"
#include "gmml/internal/structure.h"
#include "gmml/internal/tree_residue.h"
#include "utilities.h"

using std::bitset;
using std::map;
using std::stack;
using std::string;
using std::stringstream;
using std::vector;

namespace gmml {

const char *kNameMap[] = {
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
    vector<string> map(kNameMap, kNameMap + ARRAY_SIZE(kNameMap));
    for (int i = 0; i < map.size(); i += 2) {
        name_to_letter_.insert(std::make_pair(map[i], map[i + 1]));
        letter_to_name_.insert(std::make_pair(map[i + 1], map[i]));
    }
}

string GlycamCodeSet::get_code(const ParsedResidue& parsed_residue,
                               const vector<int>& open_valences) const {
    vector<int> new_valence_list = open_valences;
    map<int, string>::const_iterator it = parsed_residue.derivatives.begin();
    while (it != parsed_residue.derivatives.end()) {
        new_valence_list.push_back(it->first);
        ++it;
    }
    return get_code(parsed_residue.name, parsed_residue.isomer,
                    parsed_residue.configuration, new_valence_list);
}

string GlycamCodeSet::get_code(const string& residue_name,
                               ResidueClassification::Isomer isomer,
                               ResidueClassification::RingType ring_type,
                               ResidueClassification::Configuration config,
                               const vector<int>& open_valences) const {
    bitset<10> bs;
    for (int i = 0; i < open_valences.size(); i++) {
        if (open_valences[i] < 0 || open_valences[i] >= 10)
            throw GlycamCodeException("Invalid open valence");
        bs.set(open_valences[i]);
    }
    string code = get_first_letter(bs) +
                  get_second_letter(residue_name, isomer);
    if (code.size() < 3)
        code += get_third_letter(config, ring_type);
    return code;
}

string GlycamCodeSet::get_code(const string& residue_name,
                               ResidueClassification::Isomer isomer,
                               ResidueClassification::Configuration config,
                               const vector<int>& open_valences) const {
    string new_name = residue_name;
    ResidueClassification::RingType ring_type =
        ResidueClassification::kPyranose;
    // I think this poses a problem if the fourth letter is 'p' or 'f' and is
    // part of the residue name instead of the ring type.
    if (new_name.size() > 3) {
        char ring_letter = new_name[3];
        if (ring_letter == 'p' || ring_letter == 'f') {
            new_name.erase(new_name.begin() + 3);
            if (ring_letter == 'f') {
                ring_type = ResidueClassification::kFuranose;
            }
        }
    }
 
    return get_code(new_name, isomer, ring_type, config, open_valences);
}

string GlycamCodeSet::get_terminal_code(const string& terminal_name) const {
    if (terminal_name == "OH")
        return "ROH";
    else if (terminal_name == "OME")
        return "OME";
    else if (terminal_name == "OtBu")
        return "TBT";

    throw GlycamCodeException("Invalid aglycon " + terminal_name);
}

tree<TreeResidue*> *GlycamCodeSet::build_residue_tree(
        tree<ParsedResidue*> *parsed_tree) const {
    tree<ParsedResidue*>::pre_order_iterator src_it = parsed_tree->begin();
    tree<TreeResidue*> *residue_tree = new tree<TreeResidue*>;
    tree<TreeResidue*>::pre_order_iterator dest_it = residue_tree->begin();

    stack<tree<TreeResidue*>::pre_order_iterator> st;
    dest_it = residue_tree->begin();
    TreeResidue *terminal = new TreeResidue(get_terminal_code((*src_it)->name));
    st.push(residue_tree->insert(dest_it, terminal));
    src_it++;
    while (src_it != parsed_tree->end()) {
        ParsedResidue *parsed_residue = *src_it;
        TreeResidue *tree_residue = build_tree_residue(src_it);
        dest_it = residue_tree->append_child(st.top(), tree_residue);
        map<int, string>::const_iterator map_it;
        for (map_it = parsed_residue->derivatives.begin();
                map_it != parsed_residue->derivatives.end(); ++map_it) {
            residue_tree->append_child(
                    dest_it, get_derivative_tree_residue(map_it->second,
                                                         map_it->first));
        }
        st.push(dest_it);
        int cur_depth = parsed_tree->depth(src_it);
        ++src_it;
        int new_depth = parsed_tree->depth(src_it);
        for (int i = 0; i < cur_depth - new_depth + 1; i++)
            st.pop();
    }

    return residue_tree;
}

ArrayTree<TreeResidue*> *GlycamCodeSet::build_array_tree(
        ArrayTree<ParsedResidue*> *parsed_tree) const {
    ArrayTree<TreeResidue*> *tree = new ArrayTree<TreeResidue*>;

    vector<vector<int> > open_valences(parsed_tree->size());
    for (int i = 0; i < parsed_tree->size(); i++) {
        int parent = (*parsed_tree)[i].second;
        if ((*parsed_tree)[i].second != -1) {
            int oxygen_position = (*parsed_tree)[i].first->oxygen_position;
            open_valences[parent].push_back(oxygen_position);
        }
    }

    string terminal_name = (*parsed_tree)[0].first->name;
    tree->insert(new TreeResidue(get_terminal_code(terminal_name)));

    int cur_derivative_count = 0;
    // We keep a count of how many derivatives have appeared before each
    // residue because the parent indices will need to be modified.
    vector<int> derivatives_before(parsed_tree->size(), 0);
    for (int i = 1; i < parsed_tree->size(); i++) {
        derivatives_before[i] = cur_derivative_count;
        ParsedResidue *parsed_residue = (*parsed_tree)[i].first;
        int parent = (*parsed_tree)[i].second;
        string parent_name = (*parsed_tree)[parent].first->name;
        TreeResidue *tree_residue = build_tree_residue(*parsed_residue,
                                                       open_valences[i],
                                                       parent_name);

        // Insert the residue with possibly a different parent index,
        // accounting for derivatives that may have been added before the
        // parent.
        int residue_index = tree->insert(tree_residue, 
                                         parent + derivatives_before[parent]);

        map<int, string>::const_iterator map_it;
        for (map_it = parsed_residue->derivatives.begin();
                map_it != parsed_residue->derivatives.end(); ++map_it) {
            tree->insert(get_derivative_tree_residue(map_it->second,
                                                     map_it->first),
                         residue_index);
            cur_derivative_count++;
        }
    }

    return tree;
}

std::string GlycamCodeSet::get_name_from_code(const string& code) const {
    string uppercase_code(code);
    std::transform(uppercase_code.begin(), uppercase_code.end(),
                   uppercase_code.begin(), ::toupper);
    if (uppercase_code.substr(1) == "GL")
        return "Neu5Gc";
    // I should probably make a table for these.
    else if (uppercase_code == "SUL")
        return "sulfate";
    else if (uppercase_code == "OME")
        return "OME";
    else if (uppercase_code == "ROH")
        return "OH";
    else if (uppercase_code == "TBT")
        return "OtBu";
    string letter(1, uppercase_code[1]);
    map<string, string>::const_iterator it;
    if ((it = letter_to_name_.find(letter)) != letter_to_name_.end())
        return it->second;
    return "";
}

char GlycamCodeSet::get_first_letter(const vector<int>& positions) const {
    bitset<10> bs;
    for (int i = 0; i < positions.size(); i++) {
        if (positions[i] < 0 || positions[i] > 9) {
            throw std::invalid_argument("Invalid position.");
        }
        bs.set(positions[i]);
    }
    return get_first_letter(bs)[0];
}

string GlycamCodeSet::get_first_letter(const bitset<10>& open_valences) const {
    // This is an alias for convenience.
    const bitset<10>& bs = open_valences;
    if (bs.count() == 4) {
        if (bs[4] && bs[7] && bs[8] && bs[9])
            return "A";
        if (bs[2] && bs[3] && bs[4] && bs[6])
            return "P";
    } else if (bs.count() == 3) {
        if (bs[4] && bs[7] && bs[8])
            return "E";
        if (bs[4] && bs[7] && bs[9])
            return "D";
        if (bs[4] && bs[8] && bs[9])
            return "C";
        if (bs[7] && bs[8] && bs[9])
            return "B";
        if (bs[3] && bs[4] && bs[6])
            return "Q";
        if (bs[2] && bs[4] && bs[6])
            return "R";
        if (bs[2] && bs[3] && bs[6])
            return "S";
        if (bs[2] && bs[3] && bs[4])
            return "T";
    } else if (bs.count() == 2) {
        if (bs[4] && bs[7])
            return "K";
        if (bs[4] && bs[8])
            return "J";
        if (bs[4] && bs[9])
            return "I";
        if (bs[7] && bs[8])
            return "H";
        if (bs[7] && bs[9])
            return "G";
        if (bs[8] && bs[9])
            return "F";
        if (bs[4] && bs[6])
            return "U";
        if (bs[3] && bs[6])
            return "V";
        if (bs[3] && bs[4])
            return "W";
        if (bs[2] && bs[6])
            return "X";
        if (bs[2] && bs[4])
            return "Y";
        if (bs[2] && bs[3])
            return "Z";
    } else if (bs.count() == 1) {
        for (int i = 1; i < bs.size(); i++)
            if (bs[i])
                return to_string(i);
    } else if (bs.none()) {
        return "0";
    }

    stringstream ss;
    ss << "There is no code in the GLYCAM code set for residues with open " <<
          "valences at the given positions.";

    throw GlycamCodeException(ss.str());
}

string GlycamCodeSet::get_second_letter(
        const string& residue_name,
        ResidueClassification::Isomer isomer) const {
    map<string,string>::const_iterator it;
    if ((it = name_to_letter_.find(residue_name)) != name_to_letter_.end()) {
        string letter = it->second;
        if (isomer == ResidueClassification::kIsomerL)
            std::transform(letter.begin(), letter.end(), letter.begin(),
                           ::tolower);
        return letter;
    }
    throw GlycamCodeException(residue_name + " is not a valid residue.");
}

string GlycamCodeSet::get_third_letter(
        ResidueClassification::Configuration configuration,
        ResidueClassification::RingType ring_type) const {
    if (ring_type == ResidueClassification::kPyranose)
        return (configuration == ResidueClassification::kAlpha)?"A":"B";
    else
        return (configuration == ResidueClassification::kAlpha)?"D":"U";
}

TreeResidue *GlycamCodeSet::build_tree_residue(
        tree<ParsedResidue*>::iterator it) const {
    vector<int> open_valences;
    tree<ParsedResidue*>::sibling_iterator child_it = it.begin();
    while (child_it != it.end()) {
        open_valences.push_back((*child_it)->oxygen_position);
        ++child_it;
    }
    tree<ParsedResidue*>::iterator parent_it = tree<ParsedResidue*>::parent(it);
    string parent_name = (*parent_it)->name;
    return build_tree_residue(**it, open_valences, parent_name);
}

TreeResidue *GlycamCodeSet::build_tree_residue(
        const ParsedResidue& parsed_residue,
        const vector<int>& open_valences,
        const string& parent_name) const {
    string anomeric_carbon = "C" + to_string(parsed_residue.anomeric_carbon);
    string oxygen_position = get_oxygen_name(parent_name,
                                             parsed_residue.oxygen_position);
    return new TreeResidue(get_code(parsed_residue, open_valences),
                           anomeric_carbon, oxygen_position);

}

TreeResidue *GlycamCodeSet::get_derivative_tree_residue(
        const string& derivative, int pos) const {
    string oxygen = "O" + to_string(pos);
    if (derivative == "S")
        return new TreeResidue("SUL", "S1", oxygen);
    else if (derivative == "Me")
        return new TreeResidue("MEX", "CH3", oxygen);
    else if (derivative == "A")
        return new TreeResidue("ACE", "C1A", oxygen);

    stringstream ss;
    ss << "There is no derivative in the GLYCAM code set represented by the " <<
          "letter " << derivative;
    throw GlycamCodeException(ss.str());
}

string GlycamCodeSet::get_oxygen_name(const string& residue_code,
                                      int position) const {
    if (residue_code == "OME")
        return "O";
    return "O" + to_string(position);
}

int GlycamAttach::operator()(Structure& structure, Residue *residue,
                             const string& new_atom_name, int residue_index,
                             const string& atom_name) const {
    if (residue->name() == "SUL") {
        int oxygen = structure.get_atom_index(residue_index, atom_name);
        Atom *atom = structure.atoms(oxygen);
        atom->set_charge(atom->charge() + 0.031);
    } else if (residue->name() == "MEX" || residue->name() == "ACE") {
        int oxygen_number = atom_name[1] - '0';
        string carbon_name = "C" + to_string(oxygen_number);
        int carbon = structure.get_atom_index(residue_index, carbon_name);
        Atom *atom = structure.atoms(carbon);
        if (residue->name() == "MEX") {
            atom->set_charge(atom->charge() - 0.039);
        } else if (residue->name() == "ACE") {
            atom->set_charge(atom->charge() + 0.008);
        }
    }

    return structure.attach(residue, new_atom_name, residue_index, atom_name);
}

Structure *glycam_build(ArrayTree<TreeResidue*> *residue_tree) {
    return build_residue_tree(residue_tree, GlycamAttach());
}

Structure *glycam_build(tree<TreeResidue*> *residue_tree) {
    return build_residue_tree(residue_tree, GlycamAttach());
}

Structure *glycam_build(tree<ParsedResidue*> *parsed_tree) {
    GlycamCodeSet code_set;
    tree<TreeResidue*> *residue_tree = code_set.build_residue_tree(parsed_tree);
    Structure *structure = glycam_build(residue_tree);
    std::for_each(residue_tree->begin(), residue_tree->end(), DeletePtr());
    delete residue_tree;
    return structure;
}

Structure *glycam_build(ArrayTree<ParsedResidue*> *parsed_tree) {
    GlycamCodeSet code_set;
    ArrayTree<TreeResidue*> *tree = code_set.build_array_tree(parsed_tree);
    Structure *structure = glycam_build(tree);
    for (ArrayTree<TreeResidue*>::const_iterator it = tree->begin();
            it != tree->end(); ++it)
        delete it->first;
    delete tree;
    return structure;
}

Structure *glycam_build(const string& sequence) {
    SequenceParser *parser = new GlycamParser;
    tree<ParsedResidue*> *parsed_tree = parser->parse(sequence);
    Structure *structure = glycam_build(parsed_tree);
    std::for_each(parsed_tree->begin(), parsed_tree->end(), DeletePtr());
    delete parsed_tree;
    return structure;
}

// TODO: Add some error checking here.
Structure *glycam_build_without_aglycon(const string& sequence) {
    Structure *structure = glycam_build(sequence);
    structure->remove_residue(0);
    structure->set_head(structure->residues(0)->head());
    return structure;
}

Structure *glycam_build_with_array_tree(const string& sequence) {
    GlycamParser *parser = new GlycamParser;
    ArrayTree<ParsedResidue*> *parsed_tree = parser->get_array_tree(sequence);
    Structure *structure = glycam_build(parsed_tree);
    for (ArrayTree<ParsedResidue*>::const_iterator it = parsed_tree->begin();
            it != parsed_tree->end(); ++it)
        delete it->first;
    delete parsed_tree;
    return structure;
}

namespace carbohydrate {

bool set_phi(Structure *structure, int residue_index, double degrees) {
    const vector<Atom*>& atoms = structure->atoms();
    int carbon_index = structure->get_anomeric_index(residue_index);
    if (carbon_index == -1) {
        return false;
    }
    int atom2_index;
    if (atoms[carbon_index]->name() == "C1") {
        atom2_index = structure->get_atom_index(residue_index, "H1");
    }
    else if (atoms[carbon_index]->name()[0] == 'C') {
        int carbon_number = structure->atoms(carbon_index)->name()[1] - '0';
        string atom2_name = "C" + to_string(carbon_number - 1);
        atom2_index = structure->get_atom_index(residue_index, atom2_name);
    } else {
        return false;
    }
    int oxygen_index = kNotSet;
    {
        const Structure::AdjList& adj_list = structure->bonds(carbon_index);
        for (int i = 0; i < adj_list.size(); i++) {
            int atom_index = adj_list[i];
            if (structure->get_residue_index(atom_index) != residue_index &&
                    structure->atoms(atom_index)->name()[0] == 'O') {
                oxygen_index = atom_index;
                break;
            }
        }
        if (oxygen_index == kNotSet ||
                 structure->atoms(oxygen_index)->name().size() < 2 ||
                 !is_number(structure->atoms(oxygen_index)->name()[1])) {
            return false;
        }
    }
    int oxygen_number = char_to_number(atoms[oxygen_index]->name()[1]);
    int attaching_carbon_index = kNotSet;
    {
        const Structure::AdjList& adj_list = structure->bonds(oxygen_index);
        for (int i = 0; i < adj_list.size(); i++) {
            int atom_index = adj_list[i];
            if (structure->get_residue_index(atom_index) ==
                        structure->get_residue_index(oxygen_index) &&
                    atoms[atom_index]->name() ==
                        "C" + to_string(oxygen_number))
            attaching_carbon_index = atom_index;
        }
        if (attaching_carbon_index == kNotSet) {
            return false;
        }
    }
    structure->set_dihedral(attaching_carbon_index, oxygen_index, carbon_index,
                            atom2_index, degrees);
    return true;
}

bool set_psi(Structure *structure, int residue_index, double degrees) {
    const vector<Atom*>& atoms = structure->atoms();
    int anomeric_carbon_index = structure->get_anomeric_index(residue_index);
    if (anomeric_carbon_index == -1) {
        return false;
    }
    int oxygen_index = kNotSet;
    {
        const Structure::AdjList& adj_list =
                structure->bonds(anomeric_carbon_index);
        for (int i = 0; i < adj_list.size(); i++) {
            int atom_index = adj_list[i];
            if (structure->get_residue_index(atom_index) != residue_index &&
                    atoms[atom_index]->name()[0] == 'O') {
                oxygen_index = atom_index;
                break;
            }
        }
    }
    if (oxygen_index == kNotSet || atoms[oxygen_index]->name().size() < 2 ||
            !is_number(atoms[oxygen_index]->name()[1])) {
        return false;
    }
    int oxygen_number = char_to_number(atoms[oxygen_index]->name()[1]);
    int carbon_index = kNotSet;
    {
        const Structure::AdjList& adj_list = structure->bonds(oxygen_index);
        for (int i = 0; i < adj_list.size(); i++) {
            int atom_index = adj_list[i];
            if (structure->get_residue_index(atom_index) ==
                        structure->get_residue_index(oxygen_index) &&
                    atoms[atom_index]->name() ==
                        "C" + to_string(oxygen_number))
                carbon_index = atom_index;
        }
        if (carbon_index == kNotSet) {
            return false;
        }
    }
    string fourth_atom_name;
    if (structure->is_cyclic(carbon_index))
        fourth_atom_name = "H" + to_string(oxygen_number);
    else
        fourth_atom_name = "C" + to_string(oxygen_number - 1);
    int fourth_atom_index = kNotSet;
    {
        const Structure::AdjList& adj_list = structure->bonds(carbon_index);
        for (int i = 0; i < adj_list.size(); i++) {
            int atom_index = adj_list[i];
            if (atoms[atom_index]->name() == fourth_atom_name) {
                fourth_atom_index = atom_index;
                break;
            }
        }
        if (fourth_atom_index == kNotSet) {
            return false;
        }
    }
    structure->set_dihedral(fourth_atom_index, carbon_index, oxygen_index,
                            anomeric_carbon_index, degrees);
    return true;
}

bool set_omega(Structure *structure, int residue_index, double degrees) {
    const vector<Atom*>& atoms = structure->atoms();
    int anomeric_carbon_index = structure->get_anomeric_index(residue_index);
    if (anomeric_carbon_index == -1) {
        return false;
    }
    int oxygen_index = kNotSet;
    int adjacent_residue_index = kNotSet;
    {
        const Structure::AdjList& adj_list =
                structure->bonds(anomeric_carbon_index);
        for (int i = 0; i < adj_list.size(); i++) {
            int atom_index = adj_list[i];
            if (structure->get_residue_index(atom_index) != residue_index &&
                    atoms[atom_index]->name()[0] == 'O') {
                oxygen_index = atom_index;
                adjacent_residue_index =
                        structure->get_residue_index(atom_index);
                break;
            }
        }
        if (oxygen_index == kNotSet ||
                atoms[oxygen_index]->name().size() < 2 ||
                !is_number(atoms[oxygen_index]->name()[1])) {
            return false;
        }
    }
    int oxygen_number = char_to_number(atoms[oxygen_index]->name()[1]);
    int carbon1_index =
            structure->get_atom_index(adjacent_residue_index,
                                      "C" + to_string(oxygen_number));
    if (structure->is_cyclic(carbon1_index))
        return false;
    int carbon2_index =
            structure->get_atom_index(adjacent_residue_index,
                                      "C" + to_string(oxygen_number - 1));
    int other_oxygen_index =
            structure->get_atom_index(adjacent_residue_index,
                                      "O" + to_string(oxygen_number - 1));
    if (carbon1_index == -1 || carbon2_index == -1 || other_oxygen_index == -1)
        return false;
    structure->set_dihedral(other_oxygen_index, carbon2_index, carbon1_index,
                            oxygen_index, degrees);
    return true;
}

void set_default_torsions(Structure *structure, int new_residue_index,
                          int target_residue_index, int carbon_number,
                          int oxygen_number) {
    string oxygen = to_string(oxygen_number);
    string carbon = to_string(carbon_number);

    if (oxygen_number == 5 || oxygen_number == 6) {
        structure->set_dihedral(target_residue_index,
                                "C" + to_string(oxygen_number - 1),
                                target_residue_index, "C" + oxygen,
                                target_residue_index, "O" + oxygen,
                                new_residue_index, "C" + carbon,
                                to_degrees(kPi));
    }
    else {
        structure->set_dihedral(target_residue_index, "H" + oxygen,
                                target_residue_index, "C" + oxygen,
                                target_residue_index, "O" + oxygen,
                                new_residue_index, "C" + carbon,
                                0.0);
    }

    // Regardless of what the actual phi torsion is defined to be, setting
    // C(x+1)-Cx-O'-C' to 180.0 indirectly sets phi to the right thing.
    // So we set this torsion instead of calling set_phi(), which 
    // sets the actual phi torsion.
    structure->set_dihedral(target_residue_index, "C" + oxygen,
                            target_residue_index, "O" + oxygen,
                            new_residue_index, "C" + carbon,
                            new_residue_index,
                            "C" + to_string(carbon_number + 1),
                            180.0);

    // Should probably change this to check if the linkage is exocyclic
    if (oxygen_number == 6) {
        structure->set_dihedral(target_residue_index,
                               "O" + to_string(oxygen_number - 1),
                               target_residue_index,
                               "C" + to_string(oxygen_number - 1),
                               target_residue_index, "C" + oxygen,
                               target_residue_index, "O" + oxygen,
                               60.0);
    }
}

}  // namespace carbohydrate
}  // namespace gmml
