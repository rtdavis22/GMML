// Author: Robert Davis

#include "gmml/internal/pdb_file.h"

#include <cassert>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "gmml/internal/environment.h"
#include "gmml/internal/geometry.h"
#include "utilities.h"

namespace gmml {

using std::deque;
using std::istringstream;
using std::string;
using std::vector;

void PdbFile::read(std::istream& in) {
    string line;
    while (getline(in, line)) {
        PdbLine pdb_line(line);
        cards_.push_back(pdb_line.get_card());
    }
}

void PdbFile::write(std::ostream& out) const {
    for (int i = 0; i < cards_.size(); i++) {
        cards_[i]->write(out);
        out << std::endl;
    }
}

void PdbFile::accept(PdbCardVisitor *visitor) const {
    for (int i = 0; i < cards_.size(); i++) {
        cards_[i]->accept(visitor);
    }
}

PdbCard *PdbLine::get_card() const {
    switch (get_card_type()) {
        case ATOM:
            return new PdbAtomCard(*this);
        case TER:
            return new PdbTerCard(*this);
        case CONECT:
            return new PdbConnectCard(*this);
        case END:
            return new PdbEndCard(*this);
        case LINK:
            return new PdbLinkCard(*this);
        case UNKNOWN:
            return new PdbUnknownCard(*this);
        default:
            assert(false);
    }
}

PdbLine::CardType PdbLine::get_card_type() const {
    string first_six = data_.substr(0, 6);
    first_six.resize(6, ' ');
    if (first_six == "ATOM  " || first_six == "HETATM")
        return ATOM;
    else if (first_six == "TER   ")
        return TER;
    else if (first_six == "CONECT")
        return CONECT;
    else if (first_six == "END   ")
        return END;
    else if (first_six == "LINK  ")
        return LINK;
    else if (first_six == "ENDMDL")
        return ENDMDL;
    else
        return UNKNOWN;
}

PdbAtomCardBuilder::PdbAtomCardBuilder()
        : serial_(-1), name_(""), alt_loc_(' '), res_name_(""), chain_id_(' '),
          res_seq_(-1), i_code_(' '), coordinate_(0.0, 0.0, 0.0),
          occupancy_(1.0), temp_factor_(0.0), element_(kElementUnknown),
          charge_(kNotSet), is_hetatm_(false) {
}

void PdbAtomCardBuilder::initialize_from_atom(const Atom& atom) {
    name_ = atom.name();
    coordinate_ = atom.coordinate();
    element_ = atom.element();
}

PdbAtomCard *PdbAtomCardBuilder::build() const {
    validate();
    return new PdbAtomCard(*this);
}

void PdbAtomCardBuilder::validate() const {
    // Check here to make sure all required info is set.
}

PdbAtomCard::PdbAtomCard(const PdbAtomCardBuilder& builder)
        : serial_(builder.serial()),
          name_(builder.name()),
          alt_loc_(builder.alt_loc()),
          res_name_(builder.res_name()),
          chain_id_(builder.chain_id()),
          res_seq_(builder.res_seq()),
          i_code_(builder.i_code()),
          coordinate_(builder.coordinate()),
          occupancy_(builder.occupancy()),
          temp_factor_(builder.temp_factor()),
          element_(builder.element()),
          charge_(builder.charge()),
          is_hetatm_(builder.is_hetatm()) {
}

void PdbAtomCard::write(std::ostream& out) const {
    string str(80, ' ');
    if (!is_hetatm_) {
        set_in_string(str, "ATOM", 0, 6, 'L');
    } else {
        set_in_string(str, "HETATM", 0, 6, 'L');
    }
    set_in_string(str, serial_, 6, 5, 'R');
    if (name_.size() < 4)
        set_in_string(str, name_, 13, 3, 'L');
    else
        set_in_string(str, name_, 12, 4, 'L');
    str[16] = alt_loc_;
    set_in_string(str, res_name_, 17, 3, 'L');
    str[21] = chain_id_;
    set_in_string(str, res_seq_, 22, 4, 'R');
    str[26] = i_code_;
    set_in_string(str, coordinate_.x, 30, 8, 'R', 3);
    set_in_string(str, coordinate_.y, 38, 8, 'R', 3);
    set_in_string(str, coordinate_.z, 46, 8, 'R', 3);
    if (is_set(occupancy_))
        set_in_string(str, occupancy_, 54, 6, 'R', 2);
    if (is_set(temp_factor_))
        set_in_string(str, temp_factor_, 60, 6, 'R', 2);
    set_in_string(str, get_element_symbol(element_), 76, 2, 'R');
    if (is_set(charge_))
        set_in_string(str, charge_, 78, 2, 'R');
    trim(str);
    out << str;
}

void PdbAtomCard::read(const string& line) {
    string input(line);
    if (input.size() < 80)
	input.resize(80, ' ');
    serial_ = convert_string<int>(input.substr(6, 5));
    name_ = input.substr(12, 4);
    trim(name_);
    alt_loc_ = input[16];
    res_name_ = input.substr(17, 3);
    trim(res_name_);
    chain_id_ = input[21];
    res_seq_ = convert_string<int>(input.substr(22, 4));
    i_code_ = input[26];
    coordinate_.x = convert_string<double>(input.substr(30, 8));
    coordinate_.y = convert_string<double>(input.substr(38, 8));
    coordinate_.z = convert_string<double>(input.substr(46, 8));
    try {
        occupancy_ = convert_string<double>(input.substr(54, 6));
    } catch (const ConversionException& e) {
        occupancy_ = kNotSet;
    }
    try {
        temp_factor_ = convert_string<double>(input.substr(60, 6));
    } catch (const ConversionException e) {
        temp_factor_ = kNotSet;
    }
    string element_str = input.substr(76, 2);
    trim(element_str);
    element_ = get_element_by_symbol(element_str.c_str());
    if (element_ == kElementUnknown) {
        element_ = get_element_by_name(name_);
    }
    try {
        charge_ = convert_string<double>(input.substr(78, 2));
    } catch (const ConversionException& e) {
        charge_ = kNotSet;
    }
}

void PdbTerCard::write(std::ostream& out) const {
    out << "TER";
}

vector<PdbConnectCard*> *PdbConnectCard::create_cards(
        int source, const vector<int>& bonded_atoms) {
    vector<PdbConnectCard*> *cards = new vector<PdbConnectCard*>;
    int full_cards = bonded_atoms.size()/kMaxBondedAtoms;
    for (int i = 0; i < full_cards; i++) {
        PdbConnectCard *card = new PdbConnectCard(source);
        for (int j = 0; j < kMaxBondedAtoms; j++) {
            card->add_bonded_atom(bonded_atoms[i*kMaxBondedAtoms + j]);
        }
        cards->push_back(card);
    }
    if (bonded_atoms.size()%kMaxBondedAtoms != 0) {
        PdbConnectCard *card = new PdbConnectCard(source);
        for (int i = full_cards*kMaxBondedAtoms; i < bonded_atoms.size(); i++) {
            card->add_bonded_atom(bonded_atoms[i]);
        }
        cards->push_back(card);
    }
    return cards;
}

void PdbConnectCard::add_bonded_atom(int serial) {
    if (bonded_atoms_.size() < kMaxBondedAtoms) {
        bonded_atoms_.push_back(serial);
    } else {
        throw std::runtime_error(get_too_many_bonds_error());
    }
}

void PdbConnectCard::write(std::ostream& out) const {
    string str(80, ' ');
    set_in_string(str, "CONECT", 0, 6, 'L');

    int cur_index = 6;
    set_in_string(str, source_, cur_index, kItemWidth, 'R');

    cur_index += kItemWidth;
    for (int i = 0; i < bonded_atoms_.size(); i++) {
        set_in_string(str, bonded_atoms_[i], cur_index, kItemWidth, 'R');
        cur_index += kItemWidth;
    }

    out << str;
}

void PdbConnectCard::read(const std::string& line) {
    string input(line);
    if (input.size() < 80)
        input.resize(80, ' ');

    int cur_index = 6;
    source_ = convert_string<int>(input.substr(cur_index, kItemWidth));

    cur_index += kItemWidth;
    for (int i = 0; i < kMaxBondedAtoms; i++) {
        string data = input.substr(cur_index, kItemWidth);
        trim(data);
        if (data.empty()) {
            break;
        }
        bonded_atoms_.push_back(convert_string<int>(data));
        cur_index += kItemWidth;
    }
}

std::string PdbConnectCard::get_too_many_bonds_error() {
    return make_string() << "CONECT cards may only have " << kMaxBondedAtoms <<
                            " bonded atoms.";
}

void PdbEndCard::write(std::ostream& out) const {
    out << "END";
}

void PdbLinkCard::read(const std::string& line) {
    string input(line);
    if (input.size() < 80)
        input.resize(80, ' ');
    name1_ = input.substr(12, 4);
    trim(name1_);
    res_name1_ = input.substr(17, 3);
    trim(res_name1_);
    res_seq1_ = convert_string<int>(input.substr(22, 4));
    name2_ = input.substr(42, 4);
    trim(name2_);
    res_name2_ = input.substr(47, 3);
    trim(res_name2_);
    res_seq2_ = convert_string<int>(input.substr(52, 4));
}

void PdbLinkCard::write(std::ostream& out) const {
    string str(80, ' ');
    set_in_string(str, "LINK", 0, 6, 'L');
    if (name1_.size() < 4)
        set_in_string(str, name1_, 13, 3, 'L');
    else
        set_in_string(str, name1_, 12, 4, 'L');
    set_in_string(str, res_name1_, 17, 3, 'L');
    set_in_string(str, res_seq1_, 22, 4, 'R');
    if (name2_.size() < 4)
        set_in_string(str, name2_, 43, 3, 'L');
    else
        set_in_string(str, name2_, 42, 4, 'L');
    set_in_string(str, res_name2_, 47, 3, 'L');
    set_in_string(str, res_seq2_, 52, 4, 'R');
    trim(str);
    out << str;
}

void PdbUnknownCard::write(std::ostream& out) const {
    out << line_;
}

void PdbEndMdlCard::write(std::ostream& out) const {
    out << "ENDMDL";
}

}  // namespace gmml
