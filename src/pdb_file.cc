// Copyright (c) 2012 The University of Georgia
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// Author: Robert Davis

#include "gmml/internal/pdb_file.h"

#include <cassert>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "gmml/internal/element.h"
#include "gmml/internal/environment.h"
#include "gmml/internal/geometry.h"
#include "utilities.h"
#include "gmml/internal/pdb_file_structure.h"

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
    string first_six = data_.substr(0, 6);
    first_six.resize(6, ' ');
    if (first_six == "ATOM  " || first_six == "HETATM") {
        return new PdbAtomCard(*this);
    } else if (first_six == "TER   ") {
        return new PdbTerCard(*this);
    } else if (first_six == "CONECT") {
        return new PdbConnectCard(*this);
    } else if (first_six == "END   ") {
        return new PdbEndCard(*this);
    } else if (first_six == "LINK  ") {
        return new PdbLinkCard(*this);
    } else if (first_six == "ENDMDL") {
        return new PdbEndMdlCard();
    } else if (first_six == "SEQRES") {
        return new PdbSeqresCard(*this);
    } else if (first_six == "MODRES") {
        return new PdbModelCard(*this);
    } else if (first_six == "SSBOND") {
        return new PdbSsbondCard(*this);
    } else if (first_six == "SITE  ") {
        return new PdbSiteCard(*this);
    } else if (first_six == "MODEL ") {
        return new PdbModelCard(*this);
    } else {
        return new PdbUnknownCard(*this);
    }
}

PdbAtomCardBuilder::PdbAtomCardBuilder()
        : serial_(-1), name_(""), alt_loc_(' '), res_name_(""), chain_id_(' '),
          res_seq_(-1), i_code_(' '), coordinate_(0.0, 0.0, 0.0),
          occupancy_(1.0), temp_factor_(0.0), element_(),
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
    set_in_string(str, element_.symbol(), 76, 2, 'R');
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
    try {
        element_ = Element(element_str);
    } catch(const std::invalid_argument& e) {
        element_ = Element();
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

std::vector<PdbSeqresCard*> PdbSeqresCardBuilder::build() const {
    return PdbSeqresCard::create_cards(*this);
}

void PdbSeqresCardBuilder::add_residue(std::string residue) {
    if (residue.length() > 3)
        throw std::invalid_argument("Residue string length must be at most 3 characters.");
    else
        residues_.push_back(residue);
}

bool PdbSeqresCardBuilder::validate() const {
    /**
    Check here to see if the PdbSeqresCardBuilder is
    fully validated
    */
}

std::vector<PdbSeqresCard*> PdbSeqresCard::create_cards(const PdbSeqresCardBuilder& builder) {
    std::vector<PdbSeqresCard*> seqres_vector;
    int serial_number = 1;
    PdbSeqresCard* current_card = new PdbSeqresCard(serial_number, builder.chain_id(), builder.size());
    for (int i = 0; i < builder.size(); i++) {
        if (i % kMaxNumberOfResidues == 0) {
            seqres_vector.push_back(current_card);
            serial_number++;
            PdbSeqresCard* current_card = new PdbSeqresCard(serial_number, builder.chain_id(), builder.size());
        }
        current_card->residues_.push_back(builder.at(i));
        if (i == builder.size() - 1)
            seqres_vector.push_back(current_card);
    }
    return seqres_vector;

}


void PdbSeqresCard::write(std::ostream& out) const {
    string str(80, ' ');
    set_in_string(str, "SEQRES", 0, 6, 'L');
    set_in_string(str, serial_number_, 7, 3, 'R');
    set_in_string(str, chain_id_, 11, 1, 'L');
    set_in_string(str, number_of_chain_residues_, 13, 4, 'R');
    int start_index = 18;
    for (int i = 0; i < residues_size(); i++) {
        set_in_string(str, residues_[i], start_index, 4, 'R');
        start_index++;
    }
    trim(str);
    out << str;
}


void PdbSeqresCard::read(const string& line) {
    string input(line);
    serial_number_ = convert_string<int>(input.substr(7, 3));
    chain_id_ = input[11];
    number_of_chain_residues_ = convert_string<int>(input.substr(13, 4));
    int residue_index = 0;
    for (int i = 19; i < line.size(); i = i + 4) {
        residues_[residue_index] = input.substr(i, 4);
        residue_index++;
    }
}

void PdbModresCard::write(std::ostream& out) const {
    string str(80, ' ');
    set_in_string(str, "MODRES", 0, 6, 'L');
    set_in_string(str, id_code_, 7, 4, 'R');
    set_in_string(str, res_name_, 12, 3, 'R');
    set_in_string(str, chain_id_, 16, 1, 'L');
    set_in_string(str, seq_num_, 18, 4, 'R');
    set_in_string(str, i_code_, 22, 1, 'L');
    set_in_string(str, std_res_name_, 24, 3, 'R');
    set_in_string(str, comment_, 29, 79, 'R');
    trim(str);
    out << str;
}

void PdbModresCard::read(const string& line) {
    string input(line);
    id_code_ = input.substr(7, 4);
    trim(id_code_);
    res_name_ = input.substr(12, 3);
    trim(res_name_);
    chain_id_ = input[16];
    seq_num_ = convert_string<int>(input.substr(18, 4));
    i_code_ = input[22];
    std_res_name_ = input.substr(24, 3);
    trim(std_res_name_);
    if (input.size() > 29)
        comment_ = input.substr(29, (input.size() - 29));
    trim(comment_);
}

void PdbSsbondCard::read(const string& line) {
    string input(line);
    ser_num_ = convert_string<int>(input.substr(7, 3));
    res_name_1_ = input.substr(11, 3);
    trim(res_name_1_);
    chain_id_1_ = input[15];
    res_seq_num_1_ = convert_string<int>(input.substr(17, 4));
    i_code_1_ = input[21];
    res_name_2_ = input.substr(25, 3);
    trim(res_name_2_);
    chain_id_2_ = input[29];
    res_seq_num_2_ = convert_string<int>(input.substr(31, 4));
    i_code_2_ = input[35];
    sym_op_1_ = convert_string<int>(input.substr(59, 6));
    sym_op_2_ = convert_string<int>(input.substr(66, 6));
    length_ = convert_string<double>(input.substr(73, 5));
}

void PdbSsbondCard::write(std::ostream& out) {
    string str(80, ' ');
    set_in_string(str, "SSBOND", 0, 6, 'L');
    set_in_string(str, ser_num_, 7, 3, 'R');
    set_in_string(str, res_name_1_, 11, 3, 'L');
    set_in_string(str, chain_id_1_, 15, 1, 'L');
    set_in_string(str, res_seq_num_1_, 17, 4, 'R');
    set_in_string(str, i_code_1_, 21, 1, 'L');
    set_in_string(str, res_name_2_, 25, 3, 'L');
    set_in_string(str, chain_id_2_, 29, 1, 'L');
    set_in_string(str, res_seq_num_2_, 31, 4, 'R');
    set_in_string(str, i_code_2_, 35, 1, 'L');
    set_in_string(str, sym_op_1_, 59, 6, 'R');
    set_in_string(str, sym_op_2_, 66, 6, 'R');
    set_in_string(str, length_, 73, 5, 'R');
    trim(str);
    out << str;
}

void PdbModelCard::write(std::ostream& out) const {
    string str(80, ' ');
    set_in_string(str, "MODEL ", 0, 6, 'L');
    set_in_string(str, ser_num_, 10, 4, 'R');
    trim(str);
    out << str;
}

void PdbModelCard::read(const std::string& line) {
    string input(line);
    ser_num_ = convert_string<int>(input.substr(10, 4));
}

vector<PdbSiteCard*> PdbSiteCardBuilder::build() const {
    return PdbSiteCard::create_cards(*this);
}

void PdbSiteCardBuilder::add_residue(const NamedPdbResidueId& residue) {
    residues_.push_back(new NamedPdbResidueId(residue));
}

vector<PdbSiteCard*> PdbSiteCard::create_cards(const PdbSiteCardBuilder& builder) {
    vector<PdbSiteCard*> site_vector;
    int sequence_number = 1;
    PdbSiteCard* current_card = new PdbSiteCard(sequence_number, builder.site_name(), builder.size());
    for (int i = 0; i < builder.size(); i++) {
        if (i % kMaxNumberOfResidues == 0) {
            site_vector.push_back(current_card);
            sequence_number++;
            PdbSiteCard* current_card = new PdbSiteCard(sequence_number, builder.site_name(), builder.size());
        }
        current_card->add_residue(*builder.at(i));
        if (i == builder.size() - 1)
            site_vector.push_back(current_card);
    }
    return site_vector;
}

void PdbSiteCard::add_residue(const NamedPdbResidueId& residue) {
    residues_.push_back(new NamedPdbResidueId(residue));
}

void PdbSiteCard::write(std::ostream& out) const {

}

void PdbSiteCard::read(const std::string& line) {

}

}  // namespace gmml
