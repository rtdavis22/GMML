#include "gmml/internal/pdb_file.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "gmml/internal/environment.h"
#include "gmml/internal/utilities.h"
#include "gmml/internal/geometry.h"

namespace gmml
{

using std::istringstream;
using std::list;
using std::string;
using std::vector;

void PdbAtomCard::write(std::ostream& out) const {
    string str(80, ' ');
    set_in_string(str, "ATOM", 0, 6, 'L');
    set_in_string(str, serial, 6, 5, 'R');
    if (name.size() < 4)
        set_in_string(str, name, 13, 3, 'L');
    else
        set_in_string(str, name, 12, 4, 'L');
    str[16] = alt_loc;
    set_in_string(str, res_name, 17, 3, 'L');
    str[21] = chain_id;
    set_in_string(str, res_seq, 22, 4, 'R');
    str[26] = i_code;
    set_in_string(str, x, 30, 8, 'R', 3);
    set_in_string(str, y, 38, 8, 'R', 3);
    set_in_string(str, z, 46, 8, 'R', 3);
    if (is_set(occupancy))
        set_in_string(str, occupancy, 54, 6, 'R', 2);
    if (is_set(temp_factor))
        set_in_string(str, temp_factor, 60, 6, 'R', 2);
    set_in_string(str, element, 76, 2, 'R');
    if (is_set(charge))
        set_in_string(str, charge, 78, 2, 'R');
    trim(str);
    out << str;
}
void PdbAtomCard::read(const string& line) {
    string input(line);
    if (input.size() < 80)
	input.resize(80, ' ');
    serial = convert_string<int>(input.substr(6, 5));
    name = input.substr(12, 4);
    trim(name);
    alt_loc = input[16];
    res_name = input.substr(17, 3);
    trim(res_name);
    chain_id = input[21];
    res_seq = convert_string<int>(input.substr(22, 4));
    i_code = input[26];
    x = convert_string<double>(input.substr(30, 8));
    y = convert_string<double>(input.substr(38, 8));
    z = convert_string<double>(input.substr(46, 8));
    try {
        occupancy = convert_string<double>(input.substr(54, 6));
    } catch (const ConversionException& e) {
        occupancy = kNotSet;
    }
    try {
        temp_factor = convert_string<double>(input.substr(60, 6));
    } catch (const ConversionException e) {
        temp_factor = kNotSet;
    }
    element = input.substr(76, 2);
    trim(element);
    if (element == "" && name.size() > 0)
        element = string(1, name[0]);
    try {
        charge = convert_string<double>(input.substr(78, 2));
    } catch (const ConversionException& e) {
        charge = kNotSet;
    }
}

void PdbConnectCard::write(std::ostream& out) const {
    string str(80, ' ');
    set_in_string(str, "CONECT", 0, 6, 'L');
    if (is_set(connect1))
        set_in_string(str, connect1, 6, 5, 'R');
    if (is_set(connect2))
        set_in_string(str, connect2, 11, 5, 'R');
    if (is_set(connect3))
        set_in_string(str, connect3, 16, 5, 'R');
    if (is_set(connect4))
        set_in_string(str, connect4, 21, 5, 'R');
    if (is_set(connect5))
        set_in_string(str, connect5, 26, 5, 'R');
    trim(str);
    out << str;
}

void PdbConnectCard::read(const std::string& line) {
    string input(line);
    if (input.size() < 80)
        input.resize(80, ' ');

    connect1 = connect2 = connect3 = connect4 = connect5 = kNotSet;
    try {
        connect1 = convert_string<int>(input.substr(6, 5));
        connect2 = convert_string<int>(input.substr(11, 5));
        connect3 = convert_string<int>(input.substr(16, 5));
        connect4 = convert_string<int>(input.substr(21, 5));
        connect5 = convert_string<int>(input.substr(26, 5));
    } catch (const ConversionException& /*e*/) {
        // Do nothing. If this occurs, the remaining connects will all
        // properly be set to kNotSet.
    }
}

void PdbLinkCard::read(const std::string& line) {
    string input(line);
    if (input.size() < 80)
        input.resize(80, ' ');
    name1 = input.substr(12, 4);
    trim(name1);
    res_name1 = input.substr(17, 3);
    trim(res_name1);
    res_seq1 = convert_string<int>(input.substr(22, 4));
    name2 = input.substr(42, 4);
    trim(name2);
    res_name2 = input.substr(47, 3);
    trim(res_name2);
    res_seq2 = convert_string<int>(input.substr(52, 4));
}

void PdbLinkCard::write(std::ostream& out) const {
    string str(80, ' ');
    set_in_string(str, "LINK", 0, 6, 'L');
    if (name1.size() < 4)
        set_in_string(str, name1, 13, 3, 'L');
    else
        set_in_string(str, name1, 12, 4, 'L');
    set_in_string(str, res_name1, 17, 3, 'L');
    set_in_string(str, res_seq1, 22, 4, 'L');
    if (name2.size() < 4)
        set_in_string(str, name2, 43, 3, 'L');
    else
        set_in_string(str, name2, 42, 4, 'L');
    set_in_string(str, res_name2, 47, 3, 'L');
    set_in_string(str, res_seq2, 52, 4, 'L');
    trim(str);
    out << str;
}

void PdbFile::read(const string& file_name) {
    std::ifstream stream(find_file(file_name).c_str());
    read(stream);
    stream.close();
}

void PdbFile::print(const string& file) const {
    std::ofstream out;
    out.open(file.c_str());
    write(out);
    out.close();
}

void PdbFile::print() const {
    write(std::cout);
}

void PdbFile::read(std::istream& in) {
    string line;
    while (getline(in, line)) {
        string card_type = line.substr(0, 6);
        card_type.resize(6, ' ');
        boost::shared_ptr<PdbCard> card_ptr;
        switch (get_card_type(card_type)) {
          case ATOM:
          case HETATM: {
            boost::shared_ptr<PdbAtomCard> atom_card(new PdbAtomCard(line));
            card_ptr = atom_card;
            atom_cards_.push_back(atom_card);
            break;
          }
          case TER:
            card_ptr.reset(new PdbTerCard(line));
            break;
          case CONECT: {
            boost::shared_ptr<PdbConnectCard> connect_card(
                new PdbConnectCard(line)
            );
            card_ptr = connect_card;
            connect_cards_.push_back(connect_card);
            break;
          }
          case END:
            card_ptr.reset(new PdbEndCard(line));
            break;
          case LINK:
            card_ptr.reset(new PdbLinkCard(line));
            break;
          default:
            warning("Unknown card " + card_type);
            card_ptr.reset(new PdbUnknownCard(line));
            break;
        }
        cards_.push_back(card_ptr);
    }
}

PdbFile::CardType PdbFile::get_card_type(const string& card_name) {
    if (card_name == "ATOM  " || card_name == "HETATM")
        return PdbFile::ATOM;
    else if (card_name == "TER   ")
        return PdbFile::TER;
    else if (card_name == "CONECT")
        return PdbFile::CONECT;
    else if (card_name == "END   ")
        return PdbFile::END;
    else if(card_name == "LINK  ")
        return PdbFile::LINK;
    else
        return PdbFile::UNKNOWN;
}

void PdbFile::write(std::ostream& out) const {
    list<CardPtr>::const_iterator it;
    for (it = cards_.begin(); it != cards_.end(); ++it) {
        (*it)->write(out);
        out << std::endl;
    }
}

} //namespace gmml
