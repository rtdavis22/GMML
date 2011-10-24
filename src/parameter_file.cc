#include "gmml/internal/parameter_file.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "gmml/internal/environment.h"
#include "gmml/internal/utilities.h"

using std::istringstream;
using std::map;
using std::set;
using std::string;
using std::vector;

namespace gmml
{

ParameterFile::ParameterFile() { 
    read(std::cin); 
}

void ParameterFile::read(const string& file_name) {
    std::ifstream stream(find_file(file_name).c_str());
    read(stream);
    stream.close();
}

void ParameterFile::read(std::istream& in) {
    string line;
    int line_index = 0;

    if (!getline(in, line)) {
        throw ParameterFileProcessingException("Error reading file");
    }
    title_ = trim(line);
    line_index++;

    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_atom_type(line);
        } catch (...) {
            throw ParameterFileProcessingException(line_index,
                    "Error processing atom type");
        }
    }

    while (line_index++, getline(in, line) && line.find('-') == string::npos) {
        try {
            process_hydrophilic(line);
        } catch (...) {
            throw ParameterFileProcessingException(line_index,
                    "Error processing hydrophilic atom list");
        }
    }

    do {
        try {
            process_bond(line);
        } catch (...) {
            throw ParameterFileProcessingException(line_index, 
                    "Error processing bond");
        }
    } while (line_index++, getline(in, line) && !trim(line).empty());

    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_angle(line);
        } catch (...) {
            throw ParameterFileProcessingException(line_index,
                    "Error processing angle");
        }
    }

    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_dihedral(line, line_index, in);
        } catch (...) {
            throw ParameterFileProcessingException(line_index, 
                    "Error processing dihedral");
        }
    }

    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_improper_dihedral(line);
        } catch (...) {
            throw ParameterFileProcessingException(line_index, 
                    "Error processing improper dihedral");
        }
    }

    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_hbond_parameters(line);
        } catch (...) {
            throw ParameterFileProcessingException(line_index, 
                    "Error processing h-bond parameters");
        }
    }

    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_equivalent_symbol_list(line);
        } catch (...) {
            throw ParameterFileProcessingException(line_index, 
                    "Error processing equivalent symbol list");
        }
    }

    line_index++; getline(in, line);
    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_6_12_parameters(line);
        } catch (...) {
            throw ParameterFileProcessingException(line_index, 
                    "Error processing 6-12 parameters");
        }
    }
}

void ParameterFile::process_atom_type(const string& line) {
    double mass, polarizability;
    string type, description;

    type = line.substr(0, 2);
    trim(type);
    mass = convert_string<double>(line.substr(3, 10));
    try {
        polarizability = convert_string<double>(line.substr(14, 10));
    } catch (const ConversionException&) {
        polarizability = kNotSet;
    }
    if (line.size() > 24)
       description = line.substr(24);

    atom_types_[type] = new ParameterFileAtom(type, mass, polarizability);
}

void ParameterFile::process_hydrophilic(const string& line) {
    string type;
    istringstream in(line);
    while (in >> std::setw(4) >> type && !trim(type).empty())
        hydrophilic_types_.push_back(type);
}

void ParameterFile::process_bond(const string& line) {
    char c;
    vector<string> types(2);
    string description;
    double force_constant, length;

    istringstream in(line);
    in >> std::setw(2) >> types[0] >> c >> std::setw(2) >> types[1] >>
          std::setw(10) >> force_constant >> length;

    if (in.fail())
        throw std::exception();
    if (line.size() > 26)
        description = line.substr(26);

    bonds_.insert(new ParameterFileBond(types, force_constant, length));
}

void ParameterFile::process_angle(const string& line) {
    char c;
    vector<string> types(3);
    string description;
    double force_constant, measure;

    istringstream in(line);
    in >> std::setw(2) >> types[0] >> c >> std::setw(2) >> types[1] >> c >>
          std::setw(2) >> types[2];
    in >> std::setw(10) >> force_constant >> std::setw(10) >> measure;
    if (in.fail())
        throw std::exception();
    if (line.size() > 29)
        description = line.substr(29);

    angles_.insert(new ParameterFileAngle(types, force_constant, measure));
}

void ParameterFile::process_dihedral(string& line, int& line_index, 
                                     std::istream& in) {
    char c;
    vector<string> types(4);
    string description;
    vector<ParameterFileDihedralTerm> terms;
    ParameterFileDihedralTerm t;
    double scee, scnb;

    istringstream ss(line);
    ss >> std::setw(2) >> types[0] >> c >> std::setw(2) >> types[1] >> c >>
          std::setw(2) >> types[2] >> c >> std::setw(2) >> types[3];
    ss >> std::setw(4) >> t.factor >> std::setw(15) >> t.force_constant >> 
          std::setw(15) >> t.phase >> std::setw(15) >> t.periodicity;
    if (ss.fail())
        throw std::exception();
    terms.push_back(t);

    if (line.size() > 60)
        description = line.substr(60);

    scee = extract_kv_double(description, "SCEE");
    scnb = extract_kv_double(description, "SCNB");

    ParameterFileDihedral *dihedral = 
        new ParameterFileDihedral(types, t, scee, scnb);
    while (dihedral->terms.at(dihedral->terms.size() - 1).periodicity < 0) {
        getline(in, line); line_index++;
        istringstream ss2(line.substr(11));
        ParameterFileDihedralTerm new_term;
        // make this a separate function probably
        ss2 >> std::setw(4) >> new_term.factor >> 
               std::setw(15) >> new_term.force_constant >>
               std::setw(15) >> new_term.phase >> 
               std::setw(15) >> new_term.periodicity;
        if (ss2.fail()) {
            throw ParameterFileProcessingException(line_index, 
                    "Error processing dihedral term");
        }
        dihedral->add_term(new_term);
    }

    if (dihedral->types[0] == "X" && dihedral->types[3] == "X")
        generic_dihedrals_.insert(dihedral);
    else
        dihedrals_.insert(dihedral);
}

void ParameterFile::process_improper_dihedral(const string& line) {
    char c;
    vector<string> atom_types(4);
    string description, s;
    ParameterFileDihedralTerm t;

    istringstream ss(line);
    ss >> std::setw(2) >> atom_types[0] >> c >> 
          std::setw(2) >> atom_types[1] >> c >> 
          std::setw(2) >> atom_types[2] >> c >> 
          std::setw(2) >> atom_types[3] >>
          std::setw(19) >> t.force_constant >> 
          std::setw(15) >> t.phase >> 
          std::setw(15) >> t.periodicity;
    if (ss.fail())
        throw std::exception();

    t.factor = kNotSet;

    if (line.size() > 60)
        description = line.substr(60);

    string center = atom_types[2];
    atom_types.erase(atom_types.begin() + 2);

    improper_dihedrals_.insert(
        new ParameterFileImproperDihedral(center, atom_types, t)
    );
}

void ParameterFile::process_hbond_parameters(const string& line) {
    string t1, t2, description;
    double A, B;

    istringstream ss(line);
    ss >> t1 >> t2 >> A >> B;
    if (ss.fail())
        throw std::exception();
    if (line.size() > 58)
        description = line.substr(58);

    HbondParameter h = { t1, t2, A, B, description };
    hbond_parameters_.push_back(h);
}

void ParameterFile::process_equivalent_symbol_list(const string& line) {
    string type, t;
    istringstream in(line);
    in >> std::setw(4) >> type;
    while (in >> std::setw(4) >> t && !trim(t).empty())
        equivalent_symbol_lists_[type].push_back(t);
}

void ParameterFile::process_6_12_parameters(const string& line) {
    string type, description;
    double radius, depth;
    istringstream in(line);
    in >> type >> radius >> depth;
    if (in.fail())
        throw std::exception();
    if (line.size() > 38)
        description = line.substr(38);

    trim(type);
    if (atom_types_.find(type) == atom_types_.end()) {
        atom_types_[type] = 
            new ParameterFileAtom(type, kNotSet, kNotSet, radius, depth);
    }
    else {
        atom_types_[type]->radius = radius;
        atom_types_[type]->well_depth = depth;
    }
    vector<string>::const_iterator it;
    vector<string> equivalent_atoms;
    if (equivalent_symbol_lists_.find(type) != equivalent_symbol_lists_.end())
        equivalent_atoms = equivalent_symbol_lists_[type];
    else
        return;

    for (it = equivalent_atoms.begin(); it != equivalent_atoms.end(); ++it) {
        if (atom_types_.find(*it) == atom_types_.end()) {
            atom_types_[*it] = 
                new ParameterFileAtom(*it, kNotSet, kNotSet, radius, depth);
        }
        else {
            atom_types_[*it]->radius = radius;
            atom_types_[*it]->well_depth = depth;
        }
    }
}                                                

double ParameterFile::extract_kv_double(const string& str, const string& key) {
    double val;
    size_t pos = str.find(string(key + "="));
    if (pos == string::npos)
        return kNotSet;
    std::istringstream ss(str.substr(pos + key.size() + 1));
    if ((ss >> val).fail())
        return kNotSet;
    return val;
}

void ParameterFileSet::load(const ParameterFile& file) {
    load_atom_types(file.atom_types().begin(), file.atom_types().end());
    load_bonds(file.bonds().begin(), file.bonds().end());
    load_angles(file.angles().begin(), file.angles().end());
    load_dihedrals(file.generic_dihedrals().begin(),
                   file.generic_dihedrals().end());
    load_dihedrals(file.dihedrals().begin(), file.dihedrals().end());
    load_improper_dihedrals(file.improper_dihedrals());
}

void ParameterFileSet::load_atom_types(
        ParameterFile::AtomTypeMap::const_iterator first,
        ParameterFile::AtomTypeMap::const_iterator last) {
    ParameterFile::AtomTypeMap::iterator it;
    while (first != last) {
        string type = first->first;
        if ((it = atom_types_.find(type)) == atom_types_.end()) {
            atom_types_[type] = new ParameterFileAtom(*first->second);
            first++;
            continue;
        }
        if (is_set(first->second->mass)) {
            warning("atom type " + it->second->type + " redefined");
            if (!is_set(it->second->mass)) {
                atom_types_[type]->mass = first->second->mass;
            }
            else if (it->second->mass != first->second->mass) {
                atom_types_[type]->mass = first->second->mass;
                warning("overriding mass of atom type " + type);
            }
        }
        if (is_set(first->second->polarizability)) {
            if (!is_set(it->second->polarizability)) {
                atom_types_[type]->polarizability = 
                    first->second->polarizability;
            }
            else if (it->second->polarizability !=
                    first->second->polarizability) {
                atom_types_[type]->polarizability =
                    first->second->polarizability;
                warning("overriding polarizability of atom " + type);
            }
       }
       if (is_set(first->second->radius)) {
            if (!is_set(it->second->radius)) {
                atom_types_[type]->radius = first->second->radius;
            }
            else if (it->second->radius != first->second->radius) {
                atom_types_[type]->radius = first->second->radius;
                warning("overriding radius of atom " + type);
            }
        }
        if (is_set(first->second->well_depth)) {
            if (!is_set(it->second->well_depth)) {
                atom_types_[type]->well_depth = first->second->well_depth;
            }
            else if (it->second->well_depth != first->second->well_depth) {
                atom_types_[type]->well_depth = first->second->well_depth;
                warning("overriding well depth of atom " + type);
            }
        }
        ++first;
    }
}

void ParameterFileSet::load_bonds(ParameterFile::BondSet::const_iterator first,
                                  ParameterFile::BondSet::const_iterator last) {
    ParameterFile::BondSet::iterator it;
    while (first != last) {
        ParameterFileBond *bond = new ParameterFileBond(**first);
        if ((it = find(bond)) == bonds_.end()) {
            bonds_.insert(bond);
            first++;
            continue;
        }
        if (bond->force_constant != (*it)->force_constant)
            warning("overriding force constant of bond " + (*it)->types[0] +
                    "-" + (*it)->types[1]);
        if (bond->length != (*it)->length)
            warning("overriding length of bond " + (*it)->types[0] + "-" +
                    (*it)->types[1]);

        delete *it;
        bonds_.erase(it);
        bonds_.insert(bond);
        first++;
    }
}

void ParameterFileSet::load_angles(
        ParameterFile::AngleSet::const_iterator first,
        ParameterFile::AngleSet::const_iterator last) {
    ParameterFile::AngleSet::iterator it;
    while (first != last) {
        ParameterFileAngle *angle = new ParameterFileAngle(**first);
        if ((it = find(angle)) == angles_.end()) {
            angles_.insert(angle);
            first++;
            continue;
        }
        warning("angle " + (*it)->types[0] + "-" + (*it)->types[1] + "-" + 
                (*it)->types[2] + " redefined");
        if (angle->force_constant != (*it)->force_constant)
            warning("overriding force constant of angle " + (*it)->types[0] +
                    "-" + (*it)->types[1] + "-" + (*it)->types[2]);
        if (angle->angle != (*it)->angle)
            warning("overriding measure of angle " + (*it)->types[0] + "-" +
                    (*it)->types[1] + "-" + (*it)->types[2]);

        delete *it;
        angles_.erase(it);
        angles_.insert(angle);
        first++;
    }
}

void ParameterFileSet::load_dihedrals(
        ParameterFile::DihedralSet::const_iterator first,
        ParameterFile::DihedralSet::const_iterator last) {
    ParameterFile::DihedralSet::iterator it;
    while (first != last) {
        ParameterFileDihedral *dihedral = new ParameterFileDihedral(**first);
        ParameterFileDihedral *reverse_dihedral = 
            new ParameterFileDihedral(*dihedral);
        swap(reverse_dihedral->types[0], reverse_dihedral->types[3]);
        swap(reverse_dihedral->types[1], reverse_dihedral->types[2]);

        if (dihedral->types[0] == "X" && dihedral->types[3] == "X") {
            if (generic_dihedrals_.find(dihedral) == generic_dihedrals_.end() &&
                    generic_dihedrals_.find(reverse_dihedral) ==
                    generic_dihedrals_.end()) {
                generic_dihedrals_.insert(dihedral);
            }
            else {
                warning("redefinition of dihedral " + dihedral->types[0] + "-" +
                        dihedral->types[1] + "-" + dihedral->types[2] + "-" +
                        dihedral->types[3]);
                //overwrite it
            }
        }
        else {
            it = dihedrals_.find(dihedral);
            if (it == dihedrals_.end())
                it = dihedrals_.find(reverse_dihedral);
            //the dihedral doesn't exist
            if (it != dihedrals_.end()) {
                warning("overriding dihedral " + (*it)->types[0] + "-" +
                        (*it)->types[1] + "-" + (*it)->types[2] + "-" +
                        (*it)->types[3]);
                delete *it;
                dihedrals_.erase(it);
            }
            dihedrals_.insert(dihedral);
        }
        delete reverse_dihedral;
        ++first;
    }
}

ParameterFile::BondSet::const_iterator ParameterFileSet::find(
        const string& type1, const string& type2) const {
    ParameterFileBond *bond = new ParameterFileBond(type1, type2);
    ParameterFile::BondSet::iterator it;
    if ((it = bonds_.find(bond)) == bonds_.end()) {
        delete bond;
        bond = new ParameterFileBond(type2, type1);
        it = bonds_.find(bond);
    }
    delete bond;
    return it;
}

ParameterFile::AngleSet::const_iterator ParameterFileSet::find(
        const std::string& type1, const std::string& type2,
        const std::string& type3) const {
    ParameterFileAngle *angle = new ParameterFileAngle(type1, type2, type3);
    ParameterFile::AngleSet::iterator it;
    if ((it = angles_.find(angle)) == angles_.end()) {
        delete angle;
        angle = new ParameterFileAngle(type3, type2, type1);
        it = angles_.find(angle);
    }
    delete angle;
    return it;
}

ParameterFile::DihedralSet::const_iterator ParameterFileSet::find(
        const std::string& type1, const std::string& type2,
        const std::string& type3, const std::string& type4) const {
    ParameterFileDihedral *dihedral =
        new ParameterFileDihedral(type1, type2, type3, type4);
    ParameterFile::DihedralSet::iterator it;
    if ((it = dihedrals_.find(dihedral)) == dihedrals_.end()) {
        delete dihedral;
        dihedral = new ParameterFileDihedral(type4, type3, type2, type1);
        it = dihedrals_.find(dihedral);
    }
    delete dihedral;
    return it;
}

ParameterFile::DihedralSet::const_iterator ParameterFileSet::find_generic(
        const std::string& type1, const std::string& type2,
        const std::string& type3, const std::string& type4) const {
    ParameterFileDihedral *dihedral = 
        new ParameterFileDihedral(type1, type2, type3, type4);
    ParameterFile::DihedralSet::iterator it;
    if ((it = generic_dihedrals_.find(dihedral)) == generic_dihedrals_.end()) {
        delete dihedral;
        dihedral = new ParameterFileDihedral(type4, type3, type2, type1);
        it = generic_dihedrals_.find(dihedral);
    }
    delete dihedral;
    return it;
}

} //namespace gmml
