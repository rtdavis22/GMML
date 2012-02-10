#include "gmml/internal/parameter_file.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>

#include "gmml/internal/environment.h"
#include "gmml/internal/stubs/logging.h"
#include "utilities.h"

using std::istringstream;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

namespace gmml {

void ParameterFileImproperDihedral::init() {
    types.erase(std::remove(types.begin(), types.end(), "X"),
                types.end());
    std::sort(types.begin(), types.end());
}

struct ParameterFileBondPtrLess
        : public std::binary_function<const ParameterFileBond*,
                                      const ParameterFileBond*, bool> {
    bool operator()(const ParameterFileBond *lhs,
                    const ParameterFileBond *rhs) const {
        return lhs->types < rhs->types;
    }
};

struct ParameterFileAnglePtrLess
        : public std::binary_function<const ParameterFileAngle*,
                                      const ParameterFileAngle*, bool> {
    bool operator()(const ParameterFileAngle *lhs,
                    const ParameterFileAngle *rhs) const {
        return lhs->types < rhs->types;
    }
};

struct ParameterFileDihedralPtrLess
        : std::binary_function<const ParameterFileDihedral*,
                               const ParameterFileDihedral*, bool> {
    bool operator()(const ParameterFileDihedral *lhs,
                    const ParameterFileDihedral *rhs) const {
        return lhs->types < rhs->types;
    }
};

bool ParameterFileImproperDihedral::matches(const vector<string>& rhs) const {
    vector<string> rhs_types(rhs);
    std::sort(rhs_types.begin(), rhs_types.end());
    // No Xs
    if (types.size() == 3) {
        return rhs_types[0] == types[0] && rhs_types[1] == types[1] &&
               rhs_types[2] == types[2];
    } else if (types.size() == 2) {
        // One X
        return (rhs_types[0] == types[0] && rhs_types[1] == types[1]) ||
               (rhs_types[1] == types[0] && rhs_types[2] == types[1]) ||
               (rhs_types[0] == types[0] && rhs_types[2] == types[1]);
    } else if (types.size() == 1) {
        // Two Xs
        return rhs_types[0] == types[0] || rhs_types[1] == types[0] ||
               rhs_types[2] == types[0];
    } else {
        // All Xs
        return true;
    }
}

std::pair<ParameterFileDihedralTerm, bool> ImproperDihedralCollection::lookup(
        const string& center, const vector<string>& types) const {
    ParameterFileDihedralTerm term;
    const_iterator it = dihedrals_.find(center);
    if (it != dihedrals_.end()) {
        for (int i = it->second.size() - 1; i >= 0; i--) {
            if (it->second[i]->matches(types))
                return std::make_pair(it->second[i]->term, true);
        }
    }
    return std::make_pair(term, false);
}

void ImproperDihedralCollection::append(
        const ImproperDihedralCollection& coll) {
    vector<ParameterFileImproperDihedral*>::const_iterator it;
    const_iterator coll_it;
    for (coll_it = coll.begin(); coll_it != coll.end(); ++coll_it) {
        for (it = coll_it->second.begin(); it != coll_it->second.end(); ++it) {
            dihedrals_[coll_it->first].push_back(
                new ParameterFileImproperDihedral(**it));
        }
    }
}

const char *ParameterFileProcessingException::what() const throw() {
    what_ = "ParameterFile: " + message_;
    if (line_number_ != kNotSet)
        what_ += " (line " + to_string(line_number_) + ")";
    return what_.c_str();
}

// Private implementation
class ParameterFile::Impl {
  public:
    explicit Impl(const string& file) { read(file); }
    ~Impl();

    const AtomTypeMap& atom_types() const { return atom_types_; }
    const BondSet& bonds() const { return bonds_; }
    const AngleSet& angles() const { return angles_; }
    const DihedralSet& dihedrals() const { return dihedrals_; }
    const DihedralSet& generic_dihedrals() const { return generic_dihedrals_; }
    const ImproperDihedralCollection& improper_dihedrals() const {
        return improper_dihedrals_;
    }

  private:
    // These throw ParameterFileProcessingException.
    void read(const string& file);
    void read(std::istream&);

    void process_atom_type(const string& line);
    void process_hydrophilic(const string& line);
    void process_bond(const string& line);
    void process_angle(const string& line);
    void process_dihedral(string& line, int& line_index, std::istream&);
    void process_improper_dihedral(const string& line);
    void process_hbond_parameters(const string& line);
    void process_equivalent_symbol_list(const string& line);
    void process_6_12_parameters(const string& line);

    double extract_kv_double(const string& str, const string& key);

    string title_;
    AtomTypeMap atom_types_;
    vector<string> hydrophilic_types_;
    BondSet bonds_;
    AngleSet angles_;
    DihedralSet dihedrals_;
    DihedralSet generic_dihedrals_;
    ImproperDihedralCollection improper_dihedrals_;
    vector<HbondParameter> hbond_parameters_;
    map<string, vector<string> > equivalent_symbol_lists_;
};

ParameterFile::Impl::~Impl() {
    for (AtomTypeMap::iterator it = atom_types_.begin();
            it != atom_types_.end(); ++it)
        delete it->second;
    std::for_each(bonds_.begin(), bonds_.end(), DeletePtr());
    std::for_each(angles_.begin(), angles_.end(), DeletePtr());
    std::for_each(dihedrals_.begin(), dihedrals_.end(), DeletePtr());
    std::for_each(generic_dihedrals_.begin(), generic_dihedrals_.end(),
                  DeletePtr());
}

void ParameterFile::Impl::read(const string& file_name) {
    std::ifstream stream(find_file(file_name).c_str());
    read(stream);
    stream.close();
}

void ParameterFile::Impl::read(std::istream& in) {
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
        } catch(...) {
            throw ParameterFileProcessingException(line_index,
                    "Error processing atom type");
        }
    }

    while (line_index++, getline(in, line) && line.find('-') == string::npos) {
        try {
            process_hydrophilic(line);
        } catch(...) {
            throw ParameterFileProcessingException(line_index,
                    "Error processing hydrophilic atom list");
        }
    }

    do {
        try {
            process_bond(line);
        } catch(...) {
            throw ParameterFileProcessingException(line_index,
                    "Error processing bond");
        }
    } while (line_index++, getline(in, line) && !trim(line).empty());

    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_angle(line);
        } catch(...) {
            throw ParameterFileProcessingException(line_index,
                    "Error processing angle");
        }
    }

    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_dihedral(line, line_index, in);
        } catch(...) {
            throw ParameterFileProcessingException(line_index,
                    "Error processing dihedral");
        }
    }

    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_improper_dihedral(line);
        } catch(...) {
            throw ParameterFileProcessingException(line_index,
                    "Error processing improper dihedral");
        }
    }

    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_hbond_parameters(line);
        } catch(...) {
            throw ParameterFileProcessingException(line_index,
                    "Error processing h-bond parameters");
        }
    }

    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_equivalent_symbol_list(line);
        } catch(...) {
            throw ParameterFileProcessingException(line_index,
                    "Error processing equivalent symbol list");
        }
    }

    line_index++;
    getline(in, line);
    while (line_index++, getline(in, line) && !trim(line).empty()) {
        try {
            process_6_12_parameters(line);
        } catch(...) {
            throw ParameterFileProcessingException(line_index,
                    "Error processing 6-12 parameters");
        }
    }
}

void ParameterFile::Impl::process_atom_type(const string& line) {
    double mass, polarizability;
    string type, description;

    type = line.substr(0, 2);
    trim(type);
    mass = convert_string<double>(line.substr(3, 10));
    try {
        polarizability = convert_string<double>(line.substr(14, 10));
    } catch(const ConversionException&) {
        polarizability = kNotSet;
    }
    if (line.size() > 24)
       description = line.substr(24);

    atom_types_[type] = new ParameterFileAtom(type, mass, polarizability);
}

void ParameterFile::Impl::process_hydrophilic(const string& line) {
    string type;
    istringstream in(line);
    while (in >> std::setw(4) >> type && !trim(type).empty())
        hydrophilic_types_.push_back(type);
}

void ParameterFile::Impl::process_bond(const string& line) {
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

void ParameterFile::Impl::process_angle(const string& line) {
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

void ParameterFile::Impl::process_dihedral(string& line, int& line_index,
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
        getline(in, line);
        line_index++;
        istringstream ss2(line.substr(11));
        ParameterFileDihedralTerm new_term;
        // Make this a separate function probably.
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

void ParameterFile::Impl::process_improper_dihedral(const string& line) {
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
        new ParameterFileImproperDihedral(center, atom_types, t));
}

void ParameterFile::Impl::process_hbond_parameters(const string& line) {
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

void ParameterFile::Impl::process_equivalent_symbol_list(const string& line) {
    string type, t;
    istringstream in(line);
    in >> std::setw(4) >> type;
    while (in >> std::setw(4) >> t && !trim(t).empty())
        equivalent_symbol_lists_[type].push_back(t);
}

void ParameterFile::Impl::process_6_12_parameters(const string& line) {
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
    } else {
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
        } else {
            atom_types_[*it]->radius = radius;
            atom_types_[*it]->well_depth = depth;
        }
    }
}

double ParameterFile::Impl::extract_kv_double(const string& str,
                                              const string& key) {
    double val;
    size_t pos = str.find(string(key + "="));
    if (pos == string::npos)
        return kNotSet;
    std::istringstream ss(str.substr(pos + key.size() + 1));
    if ((ss >> val).fail())
        return kNotSet;
    return val;
}

// Public implementation
ParameterFile::ParameterFile(const string& file) : impl_(new Impl(file)) {}

ParameterFile::~ParameterFile() {}

const ParameterFile::AtomTypeMap& ParameterFile::atom_types() const {
    return impl_->atom_types();
}

const ParameterFile::BondSet& ParameterFile::bonds() const {
    return impl_->bonds();
}

const ParameterFile::AngleSet& ParameterFile::angles() const {
    return impl_->angles();
}

const ParameterFile::DihedralSet& ParameterFile::dihedrals() const {
    return impl_->dihedrals();
}

const ParameterFile::DihedralSet& ParameterFile::generic_dihedrals() const {
    return impl_->generic_dihedrals();
}

const ImproperDihedralCollection& ParameterFile::improper_dihedrals() const {
    return impl_->improper_dihedrals();
}

// Private implementation
class ParameterSet::Impl {
  public:
    Impl() {}

    void load(const ParameterFile& parameter_file);
    void load(const string& file_name) { load(ParameterFile(file_name)); }

    const ParameterFileAtom *lookup(const string& type) const;
    const ParameterFileBond *lookup(const string& type1,
                                    const string& type2) const;
    const ParameterFileAngle *lookup(const string& type1,
                                     const string& type2,
                                     const string& type3) const;
    const ParameterFileDihedral *lookup(const string& type1,
                                        const string& type2,
                                        const string& type3,
                                        const string& type4) const;
    pair<ParameterFileDihedralTerm, bool> lookup_improper_dihedral(
            const string& center_type,
            const string& type1,
            const string& type2,
            const string& type3) const;

  private:
    void load_atom_types(ParameterFile::AtomTypeMap::const_iterator first,
                         ParameterFile::AtomTypeMap::const_iterator last);
    void load_bonds(ParameterFile::BondSet::const_iterator first,
                    ParameterFile::BondSet::const_iterator last);
    void load_angles(ParameterFile::AngleSet::const_iterator first,
                     ParameterFile::AngleSet::const_iterator last);
    void load_dihedrals(ParameterFile::DihedralSet::const_iterator first,
                        ParameterFile::DihedralSet::const_iterator last);
    void load_improper_dihedrals(const ImproperDihedralCollection& coll) {
        improper_dihedrals_.append(coll);
    }

    ParameterFile::BondSet::const_iterator find(const string& type1,
                                                const string& type2) const;
    ParameterFile::BondSet::const_iterator find(
            const ParameterFileBond *bond) const {
        return find(bond->types[0], bond->types[1]);
    }
    ParameterFile::AngleSet::const_iterator find(
            const string& type1, const string& type2,
            const string& type3) const;
    ParameterFile::AngleSet::const_iterator find(
            const ParameterFileAngle *angle) const {
        return find(angle->types[0], angle->types[1], angle->types[2]);
    }
    ParameterFile::DihedralSet::const_iterator find(
            const string& type1, const string& type2,
            const string& type3, const string& type4) const;
    ParameterFile::DihedralSet::const_iterator find_generic(
            const string& type1, const string& type2,
            const string& type3, const string& type4) const;

    ParameterFile::AtomTypeMap atom_types_;
    ParameterFile::BondSet bonds_;
    ParameterFile::AngleSet angles_;
    ParameterFile::DihedralSet dihedrals_;
    ParameterFile::DihedralSet generic_dihedrals_;
    ImproperDihedralCollection improper_dihedrals_;
};

void ParameterSet::Impl::load(const ParameterFile& file) {
    load_atom_types(file.atom_types().begin(), file.atom_types().end());
    load_bonds(file.bonds().begin(), file.bonds().end());
    load_angles(file.angles().begin(), file.angles().end());
    load_dihedrals(file.generic_dihedrals().begin(),
                   file.generic_dihedrals().end());
    load_dihedrals(file.dihedrals().begin(), file.dihedrals().end());
    load_improper_dihedrals(file.improper_dihedrals());
}

const ParameterFileAtom *ParameterSet::Impl::lookup(
        const string& type) const {
    ParameterFile::AtomTypeMap::const_iterator it;
    if ((it = atom_types_.find(type)) != atom_types_.end()) {
        return it->second;
    }
    return NULL;
}

const ParameterFileBond *ParameterSet::Impl::lookup(
        const string& type1,
        const string& type2) const {
    ParameterFile::BondSet::const_iterator it;
    if ((it = find(type1, type2)) != bonds_.end()) {
        return *it;
    }
    return NULL;
}

const ParameterFileAngle *ParameterSet::Impl::lookup(
        const string& type1,
        const string& type2,
        const string& type3) const {
    ParameterFile::AngleSet::const_iterator it;
    if ((it = find(type1, type2, type3)) != angles_.end()) {
        return *it;
    }
    return NULL;
}

const ParameterFileDihedral *ParameterSet::Impl::lookup(
        const string& type1,
        const string& type2,
        const string& type3,
        const string& type4) const {
    ParameterFile::DihedralSet::const_iterator it;
    if ((it = find(type1, type2, type3, type4)) == dihedrals_.end()) {
        if ((it = find_generic("X", type2, type3, "X")) ==
                 generic_dihedrals_.end()) {
             return NULL;
        }
    }
    return *it;
}

std::pair<ParameterFileDihedralTerm, bool>
ParameterSet::Impl::lookup_improper_dihedral(
        const string& center_type,
        const string& type1,
        const string& type2,
        const string& type3) const {
    vector<string> types;
    types.reserve(3);
    types.push_back(type1);
    types.push_back(type2);
    types.push_back(type3);
    return improper_dihedrals_.lookup(center_type, types);
}

void ParameterSet::Impl::load_atom_types(
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
            LOG(WARNING) << "Atom type " + it->second->type + " redefined.";
            if (!is_set(it->second->mass)) {
                atom_types_[type]->mass = first->second->mass;
            } else if (it->second->mass != first->second->mass) {
                atom_types_[type]->mass = first->second->mass;
                LOG(WARNING) << "Overriding mass of atom type " << type << ".";
            }
        }
        if (is_set(first->second->polarizability)) {
            if (!is_set(it->second->polarizability)) {
                atom_types_[type]->polarizability =
                    first->second->polarizability;
            } else if (it->second->polarizability !=
                    first->second->polarizability) {
                atom_types_[type]->polarizability =
                    first->second->polarizability;
                LOG(WARNING) << "Overriding polarizability of atom " << type <<
                                ".";
            }
       }
       if (is_set(first->second->radius)) {
            if (!is_set(it->second->radius)) {
                atom_types_[type]->radius = first->second->radius;
            } else if (it->second->radius != first->second->radius) {
                atom_types_[type]->radius = first->second->radius;
                LOG(WARNING) << "Overriding radius of atom " << type << ".";
            }
        }
        if (is_set(first->second->well_depth)) {
            if (!is_set(it->second->well_depth)) {
                atom_types_[type]->well_depth = first->second->well_depth;
            } else if (it->second->well_depth != first->second->well_depth) {
                atom_types_[type]->well_depth = first->second->well_depth;
                LOG(WARNING) << "Overriding well depth of atom " << type << ".";
            }
        }
        ++first;
    }
}

void ParameterSet::Impl::load_bonds(
        ParameterFile::BondSet::const_iterator first,
        ParameterFile::BondSet::const_iterator last) {
    ParameterFile::BondSet::iterator it;
    while (first != last) {
        ParameterFileBond *bond = new ParameterFileBond(**first);
        if ((it = find(bond)) == bonds_.end()) {
            bonds_.insert(bond);
            first++;
            continue;
        }
        if (bond->force_constant != (*it)->force_constant) {
            LOG(WARNING) << "Overriding force constant of bond " <<
                            (*it)->types[0] << "-" << (*it)->types[1] << ".";
        }
        if (bond->length != (*it)->length) {
            LOG(WARNING) << "Overriding length of bond " << (*it)->types[0] <<
                            "-" << (*it)->types[1] << ".";
        }

        delete *it;
        bonds_.erase(it);
        bonds_.insert(bond);
        first++;
    }
}

void ParameterSet::Impl::load_angles(
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
        LOG(WARNING) << "Angle " << (*it)->types[0] << "-" << (*it)->types[1] <<
                        "-" << (*it)->types[2] << " redefined.";
        if (angle->force_constant != (*it)->force_constant) {
            LOG(WARNING) << "Overriding force constant of angle " <<
                            (*it)->types[0] << "-" << (*it)->types[1] << "-" <<
                            (*it)->types[2] << ".";
        }
        if (angle->angle != (*it)->angle) {
            LOG(WARNING) << "Overriding measure of angle " << (*it)->types[0] <<
                            "-" << (*it)->types[1] << "-" << (*it)->types[2] <<
                            ".";
        }

        delete *it;
        angles_.erase(it);
        angles_.insert(angle);
        first++;
    }
}

void ParameterSet::Impl::load_dihedrals(
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
            } else {
                LOG(WARNING) << "Redefinition of dihedral " <<
                                dihedral->types[0] << "-" <<
                                dihedral->types[1] << "-" <<
                                dihedral->types[2] << "-" <<
                                dihedral->types[3] << ".";
            }
        } else {
            it = dihedrals_.find(dihedral);
            if (it == dihedrals_.end()) {
                it = dihedrals_.find(reverse_dihedral);
            }
            // The dihedral doesn't exist.
            if (it != dihedrals_.end()) {
                LOG(WARNING) << "Overriding dihedral " << (*it)->types[0] <<
                                "-" << (*it)->types[1] << "-" <<
                               (*it)->types[2] << "-" << (*it)->types[3] << ".";
                delete *it;
                dihedrals_.erase(it);
            }
            dihedrals_.insert(dihedral);
        }
        delete reverse_dihedral;
        ++first;
    }
}

ParameterFile::BondSet::const_iterator ParameterSet::Impl::find(
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

ParameterFile::AngleSet::const_iterator ParameterSet::Impl::find(
        const string& type1, const string& type2,
        const string& type3) const {
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

ParameterFile::DihedralSet::const_iterator ParameterSet::Impl::find(
        const string& type1, const string& type2,
        const string& type3, const string& type4) const {
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

ParameterFile::DihedralSet::const_iterator ParameterSet::Impl::find_generic(
        const string& type1, const string& type2,
        const string& type3, const string& type4) const {
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

// Public implementation
ParameterSet::ParameterSet() : impl_(new Impl) {}

void ParameterSet::load(const ParameterFile& parameter_file) {
    impl_->load(parameter_file);
}

const ParameterFileAtom *ParameterSet::lookup(const string& type) const {
    return impl_->lookup(type);
}

const ParameterFileBond *ParameterSet::lookup(const string& type1,
                                              const string& type2) const {
    return impl_->lookup(type1, type2);
}

const ParameterFileAngle *ParameterSet::lookup(const string& type1,
                                               const string& type2,
                                               const string& type3) const {
    return impl_->lookup(type1, type2, type3);
}

const ParameterFileDihedral *ParameterSet::lookup(const string& type1,
                                                  const string& type2,
                                                  const string& type3,
                                                  const string& type4) const {
    return impl_->lookup(type1, type2, type3, type4);
}

std::pair<ParameterFileDihedralTerm, bool>
ParameterSet::lookup_improper_dihedral(const string& center_type,
                                           const string& type1,
                                           const string& type2,
                                           const string& type3) const {
    return impl_->lookup_improper_dihedral(center_type, type1, type2, type3);
}

}  // namespace gmml
