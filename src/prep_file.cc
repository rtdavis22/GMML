#include "gmml/internal/prep_file.h"

#include <cmath>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stack>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "gmml/internal/environment.h"
#include "gmml/internal/geometry.h"
#include "gmml/internal/utilities.h"

namespace gmml {

using std::map;
using std::stack;
using std::string;
using std::vector;

PrepFile::PrepFile() { read(std::cin); }

void PrepFile::read(const string& file_name) {
    std::ifstream stream(find_file(file_name).c_str());
    read(stream);
    stream.close();
}

void PrepFile::read(std::istream& in) {
    string line;
    getline(in, line);
    getline(in, line);
    while (process_residue(in)) {}
}

bool PrepFile::process_residue(std::istream& in) {
    ResiduePtr residue(new PrepFileResidue);
    string line;
    std::istringstream ss;

    getline(in, residue->header);
    if (trim(residue->header) == "STOP")
        return false;
    getline(in, residue->file);

    getline(in, line);
    ss.str(line);
    ss >> residue->name;
    residue->coordinate_type = extract_coordinate_type(ss);
    residue->output_format = extract_output_format(ss);

    getline(in, line);
    ss.str(line);
    residue->geometry_type = extract_geometry_type(ss);
    residue->dummy_atom_omission = extract_dummy_omission(ss);
    ss >> residue->dummy_atom_type;
    residue->dummy_atom_position = extract_dummy_position(ss);

    getline(in, line);
    residue->cutoff = convert_string<double>(line);

    while (getline(in, line) && !trim(line).empty()) {
        PrepFileAtom *atom = new PrepFileAtom;
        ss.clear();
        ss.str(line);
        ss >> atom->index >> atom->name >> atom->type;
        atom->topological_type = extract_topological_type(ss);
        ss >> atom->bond_index >> atom->angle_index >> atom->dihedral_index >>
              atom->bond_length >> atom->angle >> atom->dihedral >>
              atom->charge;
        residue->atoms.push_back(atom);
    }

    bool done = false;
    while (!done) {
        while (getline(in, line) && trim(line).empty()) {}
        switch (get_other_section(line)) {
            case kSectionLoop:
                while (getline(in, line) && !trim(line).empty()) {
                    ss.clear();
                    ss.str(line);
                    string atom_names[2];
                    ss >> atom_names[0] >> atom_names[1];
                    int from = residue->find(atom_names[0]);
                    int to = residue->find(atom_names[1]);
                    residue->loops.push_back(PrepFileResidue::Loop(from, to));
                }
                break;
            case kSectionImproper:
                while (getline(in, line) && !trim(line).empty()) {
                    PrepFileResidue::ImproperDihedral dihedral;
                    ss.clear();
                    ss.str(line);
                    ss >> dihedral.atom_names[0] >> dihedral.atom_names[1] >>
                          dihedral.atom_names[2] >> dihedral.atom_names[3];
                    residue->improper_dihedrals.push_back(dihedral);
                }
                break;
            case kSectionDone:
                done = true;
                break;
            case kSectionOther:
                warning("unrecognized section in prep file");
                break;
        }
    }
    residues_[residue->name] = residue;
    return true;
}

PrepFileResidue::CoordinateType PrepFile::extract_coordinate_type(
        std::istream& in) const {
    string s;
    in >> s;
    if (s == "XYZ")
        return PrepFileResidue::kXYZ;
    else
        return PrepFileResidue::kINT;
}

PrepFileResidue::DummyAtomOmission PrepFile::extract_dummy_omission(
        std::istream& in) const {
    string s;
    in >> s;
    if (s == "NOMIT")
        return PrepFileResidue::kNomit;
    else
        return PrepFileResidue::kOmit;
}

PrepFileResidue::DummyAtomPosition PrepFile::extract_dummy_position(
        std::istream& in) const {
    string s;
    in >> s;
    if (s == "ALL")
        return PrepFileResidue::kPositionAll;
    else
        return PrepFileResidue::kPositionBeg;
}

PrepFileResidue::OutputFormat PrepFile::extract_output_format(
        std::istream& in) const {
    int val;
    in >> val;
    if (val == 1)
        return PrepFileResidue::kBinary;
    else
        return PrepFileResidue::kFormatted;
}

PrepFileResidue::GeometryType PrepFile::extract_geometry_type(
        std::istream& in) const {
    string s;
    in >> s;
    if (s == "CHANGE")
        return PrepFileResidue::kGeometryChange;
    else
        return PrepFileResidue::kGeometryCorrect;
}

PrepFileAtom::TopologicalType PrepFile::extract_topological_type(
        std::istream& in) const {
    string s;
    in >> s;
    if (s == "M")
        return PrepFileAtom::kTopTypeM;
    else if (s == "S")
        return PrepFileAtom::kTopTypeS;
    else if (s == "B")
        return PrepFileAtom::kTopTypeB;
    else if (s == "E")
        return PrepFileAtom::kTopTypeE;
    else
        return PrepFileAtom::kTopType3;
}

PrepFile::OtherSection PrepFile::get_other_section(const string& line) const {
    if (line == "LOOP")
        return kSectionLoop;
    else if (line == "IMPROPER")
        return kSectionImproper;
    else if (line == "DONE")
        return kSectionDone;
    return kSectionOther;
}

void PrepFileSet::load(const PrepFile& prep_file) {
    PrepFile::const_iterator it;
    for (it = prep_file.begin(); it != prep_file.end(); ++it) {
        if (residues_.find(it->first) == residues_.end()) {
            residues_.insert(*it);
        } else  {
            residues_[it->first] = it->second;
            warning("PrepFile - Overwriting prep file " + it->first +
                    " in prep file set.");
        }
    }
}

Residue *BuildPrepFileResidue::operator()(
        const PrepFileResidue& prep_file_residue) const {
    const vector<PrepFileAtom*> atom_list = prep_file_residue.atoms;
    // The index of the parent of each atom
    vector<int> parent_list(atom_list.size());
    set_parent_list(atom_list, parent_list);

    vector<Coordinate*> coordinates(atom_list.size());
    set_dummy_coordinates(atom_list, coordinates);

    Graph *bonds = new Graph(atom_list.size() - 3);
    for (size_t i = 3; i < atom_list.size(); i++) {
        int parent1 = parent_list[i];
        // Grandparent of i
        int parent2 = parent_list[parent1];
        // Great-grandparent of i
        int parent3 = parent_list[parent2];
        if (parent1 > 2)
            bonds->add_edge(i - 3, parent1 - 3);

        coordinates[i] = new Coordinate(
            calculate_point(*coordinates[parent3], *coordinates[parent2],
                            *coordinates[parent1],
                            to_radians(atom_list[i]->angle),
                            to_radians(atom_list[i]->dihedral),
                            atom_list[i]->bond_length)
        );
    }
    for (size_t i = 0; i < prep_file_residue.loops.size(); i++) {
        const PrepFileResidue::Loop& loop = prep_file_residue.loops[i];
        bonds->add_edge(loop.from - 3, loop.to - 3);
    }

    vector<boost::shared_ptr<Atom> > *atoms =
        new vector<boost::shared_ptr<Atom> >;
    atoms->reserve(atom_list.size() - 3);

    for (size_t i = 3; i < atom_list.size(); i++) {
        Element element = get_element_by_char(atom_list[i]->name[0]);
        boost::shared_ptr<Atom> atom_ptr(new Atom(element, *coordinates[i],
                                                  atom_list[i]->name,
                                                  atom_list[i]->type,
                                                  atom_list[i]->charge));
        atoms->push_back(atom_ptr);
    }

    std::for_each(coordinates.begin(), coordinates.end(), DeletePtr());

    return new Residue(bonds, prep_file_residue.name, atoms);
}

// Pretty sure there's a good way to simplify this
void BuildPrepFileResidue::set_parent_list(const vector<PrepFileAtom*>& atoms,
                                           vector<int>& parent_list) {
    stack<int> st;
    for (int i = atoms.size() - 1; i >= 0; i--) {
        size_t stack_size = st.size();
        switch (atoms[i]->topological_type) {
            case PrepFileAtom::kTopTypeM:
                while (!st.empty()) {
                    parent_list[st.top()] = i;
                    st.pop();
                }
                break;
            case PrepFileAtom::kTopTypeS:
                if (stack_size > 0) {
                    parent_list[st.top()] = i;
                    st.pop();
                }
                break;
            case PrepFileAtom::kTopTypeB:
                for (size_t j = 0; j < std::min((size_t)2, stack_size); j++) {
                    parent_list[st.top()] = i;
                    st.pop();
                }
                break;
            case PrepFileAtom::kTopTypeE:
                break;
            case PrepFileAtom::kTopType3:
                for (size_t j = 0; j < std::min((size_t)3, stack_size); j++) {
                    parent_list[st.top()] = i;
                    st.pop();
                }
                break;
            default:
                warning("Unrecognized topological type ");
                break;
        }
        st.push(i);
    }
}

void BuildPrepFileResidue::set_dummy_coordinates(
        const vector<PrepFileAtom*>& atoms, vector<Coordinate*>& coordinates) {
    coordinates[0] = new Coordinate(0.0, 0.0, 0.0);

    double dist1 = atoms[1]->bond_length;
    coordinates[1] = new Coordinate(dist1, 0.0, 0.0);

    double dist2 = atoms[2]->bond_length;
    double angle = atoms[2]->angle;
    coordinates[2] = new Coordinate(dist1 - cos(to_radians(angle))*dist2,
                                    sin(to_radians(angle))*dist2,
                                    0.0);
}

}  // namespace gmml
