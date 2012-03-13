// Author: Robert Davis

#include "gmml/internal/library_file.h"

#include <cassert>
#include <cctype>

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "gmml/internal/atom.h"
#include "gmml/internal/environment.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/stubs/logging.h"
#include "utilities.h"

using std::istringstream;
using std::map;
using std::string;
using std::vector;

namespace gmml {
namespace {

void remove_quotes(string& str) {
    str.erase(std::remove(str.begin(), str.end(), '\"'), str.end());
}

void remove_spaces(string& str) {
    str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
}

}  // namespace

void LibraryFile::read(std::istream& in) {
    string line;
    std::getline(in, line);
    vector<string> names;
    while (in.peek() != '!') {
        std::getline(in, line);
        remove_quotes(line);
        remove_spaces(line);
        names.push_back(line);
    }
    std::getline(in, line);
    for (int i = 0; i < names.size(); i++) {
        boost::shared_ptr<LibraryFileStructure> structure(
            new LibraryFileStructure(in));

        iterator it = structures_.lower_bound(names[i]);
        if (it == structures_.end() ||
                structures_.key_comp()(names[i], it->first)) {
            structures_.insert(it, std::make_pair(names[i], structure));
        } else {
            it->second = structure;
            LOG(WARNING) << "LibraryFile - Overwriting structure " <<
                             names[i] << ".";
        }
    }
}

void LibraryFileSet::load(const LibraryFile& library_file) {
    LibraryFile::const_iterator it;
    for (it = library_file.begin(); it != library_file.end(); ++it) {
        LibraryFile::iterator lb = structures_.lower_bound(it->first);
        if (lb == structures_.end() ||
                structures_.key_comp()(it->first, lb->first)) {
            structures_.insert(lb, *it);
        } else {
            lb->second = it->second;
            LOG(WARNING) << "LibraryFileSet - Overwriting structure " <<
                             it->first << ".";
        }
    }
}

void LibraryFileStructure::clone_from(const LibraryFileStructure& structure) {
    BoxedStructure::clone_from(structure);
}

namespace {

struct RawAtoms {
    RawAtoms() : atoms(new vector<Atom*>) {}

    ~RawAtoms() { delete atoms; }

    void add_atom(Atom *atom, int index) {
        atoms->push_back(atom);
        indices.push_back(index);
    }

    vector<Atom*> *atoms;
    vector<int> indices;
};

}  // namespace

void LibraryFileStructure::read(std::istream& in) {
    // This is a map from the residue serial numbers found in the atom section
    // and referenced in other sections to the atoms in the residues.
    map<int, RawAtoms*> residues;
    int cur_atom = 0;
    while (in.peek() != '!') {
        std::pair<Atom*, int> atom = read_atom(in);
        // We don't use the serial field because the sections that follow seem
        // to only use the (1-based) index of the atom in the atom list. The
        // serial is the same number, however, for all the library files I've
        // seen.
        typedef std::pair<map<int, RawAtoms*>::iterator, bool> InsertRet;
        InsertRet ret = residues.insert(
                std::make_pair(atom.second, static_cast<RawAtoms*>(NULL)));
        if (ret.second)
            ret.first->second = new RawAtoms;
        ret.first->second->add_atom(atom.first, cur_atom++);
    }

    // This is a map from the (0-based) index of the atom in the atom list to
    // the index of the atom in the structure.
    map<int, int> atom_map;

    // This is a map from the library file residue numbers to the index of
    // the residues in the structure. We construct it here so we know how to
    // correspond the residue names in a later section to the residues in the
    // structure.
    map<int, int> residue_map;

    int cur_residue = 0;
    std::map<int, RawAtoms*>::const_iterator it;
    for (it = residues.begin(); it != residues.end(); ++it) {
        residue_map[it->first] = cur_residue++;
        IndexedResidue *indexed_residue =
                new IndexedResidue("", it->second->atoms, NULL);
        for (int i = 0; i < it->second->indices.size(); i++) {
            indexed_residue->set_atom_index(i, it->second->indices[i]);
        }
        // We'll go ahead and insert the residues without names. The names
        // come from a later section.
        append(atom_map, indexed_residue);
    }

    // Keep reading the stream until we bump into another residue
    string line;
    while (getline(in, line) && line.find("unit.atoms ") == string::npos) {
        if (line.find("boundbox") != string::npos) {
            read_box(in);
        } else if (line.find("connect array") != string::npos) {
            read_connect_atoms(in);
        } else if (line.find("connectivity") != string::npos) {
            read_connectivity_info(in, atom_map);
        } else if (line.find("unit.name ") != string::npos) {
            // The name here seems to be useless. We use the name at the top of
            // the file instead.
        } else if (line.find("positions") != string::npos) {
            read_positions(in, atom_map);
        } else if (line.find("unit.residues ") != string::npos) {
            read_residue_info(in, residue_map);
        }
    }
}

std::pair<Atom*, int> LibraryFileStructure::read_atom(
        std::istream& in) const {
    string line;
    getline(in, line);
    istringstream ss(line);
    string name;
    string type;
    int residue;
    int serial;
    int atomic_number;
    double charge;
    int int_t;
    ss >> name >> type >> int_t >> residue >> int_t >> serial >>
	  atomic_number >> charge;
    remove_quotes(name);
    remove_quotes(type);
    // The coordinates will be set in the sections that follow.
    Atom *atom = new Atom(static_cast<Element>(atomic_number),
                          Coordinate(), name, type, charge);
    return std::make_pair(atom, residue);
}

void LibraryFileStructure::read_box(std::istream& in) {
    string line;
    getline(in, line);
    istringstream ss(line);
    double is_box_set;
    ss >> is_box_set;
    if (is_box_set < 0) {
        for (int i = 0; i < 4; i++)
            getline(in, line);
        return;
    }

    box_ = new Box(kNotSet, kNotSet, kNotSet, kNotSet);

    getline(in, line);
    ss.str(line);
    ss.clear();
    ss >> box_->angle;

    getline(in, line);
    ss.str(line);
    ss.clear();
    ss >> box_->length;

    getline(in, line);
    ss.str(line);
    ss.clear();
    ss >> box_->width;

    getline(in, line);
    ss.str(line);
    ss.clear();
    ss >> box_->height;
}

void LibraryFileStructure::read_connect_atoms(std::istream& in) {
    string line;
    int head = -1;
    int tail = -1;
    if (in.peek() != '!') {
        getline(in, line);
        // The file numbers are 1-based.
        head = convert_string<int>(line) - 1;
    }
    if (in.peek() != '!') {
        getline(in, line);
        tail = convert_string<int>(line) - 1;
    }

    set_head(head);
    set_tail(tail);

    while (in.peek() != '!') {
        getline(in, line);
    }
}

void LibraryFileStructure::read_connectivity_info(
        std::istream& in, const map<int, int>& atom_map) {
    string line;
    while (in.peek() != '!') {
        getline(in, line);
        istringstream ss(line);
        int from, to;
        ss >> from >> to;
        // The indices are 1-based in the file so we decrement.
        from--;
        to--;
        map<int, int>::const_iterator from_it, to_it;
        from_it = atom_map.find(from);
        to_it = atom_map.find(to);
        assert(from_it != atom_map.end());
        assert(from_it->second < size());
        assert(to_it != atom_map.end());
        assert(to_it->second < size());
        this->add_bond(from_it->second, to_it->second);
    }
}

void LibraryFileStructure::read_positions(std::istream& in,
                                          const map<int, int>& atom_map) {
    string line;
    for (int i = 0; i < size(); i++) {
        getline(in, line);
        istringstream ss(line);
        double x, y, z;
        ss >> x >> y >> z;
        map<int, int>::const_iterator it = atom_map.find(i);
        assert(it != atom_map.end());
        assert(it->second < size());
        atoms(it->second)->set_coordinate(x, y, z);
    }
}

void LibraryFileStructure::read_residue_info(std::istream& in,
                                             const map<int, int>& residue_map) {
    string line;
    for (int i = 0; i < residue_count(); i++) {
        getline(in, line);
        istringstream ss(line);
        string residue_name;
        int residue_number;
        ss >> residue_name >> residue_number;
        remove_quotes(residue_name);
        remove_spaces(residue_name);
        map<int, int>::const_iterator it = residue_map.find(residue_number);
        assert(it != residue_map.end());
        assert(it->second < residue_count());
        residues(it->second)->set_name(residue_name);
    }
}

}  // namespace gmml

