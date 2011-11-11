#include "gmml/internal/structure.h"

#include <algorithm>
#include <deque>
#include <string>
#include <vector>

#include "gmml/internal/amber_top_builder.h"
#include "gmml/internal/amber_top_file.h"
#include "gmml/internal/coordinate_file.h"
#include "gmml/internal/environment.h"
#include "gmml/internal/geometry.h"
#include "gmml/internal/parameter_file.h"
#include "gmml/internal/pdb_file.h"
#include "gmml/internal/pdb_file_builder.h"
#include "gmml/internal/prep_file.h"
#include "gmml/internal/sander_minimize.h"
#include "gmml/internal/utilities.h"

using std::deque;
using std::list;
using std::map;
using std::string;
using std::vector;

namespace gmml {

Structure *Structure::build_from_pdb(const PdbFile& pdb_file) {
    // This is a map from the chain id to a map between the residue sequence
    // id and a list of atoms indexed by a pdb atom sequence id. Either or
    // both of the maps should probably be a hash table.
    map<char, map<int, IndexedResidue> > residue_map;
    const list<boost::shared_ptr<PdbAtomCard> >& atom_cards =
        pdb_file.atom_cards();

    typedef PdbFile::AtomCardPtr AtomCardPtr;
    for (list<AtomCardPtr>::const_iterator it = atom_cards.begin();
            it != atom_cards.end(); ++it) {
        AtomPtr new_atom(new Atom(get_element_by_char((*it)->element[0]),
                                  Coordinate((*it)->x, (*it)->y, (*it)->z),
                                  (*it)->name, (*it)->charge));
        residue_map[(*it)->chain_id][(*it)->res_seq].atoms.push_back(
                new IndexedAtom(new_atom, (*it)->serial));
        if (residue_map[(*it)->chain_id][(*it)->res_seq].name == "")
            residue_map[(*it)->chain_id][(*it)->res_seq].name = (*it)->res_name;
    }

    Structure *structure = new Structure;
    map<int, int> atom_map;
    map<char, map<int, IndexedResidue> >::const_iterator chain_it;
    map<int, IndexedResidue>::const_iterator residue_it;
    for (chain_it = residue_map.begin(); chain_it != residue_map.end();
            ++chain_it) {
        for (residue_it = chain_it->second.begin();
                residue_it != chain_it->second.end(); ++residue_it) {
            structure->add_indexed_residue(atom_map, residue_it->second);
            std::for_each(residue_it->second.atoms.begin(),
                          residue_it->second.atoms.end(), DeletePtr());
        }
    }

    typedef PdbFile::ConnectCardPtr ConnectCardPtr;
    const list<ConnectCardPtr>& connect_cards = pdb_file.connect_cards();
    for (list<ConnectCardPtr>::const_iterator it = connect_cards.begin();
            it != connect_cards.end(); ++it) {
        int source = atom_map[(*it)->connect1];
        if ((*it)->connect2 != kNotSet && atom_map[(*it)->connect2] < source)
             structure->add_bond(source, atom_map[(*it)->connect2]);
        if ((*it)->connect3 != kNotSet && atom_map[(*it)->connect3] < source)
             structure->add_bond(source, atom_map[(*it)->connect3]);
        if ((*it)->connect4 != kNotSet && atom_map[(*it)->connect4] < source)
             structure->add_bond(source, atom_map[(*it)->connect4]);
        if ((*it)->connect5 != kNotSet && atom_map[(*it)->connect5] < source)
             structure->add_bond(source, atom_map[(*it)->connect5]);
    }

    return structure;
}

void Structure::add_indexed_residue(map<int, int>& atom_map,
                                    const IndexedResidue& residue) {
    int cur_size = atoms_.size();
    residues_->push_back(new detail::StructureResidue(
            residue.name, cur_size, residue.atoms.size()));
    for (int i = 0; i < residue.atoms.size(); i++) {
        AtomPtr atom = residue.atoms[i]->atom;
        atoms_.push_back(atom);
        atom_map[residue.atoms[i]->index] = cur_size + i;
    }
    bonds_->append(Graph(residue.atoms.size()));
}

Structure *Structure::clone() const {
/*
    using detail::StructureResidue;

    Structure *structure = new Structure;
    structure->atoms_.reserve(atoms_.size());
    for (size_t i = 0; i < atoms_.size(); i++)
	structure->atoms_.push_back(AtomPtr(atoms_[i]->clone()));

    structure->residues_->reserve(residues_->size());
    for (size_t i = 0; i < residues_->size(); i++) {
	structure->residues_->push_back(
            new StructureResidue(*residues_->at(i)));
    }

    structure->bonds_ = bonds_->clone();

    return structure;
*/
    Structure *structure = new Structure;
    structure->clone_from(*this);
    return structure;
}

void Structure::clone_from(const Structure& structure) {
    using detail::StructureResidue;

    atoms_.reserve(structure.atoms_.size());
    for (size_t i = 0; i < structure.atoms_.size(); i++)
        atoms_.push_back(AtomPtr(structure.atoms_[i]->clone()));

    residues_->reserve(structure.residues_->size());
    for (size_t i = 0; i < structure.residues_->size(); i++) {
        residues_->push_back(new StructureResidue(*(*structure.residues_)[i]));
    }

    bonds_ = structure.bonds_->clone();
}

void Structure::append(const Structure& rhs) {
    size_t old_size = size();
    atoms_.reserve(old_size + rhs.atoms_.size());
    atoms_.insert(end(), rhs.begin(), rhs.end());

    bonds_->append(*rhs.bonds_);

    using detail::StructureResidue;
    residues_->reserve(residues_->size() + rhs.residues_->size());
    for (size_t i = 0; i < rhs.residues_->size(); i++) {
	residues_->push_back(
	    new StructureResidue(rhs.residues_->at(i)->name,
                                 rhs.residues_->at(i)->start_index + old_size,
                                 rhs.residues_->at(i)->size)
        );
    }
}

int Structure::append(const Residue *residue) {
    size_t old_size = size();
    atoms_.reserve(old_size + residue->size());
    atoms_.insert(end(), residue->begin(), residue->end());

    bonds_->append(*residue->bonds());

    using detail::StructureResidue;
    residues_->push_back(new StructureResidue(residue->name(), old_size,
                                              residue->size()));
    return residues_->size() - 1;
}

int Structure::append(const string& prep_code) {
    Residue *residue = build_prep_file(prep_code);
    return append(residue);
}

void Structure::shift(double x, double y, double z) {
    for (iterator it = begin(); it != end(); ++it)
        (*it)->translate(x, y, z);
}

int Structure::attach(Residue *new_residue, const string& new_atom_name,
                      int residue_index, const string& target_atom_name) {
    StructureAttach s;
    return s(*this, new_residue, new_atom_name, residue_index,
             target_atom_name);
}

int Structure::attach(const string& prep_code, const string& new_atom_name,
                      int residue_index, const string& target_atom_name) {
    Residue *residue = build_prep_file(prep_code);
    return attach(residue, new_atom_name, residue_index, target_atom_name);
}

MinimizationResults *Structure::minimize(const string& input_file) {
    SanderMinimize minimize;
    return minimize(*this, input_file);
}

void Structure::translate_residue(int residue_index, double x, double y,
                                  double z) {
    int residue_start = residues_->at(residue_index)->start_index;
    int residue_size = residues_->at(residue_index)->size;
    for (int i = residue_start; i < residue_size + residue_start; i++)
        atoms_[i]->translate(x, y, z);
}

void Structure::remove_residues(const std::vector<int>& residues) {
    vector<size_t> atoms;
    vector<int> new_indices(residues_->size());
    for (int i = 0; i < residues.size(); i++)
        new_indices[residues[i]] = -1;

    deque<bool> is_deleted(atoms_.size(), false);
    int cur_shift = 0;
    int i = 0;
    int k = 0;
    while (i < residues_->size()) {
        if (new_indices[k] == -1) {
            int start_index = (*residues_)[i]->start_index;
            int residue_size = (*residues_)[i]->size;
            cur_shift += residue_size;
            for (int j = start_index; j < start_index + residue_size; j++) {
                atoms.push_back(j);
                is_deleted[j] = true;
            }
            delete (*residues_)[i];
            residues_->erase(residues_->begin() + i);
        } else {
            (*residues_)[i]->start_index -= cur_shift;
            i++;
        }
        k++;
    }

    i = 0;
    k = 0;
    while (i < atoms_.size()) {
        if (is_deleted[k]) {
            atoms_.erase(atoms_.begin() + i);
        } else {
            i++;
        }
        k++;
    }

    bonds_->remove_vertex_list(atoms);
}

Graph *Structure::get_link_graph() const {
    Graph *graph = new Graph(residues_->size());
    vector<size_t> *index_table = get_residue_index_table();
    for (int i = 0; i < atoms_.size(); i++) {
        int residue_index = (*index_table)[i];
        const AdjList& edges = bonds(i);
        int j = 0;
        // The second conditions ensures that we don't add the same edge twice.
        while (j < edges.size() && edges[j] < i) {
            graph->add_edge(residue_index, (*index_table)[edges[j]]);
            j++;
        }
    }
    delete index_table;
    return graph;
}

int Structure::get_atom_index(int residue_index, int atom_index) const {
    int residue_start = residues_->at(residue_index)->start_index;
    return residue_start + atom_index;
}

int Structure::get_atom_index(int residue_index,
                              const string& atom_name) const {
    int residue_start = residues_->at(residue_index)->start_index;
    int residue_size = residues_->at(residue_index)->size;
    for (int i = residue_start; i < residue_size + residue_start; i++)
        if (atoms_[i]->name() == atom_name)
            return i;
    
    return -1;
}

int Structure::get_anomeric_index(int residue_index) const {
    int residue_start = get_residue_start(residue_index);
    int residue_size = get_residue_size(residue_index);
    int residue_end = residue_start + residue_size - 1;
    for (int i = residue_start; i <= residue_end; i++) {
        if (atoms_[i]->element() == kElementO)
            continue;
        const AdjList& adj_list = bonds(i);
        for (int j = 0; j < adj_list.size(); j++) {
            if ((adj_list[j] < residue_start || adj_list[j] > residue_end) &&
                    atoms_[adj_list[j]]->element() == kElementO)
                return i;
        }
    }
    return -1;
}

int Structure::get_parent_atom(int residue_index) const {
    int anomeric_index = get_anomeric_index(residue_index);
    if (anomeric_index != -1) {
        int residue_start = get_residue_start(residue_index);
        int residue_size = get_residue_size(residue_index);
        int residue_end = residue_start + residue_size - 1;
        const AdjList& adj_list = bonds(anomeric_index);
        for (int i = 0; i < adj_list.size(); i++) {
            if ((adj_list[i] < residue_start || adj_list[i] > residue_end) &&
                    atoms_[adj_list[i]]->element() == kElementO)
                return adj_list[i];
        }
    }
    return -1;
}

int Structure::get_parent_residue(int residue_index) const {
    int parent_atom = get_parent_atom(residue_index);
    if (parent_atom != -1)
        return get_residue_index(parent_atom);
    return -1;
}

vector<int> *Structure::get_flexible_linkages() const {
    vector<int> *residues = new vector<int>;
    for (int i = 0; i < residues_->size(); i++) {
        int parent_atom = get_parent_atom(i);
        int parent_index = get_residue_index(parent_atom);
        if (parent_atom == -1)
            continue;
        const AdjList& adj_list = bonds(parent_atom);
        bool is_flexible = false;
        for (int j = 0; j < adj_list.size(); j++) {
            if (get_residue_index(adj_list[j]) == parent_index &&
                    !is_cyclic(adj_list[j]) &&
                    atoms_[adj_list[j]]->element() == kElementC) {
                is_flexible = true;
                break;
            }
        }
        if (is_flexible)
            residues->push_back(i);
    }
    return residues;
}

bool Structure::is_cyclic(int atom_index) const {
    const AdjList& adj_list = bonds(atom_index);
    bool in_cycle = false;
    for (int i = 0; i < adj_list.size(); i++) {
        vector<size_t> *found_atoms = bonds_->edge_bfs(atom_index, adj_list[i]);
        if (std::find(found_atoms->begin(), found_atoms->end(), atom_index) !=
                found_atoms->end()) {
            in_cycle = true;
            delete found_atoms;
            break;
        }
        delete found_atoms;
    }
    return in_cycle;
}

void Structure::set_dihedral(size_t atom1, size_t atom2, size_t atom3,
                             size_t atom4, double degrees) {
    vector<size_t> *atoms = bonds_->edge_bfs(atom2, atom3);

    double current_dihedral = 
        measure(atoms_[atom1]->coordinate(), 
                atoms_[atom2]->coordinate(),
                atoms_[atom3]->coordinate(), 
                atoms_[atom4]->coordinate());

    RotationMatrix matrix(atoms_[atom2]->coordinate(),
                          Vector<3>(atoms_[atom3]->coordinate(),
                                    atoms_[atom2]->coordinate()),
                          current_dihedral - to_radians(degrees));

    for (vector<size_t>::iterator it = atoms->begin();
            it != atoms->end(); ++it)
        matrix.apply(atoms_[*it]->coordinate());

    delete atoms;
}

void Structure::set_dihedral(int residue1_index, const string& atom1_name,
                             int residue2_index, const string& atom2_name,
                             int residue3_index, const string& atom3_name,
                             int residue4_index, const string& atom4_name,
                             double degrees) {
    int atom1_index = get_atom_index(residue1_index, atom1_name);
    int atom2_index = get_atom_index(residue2_index, atom2_name);
    int atom3_index = get_atom_index(residue3_index, atom3_name);
    int atom4_index = get_atom_index(residue4_index, atom4_name);

    if (atom1_index == -1 || atom2_index == -1 || atom3_index == -1 ||
            atom4_index == -1) {
        return;
    }

    set_dihedral(atom1_index, atom2_index, atom3_index, atom4_index, degrees);
}

bool Structure::set_phi(int residue_index, double degrees) {
    int carbon_index = get_anomeric_index(residue_index);
    if (carbon_index == -1) {
        return false;
    }
    int atom2_index;
    if (atoms_[carbon_index]->name() == "C1") {
        atom2_index = get_atom_index(residue_index, "H1");
    }
    else if (atoms_[carbon_index]->name()[0] == 'C') {
        int carbon_number = atoms_[carbon_index]->name()[1] - '0';
        string atom2_name = "C" + to_string(carbon_number - 1);
        atom2_index = get_atom_index(residue_index, atom2_name);
    }
    //residue must start with Cx
    else {
        return false;
    }
    int oxygen_index = kNotSet;
    {
        const AdjList& adj_list = bonds(carbon_index);
        for (int i = 0; i < adj_list.size(); i++) {
            int atom_index = adj_list[i];
            if (get_residue_index(atom_index) != residue_index &&
                    atoms_[atom_index]->name()[0] == 'O') {
                oxygen_index = atom_index;
                break;
            }
        }
        if (oxygen_index == kNotSet ||
                 atoms_[oxygen_index]->name().size() < 2 ||
                 !is_number(atoms_[oxygen_index]->name()[1])) {
            return false;
        }
    }
    int oxygen_number = char_to_number(atoms_[oxygen_index]->name()[1]);
    int attaching_carbon_index = kNotSet;
    {
        const AdjList& adj_list = bonds(oxygen_index);
        for (int i = 0; i < adj_list.size(); i++) {
            int atom_index = adj_list[i];
            if (get_residue_index(atom_index) == 
                        get_residue_index(oxygen_index) &&
                    atoms_[atom_index]->name() ==
                        "C" + to_string(oxygen_number))
            attaching_carbon_index = atom_index;
        }
        if (attaching_carbon_index == kNotSet) {
            return false;
        }
    }
    set_dihedral(attaching_carbon_index, oxygen_index, carbon_index,
                 atom2_index, degrees);
    return true;
}

bool Structure::set_psi(int residue_index, double degrees) {
    int anomeric_carbon_index = get_anomeric_index(residue_index);
    if (anomeric_carbon_index == -1) {
        return false;
    }
    int oxygen_index = kNotSet;
    {
        const AdjList& adj_list = bonds(anomeric_carbon_index);
        for (int i = 0; i < adj_list.size(); i++) {
            int atom_index = adj_list[i];
            if (get_residue_index(atom_index) != residue_index &&
                    atoms_[atom_index]->name()[0] == 'O') {
                oxygen_index = atom_index;
                break;
            }
        }
    }
    if (oxygen_index == kNotSet || atoms_[oxygen_index]->name().size() < 2 ||
            !is_number(atoms_[oxygen_index]->name()[1])) {
        return false;
    }
    int oxygen_number = char_to_number(atoms_[oxygen_index]->name()[1]);
    int carbon_index = kNotSet;
    {
        const AdjList& adj_list = bonds(oxygen_index);
        for (int i = 0; i < adj_list.size(); i++) {
            int atom_index = adj_list[i];
            if (get_residue_index(atom_index) ==
                        get_residue_index(oxygen_index) &&
                    atoms_[atom_index]->name() ==
                        "C" + to_string(oxygen_number))
                carbon_index = atom_index;
        }
        if (carbon_index == kNotSet) {
            return false;
        }
    }
    string fourth_atom_name;
    if (is_cyclic(carbon_index))
        fourth_atom_name = "H" + to_string(oxygen_number);
    else
        fourth_atom_name = "C" + to_string(oxygen_number - 1);
    int fourth_atom_index = kNotSet;
    {
        const AdjList& adj_list = bonds(carbon_index);
        for (int i = 0; i < adj_list.size(); i++) {
            int atom_index = adj_list[i];
            if (atoms_[atom_index]->name() == fourth_atom_name) {
                fourth_atom_index = atom_index;
                break;
            }
        }
        if (fourth_atom_index == kNotSet) {
            return false;
        }
    }
    set_dihedral(fourth_atom_index, carbon_index, oxygen_index,
                 anomeric_carbon_index, degrees);
    return true;
}

bool Structure::set_omega(int residue_index, double degrees) {
    int anomeric_carbon_index = get_anomeric_index(residue_index);
    if (anomeric_carbon_index == -1) {
        return false;
    }
    int oxygen_index = kNotSet;
    int adjacent_residue_index = kNotSet;
    {
        const AdjList& adj_list = bonds(anomeric_carbon_index);
        for (int i = 0; i < adj_list.size(); i++) {
            int atom_index = adj_list[i];
            if (get_residue_index(atom_index) != residue_index &&
                    atoms_[atom_index]->name()[0] == 'O') {
                oxygen_index = atom_index;
                adjacent_residue_index = get_residue_index(atom_index);
                break;
            }
        }
        if (oxygen_index == kNotSet ||
                atoms_[oxygen_index]->name().size() < 2 ||
                !is_number(atoms_[oxygen_index]->name()[1])) {
            return false;
        }
    }
    int oxygen_number = char_to_number(atoms_[oxygen_index]->name()[1]);
    int carbon1_index = get_atom_index(adjacent_residue_index,
                                       "C" + to_string(oxygen_number));
    if (is_cyclic(carbon1_index))
        return false;
    int carbon2_index = get_atom_index(adjacent_residue_index,
                                       "C" + to_string(oxygen_number - 1));
    int other_oxygen_index = get_atom_index(adjacent_residue_index,
                                            "O" + to_string(oxygen_number - 1));
    if (carbon1_index == -1 || carbon2_index == -1 || other_oxygen_index == -1)
        return false;
    set_dihedral(other_oxygen_index, carbon2_index, carbon1_index, oxygen_index,
                 degrees);
    return true;
}

Graph *Structure::get_residue_link_table() const {
    vector<size_t> *residue_index_table = get_residue_index_table();

    Graph *links = new Graph(residues_->size());
    for (size_t i = 0; i < atoms_.size(); i++) {
        for (size_t j = 0; j < bonds_->edges(i).size(); j++)
            if (residue_index_table->at(bonds_->edges(i)[j]) !=
                    residue_index_table->at(i))
                links->add_edge(residue_index_table->at(i),
                                residue_index_table->at(bonds_->edges(i)[j]));
    }

    delete residue_index_table;

    return links;
}

void Structure::set_residue_angle(int atom1, int atom2, int atom3,
                                  int residue_index, double radians) {
    double cur_angle = measure(atoms_[atom1]->coordinate(),
                               atoms_[atom2]->coordinate(),
                               atoms_[atom3]->coordinate());

    RotationMatrix matrix(
	atoms_[atom2]->coordinate(),
	Vector<3>(atoms_[atom1]->coordinate(),
                  atoms_[atom2]->coordinate()).cross(
	              Vector<3>(atoms_[atom3]->coordinate(),
                                atoms_[atom2]->coordinate())),
        radians - cur_angle
    );

    int residue_start = residues_->at(residue_index)->start_index;
    int residue_size = residues_->at(residue_index)->size;
    for (int i = residue_start; i < residue_size + residue_start; i++)
        matrix.apply(atoms_[i]->coordinate());
}

PdbFile *Structure::build_pdb_file() const {
    PdbFileBuilder builder;
    return builder.build(*this);
}

void Structure::print_pdb_file(const string& file_name) const {
    PdbFile *file = build_pdb_file();
    file->print(file_name);
    delete file;
}

void Structure::print_pdb_file() const {
    PdbFile *file = build_pdb_file();
    file->print();
    delete file;
}

AmberTopFile *Structure::build_amber_top_file(
        const ParameterFileSet& parm_set) const {
    AmberTopBuilder builder(parm_set);
    return builder.build(*this);
}

AmberTopFile *Structure::build_amber_top_file() const {
    AmberTopBuilder builder;
    return builder.build(*this);
}

void Structure::print_amber_top_file(const string& file_name,
                                     const ParameterFileSet& parm_set) const {
    AmberTopFile *top_file = build_amber_top_file(parm_set);
    top_file->print(file_name);
    delete top_file;
}

void Structure::print_amber_top_file(const string& file_name) const {
    AmberTopFile *top_file = build_amber_top_file();
    top_file->print(file_name);
    delete top_file;
}

void Structure::print_amber_top_file() const {
    AmberTopFile *top_file = build_amber_top_file();
    top_file->print();
    delete top_file;
}

CoordinateFile *Structure::build_coordinate_file() const {
    CoordinateFile *coordinate_file = new CoordinateFile;
    for (int i = 0; i < atoms_.size(); i++)
        coordinate_file->add_coordinate(atoms_[i]->coordinate());
    return coordinate_file;
}

void Structure::print_coordinate_file(const string& file_name) const {
    CoordinateFile *coordinate_file = build_coordinate_file();
    coordinate_file->print(file_name);
    delete coordinate_file;
}

void Structure::print_coordinate_file() const {
    CoordinateFile *coordinate_file = build_coordinate_file();
    coordinate_file->print();
    delete coordinate_file;
}

void Structure::load_coordinates(const CoordinateFile& coordinate_file) {
    if (coordinate_file.size() < atoms_.size()) {
        
        warning(string("Structure: Insufficient number of coordinates") +
                " in coordinate file");
        return;
    }
    for (int i = 0; i < coordinate_file.size(); i++)
        atoms_[i]->set_coordinate(coordinate_file[i]);
}

vector<size_t> *Structure::get_residue_index_table() const {
    vector<size_t> *table = new vector<size_t>(atoms_.size());
    size_t current_index = 0;
    for (size_t i = 0; i < residues_->size(); i++)
        for (size_t j = 0; j < residues_->at(i)->size; j++)
            (*table)[current_index++] = i;
    return table;
}

int StructureAttach::operator()(Structure& structure, Residue *new_residue,
                                const string& new_atom_name, int residue_index,
                                const string& target_atom_name) const {
    const Structure::AtomList& atoms = structure.atoms();
    int new_residue_index = structure.append(new_residue);
    int new_atom_index = structure.get_atom_index(new_residue_index,
                                                  new_atom_name);
    int target_atom_index = structure.get_atom_index(residue_index,
                                                     target_atom_name);
    structure.add_bond(new_atom_index, target_atom_index);

    const ParameterFileSet *parm_set = kDefaultEnvironment.parm_set();

    const ParameterFileBond *parameter_bond =
        parm_set->lookup(atoms[target_atom_index]->type(),
                         atoms[new_atom_index]->type());
    double bond_length = parameter_bond->length;

    Vector<3> direction = get_connection_direction(structure,
                                                   new_atom_index,
                                                   target_atom_index);
    direction.normalize();
    direction *= bond_length;

    //Vector<3> oxygen_position(Vector<3>(atoms_[new_atom_index]->coordinate()) +
    //                          direction);
    direction += Vector<3>(atoms[new_atom_index]->coordinate());
    Vector<3> oxygen_position = direction;
    VectorBase<3> offset(Vector<3>(atoms[target_atom_index]->coordinate()) -
                         oxygen_position);
    structure.translate_residue(new_residue_index, offset[0], offset[1],
                                offset[2]);

    // Set the bond angles
    const Structure::AdjList& adj_atoms = structure.bonds(target_atom_index);
    for (int i = 0; i < adj_atoms.size(); i++) {
        if (adj_atoms[i] != new_atom_index){
            int third_atom = adj_atoms[i];
            const ParameterFileAngle *parameter_angle = parm_set->lookup(
                atoms[third_atom]->type(),
                atoms[target_atom_index]->type(),
                atoms[new_atom_index]->type());

            if (parameter_angle != NULL)
                structure.set_residue_angle(third_atom, target_atom_index,
                                            new_atom_index, new_residue_index,
                                            to_radians(parameter_angle->angle));
         }
    }

    int carbon_number = kNotSet;
    if (new_atom_name.size() > 1 && is_number(new_atom_name[1]))
        carbon_number = char_to_number(new_atom_name[1]);

    int oxygen_number = kNotSet;
    if (target_atom_name.size() > 1 && is_number(target_atom_name[1]))
        oxygen_number = char_to_number(target_atom_name[1]);
    else
        oxygen_number = 1;

    if (is_set(carbon_number) && is_set(oxygen_number))
        set_dihedrals(structure, new_residue_index, residue_index,
                      carbon_number, oxygen_number);

    return new_residue_index;
}

Vector<3> StructureAttach::get_connection_direction(const Structure& structure,
                                                    int source_index,
                                                    int target_index) const {
    Vector<3> direction(0.0, 0.0, 0.0);
    const Structure::AtomList& atoms = structure.atoms();
    const Structure::AdjList& adj_atoms  = structure.bonds(source_index);
    for (int i = 0; i < adj_atoms.size(); i++) {
        if (adj_atoms[i] != target_index) {
            Vector<3> v(atoms[adj_atoms[i]]->coordinate(),
                        atoms[source_index]->coordinate());
            v.normalize();
            direction += v;
            direction.normalize();
         }
    }
    return direction;
}

void StructureAttach::set_dihedrals(Structure& structure, int new_residue_index,
                                    int target_residue_index, int carbon_number,
                                    int oxygen_number) const {
    string oxygen = to_string(oxygen_number);
    string carbon = to_string(carbon_number);

    if (oxygen_number == 5 || oxygen_number == 6) {
        structure.set_dihedral(target_residue_index,
                               "C" + to_string(oxygen_number - 1),
                               target_residue_index, "C" + oxygen,
                               target_residue_index, "O" + oxygen,
                               new_residue_index, "C" + carbon,
                               to_degrees(kPi));
    }
    else {
        structure.set_dihedral(target_residue_index, "H" + oxygen_number,
                               target_residue_index, "C" + oxygen,
                               target_residue_index, "O" + oxygen_number,
                               new_residue_index, "C" + carbon_number,
                               0.0);
    }

    structure.set_dihedral(target_residue_index, "C" + oxygen,
                           target_residue_index, "O" + oxygen,
                           new_residue_index, "C" + carbon,
                           new_residue_index,
                          "C" + to_string(carbon_number + 1),
                           to_degrees(kPi));

    if (oxygen_number == 5 || oxygen_number == 6) {
        structure.set_dihedral(target_residue_index,
                              "O" + to_string(oxygen_number - 1),
                              target_residue_index,
                              "C" + to_string(oxygen_number - 1),
                              target_residue_index, "C" + oxygen,
                              target_residue_index, "O" + oxygen,
                              to_degrees(kPi/3.0));
    }
}

}
