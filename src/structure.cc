#include "gmml/internal/structure.h"

#include <algorithm>
#include <deque>
#include <string>
#include <vector>

#include "gmml/internal/amber_top_builder.h"
#include "gmml/internal/amber_top_file.h"
#include "gmml/internal/atom.h"
#include "gmml/internal/coordinate_file.h"
#include "gmml/internal/environment.h"
#include "gmml/internal/geometry.h"
#include "gmml/internal/graph.h"
#include "gmml/internal/library_file.h"
#include "gmml/internal/parameter_file.h"
#include "gmml/internal/pdb_file.h"
#include "gmml/internal/pdb_file_builder.h"
#include "gmml/internal/pdb_file_structure.h"
#include "gmml/internal/prep_file.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/sander_minimize.h"
#include "utilities.h"

using std::deque;
using std::list;
using std::map;
using std::string;
using std::vector;

namespace gmml {

struct Structure::InternalResidue : public Residue {
    InternalResidue(const Residue *residue, int start_index)
            : start_index(start_index) {
        clone_from(residue);
    }

    virtual ~InternalResidue() {}

    int start_index;
};

Structure::Structure() : atoms_(0), bonds_(new Graph), head_(-1), tail_(-1) {}

Structure::Structure(const Residue *residue) : atoms_(0), bonds_(new Graph),
                                               head_(-1), tail_(-1) {
    append(residue);
}

Structure::~Structure() {
    if (bonds_ != NULL)
        delete bonds_;
    std::for_each(residues_.begin(), residues_.end(), DeletePtr());
}

// Construction methods
Structure *Structure::clone() const {
    Structure *structure = new Structure;
    structure->clone_from(*this);
    return structure;
}

Structure *Structure::build_from_pdb(const PdbFile& pdb_file,
                                     const PdbMappingInfo& mapping_info) {
    return PdbFileStructure::build(pdb_file, mapping_info);
}

Structure *Structure::build_from_pdb(const PdbFile& pdb_file) {
    return PdbFileStructure::build(pdb_file);
}

// Bonding-related operations
void Structure::add_bond(int atom1, int atom2) {
    bonds_->add_edge(atom1, atom2);
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

Graph *Structure::get_residue_link_table() const {
    vector<size_t> *residue_index_table = get_residue_index_table();

    Graph *links = new Graph(residues_.size());
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

Graph *Structure::get_link_graph() const {
    Graph *graph = new Graph(residues_.size());
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

vector<int> *Structure::get_flexible_linkages() const {
    vector<int> *residues = new vector<int>;
    for (int i = 0; i < residues_.size(); i++) {
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

// Geometric operations
void Structure::shift(double x, double y, double z) {
    for (iterator it = begin(); it != end(); ++it)
        (*it)->translate(x, y, z);
}

void Structure::translate_residue(int residue_index, double x, double y,
                                  double z) {
    for (Residue::iterator it = residues(residue_index)->begin();
            it != residues(residue_index)->end(); ++it) {
        (*it)->translate(x, y, z);
    }
}

void Structure::translate_residues_after(int index, double x, double y,
                                         double z) {
    int num_residues = residue_count();
    for (int i = index; i < num_residues; i++) {
        translate_residue(i, x, y, z);
    }
}

void Structure::set_dihedral(size_t atom1, size_t atom2, size_t atom3,
                             size_t atom4, double degrees) {
    vector<size_t> *atoms = bonds_->edge_bfs(atom2, atom3);

    double current_dihedral = measure(atoms_[atom1]->coordinate(),
                                      atoms_[atom2]->coordinate(),
                                      atoms_[atom3]->coordinate(),
                                      atoms_[atom4]->coordinate());

    RotationMatrix matrix(atoms_[atom2]->coordinate(),
                          Vector<3>(atoms_[atom3]->coordinate(),
                                    atoms_[atom2]->coordinate()),
                          current_dihedral - to_radians(degrees));

    for (vector<size_t>::iterator it = atoms->begin();
            it != atoms->end(); ++it)
        matrix.apply(atoms_[*it]->mutable_coordinate());

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

void Structure::set_residue_angle(int atom1, int atom2, int atom3,
                                  int residue_index, double radians) {
    int start_atom = residues_[residue_index]->start_index;
    int end_atom = size() - 1;

    set_angle_in_range(start_atom, end_atom, atom1, atom2, atom3, radians);
}

void Structure::set_angle_in_range(int start_atom, int end_atom,
                                   int atom1, int atom2, int atom3,
                                   double radians) {
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

    for (int i = start_atom; i <= end_atom; i++)
        matrix.apply(atoms(i)->mutable_coordinate());
}

void Structure::set_angle_after(int residue, int atom1, int atom2, int atom3,
                                double radians) {
    int num_residues = residue_count();
    for (int i = residue; i < num_residues; i++) {
        set_residue_angle(atom1, atom2, atom3, i, radians);
    }
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

// Augmenting operations
int Structure::append(const Structure *structure) {
    if (structure->residue_count() == 0)
        return -1;
    int first_residue = residue_count();
    for (int i = 0; i < structure->residue_count(); i++) {
        // We don't want to use the bonding information in the residues.
        append(structure->residues(i), false);
    }
    
    if (structure->bonds() == NULL) {
        bonds_->append(Graph(structure->size()));
    } else if (structure->bonds()->size() != structure->size()) {
        // Error message here maybe? See below.
        bonds_->append(Graph(structure->size()));
    } else {
        bonds_->append(*structure->bonds());
    }

    return first_residue;
}

int Structure::append(const Residue *residue, bool load_bonds) {
    int prev_size = size();

    InternalResidue *new_residue = new InternalResidue(residue, size());
    residues_.push_back(new_residue);

    if (load_bonds) {
        if (new_residue->bonds() == NULL) {
            bonds_->append(Graph(new_residue->size()));
        } else if (new_residue->bonds()->size() != new_residue->size()) {
            // The user is trying to append a residue whose bond graph is a
            // different size than the number of atoms in the residue. Maybe
            // some error should go here?
            bonds_->append(Graph(new_residue->size()));
        } else {
            bonds_->append(*residue->bonds());
        }
    }

    atoms_.insert(end(), new_residue->begin(), new_residue->end());

    set_head(prev_size + new_residue->head());
    set_tail(prev_size + new_residue->tail());

    return residue_count() - 1;
}

int Structure::append(map<int, int>& atom_map, const IndexedResidue *residue) {
    int prev_size = size();
    int index = append(residue);
    for (int i = 0; i < residue->size(); i++) {
        atom_map[residue->get_atom_index(i)] = prev_size + i;
    }
    return index;
}

int Structure::append(const string& name) {
    Structure *structure = build(name);
    if (structure == NULL)
        return -1;
    int index = append(structure);
    delete structure;
    return index;
}









/*
int Structure::attach(const string& code, const string& head_name,
                      int residue, const string& tail_name) {
    int tail_index = get_atom_index(residue, tail_name);

    return attach(code, head_name, tail_index);


    //Residue *new_residue = build_prep_file(code);
    //if (new_residue == NULL)
    //    return -1;
    //return attach(new_residue, head_name, residue, tail_name);
}


int Structure::attach(const string& code, int head_index, int tail_index) {
    Residue *new_residue = build_prep_file(code);
    return attach(new_residue, head_index, tail_index);
}
*/

/*
int Structure::attach(const string& code, int tail) {
    Residue *new_residue = build_prep_file(code);
    if (new_residue == NULL)
        return -1;
    int head = new_residue->head();
    if (head == -1)
        return -1
    return attach(new_residue, head, tail);
}

int Structure::attach(const string& code) {
    int tail_index = tail();

    if (tail_index == -1)
        return -1;
    
    return attach(code, tail_index);
}
*/

int Structure::attach(const Structure *structure) {
    int tail_index = tail();
    if (tail_index == -1)
        return -1;
    return attach_from_head(structure, tail_index);
}

int Structure::attach(const Residue *residue) {
    Structure structure(residue);
    return attach(&structure);
}

int Structure::attach(const std::string& code) {
    Structure *new_structure = build(code);
    if (new_structure == NULL)
        return -1;
    int ret_val = attach(new_structure);
    delete new_structure;
    return ret_val;
}

int Structure::attach(const Structure *structure,
                      int head_residue, const string& head_name,
                      int tail_residue, const string& tail_name) {
    int head_index = structure->get_atom_index(head_residue, head_name);
    int tail_index = get_atom_index(tail_residue, tail_name);
    if (head_index == -1 || tail_index == -1)
        return -1;
    return attach(structure, head_index, tail_index);
}

int Structure::attach(const Residue *residue, const string& head_name,
                      int target_residue, const string& tail_name) {
    int head_index = residue->get_index(head_name);
    int tail_index = get_atom_index(target_residue, tail_name);
    Structure structure(residue);
    return attach(&structure, head_index, tail_index);
}

int Structure::attach(const Structure *structure, int head, int tail) {
    return StructureAttach()(*this, structure, head, tail);
}





int Structure::attach_from_head(const Structure *structure, int target_residue,
                                const string& tail_name) {
    int tail_atom = get_atom_index(target_residue, tail_name);
    if (tail_atom == -1)
        return -1;
    return attach_from_head(structure, tail_atom);
}

int Structure::attach_from_head(const Structure *structure, int tail_atom) {
    int head_atom = structure->head();
    if (head_atom == -1)
        return -1;
    return attach(structure, head_atom, tail_atom);
}



int Structure::attach_to_tail(const Structure *structure, int head_residue,
                              const string& head_name) {
    int head_index = structure->get_atom_index(head_residue, head_name);
    return attach_to_tail(structure, head_index);
}

int Structure::attach_to_tail(const Residue *residue, const string& head_name) {
    int head_index = residue->get_index(head_name);
    if (head_index == -1)
        return -1;
    Structure structure(residue);
    return attach_to_tail(&structure, head_index);
}

int Structure::attach_to_tail(const Structure *structure, int head_index) {
    int tail_index = tail();
    if (tail_index == -1)
        return -1;
    return attach(structure, head_index, tail_index);
}




// Removal operations
void Structure::remove_residues(const std::vector<int>& indices) {
    vector<size_t> atoms;
    vector<int> new_indices(residue_count());
    for (int i = 0; i < indices.size(); i++)
        new_indices[indices[i]] = -1;

    deque<bool> is_deleted(size(), false);
    int cur_shift = 0;
    int i = 0;
    int k = 0;
    while (i < residue_count()) {
        if (new_indices[k] == -1) {
            int start_index = residues_[i]->start_index;
            int residue_size = residues(i)->size();
            cur_shift += residue_size;
            for (int j = start_index; j < start_index + residue_size; j++) {
                atoms.push_back(j);
                is_deleted[j] = true;
            }
            delete residues_[i];
            residues_.erase(residues_.begin() + i);
        } else {
            residues_[i]->start_index -= cur_shift;
            i++;
        }
        k++;
    }

    i = 0;
    k = 0;
    while (i < size()) {
        if (is_deleted[k]) {
            atoms_.erase(atoms_.begin() + i);
        } else {
            i++;
        }
        k++;
    }

    bonds_->remove_vertex_list(atoms);
}

// Advanced modification operations
MinimizationResults *Structure::minimize(const string& input_file) {
    return SanderMinimize()(*this, input_file);
}

// Atom-related query operations
int Structure::get_atom_index(int residue_index, int atom_index) const {
    if (atom_index < 0 || atom_index >= residues(residue_index)->size())
        return -1;
    return residues_[residue_index]->start_index + atom_index;
}

int Structure::get_atom_index(int residue_index,
                              const string& atom_name) const {
    if (residue_index >= residue_count() || residue_index < 0)
        return -1;

    int atom_index = residues(residue_index)->get_index(atom_name);
    if (atom_index == -1)
        return -1;

    return get_atom_index(residue_index, atom_index);
}

int Structure::get_anomeric_index(int residue_index) const {
    int residue_start = residues_[residue_index]->start_index;
    int residue_size = residues(residue_index)->size();
    int residue_end = residue_start + residue_size - 1;
    for (int i = residue_start; i <= residue_end; i++) {
        if (atoms(i)->element() == kElementO)
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
        int residue_start = residues_[residue_index]->start_index;
        int residue_size = residues(residue_index)->size();
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

vector<size_t> *Structure::get_residue_index_table() const {
    vector<size_t> *table = new vector<size_t>(size());
    size_t current_index = 0;
    for (size_t i = 0; i < residue_count(); i++)
        for (size_t j = 0; j < residues(i)->size(); j++)
            (*table)[current_index++] = i;
    return table;
}

// Residue-related query operations
const Residue *Structure::residues(int index) const {
    return residues_[index];
}

Residue *Structure::residues(int index) {
    return residues_[index];
}

size_t Structure::get_residue_index(size_t atom_index) const {
    for (size_t i = 1; i < residue_count(); i++) {
        if (atom_index < residues_[i]->start_index)
            return i - 1;
    }
    return residue_count() - 1;
}

int Structure::get_parent_residue(int residue_index) const {
    int parent_atom = get_parent_atom(residue_index);
    if (parent_atom != -1)
        return get_residue_index(parent_atom);
    return -1;
}

// Accessors
const Structure::AdjList& Structure::bonds(size_t index) const {
    return bonds_->edges(index);
}

// Mutators
void Structure::set_bonds(const Graph *bonds) {
    if (bonds_ != NULL)
        delete bonds_;
    bonds_ = bonds->clone();
}

// File operations
PdbFile *Structure::build_pdb_file() const {
    return PdbFileBuilder().build(*this);
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
    if (coordinate_file.coordinate_count() < atoms_.size()) { 
        warning(string("Structure: Insufficient number of coordinates") +
                " in coordinate file");
        return;
    }
    for (int i = 0; i < atoms_.size(); i++)
        atoms_[i]->set_coordinate(coordinate_file[i]);
}

void Structure::print(std::ostream& out) const {
    out << "Total atoms: " << size() << std::endl;
    for (int i = 0; i < residue_count(); i++) {
        InternalResidue *residue = residues_[i];
        out << "Residue: " << residue->name() <<
               ", start - " << residue->start_index <<
               ", size - " << residue->size() << std::endl;
        for (int j = residue->start_index;
                j < residue->start_index + residue->size(); j++) {
            Atom *atom = atoms_[j];
            out << "Atom: name - " << atom->name() <<
                   ", type - " << atom->type() <<
                   ", element - " << static_cast<int>(atom->element()) <<
                   ", coordinate - " << atom->coordinate().x << " " <<
                   atom->coordinate().y << " " << atom->coordinate().z <<
                   std::endl;
        }
    }
}

void Structure::clone_from(const Structure& structure) {
    if (bonds_ != NULL)
        delete bonds_;
    bonds_ = new Graph;

    std::for_each(residues_.begin(), residues_.end(), DeletePtr());
    residues_.clear();
    atoms_.clear();

    head_ = structure.head();
    tail_ = structure.tail();

    append(&structure);
}

struct StructureAttach::Impl {
    static Vector<3> get_connection_direction(const Structure& structure,
                                              int source_index,
                                              int target_index);

    static void set_dihedrals(Structure& structure, int new_residue_index,
                              int target_residue_index, int carbon_number,
                              int oxygen_number);
};

/*
int StructureAttach::operator()(Structure& structure, const string& code,
                                int residue) const {
    Residue *new_residue = build_prep_file(code);
    if (new_residue == NULL)
        return -1;
    int tail = structure.tail();
    int head = new_residue->head();
    if (tail == -1 || head == -1)
        return -1;
    return operator()(structure, new_residue, head, tail);
}

int StructureAttach::operator()(Structure& structure, Residue *new_residue,
                                const string& head_atom, int residue,
                                const string& tail_atom) const {
    int head_index = new_residue->get_index(head_atom);
    int tail_index = structure.get_atom_index(residue, tail_atom);
    return operator()(structure, new_residue, head_index, tail_index);
}
*/

int StructureAttach::operator()(Structure& structure,
                                const Structure *new_structure,
                                int head, int tail) const {
    const Structure::AtomList& atoms = structure.atoms();

    int prev_size = structure.size();


    // THE FIRST RESIDUE IN THE NEW STRUCTURE MIGHT NOT BE THE ONE YOU'RE
    // ATTACHING TO. NEED TO FIX.
    int structure_head = structure.head();
    int new_residue_index = structure.append(new_structure);

    // Append modifies the head atom. This is not what we want on attach.
    structure.set_head(structure_head);
    if (new_residue_index == -1)
        return -1;

    int new_atom_index = prev_size + head;

    int target_atom_index = tail;

    int residue_index = structure.get_residue_index(tail);

    string new_atom_name = structure.atoms(new_atom_index)->name();
    
    string target_atom_name = structure.atoms(tail)->name();





    structure.add_bond(new_atom_index, target_atom_index);

    const ParameterFileSet *parm_set = kDefaultEnvironment.parm_set();

    const ParameterFileBond *parameter_bond =
        parm_set->lookup(atoms[target_atom_index]->type(),
                         atoms[new_atom_index]->type());

    // This constant should be changed and made static const in the class.
    double bond_length = 1.4;
    if (parameter_bond != NULL && parameter_bond->length != kNotSet)
        bond_length = parameter_bond->length;

    Vector<3> direction = Impl::get_connection_direction(structure,
                                                         new_atom_index,
                                                         target_atom_index);
    direction.normalize();
    direction *= bond_length;

    direction += Vector<3>(atoms[new_atom_index]->coordinate());
    Vector<3> oxygen_position = direction;
    VectorBase<3> offset(Vector<3>(atoms[target_atom_index]->coordinate()) -
                         oxygen_position);
    structure.translate_residues_after(new_residue_index, offset[0], offset[1],
                                       offset[2]);

    // Set the bond angles
    const Structure::AdjList& adj_atoms = structure.bonds(target_atom_index);
    for (int i = 0; i < adj_atoms.size(); i++) {
        if (adj_atoms[i] != new_atom_index) {
            int third_atom = adj_atoms[i];
            const ParameterFileAngle *parameter_angle = parm_set->lookup(
                atoms[third_atom]->type(),
                atoms[target_atom_index]->type(),
                atoms[new_atom_index]->type());

            if (parameter_angle != NULL)
                structure.set_angle_after(new_residue_index,
                                          third_atom, target_atom_index,
                                          new_atom_index,
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
        Impl::set_dihedrals(structure, new_residue_index, residue_index,
                            carbon_number, oxygen_number);

    return new_residue_index;
}

Vector<3> StructureAttach::Impl::get_connection_direction(
        const Structure& structure, int source_index, int target_index) {
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

void StructureAttach::Impl::set_dihedrals(Structure& structure, 
                                          int new_residue_index,
                                          int target_residue_index,
                                          int carbon_number,
                                          int oxygen_number) {
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
        structure.set_dihedral(target_residue_index, "H" + oxygen,
                               target_residue_index, "C" + oxygen,
                               target_residue_index, "O" + oxygen,
                               new_residue_index, "C" + carbon,
                               0.0);
    }

    // Regardless of what the actual phi torsion is defined to be, setting
    // C(x+1)-Cx-O'-C' to 180.0 indirectly sets phi to the right thing.
    // So we set this torsion instead of calling set_phi(), which 
    // sets the actual phi torsion.
    structure.set_dihedral(target_residue_index, "C" + oxygen,
                           target_residue_index, "O" + oxygen,
                           new_residue_index, "C" + carbon,
                           new_residue_index,
                           "C" + to_string(carbon_number + 1),
                           180.0);

    // Should probably change this to check if the linkage is exocyclic
    if (oxygen_number == 6) {
        structure.set_dihedral(target_residue_index,
                              "O" + to_string(oxygen_number - 1),
                              target_residue_index,
                              "C" + to_string(oxygen_number - 1),
                              target_residue_index, "C" + oxygen,
                              target_residue_index, "O" + oxygen,
                              60.0);
    }
}

}  // namespace gmml
