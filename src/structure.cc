// Author: Robert Davis

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
#include "gmml/internal/glycam_code_set.h"
#include "gmml/internal/graph.h"
#include "gmml/internal/library_file.h"
#include "gmml/internal/parameter_file.h"
#include "gmml/internal/pdb_file.h"
#include "gmml/internal/pdb_file_builder.h"
#include "gmml/internal/prep_file.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/sander_minimize.h"
#include "gmml/internal/stubs/logging.h"
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

// Bonding-related operations
bool Structure::is_bonded(int atom1, int atom2) const {
    return bonds_->is_adjacent(atom1, atom2);
}

bool Structure::add_bond(int atom1, int atom2) {
    return bonds_->add_edge(atom1, atom2);
}

bool Structure::remove_bond(int atom1, int atom2) {
    return bonds_->remove_edge(atom1, atom2);
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

// Augmenting operations
int Structure::append(const Structure *structure) {
    if (structure->residue_count() == 0)
        return -1;
    int first_residue = residue_count();
    int prev_size = size();
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

    // Appending the residues modifies the head and tail of the structures based
    // on the heads and tails of the residues. We want to instead use the
    // head and tail of the whole structure.
    set_head(prev_size + structure->head());
    set_tail(prev_size + structure->tail());

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

int Structure::attach_from_head(const string& code, const string& tail_atom) {
    Structure *structure = build(code);
    if (structure == NULL)
        return -1;
    int tail_residue = get_residue_index(tail());
    return attach_from_head(structure, tail_residue, tail_atom);
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

void Structure::remove_residue(int index) {
    vector<int> indices(1);
    indices[0] = index;
    remove_residues(indices);
    // Do this?
    //remove_residues(vector<int>(1, index));
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

vector<int> *Structure::get_adjacent_residues_by_atom(int atom_index) const {
    if (atom_index < 0 || atom_index >= size()) {
        return NULL;
    }
    int this_residue = get_residue_index(atom_index);
    vector<int> *residues = new vector<int>;
    const AdjList& adj_list = bonds(atom_index);
    for (int i = 0; i < adj_list.size(); i++) {
        int residue_index = get_residue_index(adj_list[i]);
        if (residue_index != this_residue) {
            residues->push_back(residue_index);
        }
    }
    return residues;
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
        const ParameterSet& parm_set) const {
    AmberTopBuilder builder(parm_set);
    return builder.build(*this);
}

AmberTopFile *Structure::build_amber_top_file() const {
    AmberTopBuilder builder;
    return builder.build(*this);
}

void Structure::print_amber_top_file(const string& file_name,
                                     const ParameterSet& parm_set) const {
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
        LOG(WARNING) << "Insufficient number of coordinates in " <<
                        "in coordinate file, cannot load.";
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

    // Append will modify the head and tail of this structure based on the
    // heads and tails of the residues of the structure. We set the head and
    // tail after appending to undo this.
    append(&structure);

    head_ = structure.head();
    tail_ = structure.tail();
}

struct StructureAttach::Impl {
    static Vector<3> get_connection_direction(const Structure& structure,
                                              int source_index,
                                              int target_index);
};

// TODO: Fix this up and figure out precisely how this should be done.
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

    const ParameterSet *parm_set = kDefaultEnvironment.parm_set();

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

            if (parameter_angle != NULL) {
                structure.set_angle_after(new_residue_index,
                                          third_atom, target_atom_index,
                                          new_atom_index,
                                          to_radians(parameter_angle->angle));
            }
        }
    }

    // This sets torsions appropriately if the residues are sugars. I'm not
    // sure how to set them otherwise.
    int carbon_number = kNotSet;
    if (new_atom_name.size() > 1 && is_number(new_atom_name[1]))
        carbon_number = char_to_number(new_atom_name[1]);

    int oxygen_number = kNotSet;
    if (target_atom_name.size() > 1 && is_number(target_atom_name[1]))
        oxygen_number = char_to_number(target_atom_name[1]);
    else
        oxygen_number = 1;

    if (is_set(carbon_number) && is_set(oxygen_number))
        carbohydrate::set_default_torsions(&structure, new_residue_index,
                                           residue_index, carbon_number,
                                           oxygen_number);

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

}  // namespace gmml
