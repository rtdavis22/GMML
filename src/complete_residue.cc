// Author: Robert Davis

#include "gmml/internal/complete_residue.h"

#include <iostream>
#include <vector>

#include "gmml/internal/environment.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/structure.h"
#include "utilities.h"

using std::string;
using std::vector;

namespace gmml {
namespace {

// Returns true only if the input coordinate isn't too close to any of the
// coordinates in the list.
bool check_candidate(const Coordinate& coordinate,
                     const vector<Coordinate>& coordinates) {
    static double threshold = 1.0;

    for (int i = 0; i < coordinates.size(); i++) {
        if (measure(coordinate, coordinates[i]) < threshold) {
            return false;
        }
    }
    return true;
}

}  // namespace

// Private implementation
struct CompleteResidue::Impl {
    typedef Residue::AtomPtr AtomPtr;

    // This function attempts to find a chain of up to three atoms with known
    // coordinates that starts with the atom at the given index.
    // Any NULL (AtomPtr()) atoms in the list indicate unknown atoms, which
    // can't be use to set the current atom's coordinate.
    static vector<int> *find_chain(int index, const vector<AtomPtr> *atoms,
                                   const Graph *bonds);

    // The function uses the chain atoms calculated in the previous function
    // to find a coordinate for the atom with the given index, using the
    // internal coordinates in the structure.
    static Coordinate place_coordinate(int index, const vector<AtomPtr> *atoms,
                                       const vector<int>& chain,
                                       const Structure *structure);

    // Calculate the coordinate based only on a distance and a reference atom.
    // Any coordinate on a sphere around the reference coordinate is a valid
    // return of this function.
    static Coordinate calculate_unknown_coordinate(const Coordinate& parent,
                                                   double distance);

    // Calculate the coordinate based on two reference atoms, a distance, and
    // an angle. The return coordinate will be |distance| units from parent1,
    // and the angle formed by the returned coordinate, parent1, and parent2
    // will be |angle| radians. As in the previous function, there are
    // infinitely many possible return values.
    static Coordinate calculate_unknown_coordinate(const Coordinate& parent1,
                                                   const Coordinate& parent2,
                                                   double distance,
                                                   double angle);

    // Calculate the coordinate that is |distance| units from parent1, that
    // forms an angle of |angle| radians with parent1 and parent2, and that
    // forms a dihedral of |dihedral| radians with parent1, parent2, and
    // parent3. Note that there is only one such coordinate. If this coordinate
    // is too close to any of the given atoms to avoid, the dihedral will
    // be adjusted 120 degrees.
    static Coordinate calculate_unknown_coordinate(
            const Coordinate& parent1, const Coordinate& parent2,
            const Coordinate& parent3, double distance, double angle,
            double dihedral, vector<Coordinate> *atoms_to_avoid);
};

vector<int> *CompleteResidue::Impl::find_chain(
        int index, const vector<Residue::AtomPtr> *atoms, const Graph *bonds) {
    typedef Structure::AdjList AdjList;

    // We try to find three atoms that will set the space for this atom.
    // If we only find two atoms, we can still use the distance and angle.
    // If we just find one atom, we can use the distance. If we don't find
    // any atoms, we can't do anything.
    int atom1 = -1;
    int atom2 = -1;
    int atom3 = -1;
    int longest_chain = 0;
    const AdjList& atom_bonds = bonds->edges(index);
    for (int j = 0; j < atom_bonds.size(); j++) {
        if ((*atoms)[atom_bonds[j]] == AtomPtr())
            continue;
        if (longest_chain == 0) {
            atom1 = atom_bonds[j];
            longest_chain = 1;
        }
        const AdjList& atom1_bonds = bonds->edges(atom_bonds[j]);
        for (int k = 0; k < atom1_bonds.size(); k++) {
            if ((*atoms)[atom1_bonds[k]] == AtomPtr())
                continue;
            if (longest_chain < 2) {
                atom1 = atom_bonds[j];
                atom2 = atom1_bonds[k];
                longest_chain = 2;
            }
            const AdjList& atom2_bonds = bonds->edges(atom1_bonds[k]);
            for (int l = 0; l < atom2_bonds.size(); l++) {
                if (atom2_bonds[l] == atom_bonds[j] ||
                        (*atoms)[atom2_bonds[l]] == AtomPtr())
                    continue;
                atom1 = atom_bonds[j];
                atom2 = atom1_bonds[k];
                atom3 = atom2_bonds[l];
                longest_chain = 3;
                break;
            }
            if (longest_chain == 3)
                break;
        }
        if (longest_chain == 3)
            break;
    }

    vector<int> *chain = new vector<int>;
    if (atom1 == -1)
        return chain;
    else
        chain->push_back(atom1);

    if (atom2 == -1)
        return chain;
    else
        chain->push_back(atom2);

    if (atom3 != -1)
        chain->push_back(atom3);

    return chain;
}

Coordinate CompleteResidue::Impl::place_coordinate(
        int index, const vector<AtomPtr> *atoms,
        const vector<int>& chain, const Structure *structure) {
    typedef Structure::AdjList AdjList;
    Coordinate coordinate;

    double distance = measure(structure->atoms(chain[0])->coordinate(),
                              structure->atoms(index)->coordinate());

    if (chain.size() >= 2) {
        double angle = measure(structure->atoms(chain[1])->coordinate(),
                               structure->atoms(chain[0])->coordinate(),
                               structure->atoms(index)->coordinate());
        if (chain.size() == 2) {
            coordinate = calculate_unknown_coordinate(
                    (*atoms)[chain[0]]->coordinate(),
                    (*atoms)[chain[1]]->coordinate(),
                    distance, angle);
        } else {
            double dihedral = measure(structure->atoms(chain[2])->coordinate(),
                                      structure->atoms(chain[1])->coordinate(),
                                      structure->atoms(chain[0])->coordinate(),
                                      structure->atoms(index)->coordinate());
            vector<Coordinate> *to_avoid = new vector<Coordinate>;
            const AdjList& atom1_bonds = structure->bonds(chain[0]);
            for (int i = 0; i < atom1_bonds.size(); i++) {
                int atom_index = atom1_bonds[i];
                if ((*atoms)[atom_index] != Residue::AtomPtr()) {
                    to_avoid->push_back((*atoms)[atom_index]->coordinate());
                }
            }
            coordinate = calculate_unknown_coordinate(
                    (*atoms)[chain[0]]->coordinate(),
                    (*atoms)[chain[1]]->coordinate(),
                    (*atoms)[chain[2]]->coordinate(),
                    distance, angle, dihedral, to_avoid);
        }

    } else {
        coordinate = calculate_unknown_coordinate(
                (*atoms)[chain[0]]->coordinate(), distance);
    }

    return coordinate;
}

Coordinate CompleteResidue::Impl::calculate_unknown_coordinate(
        const Coordinate& parent, double distance) {
    Coordinate position(parent);
    position.x += distance;
    return position;
}

Coordinate CompleteResidue::Impl::calculate_unknown_coordinate(
        const Coordinate& parent1, const Coordinate& parent2,
        double distance, double angle) {
    return calculate_point(Coordinate(0.0, 0.0, 0.0), parent2, parent1,
                           angle, 0.0, distance);
}

// I need to make sure that the chirality isn't inverted.
Coordinate CompleteResidue::Impl::calculate_unknown_coordinate(
        const Coordinate& parent1, const Coordinate& parent2,
        const Coordinate& parent3, double distance, double angle,
        double dihedral, vector<Coordinate> *atoms_to_avoid) {
    Coordinate candidate = calculate_point(parent3, parent2, parent1,
                                           angle, dihedral, distance);
    if (check_candidate(candidate, *atoms_to_avoid))
        return candidate;

    // Add 120 degrees to the dihedral.
    dihedral += 2.0*kPi/3.0;
    if (dihedral > kPi)
        dihedral -= 2.0*kPi;
    candidate = calculate_point(parent3, parent2, parent1, angle, dihedral,
                                distance);
    if (check_candidate(candidate, *atoms_to_avoid))
        return candidate;

    dihedral += 2.0*kPi/3.0;
    if (dihedral > kPi)
        dihedral -= 2.0*kPi;
    candidate = calculate_point(parent3, parent2, parent1, angle, dihedral,
                                distance);

    return candidate;
}


// Public implementation
bool CompleteResidue::operator()(Residue *residue, const string& name) const {
    const Structure *complete_structure = build(name);
    if (complete_structure == NULL)
        return false;
    return operator()(residue, complete_structure);
}

bool CompleteResidue::operator()(Residue *residue,
                                 const Structure *complete_structure) const {
    typedef Residue::AtomPtr AtomPtr;

    if (complete_structure->size() < residue->size())
        return false;

    vector<AtomPtr> *new_atoms = new vector<AtomPtr>(complete_structure->size(),
                                                     AtomPtr());

    for (int i = 0; i < residue->size(); i++) {
        bool found = false;
        for (int j = 0; j < complete_structure->size(); j++) {
            AtomPtr atom_ptr = residue->atoms(i);
            if (residue->atoms(i)->name() ==
                    complete_structure->atoms(j)->name()) {
                (*new_atoms)[j] = residue->atoms(i);
                found = true;
                break;
            }
        }
        if (!found)
            return false;
    }

    for (int i = 0; i < new_atoms->size(); i++) {
        if ((*new_atoms)[i] != AtomPtr()) {
            (*new_atoms)[i]->set_type(complete_structure->atoms(i)->type());
            (*new_atoms)[i]->set_charge(complete_structure->atoms(i)->charge());
            (*new_atoms)[i]->set_element(
                    complete_structure->atoms(i)->element());
            continue;
        }

        vector<int> *chain = Impl::find_chain(i, new_atoms, 
                                              complete_structure->bonds());
        if (chain->size() == 0)
            return false;

        Coordinate coordinate = Impl::place_coordinate(i, new_atoms, *chain,
                                                       complete_structure);

        AtomPtr new_atom(complete_structure->atoms(i)->clone());
        new_atom->set_coordinate(coordinate);
        (*new_atoms)[i] = new_atom;
    }

    residue->set_atoms(new_atoms);
    residue->set_bonds(complete_structure->bonds()->clone());
    return true;
}

}  // namespace gmml
