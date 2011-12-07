// Author: Robert Davis

#include "gmml/internal/solvated_structure.h"

#include <cassert>

#include <algorithm>
#include <vector>

#include "gmml/internal/amber_top_builder.h"
#include "gmml/internal/coordinate_file.h"
#include "gmml/internal/coordinate_grid.h"
#include "gmml/internal/geometry.h"
#include "utilities.h"

using std::map;
using std::vector;

namespace gmml {

using detail::TrimmedSolvents;

SolvatedStructure::SolvatedStructure(const Structure& structure,
                                     const Structure& solvent,
                                     double distance,
                                     double closeness) : BoxedStructure() {
    Structure::clone_from(structure);
    BoxedRegion *solvent_region = get_boxed_region(solvent);
    solvate(solvent, solvent_region, distance, closeness);
    delete solvent_region;
}


SolvatedStructure::SolvatedStructure(const Structure& structure,
                                     const BoxedStructure& solvent,
                                     double distance,
                                     double closeness)
        : BoxedStructure(), last_solute_atom_(structure.size() - 1) {
    Structure::clone_from(structure);
    BoxedRegion *solvent_region = get_boxed_region(solvent);
    const Box *box = solvent.box();
    // If we have box dimensions, we contract the solvent region so that it's
    // dimensions are the same as the box's dimensions. If the dimensions of the
    // box are larger than the solvent, the contraction is actually an
    // expansion. Currently, the angle of the box is not used and presumed
    // to be 90 degrees.
    if (box != NULL) {
        double delta_x = (solvent_region->max_x - solvent_region->min_x -
                          box->length)/2.0;
        double delta_y = (solvent_region->max_y - solvent_region->min_y -
                          box->width)/2.0;
        double delta_z = (solvent_region->max_z - solvent_region->min_z -
                          box->height)/2.0;
        solvent_region->contract(delta_x, delta_y, delta_z);
    }
    solvate(solvent, solvent_region, distance, closeness);
    delete solvent_region;
}

void SolvatedStructure::solvate(const Structure& solvent,
                                BoxedRegion *solvent_region,
                                double distance, double closeness) {
    BoxedRegion *solute_region = get_boxed_region(*this);
    solute_region->expand(distance, distance, distance);

    // The dimensions of the solvent region. These don't necessarily
    // correspond to the atoms in the solvent.
    double solvent_length = solvent_region->max_x - solvent_region->min_x;
    double solvent_width = solvent_region->max_y - solvent_region->min_y;
    double solvent_height = solvent_region->max_z - solvent_region->min_z;

    // The dimensions of the solute.
    double solute_length = solute_region->max_x - solute_region->min_x;
    double solute_width = solute_region->max_y - solute_region->min_y;
    double solute_height = solute_region->max_z - solute_region->min_z;

    // These are the total number of copies we will make in each direction,
    // including trimmed boxes.
    int x_copies = solute_length/solvent_length + 1;
    int y_copies = solute_width/solvent_width + 1;
    int z_copies = solute_height/solvent_height + 1;

    // The amounts to trim off each side of the solvent box.
    double x_to_trim = solvent_length*x_copies - solute_length;
    double y_to_trim = solvent_width*y_copies - solute_width;
    double z_to_trim = solvent_height*z_copies - solute_height;

    // We first shift the solvent and solvent box to align them on top of the
    // solute.
    double shift_x = solute_region->min_x - solvent_region->min_x;
    double shift_y = solute_region->min_y - solvent_region->min_y;
    double shift_z = solute_region->min_z - solvent_region->min_z;
    Structure *solvent_clone = solvent.clone();
    solvent_clone->shift(shift_x, shift_y, shift_z);
    solvent_region->shift(shift_x, shift_y, shift_z);

    delete solute_region;

    // We precompute the trimmed solvent boxes around the perimeter.
    TrimmedSolvents trimmed_solvents(*solvent_clone, solvent_region,
                                     x_to_trim, y_to_trim, z_to_trim);

    // This is where the solvent gets added. You can unroll the loops to get
    // rid of the boolean assignments, but it doesn't buy you much and it gets
    // a little ugly.
    for (int i = 0; i < x_copies; i++) {
        bool trim_x = (i == x_copies - 1);
        for (int j = 0; j < y_copies; j++) {
            bool trim_y = (j == y_copies - 1);
            for (int k = 0; k < z_copies; k++) {
                bool trim_z = (k == z_copies - 1);
                Structure *solvent_copy;
                if (trim_x || trim_y || trim_z) {
                    solvent_copy = trimmed_solvents.get_trimmed_solvent(
                            trim_x, trim_y, trim_z)->clone();
                } else {
                    solvent_copy = solvent_clone->clone();
                }
                solvent_copy->shift(i*solvent_length, j*solvent_width,
                                    k*solvent_height);
                append(*solvent_copy);
                delete solvent_copy;
            }
        }
    }

    delete solvent_clone;

    remove_close_solvent_residues(closeness);

    BoxedRegion *new_region = get_boxed_region(*this);
    if (box_ == NULL) {
        box_ = new Box(kNotSet, kNotSet, kNotSet, kNotSet);
    }
    box_->angle = 90.0;
    box_->length = new_region->max_x - new_region->min_x;
    box_->width = new_region->max_y - new_region->min_y;
    box_->height = new_region->max_z - new_region->min_z;
    delete new_region;
}

void SolvatedStructure::remove_close_solvent_residues(double closeness) {
    // We insert all solvent atoms into a coordinate grid.
    CoordinateGrid<int> grid(closeness);
    for (int i = last_solute_atom_ + 1; i < atoms_.size(); i++) {
        grid.insert(atoms_[i]->coordinate(), i);
    }

    vector<size_t> *residue_index_table =  get_residue_index_table();

    // The indices of all residues that are too close.
    vector<int> to_remove;

    // For each solute atom, we look for all solvent atoms that map to the same
    // cell or an adjacent cell. Since the dimension of the cells in the grid is
    // our closeness parameter, only these atoms have a chance at being too
    // close to the solute atom.
    for (int i = 0; i <= last_solute_atom_; i++) {
        const Coordinate& coordinate = atoms_[i]->coordinate();
        vector<int> *found = grid.retrieve_adjacent_cells(coordinate);
        for (int j = 0; j < found->size(); j++) {
            int solvent_atom = (*found)[j];
            double distance = measure(coordinate,
                                      atoms_[solvent_atom]->coordinate());
            // It's possible we insert the same residue twice, but that's OK.
            if (distance < closeness)
                to_remove.push_back((*residue_index_table)[solvent_atom]);
        }
        delete found;
    }
    delete residue_index_table;
    remove_residues(to_remove);
}

AmberTopFile *SolvatedStructure::build_amber_top_file() const {
    AmberTopBuilder builder;
    return builder.build(*this);
}

CoordinateFile *SolvatedStructure::build_coordinate_file() const {
    CoordinateFile *file = Structure::build_coordinate_file();
    file->add_noncoordinate(box_->length, box_->width, box_->height);
    file->add_noncoordinate(box_->angle, box_->angle, box_->angle);
    return file;
}

namespace detail {
namespace {

// This is used for sorting a list of atom indices by either the atoms
// x coordinate, y coordinate, or z coordinate.
class CoordinateLess {
  public:
    enum Dimension { kDimensionX, kDimensionY, kDimensionZ };

    CoordinateLess(const Structure::AtomList& atoms, Dimension dimension)
            : atoms_(atoms), dimension_(dimension) {}

    bool operator()(int index1, int index2) {
        const Coordinate& coordinate1 = atoms_[index1]->coordinate();
        const Coordinate& coordinate2 = atoms_[index2]->coordinate();
        switch (dimension_) {
            case kDimensionX:
                return coordinate1.x < coordinate2.x;
            case kDimensionY:
                return coordinate1.y < coordinate2.y;
            case kDimensionZ:
                return coordinate1.z < coordinate2.z;
            default:
                assert(false);
        }
    }

  private:
    const Structure::AtomList& atoms_;
    Dimension dimension_;
};

}  // namespace

TrimmedSolvents::TrimmedSolvents(const Structure& structure,
                                 const BoxedRegion *solvent_box,
                                 double trim_x, double trim_y, double trim_z)
        : trimmed_x_(structure.clone()), trimmed_y_(structure.clone()),
          trimmed_z_(structure.clone()), trimmed_xy_(structure.clone()),
          trimmed_xz_(structure.clone()), trimmed_yz_(structure.clone()),
          trimmed_xyz_(structure.clone()) {
    // We construct 3 lists of atom indices: one sorted by the atom's
    // x coordinate, one sorted by the y coordinate, and one by the z
    // coordinate.
    vector<int> sorted_by_x(structure.size());
    vector<int> sorted_by_y(structure.size());
    vector<int> sorted_by_z(structure.size());

    for (int i = 0; i < sorted_by_x.size(); i++) {
        sorted_by_x[i] = sorted_by_y[i] = sorted_by_z[i] = i;
    }

    const Structure::AtomList& atoms = structure.atoms();
    CoordinateLess XCoordinateLess(atoms, CoordinateLess::kDimensionX);
    CoordinateLess YCoordinateLess(atoms, CoordinateLess::kDimensionY);
    CoordinateLess ZCoordinateLess(atoms, CoordinateLess::kDimensionZ);

    std::sort(sorted_by_x.begin(), sorted_by_x.end(), XCoordinateLess);
    std::sort(sorted_by_y.begin(), sorted_by_y.end(), YCoordinateLess);
    std::sort(sorted_by_z.begin(), sorted_by_z.end(), ZCoordinateLess);

    vector<size_t> *residue_index_table = structure.get_residue_index_table();

    vector<int>::iterator it;

    // We find all residues that have an atom whose x, y, or z coordinate is
    // too large.
    it = sorted_by_x.end() - 1;
    vector<int> x_to_remove;
    double x_cutoff = solvent_box->max_x - trim_x;
    while (atoms[*it]->coordinate().x > x_cutoff) {
        x_to_remove.push_back((*residue_index_table)[*it]);
        --it;
    }

    it = sorted_by_y.end() - 1;
    vector<int> y_to_remove;
    double y_cutoff = solvent_box->max_y - trim_y;
    while (atoms[*it]->coordinate().y > y_cutoff) {
        y_to_remove.push_back((*residue_index_table)[*it]);
        --it;
    }

    it = sorted_by_z.end() - 1;
    vector<int> z_to_remove;
    double z_cutoff = solvent_box->max_z - trim_z;
    while (atoms[*it]->coordinate().z > z_cutoff) {
        z_to_remove.push_back((*residue_index_table)[*it]);
        --it;
    }

    delete residue_index_table;

    create_trimmed_copies(x_to_remove, y_to_remove, z_to_remove);
}

TrimmedSolvents::~TrimmedSolvents() {
    delete trimmed_x_;
    delete trimmed_y_;
    delete trimmed_z_;
    delete trimmed_xy_;
    delete trimmed_xz_;
    delete trimmed_yz_;
    delete trimmed_xyz_;
}

const Structure *TrimmedSolvents::get_trimmed_solvent(
        bool is_x_trimmed, bool is_y_trimmed, bool is_z_trimmed) const {
    // These are just aliases for readability.
    bool x = is_x_trimmed;
    bool y = is_y_trimmed;
    bool z = is_z_trimmed;
    if (x && !y && !z)
        return trimmed_x_;
    else if (!x && y && !z)
        return trimmed_y_;
    else if (!x && !y && z)
        return trimmed_z_;
    else if (x && y && !z)
        return trimmed_xy_;
    else if (x && !y && z)
        return trimmed_xz_;
    else if (!x && y && z)
        return trimmed_yz_;
    else if (x && y && z)
        return trimmed_xyz_;
    else
        return NULL;
}

void TrimmedSolvents::create_trimmed_copies(const vector<int>& x_to_remove,
                                            const vector<int>& y_to_remove,
                                            const vector<int>& z_to_remove) {
    trimmed_x_->remove_residues(x_to_remove);
    trimmed_y_->remove_residues(y_to_remove);
    trimmed_z_->remove_residues(z_to_remove);

    vector<int> xy_to_remove(x_to_remove);
    xy_to_remove.insert(xy_to_remove.end(), y_to_remove.begin(),
                        y_to_remove.end());
    trimmed_xy_->remove_residues(xy_to_remove);

    vector<int> yz_to_remove(y_to_remove);
    yz_to_remove.insert(yz_to_remove.end(), z_to_remove.begin(),
                        z_to_remove.end());
    trimmed_yz_->remove_residues(yz_to_remove);

    vector<int> xz_to_remove(x_to_remove);
    xz_to_remove.insert(xz_to_remove.end(), z_to_remove.begin(),
                        z_to_remove.end());
    trimmed_xz_->remove_residues(xz_to_remove);

    vector<int> xyz_to_remove(xy_to_remove);
    xyz_to_remove.insert(xyz_to_remove.end(), z_to_remove.begin(),
                         z_to_remove.end());
    trimmed_xyz_->remove_residues(xyz_to_remove);
}

}  // namespace detail
}  // namespace gmml
