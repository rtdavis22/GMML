// Author: Robert Davis

#ifndef GMML_INTERNAL_COORDINATE_GRID_H_
#define GMML_INTERNAL_COORDINATE_GRID_H_

#include <map>
#include <utility>
#include <vector>

#include "gmml/internal/geometry.h"
#include "gmml/internal/stubs/common.h"
#include "gmml/internal/stubs/utils.h"

namespace gmml {

// This class maps coordinates to 3-dimensional cells in space and stores
// information associated with the coordinate (such as an index, pointer,
// or iterator) in the cell. The cells are stored in such a way that
// querying neighbor cells is efficient.
template<typename T>
class CoordinateGrid {
  public:
    // This is the storage type for each cell.
    typedef typename std::vector<T>* StorageType;

    // The constructor creates a grid in which cells have a length, width,
    // and height of grid_unit.
    explicit CoordinateGrid(double grid_unit) : grid_unit_(grid_unit) {}

    virtual ~CoordinateGrid();

    // This maps the coordinate to a cell and stores the associated
    // information in the cell. The contents of the cell are returned.
    const StorageType insert(const Coordinate& coordinate, const T& info);

    // This returns the contents of the cell associated with the
    // coordinate.
    const StorageType retrieve(const Coordinate& coordinate) const;

    // This function finds the cell associated with the coordinate and returns
    // the contents of the cell along with the contents of all cells adjacent
    // to the cell (3^3 cells in total). Only coordinates mapped to one of these
    // cells can be within grid_unit distance of the coordinate parameter.
    // It's the caller's responsibility to free the memory that's returned.
    StorageType retrieve_adjacent_cells(const Coordinate& coordinate) const;

    double grid_unit() const { return grid_unit_; }

  private:
    typedef typename std::map<Triplet<int>*, StorageType,
                              TripletPtrLess<int> > GridType;
    typedef typename GridType::iterator iterator;
    typedef typename GridType::const_iterator const_iterator;

    const StorageType retrieve(Triplet<int> *triplet) const;

    GridType grid_;

    // The grid_unit is the length, width, and height of all cells.
    double grid_unit_;

    DISALLOW_COPY_AND_ASSIGN(CoordinateGrid);
};

template<typename T>
CoordinateGrid<T>::~CoordinateGrid() {
    for (iterator it = grid_.begin(); it != grid_.end(); ++it) {
        delete it->first;
        delete it->second;
    }
}

template<typename T>
const typename CoordinateGrid<T>::StorageType CoordinateGrid<T>::insert(
        const Coordinate& coordinate, const T& info) {
    Triplet<int> *triplet = new Triplet<int>(coordinate.x/grid_unit_,
                                             coordinate.y/grid_unit_,
                                             coordinate.z/grid_unit_);

    // The second element of the pair is false if an element already exists
    // with the given key;
    std::pair<iterator, bool> ret = grid_.insert(
            std::make_pair(triplet, static_cast<StorageType>(NULL)));

    if (!ret.second)
        delete triplet;
    else
        ret.first->second = new std::vector<int>;

    ret.first->second->push_back(info);

    return ret.first->second;
}

template<typename T>
const typename CoordinateGrid<T>::StorageType CoordinateGrid<T>::retrieve(
        Triplet<int> *triplet) const {
    StorageType ret = NULL;
    const_iterator it = grid_.find(triplet);
    if (it != grid_.end())
        ret = it->second;
    return ret;
}

template<typename T>
const typename CoordinateGrid<T>::StorageType CoordinateGrid<T>::retrieve(
        const Coordinate& coordinate) const {
    Triplet<int> *triplet = new Triplet<int>(coordinate.x/grid_unit_,
                                             coordinate.y/grid_unit_,
                                             coordinate.z/grid_unit_);
    StorageType info = retrieve(triplet);
    delete triplet;
    return info;
}

template<typename T>
typename CoordinateGrid<T>::StorageType
CoordinateGrid<T>::retrieve_adjacent_cells(const Coordinate& coordinate) const {
    StorageType found = new std::vector<T>;

    Triplet<int> *triplet_base = new Triplet<int>(coordinate.x/grid_unit_ - 1,
                                                  coordinate.y/grid_unit_ - 1,
                                                  coordinate.z/grid_unit_ - 1);
    Triplet<int> *triplet = new Triplet<int>(*triplet_base);
    for (int i = 0; i <= 2; i++) {
        for (int j = 0; j <= 2; j++) {
            for (int k = 0; k <= 2; k++) {
                triplet->first = triplet_base->first + i;
                triplet->second = triplet_base->second + j;
                triplet->third = triplet_base->third + k;
                const StorageType info = retrieve(triplet);
                if (info != NULL)
                    found->insert(found->end(), info->begin(), info->end());
            }
        }
    }
    delete triplet;
    delete triplet_base;

    return found;
}

}  // namespace gmml

#endif  // GMML_INTERNAL_COORDINATE_GRID_H_
