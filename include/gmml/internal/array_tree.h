// Author: Robert Davis
//
// This is a bare-bones implementation of a tree that uses an array as the
// underlying data structure. It only includes functions that are needed by
// the rest of the library.

#ifndef GMML_INTERNAL_ARRAY_TREE_H_
#define GMML_INTERNAL_ARRAY_TREE_H_

#include <stdexcept>
#include <utility>
#include <vector>

#include "gmml/internal/stubs/utils.h"

namespace gmml {
namespace internal {

template<typename T>
class ArrayTree {
  public:
    typedef typename std::vector<std::pair<T, int> >::iterator iterator;
    typedef typename std::vector<std::pair<T, int> >::const_iterator
            const_iterator;

    ArrayTree() {}

    const_iterator begin() const { return tree_.begin(); }
    const_iterator end() const { return tree_.end(); }

    size_t size() const { return tree_.size(); }

    const std::pair<T, int>& operator[](int i) const { return tree_[i]; }

    int insert(T data) { return insert(data, -1); }
    int insert(T data, int parent_id);

  private:
    typename std::vector<std::pair<T, int> > tree_;
};

template<typename T>
inline int ArrayTree<T>::insert(T data, int parent_id) {
    if (parent_id != -1 && parent_id >= tree_.size())
        throw std::invalid_argument(
                "ArrayTree::insert - invalid parent index(" +
                to_string(parent_id) + ")");
    tree_.push_back(std::make_pair(data, parent_id));
    return tree_.size() - 1;
}

}  // namespace internal

using internal::ArrayTree;

}  // namespace gmml

#endif  // GMML_INTERNAL_ARRAY_TREE_H_
