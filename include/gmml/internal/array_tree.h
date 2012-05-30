// Copyright (c) 2012 The University of Georgia
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

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

    // Creates an empty tree.
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
    if (parent_id != -1 && parent_id >= tree_.size()) {
        throw std::invalid_argument(
                "ArrayTree::insert - invalid parent index(" +
                to_string(parent_id) + ")");
    }
    tree_.push_back(std::make_pair(data, parent_id));
    return tree_.size() - 1;
}

}  // namespace internal

using internal::ArrayTree;

}  // namespace gmml

#endif  // GMML_INTERNAL_ARRAY_TREE_H_
