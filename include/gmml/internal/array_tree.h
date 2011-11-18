#ifndef ARRAY_TREE_H
#define ARRAY_TREE_H

#include <stdexcept>
#include <utility>
#include <vector>

namespace gmml {

template<typename T>
class ArrayTree {
  public:
    typedef typename std::vector<std::pair<T, int> >::iterator iterator;
    typedef typename std::vector<std::pair<T, int> >::const_iterator
            const_iterator;

    const_iterator begin() const { return tree.begin(); }
    const_iterator end() const { return tree.end(); }

    int insert(T data) { return insert(data, -1); }
    int insert(T data, int parent_id);

  private:
    typename std::vector<std::pair<T, int> > tree;
};

template<typename T>
inline int ArrayTree<T>::insert(T data, int parent_id) {
    if (parent_id >= tree.size())
        throw std::invalid_argument("ArrayTree::insert - invalid parent index");
    tree.push_back(std::make_pair(data, parent_id));
    return tree.size() - 1;
}

}  // namespace gmml

#endif  // ARRAY_TREE_H
