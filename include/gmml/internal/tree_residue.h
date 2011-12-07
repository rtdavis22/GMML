// Author: Robert Davis

#ifndef GMML_INTERNAL_TREE_RESIDUE_H_
#define GMML_INTERNAL_TREE_RESIDUE_H_

#include <exception>
#include <stack>
#include <string>

#include "gmml/internal/environment.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/structure.h"
#include "gmml/internal/stubs/common.h"
#include "gmml/internal/tree.h"

namespace gmml {

class Environment;
class Structure;
class TreeResidue;

class TreeResidue {
  public:
    TreeResidue(const std::string& name) : name_(name), has_parent_(false) {}
    TreeResidue(const std::string& name, const std::string& anomeric_carbon,
                const std::string& parent_oxygen)
            : name_(name), anomeric_carbon_(anomeric_carbon),
              parent_oxygen_(parent_oxygen), has_parent_(true) {}

    std::string name() const { return name_; }
    std::string anomeric_carbon() const { return anomeric_carbon_; }
    std::string parent_oxygen() const { return parent_oxygen_; }
    bool has_parent() const { return has_parent_; }

  private:
    std::string name_;
    std::string anomeric_carbon_;
    std::string parent_oxygen_;
    bool has_parent_;
};

class BuildResidueTreeException : public std::exception {
  public:
    explicit BuildResidueTreeException(const std::string& what) {
        what_ = what;
    }

    virtual const char *what() const throw() { return what_.c_str(); }

    virtual ~BuildResidueTreeException() throw() {}

  private:
    std::string what_;
};

template<typename AttachFunc>
Structure *build_residue_tree(tree<TreeResidue*> *residue_tree,
                              const Environment& environment,
                              AttachFunc attach_func) {
    Structure *structure = new Structure;
    std::stack<int> residue_stack;

    tree<TreeResidue*>::pre_order_iterator tree_it = residue_tree->begin();
    residue_stack.push(structure->append((*tree_it)->name()));
    ++tree_it;
    while (tree_it != residue_tree->end()) {
        Residue *residue = build_prep_file((*tree_it)->name());
        if (residue == NULL) {
            throw BuildResidueTreeException(
                "Unable to build residue " + (*tree_it)->name());
        }
        int residue_index = attach_func(*structure,
                                        residue,
                                        (*tree_it)->anomeric_carbon(),
                                        residue_stack.top(),
                                        (*tree_it)->parent_oxygen());
        delete residue;
        residue_stack.push(residue_index);
        int cur_depth = residue_tree->depth(tree_it);
        ++tree_it;
        int new_depth = residue_tree->depth(tree_it);
        for (int i = 0; i < cur_depth - new_depth + 1; i++)
            residue_stack.pop();
    }

    return structure;
}

template<typename AttachFunc>
Structure *build_residue_tree(const ArrayTree<TreeResidue*> *tree,
                              const Environment& environment,
                              AttachFunc attach_func) {
    Structure *structure = new Structure;

    ArrayTree<TreeResidue*>::const_iterator it = tree->begin();
    structure->append(it->first->name());
    ++it;

    while (it != tree->end()) {
        Residue *residue = build_prep_file(it->first->name());
        if (residue == NULL)
            throw BuildResidueTreeException(
                "Unable to build residue with code " + it->first->name());
        attach_func(*structure, residue, it->first->anomeric_carbon(),
                    it->second, it->first->parent_oxygen());
        delete residue;
        ++it;
    }

    return structure;
}

template<typename TREE, typename AttachFunc>
inline Structure *build_residue_tree(TREE *residue_tree,
                                     AttachFunc attach_func) {
    return build_residue_tree(residue_tree, kDefaultEnvironment, attach_func);
}

template<typename TREE, typename AttachFunc>
inline Structure *build_residue_tree(tree<TreeResidue*> *residue_tree,
                                     const Environment& environment) {
    return build_residue_tree(residue_tree, environment, AttachFunc());
}

template<typename TREE, typename AttachFunc>
inline Structure *build_residue_tree(TREE *residue_tree) {
    return build_residue_tree<TREE, StructureAttach>(residue_tree,
                                                     kDefaultEnvironment);
}

}  // namespace gmml

#endif  // TREE_RESIDUE_H
