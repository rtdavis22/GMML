#ifndef RESIDUE_H
#define RESIDUE_H

#include <algorithm>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "atom.h"
#include "graph.h"
#include "utilities.h"

namespace gmml {

class Residue {
  public:
    typedef boost::shared_ptr<Atom> AtomPtr;
    // try to get rid of default ctor, i think
    Residue() {}
    template<typename InputIterator>
    Residue(Graph *bonds, const std::string& name,
            InputIterator begin, InputIterator end)
            : bonds_(bonds), name_(name) {
        atoms_->assign(begin, end);
    }
    Residue(Graph *bonds, const std::string& name,
            std::vector<AtomPtr> *atoms)
            : bonds_(bonds), name_(name), atoms_(atoms) {}

    virtual ~Residue() {
        //delete bonds_;
        //std::for_each(begin(), end(), DeletePtr());
    }

    typedef std::vector<AtomPtr>::iterator iterator;
    typedef std::vector<AtomPtr>::const_iterator const_iterator;

    iterator begin() { return atoms_->begin(); }
    const_iterator begin() const { return atoms_->begin(); }

    iterator end() { return atoms_->end(); }
    const_iterator end() const { return atoms_->end(); }

    size_t size() const { return atoms_->size(); }

    const Graph *bonds() const { return bonds_; }
    std::string name() const { return name_; }

  private:
    Graph *bonds_;
    std::string name_;
    std::vector<AtomPtr> *atoms_;

    DISALLOW_COPY_AND_ASSIGN(Residue);
};

} //namespace gmml

#endif
