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

#ifndef GMML_INTERNAL_ATOM_H_
#define GMML_INTERNAL_ATOM_H_

#include <string.h>

#include <string>

#include "gmml/internal/element.h"
#include "gmml/internal/geometry.h"
#include "gmml/internal/stubs/common.h"

namespace gmml {

class Atom {
  public:
    // If a type isn't given, it's set to the empty string.
    Atom(const Element& element, const Coordinate& coordinate,
         const std::string& name, double charge);

    Atom(const Element& element, const Coordinate& coordinate,
         const std::string& name, const std::string& type, double charge); 

    virtual ~Atom() {}

    virtual Atom *clone() const {
        return new Atom(element_, coordinate_, name_, type_, charge_);
    }

    void translate(double x, double y, double z) {
        coordinate_.translate(x, y, z);
    }

    //
    // Mutators
    //
    void set_coordinate(const Coordinate& coordinate) {
        coordinate_ = coordinate;
    }

    void set_coordinate(double x, double y, double z) {
        coordinate_.x = x;
        coordinate_.y = y;
        coordinate_.z = z;
    }

    void set_element(Element element) { element_ = element; }
    void set_name(const std::string& name) { name_ = name; }
    void set_charge(double charge) { charge_ = charge; }
    void set_type(const std::string& type) { type_ = type; }

    // 
    // Accessors
    //
    const Coordinate& coordinate() const { return coordinate_; }
    Coordinate& mutable_coordinate() { return coordinate_; }

    const Element& element() const { return element_; }
    const std::string& name() const { return name_; }
    const std::string& type() const { return type_; }
    double charge() const { return charge_; }

  private:
    Coordinate coordinate_;
    Element element_;
    std::string name_;
    std::string type_;
    double charge_;

    DISALLOW_COPY_AND_ASSIGN(Atom);
};

inline Atom::Atom(const Element& element, const Coordinate& coordinate,
                  const std::string& name, double charge)
        : coordinate_(coordinate), element_(element), name_(name),
          type_(""), charge_(charge) {}

inline Atom::Atom(const Element& element, const Coordinate& coordinate,
                  const std::string& name, const std::string& type,
                  double charge)
        : coordinate_(coordinate), element_(element),  name_(name),
          type_(type), charge_(charge) {}

}  // namespace gmml

#endif  // GMML_INTERNAL_ATOM_H_
