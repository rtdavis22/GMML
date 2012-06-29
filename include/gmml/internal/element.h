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

// Author: Kyle Forrester

#ifndef GMML_INTERNAL_ELEMENT_H_
#define GMML_INTERNAL_ELEMENT_H_

#include <string>

namespace gmml {

/**
 This class represents a chemical element.
 */
class Element {
  public:
    /**
     Creates an unknown element.
     */
    Element() : atomic_number_(0) {}

    /**
     Creates an element with the given atomic number.
     */
    explicit Element(int atomic_number);

    /**
     Creates an element with the given symbol.
     */
    explicit Element(const std::string& symbol);

    /**
     Returns an element with the given full name.
     */
    static Element from_name(const std::string& name);

    bool operator==(const Element& rhs) const {
        return atomic_number_ == rhs.atomic_number_;
    }

    bool operator!=(const Element& rhs) const { return !operator==(rhs); }

    bool is_unknown() const { return atomic_number_ == 0; }

    int atomic_number() const { return atomic_number_; }
    std::string symbol() const;
    std::string name() const;
    int period() const;
    double mass() const;
    double exact_mass() const;
    double ionization_energy() const;
    double electron_affinity() const;
    double electron_negativity() const;
    double covalent_radius() const;
    double van_der_waals_radius() const;
    double boiling_point() const;
    double melting_point() const;
    
  private:
    int atomic_number_;
};

inline Element guess_element_by_letter(char letter) {
    switch (letter) {
        case 'C':
            return Element("C");
        case 'H':
            return Element("H");
        case 'N':
            return Element("N");
        case 'O':
            return Element("O");
        case 'P':
            return Element("P");
        case 'S':
            return Element("S");
        default:
            return Element();
    }
}

inline Element guess_element_by_name(const std::string& name) {
    for (int i = 0; i < name.size(); i++) {
        Element guess = guess_element_by_letter(name[i]);
        if (guess != Element())
            return guess;
    }
    return Element();
}

}  // namespace gmml

#endif // GMML_INTERNAL_ELEMENT_H_
