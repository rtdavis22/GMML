// Author: Robert Davis

#ifndef GMML_INTERNAL_ATOM_H_
#define GMML_INTERNAL_ATOM_H_

#include <string>

#include "gmml/internal/geometry.h"
#include "gmml/internal/stubs/common.h"

namespace gmml {

// kElementX is the atomic number of the element with symbol X.
enum Element {
    kElementUnknown = 0,
    kElementH = 1,
    kElementC = 6,
    kElementN = 7,
    kElementO = 8,
    kElementF = 9,
    kElementSi = 14,
    kElementP = 15,
    kElementS = 16,
    kElementCl = 17
};

// implement this and use it in place of get_element_by_char
inline Element get_element(const std::string& symbol) {
    return kElementUnknown;
}

inline std::string get_element_symbol(Element element) {
    if (element == kElementC)
        return "C";
    else if (element == kElementH)
        return "H";
    else if (element == kElementN)
        return "N";
    else if (element == kElementO)
        return "O";
    else if (element == kElementP)
        return "P";
    else if (element == kElementS)
        return "S";
    return "";
}

//may want to use a static map instead
inline Element get_element_by_char(char letter) {
    switch (letter) {
        case 'C':
            return kElementC;
        case 'H':
            return kElementH;
        case 'N':
            return kElementN;
        case 'O':
            return kElementO;
        case 'P':
            return kElementP;
        case 'S':
            return kElementS;
        default:
            return kElementUnknown;
    }
}

class Atom {
  public:
    // If a type isn't given, it's set to the empty string.
    Atom(Element element, const Coordinate& coordinate, const std::string& name,
         double charge);

    Atom(Element element, const Coordinate& coordinate, const std::string& name,
         const std::string& type, double charge); 

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

    Element element() const { return element_; }
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

inline Atom::Atom(Element element, const Coordinate& coordinate,
                  const std::string& name, double charge)
        : coordinate_(coordinate), element_(element), name_(name),
          type_(""), charge_(charge) {}

inline Atom::Atom(Element element, const Coordinate& coordinate,
                  const std::string& name, const std::string& type,
                  double charge)
        : coordinate_(coordinate), element_(element),  name_(name),
          type_(type), charge_(charge) {}

}  // namespace gmml

#endif  // GMML_INTERNAL_ATOM_H_
