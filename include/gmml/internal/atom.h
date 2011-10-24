#ifndef ATOM_H
#define ATOM_H

#include <string>

#include "geometry.h"

namespace gmml {

enum Element { kElementC, kElementN, kElementH, kElementO, kElementS,
               kElementP, kElementSi, kElementCl, kElementF, kElementUnknown };

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
    Atom(Element element, const Coordinate& coordinate, const std::string& name,
         double charge)
            : coordinate_(coordinate), element_(element), name_(name),
              charge_(charge) {}
    Atom(Element element, const Coordinate& coordinate, const std::string& name,
         const std::string& type, double charge) 
            : coordinate_(coordinate), element_(element),  name_(name), 
              type_(type), charge_(charge) {}

    Coordinate& coordinate() { return coordinate_; }
    const Coordinate& coordinate() const { return coordinate_; }
    Element element() const { return element_; }    
    const std::string& name() const { return name_; }
    const std::string& type() const { return type_; }
    double charge() const { return charge_; }

    Atom *clone() { return new Atom(*this); }

  private:
    Coordinate coordinate_;
    Element element_;
    std::string name_;
    std::string type_;
    double charge_;
};

}  // namespace gmml

#endif  // ATOM_H
