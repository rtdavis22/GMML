#ifndef GMML_INTERNAL_ELEMENT_H_
#define GMML_INTERNAL_ELEMENT_H_

#include <string>

namespace gmml {

class Element {
  public:
    // construction and destruction
    explicit Element(const std::string& symbol);
    explicit Element(int atomic_number);

    // properties
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
    
    // static methods
    static Element from_name(const std::string &name);
    
  private:
    int atomic_number_;
};

}  // namespace gmml

#endif // GMML_INTERNAL_ELEMENT_H_
