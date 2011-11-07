#ifndef TORSION_COMBINATION_BUILDER
#define TORSION_COMBINATION_BUILDER

#include <list>
#include <set>
#include <string>
#include <vector>

#include "utilities.h"

namespace gmml {

class Structure;

class TorsionCombinationBuilder {
  public:
    TorsionCombinationBuilder(const Structure& structure);

    void add_phi_value(int residue_index, double measure) {
        add_values(residue_index, measure, kNotSet, kNotSet);
    }
    void add_phi_value(const std::string& carbon_residue,
                       const std::string& oxygen_residue, double measure) {
        add_values(carbon_residue, oxygen_residue, measure, kNotSet, kNotSet);
    }
    void add_phi_value(const std::string& residue1, int carbon_number,
                       const std::string& residue2, int oxygen_number,
                       double measure) {
        add_values(residue1, carbon_number, residue2, oxygen_number, measure,
                   kNotSet, kNotSet);
    }

    void add_psi_value(int residue_index, double measure) {
        add_values(residue_index, kNotSet, measure, kNotSet);
    }
    void add_psi_value(const std::string& carbon_residue,
                      const std::string& oxygen_residue, double measure) {
        add_values(carbon_residue, oxygen_residue, kNotSet, measure, kNotSet);
    }
    void add_psi_value(const std::string& residue1, int carbon_number,
                      const std::string& residue2, int oxygen_number,
                      double measure) {
        add_values(residue1, carbon_number, residue2, oxygen_number, kNotSet,
                   measure, kNotSet);
    }

    void add_omega_value(int residue_index, double measure) {
        add_values(residue_index, kNotSet, kNotSet, measure);
    }
    void add_omega_value(const std::string& carbon_residue,
                        const std::string& oxygen_residue, double measure) {
        add_values(carbon_residue, oxygen_residue, kNotSet, kNotSet, measure);
    }
    void add_omega_value(const std::string& residue1, int carbon_number,
                        const std::string& residue2, int oxygen_number,
                        double measure) {
        add_values(residue1, carbon_number, residue2, oxygen_number, kNotSet,
                   kNotSet, measure);
    }

    void add_phi_psi_value(int residue_index, double phi, double psi) {
        add_values(residue_index, phi, psi, kNotSet);
    }
    void add_phi_psi_value(const std::string& carbon_residue,
                           const std::string& oxygen_residue, double phi,
                           double psi) {
        add_values(carbon_residue, oxygen_residue, phi, psi, kNotSet);
    }
    void add_phi_psi_value(const std::string& residue1, int carbon_number,
                           const std::string& residue2, int oxygen_number,
                           double phi, double psi) {
        add_values(residue1, carbon_number, residue2, oxygen_number, phi, psi,
                   kNotSet);
    }

    void add_phi_omega_value(int residue_index, double phi, double omega) {
        add_values(residue_index, phi, kNotSet, omega);
    }
    void add_phi_omega_value(const std::string& carbon_residue,
                             const std::string& oxygen_residue, double phi,
                             double omega) {
        add_values(carbon_residue, oxygen_residue, phi, kNotSet, omega);
    }
    void add_phi_omega_value(const std::string& residue1, int carbon_number,
                             const std::string& residue2, int oxygen_number,
                             double phi, double omega) {
        add_values(residue1, carbon_number, residue2, oxygen_number, phi, 
                   kNotSet, omega);
    }

    void add_psi_omega_value(int residue_index, double psi, double omega) {
        add_values(residue_index, kNotSet, psi, omega);
    }
    void add_psi_omega_value(const std::string& carbon_residue,
                             const std::string& oxygen_residue, double psi,
                             double omega) {
        add_values(carbon_residue, oxygen_residue, kNotSet, psi, omega);
    }
    void add_psi_omega_value(const std::string& residue1, int carbon_number,
                             const std::string& residue2, int oxygen_number,
                             double psi, double omega) {
        add_values(residue1, carbon_number, residue2, oxygen_number, kNotSet,
                   psi, omega);
    }

    void add_phi_psi_omega_value(int residue_index, double phi, double psi,
                                 double omega) {
        add_values(residue_index, phi, psi, omega);
    }
    void add_phi_psi_omega_value(const std::string& carbon_residue,
                                 const std::string& oxygen_residue, double phi,
                                 double psi, double omega) {
        add_values(carbon_residue, oxygen_residue, phi, psi, omega);
    }
    void add_phi_psi_omega_value(const std::string& residue1, int carbon_number,
                                 const std::string& residue2, int oxygen_number,
                                 double phi, double psi, double omega) {
        add_values(residue1, carbon_number, residue2, oxygen_number, phi, psi,
                   omega);
    }

    void add_likely_omega_values();

    std::list<Structure*> *build() const;
    std::vector<std::vector<std::vector<double> > > *get_build_info() const;

  private:
    class LinkageAngles {
      public:
        LinkageAngles() {
            values_.insert(std::vector<double>(3, kNotSet));
        }
        void insert(const std::vector<double>& vec);

        const std::vector<double>& at(int index) const {
            std::set<std::vector<double> >::iterator it = values_.begin();
            //yeah...
            for (int i = 0; i < index; i++)
                it++;
            return *it;
        }

        int size() const { return values_.size(); }

      private:
        std::set<std::vector<double> > values_;
    };

    void add_values(const std::string& residue1, int carbon_number,
                    const std::string& residue2, int oxygen_number,
                    double phi, double psi, double omega);
    void add_values(const std::string& residue1, const std::string& residue2, 
                    double phi, double psi, double omega);
    void add_values(int residue_index, double phi, double psi, double omega) {
        std::vector<double> vec(3);
        vec[0] = phi;
        vec[1] = psi;
        vec[2] = omega;
        warning("calling insert with residue_index " + to_string(residue_index));
        linkage_angles[residue_index].insert(vec);
    }
    int get_index(std::vector<double>&, double number) const;

    std::vector<int> *get_linkages(const std::string& carbon_residue,
                                   int carbon_number,
                                   const std::string& oxygen_residue,
                                   int oxygen_number) const;
    std::vector<int> *get_linkages(const std::string& carbon_residue,
                                   const std::string& oxygen_residue) const;

    Structure *structure;
    std::vector<LinkageAngles> linkage_angles;
};

}  // namespace gmml

#endif  // TORSION_COMBINATION_BUILDER_H
