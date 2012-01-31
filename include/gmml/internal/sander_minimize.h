// Author: Robert Davis

#ifndef GMML_INTERNAL_SANDER_MINIMIZE_H_
#define GMML_INTERNAL_SANDER_MINIMIZE_H_

#include <iosfwd>
#include <string>

namespace gmml {

class Environment;
class ParameterSet;
class Structure;

class MinimizationResults;

// This functor minimizes structures with sander, given a sander input (mdin)
// file.
class SanderMinimize {
  public:
    // If sander is not part of the user's PATH or if the the minimization
    // failed, the return is NULL.
    MinimizationResults *operator()(Structure& structure,
                                    const std::string& input_file) const;
};

// This represents the results of a minimization as found in sander's mdout
// file.
class MinimizationResults {
  public:
    // Returns the results of the minimization, given a SANDER output (mdout)
    // file. If the file doesn't exist or if there was an error parsing it,
    // NULL is returned.
    static MinimizationResults *parse(const std::string& mdout_file);

    double energy() const { return energy_; }
    double bond_energy() const { return bond_energy_; }
    double vdw_energy() const { return vdw_energy_; }

  private:
    MinimizationResults(double energy, double bond_energy, double vdw_energy)
            : energy_(energy), bond_energy_(bond_energy),
              vdw_energy_(vdw_energy) {}

    static MinimizationResults *parse(std::istream& in);

    double energy_;
    double bond_energy_;
    double vdw_energy_;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_SANDER_MINIMIZE_H_
