// Author: Robert Davis

#ifndef GMML_INTERNAL_SANDER_MINIMIZE_H_
#define GMML_INTERNAL_SANDER_MINIMIZE_H_

#include <string>

namespace gmml {

class Environment;
class ParameterSet;
class Structure;

// This represents the results of a minimization as found in sander's mdout
// file.
// TODO: Figure out what else should be included in this struct.
struct MinimizationResults {
    double energy;
};

// This functor minimizes structures with sander, given a sander input (mdin)
// file.
class SanderMinimize {
  public:
    // If sander is not part of the user's PATH or if the the minimization
    // failed, the return is NULL.
    MinimizationResults *operator()(Structure& structure,
                                    const std::string& input_file,
                                    const ParameterSet& parm_set) const;
    MinimizationResults *operator()(Structure& structure,
                                    const std::string& input_file,
                                    const Environment& environment) const;
    MinimizationResults *operator()(Structure& structure,
                                    const std::string& input_file) const;
  private:
    // This might should go in MinimizationResults.
    MinimizationResults *parse_output_file(const std::string& out_file) const;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_SANDER_MINIMIZE_H_
