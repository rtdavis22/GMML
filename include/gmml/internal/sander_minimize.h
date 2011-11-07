// Author: Robert Davis

#ifndef SANDER_MINIMIZE_H
#define SANDER_MINIMIZE_H

#include <string>

namespace gmml {

class Environment;
class ParameterFileSet;
class Structure;

// This represents the results of a minimization as found in sander's mdout
// file.
// TODO: Figure out what else should be included in this struct.
struct MinimizationResults {
    double energy;
};

class SanderMinimize {
  public:
    MinimizationResults *operator()(Structure& structure,
                                    const std::string& input_file,
                                    const ParameterFileSet& parm_set) const;
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

#endif  // SANDER_MINIMIZE_H
