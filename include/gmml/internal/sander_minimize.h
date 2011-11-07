#ifndef SANDER_MINIMIZE_H
#define SANDER_MINIMIZE_H

#include <string>

namespace gmml {

class Environment;
class ParameterFileSet;
class Structure;

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
    MinimizationResults *parse_output_file(const std::string& out_file) const;
};

}  // namespace gml

#endif  // SANDER_MINIMIZE_H
