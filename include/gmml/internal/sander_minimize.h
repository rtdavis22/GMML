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
    // Minimize the structure with the given SANDER input (mdin) file.
    // If sander is not part of the user's PATH or if the the minimization
    // failed, the return is NULL.
    MinimizationResults *operator()(Structure& structure,
                                    const std::string& mdin_file) const;
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
