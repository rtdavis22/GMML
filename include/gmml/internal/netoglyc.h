// Author: Robert Davis

#ifndef GMML_INTERNAL_NETOGLYC_H_
#define GMML_INTERNAL_NETOGLYC_H_

#include <memory>
#include <string>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

class FastaSequence;
class Structure;

class NetOGlycResults;
class OGlycosylationLocations;

// This is an interface with NetOGlyc, a program for predicting
// O-glycosylation locations.
class NetOGlycRunner {
  public:
    explicit NetOGlycRunner(const std::string& startup_script);

    ~NetOGlycRunner();

    void add_sequence(const FastaSequence& fasta_sequence);

    NetOGlycResults *run();

  private:
    struct Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(NetOGlycRunner);
};

class NetOGlycResults {
  public:
    explicit NetOGlycResults(const std::string& netoglyc_output_file);

    virtual ~NetOGlycResults();

    // Returns the predicted locations for the FASTA sequence at the given index
    // in the file.
    const OGlycosylationLocations *get_predicted_locations(int index) const;

    // Returns the number of distinct sequences in the input file.
    int sequence_count() const;

  private:
    struct Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(NetOGlycResults);
};

class OGlycosylationLocations {
  public:
    explicit OGlycosylationLocations(const std::string& location_output) { 
        init(location_output);
    }

    const std::vector<int>& get_serine_locations() const {
        return serine_locations_;
    }

    const std::vector<int>& get_threonine_locations() const {
        return threonine_locations_;
    }

  private:
    void init(const std::string& location_output);

    std::vector<int> serine_locations_;
    std::vector<int> threonine_locations_;

    DISALLOW_COPY_AND_ASSIGN(OGlycosylationLocations);
};

}  // namespace gmml

#endif  // GMML_INTERNAL_NETOGLYC_H_
