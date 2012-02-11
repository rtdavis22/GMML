// Author: Robert Davis

#ifndef GMML_INTERNAL_NETOGLYC_H_
#define GMML_INTERNAL_NETOGLYC_H_

#include <memory>
#include <string>
#include <vector>

namespace gmml {

class NetOGlycResults;
class OGlycosylationLocations;

class NetOGlycRunner {
  public:
    static void set_startup_script(const std::string& file);

    NetOGlycResults *operator()(const std::vector<std::string>& sequences);

  private:
    static std::string startup_script_;
};

class NetOGlycResults {
  public:
    NetOGlycResults(const std::string& netoglyc_file);

    virtual ~NetOGlycResults();

    const OGlycosylationLocations *get_sequence_locations(int index) const;

    int sequence_count() const;

  private:
    struct Impl;
    std::auto_ptr<Impl> impl_;
};

class OGlycosylationLocations {
  public:
    OGlycosylationLocations(const std::string& location_output) { 
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
};

}  // namespace gmml

#endif  // GMML_INTERNAL_NETOGLYC_H_
