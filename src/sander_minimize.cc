// Author: Robert Davis

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "gmml/internal/sander_minimize.h"

#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <cmath>
#include <cstdio>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "gmml/internal/amber_top_file.h"
#include "gmml/internal/coordinate_file.h"
#include "gmml/internal/environment.h"
#include "gmml/internal/parameter_file.h"
#include "gmml/internal/structure.h"
#include "gmml/internal/stubs/logging.h"
#include "utilities.h"

using std::string;
using std::vector;

namespace gmml {

MinimizationResults *SanderMinimize::operator()(Structure& structure,
                                                const string& mdin_file) const {
#ifndef HAVE_SANDER
    LOG(WARNING) << "You must have sander in your path to minimize.";
    return NULL;
#endif
    const ParameterSet& parm_set = *kDefaultEnvironment.parm_set();
    File absolute_file(mdin_file);
    set_full_pathname(&absolute_file);

    struct timeval tv;
    gettimeofday(&tv, NULL);
    string uid = to_string(static_cast<int>(tv.tv_usec));

    pid_t child = fork();
    if (child < 0) {
        LOG(WARNING) << "Fork failed.";
        return NULL;
    } else if (child == 0) {
        structure.print_amber_top_file(uid + "_temp.top", parm_set);
        structure.print_coordinate_file(uid + "_temp.rst");
        const char* args[] = { 
                "sander", "-i", absolute_file.pathname().c_str(),
                "-p", string(uid + "_temp.top").c_str(),
                "-c", string(uid + "_temp.rst").c_str(),
                "-O",
                "-r", string(uid + "_out_temp.rst").c_str(),
                "-x", string(uid + "_mdcrd").c_str(),
                "-e", string(uid + "_mden").c_str(),
                "-inf", string(uid + "_mdinfo").c_str(),
                "-o", string(uid + "_mdout").c_str(),
                (char *) NULL
        };
        execvp("sander", const_cast<char* const *>(args));
        return NULL;
        LOG(WARNING) << "Error in exec.";
    }
    
    int child_exit_status;
    waitpid(-1, &child_exit_status, 0);
    if (child_exit_status > 0) {
	LOG(WARNING) << "Error in sander, exited with status " <<
		        child_exit_status << ".";
    }
    else {
	structure.load_coordinates(CoordinateFile(uid + "_out_temp.rst"));
    }
    // Clean up after ourselves.
    remove(string(uid + "_temp.top").c_str());
    remove(string(uid + "_temp.rst").c_str());

    // Clean up after sander.
    remove(string(uid + "_mdcrd").c_str());
    remove(string(uid + "_mden").c_str());
    remove(string(uid + "_mdinfo").c_str());
    remove(string(uid + "_out_temp.rst").c_str());

    MinimizationResults *results = MinimizationResults::parse(uid + "_mdout");
    remove(string(uid + "_mdout").c_str());
    return results;
}

MinimizationResults *MinimizationResults::parse(const string& mdout_file) {
    std::ifstream file(mdout_file.c_str());
    if (file.fail()) {
        return NULL;
    }
    MinimizationResults *results = parse(file);
    file.close();
    return results;
}

MinimizationResults *MinimizationResults::parse(std::istream& in) {
    string line;
    // Skip down to the results.
    while (std::getline(in, line) &&
            (line.find("FINAL RESULTS") == string::npos))
        ;
    // The line after the one with "ENERGY" has our results.
    while (std::getline(in, line) && (line.find("ENERGY") == string::npos))
        ;
    std::getline(in, line);

    std::stringstream ss(line);
    string dummy;
    double energy;
    ss >> dummy >> energy;
    if (ss.fail()) {
        return NULL;
    }

    std::getline(in, line);
    std::getline(in, line);
    ss.str(line);
    double bond_energy;
    ss >> dummy >> dummy >> bond_energy;
    if (ss.fail()) {
        return NULL;
    }

    std::getline(in, line);
    ss.str(line);
    double vdw_energy;
    ss >> dummy >> dummy >> vdw_energy;
    if (ss.fail()) {
        return NULL;
    }

    return new MinimizationResults(energy, bond_energy, vdw_energy);
}

}  // namespace gmml
