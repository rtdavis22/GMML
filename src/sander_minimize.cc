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
#include "utilities.h"

using std::string;
using std::vector;

namespace gmml {

MinimizationResults *SanderMinimize::operator()(
        Structure& structure, const string& input_file,
        const ParameterSet& parm_set) const {
#ifndef HAVE_SANDER
    warning("SanderMinimize - You must have sander in your path to minimize.");
    return NULL;
#endif
    string absolute_file = find_file(input_file);

    struct timeval tv;
    gettimeofday(&tv, NULL);
    string uid = to_string(static_cast<int>(tv.tv_usec));

    pid_t child = fork();
    if (child < 0) {
        warning("SanderMinimize - Fork failed.");
        return NULL;
    } else if (child == 0) {
        structure.print_amber_top_file(uid + "_temp.top", parm_set);
        structure.print_coordinate_file(uid + "_temp.rst");
        // TODO: Why isn't this const? Fix this.
        char *args[] = { "sander", "-i",
                         const_cast<char*>(absolute_file.c_str()),
                         "-p",
                         const_cast<char*>(string(uid + "_temp.top").c_str()),
                         "-c",
                         const_cast<char*>(string(uid + "_temp.rst").c_str()),
                         "-O",
                         "-r",
                         const_cast<char*>(
                             string(uid + "_out_temp.rst").c_str()
                         ),
                         "-x",
                         const_cast<char*>(string(uid + "_mdcrd").c_str()),
                         "-e",
                         const_cast<char*>(string(uid + "_mden").c_str()),
                         "-inf",
                         const_cast<char*>(string(uid + "_mdinfo").c_str()),
                         "-o",
                         const_cast<char*>(string(uid + "_mdout").c_str()),
                          (char *) 0 };
        execvp("sander", args);
        return NULL;
        warning("SanderMinimize - Error in exec.");
    }
    
    int child_exit_status;
    waitpid(-1, &child_exit_status, 0);
    if (child_exit_status > 0) {
	warning("SanderMinimize - Error in sander, exited with status " +
		to_string(child_exit_status) + ".");
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

MinimizationResults *SanderMinimize::operator()(
        Structure& structure, const string& input_file,
        const Environment& environment) const {
    return operator()(structure, input_file, *environment.parm_set());
}

MinimizationResults *SanderMinimize::operator()(
        Structure& structure, const string& input_file) const {
    return operator()(structure, input_file, kDefaultEnvironment);
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
