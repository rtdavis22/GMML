// Author: Robert Davis

#include "gmml/internal/netoglyc.h"

#include <cmath>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <fstream>
#include <sstream>

#include "gmml/internal/environment.h"
#include "gmml/internal/stubs/common.h"
#include "utilities.h"

using std::endl;
using std::ifstream;
using std::istringstream;
using std::ofstream;
using std::string;
using std::vector;

namespace gmml {

string NetOGlycRunner::startup_script_("");

void NetOGlycRunner::set_startup_script(const string& file) {
    startup_script_ = find_file(file);
}

NetOGlycResults *NetOGlycRunner::operator()(const vector<string>& sequences) {
    ofstream out(".netoglyc_input");
    for (int i = 0; i < sequences.size(); i++) {
        out << ">" << endl;
        out << sequences[i] << endl;
    }
    out.close();

    pid_t child = fork();
    if (child < 0) {
        // ERROR
    } else if (child == 0) {
        int fd = open(".netoglyc_output", O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
        dup2(fd, 1);
        close(fd);
        const char* const args[] = { startup_script_.c_str(), ".netoglyc_input",
                                     (char *) NULL };
        execvp(startup_script_.c_str(), const_cast<char* const *>(args));
    }

    int exit_status;
    waitpid(-1, &exit_status, 0);
    if (exit_status > 0) {
        // ERROR
    }

    remove(".netoglyc_input");

    NetOGlycResults *results = new NetOGlycResults(".netoglyc_output");
    remove(".netoglyc_output");
    return results;
}


// Private Implementation
struct NetOGlycResults::Impl {
    ~Impl() {
        STLDeleteContainerPointers(predicted_locations.begin(),
                                   predicted_locations.end());
    }

    void init(const std::string& netoglyc_file);
    void init_from_stream(std::istream& in);
    bool advance_to_next_sequence_results(std::istream& in);
    string read_result_sequence(std::istream& in);

    vector<OGlycosylationLocations*> predicted_locations;
};

void NetOGlycResults::Impl::init(const string& netoglyc_file) {
    ifstream in(netoglyc_file.c_str());
    if (in.fail()) {
        in.close();
        throw FileNotFoundException(netoglyc_file);
    }

    init_from_stream(in);
    in.close();
}

void NetOGlycResults::Impl::init_from_stream(std::istream& in) {
    bool done = !advance_to_next_sequence_results(in);
    while (!done) {
        string sequence_result = read_result_sequence(in);
        predicted_locations.push_back(
                new OGlycosylationLocations(sequence_result));
        done = !advance_to_next_sequence_results(in);
    }
}

bool NetOGlycResults::Impl::advance_to_next_sequence_results(std::istream& in) {
    string line;
    int sequence_length = -1;
    while (getline(in, line)) {
        size_t position = line.find("Length:");
        if (position != string::npos) {
            istringstream ss(line.substr(position));
            string dummy;
            ss >> dummy >> sequence_length;
            break;
        }
    }
    if (sequence_length == -1) {
        return false;
    }
    int lines_to_skip = ceil(sequence_length/80.0);
    for (int i = 0; i < lines_to_skip; i++)
        getline(in, line);

    return in.good();
}

string NetOGlycResults::Impl::read_result_sequence(std::istream& in) {
    string result_sequence = "";
    string line;
    while (getline(in, line)) {
        trim(line);
        if (line == "") {
            break;
        } else {
            result_sequence += line;
        }
    }
    return result_sequence;
}

NetOGlycResults::NetOGlycResults(const string& netoglyc_file)
        : impl_(new Impl) {
    impl_->init(netoglyc_file);
}

NetOGlycResults::~NetOGlycResults() {
}

const OGlycosylationLocations *NetOGlycResults::get_sequence_locations(
        int index) const {
    return impl_->predicted_locations[index];
}

int NetOGlycResults::sequence_count() const {
    return impl_->predicted_locations.size();
}

void OGlycosylationLocations::init(const string& location_output) {
    int length = location_output.size();
    for (int i = 0; i < length; i++) {
        char letter = location_output[i];
        if (letter == 'S') {
            serine_locations_.push_back(i);
        } else if (letter == 'T') {
            threonine_locations_.push_back(i);
        }
    }
}


}  // namespace gmml
