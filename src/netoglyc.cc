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
#include "gmml/internal/fasta_sequence.h"
#include "gmml/internal/proteins.h"
#include "gmml/internal/structure.h"
#include "gmml/internal/stubs/file.h"
#include "gmml/internal/stubs/logging.h"
#include "utilities.h"

using std::endl;
using std::ifstream;
using std::istringstream;
using std::ofstream;
using std::string;
using std::vector;

namespace gmml {

struct NetOGlycRunner::Impl {
    static const char *kInputFileName;
    static const char *kOutputFileName;

    explicit Impl(const File& startup_script)
            : startup_script(startup_script) {
        if (!set_full_pathname(&this->startup_script)) {
            throw FileNotFoundException(this->startup_script.pathname());
        }
    }

    ~Impl() {
        STLDeleteContainerPointers(fasta_sequences.begin(),
                                   fasta_sequences.end());
    }

    void write_fasta_file();
    void exec_netoglyc();
    void cleanup_files();

    File startup_script;
    vector<FastaSequence*> fasta_sequences;
};

const char *NetOGlycRunner::Impl::kInputFileName = ".netoglyc_input";
const char *NetOGlycRunner::Impl::kOutputFileName = ".netoglyc_output";

void NetOGlycRunner::Impl::write_fasta_file() {
    ofstream out(kInputFileName);
    for (int i = 0; i < fasta_sequences.size(); i++) {
        out << ">" << endl;
        out << fasta_sequences[i]->sequence() << endl;
    }
    out.close();
}

void NetOGlycRunner::Impl::exec_netoglyc() {
    int fd = open(kOutputFileName, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
    dup2(fd, 1);
    close(fd);
    const char* const args[] = { startup_script.pathname().c_str(),
                                 kInputFileName,
                                 (char *) NULL };
    execvp(startup_script.pathname().c_str(), const_cast<char* const *>(args));
}

void NetOGlycRunner::Impl::cleanup_files() {
    remove(kInputFileName);
    remove(kOutputFileName);
}

NetOGlycRunner::NetOGlycRunner(const File& startup_script)
        : impl_(new Impl(startup_script)) {
}

NetOGlycRunner::~NetOGlycRunner() {
}

void NetOGlycRunner::add_sequence(const FastaSequence& fasta_sequence) {
    impl_->fasta_sequences.push_back(fasta_sequence.clone());
}

NetOGlycResults *NetOGlycRunner::run() {
    impl_->write_fasta_file();
    pid_t child = fork();
    if (child < 0) {
        LOG(ERROR) << "fork() failed.";
        impl_->cleanup_files();
        return NULL;
    } else if (child == 0) {
        impl_->exec_netoglyc();
    }

    int exit_status;
    waitpid(-1, &exit_status, 0);

    NetOGlycResults *results = NULL;
    if (exit_status > 0) {
        LOG(ERROR) << "Error running NetOGlyc.";
    } else {
        results = new NetOGlycResults(impl_->kOutputFileName);
    }

    impl_->cleanup_files();
    return results;
}


struct NetOGlycResults::Impl {
    ~Impl() {
        STLDeleteContainerPointers(predicted_locations.begin(),
                                   predicted_locations.end());
    }

    static const int kOutputWidth = 80;

    void init(const std::string& netoglyc_output_file);
    void init_from_stream(std::istream& in);
    bool advance_to_next_sequence_results(std::istream& in);
    string read_result_sequence(std::istream& in);

    vector<OGlycosylationLocations*> predicted_locations;
};

void NetOGlycResults::Impl::init(const string& netoglyc_output_file) {
    ifstream in(netoglyc_output_file.c_str());
    if (in.fail()) {
        in.close();
        throw FileNotFoundException(netoglyc_output_file);
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

    int lines_to_skip = ceil(static_cast<double>(sequence_length)/kOutputWidth);
    for (int i = 0; i < lines_to_skip; i++) {
        getline(in, line);
    }

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

NetOGlycResults::NetOGlycResults(const string& netoglyc_output_file)
        : impl_(new Impl) {
    impl_->init(netoglyc_output_file);
}

NetOGlycResults::~NetOGlycResults() {
}

const OGlycosylationLocations *NetOGlycResults::get_predicted_locations(
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
