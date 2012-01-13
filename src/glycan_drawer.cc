#include "gmml/internal/glycan_drawer.h"

#include <fstream>
#include <iostream>

#include "gmml/internal/glycam_code_set.h"
#include "gmml/internal/glycam_parser.h"

using std::endl;
using std::string;

namespace gmml {
namespace {

struct NameAndRingType {
    NameAndRingType(const string& parsed_name);

    string name;
    ResidueClassification::RingType ring_type;
};

NameAndRingType::NameAndRingType(const string& parsed_name)
        : name(parsed_name), ring_type(ResidueClassification::kPyranose) {
    if (name.size() > 3) {
        char ring_letter = name[3];
        if (ring_letter == 'p' || ring_letter == 'f') {
            name.erase(name.begin() + 3);
            if (ring_letter == 'f') {
                ring_type = ResidueClassification::kFuranose;
            }
        }
    }
}

}  // namespace

void GlycanDrawer::print_file(const string& glycan) const {
    write_file(glycan, std::cout);
}

void GlycanDrawer::print_file(const string& glycan, const string& file) const {
    std::ofstream out;
    out.open(file.c_str());
    write_file(glycan, out);
    out.close();
}

void GlycanDrawer::write_file(const string& glycan, std::ostream& out) const {
    GlycamParser parser;
    parser.dont_parse_derivatives();
    ArrayTree<ParsedResidue*> *residue_tree = parser.get_array_tree(glycan);
    write(residue_tree, out);
}

void GlycanDrawer::write(const ArrayTree<ParsedResidue*> *parsed_tree,
                         std::ostream& out) const {
    out << "graph G {" << endl;
    out << "  graph [splines=false dpi=" << dpi_ << "];" << endl;
    out << "  node [label=\"\" regular=true];" << endl;
    out << "  edge [labelfontsize=10 labeldistance=1];" << endl;
    out << "  rankdir=RL" << endl;

    for (int i = 1; i < parsed_tree->size(); i++) {
        out << i << " ";
        NameAndRingType info((*parsed_tree)[i].first->name);
        string name = info.name;
        if (name == "Glc")
            out << "[shape=circle style=filled fillcolor=\"#0000fa\"]";
        else if (name == "GlcNAc")
            out << "[shape=box style=filled fillcolor=\"#0000fa\"]";
        else if (name == "Gal")
            out << "[shape=circle style=filled fillcolor=\"#ffff00\"]";
        else if (name == "GalNAc")
            out << "[shape=box style=filled fillcolor=\"#ffff00\"]";
        else if (name == "Man")
            out << "[shape=circle style=filled fillcolor=\"#00c832\"]";
        else if (name == "ManNAc")
            out << "[shape=box style=filled fillcolor=\"#00c832\"]";
        else if (name == "Neu5Ac" || name == "NeuNAc")
            out << "[shape=diamond style=filled fillcolor=\"#c800c8\"]";
        else if (name == "Neu5Gc" || name == "NeuNGc")
            out << "[shape=diamond style=filled fillcolor=\"#e9ffff\"]";
        else if (name == "Fuc")
            out << "[shape=triangle style=filled fillcolor=\"#fa0000\"]";
        else if (info.ring_type == ResidueClassification::kPyranose)
            out << "[shape=hexagon label=\"" << name << "\"]";
        else
            out << "[shape=pentagon label=\"" << name << "\"]";

        out << endl;
    }

    for (int i = 2; i < parsed_tree->size(); i++) {
        out << (*parsed_tree)[i].second << "-- " << i;

        out << " [";

        if (show_config_labels_) {
            out << "headlabel=\"";
            if ((*parsed_tree)[i].first->configuration ==
                    ResidueClassification::kAlpha) {
                out << "&alpha;";
            } else {
                out << "&beta;";
            }
            out << "\" ";
        }

        if (show_position_labels_) {
            out << "taillabel=\"" << (*parsed_tree)[i].first->oxygen_position <<
                   "\" ";
        }

        if (show_edge_labels_) {
            out << "label=<<B>" << (i - 1) << "</B>>";
        }

        out << "];" << endl; // headlabel
    }
    out << "}" << endl; // graph
}

}  // namespace gmml
