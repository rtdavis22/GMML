#ifndef GMML_INTERNAL_GLYCAN_DRAWER_H_
#define GMML_INTERNAL_GLYCAN_DRAWER_H_

#include <iosfwd>
#include <string>

#include "gmml/internal/array_tree.h"

namespace gmml {

struct ParsedResidue;

class GlycanDrawer {
  public:
    GlycanDrawer() : dpi_(72), show_config_labels_(true),
                     show_edge_labels_(false), show_position_labels_(true) {}

    ~GlycanDrawer() {}

    void write_file(const std::string& glycan, std::ostream& out) const;

    void print_file(const std::string& glycan) const;
    void print_file(const std::string& glycan, const std::string& file) const;

    // Mutators
    void set_dpi(int dpi) { dpi_ = dpi; }

    void show_config_labels() { show_config_labels_ = true; }
    void hide_config_labels() { show_config_labels_ = false; }

    void show_edge_labels() { show_edge_labels_ = true; }
    void hide_edge_labels() { show_edge_labels_ = false; }

    void show_position_labels() { show_position_labels_ = true; }
    void hide_position_labels() { show_position_labels_ = false; }

    // Accessors
    int dpi() const { return dpi_; }
    bool show_edge_labels() const { return show_edge_labels_; }

  private:
    void write(const ArrayTree<ParsedResidue*> *parsed_tree,
               std::ostream& out) const;

    int dpi_;
    bool show_config_labels_;
    bool show_edge_labels_;
    bool show_position_labels_;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_GLYCAN_DRAWER_H_
