#include "utilities.h"

#include <iostream>

namespace gmml {

bool kEnableWarnings = false;

// These should be removed in favor of logging.
void warning(const std::string& message) {
    if (kEnableWarnings)
        std::cerr << "warning: " << message << std::endl;
}

void error(const std::string& message) {
    std::cerr << "error: " << message << std::endl;
    throw ExitException(1);
}

}  // namespace gmml
