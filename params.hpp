#ifndef PARAMS
#define PARAMS

#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace params {
    po::variables_map parse_options(int argc, char *argv[]);

    enum Mode {in_vivo, cell_line};
    std::istream& operator>>(std::istream& in, Mode& mode);
}

#endif
