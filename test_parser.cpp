#include "params.hpp"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    po::variables_map vm = params::parse_options(argc, argv);

    return 0;
}
