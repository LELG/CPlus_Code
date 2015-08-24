#include <iostream>
#include "parser.hpp"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    po::variables_map vm = parser::parse_options(argc, argv);

    using parser::operator<<;
    std::cout << vm << std::endl;

    return 0;
}
