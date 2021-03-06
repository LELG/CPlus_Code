#ifndef PARSER
#define PARSER

#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace parser {
    po::variables_map parse_options(int argc, char *argv[], std::ostream& out);

    enum Mode {in_vivo, cell_line};
    enum TreatmentType {single_dose, metronomic, adaptive, none};
    enum DecayType {constant, linear, exp};
    std::istream& operator>>(std::istream& in, Mode& mode);
    std::ostream& operator<<(std::ostream& out, const Mode& mode);
    std::istream& operator>>(std::istream& in, TreatmentType& treatmt_type);
    std::ostream& operator<<(std::ostream& out, const TreatmentType& treatmt_type);
    std::istream& operator>>(std::istream& in, DecayType& decay);
    std::ostream& operator<<(std::ostream& out, const DecayType& decay);
    template<class T>
    std::ostream& operator<<(std::ostream& os, const std::vector<T>& v);
    std::ostream& operator<<(std::ostream& os, const po::variables_map& vm);
}

#endif
