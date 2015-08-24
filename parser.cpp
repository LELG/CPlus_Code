#include "parser.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <boost/program_options.hpp>

namespace po = boost::program_options;


/*
 * Auxiliary function used to check that 'opt1' and 'opt2' are not specified
 * at the same time.
 * Written by Vladimir Prus, found at
 * http://www.boost.org/doc/libs/1_57_0/libs/program_options/example/real.cpp
 */
void conflicting_options(const po::variables_map& vm, 
                         const char* opt1, const char* opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted() 
        && vm.count(opt2) && !vm[opt2].defaulted())
        throw std::logic_error(std::string("Conflicting options '") 
                               + opt1 + "' and '" + opt2 + "'.");
}

/*
 * Auxiliary function to check that a multitoken option has the correct
 * number of tokens specified. This function operates on vector<string>.
 */
void validate_num_tokens( const po::variables_map& vm, const char* opt, unsigned num_tokens )
{
    std::vector<std::string> tokens;
    if (vm[opt].empty() || (tokens = vm[opt].as<std::vector<std::string> >()).size() != num_tokens) {
        throw std::logic_error(std::string("Option '") + opt + "' takes exactly " + std::to_string(num_tokens) + " arguments.");
    }
}

/*
 * Aux function to check that a probability value is reasonable
 * (i.e. is between 0 and 1)
 */
void validate_probability( const po::variables_map& vm, const char* opt)
{
    double val = vm[opt].as<double>();
    if ((val < 0.0) || (1.0 < val))
    {
        throw std::logic_error(std::string("Option '") + opt + "' is a probability; it must be a float between 0 and 1.");
    }
}

namespace parser {

/*
 * Function to parse our command line options.
 */
po::variables_map parse_options(int argc, char *argv[], std::ostream& out)
{
    std::string config_fpath;

    po::options_description general{"General Options"};
    general.add_options()
        ("help,h", "= print this help message and exit")
        ("config", po::value<std::string>(&config_fpath)->value_name("<path>"),
            "= path to optional config file")
        ("test_group", po::value<std::string>()->value_name("<name>"),
             "= name for this group of tests")
        ("param_set", po::value<int>()->value_name("<int>"),
             "= number of this param set in test group")
        ("run_number", po::value<int>()->value_name("<int>"),
             "= number of this replicate run")
        ("test_group_dir", po::value<std::string>()->value_name("<dir>"),
             "= directory for this test group")
        ("param_set_dir", po::value<std::string>()->value_name("<dir>"),
             "= directory for this particular param set")
        ("run_dir", po::value<std::string>()->value_name("<dir>")->default_value("temp"),
             "= directory for this replicate run")
        ("r_output", po::bool_switch()->default_value(false),
             "= produce output in R")
        ("auto_treatment", po::bool_switch()->default_value(false),
             "= automatically introduce treatment at size limit")
        ("no_plots", po::bool_switch()->default_value(false),
             "= do not generate plots")
        ("prune_clones", po::bool_switch()->default_value(false),
             "= remove dead clones from population")
        ("prune_muts", po::bool_switch()->default_value(false),
             "= remove dead mutations from population")
        ("mode", po::value<Mode>()->value_name("{in_vivo, cell_line}")->default_value(in_vivo, "in_vivo"),
            "= simulation mode")
    ;

    po::options_description bounds{"Simulation Bounds"};
    bounds.add_options()
        ("max_cycles", po::value<int>()->default_value(1e5, "1e5")->value_name("<int>"),
            "= maximum cycles in simulation")
        ("max_size_lim", po::value<int>()->default_value(1e9, "1e9")->value_name("<int>"),
            "= population size limit, aka carrying capacity")
        ("detectable_size", po::value<int>()->default_value(1e9, "1e9")->value_name("<int>"),
            "= population size at which tumour is detectable (triggers treatment)")
        ("mid_doublings", po::value<int>()->default_value(30)->value_name("<int>"),
            "= in cell_line mode: number of doublings at halfway point")
        ("max_doublings", po::value<int>()->default_value(60)->value_name("<int>"),
            "= in cell_line mode: number of doublings at simulation end")
        ("doublings_per_flask", po::value<int>()->default_value(3)->value_name("<int>"),
            "= cell_line mode: passage flask after this many doublings")
    ;

    po::options_description probabilities{ "Probabilities & Quantiles" };
    probabilities.add_options()
        ("prob_mut_pos", po::value<double>()->default_value(0.01, "0.01")->value_name("<float>"),
            "= prob of beneficial mutation")
        ("prob_mut_neg", po::value<double>()->default_value(0.99, "0.99")->value_name("<float>"),
            "= prob of deleterious mutation")
        ("prob_inc_mut", po::value<double>()->default_value(0.0)->value_name("<float>"),
            "= prob of increasing mutation rate")
        ("prob_dec_mut", po::value<double>()->default_value(0.0)->value_name("<float>"),
            "= prob of decreasing mutation rate")
        ("driver_quantile", po::value<double>()->value_name("<float>"),
            "= cutoff quantile for driver mutations")
        ("killer_quantile", po::value<double>()->value_name("<float>"),
            "= cutoff quantile for killer mutations")
        ("beneficial_quantile", po::value<double>()->value_name("<float>"),
            "= cutoff quantile for beneficial mutations")
        ("deleterious_quantile", po::value<double>()->value_name("<float>"),
            "= cutoff quantile for deleterious mutations")
    ;

    po::options_description tumour_params{ "Tumour characteristics" };
    tumour_params.add_options()
        ("pro", po::value<double>()->default_value(0.04, "0.04")->value_name("<float>"),
            "= initial proliferation rate")
        ("die", po::value<double>()->default_value(0.03, "0.03")->value_name("<float>"),
            "= initial death rate")
        ("mut", po::value<double>()->default_value(0.001, "0.001")->value_name("<float>"),
            "= initial mutation rate")
        ("quiesc", po::value<double>()->default_value(0.005, "0.005")->value_name("<float>"),
            "= rate of quiescence")
        ("init_size", po::value<int>()->value_name("<int>"),
            "= size of initial clone")
        ("init_diversity", po::value<std::vector<std::string> >()->multitoken()->value_name("<clone file> <mut file>"),
            "= specify a heterogeneous initial population, stored in clone and mutation files")
    ;

    po::options_description treatmt_params{ "Treatment parameters" };
    treatmt_params.add_options()
        ("treatment_type", po::value<TreatmentType>()
                           ->value_name("{single_dose, metronomic, adaptive, none}")
                           ->default_value(single_dose, "single_dose"),
            "= treatment type")
        ("decay_type", po::value<DecayType>()
                       ->value_name("{constant, linear, exp}")
                       ->default_value(constant, "constant"),
            "= treatment decay type")
        ("decay_rate", po::value<double>()->default_value(0.0)->value_name("<float>"),
            "= treatment decay rate")
        ("treatment_freq", po::value<int>()->default_value(100)->value_name("<int>"),
            "= treatment frequency")
        ("adaptive_increment", po::value<double>()->default_value(0.001, "0.001")->value_name("<float>"),
            "= adaptive treatment: size of incremental dosage changes")
        ("adaptive_threshold", po::value<double>()->default_value(0.025, "0.025")->value_name("<float>"),
            "= adaptive treatment: change in tumour size that triggers dosage change")
        ("select_time", po::value<int>()->default_value(4e5, "4e5")->value_name("<int>"),
            "= time to manually introduce selective pressure")
        ("select_pressure", po::value<double>()->default_value(0.01, "0.01")->value_name("<float>"),
            "= initial value of selective pressure")
        ("mutagenic_pressure", po::value<double>()->default_value(0.0)->value_name("<float>"),
            "= initial value of mutagenic pressure")
    ;

    po::options_description resist_params{ "Resistance parameters" };
    resist_params.add_options()
        ("resistance", po::bool_switch()->default_value(false),
            "= generate pre-existing resistance mutations")
        ("num_resist_mutns", po::value<int>()->value_name("<int>"),
            "= determine how many resistance mutations to generate")
        ("resist_strength", po::value<double>()->default_value(1.0, "1.0")->value_name("<float>"),
            "= strength of resistance mutations (between 0 and 1)")
    ;

    po::options_description saving_loading{ "Saving/Loading" };
    saving_loading.add_options()
        ("save_snapshot", po::bool_switch()->default_value(false),
            "= save a population snapshot when population reaches size limit")
        ("load_snapshot", po::value<std::string>()->value_name("<path to archive>"),
            "= load a snapshot from the specified archive")
    ;

    po::options_description all("SIMULATION OPTIONS");
    all.add(general).add(bounds).add(probabilities).add(tumour_params)
       .add(treatmt_params).add(resist_params).add(saving_loading);

    po::variables_map vm;

    try
    {
        // store options to variables map
        po::store(po::command_line_parser(argc, argv).options(all).run(), vm);

        // Display help text when requested
        if (vm.count("help"))
        {
            out << all << std::endl;
            std::exit(EXIT_SUCCESS);
        }

        // cmd line notification - may raise some errors
        po::notify(vm);

        // check for config file
        if (vm.count("config"))
        {
            std::ifstream ifs(config_fpath.c_str());
            if (!ifs) {
                throw std::runtime_error("Couldn't open config file: " + config_fpath);
            } else {
                po::store(po::parse_config_file(ifs, all), vm);
                po::notify(vm);
            }
        }

        // validate options
        conflicting_options(vm, "init_size", "init_diversity");
        conflicting_options(vm, "save_snapshot", "load_snapshot");
        if (vm.count("init_diversity"))
        {
            validate_num_tokens(vm, "init_diversity", 2);
        }
        validate_probability(vm, "prob_mut_pos");
        validate_probability(vm, "prob_mut_neg");
        validate_probability(vm, "prob_inc_mut");
        validate_probability(vm, "prob_dec_mut");
    }
    catch (std::exception &e)
    {
        out << "\nError: " << e.what() << std::endl;
        out << "See --help for usage details" << std::endl;
        throw;
    }

    return vm;
}

std::istream& operator>>(std::istream& in, Mode& mode)
{
    std::string token;
    in >> token;
    if (token == "in_vivo") {
        mode = in_vivo;
    } else if (token == "cell_line") {
        mode = cell_line;
    } else {
        throw po::invalid_option_value(token);
    }
    return in;
}

std::ostream& operator<<(std::ostream& out, const Mode& mode)
{
    if (mode == in_vivo) {
        out << "in_vivo";
    } else if (mode == cell_line) {
        out << "cell_line";
    }
    return out;
}

std::istream& operator>>(std::istream& in, TreatmentType& treatmt_type)
{
    std::string token;
    in >> token;
    if (token == "single_dose") {
        treatmt_type = single_dose;
    } else if (token == "metronomic") {
        treatmt_type = metronomic;
    } else if (token == "adaptive") {
        treatmt_type = adaptive;
    } else if (token == "none") {
        treatmt_type = none;
    } else {
        throw po::invalid_option_value(token);
    }
    return in;
}

std::ostream& operator<<(std::ostream& out, const TreatmentType& treatmt_type)
{
    if (treatmt_type == single_dose) {
        out << "single_dose";
    } else if (treatmt_type == metronomic) {
        out << "metronomic";
    } else if (treatmt_type == adaptive) {
        out << "adaptive";
    } else if (treatmt_type == none) {
        out << "none";
    }
    return out;
}

std::istream& operator>>(std::istream& in, DecayType& decay)
{
    std::string token;
    in >> token;
    if (token == "constant") {
        decay = constant;
    } else if (token == "linear") {
        decay = linear;
    } else if (token == "exp") {
        decay = exp;
    } else {
        throw po::invalid_option_value(token);
    }
    return in;
}

std::ostream& operator<<(std::ostream& out, const DecayType& decay)
{
    if (decay == constant) {
        out << "constant";
    } else if (decay == linear) {
        out << "linear";
    } else if (decay == exp) {
        out << "exp";
    }
    return out;
}

// A helper function for printing out a vector
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}

std::ostream& operator<<(std::ostream& out, const po::variables_map& vm)
{
    for (const auto& it : vm) {
        out << it.first.c_str() << " = ";
        auto& value = it.second.value();
        if (auto v = boost::any_cast<int>(&value))
            out << *v;
        else if (auto v = boost::any_cast<double>(&value))
            out << *v;
        else if (auto v = boost::any_cast<std::string>(&value))
            out << *v;
        else if (auto v = boost::any_cast<bool>(&value))
            out << std::boolalpha << *v;
        else if (auto v = boost::any_cast<Mode>(&value))
            out << *v;
        else if (auto v = boost::any_cast<DecayType>(&value))
            out << *v;
        else if (auto v = boost::any_cast<TreatmentType>(&value))
            out << *v;
        else if (auto v = boost::any_cast<std::vector<std::string> >(&value))
            out << *v;
        else
            out << "error";
        out << std::endl;
    }

    return out;
}

} // namespace parser
