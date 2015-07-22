#include "params.hpp"
#include <iostream>
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
 * Function to parse our command line options.
 */
po::variables_map params::parse_options(int argc, char *argv[])
{
    po::options_description general{"General Options"};
    general.add_options()
        ("help,h", "= print this help message and exit")
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
        ("run_dir", po::value<std::string>()->value_name("<dir>"),
             "= directory for this replicate run")
        ("r_output",
             "= produce output in R")
        ("auto_treatment",
             "= automatically introduce treatment at size limit")
        ("no_plots",
             "= do not generate plots")
        ("prune_clones",
             "= remove dead clones from population")
        ("prune_muts",
             "= remove dead mutations from population")
        ("mode", po::value<params::Mode>()->value_name("{in_vivo, cell_line}")->default_value(params::in_vivo),
            "= simulation mode - defaults to in_vivo")
    ;

    po::options_description bounds{"Simulation Bounds"};
    bounds.add_options()
        ("max_cycles", po::value<int>()->default_value(1e5)->value_name("<int>"),
            "= maximum cycles in simulation")
        ("max_size_lim", po::value<int>()->default_value(1e5)->value_name("<int>"),
            "= population size limit (triggers treatment)")
        ("mid_doublings", po::value<int>()->default_value(30)->value_name("<int>"),
            "= cell_line mode: number of doublings at halfway point")
        ("max_doublings", po::value<int>()->default_value(60)->value_name("<int>"),
            "= cell_line mode: number of doublings at simulation end")
        ("doublings_per_flask", po::value<int>()->default_value(3)->value_name("<int>"),
            "= cell_line mode: passage flask after this many doublings")
    ;

    po::options_description probabilities{ "Probabilities" };
    probabilities.add_options()
        ("prob_mut_pos", po::value<double>()->default_value(0.01)->value_name("<float>"),
            "= probability of beneficial mutation")
        ("prob_mut_neg", po::value<double>()->default_value(0.99)->value_name("<float>"),
            "= probability of deleterious mutation")
        ("prob_inc_mut", po::value<double>()->default_value(0.0)->value_name("<float>"),
            "= probability of increasing mutation rate")
        ("prob_dec_mut", po::value<double>()->default_value(0.0)->value_name("<float>"),
            "= probability of decreasing mutation rate")
    ;

    po::options_description tumour_params{ "Tumour characteristics" };
    tumour_params.add_options()
        ("pro", po::value<double>()->default_value(0.04)->value_name("<float>"),
            "= initial proliferation rate")
        ("die", po::value<double>()->default_value(0.03)->value_name("<float>"),
            "= initial death rate")
        ("mut", po::value<double>()->default_value(0.001)->value_name("<float>"),
            "= initial mutation rate")
        ("init_size", po::value<int>()->default_value(25)->value_name("<int>"),
            "= size of initial clone")
        ("init_diversity", po::value<std::string>()->value_name("<sub file>"),
            "= specify a heterogeneous initial population, stored in sub file")
    ;

    po::options_description treatmt_params{ "Treatment parameters" };
    treatmt_params.add_options()
        ("treatment_type", po::value<params::TreatmentType>()
                           ->value_name("{single_dose, metronomic, adaptive, none}")
                           ->default_value(params::single_dose),
            "= treatment type - defaults to single_dose")
        ("decay_type", po::value<params::DecayType>()
                       ->value_name("{constant, linear, exp}")
                       ->default_value(params::constant),
            "= treatment decay type - defaults to constant (no decay)")
        ("decay_rate", po::value<double>()->default_value(0.0)->value_name("<float>"),
            "= treatment decay rate")
        ("treatment_freq", po::value<int>()->default_value(100)->value_name("<int>"),
            "= treatment frequency")
        ("adaptive_increment", po::value<double>()->default_value(0.001)->value_name("<float>"),
            "= adaptive treatment: size of incremental dosage changes")
        ("adaptive_threshold", po::value<double>()->default_value(0.025)->value_name("<float>"),
            "= adaptive treatment: change in tumour size that triggers dosage change")
        ("select_time", po::value<int>()->default_value(400000)->value_name("<int>"),
            "= time to manually introduce selective pressure")
        ("select_pressure", po::value<double>()->default_value(0.01)->value_name("<float>"),
            "= initial value of selective pressure")
        ("mutagenic_pressure", po::value<double>()->default_value(0.0)->value_name("<float>"),
            "= initial value of mutagenic pressure")
    ;

    po::options_description resist_params{ "Resistance parameters" };
    resist_params.add_options()
        ("resistance", "= generate pre-existing resistance mutations")
        ("num_resist_mutns", po::value<int>()->value_name("<int>"),
            "= determine how many resistance mutations to generate")
        ("resist_strength", po::value<double>()->default_value(1.0)->value_name("<float>"),
            "= strength of resistance mutations (between 0 and 1)")
    ;

    po::options_description saving_loading{ "Saving/Loading" };
    saving_loading.add_options()
        ("save_snapshot", "= save a population snapshot when population reaches size limit")
        ("load_snapshot", po::value<std::string>()->value_name("<path to archive>"),
            "= load a snapshot from the specified archive")
    ;

    po::options_description all("SIMULATION OPTIONS");
    all.add(general).add(bounds).add(probabilities).add(tumour_params)
       .add(treatmt_params).add(resist_params).add(saving_loading);

    po::variables_map vm;

    try
    {
        po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
        // Display help text when requested
        if (vm.count("help"))
        {
            std::cout << all << std::endl;
            std::exit(EXIT_SUCCESS);
        }
        conflicting_options(vm, "init_size", "init_diversity");
        conflicting_options(vm, "save_snapshot", "load_snapshot");
        po::notify(vm);
    }
    catch (std::exception &e)
    {
        std::cout << "\nError: " << e.what() << std::endl;
        std::cout << "Run with --help for usage details" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return vm;
}

std::istream& params::operator>>(std::istream& in, params::Mode& mode)
{
    std::string token;
    in >> token;
    if (token == "in_vivo") {
        mode = params::in_vivo;
    } else if (token == "cell_line") {
        mode = params::cell_line;
    } else {
        throw po::invalid_option_value(token);
    }
    return in;
}

std::istream& params::operator>>(std::istream& in, params::TreatmentType& treatmt_type)
{
    std::string token;
    in >> token;
    if (token == "single_dose") {
        treatmt_type = params::single_dose;
    } else if (token == "metronomic") {
        treatmt_type = params::metronomic;
    } else if (token == "adaptive") {
        treatmt_type = params::adaptive;
    } else if (token == "none") {
        treatmt_type = params::none;
    } else {
        throw po::invalid_option_value(token);
    }
    return in;
}

std::istream& params::operator>>(std::istream& in, params::DecayType& decay)
{
    std::string token;
    in >> token;
    if (token == "constant") {
        decay = params::constant;
    } else if (token == "linear") {
        decay = params::linear;
    } else if (token == "exp") {
        decay = params::exp;
    } else {
        throw po::invalid_option_value(token);
    }
    return in;
}
