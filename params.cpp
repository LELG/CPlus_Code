#include "params.hpp"
#include <iostream>
#include <string>
#include <stdexcept>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

po::variables_map params::parse_options(int argc, char *argv[]) {
    po::options_description general{"General Options"};
    general.add_options()
        ("help,h", "print this help message and exit")
        ("test_group", po::value<std::string>()->value_name("NAME"),
             "name for this group of tests")
        ("param_set", po::value<unsigned>()->value_name("NUM"),
             "number of this param set in test group")
        ("run_number", po::value<unsigned>()->value_name("NUM"),
             "number of this replicate run")
        ("test_group_dir", po::value<std::string>()->value_name("DIR"),
             "directory for this test group")
        ("param_set_dir", po::value<std::string>()->value_name("DIR"),
             "directory for this particular param set")
        ("run_dir", po::value<std::string>()->value_name("DIR"),
             "directory for this replicate run")
        ("r_output",
             "produce output in R")
        ("auto_treatment",
             "automatically introduce treatment at size limit")
        ("no_plots",
             "do not generate plots")
        ("prune_clones",
             "remove dead clones from population")
        ("prune_muts",
             "remove dead mutations from population")
        ("mode", po::value<params::Mode>()->value_name("{in_vivo, cell_line}")->default_value(params::in_vivo),
            "simulation mode - defaults to in_vivo")
    ;

    po::options_description bounds{"Simulation Bounds"};
    bounds.add_options()
        ("max_cycles", po::value<unsigned int>()->default_value(1e5)->value_name("INT"),
            "maximum cycles in simulation")
        ("max_size_lim", po::value<unsigned int>()->default_value(1e5)->value_name("INT"),
            "population size limit (triggers treatment)")
        ("mid_doublings", po::value<unsigned int>()->default_value(30)->value_name("INT"),
            "cell_line mode: number of doublings at halfway point")
        ("max_doublings", po::value<unsigned int>()->default_value(60)->value_name("INT"),
            "cell_line mode: number of doublings at simulation end")
        ("doublings_per_flask", po::value<unsigned int>()->default_value(3)->value_name("INT"),
            "cell_line mode: passage flask after this many doublings")
    ;

    po::options_description probabilities{ "Probabilities" };
    probabilities.add_options()
        ("prob_mut_pos", po::value<double>()->default_value(0.01)->value_name("FLOAT"),
            "probability of beneficial mutation")
        ("prob_mut_neg", po::value<double>()->default_value(0.99)->value_name("FLOAT"),
            "probability of deleterious mutation")
        ("prob_inc_mut", po::value<double>()->default_value(0.0)->value_name("FLOAT"),
            "probability of increasing mutation rate")
        ("prob_dec_mut", po::value<double>()->default_value(0.0)->value_name("FLOAT"),
            "probability of decreasing mutation rate")
    ;

    po::options_description tumour_params{ "Tumour characteristics" };
    tumour_params.add_options()
        ("pro", po::value<double>()->default_value(0.04)->value_name("FLOAT"),
            "initial proliferation rate")
        ("die", po::value<double>()->default_value(0.03)->value_name("FLOAT"),
            "initial death rate")
        ("mut", po::value<double>()->default_value(0.001)->value_name("FLOAT"),
            "initial mutation rate")
/*  Need to deal with these mutually exclusive parameters!

    mutex_tumour = tumour_params.add_mutually_exclusive_group(required=True)
    mutex_tumour.add_argument('--init_size', type=pos_int)
    mutex_tumour.add_argument('--init_diversity',
                              action=store_flag_and_vars(['sub_file']),
                              metavar='SUB_FILE',
                              help='init_size and init_diversity are mutually exclusive')
*/
    ;

    po::options_description treatmt_params{ "Treatment parameters" };
    treatmt_params.add_options()
/*
    treatmt_params.add_argument('--treatment_type',
                                choices=('single_dose', 'metronomic',
                                         'adaptive', 'none'),
                                default='single_dose',
                                help='defaults to single_dose')
    treatmt_params.add_argument('--decay_type',
                                choices=('constant', 'linear', 'exp'),
                                default='constant',
                                help='defaults to constant')
    treatmt_params.add_argument('--decay_rate', type=float, default=0.0)
    treatmt_params.add_argument('--treatment_freq', type=pos_int, default=100)
    treatmt_params.add_argument('--adaptive_increment', type=float, default=0.001)
    treatmt_params.add_argument('--adaptive_threshold', type=float, default=0.025)
    treatmt_params.add_argument('--select_time', type=pos_int, default=400000)
    treatmt_params.add_argument('--select_pressure', type=float, default=0.01)
    treatmt_params.add_argument('--mutagenic_pressure', type=float, default=0.0)
*/
    ;

    po::options_description resist_params{ "Resistance parameters" };
    resist_params.add_options()
/*
    resist_params.add_argument('--resistance', action="store_true")
    resist_params.add_argument('--num_resist_mutns', type=nonneg_int)
    resist_params.add_argument('--resist_strength', type=float, default=1.0)
*/
    ;

    po::options_description saving_loading{ "Saving/Loading" };
    saving_loading.add_options()
/*
    mutex_saveload.add_argument('--save_snapshot', action="store_true")
    mutex_saveload.add_argument('--load_snapshot',
                                action=store_flag_and_vars(['snapshot_archive']),
                                default=False, metavar="SNAPSHOT_ARCHIVE")
*/
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
        po::notify(vm);
    }
    catch (std::exception &e)
    {
        std::cout << "\nError: " << e.what() << std::endl;
        std::cout << all << std::endl;
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
