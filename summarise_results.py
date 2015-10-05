#!/usr/bin/env python2.7
"""
A script for crawling a simulation results directory,
generating summary files, plots and reports from the data,
and, optionally, compressing the raw results.

AUTHOR
    Yoshua Wakeham
    y.wakeham@student.unimelb.edu.au
    yoshua.wakeham@petermac.org
"""
from __future__ import print_function
import argparse
import os
import sys
import shutil
import tarfile
import pandas as pd
import re
import csv

import vis_utils
from summary_utils import RunSummary, write_summaries_to_file


def main():
    args = parse_arguments()

    results_dir = args.results_dir.rstrip('/')

    if not os.path.isdir(results_dir):
        err_template = "error: directory '{}' does not exist. Aborting."
        print(err_template.format(results_dir))
        sys.exit(1)

    summarise_sim_group(results_dir)

    if args.compress:
        compress_results(results_dir)


def parse_arguments():
    """
    Parse command line arguments.

    This function is just a simple wrapper around an argparse.ArgumentParser.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("results_dir", help="results directory to summarise")
    parser.add_argument("--compress", action="store_true",
                        help="compress results to a gzipped tar archive")

    return parser.parse_args()


def compress_results(results_dir):
    """
    Compress entire results directory into an archive,
    and delete uncompressed subdirectories.
    """
    print("compressing results in '{}'".format(results_dir))

    cwd = os.getcwd()

    try:
        os.chdir(results_dir)
        arch_name = 'raw_results.tar.gz'

        param_set_dirs = get_param_set_subdirs(os.getcwd())

        arch = tarfile.open(arch_name, 'w:gz')

        for ps_dir in param_set_dirs:
            # get relative directory
            ps_dir = os.path.split(ps_dir)[1]

            print("adding dir '{}' to archive ...".format(ps_dir), end='\r')
            sys.stdout.flush()
            arch.add(ps_dir)
            print("adding dir '{}' to archive ... done".format(ps_dir))

            print("deleting uncompressed dir '{}' ...".format(ps_dir), end='\r')
            sys.stdout.flush()
            shutil.rmtree(ps_dir)
            print("deleting uncompressed dir '{}' ... done".format(ps_dir))

        arch.close()
        print("compression complete")

    finally:
        os.chdir(cwd)


def summarise_sim_group(results_dir):
    """
    Write a summary file for the simulation results stored in results_dir.
    """

    print("generating summary files for '{}'".format(results_dir))

    summary_dir = os.path.join(results_dir, 'summary')

    if not os.path.isdir(summary_dir):
        os.mkdir(summary_dir)

    param_set_dirs = get_param_set_subdirs(results_dir)

    n_psets = len(param_set_dirs)
    curr = 0

    for ps_dir in param_set_dirs:
        curr += 1
        print("summarising ps {} / {}".format(curr, n_psets))
        summarise_param_set(ps_dir, summary_dir)


def summarise_param_set(ps_dir, summary_dir):
    """
    Summarise the param set stored in ps_dir, writing results to summary_dir.
    """
    ps_id = get_param_set(ps_dir)

    print("plotting growth curves")
    growth_data = get_growth_data(ps_dir)

    fig_dir = os.path.join(summary_dir, 'fig')
    if not os.path.isdir(fig_dir):
        os.mkdir(fig_dir)

    growth_plot = vis_utils.plot_growth_curves(growth_data)
    plot_fpath = os.path.join(fig_dir, 'ps{}_growthcurves.png'.format(ps_id))
    growth_plot.savefig(plot_fpath)

    summaries = get_param_set_summaries(ps_dir)

    print("writing summaries to CSV file")
    summary_fpath = os.path.join(summary_dir, 'ps{}_summary.csv'.format(ps_id))
    write_summaries_to_file(summaries, summary_fpath)

    print("generating html report")
    vis_utils.make_html_report(ps_id, summaries, summary_dir)


def get_growth_data(ps_dir):
    """
    Store growth curve data from a set of simulations to a pandas.DataFrame.
    """
    run_dirs = get_run_subdirs(ps_dir)

    all_data = []
    for run_dir in run_dirs:
        growth_fpath = os.path.join(run_dir, 'tumour_size.txt')
        df = pd.read_csv(growth_fpath, sep='\t', header=None)
        # actual growth data is in first column
        s = df[0]
        all_data.append(s)

    growth_data = pd.concat(all_data, axis=1)
    growth_data.columns = ['run'+str(n) for n in range(1, len(run_dirs)+1)]

    return growth_data


def get_param_set_summaries(ps_dir):
    """
    Return a list of results summaries, one for each run in a param set.
    """
    summaries = []

    run_dirs = get_run_subdirs(ps_dir)

    ps_conf_fpath = get_conf_fpath(ps_dir)

    nruns = len(run_dirs)
    curr = 0

    for run_dir in run_dirs:
        curr += 1
        print("getting summary for run {} / {}".format(curr, nruns))

        run_summary = generate_run_summary(run_dir)
        run_summary.add_param_field('param_set', get_param_set(ps_dir))
        add_fields_from_conf_file(run_summary, ps_conf_fpath)
        summaries.append(run_summary)

    return summaries


def generate_run_summary(run_dir):
    """
    Create and populate a RunSummary object for an individual simulation run.
    """
    summary = RunSummary()
    summary.add_param_field('run_number', get_run_number(run_dir))
    add_fields_from_run_dir(summary, run_dir)
    return summary


def add_fields_from_run_dir(summary, run_dir):
    """
    Add fields to a run summary from files in the replicate run directory.
    """
    results_fpath = os.path.join(run_dir, 'results.txt')
    if not os.path.isfile(results_fpath):
        raise IOError("cannot find results for run_dir: " + run_dir)
    summary.add_results_from_csv(results_fpath, delim='\t')

    stats_fpath = os.path.join(run_dir, 'stats.txt')
    if not os.path.isfile(stats_fpath):
        raise IOError("cannot find clone stats for run_dir: " + run_dir)
    add_result_fields_from_clone_stats(summary, stats_fpath)


def add_fields_from_conf_file(summary, conf_fpath):
    """
    Add parameters from a conf file to the summary.
    """
    #params_to_store = ['pro', 'die', 'mut', 'qui',
    #                   'prob_mut_pos', 'prob_mut_neg',
    #                   'prob_inc_mut', 'prob_dec_mut',
    #                   'driver_quantile', 'killer_quantile',
    #                   'beneficial_quantile', 'deleterious_quantile']

    summary.add_params_from_conf(conf_fpath)

    # rename initial prolif and mutation rates
    summary.add_param_field('prolif_init', summary.param_fields.pop('pro'))
    summary.add_param_field('mut_init', summary.param_fields.pop('mut'))


def add_result_fields_from_clone_stats(summary, stats_fpath, delim='\t'):
    """
    Add summary fields computed from the clone statistics file.
    """
    with open(stats_fpath) as stats_file:
        reader = csv.DictReader(stats_file, delimiter=delim)

        tumour_size = 0
        dominant_clone_size = 0
        agg_prolif = 0.0
        agg_mut = 0.0
        nclones = 0

        for row in reader:
            clone_size = int(row['Clone_size'])
            tumour_size += clone_size
            if clone_size > dominant_clone_size:
                dominant_clone_size = clone_size

            agg_prolif += float(row['Proliferation_Rate'])
            agg_mut += float(row['Mutation_Rate'])
            nclones += 1

    summary.add_result_field('prolif_final_avg', agg_prolif/nclones)
    summary.add_result_field('mut_final_avg', agg_mut/nclones)
    summary.add_result_field('dom_clone_proportion', float(dominant_clone_size)/tumour_size)


def get_conf_fpath(sim_dir):
    """
    Search a directory for a config file, returning its filepath if found.

    Raises a ValueError if either no conf file, or
    more than one, is found.
    """
    conf_fpath = None

    for name in os.listdir(sim_dir):
        rel_path = os.path.join(sim_dir, name)
        if is_conf_fpath(rel_path):
            if not conf_fpath:
                conf_fpath = rel_path
            else:
                raise ValueError("multiple conf files in dir '" + sim_dir + "'")

    if not conf_fpath:
        raise ValueError("no conf file in dir '" + sim_dir + "'")

    return conf_fpath


def is_conf_fpath(fpath):
    has_conf_extension = os.path.splitext(fpath)[1] == '.conf'
    return os.path.isfile(fpath) and has_conf_extension


def get_simulation_id(results_dir):
    """
    Return an ID string for the simulation test group.
    """
    return os.path.split(results_dir)[1]


def get_param_set(ps_dir):
    """
    Return an ID string for a simulation parameter set.
    """
    ps_id = os.path.split(ps_dir)[1]

    if ps_id.startswith('ps'):
        ps_id = ps_id[2:]

    return ps_id


def get_run_number(run_dir):
    """
    Return the number of the replicate run (as a string).
    """
    run_number = os.path.split(run_dir)[1]

    if run_number.startswith('run'):
        run_number = run_number[3:]

    return run_number


def get_param_set_subdirs(dirpath):
    """
    Get a list of param set subdirectories of a test group directory.
    """
    patt = '^' + dirpath + r'/ps'
    return get_subdirs(dirpath, patt)


def get_run_subdirs(dirpath):
    """
    Get a list of replicate run subdirectories of a param set directory.
    """
    patt = '^' + dirpath + r'/run'
    return get_subdirs(dirpath, patt)


def get_subdirs(dirpath, pattern=None):
    """
    Get a list of subdirectories of a given directory.
    """
    all_subdirs = [os.path.join(dirpath, name)
                   for name in os.listdir(dirpath)
                   if os.path.isdir(os.path.join(dirpath, name))]

    if pattern is not None:
        matching_subdirs = [d for d in all_subdirs if re.search(pattern, d)]
        return matching_subdirs
    else:
        return all_subdirs


if __name__ == "__main__":
    main()
