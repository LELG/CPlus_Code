#!/usr/bin/env python
"""
A script for crawling a results directory,
generating a summary results file (CSV) from the data,
and compressing the raw results.

Author: Yoshua Wakeham
        yoshua.wakeham@petermac.org
        y.wakeham@student.unimelb.edu.au
"""

import argparse
import os, sys, shutil
import csv


def main():
    args = parse_arguments()

    results_dir = args.results_dir.rstrip('/')

    if not os.path.isdir(results_dir):
        err_template = "error: directory '{}' does not exist. Aborting."
        print(err_template.format(results_dir))
        sys.exit(1)

    write_summary_file(results_dir)

    if args.compress:
        compress_results(results_dir)


def parse_arguments():
    """
    Wrapper around argparse.ArgumentParser().
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("results_dir", help="results directory to summarise")
    parser.add_argument("--compress", action="store_true",
                        help="compress results to a gzipped tar archive")

    return parser.parse_args()


def write_summary_file(results_dir):
    print("generating summary for results in dir '{}'".format(results_dir))
    try:
        summ_fields, summary_dicts = get_sim_summaries(results_dir)
    except ValueError:
        # there are no simulations to summarise
        print("summary failed: no simulations to summarise")
        return

    summ_fpath = get_summary_fpath(results_dir)

    with open(summ_fpath, 'w') as summf:
        writer = csv.DictWriter(summf, fieldnames=summ_fields)
        writer.writeheader()
        for summary in summary_dicts:
            writer.writerow(summary)

    print("summary file generated")


def compress_results(results_dir):
    """
    Compress entire results directory into
    a single archive, and delete the uncompressed results.

    Note that on decompression, duplicate copies of
    any files in the main results directory (e.g. summary file)
    will be created.
    """
    print("compressing results in dir: {}".format(results_dir))

    sim_id = get_simulation_id(results_dir)
    arch_name = os.path.join(results_dir, sim_id+'_results')
    shutil.make_archive(arch_name, 'gztar', root_dir=results_dir)

    param_set_dirs = get_subdirs(results_dir)
    for ps_dir in param_set_dirs:
        shutil.rmtree(ps_dir)

    print("compression complete")


def get_sim_summaries(results_dir):
    """
    Return list of summary fields, and a list of simulation summaries (dicts).

    Walk the results directory, generating a summary dictionary
    for each simulation run. Return a 2-tuple: a list of all fields for the
    summary CSV file, and a list of all the generated summaries.

    Throws ValueError if no simulation results are found in the results dir.
    """
    summary_dicts = []

    param_set_dirs = get_subdirs(results_dir)
    for ps_dir in param_set_dirs:
        param_set = get_param_set(ps_dir)
        run_dirs = get_subdirs(ps_dir)
        for run_dir in run_dirs:
            run_summary = get_run_summary(run_dir)
            run_summary['param_set'] = param_set
            summary_dicts.append(run_summary)

    try:
        summ_fields = get_summary_fields(summary_dicts)
    except ValueError:
        raise

    return summ_fields, summary_dicts


def get_simulation_id(results_dir):
    return os.path.split(results_dir)[1]


def get_param_set(ps_dir):
    return os.path.split(ps_dir)[1]


def get_run_summary(run_dir):
    summary = {}
    summary['run_number'] = get_run_number(run_dir)
    fields_funcs_map = get_run_fields()
    for fieldname, get_field_value in fields_funcs_map.items():
        summary[fieldname] = get_field_value(run_dir)
    return summary


def get_run_fields():
    name_func_mapping = {'runtime': get_runtime,
                         'elapsed_cycles': get_elapsed_cycles}
    return name_func_mapping


def get_runtime(run_dir):
    # TODO implement this properly; currently just a stub
    # open results file, grab the 'elapsed time' value
    return ord(run_dir[-1]) * 100


def get_elapsed_cycles(run_dir):
    # TODO implement this properly
    return ord(run_dir[-1]) * 20000


def get_run_number(run_dir):
    return os.path.split(run_dir)[1]


def get_subdirs(dirpath):
    return [os.path.join(dirpath, name) for name
            in os.listdir(dirpath)
            if os.path.isdir(os.path.join(dirpath, name))]


def get_summary_fields(summaries):
    """Get a list of all fields (column names) of summary file."""
    try:
        summ_fields = summaries[0].keys()
    except IndexError:
        raise ValueError('list of simulation summaries is empty')
    # make sure that param set and run number are first columns in CSV file
    summ_fields.remove('param_set')
    summ_fields.remove('run_number')
    summ_fields = ['param_set', 'run_number'] + summ_fields
    return summ_fields


def get_summary_fpath(results_dir):
    simulation_id = get_simulation_id(results_dir)
    summary_fname = '{}_summary.csv'.format(simulation_id)
    return os.path.join(results_dir, summary_fname)


if __name__ == "__main__":
    main()
