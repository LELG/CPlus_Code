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
import re


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
    print("generating summary for '{}'".format(results_dir))
    try:
        summaries = get_sim_summaries(results_dir)
        summ_fields = summaries[0].get_field_names()
    except ValueError:
        # there are no simulations to summarise
        print("summary failed: no simulations to summarise")
        return

    summ_fpath = get_summary_fpath(results_dir)

    with open(summ_fpath, 'w') as summ_file:
        writer = csv.DictWriter(summ_file, fieldnames=summ_fields)
        writer.writeheader()
        for summary in summaries:
            writer.writerow(summary.fields)

    print("summary file generated")


def get_summary_fpath(results_dir):
    simulation_id = get_simulation_id(results_dir)
    summary_fname = '{}_summary.csv'.format(simulation_id)
    return os.path.join(results_dir, summary_fname)


def compress_results(results_dir):
    """
    Compress entire results directory into
    a single archive, and delete the uncompressed results.

    Note that on decompression, duplicate copies of
    any files in the main results directory (e.g. summary file)
    will be created.
    """
    print("compressing results in  '{}'".format(results_dir))

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
    summaries = []

    param_set_dirs = get_subdirs(results_dir)
    for ps_dir in param_set_dirs:
        ps_conf_fpath = get_conf_fpath(ps_dir)
        run_dirs = get_subdirs(ps_dir)
        for run_dir in run_dirs:
            run_summary = generate_run_summary(run_dir)
            run_summary.add_field('param_set', get_param_set(ps_dir))
            run_summary.add_fields_from_conf_file(ps_conf_fpath)
            summaries.append(run_summary)

    return summaries


def generate_run_summary(run_dir):
    summ = RunSummary()
    summ.add_field('run_number', get_run_number(run_dir))
    summ.add_fields_from_run_dir(run_dir)
    return summ


class RunSummary(object):
    def __init__(self):
        self.fields = {}

    @property
    def fields(self):
        return self.__fields

    @fields.setter
    def fields(self, val):
        if not hasattr(self, '__fields') and type(val) == dict:
            self.__fields = val
        else:
            print("warning: attempting to overwrite RunSummary fields dict.")
            print("         You probably don't want to do this. If you are")
            print("         sure you want to do this, I direct you to")
            print("         the _RunSummary__fields attribute.")

    def add_field(self, field_name, value):
        self.fields[field_name] = value

    def get_field_names(self):
        """
        Get a sorted list of field names for this summary.

        Get a list of field names for this summary, sorted
        in lexicographic order. (With the exception of 'param_set'
        and 'run_number', which are shifted to the start of the list.)
        These field names can then be passed to csv.DictWriter
        for writing the header of a CSV file.
        """
        field_names = sorted(self.fields.keys())
        # make sure that param set and run number are first columns in CSV file
        if 'run_number' in field_names:
            field_names.remove('run_number')
            field_names = ['run_number'] + field_names
        if 'param_set' in field_names:
            field_names.remove('param_set')
            field_names = ['param_set'] + field_names
        return field_names

    def add_fields_from_conf_file(self, conf_fpath):
        """
        Parse parameters from a configuration file and add them to summary.

        This method assumes the config file is a list of

            param = val

        assignments (the spaces around the '=' being optional).
        """
        params_to_store = ['pro', 'die', 'mut', 'qui',
                           'prob_mut_pos', 'prob_mut_neg',
                           'prob_inc_mut', 'prob_dec_mut',
                           'driver_quantile', 'killer_quantile',
                           'beneficial_quantile', 'deleterious_quantile']
        num_stored = 0

        with open(conf_fpath) as conf_file:
            pattern = r"(?P<param>[^\s]+)[\s]*=[\s]*(?P<val>[^\s]+)"
            for line in conf_file:
                match = re.search(pattern, line)
                line_dict = match.groupdict()
                line_param = line_dict['param']
                if line_param in params_to_store:
                    self.add_field(line_param, line_dict['val'])
                    num_stored += 1

        if num_stored != len(params_to_store):
            print("warning: config file missing some expected parameters")

    def add_fields_from_run_dir(self, run_dir):
        """
        Add summary fields from files contained in a replicate run directory.
        """
        results_fpath = os.path.join(run_dir, 'results.txt')
        self.add_fields_from_csv_file(results_fpath, delim='\t')

    def add_fields_from_csv_file(self, fpath, delim=','):
        """
        Add all fields in a csv file to the run summary.

        This assumes that every field in the results file
        should be stored to the summary. It will also only take values
        from the first row of the file; any later rows will be ignored.
        """
        with open(fpath) as results_f:
            reader = csv.DictReader(results_f, delimiter=delim)
            row = next(reader)
            for field, val in row.items():
                self.add_field(field, val)


def get_conf_fpath(sim_dir):
    """
    Get the path to the config file in this directory.
    """
    conf_files = []

    for name in os.listdir(sim_dir):
        rel_path = os.path.join(sim_dir, name)
        if os.path.isfile(rel_path) and os.path.splitext(rel_path)[1] == '.conf':
            conf_files.append(rel_path)

    if len(conf_files) != 1:
        raise IOError("no unique conf file in directory '" + sim_dir + "'")

    return conf_files[0]


def get_simulation_id(results_dir):
    return os.path.split(results_dir)[1]


def get_param_set(ps_dir):
    return os.path.split(ps_dir)[1]


def get_run_number(run_dir):
    return os.path.split(run_dir)[1]


def get_subdirs(dirpath):
    return [os.path.join(dirpath, name)
            for name
            in os.listdir(dirpath)
            if os.path.isdir(os.path.join(dirpath, name))]


if __name__ == "__main__":
    main()
