#!/usr/bin/env python2.7
"""
A script for crawling a simulation results directory,
generating a summary file (CSV format) from the data,
and, optionally, compressing the raw results.

AUTHOR
    Yoshua Wakeham
    y.wakeham@student.unimelb.edu.au
    yoshua.wakeham@petermac.org
"""
from __future__ import print_function
import argparse
import os, sys, shutil
import operator
from itertools import groupby
import csv
import re


def main():
    args = parse_arguments()

    results_dir = args.results_dir.rstrip('/')

    if not os.path.isdir(results_dir):
        err_template = "error: directory '{}' does not exist. Aborting."
        print(err_template.format(results_dir))
        sys.exit(1)

    generate_summary_files(results_dir)

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

    Note that on decompression, any files in the main results
    directory (e.g. summary file, conf file) will be duplicated.
    """
    print("compressing results in  '{}'".format(results_dir))

    sim_id = get_simulation_id(results_dir)
    arch_name = os.path.join(results_dir, sim_id+'_results')
    shutil.make_archive(arch_name, 'gztar', root_dir=results_dir)

    # delete uncompressed subdirectories
    param_set_dirs = get_subdirs(results_dir)
    for ps_dir in param_set_dirs:
        shutil.rmtree(ps_dir)

    print("compression complete")


def generate_summary_files(results_dir):
    """
    Write a summary file for the simulation results stored in results_dir.
    """
    print("generating summary files for '{}'".format(results_dir))

    summaries = get_sim_summaries(results_dir)

    summary_dir = os.path.join(results_dir, 'summary')

    sim_id = get_simulation_id(results_dir)
    tg_summ_fpath = os.path.join(summary_dir, sim_id+'_summary.csv')

    write_test_group_summary_file(summaries, tg_summ_fpath)
    print("wrote test group summary data")

    write_param_set_summary_files(tg_summ_fpath, summary_dir)
    print("wrote param set summary data")


def write_test_group_summary_file(summaries, summary_fpath):
    """
    Write a summary file for an entire group of simulations.
    """
    summ_fields = summaries[0].get_field_names()

    with open(summary_fpath, 'w') as summ_file:
        writer = csv.DictWriter(summ_file, fieldnames=summ_fields)
        writer.writeheader()
        for summary in summaries:
            writer.writerow(summary.fields)


def write_param_set_summary_files(summary_fpath, summary_dir):
    """
    Write a summary file for each param set of a simulation test group.

    This just parses rows from an existing test group summary file.
    """
    with open(summary_fpath) as tg_summ_file:
        tg_reader = csv.DictReader(tg_summ_file)
        header = tg_reader.fieldnames

        for ps, rows in groupby(tg_reader, operator.itemgetter('param_set')):
            ps_fname = "ps{}_summary.csv".format(ps)
            ps_fpath = os.path.join(summary_dir, ps_fname)

            with open(ps_fpath, "w") as ps_file:
                ps_writer = csv.DictWriter(ps_file, fieldnames=header)
                ps_writer.writeheader()
                for row in rows:
                    ps_writer.writerow(row)

def get_sim_summaries(results_dir):
    """
    Return a list of results summaries, one for each run in the test group.

    Walk the results directory, populating a RunSummary object
    for each simulation run. Return a list of all these summaries.
    """
    summaries = []

    # assume that all subdirectories of main results dir
    # are directories for individual param sets
    param_set_dirs = get_subdirs(results_dir)

    for ps_dir in param_set_dirs:
        summaries += get_param_set_summaries(ps_dir)

    return summaries


def get_param_set_summaries(ps_dir):
    """
    Return a list of results summaries, one for each run in a param set.
    """
    summaries = []

    # as above, assume that all subdirectories of
    # param set dir are directories for individual runs
    run_dirs = get_subdirs(ps_dir)

    ps_conf_fpath = get_conf_fpath(ps_dir)

    for run_dir in run_dirs:
        run_summary = generate_run_summary(run_dir)
        run_summary.add_field('param_set', get_param_set(ps_dir))
        run_summary.add_fields_from_conf_file(ps_conf_fpath)
        summaries.append(run_summary)

    return summaries


def generate_run_summary(run_dir):
    """
    Create and populate a RunSummary object for an individual simulation run.
    """
    summary = RunSummary()
    summary.add_field('run_number', get_run_number(run_dir))
    summary.add_fields_from_run_dir(run_dir)
    return summary


class RunSummary(object):
    """
    A summary of a simulation run.

    This class is a straightforward wrapper around
    a dictionary of summary fields. It provides various
    methods for scraping summary data from config files,
    results files etc.

    To add a new field or set of fields to the summary
    files generated by this script, add extra methods or
    method calls to this class.
    """
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
        Get a sorted list of the field names for this summary.

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
        Add parameters from a configuration file to the summary.

        This method assumes the config file is a list of

            param = val

        assignments (the spaces around the '=' being optional).

        Note that this method hardcodes the list of parameters
        to be added to the summary.
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
        Add fields to summary from files contained in a replicate run directory.
        """
        results_fpath = os.path.join(run_dir, 'results.txt')
        if not os.path.isfile(results_fpath):
            raise IOError("cannot find results for run_dir: " + run_dir)
        self.add_fields_from_csv_file(results_fpath, delim='\t')

    def add_fields_from_csv_file(self, fpath, delim=','):
        """
        Add all fields in a csv file to the run summary.

        This assumes that every field in the CSV file should
        be stored to the summary. It will also only take values
        from the first row of the file; any later rows will be ignored.
        """
        with open(fpath) as results_f:
            reader = csv.DictReader(results_f, delimiter=delim)
            row = next(reader)
            for field, val in row.items():
                self.add_field(field, val)


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
    return os.path.split(ps_dir)[1]


def get_run_number(run_dir):
    """
    Return the number of the replicate run (as a string).
    """
    return os.path.split(run_dir)[1]


def get_subdirs(dirpath):
    """
    Get a list of subdirectories of a given directory.
    """
    return [os.path.join(dirpath, name)
            for name
            in os.listdir(dirpath)
            if os.path.isdir(os.path.join(dirpath, name))]


if __name__ == "__main__":
    main()
