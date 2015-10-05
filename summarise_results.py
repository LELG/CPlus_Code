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

import vis_utils
from summary_utils import RunSummary, write_summaries_to_file, make_report


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
    make_report(ps_id, summaries, summary_dir)


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
        run_summary.add_fields_from_conf_file(ps_conf_fpath)
        summaries.append(run_summary)

    return summaries


def generate_run_summary(run_dir):
    """
    Create and populate a RunSummary object for an individual simulation run.
    """
    summary = RunSummary()
    summary.add_param_field('run_number', get_run_number(run_dir))
    summary.add_fields_from_run_dir(run_dir)
    return summary


# class RunSummary(object):
    # """
    # A summary of a simulation run.

    # This class is a straightforward wrapper around
    # a dictionary of summary fields. It provides various
    # methods for scraping summary data from config files,
    # results files etc.

    # To add a new field or set of fields to the summary
    # files generated by this script, add extra methods or
    # method calls to this class.
    # """
    # def __init__(self):
        # #self.fields = {}
        # self.param_fields = {}
        # self.result_fields = {}

    # @property
    # def param_fields(self):
        # return self.__param_fields

    # @param_fields.setter
    # def param_fields(self, val):
        # if not hasattr(self, '__param_fields') and type(val) == dict:
            # self.__param_fields = val
        # else:
            # print("warning: attempting to overwrite RunSummary fields dict.")
            # print("         You probably don't want to do this. If you are")
            # print("         sure you want to do this, I direct you to")
            # print("         the _RunSummary__param_fields attribute.")

    # @property
    # def result_fields(self):
        # return self.__result_fields

    # @result_fields.setter
    # def result_fields(self, val):
        # if not hasattr(self, '__result_fields') and type(val) == dict:
            # self.__result_fields = val
        # else:
            # print("warning: attempting to overwrite RunSummary fields dict.")
            # print("         You probably don't want to do this. If you are")
            # print("         sure you want to do this, I direct you to")
            # print("         the _RunSummary__result_fields attribute.")

    # @property
    # def fields(self):
        # f = self.param_fields.copy()
        # f.update(self.result_fields)
        # return f

    # def add_param_field(self, field_name, value):
        # self.param_fields[field_name] = value

    # def add_result_field(self, field_name, value):
        # self.result_fields[field_name] = value

    # def get_field_names(self):
        # """
        # Get a sorted list of the field names for this summary.

        # Get a list of field names for this summary, sorted
        # in lexicographic order. (With the exception of 'param_set'
        # and 'run_number', which are shifted to the start of the list.)
        # These field names can then be passed to csv.DictWriter
        # for writing the header of a CSV file.
        # """
        # field_names = sorted(self.param_fields.keys() + self.result_fields.keys())

        # # make sure that param set and run number are first columns in CSV file
        # if 'run_number' in field_names:
            # field_names.remove('run_number')
            # field_names = ['run_number'] + field_names
        # if 'param_set' in field_names:
            # field_names.remove('param_set')
            # field_names = ['param_set'] + field_names

        # return field_names

    # def add_fields_from_conf_file(self, conf_fpath):
        # """
        # Add parameters from a configuration file to the summary.

        # This method assumes the config file is a list of

            # param = val

        # assignments (the spaces around the '=' being optional).

        # Note that this method hardcodes the list of parameters
        # to be added to the summary.
        # """
        # params_to_store = ['pro', 'die', 'mut', 'qui',
                           # 'prob_mut_pos', 'prob_mut_neg',
                           # 'prob_inc_mut', 'prob_dec_mut',
                           # 'driver_quantile', 'killer_quantile',
                           # 'beneficial_quantile', 'deleterious_quantile']

        # with open(conf_fpath) as conf_file:
            # pattern = r"(?P<param>[^\s]+)[\s]*=[\s]*(?P<val>[^\s]+)"
            # for line in conf_file:
                # match = re.search(pattern, line)
                # line_dict = match.groupdict()
                # line_param = line_dict['param']
                # if line_param in params_to_store:
                    # self.add_param_field(line_param, line_dict['val'])

        # # for param in params_to_store:
            # # if param not in self.fields:
                # # print("warning: config file missing expected parameter '{}'".format(param))

        # # rename initial prolif and mutation rates
        # self.add_param_field('prolif_init', self.param_fields.pop('pro'))
        # self.add_param_field('mut_init', self.param_fields.pop('mut'))

    # def add_fields_from_run_dir(self, run_dir):
        # """
        # Add fields to summary from files contained in a replicate run directory.
        # """
        # results_fpath = os.path.join(run_dir, 'results.txt')
        # if not os.path.isfile(results_fpath):
            # raise IOError("cannot find results for run_dir: " + run_dir)
        # self.add_result_fields_from_csv_file(results_fpath, delim='\t')

        # stats_fpath = os.path.join(run_dir, 'stats.txt')
        # if not os.path.isfile(stats_fpath):
            # raise IOError("cannot find clone stats for run_dir: " + run_dir)
        # self.add_result_fields_from_clone_stats(stats_fpath, delim='\t')

    # def add_result_fields_from_csv_file(self, fpath, delim=','):
        # """
        # Add all fields in a csv file to the run summary.

        # This assumes that every field in the CSV file should
        # be stored to the summary. It will also only take values
        # from the first row of the file. Any later rows will be ignored!
        # """
        # with open(fpath) as results_f:
            # reader = csv.DictReader(results_f, delimiter=delim)
            # row = next(reader)
            # for field, val in row.items():
                # self.add_result_field(field, val)

    # def add_result_fields_from_clone_stats(self, stats_fpath, delim='\t'):
        # """
        # Add summary fields computed from the clone statistics file.
        # """
        # with open(stats_fpath) as stats_file:
            # reader = csv.DictReader(stats_file, delimiter=delim)

            # tumour_size = 0
            # dominant_clone_size = 0
            # agg_prolif = 0.0
            # agg_mut = 0.0
            # nclones = 0

            # for row in reader:
                # clone_size = int(row['Clone_size'])
                # tumour_size += clone_size
                # if clone_size > dominant_clone_size:
                    # dominant_clone_size = clone_size

                # agg_prolif += float(row['Proliferation_Rate'])
                # agg_mut += float(row['Mutation_Rate'])
                # nclones += 1

        # self.add_result_field('prolif_final_avg', agg_prolif/nclones)
        # self.add_result_field('mut_final_avg', agg_mut/nclones)
        # self.add_result_field('dom_clone_proportion', float(dominant_clone_size)/tumour_size)


# def write_summaries_to_file(summaries, summary_fpath):
    # """
    # Write a summary file for an entire group of simulations.
    # """
    # summ_fields = summaries[0].get_field_names()

    # with open(summary_fpath, 'w') as summ_file:
        # writer = csv.DictWriter(summ_file, fieldnames=summ_fields)
        # writer.writeheader()
        # for summary in summaries:
            # writer.writerow(summary.fields)


# def generate_html_report(ps_id, summaries, summary_dir):
    # """
    # Generate a HTML report summarising per-param set results from a test group.
    # """
    # fields_to_summarise = ['num_clones', 'prolif_final_avg', 'mut_final_avg', 'dom_clone_proportion']

    # all_params = summaries[0].param_fields.copy()
    # all_params.pop('param_set')
    # all_params.pop('run_number')

    # field_data = defaultdict(list)

    # for field_name in fields_to_summarise:
        # for summary in summaries:
            # # assume that all fields we want to summarise will be numeric
            # field_data[field_name].append(float(summary.fields[field_name]))

    # summary_output = []

    # for field_name in field_data:
        # data = pd.Series(field_data[field_name])
        # plot_fpath = os.path.join(summary_dir, 'fig', 'ps{}_{}.png'.format(ps_id, field_name))
        # plot = summary_plotting.histogram_boxplot(data, '{} (ps{})'.format(field_name, ps_id))
        # plot.savefig(plot_fpath)
        # summary_info = {'field_name': field_name,
                        # 'plot_fpath': os.path.join('fig', os.path.split(plot_fpath)[1]),
                        # 'summary_data': pd.DataFrame(data.describe()).to_html()}
        # summary_output.append(summary_info)

    # TEMPLATE = """\
    # <!DOCTYPE html>
    # <html lang="en">
    # <head>
        # <title>{{ report_title }}</title>
    # </head>
    # <body>
        # <h1>{{ report_title }}</h1>

        # <h2>Parameters</h2>

        # <table>
        # {% for key, value in params.iteritems() %}
           # <tr>
                # <th> {{ key }} </th>
                # <td> {{ value }} </td>
           # </tr>
        # {% endfor %}
        # </table>

        # <hr>

        # <h2>Growth Curves</h2>
        # <img src="{{ gc_fig }}" alt="growth curves">

        # {% for summary in summary_dicts %}
        # <hr>

            # <h2>Summary: {{ summary['field_name'] }}</h2>
            # <img src="{{ summary['plot_fpath'] }}" alt="summary plot for {{ summary['field_name'] }}">
            # <p>{{ summary['summary_data'] }}</p>

        # {% endfor %}
    # </body>
    # </html>
    # """
    # report_template = jinja2.Template(TEMPLATE)

    # results_dir = os.path.split(summary_dir)[0]
    # report_title = "Summary Report for results in \"{}\" (param set {})".format(results_dir, ps_id)

    # report_fpath = os.path.join(summary_dir, 'ps{}_report.html'.format(ps_id))
    # gc_fpath = os.path.join('fig', 'ps{}_growthcurves.png'.format(ps_id))
    # with open(report_fpath, 'w') as report_file:
        # report = report_template.render(report_title=report_title,
                                        # params=all_params,
                                        # gc_fig=gc_fpath,
                                        # summary_dicts=summary_output)
        # report_file.write(report)


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
