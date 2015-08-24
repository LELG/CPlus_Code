#!/usr/local/bin/python
"""
A script for crawling a results directory,
generating a summary results file (CSV) from the data,
and compressing the raw results.

Author: Yoshua Wakeham
        yoshua.wakeham@petermac.org
        y.wakeham@student.unimelb.edu.au
"""

import sys


def generate_summary(results_dir):
    print("generating summary for dir: {}".format(results_dir))


def compress_results(results_dir):
    print("compressing results in dir: {}".format(results_dir))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} <results dir>".format(sys.argv[0]))
        sys.exit(1)

    results_dir = sys.argv[1]
    generate_summary(results_dir)
    compress_results(results_dir)
