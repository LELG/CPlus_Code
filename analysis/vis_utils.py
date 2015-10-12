"""
Define various plotting functions for summarising simulation results.

AUTHOR
    Yoshua Wakeham
    y.wakeham@student.unimelb.edu.au
    yoshua.wakeham@petermac.org

    'histogram_boxplot' is adapted from code by Luis Lara.
"""
from collections import defaultdict
import os
import jinja2
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
plt.style.use('ggplot')


def histogram_boxplot(data, title):
    """
    Generate a combined histogram and boxplot from a vector/series of data.

    data: a pandas.Series

    Returns a matplotlib Figure.
    """
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    data.plot(ax=ax0, kind='hist', bins=20)
    ax0.set_ylabel('')

    data.plot(ax=ax1, kind='box', vert=False, widths=0.8)
    ax1.set_yticklabels([])

    title_font = {'weight': 'bold', 'size': 16}
    fig.suptitle(title, fontdict=title_font)

    return fig


def plot_growth_curves(summaries):
    """
    Plot growth curves for a collection of RunSummary objects.
    """
    fig = plt.figure()
    ax0 = fig.add_subplot(111)

    growth_data = [summ.timeseries['pop_size'] for summ in summaries]

    df = pd.concat(growth_data, axis=1)
    df.columns = [summ.run_id for summ in summaries]

    df.plot(ax=ax0)

    title_font = {'weight': 'bold', 'size': 16}
    fig.suptitle("Population Size vs. Time", fontdict=title_font)

    return fig


def make_html_report(ps_id, summaries, summary_dir):
    """
    Generate a HTML report summarising per-param set results from a test group.
    """
    fields_to_summarise = ['num_clones', 'prolif_final_avg', 'mut_final_avg', 'dom_clone_proportion']

    all_params = summaries[0].param_fields.copy()
    all_params.pop('param_set')
    all_params.pop('run_number')

    field_data = defaultdict(list)

    for field_name in fields_to_summarise:
        for summary in summaries:
            # assume that all fields we want to summarise will be numeric
            field_data[field_name].append(float(summary.fields[field_name]))

    summary_output = []

    for field_name in field_data:
        data = pd.Series(field_data[field_name])
        plot_fpath = os.path.join(summary_dir, 'fig', 'ps{}_{}.png'.format(ps_id, field_name))
        plot = histogram_boxplot(data, '{} (ps{})'.format(field_name, ps_id))
        plot.savefig(plot_fpath)
        summary_info = {'field_name': field_name,
                        'plot_fpath': os.path.join('fig', os.path.split(plot_fpath)[1]),
                        'summary_data': pd.DataFrame(data.describe()).to_html()}
        summary_output.append(summary_info)

    # guard against memory leaks
    plt.close('all')

    TEMPLATE = """\
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <title>{{ report_title }}</title>
    </head>
    <body>
        <h1>{{ report_title }}</h1>

        <h2>Parameters</h2>

        <table>
        {% for key, value in params.iteritems() %}
           <tr>
                <th> {{ key }} </th>
                <td> {{ value }} </td>
           </tr>
        {% endfor %}
        </table>

        <hr>

        <h2>Growth Curves</h2>
        <img src="{{ gc_fig }}" alt="growth curves">

        {% for summary in summary_dicts %}
        <hr>

            <h2>Summary: {{ summary['field_name'] }}</h2>
            <img src="{{ summary['plot_fpath'] }}" alt="summary plot for {{ summary['field_name'] }}">
            <p>{{ summary['summary_data'] }}</p>

        {% endfor %}
    </body>
    </html>
    """
    report_template = jinja2.Template(TEMPLATE)

    results_dir = os.path.split(summary_dir)[0]
    report_title = "Summary Report for results in \"{}\" (param set {})".format(results_dir, ps_id)

    report_fpath = os.path.join(summary_dir, 'ps{}_report.html'.format(ps_id))
    gc_fpath = os.path.join('fig', 'ps{}_growthcurves.png'.format(ps_id))
    with open(report_fpath, 'w') as report_file:
        report = report_template.render(report_title=report_title,
                                        params=all_params,
                                        gc_fig=gc_fpath,
                                        summary_dicts=summary_output)
        report_file.write(report)
