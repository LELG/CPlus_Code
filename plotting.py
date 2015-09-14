#!/usr/bin/env python2.7
"""
Define various plotting functions for summarising simulation results.

AUTHOR
    Yoshua Wakeham
    y.wakeham@student.unimelb.edu.au
    yoshua.wakeham@petermac.org

    'histogram_boxplot' is adapted from code by Luis Lara.
"""
import pandas as pd
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
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    data.plot(ax=ax0, kind='hist', bins=20)
    ax0.set_ylabel('')

    data.plot(ax=ax1, kind='box', vert=False, widths=0.8)
    ax1.set_yticklabels([])

    title_font = {'weight': 'bold', 'size': 16}
    fig.suptitle(title, fontdict=title_font);

    return fig
