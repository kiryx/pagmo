#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Generate performance scaling figures.

The figures show the scaling of the performance in terms of ERT w.r.t.
dimensionality on a log-log scale. On the y-axis, data is represented as
a number of function evaluations divided by dimension, this is in order
to compare at a glance with a linear scaling for which ERT is
proportional to the dimension and would therefore be represented by a
horizontal line in the figure.

Crosses (+) give the median number of function evaluations of successful
trials divided by dimension for the smallest *reached* target function
value.
Numbers indicate the number of successful runs for the smallest
*reached* target.
If the smallest target function value (1e-8) is not reached for a given
dimension, crosses (x) give the average number of overall conducted
function evaluations divided by the dimension.

Horizontal lines indicate linear scaling with the dimension, additional
grid lines show quadratic and cubic scaling.

**Example**

.. plot::
    :width: 50%
    
    from pylab import *
    
    import PyGMO.util._bbob_pproc as bb
    
    # Scaling figure
    ds = bb.load('<path to folder with data>')
    figure()
    bb.ppfigdim.plot(ds)
    bb.ppfigdim.beautify()
"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from pdb import set_trace
from PyGMO.util._bbob_pproc import genericsettings, toolsstats, pproc
from PyGMO.util._bbob_pproc.ppfig import saveFigure, groupByRange

values_of_interest = pproc.TargetValues((10, 1, 1e-1, 1e-2, 1e-3, 1e-5, 1e-8))  # to rename!?
xlim_max = None
ynormalize_by_dimension = True  # not at all tested yet

global_dimensions = []

styles = [  # sort of rainbow style, most difficult (red) first
          {'color': 'r', 'marker': 'o', 'markeredgecolor': 'k', 'markeredgewidth': 2, 'linewidth': 4},
          {'color': 'm', 'marker': '.', 'linewidth': 4},
          {'color': 'y', 'marker': '^', 'markeredgecolor': 'k', 'markeredgewidth': 2, 'linewidth': 4},
          {'color': 'g', 'marker': '.', 'linewidth': 4},
          {'color': 'c', 'marker': 'v', 'markeredgecolor': 'k', 'markeredgewidth': 2, 'linewidth': 4},
          {'color': 'b', 'marker': '.', 'linewidth': 4},
          {'color': 'k', 'marker': 'o', 'markeredgecolor': 'k', 'markeredgewidth': 2, 'linewidth': 4},
        ] 

refcolor = 'wheat'

def beautify(axesLabel=True):
    """Customize figure presentation.
    
    """
    # Input checking

    # Get axis handle and set scale for each axis
    axisHandle = plt.gca()
    axisHandle.set_xscale("log")
    axisHandle.set_yscale("log")

    # Grid options
    axisHandle.xaxis.grid(False, which='major')
    # axisHandle.grid(True, which='major')
    axisHandle.grid(False, which='minor')
    # axisHandle.xaxis.grid(True, linewidth=0, which='major')
    ymin, ymax = plt.ylim()

    axisHandle.yaxis.grid(True, which='major')
    # quadratic slanted "grid"
    for i in xrange(-2, 7, 1 if ymax <= 1e3 else 2):
        plt.plot((0.2, 200), (10**i, 10**(i + 3)), 'k:', linewidth=0.5)  # TODO: this should be done before the real lines are plotted? 

    dimticklist = global_dimensions
    dimannlist = global_dimensions
    # TODO: All these should depend on one given input (xlim, ylim)

    axisHandle.set_xticks(dimticklist)
    axisHandle.set_xticklabels([str(n) for n in dimannlist])
    logyend = 11  # int(1 + np.log10(plt.ylim()[1]))
    axisHandle.set_yticks([10.**i for i in xrange(0, logyend)])
    axisHandle.set_yticklabels(range(0, logyend))
    if 11 < 3:
        tmp = axisHandle.get_yticks()
        tmp2 = []
        for i in tmp:
            tmp2.append('%d' % round(np.log10(i)))
        axisHandle.set_yticklabels(tmp2)
    if 11 < 3:
        # ticklabels = 10**np.arange(int(np.log10(plt.ylim()[0])), int(np.log10(1 + plt.ylim()[1])))
        ticks = []
        for i in xrange(int(np.log10(plt.ylim()[0])), int(np.log10(1 + plt.ylim()[1]))):
            ticks += [10 ** i, 2 * 10 ** i, 5 * 10 ** i]
        axisHandle.set_yticks(ticks)
        # axisHandle.set_yticklabels(ticklabels)
    # axes limites
    plt.xlim(0.9 * global_dimensions[0], 1.125 * global_dimensions[-1]) 
    plt.ylim(ymin=np.max((ymin, 10**-0.2)), ymax=int(ymax + 1))  # Set back the default maximum.

    if axesLabel:
        plt.xlabel('Dimension')
        if ynormalize_by_dimension:
            plt.ylabel('Run Lengths / Dimension')
        else:
            plt.ylabel('Run Lengths')
            

def generateData(dataSet, targetFuncValue):
    """Computes an array of results to be plotted.
    
    :returns: (ert, success rate, number of success, total number of
               function evaluations, median of successful runs).

    """

    it = iter(reversed(dataSet.evals))
    i = it.next()
    prev = np.array([np.nan] * len(i))

    while i[0] <= targetFuncValue:
        prev = i
        try:
            i = it.next()
        except StopIteration:
            break

    data = prev[1:].copy()  # keep only the number of function evaluations.
    succ = (np.isnan(data) == False)
    if succ.any():
        med = toolsstats.prctile(data[succ], 50)[0]
        # Line above was modified at rev 3050 to make sure that we consider only
        # successful trials in the median
    else:
        med = np.nan

    data[np.isnan(data)] = dataSet.maxevals[np.isnan(data)]

    res = []
    res.extend(toolsstats.sp(data, issuccessful=succ, allowinf=False))
    res.append(np.mean(data))  # mean(FE)
    res.append(med)

    return np.array(res)  

def plot_a_bar(x, y,
               plot_cmd=plt.loglog,
               rec_width=0.1, # box ("rectangle") width, log scale 
               rec_taille_fac=0.3,  # notch width parameter
               styles={'color': 'b'},
               linewidth=1,
               fill_color=None, # None means no fill
               fill_transparency=0.7  # 1 should be invisible
               ):
    """plot/draw a notched error bar, x is the x-position,
    y[0,1,2] are lower, median and upper percentile respectively. 
    
    hold(True) to see everything.

    TODO: with linewidth=0, inf is not visible
    
    """
    if not np.isfinite(y[2]):
        y[2] = y[1] + 100 * (y[1] - y[0])
        if plot_cmd in (plt.loglog, plt.semilogy):
            y[2] = (1 + y[1]) * (1 + y[1] / y[0])**10
    if not np.isfinite(y[0]):
        y[0] = y[1] - 100 * (y[2] - y[1])
        if plot_cmd in (plt.loglog, plt.semilogy):
            y[0] = y[1] / (1 + y[2] / y[1])**10
    styles2 = {}
    for s in styles:
        styles2[s] = styles[s]
    styles2['linewidth'] = linewidth
    styles2['markeredgecolor'] = styles2['color']
    dim = 1  # to remove error
    x0 = x
    if plot_cmd in (plt.loglog, plt.semilogx):
        r = np.exp(rec_width) # ** ((1. + i_target / 3.) / 4)  # more difficult targets get a wider box
        x = [x0 * dim / r, x0 * r * dim] # assumes log-scale of x-axis
        xm = [x0 * dim / (r**rec_taille_fac), x0 * dim * (r**rec_taille_fac)]
    else:
        r = rec_width
        x = [x0 * dim - r, x0 * dim + r]
        xm = [x0 * dim - (r * rec_taille_fac), x0 * dim + (r * rec_taille_fac)]

    y = np.array(y) / dim
    if fill_color is not None:
        plt.fill_between([x[0], xm[0], x[0], x[1], xm[1], x[1], x[0]],
                         [y[0], y[1],  y[2], y[2], y[1],  y[0], y[0]], 
                         color=fill_color, alpha=1-fill_transparency)
    plot_cmd([x[0], xm[0], x[0], x[1], xm[1], x[1], x[0]],
             [y[0], y[1],  y[2], y[2], y[1],  y[0], y[0]],
             markersize=0, **styles2)
    styles2['linewidth'] = 0
    plot_cmd([x[0], x[1], x[1], x[0], x[0]],
             [y[0], y[0], y[2], y[2], y[0]],
             **styles2)
    styles2['linewidth'] = 2  # median
    plot_cmd([x[0], x[1]], [y[1], y[1]],
             markersize=0, **styles2)

    
def plot(dsList, valuesOfInterest=values_of_interest, styles=styles):
    """From a DataSetList, plot a figure of ERT/dim vs dim.
    
    There will be one set of graphs per function represented in the
    input data sets. Most usually the data sets of different functions
    will be represented separately.
    
    :param DataSetList dsList: data sets
    :param seq valuesOfInterest: 
        target precisions via class TargetValues, there might 
        be as many graphs as there are elements in
        this input. Can be different for each
        function (a dictionary indexed by ifun). 
    
    :returns: handles

    """
    global global_dimensions
    valuesOfInterest = pproc.TargetValues.cast(valuesOfInterest)
    styles = list(reversed(styles[:len(valuesOfInterest)]))
    dsList = pproc.DataSetList(dsList)
    dictFunc = dsList.dictByFunc()
    res = []

    global_dimensions = sorted(dsList.dictByDim().keys())

    for func in dictFunc:
        dictFunc[func] = dictFunc[func].dictByDim()
        dimensions = sorted(dictFunc[func])
        # legend = []
        line = []
        mediandata = {}
        displaynumber = {}
        for i_target in range(len(valuesOfInterest)):
            succ = []
            unsucc = []
            # data = []
            maxevals = np.ones(len(dimensions))
            maxevals_succ = np.ones(len(dimensions)) 
            # Collect data that have the same function and different dimension.
            for idim, dim in enumerate(dimensions):
                assert len(dictFunc[func][dim]) == 1
                # (ert, success rate, number of success, total number of
                #        function evaluations, median of successful runs)
                tmp = generateData(dictFunc[func][dim][0], valuesOfInterest((func, dim))[i_target])
                maxevals[idim] = max(dictFunc[func][dim][0].maxevals)
                # data.append(np.append(dim, tmp))
                if tmp[2] > 0:  # Number of success is larger than 0
                    succ.append(np.append(dim, tmp))
                    if tmp[2] < dictFunc[func][dim][0].nbRuns():
                        displaynumber[dim] = ((dim, tmp[0], tmp[2]))
                    mediandata[dim] = (i_target, tmp[-1])
                    unsucc.append(np.append(dim, np.nan))
                else:
                    unsucc.append(np.append(dim, tmp[-2]))  # total number of fevals

            if len(succ) > 0:
                tmp = np.vstack(succ)
                # ERT
                if genericsettings.scaling_figures_with_boxes:
                    for dim in dimensions: 
                        # to find finite simulated runlengths we need to have at least one successful run
                        if dictFunc[func][dim][0].detSuccesses([valuesOfInterest((func, dim))[i_target]])[0]:
                            # make a box-plot
                            y = toolsstats.drawSP_from_dataset(
                                                dictFunc[func][dim][0],
                                                valuesOfInterest((func, dim))[i_target],
                                                [25, 50, 75], 
                                                genericsettings.simulated_runlength_bootstrap_sample_size)[0]
                            rec_width = 1.1 # box ("rectangle") width
                            rec_taille_fac = 0.3  # notch width parameter
                            r = rec_width ** ((1. + i_target / 3.) / 4)  # more difficult targets get a wider box
                            styles2 = {}
                            for s in styles[i_target]:
                                styles2[s] = styles[i_target][s]
                            styles2['linewidth'] = 1
                            styles2['markeredgecolor'] = styles2['color'] 
                            x = [dim / r, r * dim]
                            xm = [dim / (r**rec_taille_fac), dim * (r**rec_taille_fac)]
                            y = np.array(y) / dim
                            plt.plot([x[0], xm[0], x[0], x[1], xm[1], x[1], x[0]],
                                     [y[0], y[1],  y[2], y[2], y[1],  y[0], y[0]],
                                     markersize=0, **styles2)
                            styles2['linewidth'] = 0
                            plt.plot([x[0], x[1], x[1], x[0], x[0]],
                                     [y[0], y[0], y[2], y[2], y[0]],
                                     **styles2)
                            styles2['linewidth'] = 2  # median
                            plt.plot([x[0], x[1]], [y[1], y[1]],
                                     markersize=0, **styles2)
                # plot lines, we have to be smart to connect only adjacent dimensions
                for i, n in enumerate(tmp[:, 0]):
                    j = list(dimensions).index(n)
                    if i == len(tmp[:, 0]) - 1 or j == len(dimensions) - 1: 
                        break
                    if dimensions[j+1] == tmp[i+1, 0]:
                        res.extend(plt.plot(tmp[i:i+2, 0], tmp[i:i+2, 1] / tmp[i:i+2, 0]**ynormalize_by_dimension,
                                            markersize=0, clip_on=True, **styles[i_target]))
                # plot only marker
                lw = styles[i_target].get('linewidth', None) 
                styles[i_target]['linewidth'] = 0
                res.extend(plt.plot(tmp[:, 0], tmp[:, 1] / tmp[:, 0]**ynormalize_by_dimension,
                           markersize=20, clip_on=True, **styles[i_target]))
                # restore linewidth
                if lw:
                    styles[i_target]['linewidth'] = lw
                else:
                    del styles[i_target]['linewidth']

        # To have the legend displayed whatever happens with the data.
        for i in reversed(range(len(valuesOfInterest))):
            res.extend(plt.plot([], [], markersize=10,
                                label=valuesOfInterest.loglabel(i),
                                **styles[i]))
        # Only for the last target function value
        if unsucc:  # obsolete
            tmp = np.vstack(unsucc)  # tmp[:, 0] needs to be sorted!
            # res.extend(plt.plot(tmp[:, 0], tmp[:, 1]/tmp[:, 0],
            #            color=styles[len(valuesOfInterest)-1]['color'],
            #            marker='x', markersize=20))
        if 1 < 3: # maxevals
            ylim = plt.ylim()
            res.extend(plt.plot(tmp[:, 0], maxevals / tmp[:, 0]**ynormalize_by_dimension,
                       color=styles[len(valuesOfInterest) - 1]['color'],
                       ls='', marker='x', markersize=20))
            plt.ylim(ylim)
        # median
        if mediandata:
            # for i, tm in mediandata.iteritems():
            for i in displaynumber:  # display median where success prob is smaller than one
                tm = mediandata[i]
                plt.plot((i,), (tm[1] / i**ynormalize_by_dimension,), 
                         color=styles[tm[0]]['color'],
                         linestyle='', marker='+', markersize=30,
                         markeredgewidth=5, zorder= -1)

        a = plt.gca()
        # the displaynumber is emptied for each new target precision
        # therefore the displaynumber displayed below correspond to the
        # last target (must be the hardest)
        if displaynumber:  # displayed only for the smallest valuesOfInterest
            for _k, j in displaynumber.iteritems():
                # the 1.5 factor is a shift up for the digits 
                plt.text(j[0], 1.5 * j[1] / j[0]**ynormalize_by_dimension, 
                         "%.0f" % j[2], axes=a,
                         horizontalalignment="center",
                         verticalalignment="bottom", fontsize=plt.rcParams['font.size'] * 0.85)
        # if later the ylim[0] becomes >> 1, this might be a problem
    return res

def main(dsList, _valuesOfInterest=values_of_interest, outputdir="", verbose=True):
    """From a DataSetList, returns a convergence and ERT/dim figure vs dim.
    
    :param DataSetList dsList: data sets
    :param seq _valuesOfInterest: target precisions, either as list or as
                                  ``pproc.TargetValues`` class instance. 
                                  There will be as many graphs as there are 
                                  elements in this input. 
    :param string outputdir: output directory
    :param bool verbose: controls verbosity
    
    """

    # plt.rc("axes", labelsize=20, titlesize=24)
    # plt.rc("xtick", labelsize=20)
    # plt.rc("ytick", labelsize=20)
    # plt.rc("font", size=20)
    # plt.rc("legend", fontsize=20)

    _valuesOfInterest = pproc.TargetValues.cast(_valuesOfInterest)

    dictFunc = dsList.dictByFunc()

    for func in dictFunc:
        plot(dictFunc[func], _valuesOfInterest, styles=styles)  # styles might have changed via config
        beautify(axesLabel=True)
        plt.text(plt.xlim()[0], plt.ylim()[0], _valuesOfInterest.short_info, fontsize=14)
        plt.gca().set_title(str(func))
        filename = os.path.join(outputdir, 'ppfigdim_f%s' % str(func))
        saveFigure(filename, verbose=verbose)
        plt.close()
