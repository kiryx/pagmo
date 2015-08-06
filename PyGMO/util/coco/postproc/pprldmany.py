#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Generates figure of the bootstrap distribution of ERT.
    
The main method in this module generates figures of Empirical
Cumulative Distribution Functions of the bootstrap distribution of
the Expected Running Time (ERT) divided by the dimension for many
algorithms.

The outputs show the ECDFs of the running times of the simulated runs
divided by dimension for 50 different targets logarithmically uniformly
distributed in [1e−8, 1e2]. The crosses (×) give the median number of
function evaluations of unsuccessful runs divided by dimension.

**Example**

.. plot::
    :width: 50%

    import urllib
    import tarfile
    import glob
    from pylab import *
    
    import PyGMO.util.coco.postproc as bb
      
    # Empirical cumulative distribution function of bootstrapped ERT figure
    bb.pprldmany.plot(ds) # must rather call main instead of plot?
    bb.pprldmany.beautify()

"""

from __future__ import absolute_import

import os
import warnings
from pdb import set_trace
import numpy as np
import matplotlib.pyplot as plt
from PyGMO.util.coco.postproc import toolsstats, genericsettings
from PyGMO.util.coco.postproc import pproc as pp # import dictAlgByDim, dictAlgByFun 
from PyGMO.util.coco.postproc import toolsdivers  # strip_pathname, str_to_latex
from PyGMO.util.coco.postproc import pprldistr # plotECDF, beautifyECDF
from PyGMO.util.coco.postproc import ppfig  # consecutiveNumbers, saveFigure, plotUnifLogXMarkers, logxticks
from PyGMO.util.coco.postproc import pptex

target_values = pp.TargetValues(10**np.arange(2, -8, -0.2)) # possibly changed in config
x_limit = None  # not sure whether this is necessary/useful
x_limit_default = 1e7 # better: 10 * genericsettings.evaluation_setting[1], noisy: 1e8, otherwise: 1e7. maximal run length shown
divide_by_dimension = True
annotation_line_end_relative = 1.11  # lines between graph and annotation
annotation_space_end_relative = 1.24  # figure space end relative to x_limit
save_zoom = False  # save zoom into left and right part of the figures
perfprofsamplesize = genericsettings.simulated_runlength_bootstrap_sample_size_rld  # number of bootstrap samples drawn for each fct+target in the performance profile
dpi_global_var = 100  # 100 ==> 800x600 (~160KB), 120 ==> 960x720 (~200KB), 150 ==> 1200x900 (~300KB) looks ugly in latex
nbperdecade = 1
median_max_evals_marker_format = ['x', 24, 3]

styles = [d.copy() for d in genericsettings.line_styles]  # deep copy

show_algorithms = []

refcolor = 'wheat'

def plt_plot(*args, **kwargs):
    return plt.plot(*args, clip_on=False, **kwargs)
    
def beautify():
    """Customize figure presentation."""

    #plt.xscale('log') # Does not work with matplotlib 0.91.2
    a = plt.gca()
    a.set_xscale('log')
    #Tick label handling
    plt.xlim(xmin=1e-0)

    global divide_by_dimension
    if divide_by_dimension:
        plt.xlabel('log10 of (# f-evals / dimension)')
    else:
        plt.xlabel('log10 of # f-evals')
    plt.ylabel('Proportion of function+target pairs')
    ppfig.logxticks()
    pprldistr.beautifyECDF()

def plotdata(data, maxval=None, maxevals=None, CrE=0., **kwargs):
    """Draw a normalized ECDF. What means normalized?
    
    :param seq data: data set, a 1-D ndarray of runlengths
    :param float maxval: right-most value to be displayed, will use the
                         largest non-inf, non-nan value in data if not
                         provided
    :param seq maxevals: if provided, will plot the median of this
                         sequence as a single cross marker
    :param float CrE: Crafting effort the data will be multiplied by
                      the exponential of this value.
    :param kwargs: optional arguments provided to plot function.
    
    """

    #Expect data to be a ndarray.
    x = data[np.isnan(data)==False] # Take away the nans
    nn = len(x)

    x = x[np.isinf(x)==False] # Take away the infs
    n = len(x)

    x = np.exp(CrE) * x  # correction by crafting effort CrE

    if n == 0:
        #res = plt.plot((1., ), (0., ), **kwargs)
        res = pprldistr.plotECDF(np.array((1., )), n=np.inf, **kwargs)
    else:
        dictx = {} # number of appearances of each value in x
        for i in x: 
            dictx[i] = dictx.get(i, 0) + 1  

        x = np.array(sorted(dictx))  # x is not a multiset anymore
        y = np.cumsum(list(dictx[i] for i in x)) # cumsum of size of y-steps (nb of appearences)
        idx = sum(x <= x_limit**annotation_space_end_relative) - 1 
        y_last, x_last = y[idx] / float(nn), x[idx] 
        if maxval is None:
            maxval = max(x)
        end = np.sum(x <= maxval)
        x = x[:end]
        y = y[:end]
        
        try:  # plot the very last point outside of the "normal" plotting area
            c = kwargs['color']
            plt_plot([x_last] * 2, [y_last] * 2, '.', color=c, markeredgecolor=c) 
        except:
            pass
        x2 = np.hstack([np.repeat(x, 2), maxval]) # repeat x-values for each step in the cdf
        y2 = np.hstack([0.0, np.repeat(y / float(nn), 2)])

        res = ppfig.plotUnifLogXMarkers(x2, y2, nbperdecade * 3 / np.log10(maxval), 
                                  logscale=False, clip_on=False, **kwargs)
        # res = plotUnifLogXMarkers(x2, y2, nbperdecade, logscale=False, **kwargs)

        if maxevals: # Should cover the case where maxevals is None or empty
            x3 = np.median(maxevals)
            if (x3 <= maxval):
                # np.any(x2 <= x3) and   # maxval < median(maxevals)
                try:
                    y3 = y2[x2<=x3][-1]  # find right y-value for x3==median(maxevals)
                except IndexError:  # median(maxevals) is smaller than any data, can only happen because of CrE?
                    y3 = y2[0]
                h = plt.plot((x3,), (y3,), 
                             marker=median_max_evals_marker_format[0],
                             markersize=median_max_evals_marker_format[1],
                             markeredgewidth=median_max_evals_marker_format[2],
                             # marker='x', markersize=24, markeredgewidth=3, 
                             markeredgecolor=plt.getp(res[0], 'color'),
                             ls=plt.getp(res[0], 'ls'),
                             color=plt.getp(res[0], 'color'))
                h.extend(res)
                res = h # so the last element in res still has the label.
                # Only take sequences for x and y!

    return res

def plotLegend(handles, maxval):
    """Display right-side legend.
    
    :param float maxval: rightmost x boundary
    :returns: list of (ordered) labels and handles.

    The figure is stopped at maxval (upper x-bound), and the graphs in
    the figure are prolonged with straight lines to the right to connect
    with labels of the graphs (uniformly spread out vertically). The
    order of the graphs at the upper x-bound line give the order of the
    labels, in case of ties, the best is the graph for which the x-value
    of the first step (from the right) is smallest.
    
    The annotation string is stripped from preceeding pathnames. 

    """
    reslabels = []
    reshandles = []
    ys = {}
    lh = 0
    for h in handles:
        x2 = []
        y2 = []
        for i in h:
            x2.append(plt.getp(i, "xdata"))
            y2.append(plt.getp(i, "ydata"))

        x2 = np.array(np.hstack(x2))
        y2 = np.array(np.hstack(y2))
        tmp = np.argsort(x2)
        x2 = x2[tmp]
        y2 = y2[tmp]

        h = h[-1] # we expect the label to be in the last element of h
        tmp = (x2 <= maxval)
        try:
            x2bis = x2[y2 < y2[tmp][-1]][-1]
        except IndexError: # there is no data with a y smaller than max(y)
            x2bis = 0.
        ys.setdefault(y2[tmp][-1], {}).setdefault(x2bis, []).append(h)
        lh += 1

    if len(show_algorithms) > 0:
        lh = min(lh, len(show_algorithms))
    if lh <= 1:
        lh = 2
    fontsize = genericsettings.minmax_algorithm_fontsize[0] + np.min((1, np.exp(9-lh))) * (
        genericsettings.minmax_algorithm_fontsize[-1] - genericsettings.minmax_algorithm_fontsize[0])
    i = 0 # loop over the elements of ys
    for j in sorted(ys.keys()):
        for k in reversed(sorted(ys[j].keys())):
            #enforce best ever comes last in case of equality
            tmp = []
            for h in ys[j][k]:
                tmp.append(h)
            tmp.reverse()
            ys[j][k] = tmp

            for h in ys[j][k]:
                if (not plt.getp(h, 'label').startswith('_line') and
                    (len(show_algorithms) == 0 or
                     plt.getp(h, 'label') in show_algorithms)):
                    y = 0.02 + i * 0.96/(lh-1)
                    tmp = {}
                    for attr in ('lw', 'ls', 'marker',
                                 'markeredgewidth', 'markerfacecolor',
                                 'markeredgecolor', 'markersize', 'zorder'):
                        tmp[attr] = plt.getp(h, attr)
                    legx = maxval**annotation_line_end_relative
                    if 'marker' in attr:
                        legx = maxval**annotation_line_end_relative
                    # reshandles.extend(plt_plot((maxval, legx), (j, y),
                    reshandles.extend(plt_plot((maxval, legx), (j, y),
                                      color=plt.getp(h, 'markeredgecolor'), **tmp))
                    reshandles.append(plt.text(maxval**(0.02 + annotation_line_end_relative), y,
                                               plt.getp(h, 'label').split(os.sep)[-1],
                                               horizontalalignment="left",
                                               verticalalignment="center", size=fontsize))
                    reslabels.append(plt.getp(h, 'label'))
                    #set_trace()
                    i += 1

    #plt.axvline(x=maxval, color='k') # Not as efficient?
    reshandles.append(plt_plot((maxval, maxval), (0., 1.), color='k'))
    reslabels.reverse()
    plt.xlim(xmax=maxval**annotation_space_end_relative)
    return reslabels, reshandles

def plot(dsList, targets=None, craftingeffort=0., **kwargs):
    """This function is obsolete?
    Generates a graph of the run length distribution of an algorithm.

    We display the empirical cumulative distribution function ECDF of
    the bootstrapped distribution of the runlength for an algorithm
    (in number of function evaluations) to reach the target functions 
    value :py:data:`targets`.

    :param DataSetList dsList: data set for one algorithm
    :param seq targets: target function values
    :param float crafting effort: the data will be multiplied by the
                                  exponential of this value
    :param dict kwargs: additional parameters provided to plot function.
    
    :returns: handles

    """
    if targets is None:
        targets = target_values  # set above or in config.py
    try:
        if np.min(targets) >= 1:
            ValueError('smallest target f-value is not smaller than one, use ``pproc.TargetValues(targets)`` to prevent this error')
        targets = pp.TargetValues(targets)
    except TypeError:
        pass
    res = []
    assert len(pp.DataSetList(dsList).dictByDim()) == 1 # We never integrate over dimensions...
    data = []
    maxevals = []
    for entry in dsList:
        for t in targets((entry.funcId, entry.dim)):
            divisor = entry.dim if divide_by_dimension else 1
            x = [np.inf] * perfprofsamplesize
            runlengthunsucc = []
            evals = entry.detEvals([t])[0]
            runlengthsucc = evals[np.isnan(evals) == False] / divisor
            runlengthunsucc = entry.maxevals[np.isnan(evals)] / divisor
            if len(runlengthsucc) > 0:
                x = toolsstats.drawSP(runlengthsucc, runlengthunsucc,
                                     percentiles=[50],
                                     samplesize=perfprofsamplesize)[1]
            data.extend(x)
            maxevals.extend(runlengthunsucc)

    # Display data
    data = np.array(data)
    data = data[np.isnan(data)==False] # Take away the nans
    n = len(data)
    data = data[np.isinf(data)==False] # Take away the infs
    # data = data[data <= maxval] # Take away rightmost data
    data = np.exp(craftingeffort) * data  # correction by crafting effort CrE
    if len(data) == 0: # data is empty.
        res = pprldistr.plotECDF(np.array((1., )), n=np.inf, **kwargs)
    else:
        res = pprldistr.plotECDF(np.array(data), n=n, **kwargs)

    if maxevals: # Should cover the case where maxevals is None or empty
        x3 = np.median(maxevals)
        if np.any(data > x3):
            y3 = float(np.sum(data <= x3)) / n
            h = plt_plot((x3,), (y3,), marker='x', markersize=24, markeredgewidth=3,
                         markeredgecolor=plt.getp(res[0], 'color'),
                         ls='', color=plt.getp(res[0], 'color'))
            h.extend(res)
            res = h # so the last element in res still has the label.
    return res

def main(dictAlg, order=None, outputdir='.', info='default',
         dimension=None, verbose=True):
    """Generates a figure showing the performance of algorithms.

    From a dictionary of :py:class:`DataSetList` sorted by algorithms,
    generates the cumulative distribution function of the bootstrap
    distribution of ERT for algorithms on multiple functions for
    multiple targets altogether.

    :param dict dictAlg: dictionary of :py:class:`DataSetList` instances
                         one instance is equivalent to one algorithm,
    :param list targets: target function values
    :param list order: sorted list of keys to dictAlg for plotting order
    :param str outputdir: output directory
    :param str info: output file name suffix
    :param bool verbose: controls verbosity

    """
    global x_limit  # late assignment of default, because it can be set to None in config 
    global divide_by_dimension  # not fully implemented/tested yet
    if 'x_limit' not in globals() or x_limit is None:
        x_limit = x_limit_default

    tmp = pp.dictAlgByDim(dictAlg)
    # tmp = pp.DictAlg(dictAlg).by_dim()

    if len(tmp) != 1 and dimension is None:
        raise ValueError('We never integrate over dimension.')
    if dimension is not None:
        if dimension not in tmp.keys():
            raise ValueError('dimension %d not in dictAlg dimensions %s'
                             % (dimension, str(tmp.keys())))
        tmp = {dimension: tmp[dimension]}
    dim = tmp.keys()[0]
    divisor = dim if divide_by_dimension else 1

    algorithms_with_data = [a for a in dictAlg.keys() if dictAlg[a] != []]

    dictFunc = pp.dictAlgByFun(dictAlg)

    dictData = {} # list of (ert per function) per algorithm
    dictMaxEvals = {} # list of (maxevals per function) per algorithm
    bestERT = [] # best ert per function
    # funcsolved = [set()] * len(targets) # number of functions solved per target
    for f, dictAlgperFunc in dictFunc.iteritems():
        # print target_values((f, dim))
        for j, t in enumerate(target_values((f, dim))):
        # for j, t in enumerate(genericsettings.current_testbed.ecdf_target_values(1e2, f)):
            # funcsolved[j].add(f)

            for alg in algorithms_with_data:
                x = [np.inf] * perfprofsamplesize
                runlengthunsucc = []
                try:
                    entry = dictAlgperFunc[alg][0] # one element per fun and per dim.
                    evals = entry.detEvals([t])[0]
                    assert entry.dim == dim
                    runlengthsucc = evals[np.isnan(evals) == False] / divisor
                    runlengthunsucc = entry.maxevals[np.isnan(evals)] / divisor
                    if len(runlengthsucc) > 0:
                        x = toolsstats.drawSP(runlengthsucc, runlengthunsucc,
                                             percentiles=[50],
                                             samplesize=perfprofsamplesize)[1]
                except (KeyError, IndexError):
                    #set_trace()
                    warntxt = ('Data for algorithm %s on function %d in %d-D '
                           % (alg, f, dim)
                           + 'are missing.\n')
                    warnings.warn(warntxt)

                dictData.setdefault(alg, []).extend(x)
                dictMaxEvals.setdefault(alg, []).extend(runlengthunsucc)
                
    if order is None:
        order = dictData.keys()

    # Display data
    lines = []

    for i, alg in enumerate(order):
        try:
            data = dictData[alg]
            maxevals = dictMaxEvals[alg]
        except KeyError:
            continue

        args = styles[(i) % len(styles)]
        args['linewidth'] = 1.5
        args['markersize'] = 12.
        args['markeredgewidth'] = 1.5
        args['markerfacecolor'] = 'None'
        args['markeredgecolor'] = args['color']
        args['label'] = alg
        #args['markevery'] = perfprofsamplesize # option available in latest version of matplotlib
        #elif len(show_algorithms) > 0:
            #args['color'] = 'wheat'
            #args['ls'] = '-'
            #args['zorder'] = -1
        lines.append(plotdata(np.array(data), x_limit, maxevals,
                                **args))

    labels, handles = plotLegend(lines, x_limit)
    if True:  # isLateXLeg:
        fileName = os.path.join(outputdir,'pprldmany_%s.tex' % (info))
        try:
            f = open(fileName, 'w')
            f.write(r'\providecommand{\nperfprof}{7}')
            algtocommand = {}
            for i, alg in enumerate(order):
                tmp = r'\alg%sperfprof' % pptex.numtotext(i)
                f.write(r'\providecommand{%s}{\StrLeft{%s}{\nperfprof}}' % (tmp, toolsdivers.str_to_latex(toolsdivers.strip_pathname2(str(alg)))))
                algtocommand[str(alg)] = tmp
            commandnames = []

            for l in labels:
                commandnames.append(algtocommand[l])
            # f.write(headleg)
            f.write(r'\providecommand{\perfprofsidepanel}{\mbox{%s}' % commandnames[0]) # TODO: check len(labels) > 0
            for i in range(1, len(labels)):
                f.write('\n' + r'\vfill \mbox{%s}' % commandnames[i])
            f.write('}\n')
            # f.write(footleg)
            if verbose:
                print 'Wrote right-hand legend in %s' % fileName
        except:
            raise # TODO: Does this make sense?
        else:
            f.close()

    figureName = os.path.join(outputdir,'pprldmany_%s' % (info))
    #beautify(figureName, funcsolved, x_limit*x_annote_factor, False, fileFormat=figformat)
    beautify()

    text = '%s' % (','.join(dictFunc.keys()))
    text += ',%d-D' % dim  # TODO: this is strange when different dimensions are plotted
    plt.text(0.01, 0.98, text, horizontalalignment="left",
             verticalalignment="top", transform=plt.gca().transAxes)
    if len(dictFunc) == 1:
        plt.title(str(dictFunc.keys()[0]))
    a = plt.gca()

    plt.xlim(xmin=1e-0, xmax=x_limit**annotation_space_end_relative)
    xticks, labels = plt.xticks()
    tmp = []
    for i in xticks:
        tmp.append('%d' % round(np.log10(i)))
    a.set_xticklabels(tmp)

    # print 'in_a_hurry ==', genericsettings.in_a_hurry
    if 1 < 3:
        ppfig.saveFigure(figureName, verbose=verbose)
        plt.close()

    # TODO: should return status or sthg

if __name__ == "__main__":
    # should become a test case
    import sys
    import PyGMO.util.coco.postproc
    sys.path.append('.')
    
