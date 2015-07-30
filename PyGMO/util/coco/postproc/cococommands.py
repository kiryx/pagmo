#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module for using COCO from the (i)Python interpreter.

For all operations in the Python interpreter, it will be assumed that
the package has been imported as bb, just like it is done in the first
line of the examples below.

The main data structures used in COCO are :py:class:`DataSet`, which
corresponds to data of one algorithm on one problem, and
:py:class:`DataSetList`, which is for collections of :py:class:`DataSet`
instances. Both classes are implemented in :py:mod:`PyGMO.util.coco.postproc.pproc`.

Examples:

* Start by importing :py:mod:`PyGMO.util.coco.postproc`::

    >>> import PyGMO.util.coco.postproc as bb # load PyGMO.util.coco.postproc
    >>> import os
    >>> os.chdir(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

* Load a data set, assign to variable :py:data:`ds`::

      >>> ds = bb.load('./test') # Test data stored in ./test folder
      Processing ./test

* Get some information on a :py:class:`DataSetList` instance::

      >>> print ds #
      >>> bb.info(ds)
      1 data set
      Algorithm: PSO
      Dimension: 10
      Function: 1
      Max evals: 75017
      Df      |     min       10      med       90      max
      --------|--------------------------------------------
      1.0e+01 |      55      151     2182    49207    55065
      1.0e+00 |     124      396     2820    56879    59765
      1.0e-01 |     309      466     2972    61036    66182
      1.0e-03 |     386      519     3401    67530    72091
      1.0e-05 |     446      601     3685    70739    73472
      1.0e-08 |     538      688     4052    72540    75010

"""

from __future__ import absolute_import

#from PyGMO.util.coco.postproc import ppsingle, ppfigdim, dataoutput
# from PyGMO.util.coco.postproc.pproc import DataSetList, DataSet
from PyGMO.util.coco.postproc import pproc

#__all__ = ['load', 'info', 'pickle', 'systeminfo', 'DataSetList', 'DataSet']

def load(filename):
    """Create a :py:class:`DataSetList` instance from a file or folder.
    
    Input argument filename can be a single :file:`info` file name, a
    single pickle filename or a folder name. In the latter case, the
    folder is browsed recursively for :file:`info` or :file:`pickle`
    files.

    """
    return pproc.DataSetList(filename)

# info on the DataSetList: algId, function, dim

def info(dsList):
    """Display more info on an instance of DatasetList."""
    dsList.info()

# TODO: method for pickling data in the current folder!
def pickle(dsList):
    """Pickle a DataSetList."""
    dsList.pickle(verbose=True)
    # TODO this will create a folder with suffix -pickle from anywhere:
    # make sure the output folder is created at the right location

def systeminfo():
    """Display information on the system."""
    import sys
    print sys.version
    import numpy
    print 'Numpy %s' % numpy.__version__
    import matplotlib
    print 'Matplotlib %s' % matplotlib.__version__
    import PyGMO.util.coco.postproc
    print 'PyGMO.util.coco.postproc %s' % PyGMO.util.coco.postproc.__version__
