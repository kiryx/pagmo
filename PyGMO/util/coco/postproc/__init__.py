#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""COmparing Continuous Optimisers (COCO) post-processing software

This package is meant to generate output figures and tables for the
benchmarking of continuous optimisers in the case of black-box
optimisation.
The post-processing tool takes as input data from experiments and
generates outputs that will be used in the generation of the LateX-
formatted article summarizing the experiments.

The main method of this package is :py:func:`PyGMO.util.coco.postproc.rungeneric.main`
This method allows to use the post-processing through a command-line
interface.

To obtain more information on the use of this package from the python
interpreter, assuming this package has been imported as ``bb``, type:
``help(bb.cococommands)``

"""

from __future__ import absolute_import

import sys

from PyGMO.util.coco.postproc.cococommands import *
from PyGMO.util.coco.postproc import pprldistr, ppfigdim, ppsingle, ppfigparam

__all__  = ['pprldistr', 'pproc', 'ppfig', 'ppfigdim', 'ppsingle', 'ppfigparam']

__version__ = '15.00'
