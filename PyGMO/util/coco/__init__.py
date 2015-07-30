# -*- coding: utf-8 -*-
from PyGMO.util.coco._benchmark import *

__all__ = ['benchmark']

def _benchmark_ctor(self, problem, datapath="./", algname="", instanceId=1, comments = ""):
    """
    Constructs a meta-problem for COCO benchmarking.
    USAGE: benchmark(problem, datapath, algname, instanceId)
    
    * problem: The problem to be used as benchmark
    * datapath: The data path where log files will be stored
    * algname: The name of the optimizer used like PSO
    * instanceId: Id of the instance
    """
    arg_list = []
    arg_list.append(problem)
    arg_list.append(datapath)
    arg_list.append(algname)
    arg_list.append(instanceId)
    arg_list.append(comments)
    self._orig_init(*arg_list)

benchmark._orig_init = benchmark.__init__
benchmark.__init__ = _benchmark_ctor
