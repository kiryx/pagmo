ADD_LIBRARY(_util SHARED util.cpp)
PYTHON_CPPMODULE(_util)
INSTALL(FILES __init__.py DESTINATION ${PYGMO_INSTALL_PATH}/util)
INSTALL(FILES _analysis.py DESTINATION ${PYGMO_INSTALL_PATH}/util)
INSTALL(FILES _tsp.py DESTINATION ${PYGMO_INSTALL_PATH}/util)
INSTALL(DIRECTORY _bbob_pproc/ DESTINATION ${PYGMO_INSTALL_PATH}/util/_bbob_pproc FILES_MATCHING PATTERN "*.py")
# consider adding *.py for all python files in directory, 
# also add it to __all__ in __init__.py to be able to import it with from PyGMO.util import *

ADD_SUBDIRECTORY(coco)