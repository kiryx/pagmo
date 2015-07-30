ADD_LIBRARY(_benchmark SHARED benchmark.cpp)
PYTHON_CPPMODULE_ADDITION(_benchmark /util/coco)
INSTALL(FILES __init__.py DESTINATION ${PYGMO_INSTALL_PATH}/util/coco)
# consider adding *.py for all python files in directory, 
# also add it to __all__ in __init__.py to be able to import it with from PyGMO.util import *
