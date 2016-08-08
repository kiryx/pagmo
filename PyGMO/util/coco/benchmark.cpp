/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#include <Python.h>
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>

#include "../../../src/util/benchmark.h"
#include "../../utils.h"

using namespace boost::python;
using namespace pagmo;
using namespace pagmo::util::coco;

BOOST_PYTHON_MODULE(_benchmark) {

	//COCO meta-problem for benchmarking.
	class_<util::coco::benchmark, bases<problem::base> >("benchmark","COCO meta-problem for using a pagmo::problem as benchmark", init<const util::coco::benchmark &>())
		.def(init<>())
		.def("__copy__", &Py_copy_from_ctor<util::coco::benchmark>)
		.def("__deepcopy__", &Py_deepcopy_from_ctor<util::coco::benchmark>)
		.def_pickle(python_class_pickle_suite<util::coco::benchmark>())
		.def("cpp_loads", &py_cpp_loads<util::coco::benchmark>)
		.def("cpp_dumps", &py_cpp_dumps<util::coco::benchmark>)
		.def("finalize", &util::coco::benchmark::finalize)
		.def(init<pagmo::problem::base &, std::string, std::string, unsigned int, std::string>());

}
