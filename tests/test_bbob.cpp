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
 *   the Free Software Foundation; either version 2 of the License, or       *
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

#include <iostream>
#include <iomanip>
#include "../src/pagmo.h"

using namespace pagmo;

//Test BBOB benchmarking.
int main()
{
	//We instantiate the problem Schwefel with dimension 10
	pagmo::problem::schwefel prob(10);
	std::vector<decision_vector> bestx;
	decision_vector x(10, 418.9828872724338);
	bestx.push_back(x);
	prob.set_best_x(bestx);
	//We instantiate the algorithm PSO with 10 generation
	pagmo::algorithm::pso algo(10);

	pagmo::util::bbob benchmarking(prob, "./", "PSO");
	//1 - Evolution takes place on the same thread as main
	//We instantiate a population containing 20 candidate solutions to the Schwefel problem
	pagmo::population pop(benchmarking,20, 1);
	algo.evolve(pop);

	benchmarking.finalize(pop);

	return 0;
}
