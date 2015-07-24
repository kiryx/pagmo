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
	//We instantiate the Ackley problem with dimension 8
	pagmo::problem::ackley prob(8);
	for(unsigned int instance=1; instance <= 100; instance++)
	{
		//We instantiate the algorithm PSO with 5000 generation
		pagmo::algorithm::pso algo(5000);

		pagmo::util::bbob benchmarking(prob, "./", "PSO", instance);
		//1 - Evolution takes place on the same thread as main
		//We instantiate a population containing 50 candidate solutions to the Ackley problem
		pagmo::population pop(benchmarking, 50);
		algo.evolve(pop);

		benchmarking.finalize(pop);
	}
	return 0;
}
