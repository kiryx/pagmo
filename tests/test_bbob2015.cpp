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
#include "../src/pagmo.h"
#include "test.h"

using namespace pagmo;

//Test BBOB2015 benchmark functions
int test_bbob2015(int dimension)
{

	//Test non-noisy benchmark functions
	for(unsigned int i=1; i<25; i++)
	{
		std::cout<<"Problem Number:"<<i<<std::endl;
		problem::bbob2015 prob(i, dimension);
		prob.set_seed(1);
		
		decision_vector x(dimension);

		for(int j = 0; j<dimension; j++)
		{
			x[j] = 2;
		}

		//get fitness
		fitness_vector fitness = prob.objfun(x);
		
		std::cout<<fitness[0];
		
		std::cout<<std::endl;
	}

	//Test noisy benchmark functions
	for(unsigned int i=101; i<131; i++)
	{
		std::cout<<"Problem Number:"<<i<<std::endl;
		problem::bbob2015 prob(i, dimension);
		prob.set_seed(1);
		
		decision_vector x(dimension);

		for(int j = 0; j<dimension; j++)
		{
			x[j] = 2;
		}

		//get fitness
		fitness_vector fitness = prob.objfun(x);
		
		std::cout<<fitness[0];
		
		std::cout<<std::endl;
	}

	return 0;

}

int main()
{

	int dimension = 5;
	return test_bbob2015(dimension);
}
