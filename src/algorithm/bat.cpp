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

#include <boost/random/uniform_real.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <vector>
#include <cmath>
#include <math.h>
#include <iostream>

#include "../rng.h"
#include "bat.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen number of generations
 * @param[in] qmax maximum frequency of bats must be nonnegative
 * @param[in] qmin minimum frequency of bats must be less than qmax
 * @param[in] alpha rate of decrease in loudness (1 for constant loudness)
 * @param[in] gam rate of increase in pulse rate should be greater than 0
 * @param[in] loudness initial loudness of all bats.
 * @param[in] pulserate initial pulse rate of all bats.
 */
bat::bat(int gen, double qmax, double qmin, double alpha, double gam, double loudness, double pulserate):base(),
		m_gen(gen), m_qmax(qmax), m_qmin(qmin), m_pulserate(pulserate),
		m_loudness(loudness), m_alpha(alpha), m_gamma(gam)
{
	if (gen < 0)
		pagmo_throw(value_error,"number of generations must be nonnegative");

	if(pulserate < 0 || pulserate > 1)
		pagmo_throw(value_error,"pulse rate must be in [0, 1]");

	if(loudness <= 0)
		pagmo_throw(value_error,"initial loudness must be greater than 0");

	if(qmin < 0)
		pagmo_throw(value_error,"minimum frequency must be nonnegative");

	if(qmax < 0)
		pagmo_throw(value_error,"minimum frequency must be nonnegative");

	if(qmin > qmax)
		pagmo_throw(value_error,"maximum frequency must be greater than minimum frequency");

	if(alpha <= 0 || alpha > 1)
		pagmo_throw(value_error,"alpha must be in (0, 1]");

	if(gam <= 0)
		pagmo_throw(value_error,"gamma must be greater than 0");
}

/// Clone method.
base_ptr bat::clone() const
{
	return base_ptr(new bat(*this));
}


/// Evolve implementation.
/**
 * Run the BA algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void bat::evolve(population &pop) const
{
	rng_uint32 rng = rng_generator::get<rng_uint32>();

	boost::normal_distribution<> nd(0.0, 1.0);

	boost::variate_generator<rng_uint32&, boost::normal_distribution<> > ndrng(rng, nd);

	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(),
		prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();

	const problem::base::size_type Dc = D - prob_i_dimension;
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type swarm_size = pop.size();

	//We perform some checks to determine wether the problem/population are suitable for BA
	if(Dc == 0)
	{
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for BA to optimise");
	}

	if(prob_c_dimension != 0)
	{
		pagmo_throw(value_error,"The problem is not box constrained and BA is not suitable to solve it");
	}

	if(prob_f_dimension != 1)
	{
		pagmo_throw(value_error,"The problem is not single objective and BA is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if (swarm_size == 0 || m_gen == 0)
		return;

	// Some vectors used during evolution are allocated here.
	std::vector<double> dummy(Dc, 0);	// used for initialisation purposes
	std::vector<double> freq(swarm_size, 0);   // bats' frequencies
	std::vector<double> pulserate(swarm_size, m_pulserate);   // bats' frequencies
	std::vector<double> loudness(swarm_size, m_loudness);   // bats' loudness

	std::vector<decision_vector> X(swarm_size, dummy);  // bats' current positions

	decision_vector Sol;  // Temporary position
	fitness_vector Fnew;	// Temporary fitness

	std::vector<fitness_vector> lbfit(swarm_size);	// bats' current fitness values

	std::vector<decision_vector> V(swarm_size, dummy);   // bats' velocities

	decision_vector gbest_x;   // bats' previous best position among all bats

	fitness_vector gbest_fit;	// fitness at the best found search space position

	decision_vector minv(Dc), maxv(Dc); // Maximum and minimum velocity allowed

	double new_x;   // Temporary variable

	population::size_type p;	// for iterating over bats
	problem::base::size_type d; // for iterating over problem dimensions

	// Copy the bat positions, their velocities and their fitness
	for(p = 0; p < swarm_size; p++)
	{
		X[p] = pop.get_individual(p).cur_x;
		V[p] = pop.get_individual(p).cur_v;
		lbfit[p] = pop.get_individual(p).cur_f;
	}

	//initialize temporary variables
	Fnew = lbfit[0];
	Sol = X[0];

	gbest_fit = pop.champion().f;
	gbest_x = pop.champion().x;

	/* --- Main BA loop ---*/
	// For each generation
	for(int g = 0; g < m_gen; ++g)
	{
		// For each bat in the swarm
		for(p = 0; p < swarm_size; p++)
		{
			/*-------Bat Algorithm---------------------------------------------*/
			freq[p] = m_qmin + ((m_qmin - m_qmax) * m_drng());

			//Calculate velocity
			for(d = 0; d < Dc; d++)
			{
				V[p][d] = V[p][d] + ((X[p][d] - gbest_x[d]) * freq[p]);
			}

			// and we perform the position update and the feasibility correction
			for(d = 0; d < Dc; d++)
			{

				// update position
				new_x = X[p][d] + V[p][d];

				// feasibility correction
				// (velocity updated to that which would have taken the previous position
				// to the newly corrected feasible position)
				if(new_x < lb[d])
				{
					new_x = lb[d];
				}
				else if(new_x > ub[d])
				{
					new_x = ub[d];
				}

				Sol[d] = new_x;
			}

			// pulse rate
			if(m_drng() > pulserate[p])
			{
				//Random walk, the factor 0.001 limits the step sizes of random walks
				for(d = 0; d < Dc; d++)
				{
					double r = ndrng();
					Sol[d] = gbest_x[d] + (0.001 * r);
				}
			}

			//feasibility correction
			for(d = 0; d < Dc; d++)
			{
				if(Sol[d] < lb[d])
				{
					Sol[d] = lb[d];
				}
				else if(Sol[d] > ub[d])
				{
					Sol[d] = ub[d];
				}
			}

			// We evaluate here the new individual fitness as to be able to update the global best in real time
			prob.objfun(Fnew, Sol);
			m_fevals++;

			// update if solution is improved, and not too loud
			if(prob.compare_fitness(Fnew, lbfit[p]) && m_drng() < loudness[p])
			{
				// update the bat's previous best position and fitness
				X[p] = Sol;
				lbfit[p] = Fnew;

				//decrease loudness
				loudness[p] = m_alpha * loudness[p];

				//increase pulse rate
				pulserate[p] = pulserate[p] * (1 - exp(-1 * m_gamma * g));
			}

			// update the best position observed so far by any bat in the swarm
			if(prob.compare_fitness(Fnew, gbest_fit))
			{
				gbest_x = Sol;
				gbest_fit = Fnew;
			}
		} //End of loop on the population members

	} // end of main BA loop

	// copy bats' positions & velocities back to the main population
	for(p = 0; p < swarm_size; p++)
	{
		pop.set_x(p, X[p]); // sets: cur_x, cur_f, best_x, best_f
		pop.set_v(p, V[p]); // sets: cur_v
	}
}

/// Algorithm name
std::string bat::get_name() const
{
	return "Bat algorithm";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string bat::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "qmax:" <<  m_qmax << ' ';
	s << "qmin:" <<  m_qmin << ' ';
	s << "pulserate:" <<  m_pulserate << ' ';
	s << "loudness:" <<  m_loudness << ' ';
	s << "alpha:" <<  m_alpha << ' ';
	s << "gamma:" <<  m_gamma << ' ';

	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::bat)
