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
#include <algorithm>
#include <cmath>
#include <math.h>
#include <iostream>

#include "../rng.h"
#include "cuckoosearch.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen number of generations
 * @param[in] pa Discovery rate of alien eggs/solutions
 */
cuckoosearch::cuckoosearch(int gen, double pa):base(), m_gen(gen), m_pa(pa)
{
    if (gen < 0)
        pagmo_throw(value_error,"number of generations must be nonnegative");

    if(pa <= 0 || pa >= 1)
        pagmo_throw(value_error,"pulse rate must be in (0, 1)");
}

/// Clone method.
base_ptr cuckoosearch::clone() const
{
    return base_ptr(new cuckoosearch(*this));
}


/// Evolve implementation.
/**
 * Run the Cuckoo Search algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void cuckoosearch::evolve(population &pop) const
{
    // Let's store some useful variables.
    const problem::base &prob = pop.problem();
    const problem::base::size_type D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(),
        prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();

    const problem::base::size_type Dc = D - prob_i_dimension;
    const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
    const population::size_type swarm_size = pop.size();

    //We perform some checks to determine wether the problem/population are suitable for Cuckoo Search
    if(Dc == 0)
    {
        pagmo_throw(value_error,"There is no continuous part in the problem decision vector for Cuckoo Search to optimise");
    }

    if(prob_c_dimension != 0)
    {
        pagmo_throw(value_error,"The problem is not box constrained and Cuckoo Search is not suitable to solve it");
    }

    if(prob_f_dimension != 1)
    {
        pagmo_throw(value_error,"The problem is not single objective and Cuckoo Search is not suitable to solve it");
    }

    // Get out if there is nothing to do.
    if (swarm_size == 0 || m_gen == 0)
        return;

    // Some vectors used during evolution are allocated here.
    std::vector<double> dummy(Dc, 0);    // used for initialisation purposes

    std::vector<decision_vector> nest(swarm_size, dummy), newnest(swarm_size, dummy);  // nests/solutions
    std::vector<fitness_vector> fitness(swarm_size);    //corresponding fitness values

    population::size_type p;    // for iterating over cuckoosearchs
    problem::base::size_type d; // for iterating over problem dimensions

    decision_vector temp_x, gbest_x;
    fitness_vector temp_fit, gbest_fit;

    // Copy the cuckoosearch positions, their velocities and their fitness
    for(p = 0; p < swarm_size; p++)
    {
        nest[p] = pop.get_individual(p).cur_x;
        fitness[p] = pop.get_individual(p).cur_f;
    }

    gbest_fit = pop.champion().f;
    gbest_x = pop.champion().x;

    temp_x = gbest_x;
    temp_fit = gbest_fit;

    /* --- Main Cuckoo Search loop ---*/
    // For each iteration
    for(int g = 0; g < m_gen; g++)
    {
        //Generate new solutions (but keep the current best)
        get_cuckoos(prob, newnest, nest, gbest_x);
        get_best_nest(prob, temp_fit, temp_x, nest, newnest, fitness);

        //Discovery and randomization
        empty_nests(prob, newnest, nest);

        //Evaluate this set of solutions
        get_best_nest(prob, temp_fit, temp_x, nest, newnest, fitness);

        if(prob.compare_fitness(temp_fit, gbest_fit))
        {
            gbest_x = temp_x;
            gbest_fit = temp_fit;
        }
    } // end of main Cuckoo Search loop

    // copy nests/solutions back to the main population
    for(p = 0; p < swarm_size; p++)
    {
        pop.set_x(p, nest[p]); // sets: cur_x, cur_f, best_x, best_f
    }
}

///Levy flights
///Levy exponent and coefficient
///For details, see equation (2.21), Page 16 (chapter 2) of the book
///X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
void cuckoosearch::get_cuckoos(const problem::base &prob, std::vector<decision_vector> &newnest,
    const std::vector<decision_vector> &nest, const decision_vector best) const
{
    rng_uint32 rng = rng_generator::get<rng_uint32>();
    boost::normal_distribution<> nd(0.0, 1.0);

    boost::variate_generator<rng_uint32&, boost::normal_distribution<> > ndrng(rng, nd);
    double beta = 3.0/2.0;

    //sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta)
    double sigma = pow((tgamma(1+beta)*sin(M_PI*beta/2)/(tgamma((1+beta)/2)*beta*pow(2, (beta-1)/2))), (1/beta));
    std::vector<decision_vector>::size_type n;

    decision_vector::size_type dim = nest[0].size();

    decision_vector step(dim, 0);

    //In the next equation, the difference factor (s-best) means that
    //when the solution is the best solution, it remains unchanged.
    //Here the factor 0.01 comes from the fact that L/100 should the typical
    //step size of walks/flights where L is the typical lenghtscale;
    //otherwise, Levy flights may become too aggresive/efficient,
    //which makes new solutions (even) jump out side of the design domain
    //(and thus wasting evaluations).
    //Now the actual random walks or flights
    for(n=0; n<nest.size(); n++)
    {
        for(decision_vector::size_type i=0; i<dim; i++)
        {
            step[i] = ndrng()*sigma/pow(std::abs(ndrng()), 1/beta);
            newnest[n][i] = nest[n][i] + 0.01*step[i]*(nest[n][i] - best[i])*ndrng();
        }
    }
    for(n=0; n<newnest.size(); n++)
        applysimplebounds(prob, newnest[n]);
}

/// Find the current best nest
void cuckoosearch::get_best_nest(const problem::base &prob, fitness_vector &best_f, decision_vector &best_x,
    std::vector<decision_vector> &nests, const std::vector<decision_vector> &newnest,
    std::vector<fitness_vector> &fitness) const
{
    //Evaluating all new solutions
    if(nests.size() != newnest.size())
        pagmo_throw(value_error, "Dimensions of old nests and new nests should be same");

    else if(nests.size() != fitness.size())
        pagmo_throw(value_error, "Dimensions of nests and vector of fitness vectors should be same");

    fitness_vector fnew = best_f;

    std::vector<decision_vector>::size_type d;
    for(d = 0; d < nests.size(); d++)
    {
        prob.objfun(fnew, newnest[d]);

        if(prob.compare_fitness(fnew, fitness[d]))
        {
            fitness[d] = fnew;
            nests[d] = newnest[d];
        }
    }

    //find best and return
    best_f = fitness[0];
    best_x = nests[0];
    for(d = 0; d < nests.size(); d++)
    {
        if(prob.compare_fitness(fitness[d], best_f))
        {
            best_f = fitness[d];
            best_x = nests[d];
        }
    }
}

///Replace some nests by constructing new solutions/nests
///A fraction of worse nests are discovered with a probability pa
void cuckoosearch::empty_nests(const problem::base &prob, std::vector<decision_vector> &newnest,
    const std::vector<decision_vector> &nest) const
{
    //In the real world, if a cuckoo's egg is very similar to a host's eggs, then
    //this cuckoo's egg is less likely to be discovered, thus the fitness should
    //be related to the difference in solutions.  Therefore, it is a good idea
    //to do a random walk in a biased way with some random step sizes.
    //New solution by biased/selective random walks

    std::vector<decision_vector>::size_type n=nest.size(), i;

    newnest = nest;

    decision_vector::size_type dim;

    dim = nest[0].size();

    decision_vector::size_type d;

    //Discovered or not -- a status vector
    std::vector<std::vector<int> > K(n, std::vector<int>(dim, 0));
    for(std::vector<std::vector<int> >::size_type j=0; j < K.size(); j++)
    {
        for(d=0; d<dim; d++)
            K[j][d] = m_drng() > m_pa;
    }

    std::vector<std::vector<decision_vector>::size_type> indexes1, indexes2;
    indexes1.reserve(nest.size());
    indexes2.reserve(nest.size());

    for(std::vector<decision_vector>::size_type t=0; t<n; t++)
    {
        indexes1[t] = t;
        indexes2[t] = t;
    }

    //Generate random permutation of indexes
    for(i=n-1; i>0; --i)
    {
        std::swap(indexes1[i], indexes1[m_urng()%i]);
        std::swap(indexes2[i], indexes2[m_urng()%i]);
    }

    double r = m_drng();

    //New solution by biased/selective random walks
    for(i=0; i<n; i++)
    {
        for(d=0; d<dim; d++)
        {
            newnest[i][d] = nest[i][d] + ((r * (nest[indexes1[i]][d] - nest[indexes2[i]][d])) * K[i][d]);
        }
    }

    for(i=0; i<n; i++)
        applysimplebounds(prob, newnest[i]);
}

///Apply simple bounds
void cuckoosearch::applysimplebounds(const problem::base &prob, decision_vector &x) const
{
    const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();

    decision_vector::size_type d, len;

    len = lb.size();

    for(d = 0; d < len; d++)
    {
        if(x[d] < lb[d])
        {
            x[d] = lb[d];
        }
        else if(x[d] > ub[d])
        {
            x[d] = ub[d];
        }
    }
}

/// Algorithm name
std::string cuckoosearch::get_name() const
{
    return "Cuckoo Search algorithm";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string cuckoosearch::human_readable_extra() const
{
    std::ostringstream s;
    s << "gen:" << m_gen;
    s << "pa:" << m_pa;

    return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::cuckoosearch)
