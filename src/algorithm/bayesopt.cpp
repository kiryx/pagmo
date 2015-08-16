/*****************************************************************************
 *   Copyright (C) 2015 The PaGMO development team,                          *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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

#include "bayesopt.h"
#include "./bayesopt_cpp_wrapper/bayesoptproblem_pagmo.h"

namespace pagmo { namespace algorithm {

/// Constructor
 /**
 * Constructs a BAYESOPT algorithm
 */
bayesopt::bayesopt(const int n_iterations, const int n_inner_iterations, const int n_init_samples, 
    const int n_iter_relearn, const int init_method, const int verbose_level)
{
    m_params = initialize_parameters_to_default();
    m_params.n_iterations = n_iterations;
    m_params.n_inner_iterations = n_inner_iterations;
    m_params.n_init_samples = n_init_samples;
    m_params.n_iter_relearn = n_iter_relearn;
    m_params.init_method = init_method;
    m_params.verbose_level = verbose_level;
    m_params.random_seed = m_urng();
}

/// Evolve implementation.
/**
 * Run the BAYESOPT algorithm.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void bayesopt::evolve(pagmo::population& pop) const
{
    if (pop.size() == 0)
    {
        return;
    }

    const auto& prob = pop.problem();

    if (prob.get_i_dimension() != 0)
    {
        pagmo_throw(value_error,
                    "The problem has an integer part and BAYESOPT is not suitable to solve it.");
    }

    if (prob.get_f_dimension() != 1)
    {
        pagmo_throw(value_error,
                    "The problem is not single objective and BAYESOPT is not suitable to solve it");
    }

    // Initialization of variables
    const auto best_idx = pop.get_best_idx();
    decision_vector x = pop.get_individual(best_idx).cur_x;

    bayesoptproblem_pagmo optimizer(m_params, prob);
    boost::numeric::ublas::vector<double> bestPoint(prob.get_dimension());

    //optimize
    optimizer.optimize(bestPoint);

    std::copy(bestPoint.begin(), bestPoint.end(), x.begin());

    pop.set_x(best_idx, x);
}

void bayesopt::set_kernel(std::string name) const
{
    strcpy(m_params.kernel.name, name.c_str());
}

void bayesopt::set_mean(std::string name) const
{
    strcpy(m_params.mean.name, name.c_str());
}

void bayesopt::set_criteria(std::string name) const
{
    strcpy(m_params.crit_name, name.c_str());
}

void bayesopt::set_surrogate(std::string name) const
{
    strcpy(m_params.surr_name, name.c_str());
}

void bayesopt::set_log_file(std::string name) const
{
    strcpy(m_params.log_filename, name.c_str());
}

void bayesopt::set_load_file(std::string name) const
{
    strcpy(m_params.load_filename, name.c_str());
}

void bayesopt::set_save_file(std::string name) const
{
    strcpy(m_params.save_filename, name.c_str());
}

void bayesopt::set_learning(std::string name) const
{
    m_params.l_type = str2learn(name.c_str());
}

void bayesopt::set_score(std::string name) const
{
    m_params.sc_type = str2score(name.c_str());
}

/// Clone method.
base_ptr bayesopt::clone() const
{
    return pagmo::algorithm::base_ptr(new bayesopt(*this));
}

/// Algorithm name
std::string bayesopt::get_name() const
{
    return "BAYESOPT";
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string bayesopt::human_readable_extra() const
{
    std::ostringstream s;
    s << "n_iterations:" << m_params.n_iterations << " n_inner_iterations:" << m_params.n_inner_iterations;
    s << " n_init_samples:"<< m_params.n_init_samples << " n_iter_relearn:" << m_params.n_iter_relearn;
    s << " init_method:" << m_params.init_method << std::endl;
    return s.str();
}

}} // namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::bayesopt)
