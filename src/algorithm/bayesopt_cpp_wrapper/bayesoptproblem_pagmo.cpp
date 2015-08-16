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

#include "bayesoptproblem_pagmo.h"


bayesoptproblem_pagmo::bayesoptproblem_pagmo(bopt_params param, const pagmo::problem::base &prob):
    bayesopt::ContinuousModel(prob.get_dimension(), param),
    m_prob(prob.clone()),
    m_dim(prob.get_dimension()),
    m_params(param){}

double bayesoptproblem_pagmo::evaluateSample(const boost::numeric::ublas::vector<double> &query)
{
    pagmo::fitness_vector f(1,0);
    pagmo::decision_vector d(m_dim, 0);
    std::copy(query.begin(), query.end(), d.begin());

    m_prob->objfun(f, d);

    return f[0];
}

bool bayesoptproblem_pagmo::checkReachability(const boost::numeric::ublas::vector<double> &query)
{
    pagmo::decision_vector d(m_dim, 0);
    std::copy(query.begin(), query.end(), d.begin());

    return m_prob->feasibility_x(d);
}
