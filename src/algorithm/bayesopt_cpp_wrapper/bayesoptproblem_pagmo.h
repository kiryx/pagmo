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

#ifndef BAYESOPTPROBLEM_PAGMO_H
#define BAYESOPTPROBLEM_PAGMO_H

#include "../../serialization.h"
#include "../../types.h"
#include "../../problem/base.h"

#include "bayesopt.hpp"

class bayesoptproblem_pagmo : public bayesopt::ContinuousModel
{
    public:
        bayesoptproblem_pagmo(bopt_params, const pagmo::problem::base &);
        double evaluateSample(const boost::numeric::ublas::vector<double> &);
        bool checkReachability(const boost::numeric::ublas::vector<double> &);
        const pagmo::problem::base_ptr m_prob;


    private:
        //problem
        //dimensions
        pagmo::decision_vector::size_type m_dim;
        // Structure containing all the BAYESOPT parameters, see BAYESOPT documentation
        bopt_params m_params;
};

#endif // BAYESOPTPROBLEM_PAGMO_H
