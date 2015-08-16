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

#ifndef PAGMO_ALGORITHM_BAYESOPT_H
#define PAGMO_ALGORITHM_BAYESOPT_H

#include "./bayesopt_cpp_wrapper/bayesoptproblem_pagmo.h"
#include <string>

#include "base.h"

namespace pagmo {

// Forward declaration
class population;

namespace algorithm {

class __PAGMO_VISIBLE bayesopt : public base
{
public:
    bayesopt(const int = 190, const int = 500, const int = 10, const int = 50, const int = 1, const int = 0);
    void evolve(pagmo::population&) const;
    pagmo::algorithm::base_ptr clone() const;
    std::string get_name() const;

    //parameter settings
    //see http://rmcantin.bitbucket.org/html/usemanual.html#params
    void set_kernel(std::string) const;
    void set_mean(std::string) const;
    void set_criteria(std::string) const;
    void set_surrogate(std::string) const;
    void set_log_file(std::string) const;
    void set_load_file(std::string) const;
    void set_save_file(std::string) const;
    void set_learning(std::string) const;
    void set_score(std::string) const;

protected:
    std::string human_readable_extra() const;

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar & boost::serialization::base_object<base>(*this);
    }
    // Structure containing all the BAYESOPT parameters, see BAYESOPT documentation
    mutable bopt_params m_params;
};

}} // namespaces

#endif // PAGMO_ALGORITHM_BAYESOPT_H

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::bayesopt)
