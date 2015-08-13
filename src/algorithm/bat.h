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

#ifndef PAGMO_ALGORITHM_BAT_H
#define PAGMO_ALGORITHM_BAT_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"


namespace pagmo { namespace algorithm {

/// Bat Algorithm
/**
 *
 * Bat Algorithm (BA) is a population based algorithm that has been proposed in the year 2010.
 * Bat algorithm is based on the echolocation behaviour of microbats with varying pulse rates of
 * emission and loudness. In BA each bat has each bat is associated with velocity \f$ \mathbf v_i^t \f$ and
 * location \f$ \mathbf x_i^t \f$
 * memory of the position where it achieved the best performance \f$\mathbf x^l_i\f$ and of the swarm
 * 'champion' position \f$ \mathbf x^g \f$ and uses this information to update its parameters using the equations:
 * \f[
 * \mathbf f_{i} = \mathbf f_{min} + \left( \mathbf f_{max} - \mathbf f_{min} \right) \cdot \mathbf \beta
 * \f]
 * \f[
 * \mathbf v_i^t = \mathbf v_i^{t-1} + \left( \mathbf x_i^{t-1} - \mathbf x^g \right) \cdot \mathbf f_i
 * \f]
 * \f[
 * \mathbf x_i^t = \mathbf x_i^{t-1} + \mathbf v_i^t
 * \f]
 * where \mathbf \beta is a uniformly drawn random number in range [0, 1]
 * The algorithm is suitable for box-constrained single-objective continuous optimization.
 *
 * @see http://arxiv.org/pdf/1004.4170v1.pdf for the paper on this algorithm
 *
 * @author Sunil K. Mahendrakar (sunil.mahendrakar19@gmail.com)
 */

class __PAGMO_VISIBLE bat: public base
{
public:
    bat(int gen=1, double qmax=2, double qmin=0, double alpha=0.9, double gam=0.9, double loudness=0.5, double pulserate=0.5);
    base_ptr clone() const;
    void evolve(population &) const;
    std::string get_name() const;

protected:
    std::string human_readable_extra() const;
private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar & boost::serialization::base_object<base>(*this);
        ar & const_cast<int &>(m_gen);
        ar & const_cast<double &>(m_alpha);
        ar & const_cast<double &>(m_gamma);
        ar & const_cast<double &>(m_qmax);
        ar & const_cast<double &>(m_qmin);
        ar & const_cast<double &>(m_loudness);
        ar & const_cast<double &>(m_pulserate);
    }
    // Number of generations
    const int m_gen;
    // Maximum frequency
    const double m_qmax;
    // Minimum frequency
    const double m_qmin;
    // initial pulserate
    const double m_pulserate;
    // initial loudness
    const double m_loudness;
    // rate of decrease in loadness
    const double m_alpha;
    // rate of increase in pulserate
    const double m_gamma;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::bat)

#endif // BAT_H
