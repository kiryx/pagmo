/*****************************************************************************
*   Copyright (C) 2004-2009 The PaGMO development team,                     *
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

#ifndef PAGMO_PROBLEM_PL2PL_FIXED_TIME_H
#define PAGMO_PROBLEM_PL2PL_FIXED_TIME_H

#include <vector>
#include <string>
#include <keplerian_toolbox/epoch.h>
#include <keplerian_toolbox/planet/jpl_low_precision.h>
#include <keplerian_toolbox/sims_flanagan/leg.h>
#include <keplerian_toolbox/sims_flanagan/spacecraft.h>
#include <keplerian_toolbox/sims_flanagan/throttle.h>


#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"


namespace pagmo { namespace problem {


/// GTOC_2 Low-Thrust Multiple Asteroid Randezvous Problem
/**
 * This is the problem given by Jet Propulsion Laboratories as the 2nd Global Trajectory
 * Optimization Competition. It is transcribed assembling 4 Sims-Flanagan trajectories legs
 *
 * @see http://www.esa.int/gsp/ACT/inf/op/globopt/evevejsa.htm
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE pl2pl_fixed_time: public base
{
	public:
		/// Constructor
		pl2pl_fixed_time(const kep_toolbox::planet::planet_ptr = kep_toolbox::planet::jpl_lp("earth").clone(),const kep_toolbox::planet::planet_ptr = kep_toolbox::planet::jpl_lp("mars").clone(),const kep_toolbox::epoch = kep_toolbox::epoch(0),const kep_toolbox::epoch = kep_toolbox::epoch(1000),const
 kep_toolbox::sims_flanagan::spacecraft = kep_toolbox::sims_flanagan::spacecraft(2000,0.3,3000),const int = 5);
		base_ptr clone() const;
		std::string get_name() const;
		std::string pretty(const std::vector<double> &x) const;
		bool get_high_fidelity() const;
		void set_high_fidelity(bool hf_option);
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
	private:
		template <class Iterator>
		kep_toolbox::sims_flanagan::throttle get_nth_throttle(int n, Iterator it, const kep_toolbox::epoch &start, const kep_toolbox::epoch &end) const
		{
			Iterator n_it = it + 3 * n;
			kep_toolbox::array3D tmp = {{ *n_it, *(n_it + 1), *(n_it + 2)}};
			const double seg_duration = (end.mjd() - start.mjd()) / m_n_seg;
			return kep_toolbox::sims_flanagan::throttle( kep_toolbox::epoch(start.mjd() + seg_duration * n,kep_toolbox::epoch::MJD),
				kep_toolbox::epoch(start.mjd() + seg_duration * (n + 1),kep_toolbox::epoch::MJD),
				tmp);
		}
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			(void)ar;
			/*ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<int &>(m_n_seg);
			ar & m_asteroids;
			ar & m_leg;
			ar & m_t0;
			ar & m_t1;
			ar & const_cast< kep_toolbox::sims_flanagan::spacecraft &>(m_spacecraft);*/
		}
		const int						m_n_seg;
		const kep_toolbox::planet::planet_ptr		        m_ast0;
		const kep_toolbox::planet::planet_ptr		        m_ast1;
		mutable kep_toolbox::sims_flanagan::leg			m_leg;
		const kep_toolbox::epoch				m_t0;
		const kep_toolbox::epoch				m_t1;		
		const kep_toolbox::sims_flanagan::spacecraft		m_spacecraft;
};

} } // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::pl2pl_fixed_time)

#endif

