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

#include <algorithm>
#include <functional>
#include <vector>
#include <string>
#include <sstream>
#include <keplerian_toolbox/astro_constants.h>

#include "pl2pl_fixed_time.h"
#include "../exceptions.h"
#include "../types.h"
#include "base.h"


using namespace kep_toolbox;
using namespace kep_toolbox::sims_flanagan;

namespace pagmo { namespace problem {

/// Problem Constructor
/**
 * Constructs a pl2pl_fixed_time problem
 *
 * @param[in] ast0 id of the first asteroid to visit
 * @param[in] ast1 id of the second asteroid to visit
 * @param[in] t0 epoch at departure from ast0
 * @param[in] t1 epoch at arrival at ast1
 * @param[in] S object containing specific impulse, max thrust and mass of the spacecraft 
 * @param[in] n_seg number of segments to be used for modeling the transfer
 *
 *
 * @see problem::base constructors.
 */

pl2pl_fixed_time::pl2pl_fixed_time(const planet::planet_ptr ast0,const planet::planet_ptr ast1,const epoch t0,const epoch t1,const spacecraft S,const int n_seg):
	base(3 * n_seg + 1,0,1,n_seg + 7, n_seg, 1E-5),
	m_ast0(ast0), m_ast1(ast1), m_t0(t0), m_t1(t1), m_spacecraft(S), m_n_seg(n_seg)
{
	if (n_seg <= 0) {
		pagmo_throw(value_error,"invalid number of segments");
	}
	
	array3D r, v;
	double initial_mass = m_spacecraft.get_mass();
	// Build leg.
	m_leg.set_spacecraft(m_spacecraft);
	m_leg.set_mu(ASTRO_MU_SUN);
	m_leg.set_throttles_size(n_seg);
	m_leg.set_t_i(t0);
	m_leg.set_t_f(t1);
	m_leg.set_high_fidelity(true);
	// Initial state.
	m_ast0->eph(m_t0,r,v);
	m_leg.set_x_i(sc_state(r,v,initial_mass));

	decision_vector lb_v, ub_v;
	// Mass.
	//lb_v.push_back(std::max(initial_mass-(t1-t0)*ASTRO_DAY2SEC/m_spacecraft.get_isp(),5));
	lb_v.push_back(800);
	ub_v.push_back(initial_mass);
	// Throttles.
	for (int i = 0; i < 3 * n_seg; ++i) {
		lb_v.push_back(-1);
		ub_v.push_back(1);
	}
	set_bounds(lb_v,ub_v);
}

base_ptr pl2pl_fixed_time::clone() const
{
	return base_ptr(new pl2pl_fixed_time(*this));
}

void pl2pl_fixed_time::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	f[0] = -x[0];
}

std::vector<kep_toolbox::planet::planet_ptr> pl2pl_fixed_time::get_sequence() const {
	std::vector<kep_toolbox::planet::planet_ptr> seq;
	seq.push_back(m_ast0);
	seq.push_back(m_ast1);
	return seq;
}

kep_toolbox::sims_flanagan::leg pl2pl_fixed_time::get_leg() const {
	return m_leg;
}

kep_toolbox::epoch pl2pl_fixed_time::get_t0() const {
	return m_t0;
}
kep_toolbox::epoch pl2pl_fixed_time::get_t1() const {
	return m_t1;
}

void pl2pl_fixed_time::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	// Cached values.
	array3D r, v;
	// Build leg.
	for (int j = 0; j < m_n_seg; ++j) {
		m_leg.set_throttles(j,get_nth_throttle(j,x.begin() + 1,m_t0,m_t1));
	}
	// Final state.
	double final_mass = x[0];
	m_ast1->eph(m_t1,r,v);
	m_leg.set_x_f(sc_state(r,v,final_mass));

	//Load state mismatches into constraints vector.
	m_leg.get_mismatch_con(c.begin(), c.begin() + 7);
	// Passing non-dimensional units to the solver.
	for (int j = 0; j < 3; ++j) {
		c[j] /= ASTRO_AU;
		c[j + 3] /= ASTRO_EARTH_VELOCITY;
	}
	c[6] /= m_spacecraft.get_mass();
	//Throttles constraints.
	m_leg.get_throttles_con(c.begin() + 7, c.begin() + 7 + m_n_seg);
}


std::string pl2pl_fixed_time::get_name() const
{
	return "pl2pl_fixed_time";
}

bool pl2pl_fixed_time::get_high_fidelity() const
{
	return m_leg.get_high_fidelity();
}

void pl2pl_fixed_time::set_high_fidelity(bool hf_option)
{
	reset_caches();
	m_leg.set_high_fidelity(hf_option);
}

/// A pretty description of the chromosome
/**
 * @return a string containing the human readable decision vector description
 */
std::string pl2pl_fixed_time::pretty(const std::vector<double> &x) const
{
	std::ostringstream s;
	//We start by filling up the m_legs with the correct information
	constraint_vector c(this->get_c_dimension());
	this->compute_constraints_impl(c,x);
	s << "Final Mass: " << x[0] << " Kg" << std::endl;
	s << this->m_leg << std::endl;
	s << std::endl;

	return s.str();
}

}} // namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::pl2pl_fixed_time)
