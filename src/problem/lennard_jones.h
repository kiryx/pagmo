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

#ifndef PAGMO_PROBLEM_LENNARD_JONES_H
#define PAGMO_PROBLEM_LENNARD_JONES_H

#include <string>
#include <vector>

#include "../types.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// The Lennard-Jones problem.
/**
 * \image html lennardjones.jpg "Minimum energy onfiguration with 38 atoms."
 * \image latex lennardjones.jpg "Minimum energy onfiguration with 38 atoms." width=5cm
 *
 * This is a box-constrained continuous single-objecive problem. Depending on the number of
 * atoms, the global optima will be different. In the link below a database containing all
 * putative global optima is given.
 *
 * @see http://physchem.ox.ac.uk/~doye/jon/structures/LJ/tables.150.html
 * @author Dario Izzo (dario.izzo@esa.int)
 */

class __PAGMO_VISIBLE lennard_jones : public base
{
	public:
		lennard_jones(int);
		base_ptr clone() const;
		std::string get_name() const;
	private:
		static double r(const int& atom, const int& coord, const std::vector <double>& x);
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
};

}} //namespaces

#endif // PAGMO_PROBLEM_LENNARD_JONES_H
