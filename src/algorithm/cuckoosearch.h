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

#ifndef PAGMO_ALGORITHM_CUCKOOSEARCH_H
#define PAGMO_ALGORITHM_CUCKOOSEARCH_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"


namespace pagmo { namespace algorithm {

/// Cuckoo Search Algorithm
/**
 *
 * Cuckoo Search (CS) is a population based algorithm that has been proposed in the year 2009 by Yang and Deb.
 * It was inspired by the obligate brood parasitism of some cuckoo species by laying their eggs in the
 * nests of other host birds (of other species). Some host birds can engage direct conflict with the intruding cuckoos.
 * @see http://arxiv.org/pdf/1003.1594.pdf for the paper on this algorithm
 *
 * @author Sunil K. Mahendrakar (sunil.mahendrakar19@gmail.com)
 */

class __PAGMO_VISIBLE cuckoosearch: public base
{
public:
	cuckoosearch(int gen=1, double pa=0.25);
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
		ar & const_cast<double &>(m_pa);
	}

	void get_cuckoos(const problem::base &, std::vector<decision_vector>&, const std::vector<decision_vector>&,
		const decision_vector) const;

	void empty_nests(const problem::base &, std::vector<decision_vector> &, const std::vector<decision_vector> &) const;

	void get_best_nest(const problem::base &, fitness_vector &, decision_vector &, std::vector<decision_vector> &,
		const std::vector<decision_vector> &, std::vector<fitness_vector> &) const;

	void applysimplebounds(const problem::base &, decision_vector &) const;

	// Number of generations
	const int m_gen;
	// Discovery rate of alien eggs/solutions
	const double m_pa;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::cuckoosearch)

#endif // CUCKOOSEARCH_H
