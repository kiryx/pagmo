/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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

#ifndef PAGMO_TOPOLOGY_BASE_H
#define PAGMO_TOPOLOGY_BASE_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../config.h"
#include "../serialization.h"

namespace pagmo {

/// Topology namespace.
/**
 * This namespace contains all the topologies implemented in PaGMO.
 */
namespace topology {

/// Base topology class.
class base;

/// Alias for shared pointer to base topology.
typedef boost::shared_ptr<base> base_ptr;

/// Base topology class.
/**
 * This class represents a topology connecting islands in an archipelago
 * using a directed graph in which the vertices are uniquely identified by an integer index representing:
 * - the order of creation of the vertex (i.e., the first vertex added has index 0, the second one has integer 1 and so on),
 * - the positional index of the island in the archipelago (i.e., the first island added to the archipelago corresponds
 *   to vertex 0, the second one to vertex 1 and so on).
 *
 * The internal implementation of the graph uses the Boost graph library, and most methods of this class are thin wrappers around the corresponding
 * Boost graph functions.
 *
 * The user is required to implement the connect() method, which will be called upon island insertion in an archipelago to establish the connection(s)
 * between the newly-added island and the islands already present in the archipelago.
 *
 * The re-implementable methods are:
 * - get_name(), to give a user-friendly name to the derived topology (otherwise the mangled C++ name will be used),
 * - human_readable_extra(), to provide extra, topology-specific information when printing the topology to a stream.
 *
 * @see http://www.boost.org/doc/libs/release/libs/graph/doc/
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 * @author Marek Ruciński (marek.rucinski@gmail.com)
 */
class __PAGMO_VISIBLE base
{
	protected:
		// Useful shortcut typedefs for graph-related types. The graph is directed and bidirectional,
		// since we need access to both in and out edges.
		/// Underlying graph type for the representation of the topology.
		/**
		 * The graph is a directed and bidirectional adjacency list.
		 */
		typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::bidirectionalS,boost::no_property,boost::property<boost::edge_weight_t, int, boost::property<boost::edge_name_t, double> >,boost::vecS> graph_type;
		/// Iterator over the vertices.
		typedef boost::graph_traits<graph_type>::vertex_iterator v_iterator;
		/// Iterator over adjacent vertices.
		typedef boost::graph_traits<graph_type>::adjacency_iterator a_iterator;
		/// Iterator over inversely adjacent vertices.
		typedef graph_type::inv_adjacency_iterator ia_iterator;
		/// Vertex descriptor.
		typedef boost::graph_traits<graph_type>::vertex_descriptor v_descriptor;
		/// Edge descriptor.
		typedef boost::graph_traits<graph_type>::edge_descriptor e_descriptor;
	public:
		/// Vertices size type.
		typedef graph_type::vertices_size_type vertices_size_type;
		/// Edges size type.
		typedef graph_type::edges_size_type edges_size_type;
		base();
		base(const base &);
		base &operator=(const base &);
		/// Clone method.
		/**
		 * Provided that the derived topology implements properly the copy constructor, virtually all implementations of this method will
		 * look like this:
		 * @code
		 * return base_ptr(new derived_topology(*this));
		 * @endcode
		 *
		 * @return topology::base_ptr to a copy of this.
		 */
		virtual base_ptr clone() const = 0;
		virtual ~base();
		virtual std::string get_name() const;
		std::string human_readable_terse() const;
		std::string human_readable() const;
		/** @name High-level graph access and manipulation methods. */
		//@{
		vertices_size_type get_number_of_vertices() const;
		edges_size_type get_number_of_edges() const;
		void push_back();
		double get_average_shortest_path_length() const;
                double get_clustering_coefficient() const;
                std::vector<double> get_degree_distribution();
		bool are_adjacent(const vertices_size_type &, const vertices_size_type &) const;
		bool are_inv_adjacent(const vertices_size_type &,const vertices_size_type &) const;
		std::vector<vertices_size_type> get_v_adjacent_vertices(const vertices_size_type &) const;
		std::vector<vertices_size_type> get_v_inv_adjacent_vertices(const vertices_size_type &) const;
		edges_size_type get_num_adjacent_vertices(const vertices_size_type &) const;
		edges_size_type get_num_inv_adjacent_vertices(const vertices_size_type &) const;
		void set_weight(const vertices_size_type &, const vertices_size_type &, double);
		double get_weight(const vertices_size_type &, const vertices_size_type &);
		//@}
	protected:
		/** @name Low-level graph access and manipulation methods. */
		//@{
		void add_vertex();
		std::pair<a_iterator,a_iterator> get_adjacent_vertices(const vertices_size_type &) const;
		std::pair<ia_iterator,ia_iterator> get_inv_adjacent_vertices(const vertices_size_type &) const;
		void add_edge(const vertices_size_type &, const vertices_size_type &);
		void remove_edge(const vertices_size_type &, const vertices_size_type &);
		void remove_all_edges();
		std::pair<v_iterator,v_iterator> get_vertices() const;
		/// Establish connections between islands during a push_back() operation.
		/**
		 * This method will be called by push_back() after a vertex has been added to the graph.
		 * The purpose of this method is to connect the newly-added vertex to other vertices according to the properties of the topology.
		 *
		 * @param[in] idx index of the newly-added vertex.
		 */
		virtual void connect(const vertices_size_type &idx) = 0;
		//@}
		virtual std::string human_readable_extra() const;
	private:
		void check_vertex_index(const vertices_size_type &) const;
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & m_graph;
		}
	private:
		graph_type m_graph;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base &);

}}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::topology::base)

#endif
