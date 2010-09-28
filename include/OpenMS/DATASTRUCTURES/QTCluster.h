// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Steffen Sass $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_QTCLUSTER_H
#define OPENMS_DATASTRUCTURES_QTCLUSTER_H

#include<OpenMS/DATASTRUCTURES/GridFeature.h>

namespace OpenMS {

/**
	 @brief A representation of a QT cluster used for feature grouping.

	 Ultimately, a cluster represents a group of corresponding features (or consensus features) from different input maps (feature maps or consensus maps).

	 Clusters are defined by their center points (one feature each). A cluster also stores a number of potential cluster elements (other features) from different input maps, together with their distances to the cluster center.
	 Every feature that satisfies certain constraints with respect to the cluster center is a @e potential cluster element. However, since a feature group can only contain one feature from each input map, only the "best" (i.e. closest to the cluster center) such feature is considered a true cluster element.

	 The QT clustering algorithm has the characteristic of initially producing all possible, overlapping clusters. Iteratively, the best cluster is then extracted and the clustering is recomputed for the remaining points.

	 In our implementation, multiple rounds of clustering are not necessary. Instead, the clustering is updated in each iteration. This is the reason for storing all potential cluster elements: When a certain cluster is finalized, its elements have to be removed from the remaining clusters, and affected clusters change their composition. (Note that clusters can also be invalidated by this, if the cluster center is being removed.)
	 
	 The quality of a cluster is the normalized average distance to the cluster center for present and missing cluster elements. The distance value for missing elements (if the cluster contains no feature from a certain input map) is the user-defined threshold that marks the maximum allowed radius of a cluster.

	 @see QTClusterFinder

	 @ingroup Datastructures
*/

	class OPENMS_DLLAPI QTCluster
	{
	private:
		/**
		 * @brief Mapping: input map -> distance to center -> neighboring point
		 * @note There should never be an empty sub-map! (When a sub-map becomes empty, it should be removed from the overall map.)
		 */
		typedef std::map<Size, std::multimap<DoubleReal, GridFeature*> >
			NeighborMap;
		
		/// Pointer to the cluster center
		GridFeature* center_point_;

		/**
		 * @brief Neighbors of the cluster center, sorted by distance, for different input maps.
		 *
		 * The first (best) point in each sub-map is considered a cluster element.
		 */
		NeighborMap neighbors_;

		/// Maximum distance of a point that can still belong to the cluster
		DoubleReal max_distance_;

		/// Number of input maps
		Size num_maps_;

		/// Quality of the cluster
		Size quality_;

		/// Has the cluster changed (if yes, quality needs to be recomputed)?
		bool changed_;

		/// Base constructor (not accessible)
		QTCluster();

		/// Computes the quality of the cluster
		void computeQuality_();

	public:
		/**
		 * @brief Detailed constructor
		 * @param center_point Pointer to the center point
		 * @param num_maps Number of input maps
		 * @param max_distance Maximum allowed distance of two points
		 */
		QTCluster(GridFeature* center_point, Size num_maps, 
							DoubleReal max_distance);

		/// Destructor
		virtual ~QTCluster();

		/// Assignment operator
		QTCluster& operator=(const QTCluster& rhs);

		/// Returns the RT value of the cluster
		DoubleReal getCenterRT();

		/// Returns the m/z value of the cluster center
		DoubleReal getCenterMZ();

		/// Returns the size of the cluster (number of elements, incl. center)
		Size size() const;

		/// Compare by quality
		bool operator<(QTCluster &cluster);

		/**
		 * @brief Adds a new element/neighbor to the cluster
		 * @note There is no check whether the element/neighbor already exists in the cluster!
		 * @param element The element to be added
		 * @param distance Distance of the element to the center point
		 */
		void add(GridFeature* element, DoubleReal distance);

		/// Non-mutable access to the clustered elements
		void getElements(std::map<Size, GridFeature*>& elements) const;

		/**
		 * @brief Updates the cluster after data points were removed
		 * @return Whether the cluster is still valid (it's not if the cluster center is among the removed points).
		 */
		bool update(const std::map<Size, GridFeature*>& removed);

		/// Returns the cluster quality
		DoubleReal getQuality();
	};
}

#endif // OPENMS_DATASTRUCTURES_QTCLUSTER_H
