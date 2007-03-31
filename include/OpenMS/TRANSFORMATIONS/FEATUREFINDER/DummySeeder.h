// -*- C++: make; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYSEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYSEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/DSignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>

namespace OpenMS
{
	/** 
		@brief Seeding module which selects single peaks based on their s/n ratio.
		
		Groups of peaks a clustered within a certain distance and traced over consecutive scans.
		
		<table>
			<tr>
				<td>min_snratio</td>
				<td>Minimum signal to noise ratio of a peak</td>
			</tr>
			<tr>
				<td>tolerance_mz</td>
				<td>tolerance in m/z for a peak in the previous scan.</td>
			</tr>
			<tr>
				<td>min_number_scans</td>
				<td>minimum number of scan per isotopic cluster.</td>
			</tr>
				<td>min_number_peaks</td>
				<td>minimum number of peaks per cluster</td>
			</tr>
			<td>max_peak_distance</td>
				<td>minimum number of peaks per cluster</td>
			</tr>
		</table>
		
		@ingroup FeatureFinder
	*/ 
  class DummySeeder 
    : public BaseSeeder
  {
  	public:	
			typedef FeaFiTraits::IntensityType IntensityType;
		  typedef FeaFiTraits::CoordinateType CoordinateType;
		  typedef DoubleReal ProbabilityType;	
	
			typedef FeaFiTraits::MapType MapType;
			typedef MapType::SpectrumType SpectrumType;
			typedef MapType::PeakType PeakType;
			typedef std::multimap<CoordinateType,IsotopeCluster> TableType;
	
		  /// Default constructor
	    DummySeeder();
	
	    /// destructor 
	    virtual ~DummySeeder();

	    /// Copy constructor
	    DummySeeder(const DummySeeder& rhs);
	    
	    /// Assignment operator
	    DummySeeder& operator= (const DummySeeder& rhs);
	
	    /// return next seed 
	    ChargedIndexSet nextSeed() throw (NoSuccessor);
	
	    static BaseSeeder* create()
	    {
	      return new DummySeeder();
	    }
	
	    static const String getProductName()
	    {
	      return "DummySeeder";
	    }
	
	  protected:
	 
	    /// Finds the neighbour of the peak denoted by @p current_mz in the previous scan
	    std::vector<double>::const_iterator searchInScan_(const std::vector<CoordinateType>::const_iterator& scan_begin, const std::vector<CoordinateType>::const_iterator& scan_end, CoordinateType current_mz)
	    {
	      // perform binary search to find the neighbour in rt dimension
	      // 	lower_bound finds the peak with m/z current_mz or the next larger peak if this peak does not exist.
	      std::vector<CoordinateType>::const_iterator insert_iter = lower_bound(scan_begin,scan_end,current_mz);
	
	      // the peak found by lower_bound does not have to be the closest one, therefore we have
	      // to check both neighbours
	      if ( insert_iter == scan_end ) // we are at the and have only one choice
	      {
	      	--insert_iter;
	      }
        // if the found peak is at the beginning of the spectrum,
        // there is not much we can do.
        else if ( insert_iter != scan_begin )
        {
          if ( *insert_iter - current_mz < current_mz - *(--insert_iter) )
          {
          	++insert_iter;    // peak to the right is closer
          }
	      }

				return insert_iter;
	    }

			/// Finds local maxima in the cwt
			void filterAndComputeLocalMax_(const SpectrumType & vec, 
														 						 							std::vector<int>& localmax
																											#ifdef DEBUG_FEATUREFINDER
														 													,UInt currscan_index
																											#endif
														 													);
																											
			/// Retrieves the iterator for a peak cluster at mz @p  curr_mz.
			TableType::iterator retrieveHashIter_(CoordinateType curr_mz, 
																													 CoordinateType& mz_in_hash, 
																													 const std::vector<CoordinateType>& iso_last_scan,
																													 UInt currscan_index );
	
		  /// Sweeps through scans and detects isotopic patterns
		  void sweep_();
		
		  /// stores the retention time of each isotopic cluster
		  TableType iso_map_;
		
		  /// Pointer to the current region
		  TableType::const_iterator curr_region_;
		
		  /// indicates whether the extender has been initialized
		  bool is_initialized_;
		 
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYSEEDER_H
