// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASESUPERIMPOSER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASESUPERIMPOSER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <utility>
#include <fstream>

namespace OpenMS
{

  /**
    @brief The base class of all superimposer algorithms.

   This class defines the basic interface for all superimposer algorithms. It works on several element maps and
   computes transformations, that map the elements of the maps as near as possible to each other.
   
   The element map must be a random access container (e.g. vector, DPeakArray, FeatureMap)
   of elements that have the same interface as RawDataPoint2D.
  */
  template <typename MapT>
  class BaseSuperimposer 
  	: public FactoryProduct
  {
  	
  	public:

	    /// Container for input elements
	    typedef MapT ElementMapType;
	
	    /// Constructor
	    BaseSuperimposer()
	    	: FactoryProduct("BaseSuperimposer")
	    {
	    }
	
	    /// Destructor
	    virtual ~BaseSuperimposer()
	  	{
	  	}

	    /**
	    	@brief Estimates the transformation and fills the given mapping function
	    	
	    	@exception IllegalArgument is thrown if the input maps are invalid.
	    */
	    virtual void run(const std::vector<ElementMapType>& maps, LinearMapping& mapping) = 0;
	
	    /// Register all derived classes here
	    static void registerChildren();
  
  }; // BaseSuperimposer

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASESUPERIMPOSER_H
