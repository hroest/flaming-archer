// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_PEAKMARKER_H
#define OPENMS_FILTERING_TRANSFORMERS_PEAKMARKER_H

#include <vector>
#include <map>
#include <string>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

namespace OpenMS
{
  /**
		@defgroup PeakMarker Peak Marker

		@ingroup Filtering
	*/

	/**
		@ingroup PeakMarker
	
  	@brief PeakMarker marks peaks that seem to fulfil some criterion
  */
  class PeakMarker
    : public FactoryProduct
  {
  public:
    /// standard constructor
    PeakMarker() ;

    /// copy constructor
    PeakMarker(const PeakMarker& source);

    /// destructor
    virtual ~PeakMarker();

    /// assignment operator
    PeakMarker& operator = (const PeakMarker& source);

		/// method to mark peaks, which are with their positions in @param marked afterwards
		template <typename SpectrumType> void apply(std::map<double, bool>& marked, SpectrumType& spectrum) = 0;
  };

}
#endif // OPENMS_FILTERING_TRANSFORMERS_PEAKMARKER_H
