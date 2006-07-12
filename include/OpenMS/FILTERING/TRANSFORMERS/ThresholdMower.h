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
#ifndef OPENMS_FILTERING_TRANSFORMERS_THRESHOLDMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_THRESHOLDMOWER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>

namespace OpenMS
{
  /**
  	@brief ThresholdMover removes all peaks below a Threshold
  
	 	@param threshold: the threshold
	  for comparable results we suggest normalizing (for example with Normalizer) all
	  Spectra first

		@ingroup SpectraPreprocessing
  */
  class ThresholdMower
    :	public PreprocessingFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    ThresholdMower();

    /// copy constructor
    ThresholdMower(const ThresholdMower& source);

    /// destructor
    virtual ~ThresholdMower();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    ThresholdMower& operator=(const ThresholdMower& source);
		// @}

		// @name Accessors
		// @{
		///
    static FactoryProduct* create() { return new ThresholdMower(); }

		///
		template <typename SpectrumType> void apply(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::Iterator Iterator;

			double threshold = (double)param_.getValue("threshold");

			for (Iterator it = spectrum.begin(); it != spectrum.end(); )
			{
				if (it->getIntensity() < threshold)
				{
					it = spectrum.getContainer().erase(it);
				}
				else
				{
					++it;
				}
			}
		}

		/// 
		static const String getName()
		{
			return "ThresholdMower";
		}
		// @}
  };

}

#endif //OPENMS_FILTERING_TRANSFORMERS_THRESHOLDMOWER_H
