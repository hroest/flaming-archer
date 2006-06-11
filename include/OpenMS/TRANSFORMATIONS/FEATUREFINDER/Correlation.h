// -*- C++: make; tab-width: 2; -*-
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
// $Id: Correlation.h,v 1.10 2006/05/01 09:46:23 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_CORRELATION_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_CORRELATION_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseQuality.h>

namespace OpenMS
{
	/** @brief Measures the quality of a modelfit to some realworld data. 
	  
	 		Implementation of class BaseQuality. The quality is measured as
	 		the cross correlation of data and model.
	 		
	 		@ingroup FeatureFinder
	 	*/
  class Correlation 
    : public BaseQuality
  {

  public:
    /// standard constructor
    Correlation();

    /// destructor 
    virtual ~Correlation();

    /// evaluates the quality of the fit of @p model to @p set
    double evaluate(const IndexSet& set, const BaseModel<2>& model);
    
    /// evaluates the quality of the fit of @p model to @p set along dimension @p dim
    double evaluate(const IndexSet& set, const BaseModel<1>& model, UnsignedInt dim);

		/// creates instance of this class (function is called by factory).
    static BaseQuality* create()
    {
      return new Correlation();
    }

    static const String getName()
    {
      return "Correlation";
    }
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_CORRELATION_H
