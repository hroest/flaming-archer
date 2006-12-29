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
#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDREPCOMPAREFUNCTOR_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDREPCOMPAREFUNCTOR_H

#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

namespace OpenMS
{

  /**
		@defgroup Comparison Comparison
	*/

	/**
		@defgroup SpectraComparison Spectra Comparison
	*/

	/**
	
		@brief Base class for compare functors of spectra; compare functors returns a similiarity value of two spectra
	
  	BinnedRepCompareFunctor classes return a value for a pair of BinnedRep objects
  	ideally the value should reflect the similarity of the pair
  	similarities of spectra should be > 0
		
  	@param filterwindow
		
    maximum mass difference for spectra that get similarity > 0

		@ingroup Comparison
  */
  class BinnedRepCompareFunctor : public FactoryProduct
  {

  public:

    /// default constructor
    BinnedRepCompareFunctor();

    /// copy constructor
    BinnedRepCompareFunctor(const BinnedRepCompareFunctor& source);

    /// destructor
    virtual ~BinnedRepCompareFunctor();

    /// assignment operator
    BinnedRepCompareFunctor& operator = (const BinnedRepCompareFunctor& source);

    /// function call operator, calculates the similarity
    virtual double operator () (const BinnedRep&, const BinnedRep&) const = 0;

    /// function call operator, calculates the self similarity
    virtual double operator () (const BinnedRep& a) const = 0;

		/// registers all derived products 
		static void registerChildren();

  protected:

  };

}
#endif // OPENMS_COMPARISON_SPECTRA_COMPAREFUNCTOR_H
