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
#ifndef OPENMS_COMPARISON_SPECTRA_ZHANGSIMILARITYSCORE_H
#define OPENMS_COMPARISON_SPECTRA_ZHANGSIMILARITYSCORE_H

#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{

  /**
	  @brief Similarity score of Zhang

		The details of the score can be found in:
		Z. Zhang, Prediction of Low-Energy Collision-Induced Dissociation Spectra of Peptides,
		Anal. Chem., 76 (14), 3908 -3922, 2004

		@param epsilon - defines the absolut error of the mass spectrometer; default value is 0.2 Th
		
		@ingroup SpectraComparison
  */
	
  class ZhangSimilarityScore : public CompareFunctor
  {
  public:
	
		// @name Constructors and Destructors
		// @{
    /// default constructor
    ZhangSimilarityScore();

    /// copy constructor
    ZhangSimilarityScore(const ZhangSimilarityScore& source);

    /// destructor
    virtual ~ZhangSimilarityScore();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    ZhangSimilarityScore& operator = (const ZhangSimilarityScore& source);
	
		/// 
		double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const;
		// @}

		// @name Accessors
		// @{
		///
    static CompareFunctor* create() { return new ZhangSimilarityScore(); }

		///
		static const String getName()
		{
			return "ZhangSimilarityScore";
		}

		// @}

  };

}
#endif //OPENMS_COMPARISON_SPECTRA_ZHANGSIMILARTIYSCORE_H
