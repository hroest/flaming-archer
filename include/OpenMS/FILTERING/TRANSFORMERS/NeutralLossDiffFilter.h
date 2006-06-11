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
// $Id: NeutralLossDiffFilter.h,v 1.4 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSDIFFFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSDIFFFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <string>

namespace OpenMS
{
  /**
  NeutralLossDiffFilter returns the total intensity ob peak pairs whose m/z difference can be explained by a neutral loss <br>
  
  \param tolerance m/z tolerance
  */
  class NeutralLossDiffFilter : public FilterFunctor
  {
  public:
    /** @brief standard constructor <br> */
    NeutralLossDiffFilter();

    /** @brief copy constructor <br> */
    NeutralLossDiffFilter(const NeutralLossDiffFilter& source );

    /** @brief assignment operator <br> */
    NeutralLossDiffFilter& operator=(const NeutralLossDiffFilter& source);

    /** @brief destructor <br> */
    ~NeutralLossDiffFilter();

    static FactoryProduct* create() { return new NeutralLossDiffFilter();}

    //std::vector<double> operator()(const ClusterSpectrum& spec);

		template <typename SpectrumType> double apply(SpectrumType& spectrum)
		{
    	double tolerance = (double)param_.getValue("tolerance");
    	double isodiff = 0;
    	//iterate over all peaks
    	for (int i = 0; i < (int)spectrum.size(); ++i)
    	{
      	for (int j = 1; i-j >= 0; ++j)
      	{
        	if (fabs(spectrum.getContainer()[i-j].getPosition()[0] - spectrum.getContainer()[i].getPosition()[0] - 18) < tolerance || // water
						  fabs(spectrum.getContainer()[i-j].getPosition()[0] - spectrum.getContainer()[i].getPosition()[0] - 17) < tolerance) //ammoniom
        	{
          	isodiff += spectrum.getContainer()[i-j].getIntensity() + spectrum.getContainer()[i].getIntensity();
        	}
        	else if (fabs(spectrum.getContainer()[i-j].getPosition()[0] - spectrum.getContainer()[i].getPosition()[0]  ) > 18 + tolerance)
        	{
          	break;
        	}
      	}
    	}
    	//vector<double> result;
    	//result.push_back(isodiff);
    	return isodiff;
		}

    //String info() const;

		static const String getName()
		{
			return "NeutralLossDiffFilter";
		}

  private:
    //static const String info_;
  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSDIFFFILTER_H
