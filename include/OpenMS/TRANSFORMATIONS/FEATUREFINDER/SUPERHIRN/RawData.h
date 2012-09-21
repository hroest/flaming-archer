// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
/*
 *  RawData.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#ifndef _RAWDATA_h_
#define _RAWDATA_h_
#include <vector>
#include <ostream>

#include <OpenMS/config.h>

namespace OpenMS
{

// Class for the storage of raw MS data
	class OPENMS_DLLAPI RawData
	{

		public:

			RawData() {};
			RawData(std::vector<double>&, std::vector<double>&);
			virtual ~RawData();

			friend std::ostream& operator<<(std::ostream&, RawData&);

			/*
			 * @brief Retrieve raw data as mass and intensity vectors. First argument: Mass values in profile mode
			 * Second argument: Intensity values in profile mode
			 */
			void get(std::vector<double>&, std::vector<double>&);

			/*
			 * @brief Set raw data as mass and intensity vectors. First argument: Mass values in profile mode
			 * Second argument: Intensity values in profile mode
			 */
			void set(std::vector<double>&, std::vector<double>&);

			// Virtual functions
			virtual void smooth()
			{
			}

		protected:
			std::vector<double> profileMasses_;
			std::vector<double> profileIntensities_;
	};

	std::ostream& operator<<(std::ostream& out, RawData& data);

	inline RawData::RawData(std::vector<double>& masses,	std::vector<double>& intensities)
	{
		profileMasses_ = masses;
		profileIntensities_ = intensities;
	}


	inline void RawData::get(std::vector<double> &masses, std::vector<double> &intensities)
	{
		masses = profileMasses_;
		intensities = profileIntensities_;
	}

	inline void RawData::set(std::vector<double> &masses, std::vector<double> &intensities)
	{
		profileMasses_ = masses;
		profileIntensities_ = intensities;
	}


} // ns

#endif
