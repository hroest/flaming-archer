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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

// author: Marcel Grunert
// date: 7.7.2006
// Implementation in context of the bachelor thesis.

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_LMAGAUSSMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_LMAGAUSSMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>


namespace OpenMS
{

	class LmaGaussModel
		: public InterpolationModel<>
	{

	 public:
		typedef InterpolationModel<>::CoordinateType CoordinateType;
		typedef Math::BasicStatistics<CoordinateType > BasicStatistics;
		typedef InterpolationModel<> InterpolationModel;

		/// standard constructor
		LmaGaussModel();

		/// copy constructor
		LmaGaussModel(const LmaGaussModel& source);

		/// destructor
		virtual ~LmaGaussModel();

		/// assignment operator
		virtual LmaGaussModel& operator = (const LmaGaussModel& source);

		void setParam(const Param& param);

		void setParam(const BasicStatistics& statistics, CoordinateType scale_factor, CoordinateType standard_deviation, CoordinateType expected_value, CoordinateType min, CoordinateType max);

		/// get parameters (const access)
		const Param& getParam() const;

		/// get parameters
		Param& getParam();

		/// create new EmgModel object (needed by Factory)
		static BaseModel<1>* create()
		{
			return new LmaGaussModel();
  	}

		/// name of the model (needed by Factory)
		static const String getName()
		{
			return "LmaGaussModel";
		}

		/// set offset without being computing all over and without any discrepancy
		void setOffset(CoordinateType  offset);

		/// set sample/supporting points of interpolation
		void setSamples();
		
		/// get the center of the Gaussian model i.e. the position of the maximum
		const CoordinateType getCenter() const;

	 protected:
		CoordinateType  min_;
		CoordinateType  max_;
		BasicStatistics statistics_;
		CoordinateType scale_factor_;
		CoordinateType standard_deviation_;
		CoordinateType expected_value_;
	
	};
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_LMAGAUSSMODEL_H
