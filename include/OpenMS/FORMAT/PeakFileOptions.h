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
// $Maintainer: Thomas Kadauke $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PEAKFILEOPTIONS_H
#define OPENMS_FORMAT_PEAKFILEOPTIONS_H

#include <OpenMS/DATASTRUCTURES/DRange.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief Options for loading files containing peak data.

		@ingroup FileIO
	*/
	class PeakFileOptions
	{
	public:
		typedef std::vector<int> MSLevels;
		
		///Default constructor
		PeakFileOptions();
		///Destructor
		~PeakFileOptions();

		///sets whether or not to load only meta data
		void setMetadataOnly(bool only);
		///returns whether or not to load only meta data
		bool getMetadataOnly() const;
		
		///sets whether or not to write supplemental peak data in MzData files
		void setWriteSupplementalData(bool write);
		///returns whether or not to write supplemental peak data in MzData files
		bool getWriteSupplementalData() const;
		
		///restricts the range of RT values for peaks to load
		void setRTRange(const DRange<1>& range);
		///returns @c true if an RT range has been set
		bool hasRTRange() const;
		///returns the RT range
		const DRange<1>& getRTRange() const;
		
		///restricts the range of MZ values for peaks to load
		void setMZRange(const DRange<1>& range);
		///returns @c true if an MZ range has been set
		bool hasMZRange() const;
		///returns the MZ range
		const DRange<1>& getMZRange() const;
		
		///restricts the range of intensity values for peaks to load
		void setIntensityRange(const DRange<1>& range);
		///returns @c true if an intensity range has been set
		bool hasIntensityRange() const;
		///returns the intensity range
		const DRange<1>& getIntensityRange() const;
		
		///sets the desired MS levels for peaks to load
		void setMSLevels(const MSLevels& levels);
		///adds a desired MS level for peaks to load
		void addMSLevel(int level);
		///clears the MS levels
		void clearMSLevels();
		///returns @c true, if MS levels have been set
		bool hasMSLevels() const;
		///returns @c true, if MS level @p level has been set
		bool containsMSLevel(int level) const;
		///returns the set MS levels
		const MSLevels& getMSLevels() const;

	private:
		bool metadata_only_;
		bool write_supplemental_data_;
		bool has_rt_range_;
		bool has_mz_range_;
		bool has_intensity_range_;
		DRange<1> rt_range_;
		DRange<1> mz_range_;
		DRange<1> intensity_range_;
		MSLevels ms_levels_;
	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_PEAKFILEOPTIONS_H
