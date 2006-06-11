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
// $Id: SimpleFeaFiTraits.h,v 1.31 2006/04/25 12:41:44 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEFEAFITRAITS_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEFEAFITRAITS_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseFeaFiTraits.h>
#include <OpenMS/DATASTRUCTURES/ScanIndex.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

	/**
	  @brief Simple (as one might guess) implementation of the BaseFeaFiTraits class.
	  
		<table>
			<tr><td>min_intensity</td>
					<td>minimum intensity of a data point i.e. data with an intensity
					below this value are removed from the data set before the algorithm
					starts.</td></tr>
		</table>	 
			  	  
	 */
  class SimpleFeaFiTraits 
    : public BaseFeaFiTraits 
  {

  public:
    /// standard constructor
    SimpleFeaFiTraits();

    /// destructor 
    virtual ~SimpleFeaFiTraits();

    /// add peak to internal datastructure 
    void addSinglePeak(const DRawDataPoint<2>& peak);
    
    /// set fill the internal data structure using an instance of MSExperiment (not implemented in this traits class).
    void setData(MSExperiment<DPeak<1> >& exp);

	   /// Non-mutable acess flag with Index @p index 
    const Flag& getPeakFlag(const UnsignedInt index) const throw (Exception::IndexOverflow);
    /// Mutable acess flag with Index @p index
    Flag& getPeakFlag(const UnsignedInt index) throw (Exception::IndexOverflow);

    /// acess range of flags through pointers 
    const FlagRefVector& getFlags(const IndexSet& range) throw (Exception::IndexOverflow);

    /// Non-mutable acess to all flags 
    const FlagVector& getAllFlags() const;
    /// Mutable acess to all flags 
    FlagVector& getAllFlags();

    /// acess peak with Index @p index 
    const PeakType& getPeak(const UnsignedInt index) const throw (Exception::IndexOverflow);
    /// acess range of peaks 
    const PeakRefVector& getPeaks(const IndexSet& range) throw (Exception::IndexOverflow);
    /// acess all peaks 
    const PeakVector& getAllPeaks();
    /// retrieve the number of peaks 
    const UnsignedInt getNumberOfPeaks();

		 /// acess intensity of peak with Index @p index
    const IntensityType& getPeakIntensity(const UnsignedInt index) const throw (Exception::IndexOverflow);
    /// acess m/z of peak with Index @p index 
    const CoordinateType& getPeakMz(const UnsignedInt index) const throw (Exception::IndexOverflow);
    /// acess retention time of peak with Index @p index
    const CoordinateType& getPeakRt(const UnsignedInt index) const throw (Exception::IndexOverflow);
    /// acess scan number of peak with Index @p index 
    const UnsignedInt getPeakScanNr(const UnsignedInt index) const throw (Exception::IndexOverflow);
   
    /** @brief get index of next peak in m/z dimension 
    
       \param index of the peak whose successor is requested
       \return index of the next peak 
    */
    UnsignedInt getNextMz(const UnsignedInt index) const throw (Exception::IndexOverflow, NoSuccessor);

    /** @brief get index of previous peak in m/z dimension 
    
       \param index of the peak whose predecessor is requested
       \return index of the previous peak
    */
    UnsignedInt getPrevMz(const UnsignedInt index) const throw (Exception::IndexOverflow, NoSuccessor);

    /** @brief get index of next peak in retention time dimension 
    
       \param index of the peak whose successor is requested
       \return index of the next peak
    */
    UnsignedInt getNextRt(const UnsignedInt index) const throw (Exception::IndexOverflow, NoSuccessor);

    /** @brief get index of next peak in retiontion time dimension 
    
       \param index of the peak whose predecessor is requested
       \return index of the previous peak
    */
    UnsignedInt getPrevRt(const UnsignedInt index) const throw (Exception::IndexOverflow, NoSuccessor);

    /// run main loop using set seeders, extenders and fitters
    const FeatureVector& run();
 	    
 	  /// returns an instance of this traits class 
    static BaseFeaFiTraits* create()
    {
      return new SimpleFeaFiTraits();
    }

    static const String getName()
    {
      return "SimpleFeaFiTraits";
    }
   
  protected:
  
  	/** @brief Sorts the peaks according by position.
  	 
	       In 1D m/z, in the 2D case m/z and rt. That is,
	 			 the peaks are first sorted by their rt value
	 			 and peaks with equal rt (i.e. scan index) are 
	 			 then sorted by m/z. In addition,
	       we initialise the vector of scan indizes
	 			 to quickly retrieve the scan number of a peak.
	 	*/
  	void sortData_();

    FlagVector flags_;
    FlagRefVector selected_flags_;

    PeakVector peaks_;
    PeakRefVector selected_peaks_;
    
    /// stores reference to the scan numbers for each peak
    ScanIndex<PeakVector> scan_index_;
    
     /// Threshold for the smallest accepted intensity.
    IntensityType min_intensity_;
         
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEFEAFITRAITS_H
