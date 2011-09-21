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
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MASSTRACE_H
#define OPENMS_KERNEL_MASSTRACE_H

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/FeatureMap.h>


#include <vector>
#include <list>
#include <map>


namespace OpenMS
{
    typedef Peak2D PeakType;



    /**@brief MassTrace kernel class
                @ingroup Kernel
        */
    class OPENMS_DLLAPI MassTrace
        : public std::list<PeakType>
    {
    public:

        /// Default constructor
        MassTrace();
        /// Default destructor
        ~MassTrace();
        /// Copy constructor
        MassTrace(const MassTrace &tr_obj);

        MassTrace& operator= (const MassTrace& rhs);

        /// getter & setter

        inline String getLabel()
        {
            return label_;
        }

        inline void setLabel(const String& label)
        {
            label_ = label;
        }

        inline DoubleReal getCentroidMZ()
        {
            return centroid_mz_;
        }

        inline DoubleReal getCentroidRT()
        {
            return centroid_rt_;
        }

        inline std::vector<DoubleReal> getSmoothedIntensities()
        {
            return smoothed_intensities_;
        }

        inline std::vector<DoubleReal> getSmoothedIntensities() const
        {
            return smoothed_intensities_;
        }

        inline void setSmoothedIntensities(const std::vector<DoubleReal>& db_vec)
        {
            smoothed_intensities_ = db_vec;
        }

        inline DoubleReal getSmoothedMaxRT()
        {
            Size max_idx(this->findSmoothedMaxIdx());

            MassTrace::const_iterator c_it = this->begin();
            std::advance(c_it, max_idx);

            return c_it->getRT();
        }

        //        inline Size getMinimumPeakWidth()
        //        {
        //            // ????
        //        }
        /// prepend & append peaks, update centroid mz
        void prependPeak(PeakType);
        void appendPeak(PeakType);

        DoubleReal computeWeightedMeanMZ();
        DoubleReal computePeakArea();
        Size findSmoothedMaxIdx();
        DoubleReal findMaxPeakRT();
        DoubleReal estimateFWHM();
        void findLocalExtrema(const Size&, std::vector<Size>&, std::vector<Size>&);
        ConvexHull2D getConvexhull() const;

    private:
        void updateMedianMZ_();
        void updateMeanMZ_();
        void updateIterativeWeightedMeanMZ_(const PeakType&);

        DoubleReal centroid_mz_;
        DoubleReal centroid_rt_;

        String label_;

        std::vector<DoubleReal> smoothed_intensities_;

        DoubleReal prev_counter_;
        DoubleReal prev_denom_;


    };

}

#endif // OPENMS_KERNEL_MASSTRACE_H
