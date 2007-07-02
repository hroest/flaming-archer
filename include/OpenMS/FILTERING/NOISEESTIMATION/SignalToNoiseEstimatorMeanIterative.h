// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_FILTERING_NOISEESTIMATION_SIGNALTONOISEESTIMATORMEANITERATIVE_H
#define OPENMS_FILTERING_NOISEESTIMATION_SIGNALTONOISEESTIMATORMEANITERATIVE_H

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimator.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>

namespace OpenMS
{
  /**
    @brief Estimates the signal/noise (S/N) ratio of each data point in a scan
           based on an iterative scheme which discards high intensities
   
    For each datapoint in the given scan, we collect a range of data points around it (param: <i>WinLen</i>).
    The noise for a datapoint is estimated iteratively by discarding peaks which are more than
    (<i>StdevMP</i> * StDev) above the mean value. After three iterations, the mean value is 
    considered to be the noise level. If the number of elements in the current window is not sufficient (param: <i>MinRequiredElements</i>),
    the noise level is set to a default value (param: <i>NoiseForEmptyWindow</i>).
    
    The whole computation is histogram based, so the user will need to supply a number of bins (param: <i>BinCount</i>), which determines
    the level of error and runtime. The maximal intensity for a datapoint to be included in the histogram can be either determined 
    automatically (param: <i>AutoMode</i>) by two different methods or can be set directly by the user (param: <i>MaxIntensity</i>).
    
    Changing any of the parameters will invalidate the S/N values (which will invoke a recomputation on the next request).


    @note 
    Warning to *stderr* if sparse_window_percent > 20
            - percent of windows that have less than <i>MinRequiredElements</i> of elements
              (noise estimates in those windows are simply a constant <i>NoiseForEmptyWindow</i>).   
		 
    @ref SignalToNoiseEstimatorMeanIterative_Parameters are explained on a separate page.
    
    @ingroup Filtering
  */
  template < typename Container = MSSpectrum< > >
  class SignalToNoiseEstimatorMeanIterative : public SignalToNoiseEstimator< Container >
  {

    public:

      /// method to use for estimating the maximal intensity that is used for histogram calculation
      enum IntensityThresholdCalculation { MANUAL=-1, AUTOMAXBYSTDEV=0, AUTOMAXBYPERCENT=1 };

      using SignalToNoiseEstimator< Container >::stn_estimates_;
      using SignalToNoiseEstimator< Container >::first_;
      using SignalToNoiseEstimator< Container >::last_;
      using SignalToNoiseEstimator< Container >::is_result_valid_;
      using SignalToNoiseEstimator< Container >::defaults_;
      using SignalToNoiseEstimator< Container >::param_;
      
      typedef typename SignalToNoiseEstimator< Container >::PeakIterator PeakIterator;
      typedef typename SignalToNoiseEstimator< Container >::PeakType PeakType;
      
      typedef typename SignalToNoiseEstimator< Container >::GaussianEstimate GaussianEstimate;
      

      /// default constructor
      inline SignalToNoiseEstimatorMeanIterative()
      {
	    	//set the name for DefaultParamHandler error messages
	    	this->setName("SignalToNoiseEstimatorMeanIterative");	
    	
        defaults_.setValue("MaxIntensity", -1, "maximal intensity considered for histogram construction. By default, it will be calculated automatically (see AutoMode)."\
" Only provide this parameter if you know what you are doing (and change 'AutoMode' to '-1')!"\
" All intensities EQUAL/ABOVE 'MaxIntensity' will not be added to the histogram."\
" If you choose 'MaxIntensity' too small, the noise estimate might be too small as well."\
" If chosen too big, the bins become quite large (which you could counter by increasing 'BinCount', which increases runtime)."); 
        defaults_.setValue("AutoMaxStdevFactor", 3.0, "parameter for 'MaxIntensity' estimation (if 'AutoMode' == 0): mean + 'AutoMaxStdevFactor' * stdev"); 
        defaults_.setValue("AutoMaxPercentile", 95, "parameter for 'MaxIntensity' estimation (if 'AutoMode' == 1): AutoMaxPercentile th percentile"); 
        defaults_.setValue("AutoMode", 0, "method to use to determine maximal intensity: -1 --> use 'MaxIntensity'; 0 --> 'AutoMaxStdevFactor' method (default); 1 --> 'AutoMaxPercentile' method"); 
        defaults_.setValue("WinLen", 200.0, "window length in Thomson"); 
        defaults_.setValue("BinCount", 30, "number of bins used for histogram"); 
        defaults_.setValue("StdevMP", 3.0, "multiplier for stdev"); 
        defaults_.setValue("MinRequiredElements", 10, "minimum number of elements required in a window (otherwise it is considered sparse)"); 
        defaults_.setValue("NoiseForEmptyWindow", std::pow(10.0,20), "noise value used for sparse windows"); 

        SignalToNoiseEstimator< Container >::defaultsToParam_();
      }


      /// Copy Constructor
      inline SignalToNoiseEstimatorMeanIterative(const SignalToNoiseEstimatorMeanIterative&  source)
          : SignalToNoiseEstimator< Container >(source)
      {
        updateMembers_();
      }


      /** @name Assignment
       */
      //@{
      ///
      inline SignalToNoiseEstimatorMeanIterative& operator=(const SignalToNoiseEstimatorMeanIterative& source)
      {
        if(&source == this) return *this; 
        SignalToNoiseEstimator< Container >::operator=(source);
        updateMembers_();
        return *this;
      }
      //@}


      /// Destructor
      virtual ~SignalToNoiseEstimatorMeanIterative()
      {}
      
      /** @name Accessors
       */

      //@{
      ///

      /// Non-mutable access to the maximal intensity that is included in the histogram (higher values get discarded)
      inline DoubleReal getMaxIntensity() const     {    return max_intensity_;   }
      /// Mutable access to the maximal intensity that is included in the histogram (higher values get discarded)
      inline void setMaxIntensity(DoubleReal max_intensity)
      {
        max_intensity_ = max_intensity;
        param_.setValue("MaxIntensity", max_intensity_);
      }
//

      /// Non-Mutable access to the AutoMaxStdevFactor-Param, which holds a factor for stddev (only used if autoMode=1)
      inline DoubleReal getAutoMaxStdevFactor() const  {  return auto_max_stdev_Factor_;    }
      /// Mutable access to the AutoMaxStdevFactor-Param, which holds a factor for stddev (only used if autoMode=1)
      inline void setAutoMaxStdevFactor(DoubleReal value)
      {
        auto_max_stdev_Factor_ = value;
        param_.setValue("AutoMaxStdevFactor", auto_max_stdev_Factor_);
      }
//      

      /// get the AutoMaxPercentile-Param, which holds a percentile (only used if autoMode=2)
      inline DoubleReal getAutoMaxPercentile() const  {  return auto_max_percentile_;      }
      /// Mutable access to the AutoMaxPercentile-Param, which holds a percentile (only used if autoMode=2)
      inline void setAutoMaxPercentile(DoubleReal value)
      {
        auto_max_percentile_ = value;
        param_.setValue("AutoMaxPercentile", auto_max_percentile_);
      }      
//

      /// @brief -1 will disable it. 0 is default. 1 is alternative method
      /// Non-mutable access to AutoMode, which determines the heuristic to find MaxIntensity. See Class description.
      inline Int getAutoMode() const      {    return auto_mode_;     }
      /// @brief -1 will disable it. 0 is default. 1 is alternative method
      /// Mutable access to AutoMode, which determines the heuristic to find MaxIntensity. See Class description.
      inline void setAutoMode(Int auto_mode)
      {
        auto_mode_ = auto_mode;
        param_.setValue("AutoMode", auto_mode_);
      }
//

      /// Non-mutable access to the window length (in Thomson)
      inline DoubleReal getWinLen() const   {   return win_len_;      }
      /// Mutable access to the window length (in Thomson)
      inline void setWinLen(DoubleReal win_len)
      {
        win_len_ = win_len;
        param_.setValue("WinLen", win_len_);
      }

//
      /// Non-mutable access to the number of bins used for the histogram (the more bins, the better the approximation, but longer runtime)
      inline Int getBinCount() const       {   return bin_count_;      }
      /// Mutable access to the number of bins used for the histogram
      inline void setBinCount(Int bin_count)
      {
        bin_count_ = bin_count;
        param_.setValue("BinCount", bin_count_);
      }
//

      /// Non-mutable access to the multiplier of the stdev used during iterative threshold calculation
      inline DoubleReal getSTDEVMultiplier() const    {    return stdev_;     }
      /// Mutable access to the multiplier of the stdev used during iterative threshold calculation
      inline void setSTDEVMultiplier(DoubleReal stdev)
      {
        stdev_ = stdev;
        param_.setValue("StdevMP", stdev_);
      }

//
      /// Non-mutable access to the minimum required elements in a window, to be evaluated.
      inline Int getMinReqElements() const          {    return min_required_elements_;     }
      /// Mutable access to the minimum required elements in a window, to be evaluated.
      inline void setMinReqElements(Int min_required_elements)
      {
        min_required_elements_ = min_required_elements;
        param_.setValue("MinRequiredElements", min_required_elements_);
      }

//
      /// Non-mutable access to the noise value that is used if a window contains not enough elements
      inline DoubleReal getNoiseForEmtpyWindow() const     {   return noise_for_empty_window_;   }
      /// Mutable access to the noise value that is used if a window contains not enough elements
      inline void setNoiseForEmtpyWindow(DoubleReal noise_for_empty_window)
      {
        noise_for_empty_window_ = noise_for_empty_window;
        param_.setValue("NoiseForEmptyWindow", noise_for_empty_window_);
      }
      //@}



    protected:


      /// calculate StN values for all datapoints given, by using a sliding window approach
      /// @param scan_first_ first element in the scan
      /// @param scan_last_ last element in the scan (disregarded)
      virtual void computeSTN_(const PeakIterator& scan_first_, const PeakIterator& scan_last_) 
      throw(Exception::InvalidValue)
      {
        // reset counter for sparse windows
        double sparse_window_percent = 0;

        // reset the results
        stn_estimates_.clear();

        // maximal range of histogram needs to be calculated first
        if (auto_mode_ == AUTOMAXBYSTDEV)
        {
          // use MEAN+auto_max_intensity_*STDEV as threshold
          GaussianEstimate gauss_global = SignalToNoiseEstimator< Container >::estimate_(scan_first_, scan_last_);
          max_intensity_ = gauss_global.mean + std::sqrt(gauss_global.variance)*auto_max_stdev_Factor_;
        }
        else if (auto_mode_ == AUTOMAXBYPERCENT)
        {
          // get value at "auto_max_percentile_"th percentile
          // we use a histogram approach here as well.
          if ((auto_max_percentile_ < 0) || (auto_max_percentile_ > 100))
          {
            String s = auto_max_percentile_;
            throw Exception::InvalidValue(__FILE__, 
                                           __LINE__, 
                                           __PRETTY_FUNCTION__, 
                                           "AutoMode is on AUTOMAXBYPERCENT! AutoMaxPercentile is not in [0,100]. Use setAutoMaxPercentile(<value>) to change it!", 
                                           s);
          }

          std::vector <int> histogram_auto(100, 0);

          // find maximum of current scan
          int size = 0;
          typename PeakType::IntensityType maxInt = 0;
          PeakIterator run = scan_first_;
          while (run != scan_last_)
          {
            maxInt = std::max(maxInt, (*run).getIntensity());
            ++size;
            ++run;
          }

          double bin_size = maxInt / 100;

          // fill histogram
          run = scan_first_;
          while (run != scan_last_)
          {
            ++histogram_auto[(int) (((*run).getIntensity()-1) / bin_size)];
            ++run;
          }

          // add up element counts in histogram until ?th percentile is reached
          int elements_below_percentile = (int) (auto_max_percentile_ * size / 100);
          int elements_seen = 0;
          int i = -1;
          run = scan_first_;

          while (run != scan_last_ && elements_seen < elements_below_percentile)
          {
            ++i;
            elements_seen += histogram_auto[i];
            ++run;
          }

          max_intensity_ = (((double)i) + 0.5) * bin_size;
        }
        else //if (auto_mode_ == MANUAL)
        {
          if (max_intensity_<=0) 
          {
            String s = max_intensity_;
            throw Exception::InvalidValue(__FILE__, 
                                           __LINE__, 
                                           __PRETTY_FUNCTION__, 
                                           "AutoMode is on MANUAL! MaxIntensity is <=0. Needs to be positive! Use setMaxIntensity(<value>) or enable AutoMode!", 
                                           s);
          }
        }

        PeakIterator window_pos_center  = scan_first_;
        PeakIterator window_pos_borderleft = scan_first_;
        PeakIterator window_pos_borderright = scan_first_;

        double window_half_size = win_len_ / 2;
        double bin_size = max_intensity_ / bin_count_;

        std::vector <int> histogram(bin_count_, 0);
        std::vector <double> bin_value(bin_count_, 0);
        // calculate average intensity that is represented by a bin
        for (int bin=0; bin<bin_count_; bin++)
        {
          histogram[bin] = 0;
          bin_value[bin] = (bin + 0.5) * bin_size;
        }
        // index of last valid bin during iteration
        int hist_rightmost_bin;
        // bin in which a datapoint would fall
        int to_bin;
        // mean & stdev of the histogram
        double hist_mean;
        double hist_stdev;

        // tracks elements in current window, which may vary because of uneven spaced data
        int elements_in_window = 0;
        int window_count = 0;

        double noise;    // noise value of a datapoint

        // determine how many elements we need to estimate (for progress estimation)
        int windows_overall = 0;
        PeakIterator run = scan_first_;
        while (run != scan_last_)
        {
          ++windows_overall;
          ++run;
        }
        SignalToNoiseEstimator< Container >::startProgress(0,windows_overall,"noise estimation of data");

        // MAIN LOOP
        while (window_pos_center != scan_last_)
        {
          // erase all elements from histogram that will leave the window on the LEFT side
          while ( (*window_pos_borderleft).getMZ() <  (*window_pos_center).getMZ() - window_half_size )
          {
            //std::cout << "S: " << (*window_pos_borderleft).getMZ()  <<  " " << ( (*window_pos_center).getMZ() - window_half_size ) << "\n";
            to_bin = (int) (((*window_pos_borderleft).getIntensity()) / bin_size);
            if (to_bin < bin_count_)
            {
              --histogram[to_bin];
              --elements_in_window;
            }
            ++window_pos_borderleft;
          }
          
          //std::printf("S1: %E %E\n", (*window_pos_borderright).getMZ(), (*window_pos_center).getMZ() + window_half_size);
            
 
          // add all elements to histogram that will enter the window on the RIGHT side
          while (     (window_pos_borderright != scan_last_)
                      && ((*window_pos_borderright).getMZ() < (*window_pos_center).getMZ() + window_half_size )                     )
          {
            //std::printf("Sb: %E %E %E\n", (*window_pos_borderright).getMZ(), (*window_pos_center).getMZ() + window_half_size, (*window_pos_borderright).getMZ() - ((*window_pos_center).getMZ() + window_half_size));
            
            to_bin = (int) (((*window_pos_borderright).getIntensity()) / bin_size);
            if (to_bin < bin_count_)
            {
              ++histogram[to_bin];
              ++elements_in_window;
            }
            ++window_pos_borderright;
          }

          if (elements_in_window < min_required_elements_)
          {
            noise = noise_for_empty_window_;
            ++sparse_window_percent;
          }
          else
          {

            hist_rightmost_bin = bin_count_;

            // do iteration on histogram and find threshold
            for (int i=0;i<3;++i)
            {
              // mean
              hist_mean = 0;
              for (int bin = 0; bin < hist_rightmost_bin; ++bin)
              {
                //std::cout << "V: " << bin << " " << hist_mean << " " << histogram[bin] << " " << elements_in_window << " " << bin_value[bin] << "\n";
                // immediate division is numerically more stable
                hist_mean += histogram[bin] / (double) elements_in_window * bin_value[bin] ;
              }
              //hist_mean = hist_mean / elements_in_window;

              // stdev
              hist_stdev = 0;
              for (int bin = 0; bin < hist_rightmost_bin; ++bin)
              {
                hist_stdev += histogram[bin]/ (double) elements_in_window * std::pow(bin_value[bin]-hist_mean, 2);
              }
              hist_stdev = std::sqrt(hist_stdev);

              //determine new threshold (i.e. the rightmost bin we consider)
              int estimate = (int) ((hist_mean + hist_stdev * stdev_ - 1) / bin_size + 1);
              //std::cout << "E: " << hist_mean << " " << hist_stdev << " " << stdev_ << " " << bin_size<< " " << estimate << "\n";
              hist_rightmost_bin = std::min(estimate, bin_count_);
            }

            // just avoid division by 0
            noise = std::max(1.0, hist_mean);
          }

          // store result
          stn_estimates_[*window_pos_center] = (*window_pos_center).getIntensity() / noise;



          // advance the window center by one datapoint
          ++window_pos_center;
          ++window_count;
          // update progress 
          SignalToNoiseEstimator< Container >::setProgress(window_count);
                    
        } // end while

        SignalToNoiseEstimator< Container >::endProgress();
        
        sparse_window_percent = sparse_window_percent *100 / window_count;
        // warn if percentage of sparse windows is above 20%
        if (sparse_window_percent > 20)
        {
          std::cerr << "WARNING in SignalToNoiseEstimatorMeanIterative: "
          << sparse_window_percent
          << "% of all windows were sparse. You should consider increasing WinLen or increasing MinReqElementsInWindow"
          << " You should also check the MaximalIntensity value (or the parameters for its heuristic estimation)"
          << " If it is too low, then too many high intensity peaks will be discarded, which leads to a sparse window!"
          << std::endl;
        }

        return;
        
      } // end of shiftWindow_


      /// overridden function from DefaultParamHandler to keep members up to date, when a parameter is changed
      void updateMembers_()
      {
        max_intensity_         = (double)param_.getValue("MaxIntensity"); 
        auto_max_stdev_Factor_ = (double)param_.getValue("AutoMaxStdevFactor"); 
        auto_max_percentile_   = param_.getValue("AutoMaxPercentile"); 
        auto_mode_             = param_.getValue("AutoMode"); 
        win_len_               = (double)param_.getValue("WinLen"); 
        bin_count_             = param_.getValue("BinCount"); 
        stdev_                 = (double)param_.getValue("StdevMP"); 
        min_required_elements_ = param_.getValue("MinRequiredElements"); 
        noise_for_empty_window_= (double)param_.getValue("NoiseForEmptyWindow"); 
        is_result_valid_ = false;
      }
  
      /// maximal intensity considered during binning (values above get discarded)
      double max_intensity_;
      /// parameter for initial automatic estimation of "max_intensity_": a stdev multiplier
      double auto_max_stdev_Factor_;
      /// parameter for initial automatic estimation of "max_intensity_" percentile or a stdev
      double auto_max_percentile_;
      /// determines which method shall be used for estimating "max_intensity_". valid are MANUAL=-1, AUTOMAXBYSTDEV=0 or AUTOMAXBYPERCENT=1
      int    auto_mode_;
      /// range of data points which belong to a window in Thomson
      double win_len_;
      /// number of bins in the histogram
      int    bin_count_;
      /// multiplier for the stdev of intensities
      double stdev_;
      /// minimal number of elements a window needs to cover to be used
      int min_required_elements_;
      /// used as noise value for windows which cover less than "min_required_elements_" 
      /// use a very high value if you want to get a low S/N result
      double noise_for_empty_window_;




  };

}// namespace OpenMS

#endif //OPENMS_FILTERING_NOISEESTIMATION_SIGNALTONOISEESTIMATORMEANITERATIVE_H
