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

#ifndef OPENMS_FILTERING_NOISEESTIMATION_SIGNALTONOISEESTIMATORMEDIAN_H
#define OPENMS_FILTERING_NOISEESTIMATION_SIGNALTONOISEESTIMATORMEDIAN_H


#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimator.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>

namespace OpenMS
{
  /**
    @brief Estimates the signal/noise (S/N) ratio of each data point in a scan by using the median (histogram based)
   
    For each datapoint in the given scan, we collect a range of data points around it (param: WindowLength).
    The noise for a datapoint is estimated to be the median of the intensities of the current window.
    If the number of elements in the current window is not sufficient (param: MinReqElementsInWindow),
    the noise level is set to a default value (param: NoiseEmptyWindow).
    The whole computation is histogram based, so the user will need to supply a number of bins (param: BinCount), which determines
    the level of error and runtime. The maximal intensity for a datapoint to be included in the histogram can be either determined 
    automatically (params: AutoMaxIntensity, AutoMode) by two different methods or can be set directly by the user (param: MaxIntensity).
    If provided, the MaxIntensity param is always favoured over AutoMaxIntensity&AutoMode. If the (estimated) MaxIntensity value is too
    low and the median is found to be in the last (&highest) bin, a warning to std:err will be given. In this case you should increase
    MaxIntensity (and optionally the BinCount).
    
    The class uses a lazy evaluation approach and will compute ALL S/N values when the first request is made. 
    Changing any of the parameters will invalidate the S/N values (which will invoke a recomputation on the next request).

    Accepted Ini-File Parameters:
    <table>
    <tr><th>Parameter                                               </th><th>  Description                                </th></tr>
    <tr><td> WindowLength          </td><td>  window length in Thomson                   </td></tr>
    <tr><td> BinCount              </td><td>  number of bins used                        </td></tr>
    <tr><td> MinReqElementsInWindow</td><td>  minimum number of elements required in a window   </td></tr>
    <tr><td> NoiseEmptyWindow      </td><td>  noise value used for sparse windows        </td></tr>   
      
    <tr><td> MaxIntensity          </td><td> maximal intensity used for histogram construction. By default, it will be calculated automatically (see AutoMode)<br>
                                                 If this parameter is given, the "AutoMode" is ignored 
                                                 Be aware though that choosing an adverse "MaxIntensity" might lead to bad results
                                                 All intensities EQUAL/ABOVE "MaxIntensity" will not be added to the histogram. If you 
                                                 choose "MaxIntensity" too small, the noise estimate might be too small as well.
                                                 If chosen too big, the bins become quite large (which you could counter by increasing
                                                 "BinCount", which increases runtime). 
                                                                                                        </td></tr>
    <tr><td colspan=2> [The following 3 parameters belong together and should both be provided if <b>MaxIntensity</b> was NOT given]</td></tr>

    <tr><td> AutoMode                 </td><td> method to use to determine "MaxIntensity": 
                                                             <pre>-1           disable and use <u>MaxIntensity</u>
                                                             <pre>0 [default]  use <i>mean</i> + <u>AutoMaxStdevFactor</u> * <i>stdev</i></pre> 
                                                             <pre>1            <u>AutoMaxPercentile</u>th percentile</pre>
                                                                                                                          </td></tr>
    <tr><td> AutoMaxStdevFactor       </td><td> parameter for "MaxIntensity" estimation (if "AutoMode" == 0);<br>default: 3  </td></tr>
                                                              
    <tr><td> AutoMaxPercentile        </td><td> parameter for "MaxIntensity" estimation (if "AutoMode" == 1);<br>default: 95 (do NOT exceed 100)</td></tr>

    </table>   
    
        
   @note 
    Warning to *stderr* if sparse_window_percent > 20
            - percent of windows that have less than "MinReqElementsInWindow" of elements
              (noise estimates in those windows are simply a constant "NoiseEmptyWindow")
            .             
    Warning to *stderr* if histogram_oob_percent (oob=_out_of_bounds) > 1
            - percentage of median estimations that had to rely on the last(=rightmost) bin
              which gives an unreliable result
            .  
    
    
    @ingroup Filtering
    
  */
  
  
  class SignalToNoiseEstimatorMedian : public SignalToNoiseEstimator
  {

  public:

    /// method to use for estimating the maximal intensity that is used for histogram calculation
    enum IntensityThresholdCalculation { MANUAL=-1, AUTOMAXBYSTDEV=0, AUTOMAXBYPERCENT=1 };

    using SignalToNoiseEstimator::stn_estimates_;
    using SignalToNoiseEstimator::first_;
    using SignalToNoiseEstimator::last_;
    using SignalToNoiseEstimator::is_result_valid_;

    /// default constructor
    inline SignalToNoiseEstimatorMedian()
    {
      defaults_.setValue("MaxIntensity", -1); 
      defaults_.setValue("AutoMaxStdevFactor", 3.0); 
      defaults_.setValue("AutoMaxPercentile", 95); 
      defaults_.setValue("AutoMode", 0); 
      defaults_.setValue("WinLen", 200.0); 
      defaults_.setValue("BinCount", 30); 
      defaults_.setValue("MinRequiredElements", 10); 
      defaults_.setValue("NoiseForEmptyWindow", 2.0); 

      defaultsToParam_();
    }


    /// Copy Constructor
    inline SignalToNoiseEstimatorMedian(const SignalToNoiseEstimatorMedian&  source)
        : SignalToNoiseEstimator(source)
    {
      updateMembers_();
    }


    /** @name Assignment
     */
    //@{
    ///
    inline SignalToNoiseEstimatorMedian& operator=(const SignalToNoiseEstimatorMedian& source)
    {
      if(&source == this) return *this; 
      SignalToNoiseEstimator::operator=(source);
      updateMembers_();
      return *this;
    }
    //@}


    /// Destructor
    virtual ~SignalToNoiseEstimatorMedian()
    {}

    /** Accessors
     */
    
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

    /// Return to signal/noise estimate for data point @p data_point
    /// @note the first query to this function will take longer, as
    ///        all SignalToNoise values are calculated
    /// @note you will get a warning to stderr if more than 20% of the 
    ///        noise estimates used sparse windows
    /// @note you will get a warning to stderr if more than 1% of the 
    ///        median estimates was in the rightmost histogram bin
    ///        (you should consider increasing MaxIntensity)
    virtual double getSignalToNoise(PeakIterator data_point)
    {
      double sparse_window = 0;
      double histogram_oob = 0;
      
      if (!is_result_valid_)
      {
        shiftWindow_(first_, last_, sparse_window, histogram_oob);
        
        is_result_valid_ = true;

        // warn if percentage of sparse windows is above 20%
        if (sparse_window > 20) 
        {
          std::cerr << "WARNING in DSignalToNoiseEstimatorMedian: " 
                   << sparse_window 
                   << "% of all windows were sparse. You should consider decreasing WindowLength and/or MinReqElementsInWindow" 
                   << std::endl;
        }
        
        // warn if percentage of possibly wrong median estimates is above 1%
        if (histogram_oob > 1) 
        {
          std::cerr << "WARNING in DSignalToNoiseEstimatorMedian: " 
                   << histogram_oob 
                   << "% of all Signal-to-Noise estimates are too high, because the median was found in the rightmost histogram-bin. " 
                   << "You should consider increasing MaxIntensity (and maybe BinCount with it, to keep bin width reasonable)" 
                   << std::endl;
        }
        
      }  
      
      return stn_estimates_[*data_point];
    }

  protected:


    /// calculate StN values for all datapoints given, by using a sliding window approach
    /// @param scan_first_ first element in the scan
    /// @param scan_last_ last element in the scan (disregarded)
    /// @param sparse_window_percent percent of windows that have less than "min_required_elements_" of elements
    ///          (noise estimates in those windows are simply a constant "noise_for_empty_window_")           
    /// @param histogram_oob_percent (oob=_out_of_bounds) percentage of median estimations that had to rely on the last(=rightmost) bin
    ///          which gives an unreliable result
    void shiftWindow_(const PeakIterator& scan_first_, const PeakIterator& scan_last_, double &sparse_window_percent, double &histogram_oob_percent)
    throw(Exception::InvalidValue)
    {
      // reset counter for sparse windows
      sparse_window_percent = 0;
      // reset counter for histogram overflow
      histogram_oob_percent = 0;
      
      // reset the results
      stn_estimates_.clear();
      
      // maximal range of histogram needs to be calculated first
      if (auto_mode_ == AUTOMAXBYSTDEV)
      {
        // use MEAN+auto_max_intensity_*STDEV as threshold
        GaussianEstimate gauss_global = estimate(scan_first_, scan_last_);
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
        double maxInt = 0;
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
      int bin_count_minus_1 = bin_count_ - 1;
      
      std::vector <int> histogram(bin_count_, 0);
      std::vector <double> bin_value(bin_count_, 0);
      // calculate average intensity that is represented by a bin
      for (int bin=0; bin<bin_count_; bin++)
      {
         histogram[bin] = 0;
         bin_value[bin] = (bin + 0.5) * bin_size;           
      }
      // bin in which a datapoint would fall
      int to_bin;

      // index of bin where the median is located
      int median_bin;
      // additive number of elements from left to x in histogram
      int element_inc_count;
      
      // tracks elements in current window, which may vary because of uneven spaced data
      int elements_in_window = 0;
      // number of windows
      int window_count = 0;
      
      // number of elements where we find the median
      int element_in_window_half;
      
      double noise;    // noise value of a datapoint      

      while (window_pos_center != scan_last_)
      {
        
        // erase all elements from histogram that will leave the window on the LEFT side
        while ( (*window_pos_borderleft).getMZ() <  (*window_pos_center).getMZ() - window_half_size )
        {
          to_bin = std::min((int) (((*window_pos_borderleft).getIntensity()) / bin_size), bin_count_minus_1);
          --histogram[to_bin];
          --elements_in_window;
          ++window_pos_borderleft;
        }
        
        // add all elements to histogram that will enter the window on the RIGHT side
        while (    (window_pos_borderright != scan_last_)
                &&((*window_pos_borderright).getMZ() <= (*window_pos_center).getMZ() + window_half_size ) )
        {
          to_bin = std::min((int) (((*window_pos_borderright).getIntensity()) / bin_size), bin_count_minus_1);
          ++histogram[to_bin];
          ++elements_in_window;
          ++window_pos_borderright;
        }

        if (elements_in_window < min_required_elements_)
        {
          noise = noise_for_empty_window_;
          ++sparse_window_percent;
        }
        else
        {
          // find bin i where ceil[elements_in_window/2] <= sum_c(0..i){ histogram[c] }
          median_bin = -1;
          element_inc_count = 0;
          element_in_window_half = (elements_in_window+1) / 2;
          while (median_bin < bin_count_minus_1 && element_inc_count < element_in_window_half) {
            ++median_bin;
            element_inc_count += histogram[median_bin];
          }

          // increase the error count
          if (median_bin == bin_count_minus_1) {++histogram_oob_percent;}
          
          // just avoid division by 0
          noise = std::max(1.0, bin_value[median_bin]);
        }
        
        // store result
        stn_estimates_[*window_pos_center] = (*window_pos_center).getIntensity() / noise;
        
        
        // advance the window center by one datapoint
        ++window_pos_center;
        ++window_count;  
      } // end while

      sparse_window_percent = sparse_window_percent *100 / window_count;
      histogram_oob_percent = histogram_oob_percent *100 / window_count;
      
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
      min_required_elements_ = param_.getValue("MinRequiredElements"); 
      noise_for_empty_window_= (double)param_.getValue("NoiseForEmptyWindow"); 
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
    /// minimal number of elements a window needs to cover to be used
    int min_required_elements_;
    /// used as noise value for windows which cover less than "min_required_elements_" 
    /// use a very high value if you want to get a low S/N result
    double noise_for_empty_window_;



  };

}// namespace OpenMS

#endif //OPENMS_FILTERING_NOISEESTIMATION_DSIGNALTONOISEESTIMATORMEDIAN_H
