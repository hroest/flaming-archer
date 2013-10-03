// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_CHROMATOGRAMEXTRACTORALGORITHM_H
#define OPENMS_ANALYSIS_OPENSWATH_CHROMATOGRAMEXTRACTORALGORITHM_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

namespace OpenMS
{

  /**
  @brief The ChromatogramExtractor extracts chromatograms from a mzML file.

  It will take as input a set of (TraML) transitions and will extract the
  signal of the provided map at the product ion m/z values specified by the
  transitions. The map is thus assumed to be an MS2 map from a SWATH / DIA
  experiment.

  */
  class OPENMS_DLLAPI ChromatogramExtractorAlgorithm :
    public ProgressLogger
  {

public:

    struct ExtractionCoordinates
    {
      double mz; /// mz around which should be extracted
      double rt; /// rt around which should be extracted
      std::string id; /// identifier

      static bool SortExtractionCoordinatesByMZ(
          const ChromatogramExtractorAlgorithm::ExtractionCoordinates& left,
          const ChromatogramExtractorAlgorithm::ExtractionCoordinates& right)
      {
        return left.mz < right.mz;
      }
      static bool SortExtractionCoordinatesReverseByMZ(
          const ChromatogramExtractorAlgorithm::ExtractionCoordinates& left,
          const ChromatogramExtractorAlgorithm::ExtractionCoordinates& right)
      {
        return left.mz > right.mz;
      }
    };

    /**
     * @brief Extract chromatograms at the m/z and RT defined by the ExtractionCoordinates.
     *
     * @param mz_extraction_window Extracts a window of this size in m/z
     * dimension (e.g. a window of 50 ppm means an extraction of 25 ppm on
     * either side)
     * @param rt_extraction_window Extracts a window of this size in RT
     * dimension (e.g. a window of 600 seconds means an extraction of 300
     * seconds on either side)
     *
    */
    void extractChromatograms(const OpenSwath::SpectrumAccessPtr input, 
        std::vector< OpenSwath::ChromatogramPtr >& output, 
        std::vector<ExtractionCoordinates> extraction_coordinates, double& mz_extraction_window,
        bool ppm, double rt_extraction_window, String filter);

public:

    /**
     * @brief Extract the next mz value and add the integrated intensity to integrated_intensity. 
     *
     * This function will sum up all intensities within mz +/-
     * mz_extract_window / 2.0 and add the result to integrated_intensity.
     *
     * @param mz_extraction_window Extracts a window of this size in m/z
     * dimension (e.g. a window of 50 ppm means an extraction of 25 ppm on
     * either side)
     *
     * @note It will change the position of the iterators mz_it and int_it and
     * it can *not* extract any data if the mz-iterator is already passed the
     * mz value given. It is thus critically important to provide all mz values
     * to be extracted in ascending order!
     *
    */
    void extract_value_tophat(const std::vector<double>::const_iterator& mz_start, std::vector<double>::const_iterator& mz_it,
                              const std::vector<double>::const_iterator& mz_end, std::vector<double>::const_iterator& int_it,
                              const double& mz, double& integrated_intensity, double& mz_extraction_window, bool ppm);

private:

    int get_filter_nr(String filter);

  };

}

#endif
