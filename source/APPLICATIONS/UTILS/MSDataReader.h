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

#ifndef OPENMS_LOWMEM_COMMON
#define OPENMS_LOWMEM_COMMON

#include <boost/shared_ptr.hpp>
#include "MSDataIOInterface.h"

/*
 * MSDataReaders contain classes that have the same interface as MSExperiment and
 * can be handed into mzMLHandler.
 *
 * These classes are derived from MSDataIOInterface
 *
*/

namespace OpenMS
{

  namespace Internal
  {

    /*
     * Provides the same interface as MSExperiment in order to be passed into mzMLHandler (for example).
     *
     * Only the interface used by the mzMLHandler is implemented here and this
     * class should be used as base class. The child class then implements those
     * functions which it actually expects to be used.
     *
     * Example usage:
     *
          MSDataReader<MzMLConsumer, Peak1D, ChromatogramPeak> exp_reader;
          MzMLFile mz_data_file;
          exp_reader.setConsumer(cacher);
          mz_data_file.load(in, exp_reader);
     *
    */
    template <typename ConsumerT, typename PeakT = Peak1D, typename ChromatogramPeakT = ChromatogramPeak>
    class OPENMS_DLLAPI MSDataReader :
      public Interfaces::MSDataIOInterface<PeakT, ChromatogramPeakT>
    {

  public:

      /// Constructor
      MSDataReader() {}

      // Need implementation, but we do not care about them
      inline void reserve(Size /* s */) {}
      inline void reserveSpaceChromatograms(Size /* s */) {}
      inline void reserveSpaceSpectra(Size /* s */) {}

      /// Resets all internal values
      void reset()
      {
        ExperimentalSettings::operator=(ExperimentalSettings()); //reset meta info
      }

      /// adds a spectrum to the consumer and keeps the meta-data (SpectrumSettings)
      void addSpectrum(const MSSpectrum<PeakT> & spectrum)
      {
        consumer->consumeSpectrum(spectrum);

        // We copy the meta-data of the spectrum
        MSSpectrum<PeakT> cpy = spectrum;
        cpy.clear(false);
        this->spectra_.push_back(cpy);
      }

      /// adds a chromatogram to the consumer and keeps the meta-data (ChromatogramSettings)
      void addChromatogram(const MSChromatogram<ChromatogramPeakT> & chromatogram)
      {
        consumer->consumeChromatogram(chromatogram);

        // We copy the meta-data of the chromatogram
        MSChromatogram<ChromatogramPeakT> cpy = chromatogram;
        cpy.clear(false);
        this->chromatograms_.push_back(cpy);
      }

      /// returns the list with the spectra settings
      const std::vector<MSSpectrum<PeakT> > & getSpectraSettings() const
      {
        return this->spectra_;
      }

      /// returns the list with the chromatogram settings
      const std::vector<MSChromatogram<ChromatogramPeakT> > & getChromatogramSettings() const
      {
        return this->chromatograms_;
      }

      inline void setConsumer(boost::shared_ptr<ConsumerT> c) { consumer = c; }

    protected:
      boost::shared_ptr<ConsumerT> consumer;
    };

    /*
     * Can count how many spectra and chromatograms are present in an mzML file. 
     *
     * This class relies on the fact that each spectrum and chromatogram will be
     * added to a writer through the addSpectrum and addChromatogram calls. The
     * drawback of this method (compared to the MSDataCounterReserve) is
     * that it is much slower. 
     *
     * Example usage:
     *
          MzMLFile f;
          MSDataCounter<Peak1D> exp_cnt;
          f.load(in, exp_cnt);
          // Result contained in exp_cnt.spectraCounts and exp_cnt.chromatogramCounts
     *
    */
    template <typename PeakT = Peak1D, typename ChromatogramPeakT = ChromatogramPeak>
    class OPENMS_DLLAPI MSDataCounter :
      public Interfaces::MSDataIOInterface<PeakT, ChromatogramPeakT>
    {

  public:

      /// Constructor
      MSDataCounter():
          spectraCounts(0),
          chromatogramCounts(0)
      {}

      inline void reserve(Size /* s */) {}
      inline void reserveSpaceChromatograms(Size /* s */) {}
      inline void reserveSpaceSpectra(Size /* s */) {}
      void reset() {}

      void addSpectrum(const MSSpectrum<PeakT> & /* spectrum */) { spectraCounts++; }
      void addChromatogram(const MSChromatogram<ChromatogramPeakT> & /* chromatogram */) { chromatogramCounts++; }

      Size spectraCounts;
      Size chromatogramCounts;
    };

    /*
     * Can count how many spectra and chromatograms are present in an mzML file. 
     *
     * This class relies on the fact that the spectrumList and chromatogramList
     * count attributes are accurate and that the MzMLHandler will try to reserve
     * appropriate space for them. 
     *
     * Example usage:
     *
          MzMLFile f;
          MSDataCounterReserve<Peak1D> exp_cnt;
          f.getOptions().addMSLevel(-1);
          f.load(in, exp_cnt);
          // Result contained in exp_cnt.spectraCounts and exp_cnt.chromatogramCounts
     *
    */
    template <typename PeakT = Peak1D, typename ChromatogramPeakT = ChromatogramPeak>
    class OPENMS_DLLAPI MSDataCounterReserve :
      public Interfaces::MSDataIOInterface<PeakT, ChromatogramPeakT>
    {

  public:

      /// Constructor
      MSDataCounterReserve():
          spectraCounts(0),
          chromatogramCounts(0)
      {}

      // grab the size of the chromatogram/spectra vector from the reserve calls
      inline void reserve(Size s) { spectraCounts = s; }
      inline void reserveSpaceSpectra(Size s) { spectraCounts = s; }
      inline void reserveSpaceChromatograms(Size s) { chromatogramCounts = s; }

      void reset() {}
      void addSpectrum(const MSSpectrum<PeakT> & /* spectrum */) {}
      void addChromatogram(const MSChromatogram<ChromatogramPeakT> & /* chromatogram */) {}

      Size spectraCounts;
      Size chromatogramCounts;
    };

  } // Internal

} // NS OpenMS

#endif

