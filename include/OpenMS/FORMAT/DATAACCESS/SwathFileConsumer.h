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

#ifndef OPENMS_FORMAT_DATAACCESS_SWATHFILECONSUMER_H
#define OPENMS_FORMAT_DATAACCESS_SWATHFILECONSUMER_H

#include <boost/cast.hpp>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// Consumers
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
/*
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
*/

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>


#include <OpenMS/ANALYSIS/OPENSWATH/CachedmzML.h>

namespace OpenMS 
{



  ///////////////////////////////////
  // Cached / Regular Swath Map Loaders
  ///////////////////////////////////

  /**
   * @brief Abstract base class which can consume spectra coming from SWATH experiment stored in a single file.
   *
   * The class consumes spectra which are coming from a complete SWATH experiment.
   * It expects each set of SWATH spectra to be separated by an MS1 spectrum and
   * the order of the SWATH spectra to be preserved. For example, the spectra
   * could be arranged in the following fashion:
   *
   * - MS1 Spectrum Precursor = [0,0]
   * - MS2 Spectrum Precursor = [400,425]
   * - MS2 Spectrum Precursor = [425,450]
   * [...]
   * - MS2 Spectrum Precursor = [1175,1200]
   * - MS1 Spectrum Precursor = [0,0]
   * - MS2 Spectrum Precursor = [400,425]
   * - MS2 Spectrum Precursor = [425,450]
   * [...]
   *
   * Base classes are expected to implement functions consuming a spectrum coming
   * from a specific SWATH or an MS1 spectrum and a final function
   * ensureMapsAreFilled_ after which the swath_maps_ vector needs to contain
   * valid pointers to MSExperiment.
   *
   * Usage:
   *
   * FullSwathFileLoader * dataConsumer;
   * // assign dataConsumer to an implementation of FullSwathFileLoader
   * MzMLFile().transform(file, dataConsumer);
   *
   */
  class OPENMS_DLLAPI FullSwathFileLoader :
    public Interfaces::IMSDataConsumer<> 
  {

  public:
    typedef MSExperiment<> MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef MapType::ChromatogramType ChromatogramType;

    FullSwathFileLoader() :
      ms1_counter_(0),
      ms2_counter_(0)
    {}

    ~FullSwathFileLoader() { }

    void setExpectedSize(Size, Size) {}
    void setExperimentalSettings(ExperimentalSettings& exp) {settings_ = exp;}

    /**
     * @brief Populate the vector of swath maps after consuming all spectra. 
     *
     * Will populate the input vector with SwathMap objects which correspond to
     * the MS1 map (if present) and the MS2 maps (SWATH maps). This should be
     * called after all spectra are consumed.
     *
     */
    void retrieveSwathMaps(std::vector< OpenSwath::SwathMap > & maps)
    {
      ensureMapsAreFilled_();
      if (ms1_map_)
      {
        OpenSwath::SwathMap map;
        map.sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(ms1_map_);
        map.lower = -1;
        map.upper = -1;
        map.ms1 = true;
        maps.push_back(map);
      }

      // TODO handle these cases...
      assert(swath_prec_lower_.size() == swath_maps_.size());
      assert(swath_prec_upper_.size() == swath_maps_.size());

      for (Size i = 0; i < swath_maps_.size(); i++)
      {
        OpenSwath::SwathMap map;
        map.sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_maps_[i]);
        map.lower = swath_prec_lower_[i];
        map.upper = swath_prec_upper_[i];
        map.ms1 = false;
        maps.push_back(map);
      }
    }

    /// Consume a chromatogram -> should not happen when dealing with SWATH maps
    void consumeChromatogram(MapType::ChromatogramType &) 
    {
      std::cerr << "Read spectrum while reading SWATH files, did not expect that!" << std::endl;
    }

    /// Consume a spectrum which may belong either to an MS1 scan or one of n MS2 (SWATH) scans
    void consumeSpectrum(MapType::SpectrumType & s)
    {
      if (s.getMSLevel() == 1)
      {
        // append a new MS1 scan, set the ms2 counter to zero and proceed
        consumeMS1Spectrum_(s);
        ms2_counter_ = 0;
        ms1_counter_++;
      }
      else 
      {
        // If this is the first encounter of this SWATH map, try to read the isolation windows
        if (ms2_counter_ == swath_maps_.size())
        {
          if (!s.getPrecursors().empty())
          {
            const std::vector<Precursor> prec = s.getPrecursors();
            double lower = prec[0].getIsolationWindowLowerOffset();
            double upper = prec[0].getIsolationWindowUpperOffset();
            if ( prec[0].getIsolationWindowLowerOffset() > 0.0) swath_prec_lower_.push_back(lower);
            if ( prec[0].getIsolationWindowUpperOffset() > 0.0) swath_prec_upper_.push_back(upper);
            swath_prec_center_.push_back( prec[0].getMZ() );
          }
        }
        else if (ms2_counter_ > swath_prec_center_.size() && ms2_counter_ > swath_prec_lower_.size())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
              "FullSwathFileLoader: MS2 counter is larger than size of swath maps! Are the swath_maps representing the number of read in maps?");
        }
        consumeSwathSpectrum_(s, ms2_counter_);
        ms2_counter_++;
      }
    }

  protected:
    /**
     * @brief Consume an MS2 spectrum belonging to SWATH "swath_nr"
     *
     * This function should handle a Spectrum belonging to a specific SWATH
     * (indicated by swath_nr).
     *
     * @note after this call, swath_maps_.size() _must_ increase by one if
     * ms2_counter_ == swath_maps_.size() (i.e. if a new swath was encountered
     * the first time)
     */
    virtual void consumeSwathSpectrum_(MapType::SpectrumType & s, int swath_nr) = 0;
    /// @brief Consume an MS1 spectrum
    virtual void consumeMS1Spectrum_(MapType::SpectrumType & s) = 0;
    /**
     * @brief Callback function after the reading is complete 
     *
     * Has to ensure that swath_maps_ and ms1_map_ are correctly populated.
     */
    virtual void ensureMapsAreFilled_() = 0;

    size_t ms1_counter_;
    size_t ms2_counter_;

    /// A list of SWATH maps and the MS1 map
    std::vector< boost::shared_ptr<MSExperiment<> > > swath_maps_;
    boost::shared_ptr<MSExperiment<> > ms1_map_;

    /// Values of lower limit, center and upper limit of the isolation windows
    std::vector<double> swath_prec_center_;
    std::vector<double> swath_prec_lower_;
    std::vector<double> swath_prec_upper_;
    /// The Experimental settings 
    // (MSExperiment has no constructor using ExperimentalSettings)
    MSExperiment<> settings_;

  };

  /**
   * @brief In-memory implementation of FullSwathFileLoader 
   *
   * Keeps all the spectra in memory by just appending them to an MSExperiment.
   *
   */
  class OPENMS_DLLAPI RegularSwathFileLoader :
    public FullSwathFileLoader
  {

  public:
    typedef MSExperiment<> MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef MapType::ChromatogramType ChromatogramType;

  protected:
    void addNewSwathMap_()
    {
      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>(settings_));
      swath_maps_.push_back(exp);
    }
    void consumeSwathSpectrum_(MapType::SpectrumType & s, int swath_nr)
    {
      if (swath_nr == (int)swath_maps_.size() )
        addNewSwathMap_();
      swath_maps_[swath_nr]->addSpectrum(s);
    }

    void addMS1Map_()
    {
      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>(settings_));
      ms1_map_ = exp;
    }
    void consumeMS1Spectrum_(MapType::SpectrumType & s)
    {
      if (! ms1_map_ )
        addMS1Map_();
      ms1_map_->addSpectrum(s);
    }

    void ensureMapsAreFilled_() {}
  };

  /**
   * @brief On-disked cached implementation of FullSwathFileLoader 
   *
   * Writes all spectra immediately to disk in a user-specified caching
   * location using the CachedMzMLConsumer. Internally, it handles 
   * n+1 (n SWATH + 1 MS1 map) CachedMzMLConsumers which can consume the
   * spectra and write them to disk immediately.
   *
   */
  class OPENMS_DLLAPI CachedSwathFileLoader :
    public FullSwathFileLoader
  {

  public:
    typedef MSExperiment<> MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef MapType::ChromatogramType ChromatogramType;

    CachedSwathFileLoader(String cachedir, String basename, Size nr_ms1_spectra, std::vector<int> nr_ms2_spectra) :
      ms1_consumer_(NULL),
      swath_consumers_(),
      cachedir_(cachedir),
      basename_(basename),
      nr_ms1_spectra_(nr_ms1_spectra),
      nr_ms2_spectra_(nr_ms2_spectra)
    {}

    ~CachedSwathFileLoader() 
    { 
      // Properly delete the CachedMzMLConsumers -> free memory and _close_ filestream
      while(!swath_consumers_.empty()) {delete swath_consumers_.back(); swath_consumers_.pop_back();}
      if (ms1_consumer_ != NULL) delete ms1_consumer_;
    }

  protected:
    void addNewSwathMap_()
    {
      String meta_file = cachedir_ + basename_ + "_" + String(swath_consumers_.size()) +  ".mzML";
      String cached_file = meta_file + ".cached";
      CachedMzMLConsumer * consumer = new CachedMzMLConsumer(cached_file, true);
      consumer->setExpectedSize(nr_ms2_spectra_[swath_consumers_.size()], 0);
      swath_consumers_.push_back(consumer);

      // maps for meta data
      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>(settings_));
      swath_maps_.push_back(exp);
    }
    void consumeSwathSpectrum_(MapType::SpectrumType & s, int swath_nr)
    {
      if (swath_nr == (int)swath_consumers_.size() )
        addNewSwathMap_();

      swath_consumers_[swath_nr]->consumeSpectrum(s);
      swath_maps_[swath_nr]->addSpectrum(s); // append for the metadata (actual data is deleted)
    }

    void addMS1Map_()
    {
      String meta_file = cachedir_ + basename_ + "_ms1.mzML";
      String cached_file = meta_file + ".cached";
      ms1_consumer_ = new CachedMzMLConsumer(cached_file, true);
      ms1_consumer_->setExpectedSize(nr_ms1_spectra_, 0);
      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>(settings_));
      ms1_map_ = exp;
    }
    void consumeMS1Spectrum_(MapType::SpectrumType & s)
    {
      if (ms1_consumer_ == NULL)
        addMS1Map_();
      ms1_consumer_->consumeSpectrum(s);
      ms1_map_->addSpectrum(s); // append for the metadata (actual data is deleted)
    }

    void ensureMapsAreFilled_() 
    {
      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
      String meta_file = cachedir_ + basename_ + "_ms1.mzML";
      // write metadata to disk and store the correct data processing tag
      CachedmzML().writeMetadata(*ms1_map_, meta_file, true);
      MzMLFile().load(meta_file, *exp.get());
      ms1_map_ = exp;

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (Size i = 0; i < swath_consumers_.size(); i++)
      {
        boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
        String meta_file = cachedir_ + basename_ + "_" + String(i) +  ".mzML";
        // write metadata to disk and store the correct data processing tag
        CachedmzML().writeMetadata(*swath_maps_[i], meta_file, true);
        MzMLFile().load(meta_file, *exp.get());
        swath_maps_[i] = exp;
      }
    }

    CachedMzMLConsumer * ms1_consumer_;
    std::vector< CachedMzMLConsumer * > swath_consumers_;

    String cachedir_;
    String basename_;
    int nr_ms1_spectra_;
    std::vector<int> nr_ms2_spectra_;
  };

  /*
  ///////////////////////////////////
  // Data Reducers
  ///////////////////////////////////
  class OPENMS_DLLAPI DataReducer :
    public MSDataTransformingConsumer 
  {

  public:
    DataReducer(GaussFilter nf, PeakPickerHiRes pp) :
      pp_(pp), nf_(nf) {}

    void consumeSpectrum(MapType::SpectrumType & s)
    {
      MapType::SpectrumType sout;
      nf_.filter(s);
      pp_.pick(s, sout);
      s = sout;
    }

    PeakPickerHiRes pp_;
    GaussFilter nf_;
  };

  class OPENMS_DLLAPI DataReducerIterative :
    public MSDataTransformingConsumer 
  {

  public:
    DataReducerIterative(GaussFilter nf, PeakPickerIterative pp) :
      pp_(pp), nf_(nf) {}

    void consumeSpectrum(MapType::SpectrumType & s)
    {
      MapType::SpectrumType sout;
      nf_.filter(s);
      pp_.pick(s, sout);
      s = sout;
    }

    PeakPickerIterative pp_;
    GaussFilter nf_;
  };
  */

  ///////////////////////////////////
  // Map loader
  ///////////////////////////////////
  class OPENMS_DLLAPI SwathMapLoader :
    public ProgressLogger
  {
    public:

    /// Cache a file to disk
    OpenSwath::SpectrumAccessPtr doCacheFile(String in, String tmp, String tmp_fname, 
        boost::shared_ptr<MSExperiment<Peak1D> > experiment_metadata )
    {
      String cached_file = tmp + tmp_fname + ".cached";
      String meta_file = tmp + tmp_fname;

      // Create new consumer, transform infile, write out metadata
      CachedMzMLConsumer cachedConsumer(cached_file, true);
      MzMLFile().transform(in, &cachedConsumer, *experiment_metadata.get());
      CachedmzML().writeMetadata(*experiment_metadata.get(), meta_file, true);

      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
      MzMLFile().load(meta_file, *exp.get());
      return SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
  }

    /// only read the meta data from a file and use it to populate exp_meta
    void populateMetaData_(String file, boost::shared_ptr<ExperimentalSettings>& exp_meta)
    {
      MSExperiment<Peak1D> tmp;
      MSDataTransformingConsumer c;
      MzMLFile().transform(file, &c, tmp);
      *exp_meta = tmp;
    }

    std::vector< OpenSwath::SwathMap > load_files(StringList file_list, String tmp, 
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions="normal")
    {
      int progress = 0;
      startProgress(0, file_list.size(), "Loading data");

      std::vector< OpenSwath::SwathMap > swath_maps;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(file_list.size()); ++i)
      {
        std::cout << "Loading file " << file_list[i] << std::endl;
        String tmp_fname = "openswath_tmpfile_" + String(i) + ".mzML";

        boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
        OpenSwath::SpectrumAccessPtr spectra_ptr;

        // Populate meta-data
        if (i == 0) { populateMetaData_(file_list[i], exp_meta); }

        if (readoptions == "normal")
        {
          MzMLFile().load(file_list[i], *exp.get());
          spectra_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
        }
        else if (readoptions == "cache")
        {
          // Cache and load the exp (metadata only) file again
          spectra_ptr = doCacheFile(file_list[i], tmp, tmp_fname, exp);
        }
        /*
        else if (readoptions == "reduce")
        {
          GaussFilter gf;
          PeakPickerHiRes pp;

          Param p = gf.getParameters();
          p.setValue("use_ppm_tolerance", "true");
          p.setValue("ppm_tolerance", 50.0);
          gf.setParameters(p);

          // using the consumer to reduce the input data
          DataReducer dataConsumer(gf, pp);
          MzMLFile().transform(file_list[i], &dataConsumer, *exp.get());
          // ownership is transferred to AccessPtr
          spectra_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
        }
        else if (readoptions == "reduce_iterative")
        {
          GaussFilter gf;
          PeakPickerIterative pp;

          Param p = gf.getParameters();
          p.setValue("use_ppm_tolerance", "true");
          p.setValue("ppm_tolerance", 10.0);
          gf.setParameters(p);

          p = pp.getParameters();
          p.setValue("peak_width", 0.04);
          p.setValue("spacing_difference", 2.5);
          p.setValue("signal_to_noise_", 0.0);
          p.setValue("check_width_internally", "true");
          p.setValue("clear_meta_data", "true");
          pp.setParameters(p);

          // using the consumer to reduce the input data
          DataReducerIterative dataConsumer(gf, pp);
          MzMLFile().transform(file_list[i], &dataConsumer, *exp.get());
          // ownership is transferred to AccessPtr
          spectra_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
        }
        */
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
              "Unknown option " + readoptions);
        }

        OpenSwath::SwathMap swath_map;

        bool ms1 = false;
        double upper = -1, lower = -1;
        if (exp->size() == 0)
        {
          std::cerr << "WARNING: File " << file_list[i] << "\n does not have any scans - I will skip it" << std::endl;
          continue;
        }
        if (exp->getSpectra()[0].getPrecursors().size() == 0)
        {
          std::cout << "NOTE: File " << file_list[i] << "\n does not have any precursors - I will assume it is the MS1 scan." << std::endl;
          ms1 = true;
        }
        else
        {
          // Checks that this is really a SWATH map and extracts upper/lower window
          OpenSwathHelper::checkSwathMap(*exp.get(), lower, upper);
        }

        swath_map.sptr = spectra_ptr;
        swath_map.lower = lower;
        swath_map.upper = upper;
        swath_map.ms1 = ms1;
#ifdef _OPENMP
#pragma omp critical (load_files)
#endif
        {
          swath_maps.push_back( swath_map );
          setProgress(progress++);
        }
      }
      endProgress();
      return swath_maps;
    }

    std::vector< OpenSwath::SwathMap > load_files_from_single(String file, String tmp, 
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions="normal")
    {
      startProgress(0, 1, "Loading data file " + file);
      std::vector< OpenSwath::SwathMap > swath_maps;
      FullSwathFileLoader * dataConsumer;
      String tmp_fname = "openswath_tmpfile";

      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
      OpenSwath::SpectrumAccessPtr spectra_ptr;
      OpenSwath::SwathMap swath_map;

      populateMetaData_(file, exp_meta); 

      if (readoptions == "normal")
      {
        dataConsumer = new RegularSwathFileLoader();
        MzMLFile().transform(file, dataConsumer, *exp.get());
      }
      else if (readoptions == "cache")
      {

        std::cout << "Will analyze the metadata first to determine the number of SWATH windows and the window sizes." << std::endl;
        boost::shared_ptr<MSExperiment<Peak1D> > experiment_metadata(new MSExperiment<Peak1D>);
        // First pass through the file -> get the meta data
        {
          MzMLFile f;
          f.getOptions().setAlwaysAppendData(true);
          f.getOptions().setFillData(false);
          f.load(file, *experiment_metadata);
        }

        std::vector<int> swath_counter;
        int nr_ms1_spectra;
        analyzeFullSwath(experiment_metadata->getSpectra(), swath_counter, nr_ms1_spectra);

        std::cout << "Determined there to be " << swath_counter.size() << " SWATH windows and in total " << nr_ms1_spectra << " MS1 spectra" << std::endl;
        dataConsumer = new CachedSwathFileLoader(tmp, tmp_fname, nr_ms1_spectra, swath_counter);
        MzMLFile().transform(file, dataConsumer, *exp.get());
      }
      //else if (readoptions == "reduce") { }
      //else if (readoptions == "reduce_iterative") { }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "Unknown or unsupported option " + readoptions);
      }
      dataConsumer->retrieveSwathMaps(swath_maps);
      delete dataConsumer;

      endProgress();
      return swath_maps;
    }

#ifndef MZXMLSUPPORT
    std::vector< OpenSwath::SwathMap > load_files_from_single_mzxml(String, String, 
      boost::shared_ptr<ExperimentalSettings>&, String)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
          "MzXML not supported");
    }
#else
    std::vector< OpenSwath::SwathMap > load_files_from_single_mzxml(String file, String tmp, 
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions="normal")
    {
      startProgress(0, 1, "Loading data file " + file);
      std::vector< OpenSwath::SwathMap > swath_maps;
      boost::shared_ptr<FullSwathFileLoader> dataConsumer;
      String tmp_fname = "openswath_tmpfile";

      boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
      OpenSwath::SpectrumAccessPtr spectra_ptr;
      OpenSwath::SwathMap swath_map;

      if (readoptions == "normal")
      {
        dataConsumer = boost::shared_ptr<RegularSwathFileLoader>( new RegularSwathFileLoader() ) ; 
        Internal::MSMzXMLDataReader<FullSwathFileLoader> datareader;
        datareader.setConsumer(dataConsumer);
        MzXMLFile().load(file, datareader);
        *exp_meta = datareader;
      }
      else if (readoptions == "cache")
      {

        // First pass through the file -> get the meta data
        std::cout << "Will analyze the metadata first to determine the number of SWATH windows and the window sizes." << std::endl;
        boost::shared_ptr<MSExperiment<Peak1D> > experiment_metadata(new MSExperiment<Peak1D>);
        std::vector<int> swath_counter;
        int nr_ms1_spectra;
        {
          boost::shared_ptr<MSDataTransformingConsumer> noopConsumer = 
            boost::shared_ptr<MSDataTransformingConsumer>( new MSDataTransformingConsumer() ) ; 
          Internal::MSMzXMLDataReader<MSDataTransformingConsumer> datareader;
          datareader.setConsumer(noopConsumer);
          MzXMLFile f;
          f.getOptions().setFillData(false);
          f.load(file, datareader);
          analyzeFullSwath(datareader.getRealSpectra(), swath_counter, nr_ms1_spectra);
          *exp_meta = datareader;
        }

        std::cout << "Determined there to be " << swath_counter.size() << " SWATH windows and in total " << nr_ms1_spectra << " MS1 spectra" << std::endl;
        dataConsumer = boost::shared_ptr<CachedSwathFileLoader>( new CachedSwathFileLoader(tmp, tmp_fname, nr_ms1_spectra, swath_counter) ) ; 
        Internal::MSMzXMLDataReader<FullSwathFileLoader> datareader;
        datareader.setConsumer(dataConsumer);
        MzXMLFile().load(file, datareader);
      }
      //else if (readoptions == "reduce") { }
      //else if (readoptions == "reduce_iterative") { }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "Unknown or unsupported option " + readoptions);
      }
      dataConsumer->retrieveSwathMaps(swath_maps);

      endProgress();
      return swath_maps;
    }
#endif

    void analyzeFullSwath(const std::vector<MSSpectrum<> > exp, std::vector<int> & swath_counter_, int & nr_ms1_spectra)
    {
      int ms1_counter_ = 0;
      int ms2_counter_ = 0;
      for (Size i = 0; i < exp.size(); i++)
      {
        const MSSpectrum<> & s = exp[i];
        {
          if (s.getMSLevel() == 1)
          {
            ms2_counter_ = 0;
            ms1_counter_++;
          }
          else 
          {
            if (ms2_counter_ == (int)swath_counter_.size())
            {
              swath_counter_.push_back(0);
            }
            swath_counter_[ms2_counter_]++;
            ms2_counter_++;
          }
        }
      }
      nr_ms1_spectra = ms1_counter_;
    }

  }; // SwathMapLoader

  void selectChrom_(const MSChromatogram<ChromatogramPeak>& chromatogram_old, 
    MSSpectrum<ChromatogramPeak>& chromatogram, double rt_extraction_window, double center_rt)
  {
    double rt_max = center_rt + rt_extraction_window;
    double rt_min = center_rt - rt_extraction_window;
    for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
    {
      if (rt_extraction_window >= 0 && (it->getRT() < rt_min || it->getRT() > rt_max))
      {
        continue;
      }
      ChromatogramPeak peak;
      peak.setMZ(it->getRT());
      peak.setIntensity(it->getIntensity());
      chromatogram.push_back(peak);
    }
    // TODO should we always warn?
    /*
    if (chromatogram.empty())
    {
      std::cerr << "Error: Could not find any points for chromatogram " + chromatogram.getNativeID() + \
      ". Maybe your retention time transformation is off?" << std::endl;
    }
    */
  }

  /*
  // TODO shared code!! -> OpenSwathRTNormalizer...
  void simple_find_best_feature(OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, 
      std::vector<std::pair<double, double> > & pairs, std::map<OpenMS::String, double> PeptideRTMap)
  {
    for (OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType::iterator trgroup_it = transition_group_map.begin();
        trgroup_it != transition_group_map.end(); trgroup_it++)
    {
      // we need at least one feature to find the best one
      OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType * transition_group = &trgroup_it->second;
      if (transition_group->getFeatures().size() == 0)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "RT normalization: did not find any features for group " + transition_group->getTransitionGroupID());
      }

      // Find the feature with the highest score
      double bestRT = -1;
      double highest_score = -1000;
      for (std::vector<MRMFeature>::iterator mrmfeature = transition_group->getFeaturesMuteable().begin();
           mrmfeature != transition_group->getFeaturesMuteable().end(); mrmfeature++)
      {
        if (mrmfeature->getOverallQuality() > highest_score)
        {
          bestRT = mrmfeature->getRT();
          highest_score = mrmfeature->getOverallQuality();
        }
      }
      String pepref = trgroup_it->second.getTransitions()[0].getPeptideRef();
      pairs.push_back(std::make_pair(bestRT, PeptideRTMap[pepref]));
    }
  }
  */

}

#endif

