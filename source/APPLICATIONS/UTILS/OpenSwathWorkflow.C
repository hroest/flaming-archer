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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/ANALYSIS/OPENSWATH/CachedmzML.h>

// files
#include <OpenMS/FORMAT/TraMLFile.h>
// #include <OpenMS/FORMAT/CachedMzMLFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>

// helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

// interfaces
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSCached.h>

// algos
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerMRM.h>

#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerIterative.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

using namespace std;
using namespace OpenMS;

#define DEBUG_OPENSWATHWORKFLOW

#ifdef _OPENMP
  #define IF_MASTERTHREAD if (omp_get_thread_num() ==0)  
#else
  #define IF_MASTERTHREAD 
#endif    

namespace OpenMS 
{

  class OPENMS_DLLAPI DataReducer :
    public MSDataTransformingConsumer 
  {

  public:
    DataReducer(GaussFilter nf, PeakPickerHiRes pp) :
      pp_(pp), nf_(nf) {}

    void consumeSpectrum(typename MapType::SpectrumType & s)
    {
      typename MapType::SpectrumType sout;
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

    void consumeSpectrum(typename MapType::SpectrumType & s)
    {
      typename MapType::SpectrumType sout;
      nf_.filter(s);
      pp_.pick(s, sout);
      s = sout;
    }

    PeakPickerIterative pp_;
    GaussFilter nf_;
  };

  class OPENMS_DLLAPI SwathMapLoader :
    public ProgressLogger
  {
    public:
    struct SwathMap
    {
        OpenSwath::SpectrumAccessPtr sptr;
        double lower;
        double upper;
        bool ms1;
    };

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

    std::vector< SwathMap > load_files(StringList file_list, String tmp, String readoptions="normal")
    {
      int progress = 0;
      startProgress(0, file_list.size(), "Loading data");

      std::vector< SwathMap > swath_maps;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(file_list.size()); ++i)
      {
        std::cout << "Loading file " << file_list[i] << std::endl;
        String tmp_fname = "openswath_tmpfile_" + String(i) + ".mzML";

        boost::shared_ptr<MSExperiment<Peak1D> > exp(new MSExperiment<Peak1D>);
        OpenSwath::SpectrumAccessPtr spectra_ptr;

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
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
              "Unknown option " + readoptions);
        }

        SwathMap swath_map;

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
  };

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
      if (chromatogram.empty())
      {
        std::cerr << "Error: Could not find any points for chromatogram " + chromatogram.getNativeID() + \
        ". Maybe your retention time transformation is off?" << std::endl;
      }
  }

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
}

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_OpenSwathWorkflow Workflow

  @brief Complete workflow to run OpenSWATH

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathWorkflow 
  : public TOPPBase
{
public:

  TOPPOpenSwathWorkflow() 
    : TOPPBase("OpenSwathWorkflow", "Complete workflow to run OpenSWATH", true)
  {
  }

protected:

    typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram; // this is the type in which we store the chromatograms for this analysis
    typedef OpenSwath::LightTransition TransitionType;
    typedef SwathMapLoader::SwathMap SwathMap;
    typedef OpenSwath::LightTargetedExperiment TargetedExpType;
    typedef OpenSwath::LightPeptide PeptideType;
    typedef OpenSwath::LightProtein ProteinType;
    typedef OpenSwath::LightModification ModificationType;
    typedef MRMTransitionGroup<MSSpectrum <ChromatogramPeak>, TransitionType> MRMTransitionGroupType; // a transition group holds the MSSpectra with the Chromatogram peaks from above
    typedef std::map<String, MRMTransitionGroupType> TransitionGroupMapType;

    typedef OpenSwath::ChromatogramPtr SmallChromatogram;
    typedef std::map<String, std::vector<const ReactionMonitoringTransition*> > PeptideTransitionMapType;
    // TODO duplicate code!
    // prepare the extraction coordinates for extraction
    void prepare_coordinates(std::vector< SmallChromatogram > & small_chromatograms,
      std::vector< ChromatogramExtractor::ExtractionCoordinates > & coordinates,
      OpenMS::TargetedExperiment & transition_exp_used,
      TransformationDescription trafo, DoubleReal rt_extraction_window, bool ms1)
    {
      trafo.invert(); // copy

      // hash of the peptide reference containing all transitions
      PeptideTransitionMapType peptide_trans_map;
      for (Size i = 0; i < transition_exp_used.getTransitions().size(); i++)
      {
        peptide_trans_map[transition_exp_used.getTransitions()[i].getPeptideRef()].push_back(&transition_exp_used.getTransitions()[i]);
      }

      Size itersize;
      if (ms1) {itersize = transition_exp_used.getPeptides().size();}
      else     {itersize = transition_exp_used.getTransitions().size();}

      for (Size i = 0; i < itersize; i++)
      {
        OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
        small_chromatograms.push_back(s);

        ChromatogramExtractor::ExtractionCoordinates coord;
        TargetedExperiment::Peptide pep;
        OpenMS::ReactionMonitoringTransition transition;

        if (ms1) 
        {
          pep = transition_exp_used.getPeptides()[i];
          transition = (*peptide_trans_map[pep.id][0]);
          coord.mz = transition.getPrecursorMZ();
          coord.id = pep.id;
        }
        else 
        {
          transition = transition_exp_used.getTransitions()[i];
          pep = transition_exp_used.getPeptideByRef(transition.getPeptideRef()); 
          coord.mz = transition.getProductMZ();
          coord.id = transition.getNativeID();
        }

        if (pep.rts.empty() || pep.rts[0].getCVTerms()["MS:1000896"].empty())
        {
          // we dont have retention times -> this is only a problem if we actually
          // wanted to use the RT limit feature.
          if (rt_extraction_window >= 0)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
              "Error: Peptide " + pep.id + " does not have normalized retention times (term 1000896) which are necessary to perform an RT-limited extraction");
          }
          coord.rt = -1;
        }
        else
        {
          coord.rt = pep.rts[0].getCVTerms()["MS:1000896"][0].getValue().toString().toDouble();
          coord.rt = trafo.apply(coord.rt); // apply RT transformation if necessary
        }
        coordinates.push_back(coord);
      }
      std::sort(coordinates.begin(), coordinates.end(), ChromatogramExtractor::ExtractionCoordinates::SortExtractionCoordinatesByMZ);
    }

  void scoreAll_(OpenSwath::SpectrumAccessPtr input,
         OpenSwath::SpectrumAccessPtr swath_map,
         TargetedExpType& transition_exp, 
         TransformationDescription trafo, double rt_extraction_window, 
         FeatureMap<Feature>& output, Param& feature_finder_param)
  {
    double expected_rt;
    TransformationDescription trafo_inv = trafo;
    trafo_inv.invert();

    MRMFeatureFinderScoring featureFinder;
    MRMTransitionGroupPicker trgroup_picker;

    trgroup_picker.setParameters(feature_finder_param.copy("TransitionGroupPicker:", true));
    featureFinder.setParameters(feature_finder_param);
    featureFinder.prepareProteinPeptideMaps_(transition_exp);

    std::map<String, int> chromatogram_map;
    //Size nr_chromatograms = input->getNrChromatograms();
    for (Size i = 0; i < input->getNrChromatograms(); i++)
    {
      chromatogram_map[input->getChromatogramNativeID(i)] = boost::numeric_cast<int>(i);
    }

    // map peptides
    std::map<String, int> assay_peptide_map;
    for (Size i = 0; i < transition_exp.getPeptides().size(); i++)
    {
      // Map peptide id
      assay_peptide_map[transition_exp.getPeptides()[i].id] = boost::numeric_cast<int>(i);
    }

    // Group transitions
    typedef std::map<String, std::vector< const TransitionType* > > AssayMapT;
    AssayMapT assay_map;
    for (Size i = 0; i < transition_exp.getTransitions().size(); i++)
    {
      assay_map[transition_exp.getTransitions()[i].getPeptideRef()].push_back(&transition_exp.getTransitions()[i]);
    }

    // Iterating over all the assays
    for (AssayMapT::iterator assay_it = assay_map.begin(); assay_it != assay_map.end(); assay_it++)
    {
      String id = assay_it->first;

      // Create new transition group if there is none for this peptide
      MRMTransitionGroupType transition_group;
      transition_group.setTransitionGroupID(id);

      expected_rt = transition_exp.getPeptides()[ assay_peptide_map[id] ].rt;

      // Go through all transitions
      for (Size i = 0; i < assay_it->second.size(); i++)
      {
        const TransitionType* transition = assay_it->second[i];

        if (chromatogram_map.find(transition->getNativeID()) == chromatogram_map.end() )
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
              "Error, did not find chromaotgram for transitions" + transition->getNativeID() );
        }

        OpenSwath::ChromatogramPtr cptr = input->getChromatogramById(chromatogram_map[transition->getNativeID()]);
        MSChromatogram<ChromatogramPeak> chromatogram_old;
        OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chromatogram_old, cptr);
        RichPeakChromatogram chromatogram;

        // expected_rt = PeptideRefMap_[transition->getPeptideRef()]->rt;
        chromatogram.setMetaValue("product_mz", transition->getProductMZ());
        chromatogram.setMetaValue("precursor_mz", transition->getPrecursorMZ());
        chromatogram.setNativeID(transition->getNativeID());
        double de_normalized_experimental_rt = trafo_inv.apply(expected_rt);
        selectChrom_(chromatogram_old, chromatogram, rt_extraction_window, de_normalized_experimental_rt);

        // Now add the transition and the chromatogram to the group
        transition_group.addTransition(*transition, transition->getNativeID());
        transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());
      }

      // Process the transition_group
      trgroup_picker.pickTransitionGroup(transition_group);
      featureFinder.scorePeakgroups(transition_group, trafo, swath_map, output);
    }
  }
      
  TransformationDescription RTNormalization(TargetedExperiment transition_exp_,
          std::vector< OpenMS::MSChromatogram<> > chromatograms, double min_rsq, double min_coverage, 
          Param& feature_finder_param)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    int progress = 0;
    progresslogger.startProgress(0, 1, "Retention time normalization");

    OpenSwath::LightTargetedExperiment targeted_exp;
    OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, targeted_exp);

    std::vector<std::pair<double, double> > pairs; // store the RT pairs to write the output trafoXML

    // Store the peptide retention times in an intermediate map
    std::map<OpenMS::String, double> PeptideRTMap;
    for (Size i = 0; i < targeted_exp.getPeptides().size(); i++)
    {
      PeptideRTMap[targeted_exp.getPeptides()[i].id] = targeted_exp.getPeptides()[i].rt; 
    }

    OpenSwath::LightTargetedExperiment transition_exp_used = targeted_exp;

    MRMFeatureFinderScoring featureFinder;
    feature_finder_param.setValue("Scores:use_rt_score", "false");
    feature_finder_param.setValue("Scores:use_elution_model_score", "false");
    feature_finder_param.setValue("Scores:use_elution_model_score", "false");
    feature_finder_param.setValue("rt_extraction_window", -1.0);
    featureFinder.setParameters(feature_finder_param);
    
    FeatureMap<> featureFile; // also for results
    OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map; // for results
    boost::shared_ptr<MSExperiment<Peak1D> > swath_map(new MSExperiment<Peak1D>);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

    boost::shared_ptr<MSExperiment<Peak1D> > xic_map(new MSExperiment<Peak1D>); // the map with the extracted ion chromatograms
    xic_map->setChromatograms(chromatograms);
    OpenSwath::SpectrumAccessPtr chromatogram_ptr = OpenSwath::SpectrumAccessPtr(new OpenMS::SpectrumAccessOpenMS(xic_map));
    TransformationDescription empty_trafo;

    featureFinder.setStrictFlag(false); // TODO remove this, it should be strict (e.g. all transitions need to be present for RT norm)
    featureFinder.pickExperiment(chromatogram_ptr, featureFile, transition_exp_used, empty_trafo, swath_ptr, transition_group_map);

    // find best feature, compute pairs of iRT and real RT
    simple_find_best_feature(transition_group_map, pairs, PeptideRTMap);

    std::vector<std::pair<double, double> > pairs_corrected;
    pairs_corrected = MRMRTNormalizer::rm_outliers(pairs, min_rsq, min_coverage);

    // store transformation, using a linear model as default
    TransformationDescription trafo_out;
    trafo_out.setDataPoints(pairs_corrected);
    Param model_params;
    model_params.setValue("symmetric_regression", "false");
    String model_type = "linear";
    trafo_out.fitModel(model_type, model_params);

    progresslogger.endProgress();
    return trafo_out;
  }

   struct ChromExtractParams 
   {
     double min_upper_edge_dist;
     double extraction_window;
     bool ppm; 
     double rt_extraction_window; 
     String extraction_function;
   };
    
  void simpleExtractChromatograms(const std::vector< SwathMap > & swath_maps,
    const OpenMS::TargetedExperiment & irt_transitions, 
    std::vector< OpenMS::MSChromatogram<> > & chromatograms, ChromExtractParams cp)
  {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(Size i = 0; i < swath_maps.size(); i++)
    {
      if (swath_maps[i].ms1) {continue;}
      TargetedExperiment transition_exp_used;
      OpenSwathHelper::selectSwathTransitions(irt_transitions, transition_exp_used,
          cp.min_upper_edge_dist, swath_maps[i].lower, swath_maps[i].upper);
      if (transition_exp_used.getTransitions().size() == 0) { continue;}

      std::vector< OpenSwath::ChromatogramPtr > tmp_out;
      std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
      prepare_coordinates(tmp_out, coordinates, transition_exp_used,  TransformationDescription(), cp.rt_extraction_window, false);

      ChromatogramExtractor().extractChromatograms(swath_maps[i].sptr, tmp_out, coordinates, cp.extraction_window,
          cp.ppm, cp.rt_extraction_window, cp.extraction_function);
#ifdef _OPENMP
#pragma omp critical (featureFinder)
#endif
      {
        for (Size i = 0; i < tmp_out.size(); i++)
        { 
          OpenMS::MSChromatogram<> chrom;
          OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chrom, tmp_out[i]);
          chrom.setNativeID(coordinates[i].id);
          chromatograms.push_back(chrom);
        }
      }
    }
  }

  void extractAndScore(const std::vector< SwathMap > & swath_maps,
    const TransformationDescription trafo,
    ChromExtractParams cp, String tr_file, String out, 
    Param& feature_finder_param)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    int progress = 0;
    progresslogger.startProgress(0, swath_maps.size(), "Extracting and scoring transitions");

    // Load the OpenMS TransitionExperiment (since ChromExtractor wants an old TransitionExp)
    // OpenSwath::LightTargetedExperiment transition_exp;
#ifdef DEBUG_OPENSWATHWORKFLOW
    std::cout << " Loading TraML file " << std::endl;
#endif
    TargetedExperiment transition_exp;
    {
      // TargetedExperiment *transition_exp__ = new TargetedExperiment();
      // TargetedExperiment &transition_exp_ = *transition_exp__;
      {
        TraMLFile *t = new TraMLFile;
        t->load(tr_file, transition_exp);
        delete t;
      }
      // OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transition_exp);
      // delete transition_exp__;
    }

    FeatureMap<> out_featureFile;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(Size i = 0; i < swath_maps.size(); i++)
    {
      if (swath_maps[i].ms1) {continue;}

      // Step 1: select transitions
      TargetedExperiment transition_exp_used;
      OpenSwathHelper::selectSwathTransitions(transition_exp, transition_exp_used,
          cp.min_upper_edge_dist, swath_maps[i].lower, swath_maps[i].upper);
      if (transition_exp_used.getTransitions().size() == 0) { continue;}

      // Step 2: extract these transitions
      ChromatogramExtractor extractor;
      boost::shared_ptr<MSExperiment<Peak1D> > chrom_tmp(new MSExperiment<Peak1D>);

      std::vector< OpenSwath::ChromatogramPtr > tmp_out; // chrom_tmp
      std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
      prepare_coordinates(tmp_out, coordinates, transition_exp_used, trafo, cp.rt_extraction_window, false);

      extractor.extractChromatograms(swath_maps[i].sptr, tmp_out, coordinates, cp.extraction_window,
          cp.ppm, cp.rt_extraction_window, cp.extraction_function);

      std::vector< OpenMS::MSChromatogram<> > chromatograms;
      for (Size j = 0; j < tmp_out.size(); j++)
      { 
        OpenMS::MSChromatogram<> chrom;
        OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chrom, tmp_out[j]);
        chrom.setNativeID(coordinates[j].id);
        chromatograms.push_back(chrom);
      }
      chrom_tmp->setChromatograms(chromatograms);

      // Convert to LightTargetedExperiment as FeatureFinder requests and convert chroms to ptr
      OpenSwath::LightTargetedExperiment light_transition_exp_used;
      OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_used, light_transition_exp_used);
      OpenSwath::SpectrumAccessPtr chromatogram_ptr = OpenSwath::SpectrumAccessPtr(new OpenMS::SpectrumAccessOpenMS(chrom_tmp));

      // Step 3: score these extracted transitions
      FeatureMap<> featureFile;
      scoreAll_(chromatogram_ptr, swath_maps[i].sptr, light_transition_exp_used, trafo, cp.rt_extraction_window, featureFile, feature_finder_param);

      // write all features and the protein identifications from tmp_featureFile into featureFile
#ifdef _OPENMP
#pragma omp critical (featureFinder)
#endif
      {
        for (FeatureMap<Feature>::iterator feature_it = featureFile.begin();
             feature_it != featureFile.end(); feature_it++)
        {
          out_featureFile.push_back(*feature_it);
        }
        for (std::vector<ProteinIdentification>::iterator protid_it =
               featureFile.getProteinIdentifications().begin();
             protid_it != featureFile.getProteinIdentifications().end();
             protid_it++)
        {
          out_featureFile.getProteinIdentifications().push_back(*protid_it);
        }
        progresslogger.setProgress(progress++);
      }

    }
    addDataProcessing_(out_featureFile, getProcessingInfo_(DataProcessing::QUANTITATION));
    out_featureFile.ensureUniqueId();
    FeatureXMLFile().store(out, out_featureFile);
    progresslogger.endProgress();
  }

  void registerOptionsAndFlags_()
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", StringList::create("mzML"));

    registerInputFile_("tr_irt", "<file>", "", "transition file ('TraML' or 'csv')");
    setValidFormats_("tr_irt", StringList::create("csv,traML"));

    registerInputFile_("tr", "<file>", "", "transition file ('TraML' or 'csv')");
    setValidFormats_("tr", StringList::create("csv,traML"));

    /*
    registerInputFile_("rt_norm", "<file>", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library)", false);
    setValidFormats_("rt_norm", StringList::create("trafoXML"));
    */

    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", StringList::create("featureXML"));


    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0, "Minimal distance to the edge to still consider a precursor, in Thomson", false);
    registerDoubleOption_("extraction_window", "<double>", 0.05, "Extraction window used (in Thomson, to use ppm see -ppm flag)", false);
    registerDoubleOption_("rt_extraction_window", "<double>", -1, "Only extract RT around this value (-1 means extract over the whole range, a value of 500 means to extract around +/- 500 s of the expected elution).", false);
    setMinFloat_("extraction_window", 0.0);

    registerDoubleOption_("min_rsq", "<double>", 0.95, "Minimum r-squared of RT peptides regression", false);
    registerDoubleOption_("min_coverage", "<double>", 0.6, "Minimum relative amount of RT peptides to keep", false);

    registerFlag_("is_swath", "Set this flag if the data is SWATH data");
    registerFlag_("ppm", "extraction_window is in ppm");

    registerStringOption_("readOptions", "<name>", "normal", "Whether to run OpenSWATH directly on the input data, cache data to disk first or to perform a datareduction step first", false);
    setValidStrings_("readOptions", StringList::create("normal,cache,reduce,reduce_iterative"));

    registerStringOption_("extraction_function", "<name>", "tophat", "Function used to extract the signal", false);
    setValidStrings_("extraction_function", StringList::create("tophat,bartlett"));

    registerSubsection_("Scoring", "Scoring parameters section");
  }

  Param getSubsectionDefaults_(const String & name) const
  {
    if (name == "Scoring")
    {
      // set sensible default parameters
      Param feature_finder_param = MRMFeatureFinderScoring().getDefaults();
      feature_finder_param.setValue("rt_extraction_window", 300.0); // 10 minutes (+/- 5 minutes)
      feature_finder_param.setValue("rt_normalization_factor", 100.0);

      feature_finder_param.setValue("TransitionGroupPicker:min_peak_width", 12.0);
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "true");
      feature_finder_param.setValue("TransitionGroupPicker:minimal_quality", -1.5);
      feature_finder_param.remove("TransitionGroupPicker:background_subtraction");
      feature_finder_param.remove("TransitionGroupPicker:stop_after_intensity_ratio");

      // Peak Picker
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:use_gauss", "false");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:sgolay_polynomial_order", 3);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length", 9);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:remove_overlapping_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:method", "corrected");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:gauss_width");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:sn_win_len");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:sn_bin_count");

      // EMG Scoring
      feature_finder_param.setValue("EMGScoring:max_iteration", 10);
      feature_finder_param.setValue("EMGScoring:deltaRelError", 0.1);
      feature_finder_param.remove("EMGScoring:interpolation_step");
      feature_finder_param.remove("EMGScoring:tolerance_stdev_bounding_box");
      feature_finder_param.remove("EMGScoring:deltaAbsError");
         
      // remove these params
      feature_finder_param.remove("stop_report_after_feature");
      feature_finder_param.remove("add_up_spectra");
      feature_finder_param.remove("spacing_for_spectra_resampling");
      feature_finder_param.remove("EMGScoring:statistics:mean");
      feature_finder_param.remove("EMGScoring:statistics:variance");
      return feature_finder_param;
    }
  }

  ExitCodes main_(int, const char **)
  {
    StringList file_list = getStringList_("in");
    String irt_tr_file = getStringOption_("tr_irt");
    String tr_file = getStringOption_("tr");
    String out = getStringOption_("out");
    bool ppm = getFlag_("ppm");
    DoubleReal min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    DoubleReal extraction_window = getDoubleOption_("extraction_window");
    DoubleReal rt_extraction_window = getDoubleOption_("rt_extraction_window");
    String extraction_function = getStringOption_("extraction_function");

    String readoptions = getStringOption_("readOptions");

    String tmp = "/tmp/";
    double irt_extraction_window = -1;

    ChromExtractParams cp;
    cp.min_upper_edge_dist   = min_upper_edge_dist;
    cp.extraction_window     = extraction_window;
    cp.ppm                   = ppm;
    cp.rt_extraction_window  = rt_extraction_window, 
    cp.extraction_function   = extraction_function;

    ChromExtractParams cp_irt = cp;
    cp_irt.rt_extraction_window = irt_extraction_window;

    // Load the SWATH files
    SwathMapLoader sml;
    sml.setLogType(log_type_);
    std::vector< SwathMap > swath_maps = sml.load_files(file_list, tmp, readoptions);

    // Loading iRT file
    TraMLFile traml;
    OpenMS::TargetedExperiment irt_transitions;
    traml.load(irt_tr_file, irt_transitions);
#ifdef DEBUG_OPENSWATHWORKFLOW
    std::cout << "Loaded iRT files" << std::endl;
#endif

    // Extracting the iRT file
    std::vector< OpenMS::MSChromatogram<> > irt_chromatograms;
    simpleExtractChromatograms(swath_maps, irt_transitions, irt_chromatograms, cp_irt);
#ifdef DEBUG_OPENSWATHWORKFLOW
    std::cout << "Extracted iRT files: " << irt_chromatograms.size() <<  std::endl;
#endif

    Param feature_finder_param = getParam_().copy("Scoring:", true);
    DoubleReal min_rsq = getDoubleOption_("min_rsq");
    DoubleReal min_coverage = getDoubleOption_("min_coverage");
    TransformationDescription trafo_rtnorm = RTNormalization(irt_transitions,
            irt_chromatograms, min_rsq, min_coverage, feature_finder_param);

    extractAndScore(swath_maps, trafo_rtnorm, cp, tr_file, out, feature_finder_param);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPOpenSwathWorkflow tool;
  return tool.main(argc, argv);
}

/// @endcond


// Speed analysis
/*
 *_         


-- extraction takes 0.39 seconds or 1.67 %

will score 166 chromatograms
7.15208961654 chromatograms / second [ 23.3 seconds] with all scores
8.4           chromatograms / second [ 19.76 seconds] without the DIA scores
12.7496159754 chromatograms / second [ 13.02 seconds] without the model fit
13.0708661417 chrom / second [12 seconds] without the dia and model fit

504.0 / minute

- copying the chromatograms around costs 0.62 seconds 
- picking the peaks costs [without signal to noise] 2.64 seconds (8%)
- picking the peaks costs [with signal to noise] 7.5 seconds (30%)
- creating the SignalToNoise estimators alone costs 7.82 seconds (34%)
- calculateChromatographicScores (without elution model) 8.14 seconds (alone 0.3 seconds)
- calculateLibraryScores 8.17 seconds
- calculateDIAScores 11.09 seconds (3 seconds or 13%)
- calculateDIAScores 21.10 seconds (10 seconds or 42%) with the EMG scores


print (7.54-0.62)/ 23.3


0.296995708155

0.08669527897
0.11330472103

# 1 Million transitions
print 1000000 / 12.7 *1/(3600.0), "hours"
21.8722659668
print 1000000 / 7.15 *1/(3600.0), "hours"
38.85003885 hours

 will score 517 chromatograms
-- done [took 59.37 s (CPU), 59.69 s (Wall)] -- 
print 517 / 59.3
8.71838111298 chromatograms / second


  File "<string>", line 3
    will score 517 chromatograms
    ^
IndentationError: unexpected indent

*/

// NoiseEstimator Rapid vs standard -> speedtest
/*
 *
 * NoiseEstimator: regular vs rapid -> when just measuring the performance of doing the same estimation 4000 times:
 *  ca 10 seconds for old one, ca 1.2 seconds for Rapid -> 8x increase in performance ... 
 *
 *
 * When extracting real chromatograms and running them through :
 *
 * without any S/N
 *
 * real    1m6.057s
 * user    1m4.592s
 * sys     0m1.116s
 *
 *
 *  OpenMS::SignalToNoiseEstimatorMedian< MSChromatogram<ChromatogramPeak> >().init(chromatogram_old);
 * real    1m46.083s
 * user    1m44.307s
 * sys     0m1.240s
 *
 * delta = 40 seconds
 *
 * SignalToNoiseEstimatorMedianRapid(200).compute(cptr->getTimeArray()->data, cptr->getIntensityArray()->data);
 * real    1m10.709s
 * user    1m9.132s
 * sys     0m1.196s
 *
 * delta = 4 seconds
 *
 *      std::vector<double> mz(chromatogram_old.size()), intensity(chromatogram_old.size());
 *      for (Size p = 0; p < chromatogram_old.size(); ++p)
 *      {
 *        mz[p] = chromatogram_old[p].getMZ();
 *        intensity[p] = chromatogram_old[p].getIntensity();
 *      }
 *      SignalToNoiseEstimatorMedianRapid(200).estimateNoise(mz, intensity);
 *
 * with copying of all data: 
 * real    1m11.566s
 * user    1m9.948s
 * sys     0m1.248s
 * hr@hr-Latitude-E6410
 *
 * delta = 5 seconds
 *

*/

// PeakPicker Maxima vs standard -> speedtest
/* 
  std::vector<PeakPickerMaxima::PeakCandidate> pc;
  PeakPickerMaxima(1.0).findMaxima(cptr->getTimeArray()->data, cptr->getIntensityArray()->data, pc);
  pks += pc.size();
 * Use // -rt_extraction_window -1
 *  peaks 0
 *  -- done [took 17.52 s (CPU), 18.47 s (Wall)] -- 
 *
    MSChromatogram<ChromatogramPeak> other;
    PeakPickerHiRes().pick(chromatogram_old, other);
    pks += other.size();

 peaks 354681
   -- done [took 36.66 s (CPU), 37.70 s (Wall)] -- 
 *
 * result = ca 9.1 seconds time
 *
 *
  std::vector<PeakPickerMaxima::PeakCandidate> pc;
  PeakPickerMaxima(1.0).pick(cptr->getTimeArray()->data, cptr->getIntensityArray()->data, pc);
  pks += pc.size();
  peaks 360842
   -- done [took 22.42 s (CPU), 23.41 s (Wall)] -- 
 *
 * result = ca 4.9 seconds time (of which 1.3 seconds are GSL time)
 *
 *



 * empty to do 14310 pickings: 
 * -- done [took 4.29 s (CPU), 4.38 s (Wall)] --
 *
-- done [took 8.19 s (CPU), 8.33 s (Wall)] -- 
--> time for picking is 3.9 seconds

std::vector<PeakPickerMaxima::PeakCandidate> pc;
PeakPickerMaxima(1.0).pick(cptr->getTimeArray()->data, cptr->getIntensityArray()->data, pc);
pks += pc.size();


--> time for picking is 14.0 seconds
-- done [took 17.99 s (CPU), 18.20 s (Wall)] -- 

MSChromatogram<ChromatogramPeak> other;
PeakPickerHiRes().pick(chromatogram_old, other);
pks += other.size();

*/

// Full analysis using a large library -> speed is up to 30k transitions / minute
/*
 *
 *
 * using strep, on the debug build I have  using the nf values
 *  will score 14310 chromatograms
 *  -- done [took 02:01 m (CPU), 02:02 m (Wall)] --
 * or 7k transitions / minute
 *

 ./bin/OpenSwathWorkflow -in data/split_*_[8].nf.pp.mzML.gz -rt_norm split_napedro_L120420_010_SW-400AQUA.rtnorm.trafoXML  
  -tr /media/data/tmp/oge_plus_shotgun_only_300_fixed_decoy.csv -out_tsv legacy_smallstrep.csv
 real    2m24.533s
 user    2m23.277s
 sys     0m0.672s
 $ wc legacy_smallstrep.csv
    2422  121100 2337568 legacy_smallstrep.csv



    As soon as we go to full swath files, it drops to 1.6k transitions / minute
    but also there are many more features (7 per peakgroup)

    ./bin/OpenSwathWorkflow -in /media/data/tmp/split_napedro/split_napedro_L120420_010_SW-400AQUA__human_2ul_dilution_1_13.mzML.gz -rt_norm split_napedro_L120420_010_SW-400AQUA.rtnorm.trafoXML   -tr /media/data/tmp/oge_plus_shotgun_only_300_fixed_decoy.csv -out_tsv legacy_smallstrep.csv  -readOptions cache

 $ wc legacy_smallstrep.csv
   11620   581000 11504158 legacy_smallstrep.csv


   load file 1.3 GB
   will extract 9718 chromatograms
   will score 9718 chromatograms
   -- done [took 05:47 m (CPU), 05:49 m (Wall)] -- 
   OpenSwathWorkflow took 07:40 m (wall), 07:35 m (CPU), 0.00 s (system), 07:35 m (user).

   real    7m40.504s
   user    7m35.332s

   full strep library, full swathes
   with cache readoptions -> only 48 MB memory while caching, 180 MB while loading TraML, 212 MB while scoring
   -- done [took 05:18 m (CPU), 05:24 m (Wall)] --
   OpenSwathWorkflow took 07:40 m (wall), 07:13 m (CPU), 0.00 s (system), 07:13 m (user).

   real    7m40.630s
   user    7m13.579s
   sys     0m5.548s

   -- optimized
    will score 9718 chromatograms
    -- done [took 01:58 m (CPU), 02:03 m (Wall)] -- 
    OpenSwathWorkflow took 04:13 m (wall), 02:55 m (CPU), 0.00 s (system), 02:55 m (user).

    without DIA scores
    -- done [took 01:40 m (CPU), 01:41 m (Wall)] -- 
    OpenSwathWorkflow took 02:54 m (wall), 02:35 m (CPU), 0.00 s (system), 02:35 m (user).

    without EMG scoring and w/o DIA scoring -> 30k transitions / minute
    will extract 9718 chromatograms
    will score 9718 chromatograms
    -- done [took 21.47 s (CPU), 22.03 s (Wall)] -- 
    OpenSwathWorkflow took 01:43 m (wall), 01:17 m (CPU), 0.00 s (system), 01:17 m (user).

    Datareduction:
      - takes about 03:29 minutes for the regular reduce
      - takes about 2:51 minutes for the iterative reduce

*/

