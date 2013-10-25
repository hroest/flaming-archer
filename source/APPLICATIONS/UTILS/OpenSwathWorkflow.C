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

// Interfaces
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

// Consumers
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>

// Files
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVReader.h>
#include <OpenMS/ANALYSIS/OPENSWATH/CachedmzML.h>
#ifdef OPENMS_FORMAT_SWATHFILE_MZXMLSUPPORT
#include "MSDataReader.h"
#endif
#include <OpenMS/FORMAT/SwathFile.h>

// Kernel and implementations
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

// Algorithms
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerIterative.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

// OpenMS base classes
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <assert.h>

using namespace OpenMS;

// Some extra code for the data reducers
namespace OpenMS 
{

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

  void populateMetaData_(String file, boost::shared_ptr<ExperimentalSettings>& exp_meta)
  {
    MSExperiment<Peak1D> tmp;
    MSDataTransformingConsumer c;
    MzMLFile().transform(file, &c, tmp);
    *exp_meta = tmp;
  }

  std::vector< OpenSwath::SwathMap > load_files_reduce(StringList file_list, String /* tmp */, 
    boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions="normal")
  {
    //int progress = 0;
    //startProgress(0, file_list.size(), "Loading data");

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

      if (readoptions == "reduce")
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
        //setProgress(progress++);
      }
    }
    //endProgress();
    return swath_maps;
  }

}

// The workflow class and the TSV writer
namespace OpenMS 
{

  /**
   * @brief Class to write out an OpenSwath TSV output (mProphet input)
   *
   */
  class OPENMS_DLLAPI OpenSwathTSVWriter
  {
    std::ofstream ofs;
    String input_filename_;
    bool doWrite_;

  public:

    OpenSwathTSVWriter(String output_filename, String input_filename = "inputfile") :
      ofs(output_filename.c_str()), 
      input_filename_(input_filename),
      doWrite_(!output_filename.empty())
      {}

    bool isActive() {return doWrite_;}

    void writeHeader()
    {
      ofs << "transition_group_id\trun_id\tfilename\tRT\tid\tSequence\tFullPeptideName\tCharge\tm/z\tIntensity\tProteinName\tdecoy\tassay_rt\tdelta_rt\tleftWidth\tmain_var_xx_swath_prelim_score\tnorm_RT\ttnr_peaks\tpeak_apices_sum\tpotentialOutlier\trightWidth\trt_score\tsn_ratio\ttotal_xic\tvar_bseries_score\tvar_dotprod_score\tvar_intensity_score\tvar_isotope_correlation_score\tvar_isotope_overlap_score\tvar_library_corr\tvar_library_dotprod\tvar_library_manhattan\tvar_library_rmsd\tvar_library_rootmeansquare\tvar_library_sangle\tvar_log_sn_score\tvar_manhatt_score\tvar_massdev_score\tvar_massdev_score_weighted\tvat_norm_rt_score\tvar_xcorr_coelution\tvar_xcorr_coelution_weighted\tvar_xcorr_shape\tvar_xcorr_shape_weighted\tvar_yseries_score\txx_lda_prelim_score\txx_swath_prelim_score\taggr_Peak_Area\taggr_Peak_Apex\taggr_Fragment_Annotation\n";
    }

    String prepareLine(const OpenSwath::LightPeptide & pep,
        const OpenSwath::LightTransition* transition,
        FeatureMap<>& output, String id)
    {
        String result = "";
        String decoy = "0"; // 0 = false
        if (transition->decoy) decoy = "1";
        for (FeatureMap<>::iterator feature_it = output.begin(); feature_it != output.end(); feature_it++)
        {

          char intensity_char[40];
          String aggr_Peak_Area = "";
          String aggr_Peak_Apex = "";
          String aggr_Fragment_Annotation = "";
          for (std::vector<Feature>::iterator sub_it = feature_it->getSubordinates().begin(); sub_it != feature_it->getSubordinates().end(); ++sub_it)
          {
            sprintf(intensity_char, "%f", sub_it->getIntensity());
            aggr_Peak_Area += (String)intensity_char + ";";
            aggr_Peak_Apex +=  "NA;";
            aggr_Fragment_Annotation += (String)sub_it->getMetaValue("native_id") + ";";
          }
          if (!feature_it->getSubordinates().empty())
          {
            aggr_Peak_Area = aggr_Peak_Area.substr(0, aggr_Peak_Area.size() - 1);
            aggr_Peak_Apex = aggr_Peak_Apex.substr(0, aggr_Peak_Apex.size() - 1);
            aggr_Fragment_Annotation = aggr_Fragment_Annotation.substr(0, aggr_Fragment_Annotation.size() - 1);
          }

          String full_peptide_name = "";
          for (int loc = -1; loc <= (int)pep.sequence.size(); loc++)
          {
            if (loc > -1 && loc < (int)pep.sequence.size())
            {
              full_peptide_name += pep.sequence[loc];
            }
            // C-terminal and N-terminal modifications may be at positions -1 or pep.sequence
            for (Size modloc = 0; modloc < pep.modifications.size(); modloc++)
            {
              if (pep.modifications[modloc].location == loc)
              {
                full_peptide_name += "(" + pep.modifications[modloc].unimod_id + ")";
              }
            }
          }

          String line = "";
          line += id + "_run0"
            + "\t" + "0" 
            + "\t" + input_filename_
            + "\t" + (String)feature_it->getRT() 
            + "\t" + "f_" + feature_it->getUniqueId()  // TODO might not be unique!!! 
            + "\t" + pep.sequence
            + "\t" + full_peptide_name
            + "\t" + (String)pep.charge
            + "\t" + (String)transition->precursor_mz
            + "\t" + (String)feature_it->getIntensity() 
            + "\t" + pep.protein_ref
            + "\t" + decoy 
            // Note: missing MetaValues will just produce a DataValue::EMPTY which lead to an empty column
            + "\t" + (String)feature_it->getMetaValue("assay_rt") 
            + "\t" + (String)feature_it->getMetaValue("delta_rt") 
            + "\t" + (String)feature_it->getMetaValue("leftWidth") 
            + "\t" + (String)feature_it->getMetaValue("main_var_xx_swath_prelim_score") 
            + "\t" + (String)feature_it->getMetaValue("norm_RT") 
            + "\t" + (String)feature_it->getMetaValue("nr_peaks") 
            + "\t" + (String)feature_it->getMetaValue("peak_apices_sum") 
            + "\t" + (String)feature_it->getMetaValue("potentialOutlier") 
            + "\t" + (String)feature_it->getMetaValue("rightWidth") 
            + "\t" + (String)feature_it->getMetaValue("rt_score") 
            + "\t" + (String)feature_it->getMetaValue("sn_ratio") 
            + "\t" + (String)feature_it->getMetaValue("total_xic") 
            + "\t" + (String)feature_it->getMetaValue("var_bseries_score") 
            + "\t" + (String)feature_it->getMetaValue("var_dotprod_score") 
            + "\t" + (String)feature_it->getMetaValue("var_intensity_score") 
            + "\t" + (String)feature_it->getMetaValue("var_isotope_correlation_score") 
            + "\t" + (String)feature_it->getMetaValue("var_isotope_overlap_score") 
            + "\t" + (String)feature_it->getMetaValue("var_library_corr") 
            + "\t" + (String)feature_it->getMetaValue("var_library_dotprod") 
            + "\t" + (String)feature_it->getMetaValue("var_library_manhattan") 
            + "\t" + (String)feature_it->getMetaValue("var_library_rmsd") 
            + "\t" + (String)feature_it->getMetaValue("var_library_rootmeansquare") 
            + "\t" + (String)feature_it->getMetaValue("var_library_sangle") 
            + "\t" + (String)feature_it->getMetaValue("var_log_sn_score") 
            + "\t" + (String)feature_it->getMetaValue("var_manhatt_score") 
            + "\t" + (String)feature_it->getMetaValue("var_massdev_score") 
            + "\t" + (String)feature_it->getMetaValue("var_massdev_score_weighted") 
            + "\t" + (String)feature_it->getMetaValue("var_norm_rt_score") 
            + "\t" + (String)feature_it->getMetaValue("var_xcorr_coelution") 
            + "\t" + (String)feature_it->getMetaValue("var_xcorr_coelution_weighted") 
            + "\t" + (String)feature_it->getMetaValue("var_xcorr_shape") 
            + "\t" + (String)feature_it->getMetaValue("var_xcorr_shape_weighted") 
            + "\t" + (String)feature_it->getMetaValue("var_yseries_score") 
            + "\t" + (String)feature_it->getMetaValue("xx_lda_prelim_score") 
            + "\t" + (String)feature_it->getMetaValue("xx_swath_prelim_score") 
            + "\t" + aggr_Peak_Area + "\t" + aggr_Peak_Apex + "\t" + aggr_Fragment_Annotation + "\n";
          result += line;
        } // end of iteration
      return result;
    }

    void writeLines(std::vector<String> to_output)
    {
      for (Size i = 0; i < to_output.size(); i++) { ofs << to_output[i]; }
    }

  };

  /**
   * @brief Class to execute an OpenSwath Workflow
   *
   * performExtraction will perform the OpenSWATH analysis. Optionally, an RT
   * transformation (mapping peptides to normalized space) can be obtained
   * beforehand using performRTNormalization.
   *
   */
  class OPENMS_DLLAPI OpenSwathWorkflow :
    public ProgressLogger
  {
  public:

    /** @brief ChromatogramExtractor parameters
     *
    */
    struct ChromExtractParams 
    {
      double min_upper_edge_dist;
      double extraction_window;
      bool ppm; 
      double rt_extraction_window; 
      String extraction_function;
    };

    /** @brief Compute the alignment against a set of RT-normalization peptides
     *
    */
    TransformationDescription performRTNormalization(const OpenMS::TargetedExperiment & irt_transitions, 
            const std::vector< OpenSwath::SwathMap > & swath_maps, double min_rsq, double min_coverage, 
            const Param& feature_finder_param, const ChromExtractParams cp_irt)
    {
      std::vector< OpenMS::MSChromatogram<> > irt_chromatograms;
      simpleExtractChromatograms(swath_maps, irt_transitions, irt_chromatograms, cp_irt);
      std::cout << "Extracted iRT files: " << irt_chromatograms.size() <<  std::endl;
      // get RT normalization from data
      return RTNormalization(irt_transitions,
              irt_chromatograms, min_rsq, min_coverage, feature_finder_param);
    }

    /** @brief Execute the OpenSWATH workflow on a set of SwathMaps and transitions.
     *
     * Executes the following operations on the given input:
     * 
     * 1. OpenSwathHelper::selectSwathTransitions
     * 2. ChromatogramExtractor prepare, extract
     * 3. scoreAllChromatograms
     * 4. Write out chromatograms and found features
     *
    */
    void performExtraction(const std::vector< OpenSwath::SwathMap > & swath_maps,
      const TransformationDescription trafo,
      ChromExtractParams cp, OpenSwath::LightTargetedExperiment& transition_exp, 
      FeatureMap<>& out_featureFile, String out,
      Param& feature_finder_param, OpenSwathTSVWriter & tsv_writer, 
      MSDataWritingConsumer * chromConsumer, int batchSize)
    {
      tsv_writer.writeHeader();

      TransformationDescription trafo_inverse = trafo;
      trafo_inverse.invert();


      std::cout << "Will analyze " << transition_exp.getTransitions().size() << " transitions in total." << std::endl;
      int progress = 0;
      this->startProgress(0, swath_maps.size(), "Extracting and scoring transitions");
      
      // We set dynamic scheduling such that the maps are worked on in the order
      // in which they were given to the program / acquired. This gives much
      // better load balancing than static allocation.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
      for(Size i = 0; i < swath_maps.size(); i++)
      {
        if (swath_maps[i].ms1) {continue;}

        // Step 1: select transitions
        OpenSwath::LightTargetedExperiment transition_exp_used_all;
        OpenSwathHelper::selectSwathTransitions(transition_exp, transition_exp_used_all,
            cp.min_upper_edge_dist, swath_maps[i].lower, swath_maps[i].upper);
        if (transition_exp_used_all.getTransitions().size() == 0) { continue;}

        int batch_size;
        if (batchSize <= 0 || batchSize >= (int)transition_exp_used_all.getPeptides().size()) 
        {
          batch_size = transition_exp_used_all.getPeptides().size();
        }
        else {batch_size = batchSize;}
#ifdef _OPENMP
#pragma omp critical (featureFinder)
#endif
        { std::cout << "Thread " << 
#ifdef _OPENMP
          omp_get_thread_num() << " " <<
#endif
          "will analyze " << transition_exp_used_all.getPeptides().size() <<  " peptides and "
          << transition_exp_used_all.getTransitions().size() <<  " transitions "
          "from SWATH " << i << " in batches of " << batch_size << std::endl; }
        for (size_t j = 0; j <= (transition_exp_used_all.getPeptides().size() / batch_size) ; j++)
        {
          // Create the new, batch-size transition experiment
          OpenSwath::LightTargetedExperiment transition_exp_used;
          selectPeptidesForBatch_(transition_exp_used_all, transition_exp_used, batch_size, j);

          // Step 2: extract these transitions
          ChromatogramExtractor extractor;
          boost::shared_ptr<MSExperiment<Peak1D> > chrom_exp(new MSExperiment<Peak1D>);

          std::vector< OpenSwath::ChromatogramPtr > chrom_list;
          std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;

          // Step 2.1: prepare the extraction coordinates
          if (cp.rt_extraction_window < 0)
          {
            prepare_coordinates(chrom_list, coordinates, transition_exp_used, cp.rt_extraction_window, false);
          }
          else
          {
            // Use an rt extraction window of 0.0 which will just write the retention time in start / end positions
            prepare_coordinates(chrom_list, coordinates, transition_exp_used, 0.0, false);
            for (std::vector< ChromatogramExtractor::ExtractionCoordinates >::iterator it = coordinates.begin(); it != coordinates.end(); it++)
            {
              it->rt_start = trafo_inverse.apply(it->rt_start) - cp.rt_extraction_window / 2.0;
              it->rt_end = trafo_inverse.apply(it->rt_end) + cp.rt_extraction_window / 2.0;
            }
          }

          // Step 2.2: extract chromatograms
          extractor.extractChromatograms(swath_maps[i].sptr, chrom_list, coordinates, cp.extraction_window,
              cp.ppm, cp.extraction_function);

          // Step 2.3: convert chromatograms back and write to output
          std::vector< OpenMS::MSChromatogram<> > chromatograms;
          extractor.return_chromatogram(chrom_list, coordinates, transition_exp_used,  SpectrumSettings(), chromatograms, false);
          chrom_exp->setChromatograms(chromatograms);
          OpenSwath::SpectrumAccessPtr chromatogram_ptr = OpenSwath::SpectrumAccessPtr(new OpenMS::SpectrumAccessOpenMS(chrom_exp));

          // Step 3: score these extracted transitions
          FeatureMap<> featureFile;
          scoreAllChromatograms(chromatogram_ptr, swath_maps[i].sptr, transition_exp_used, trafo,
              cp.rt_extraction_window, featureFile, feature_finder_param, tsv_writer);

          // Step 4: write all chromatograms and features out into an output object / file 
          // (this needs to be done in a critical section since we only have one
          // output file and one output map).
#ifdef _OPENMP
#pragma omp critical (featureFinder)
#endif
          {
            // write chromatograms to output if so desired
            for (Size j = 0; j < chromatograms.size(); j++)
            {
              chromConsumer->consumeChromatogram(chromatograms[j]); 
            }
            // write features to output if so desired
            if (!out.empty())
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
              this->setProgress(progress++);
            }
          }
        }

      }
      this->endProgress();
    }

  private:

    void selectPeptidesForBatch_(const OpenSwath::LightTargetedExperiment& transition_exp_used_all, 
      OpenSwath::LightTargetedExperiment& transition_exp_used, int batch_size, size_t j)
    {
      // compute batch start/end
      size_t start = j*batch_size;
      size_t end = j*batch_size+batch_size;
      if (end > transition_exp_used_all.peptides.size() ) {end = transition_exp_used_all.peptides.size();}

      // Create the new, batch-size transition experiment
      transition_exp_used.proteins = transition_exp_used_all.proteins;
      transition_exp_used.peptides.insert(transition_exp_used.peptides.end(), 
          transition_exp_used_all.peptides.begin() + start, transition_exp_used_all.peptides.begin() + end);
      copyBatchTransitions_(transition_exp_used.peptides, transition_exp_used_all.transitions, transition_exp_used.transitions);
    }

    void copyBatchTransitions_(const std::vector<OpenSwath::LightPeptide>& used_peptides, 
        const std::vector<OpenSwath::LightTransition>& all_transitions, std::vector<OpenSwath::LightTransition>& output)
    {
      std::set<std::string> selected_peptides;
      for (Size i = 0; i < used_peptides.size(); i++)
      {
        selected_peptides.insert(used_peptides[i].id);
      }

      for (Size i = 0; i < all_transitions.size(); i++)
      {
        if (selected_peptides.find(all_transitions[i].peptide_ref) != selected_peptides.end())
        {
          output.push_back(all_transitions[i]);
        }
      }
    }

    void simpleExtractChromatograms(const std::vector< OpenSwath::SwathMap > & swath_maps,
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
        ChromatogramExtractor extractor;
        // TODO for lrage rt extraction windows!
        extractor.prepare_coordinates(tmp_out, coordinates, transition_exp_used,  cp.rt_extraction_window, false);
        extractor.extractChromatograms(swath_maps[i].sptr, tmp_out, coordinates, cp.extraction_window,
            cp.ppm, cp.extraction_function);

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

    /// @note: feature_finder_param are copied because they are changed here.
    TransformationDescription RTNormalization(TargetedExperiment transition_exp_,
            std::vector< OpenMS::MSChromatogram<> > chromatograms, double min_rsq, double min_coverage, 
            Param feature_finder_param)
    {
      this->startProgress(0, 1, "Retention time normalization");

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
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise", 1.0); // set to 1.0 in all cases

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

      this->endProgress();
      return trafo_out;
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

    /// Helper function to score a set of chromatograms
    void scoreAllChromatograms(OpenSwath::SpectrumAccessPtr input,
           OpenSwath::SpectrumAccessPtr swath_map,
           OpenSwath::LightTargetedExperiment& transition_exp, 
           TransformationDescription trafo, double rt_extraction_window, 
           FeatureMap<Feature>& output, Param& feature_finder_param, OpenSwathTSVWriter & tsv_writer)
    {
      typedef OpenSwath::LightTransition TransitionType;
      // a transition group holds the MSSpectra with the Chromatogram peaks from above
      typedef MRMTransitionGroup<MSSpectrum <ChromatogramPeak>, TransitionType> MRMTransitionGroupType; 
      typedef std::map<String, MRMTransitionGroupType> TransitionGroupMapType;
      // this is the type in which we store the chromatograms for this analysis
      typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram; 

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

      std::vector<String> to_output;
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
        
        if (tsv_writer.isActive()) output.clear();
        featureFinder.scorePeakgroups(transition_group, trafo, swath_map, output);

        if (tsv_writer.isActive())
        {
          const OpenSwath::LightPeptide pep = transition_exp.getPeptides()[ assay_peptide_map[id] ];
          const TransitionType* transition = assay_it->second[0];
          to_output.push_back(tsv_writer.prepareLine(pep, transition, output, id));
        }
      }

      if(tsv_writer.isActive())
      {
#ifdef _OPENMP
#pragma omp critical (scoreAll)
#endif
        {
          tsv_writer.writeLines(to_output);
        }
      }
    }

    void prepare_coordinates(std::vector< OpenSwath::ChromatogramPtr > & output_chromatograms,
      std::vector< ChromatogramExtractorAlgorithm::ExtractionCoordinates > & coordinates,
      OpenSwath::LightTargetedExperiment & transition_exp_used,
      const double rt_extraction_window, const bool ms1) const
    {
      // hash of the peptide reference containing all transitions
      std::map<String, std::vector<OpenSwath::LightTransition*> > peptide_trans_map;
      for (Size i = 0; i < transition_exp_used.getTransitions().size(); i++)
      {
        peptide_trans_map[transition_exp_used.getTransitions()[i].getPeptideRef()].push_back(&transition_exp_used.getTransitions()[i]);
      }
      std::map<String, OpenSwath::LightPeptide* > trans_peptide_map;
      for (Size i = 0; i < transition_exp_used.getPeptides().size(); i++)
      {
        trans_peptide_map[transition_exp_used.getPeptides()[i].id] = &transition_exp_used.getPeptides()[i];
      }

      // Determine iteration size (nr peptides or nr transitions)
      Size itersize;
      if (ms1) {itersize = transition_exp_used.getPeptides().size();}
      else     {itersize = transition_exp_used.getTransitions().size();}

      for (Size i = 0; i < itersize; i++)
      {
        OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
        output_chromatograms.push_back(s);

        ChromatogramExtractor::ExtractionCoordinates coord;
        OpenSwath::LightPeptide pep; // TargetedExperiment::Peptide pep;
        OpenSwath::LightTransition transition;

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
          pep = (*trans_peptide_map[transition.getPeptideRef()]);
          coord.mz = transition.getProductMZ();
          coord.id = transition.getNativeID();
        }

        double rt = pep.rt;
        coord.rt_start = rt - rt_extraction_window / 2.0;
        coord.rt_end = rt + rt_extraction_window / 2.0;
        coordinates.push_back(coord);
      }

      // sort result
      std::sort(coordinates.begin(), coordinates.end(), ChromatogramExtractor::ExtractionCoordinates::SortExtractionCoordinatesByMZ);
    }

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

  };

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

  void registerOptionsAndFlags_()
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", StringList::create("mzML,mzXML"));

    registerInputFile_("tr", "<file>", "", "transition file ('TraML' or 'csv')");
    setValidFormats_("tr", StringList::create("csv,traML"));

    // one of the following two needs to be set
    registerInputFile_("tr_irt", "<file>", "", "transition file ('TraML' or 'csv')", false);
    setValidFormats_("tr_irt", StringList::create("csv,traML"));

    registerInputFile_("rt_norm", "<file>", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library). If set, tr_irt may be omitted.", false);
    setValidFormats_("rt_norm", StringList::create("trafoXML"));

    registerStringOption_("swath_windows_file", "<file>", "", "Optional, tab separated file containing the SWATH windows: lower_offset upper_offset \\newline 400 425 \\newline ... ", false);

    // one of the following two needs to be set
    registerOutputFile_("out_features", "<file>", "", "output file", false);
    setValidFormats_("out_features", StringList::create("featureXML"));

    registerStringOption_("out_tsv", "<file>", "", "tsv output file (mProphet compatible)", false);

    registerOutputFile_("out_chrom", "<file>", "", ".chrom.mzML output (all chromatograms)", false);
    setValidFormats_("out_chrom", StringList::create("mzML"));

    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0, "Minimal distance to the edge to still consider a precursor, in Thomson", false);
    registerDoubleOption_("extraction_window", "<double>", 0.05, "Extraction window used (in Thomson, to use ppm see -ppm flag)", false);
    registerDoubleOption_("rt_extraction_window", "<double>", 300.0, "Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution).", false);
    setMinFloat_("extraction_window", 0.0);

    registerDoubleOption_("min_rsq", "<double>", 0.95, "Minimum r-squared of RT peptides regression", false);
    registerDoubleOption_("min_coverage", "<double>", 0.6, "Minimum relative amount of RT peptides to keep", false);

    registerFlag_("ppm", "extraction_window is in ppm");
    registerFlag_("split_file_input", "The input files each contain one single SWATH (alternatively: all SWATH are in separate files)");

    registerStringOption_("readOptions", "<name>", "normal", "Whether to run OpenSWATH directly on the input data, cache data to disk first or to perform a datareduction step first", false);
    setValidStrings_("readOptions", StringList::create("normal,cache,reduce,reduce_iterative"));

    registerStringOption_("tempDirectory", "<tmp>", "/tmp/", "Temporary directory to store cached files for example", false);

    registerStringOption_("extraction_function", "<name>", "tophat", "Function used to extract the signal", false);
    setValidStrings_("extraction_function", StringList::create("tophat,bartlett"));

    registerIntOption_("batchSize", "<number>", 0, "The batch size of chromatograms to process (0 means to only have one batch, sensible values are around 500-1000)", false);
    setMinInt_("batchSize", 0);

    registerSubsection_("Scoring", "Scoring parameters section");
  }

  Param getSubsectionDefaults_(const String & name) const
  {
    if (name == "Scoring")
    {
      // set sensible default parameters
      Param feature_finder_param = MRMFeatureFinderScoring().getDefaults();
      feature_finder_param.remove("rt_extraction_window");
      feature_finder_param.setValue("rt_normalization_factor", 100.0); // for iRT peptides between 0 and 100 (more or less)

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
      // TODO it seems that the legacy method produces slightly larger peaks, e.g. it will not cut off peaks too early
      // however the same can be achieved by using a relatively low SN cutoff in the -Scoring:TransitionGroupPicker:PeakPickerMRM:signal_to_noise 0.5
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:method", "corrected");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise", 0.1);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:gauss_width", 30);
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
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Unknown subsection", name);
    }
  }

  void readSwathWindows(String filename, std::vector<double> & swath_prec_lower_,
    std::vector<double> & swath_prec_upper_ )
  {
    std::ifstream data(filename.c_str());
    std::string   line;
    std::string   tmp;
    std::getline(data, line); //skip header
    double lower, upper;
    while (std::getline(data, line))
    {
      std::stringstream lineStream(line);

      lineStream >> lower;
      lineStream >> upper;

      swath_prec_lower_.push_back(lower);
      swath_prec_upper_.push_back(upper);
    }
    assert(swath_prec_lower_.size() == swath_prec_upper_.size());
  }

  void annotateSwathMapsFromFile(String filename,
    std::vector< OpenSwath::SwathMap >& swath_maps)
  {
    std::vector<double> swath_prec_lower_, swath_prec_upper_;
    readSwathWindows(filename, swath_prec_lower_, swath_prec_upper_);
    assert(swath_prec_lower_.size() == swath_maps.size());
    for (Size i = 0; i < swath_maps.size(); i++)
    {
      swath_maps[i].lower = swath_prec_lower_[i];
      swath_maps[i].upper = swath_prec_upper_[i];
    }
  }
      
  void loadSwathFiles(StringList& file_list, bool split_file, String tmp, String readoptions,
    boost::shared_ptr<ExperimentalSettings > & exp_meta,
    std::vector< OpenSwath::SwathMap > & swath_maps)
  {
    SwathFile swath_file;
    swath_file.setLogType(log_type_);

    if (split_file || file_list.size() > 1)
    {
      // TODO cannot use data reduction here any more ...
      swath_maps = swath_file.loadSplit(file_list, tmp, exp_meta, readoptions);
    }
    else 
    {
      FileTypes::Type in_file_type = FileTypes::nameToType(file_list[0]);
      if (in_file_type == FileTypes::MZML || file_list[0].suffix(4).toLower() == "mzml"  
        || file_list[0].suffix(7).toLower() == "mzml.gz"  )
      {
        swath_maps = swath_file.loadMzML(file_list[0], tmp, exp_meta, readoptions);
      }
      else if (in_file_type == FileTypes::MZXML || file_list[0].suffix(5).toLower() == "mzxml"  
        || file_list[0].suffix(8).toLower() == "mzxml.gz"  )
      {
        swath_maps = swath_file.loadMzXML(file_list[0], tmp, exp_meta, readoptions);
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "Input file needs to have ending mzML or mzXML");
      }
    }
  }

  TransformationDescription loadTrafoFile(String trafo_in, String irt_tr_file,
    const std::vector< OpenSwath::SwathMap > & swath_maps, double min_rsq, double min_coverage, 
    const Param& feature_finder_param, const OpenSwathWorkflow::ChromExtractParams& cp_irt)
  {
    TransformationDescription trafo_rtnorm;
    if (trafo_in.size() > 0) 
    {
      // get read RT normalization file
      TransformationXMLFile trafoxml;
      trafoxml.load(trafo_in, trafo_rtnorm);
      Param model_params = getParam_().copy("model:", true);
      model_params.setValue("symmetric_regression", "false");
      String model_type = "linear";
      trafo_rtnorm.fitModel(model_type, model_params);
    }
    else
    {
      OpenSwathWorkflow wf;
      wf.setLogType(log_type_);
      // Loading iRT file
      std::cout << "Will load iRT transitions and try to find iRT peptides" << std::endl;
      TraMLFile traml;
      OpenMS::TargetedExperiment irt_transitions;
      traml.load(irt_tr_file, irt_transitions);
      trafo_rtnorm = wf.performRTNormalization(irt_transitions, swath_maps, min_rsq, min_coverage, 
          feature_finder_param, cp_irt);
    }
    return trafo_rtnorm;
  }

  ExitCodes main_(int, const char **)
  {
    ///////////////////////////////////
    // Prepare Parameters
    ///////////////////////////////////
    StringList file_list = getStringList_("in");
    String tr_file = getStringOption_("tr");

    String out = getStringOption_("out_features");
    String out_tsv = getStringOption_("out_tsv");

    String irt_tr_file = getStringOption_("tr_irt");
    String trafo_in = getStringOption_("rt_norm");

    String out_chrom = getStringOption_("out_chrom");
    bool ppm = getFlag_("ppm");
    bool split_file = getFlag_("split_file_input");
    DoubleReal min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    DoubleReal extraction_window = getDoubleOption_("extraction_window");
    DoubleReal rt_extraction_window = getDoubleOption_("rt_extraction_window");
    String extraction_function = getStringOption_("extraction_function");
    String swath_windows_file = getStringOption_("swath_windows_file");
    int batchSize = (int)getIntOption_("batchSize");

    DoubleReal min_rsq = getDoubleOption_("min_rsq");
    DoubleReal min_coverage = getDoubleOption_("min_coverage");

    String readoptions = getStringOption_("readOptions");
    String tmp = getStringOption_("tempDirectory");

    if (trafo_in.empty() && irt_tr_file.empty()) 
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
              "Either rt_norm or tr_irt needs to be set");
    if (out.empty() && out_tsv.empty()) 
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
              "Either out_features or out_tsv needs to be set");

    OpenSwathWorkflow::ChromExtractParams cp;
    cp.min_upper_edge_dist   = min_upper_edge_dist;
    cp.extraction_window     = extraction_window;
    cp.ppm                   = ppm;
    cp.rt_extraction_window  = rt_extraction_window, 
    cp.extraction_function   = extraction_function;

    OpenSwathWorkflow::ChromExtractParams cp_irt = cp;
    cp_irt.rt_extraction_window = -1; // extract the whole RT range

    Param feature_finder_param = getParam_().copy("Scoring:", true);

    ///////////////////////////////////
    // Load the SWATH files
    ///////////////////////////////////
    boost::shared_ptr<ExperimentalSettings > exp_meta(new ExperimentalSettings);
    std::vector< OpenSwath::SwathMap > swath_maps;
    loadSwathFiles(file_list, split_file, tmp, readoptions, exp_meta, swath_maps);

    if (!swath_windows_file.empty())
      annotateSwathMapsFromFile(swath_windows_file, swath_maps);

    ///////////////////////////////////
    // Get the transformation information (using iRT peptides)
    ///////////////////////////////////
    
    TransformationDescription trafo_rtnorm = loadTrafoFile(trafo_in, irt_tr_file,
        swath_maps, min_rsq, min_coverage, feature_finder_param, cp_irt);

    ///////////////////////////////////
    // Load the transitions
    ///////////////////////////////////
    OpenSwath::LightTargetedExperiment transition_exp;
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    progresslogger.startProgress(0, swath_maps.size(), "Load TraML file");
    FileTypes::Type tr_file_type = FileTypes::nameToType(tr_file);
    if (tr_file_type == FileTypes::TRAML || tr_file.suffix(5).toLower() == "traml"  )
    {
      TargetedExperiment targeted_exp;
      TraMLFile().load(tr_file, targeted_exp);
      OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, transition_exp);
    }
    else
    {
      TransitionTSVReader().convertTSVToTargetedExperiment(tr_file.c_str(), transition_exp);
    }
    progresslogger.endProgress();

    ///////////////////////////////////
    // Set up chrom.mzML output
    ///////////////////////////////////
    MSDataWritingConsumer * chromConsumer;
    if (!out_chrom.empty())
    {
      chromConsumer = new PlainMSDataWritingConsumer(out_chrom);
      chromConsumer->setExpectedSize(0, transition_exp.transitions.size());
      chromConsumer->setExperimentalSettings(*exp_meta);
      chromConsumer->addDataProcessing(getProcessingInfo_(DataProcessing::SMOOTHING));
    }
    else
    {
      chromConsumer = new NoopMSDataWritingConsumer(out_chrom);
    }

    ///////////////////////////////////
    // Extract and score
    ///////////////////////////////////
    FeatureMap<> out_featureFile;

    OpenSwathTSVWriter tsvwriter(out_tsv, "/tmp/out.featureXML");
    OpenSwathWorkflow wf;
    wf.setLogType(log_type_);
    wf.performExtraction(swath_maps, trafo_rtnorm, cp, transition_exp, out_featureFile, out,
        feature_finder_param, tsvwriter, chromConsumer, batchSize);
    if (!out.empty())
    {
      addDataProcessing_(out_featureFile, getProcessingInfo_(DataProcessing::QUANTITATION));
      out_featureFile.ensureUniqueId();
      FeatureXMLFile().store(out, out_featureFile);
    }

    delete chromConsumer;

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

