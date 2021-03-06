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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
        @page TOPP_FeatureFinderIdentification FeatureFinderIdentification

        @brief Detects features in MS1 data based on peptide identifications.

        <CENTER>
        <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderIdentification \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_MapAlignerIdentification</td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
        </tr>
        </table>
        </CENTER>

        This tool uses algorithms for targeted data analysis from the OpenSWATH pipeline.

        <B>The command line parameters of this tool are:</B>
        @verbinclude TOPP_FeatureFinderIdentification.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderIdentification :
  public TOPPBase
{
public:
  TOPPFeatureFinderIdentification() :
    TOPPBase("FeatureFinderIdentification", "Detects features in MS1 data based on peptide identifications.", false)
  {
    rt_term_.setCVIdentifierRef("MS");
    rt_term_.setAccession("MS:1000896");
    rt_term_.setName("normalized retention time");
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file (LC-MS raw data)");
    setValidFormats_("in", StringList::create("mzML"));
    registerInputFile_("id", "<file>", "", 
                       "input file (peptide identifications)");
    setValidFormats_("id", StringList::create("idXML"));
    registerOutputFile_("out", "<file>", "", "output file (features)");
    setValidFormats_("out", StringList::create("featureXML"));
    registerOutputFile_("lib_out","<file>", "", "output file (assay library)",
                        false);
    setValidFormats_("lib_out", StringList::create("traML"));
    registerOutputFile_("chrom_out","<file>", "", "output file (chromatograms)",
                        false);
    setValidFormats_("chrom_out", StringList::create("mzML"));
    registerOutputFile_("trafo_out","<file>", "", 
                        "output file (RT transformation)", false);
    setValidFormats_("trafo_out", StringList::create("trafoXML"));

    addEmptyLine_();
    registerStringOption_("reference_rt", "<choice>", "score", "Method for selecting the reference RT, if there are multiple IDs for a peptide and charge ('score': RT of the best-scoring ID; 'intensity': RT of the ID with the most intense precursor; 'median': median RT of all IDs; 'all': no single reference, use RTs of all IDs, 'adapt': adapt RT windows based on IDs)", false);
    setValidStrings_("reference_rt", 
                     StringList::create("score,intensity,median,all,adapt"));
    registerDoubleOption_("rt_window", "<value>", 180, "RT window size (in sec.) for chromatogram extraction.", false);
    setMinFloat_("rt_window", 0);
    registerDoubleOption_("mz_window", "<value>", 0.03, "m/z window size for chromatogram extraction (in Th or ppm, see 'mz_window_ppm').", false);
    setMinFloat_("mz_window", 0);
    registerFlag_("mz_window_ppm", "Interpret 'mz_window' parameter as a ppm value (default: Th)");
    registerDoubleOption_("isotope_pmin", "<value>", 0.01, "Minimum probability for an isotope to be included in the assay for a peptide.", false);
    setMinFloat_("isotope_pmin", 0);
    setMaxFloat_("isotope_pmin", 1);

    // addEmptyLine_();
    // registerSubsection_("algorithm", "Algorithm parameters section");
  }


  // Param getSubsectionDefaults_(const String& /*section*/) const
  // {
  //   Param combined;
  //   return combined;
  // }


  typedef MSExperiment<Peak1D> PeakMap;

  // mapping: charge -> iterator to peptide
  typedef Map<Int, vector<vector<PeptideIdentification>::iterator> > ChargeMap;
  // mapping: sequence -> charge -> iterator to peptide
  typedef Map<AASequence, ChargeMap> PeptideMap;
  // mapping: assay ID -> RT begin/end

  PeakMap ms_data_; // input LC-MS data
  TargetedExperiment library_; // assay library
  CVTerm rt_term_; // controlled vocabulary term for reference RT
  IsotopeDistribution iso_dist_; // isotope distribution for current peptide
  TransformationDescription trafo_; // RT transformation (to range 0-1)
  String reference_rt_; // value of "reference_rt" parameter
  DoubleReal rt_tolerance_; // half the RT window width

  // comparator for spectrum and RT (for "lower_bound"):
  struct RTLess: public binary_function<MSSpectrum<>, DoubleReal, bool>
  {
    inline bool operator()(const MSSpectrum<>& spec, DoubleReal rt) const
    {
      return spec.getRT() < rt;
    }
    // this overload shouldn't be necessary according to the C++ standard, but
    // may be required in older versions of Microsoft Visual C++:
    inline bool operator()(DoubleReal rt, const MSSpectrum<>& spec) const
    {
      return rt < spec.getRT();
    }
  };

  // comparator for peak and m/z (for "lower_bound"):
  struct MZLess: public binary_function<Peak1D, DoubleReal, bool>
  {
    inline bool operator()(const Peak1D& peak, DoubleReal mz) const
    {
      return peak.getMZ() < mz;
    }
    // this overload shouldn't be necessary according to the C++ standard, but
    // may be required in older versions of Microsoft Visual C++:
    inline bool operator()(DoubleReal mz, const Peak1D& peak) const
    {
      return mz < peak.getMZ();
    }
  };


  // add transitions for a peptide ion to the library:
  void addTransitions_(const String& peptide_id, DoubleReal mz, Int charge)
  {
    // go through different isotopes:
    Size counter = 0;
    for (IsotopeDistribution::ConstIterator iso_it = iso_dist_.begin();
         iso_it != iso_dist_.end(); ++iso_it, ++counter)
    {
      String annotation = "i" + String(counter);
      String transition_name = peptide_id + "_" + annotation;
          
      ReactionMonitoringTransition transition;
      transition.setNativeID(transition_name);
      transition.setPrecursorMZ(mz);
      transition.setProductMZ(mz + Constants::C13C12_MASSDIFF_U * 
                              float(counter) / charge);
      transition.setLibraryIntensity(iso_it->second * 100);
      transition.setMetaValue("annotation", annotation);
      transition.setPeptideRef(peptide_id);
      library_.addTransition(transition);
    }
  }


  // add an assay (peptide and transitions) to the library:
  void addAssay_(TargetedExperiment::Peptide& peptide, const AASequence& seq, 
                 const ChargeMap::value_type& charge_data)
  {
    // get reference RT(s):
    DoubleList rts;
    if (charge_data.second.size() == 1) // only one peptide ID
    {
      rts.push_back(charge_data.second[0]->getMetaValue("RT"));
    }
    else if (reference_rt_ == "score")
    {
      rts.resize(1);
      DoubleReal best_score;
      for (ChargeMap::mapped_type::const_iterator pi_it = 
             charge_data.second.begin(); pi_it != charge_data.second.end();
           ++pi_it)
      {
        const PeptideHit& hit = (*pi_it)->getHits()[0];
        bool higher_better = (*pi_it)->isHigherScoreBetter();
        if ((pi_it == charge_data.second.begin()) || // initial case
            (higher_better && (hit.getScore() > best_score)) ||
            (!higher_better && (hit.getScore() < best_score)))
        {
          best_score = hit.getScore();
          rts[0] = (*pi_it)->getMetaValue("RT");
        }
      }
    }
    else if (reference_rt_ == "intensity")
    {
      rts.resize(1);
      DoubleReal highest_intensity = -1;
      for (ChargeMap::mapped_type::const_iterator pi_it = 
             charge_data.second.begin(); pi_it != charge_data.second.end(); 
           ++pi_it)
      {
        // find precursor:
        DoubleReal ms2_rt = (*pi_it)->getMetaValue("RT");
        DoubleReal prec_mz = (*pi_it)->getMetaValue("MZ");
        // "lower_bound" gives the MS1 spectrum _after_ the MS2 of the ID:
        PeakMap::ConstIterator ms1_it = 
          lower_bound(ms_data_.begin(), ms_data_.end(), ms2_rt, RTLess());
        // the following shouldn't happen, but might if input is combined IDs
        // from different samples - use the current ID only if we have to:
        if ((ms1_it == ms_data_.begin()) && (highest_intensity < 0))
        {
          rts[0] = ms2_rt;
          continue;
        }
        --ms1_it;
        MSSpectrum<>::ConstIterator peak_it = 
          lower_bound(ms1_it->begin(), ms1_it->end(), prec_mz, MZLess());
        if (peak_it == ms1_it->end())
        {
          --peak_it; // assuming the spectrum isn't empty, which it shouldn't be
        }
        // check if previous peak is closer to the precursor in m/z:
        else if ((peak_it != ms1_it->begin()) &&
                 (fabs(peak_it->getMZ() - prec_mz) < 
                  fabs((--peak_it)->getMZ() - prec_mz)))
        {
          ++peak_it; 
        }
        if (peak_it->getIntensity() > highest_intensity)
        {
          highest_intensity = peak_it->getIntensity();
          rts[0] = ms2_rt;
        }
      }    
    }
    else // "median", "all", or "adapt"
    {
      for (ChargeMap::mapped_type::const_iterator pi_it = 
             charge_data.second.begin(); pi_it != charge_data.second.end();
           ++pi_it)
      {
        rts << (*pi_it)->getMetaValue("RT");
      }
      if (reference_rt_ != "all")
      {
        sort(rts.begin(), rts.end());
        bool median_fallback = false; // use "median" to resolve ties in "adapt"
      
        if (reference_rt_ == "adapt")
        {
          // store RT region as pair (length, start point) for easier sorting:
          vector<pair<DoubleReal, DoubleReal> > rt_regions;
          rt_regions.push_back(make_pair(rt_tolerance_ * 2.0, 
                                         rts[0] - rt_tolerance_));
          for (DoubleList::iterator rt_it = ++rts.begin(); rt_it != rts.end();
               ++rt_it)
          {
            pair<DoubleReal, DoubleReal>& rt_region = rt_regions.back();
            if (rt_region.second + rt_region.first >= *rt_it - rt_tolerance_)
            { // regions overlap, join them (same start point, new length):
              rt_region.first = *rt_it + rt_tolerance_ - rt_region.second;
            }
            else // no overlap, start new region:
            {
              rt_regions.push_back(make_pair(rt_tolerance_ * 2.0,
                                             *rt_it - rt_tolerance_));
            }
          }
          sort(rt_regions.begin(), rt_regions.end()); // sort regions by size
          DoubleReal rt_window = rt_regions.back().first;
          peptide.setMetaValue("rt_window", rt_window);
          // are there multiple regions of maximal size?
          Int n = rt_regions.size() - 2; // second to last, counting from zero
          while ((n >= 0) && (rt_regions[n].first == rt_window)) --n;
          if (n == Int(rt_regions.size()) - 2) // only one longest region
          {
            DoubleReal rt_start = rt_regions.back().second;
            peptide.setMetaValue("rt_start", rt_start);
            peptide.setMetaValue("rt_end", rt_start + rt_window);
            rts.resize(1);
            rts[0] = rt_start + rt_window / 2.0;
          }
          else // multiple longest regions -> resolve using median below
          {
            median_fallback = true;
            rts.clear();
            for (++n; n < Int(rt_regions.size()); ++n)
            {
              rts << rt_regions[n].second + rt_regions[n].first / 2.0;
            }
          }
        }
        if ((reference_rt_ == "median") || median_fallback)
        {
          DoubleList::iterator start = rts.begin();
          // even number of IDs? don't take the RT _between_ the middle ones!
          if (rts.size() % 2 == 0) ++start;
          rts[0] = Math::median(start, rts.end(), true);
          rts.resize(1);
        }
      }
    }

    // complete peptide information:
    Int charge = charge_data.first;
    peptide.setChargeState(charge);
    peptide.id = peptide.sequence + "/" + String(charge);
    DoubleReal mz = seq.getMonoWeight(Residue::Full, charge) / charge;

    TargetedExperiment::Peptide copy = peptide;
    for (Size i = 0; i < rts.size(); ++i)
    {
      rt_term_.setValue(trafo_.apply(rts[i]));
      TargetedExperiment::RetentionTime rt;
      rt.addCVTerm(rt_term_);
      peptide.rts.push_back(rt);
      if (rts.size() > 1) peptide.id += ":" + String(i + 1); // use multiple IDs
      library_.addPeptide(peptide);
      addTransitions_(peptide.id, mz, charge);
      peptide = copy; // reset
    }
  }


  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String id = getStringOption_("id");
    String out = getStringOption_("out");
    String lib_out = getStringOption_("lib_out");
    String chrom_out = getStringOption_("chrom_out");
    String trafo_out = getStringOption_("trafo_out");
    reference_rt_ = getStringOption_("reference_rt");
    DoubleReal rt_window = getDoubleOption_("rt_window");
    rt_tolerance_ = rt_window / 2.0;
    DoubleReal mz_window = getDoubleOption_("mz_window");
    bool mz_window_ppm = getFlag_("mz_window_ppm");
    DoubleReal isotope_pmin = getDoubleOption_("isotope_pmin");

    //-------------------------------------------------------------
    // load input
    //-------------------------------------------------------------
    LOG_INFO << "Loading input data..." << endl;
    MzMLFile mzml;
    mzml.setLogType(log_type_);
    mzml.getOptions().addMSLevel(1);
    mzml.load(in, ms_data_);
    if (reference_rt_ == "intensity") ms_data_.sortSpectra(true);

    // RT transformation to range 0-1:
    ms_data_.updateRanges();
    DoubleReal min_rt = ms_data_.getMinRT(), max_rt = ms_data_.getMaxRT();
    TransformationDescription::DataPoints points;
    points.push_back(make_pair(min_rt, 0.0));
    points.push_back(make_pair(max_rt, 1.0));
    trafo_.setDataPoints(points);
    trafo_.fitModel("linear");
    if (!trafo_out.empty())
    {
      TransformationXMLFile().store(trafo_out, trafo_);
    }
    
    vector<PeptideIdentification> peptides;
    vector<ProteinIdentification> proteins;
    IdXMLFile().load(id, proteins, peptides);

    //-------------------------------------------------------------
    // prepare peptide map
    //-------------------------------------------------------------
    LOG_INFO << "Preparing mapping of peptide data..." << endl;
    PeptideMap peptide_map;
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin(); 
         pep_it != peptides.end(); ++pep_it)
    {
      if (pep_it->getHits().empty()) continue;
      pep_it->sort();
      PeptideHit& hit = pep_it->getHits()[0];
      peptide_map[hit.getSequence()][hit.getCharge()].push_back(pep_it);
    }

    //-------------------------------------------------------------
    // create assay library from peptides
    //-------------------------------------------------------------
    LOG_INFO << "Creating assay library..." << endl;
    set<String> protein_accessions;

    for (PeptideMap::iterator pm_it = peptide_map.begin(); 
         pm_it != peptide_map.end(); ++pm_it)
    {
      const AASequence& seq = pm_it->first;
      // LOG_DEBUG << "Peptide: " << seq.toString() << endl;

      // keep track of protein accessions:
      const PeptideHit& hit = pm_it->second.begin()->second[0]->getHits()[0];
      vector<String> current_accessions = hit.getProteinAccessions();
      // missing protein accession would crash OpenSwath algorithms:
      if (current_accessions.empty())
      {
        current_accessions.push_back("not_available");
      }
      protein_accessions.insert(current_accessions.begin(), 
                                current_accessions.end());

      // get isotope distribution for peptide:
      iso_dist_ = seq.getFormula(Residue::Full, 0).getIsotopeDistribution(10);
      iso_dist_.trimLeft(isotope_pmin);
      iso_dist_.trimRight(isotope_pmin);
      iso_dist_.renormalize();

      // create assay for current peptide (fill in charge etc. later):
      TargetedExperiment::Peptide peptide;
      peptide.sequence = seq.toString();
      peptide.protein_refs = current_accessions;

      // go through different charge states:
      for (ChargeMap::iterator cm_it = pm_it->second.begin(); 
           cm_it != pm_it->second.end(); ++cm_it)
      {
        addAssay_(peptide, seq, *cm_it);
      }      
    }
    // add protein references:
    for (set<String>::iterator acc_it = protein_accessions.begin();
         acc_it != protein_accessions.end(); ++acc_it)
    {
      TargetedExperiment::Protein protein;
      protein.id = *acc_it;
      library_.addProtein(protein);
    }

    if (!lib_out.empty())
    {
      TraMLFile().store(lib_out, library_);
    }

    //-------------------------------------------------------------
    // extract chromatograms
    //-------------------------------------------------------------
    LOG_INFO << "Extracting chromatograms..." << endl;
    ChromatogramExtractor extractor;
    PeakMap chrom_data;
    extractor.setLogType(log_type_);
    if (reference_rt_ != "adapt")
    {
      extractor.extractChromatograms(ms_data_, chrom_data, library_, mz_window,
                                     mz_window_ppm, trafo_, rt_window, 
                                     "tophat");
    }
    else
    {
      trafo_.invert(); // needed to reverse RT transformation below
      vector<ChromatogramExtractor::ExtractionCoordinates> coords;
      for (vector<ReactionMonitoringTransition>::const_iterator trans_it =
             library_.getTransitions().begin(); trans_it != 
             library_.getTransitions().end(); ++trans_it)
      {
        const TargetedExperiment::Peptide& peptide = 
          library_.getPeptideByRef(trans_it->getPeptideRef());
        ChromatogramExtractor::ExtractionCoordinates current;
        current.id = trans_it->getNativeID();
        current.mz = trans_it->getProductMZ();
        if (peptide.metaValueExists("rt_start"))
        {
          current.rt_start = peptide.getMetaValue("rt_start");
          current.rt_end = peptide.getMetaValue("rt_end");
        }
        else
        {
          // is this an intuitive way to store/access the RT?!
          DoubleReal rt = peptide.rts[0].getCVTerms()["MS:1000896"][0].
            getValue().toString().toDouble();
          rt = trafo_.apply(rt); // reverse RT transformation
          DoubleReal rt_win = rt_window;
          if (peptide.metaValueExists("rt_window"))
          {
            rt_win = peptide.getMetaValue("rt_window");
          }
          current.rt_start = rt - rt_win / 2.0;
          current.rt_end = rt + rt_win / 2.0;
        }
        coords.push_back(current);
      }
      sort(coords.begin(), coords.end(), ChromatogramExtractor::
           ExtractionCoordinates::SortExtractionCoordinatesByMZ);

      boost::shared_ptr<PeakMap> shared = boost::make_shared<PeakMap>(ms_data_);
      OpenSwath::SpectrumAccessPtr input = 
        SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(shared);
      vector<OpenSwath::ChromatogramPtr> output;
      for (Size i = 0; i < coords.size(); ++i)
      {
        OpenSwath::ChromatogramPtr cp(new OpenSwath::Chromatogram);
        output.push_back(cp);
      }
      vector<MSChromatogram<> > chromatograms;
      extractor.extractChromatograms(input, output, coords, mz_window,
                                     mz_window_ppm, "tophat");
      extractor.return_chromatogram(output, coords, library_, (*shared)[0],
                                    chromatograms, false);
      chrom_data.setChromatograms(chromatograms);
      trafo_.invert(); // revert the "invert" above!
    }
    ms_data_.reset(); // not needed anymore, free up the memory
    if (!chrom_out.empty())
    {
      MzMLFile().store(chrom_out, chrom_data);
    }

    //-------------------------------------------------------------
    // find chromatographic peaks
    //-------------------------------------------------------------
    LOG_INFO << "Finding chromatographic peaks..." << endl;
    FeatureMap<> features;
    MRMFeatureFinderScoring mrm_finder;
    Param params = mrm_finder.getParameters();
    params.setValue("stop_report_after_feature", 1);
    params.setValue("TransitionGroupPicker:PeakPickerMRM:use_gauss", "false");
    params.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
    params.setValue("TransitionGroupPicker:PeakPickerMRM:method", "corrected");
    mrm_finder.setParameters(params);
    mrm_finder.setLogType(log_type_);
    mrm_finder.setStrictFlag(false);
    mrm_finder.pickExperiment(chrom_data, features, library_, trafo_, ms_data_);

    // @TODO add method for resolving overlaps if "reference_rt" is "all"

    //-------------------------------------------------------------
    // fill in missing feature data
    //-------------------------------------------------------------
    LOG_INFO << "Adapting feature data..." << endl;
    for (FeatureMap<>::Iterator feat_it = features.begin(); 
         feat_it != features.end(); ++feat_it)
    {
      feat_it->setMZ(feat_it->getMetaValue("PrecursorMZ"));
      feat_it->setCharge(feat_it->getPeptideIdentifications()[0].getHits()[0].
                         getCharge());
      DoubleReal rt_min = feat_it->getMetaValue("leftWidth");
      DoubleReal rt_max = feat_it->getMetaValue("rightWidth");
      if (feat_it->getConvexHulls().empty()) // add hulls for mass traces
      {
        for (vector<Feature>::iterator sub_it = 
               feat_it->getSubordinates().begin(); sub_it !=
               feat_it->getSubordinates().end(); ++sub_it)
        {
          DoubleReal abs_mz_tol = mz_window / 2.0;
          if (mz_window_ppm) abs_mz_tol = sub_it->getMZ() * abs_mz_tol * 1.0e-6;
          ConvexHull2D hull;
          hull.addPoint(DPosition<2>(rt_min, sub_it->getMZ() - abs_mz_tol));
          hull.addPoint(DPosition<2>(rt_min, sub_it->getMZ() + abs_mz_tol));
          hull.addPoint(DPosition<2>(rt_max, sub_it->getMZ() - abs_mz_tol));
          hull.addPoint(DPosition<2>(rt_max, sub_it->getMZ() + abs_mz_tol));
          feat_it->getConvexHulls().push_back(hull);
        }
      }
    }

    //-------------------------------------------------------------
    // write output
    //-------------------------------------------------------------
    LOG_INFO << "Writing results..." << endl;
    features.ensureUniqueId();
    addDataProcessing_(features, 
                       getProcessingInfo_(DataProcessing::QUANTITATION));
    FeatureXMLFile().store(out, features);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderIdentification tool;
  return tool.main(argc, argv);
}

/// @endcond
