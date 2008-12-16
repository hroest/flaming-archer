// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_ANALYSIS_ID_PRECURSORIONSELECTION_H
#define OPENMS_ANALYSIS_ID_PRECURSORIONSELECTION_H

#include <OpenMS/ANALYSIS/ID/PrecursorIonSelectionPreprocessing.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>

//#include <cmath>
#include <set>
#include <fstream>
namespace OpenMS
{
  
  class PrecursorIonSelection : public DefaultParamHandler
  {
  public:

    enum Type
			{
				IPS,
				SPS,
				UPSHIFT,
				DOWNSHIFT,
				DEX
      };
      
		PrecursorIonSelection();
		PrecursorIonSelection(const PrecursorIonSelection& source);
		~PrecursorIonSelection();
		
		const DoubleReal& getMaxScore() const;
		void setMaxScore(const DoubleReal& max_score);

		/// returns a const reference to the PeptideIdentification vector
		inline const std::vector<PeptideIdentification>& getPeptideIdentifications() const
		{
			return peptide_ids_;
		};
		
		/// returns a mutable reference to the PeptideIdentification vector
		inline std::vector<PeptideIdentification>& getPeptideIdentifications()
		{
			return peptide_ids_;
		};
		
		/// sets the PeptideIdentification vector
		inline void setPeptideIdentifications( const std::vector<PeptideIdentification>& peptide_ids )
		{
			peptide_ids_ = peptide_ids;
		};
		
		/// Compare by score
		struct TotalScoreMore
			: std::binary_function < Feature, Feature, bool >
		{
			inline bool operator () ( Feature const & left, Feature const & right ) const
			{
				return ( (DoubleReal)left.getMetaValue("msms_score") > (DoubleReal)right.getMetaValue("msms_score") );
			}
		};

  	/** @brief Sort features by total score. */
		void sortByTotalScore(FeatureMap<>& features) 
		{ 
			FeatureMap<>::Iterator beg = features.begin();
			FeatureMap<>::Iterator end  = features.end();
			std::sort(beg,end,TotalScoreMore()); 
		}

		/**
		 *	@brief Returns features with highest score for MS/MS
		 *	
		 *	@param features FeatureMap with all possible precursors
		 *	@param next_features FeatureMap with next precursors
		 *	@param number Number of features to be reported
		 *
		 */
    void getNextPrecursors(FeatureMap<>& features,FeatureMap<>& next_features,UInt number);
    
// 		/**
// 		 *	@brief Change scoring of features using peptide identifications only from spectra of the last
// 		 *  iteration
// 		 *	
// 		 *	@param features FeatureMap with all possible precursors
// 		 *	@param new_pep_ids Peptide identifications
// 		 *	@param preprocessed_db Information from preprocessed database
// 		 *
// 		 */
//     void rescoreIncremental(FeatureMap<>& features,std::vector<PeptideIdentification>& new_pep_ids,
// 														std::vector<ProteinIdentification>& prot_ids,
// 														PrecursorIonSelectionPreprocessing& preprocessed_db);

		
		/**
		 *	@brief Change scoring of features using peptide identifications from all spectra.
		 *	
		 *	@param features FeatureMap with all possible precursors
		 *	@param new_pep_ids Peptide identifications
		 *	@param preprocessed_db Information from preprocessed database
		 *
		 */
    void rescore(FeatureMap<>& features,std::vector<PeptideIdentification>& new_pep_ids,
								 std::vector<ProteinIdentification>& prot_ids,
								 PrecursorIonSelectionPreprocessing& preprocessed_db, bool check_meta_values=true);


		/**
		 *	@brief Simulate the iterative precursor ion selection. 
		 *	
		 *	@param features FeatureMap with all possible precursors
		 *	@param pep_ids Peptide identifications
		 *	@param preprocessed_db Information from preprocessed database
		 *  @param step_size Number of MS/MS spectra considered per iteration
		 *
		 */
    void simulateRun(FeatureMap<>& features,std::vector<PeptideIdentification>& pep_ids,
										 std::vector<ProteinIdentification>& prot_ids,
										 PrecursorIonSelectionPreprocessing& preprocessed_db,
										 UInt step_size, String path);
    

		
  private:

		void shiftDown_(FeatureMap<>& features,PrecursorIonSelectionPreprocessing& preprocessed_db,String protein_acc);

		void shiftUp_(FeatureMap<>& features,PrecursorIonSelectionPreprocessing& preprocessed_db,String protein_acc);

		void inferProteinIds_(std::vector<PeptideIdentification>& new_pep_ids,std::vector<ProteinIdentification>& prot_ids);
		
		/// update members method from DefaultParamHandler to update the members 
		void updateMembers_();

		void rescore_(FeatureMap<>& features,std::vector<PeptideIdentification>& new_pep_ids,
									PrecursorIonSelectionPreprocessing& preprocessed_db);

		/**
		 *	@brief Adds user params, required for the use of IPS, to a feature map using default values.
		 *	
		 *	@param features FeatureMap with all possible precursors
		 *
		 */
		void checkForRequiredUserParams_(FeatureMap<>& features);

		/**
		 *	@brief Groups protein identifications that cannot be distinguished by their peptide identifications.
		 *	
		 *	@param prot_ids All protein identifications.
		 *
		 */
		UInt filterProtIds_(std::vector<ProteinIdentification>& prot_ids);
		
		std::vector<PeptideIdentification> filterPeptideIds_(std::vector<PeptideIdentification>& pep_ids);
		
		/// minimal number of peptides identified for a protein to be declared identified
		UInt min_pep_ids_;
    /// maximal score in the FeatureMap
    DoubleReal max_score_;
    /// precursor ion selection strategy
    Type type_;
		/// peptide identications
		std::vector<PeptideIdentification> peptide_ids_;
		/// stores the peptide sequences for all protein identifications
		std::map<String,std::set<String> > prot_id_counter_;
		/// precursor ion error tolerance 
		DoubleReal mz_tolerance_;
		/// precursor ion error tolerance unit (ppm or Da)
		String mz_tolerance_unit_;
		
  };

}


#endif // #ifndef OPENMS_ANALYSIS_ID_PRECURSORIONSELECTION_H
