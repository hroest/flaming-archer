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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_SEQUESTINFILE_H
#define OPENMS_FORMAT_SEQUESTINFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>


namespace OpenMS
{
	/**
		@brief Sequest input file adapter.
		
		Creates a sequest.params file for Sequest search from a peak list.
  	  	
  	@ingroup FileIO
	*/
  class SequestInfile
  {
		public:
			/// default constructor
			SequestInfile();

			/// copy constructor
			SequestInfile(const SequestInfile& sequest_infile);

			/// destructor
			virtual ~SequestInfile();

			/// assignment operator
			SequestInfile& operator=(const SequestInfile& sequest_infile);

			/// stores the experiment data in a Sequest input file that can be used as input for Sequest shell execution
			void store(const String& filename) throw (Exception::UnableToCreateFile);

			/// returns the enzyme list
			const std::map< String, std::vector< String > >& getEnzymeInfo() const;
			/// returns the enzyme list as a string
			const String getEnzymeInfoAsString() const;

			/// returns the used database
			const String& getDatabase() const;
			/// sets the used database
			void setDatabase(const String& database);

			/// returns whether neutral losses are considered for the a-, b- and y-ions
			const String& getNeutralLossesForIons() const;
			/// sets whether neutral losses are considered for the a-, b- and y-ions
			void setNeutralLossesForIons(const String& neutral_losses_for_ions);

			/// returns the weights for the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series
			const String& getIonSeriesWeights() const;
			/// sets the weights for the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series
			void setIonSeriesWeights(const String& ion_series_weights);

			/// returns the partial sequences (space delimited) that have to occur in the theortical spectra
			const String& getPartialSequence() const;
			/// sets the partial sequences (space delimited) that have to occur in the theortical spectra
			void setPartialSequence(const String& partial_sequence);

			/// returns the sequences (space delimited) that have to occur, or be absent (preceeded by a tilde) in the header of a protein to be considered
			const String& getSequenceHeaderFilter() const;
			/// sets the sequences (space delimited) that have to occur, or be absent (preceeded by a tilde) in the header of a protein to be considered
			void setSequenceHeaderFilter(const String& sequence_header_filter);

			/// returns the protein mass filter (either min and max mass, or mass and tolerance value in percent)
			const String& getProteinMassFilter() const;
			/// sets the protein mass filter (either min and max mass, or mass and tolerance value in percent)
			void setProteinMassFilter(const String& protein_mass_filter);


			/// returns the peak mass tolerance
			Real getPeakMassTolerance() const;
			/// sets the peak mass tolerance
			void setPeakMassTolerance(Real peak_mass_tolerance);

			/// returns the precursor mass tolerance
			Real getPrecursorMassTolerance() const;
			/// sets the precursor mass tolerance
			void setPrecursorMassTolerance(Real precursor_mass_tolerance);

			/// returns the match peak tolerance
			Real getMatchPeakTolerance() const;
			/// sets the match peak tolerance
			void setMatchPeakTolerance(Real match_peak_tolerance);

			/// returns the the cutoff of the ratio matching theoretical peaks/theoretical peaks
			Real getIonCutoffPercentage() const;
			/// sets the ion cutoff of the ratio matching theoretical peaks/theoretical peaks
			void setIonCutoffPercentage(Real cutoff_percentage);

			/// returns the peptide mass unit
			Int getPeptideMassUnit() const;
			/// sets the peptide mass unit
			void setPeptideMassUnit(Int peptide_mass_unit);

			/// return the number of peptides to be displayed
			Int getOutputLines() const;
			/// sets the number of peptides to be displayed
			void setOutputLines(Int output_lines);

			/// returns the enzyme used for cleavage (by means of the number from a list of enzymes)
			Int getEnzymeNumber() const;
			/// returns the enzyme used for cleavage
			String getEnzyme() const;
			/// sets the enzyme used for cleavage (by means of the number from a list of enzymes)
			Int setEnzyme(String enzyme);

			/// returns the maximum number of amino acids containing the same modification in a peptide
			Int getMaxAAPerModPerPeptide() const;
			/// sets the maximum number of amino acids containing the same modification in a peptide
			void setMaxAAPerModPerPeptide(Int max_aa_per_mod_per_peptide);

			/// returns the maximum number of modifications that are allowed in a peptide
			Int getMaxModsPerPeptide() const;
			/// set the maximum number of modifications that are allowed in a peptide
			void setMaxModsPerPeptide(Int max_mods_per_peptide);

			/// returns the nucleotide reading frame
			Int getNucleotideReadingFrame() const;
			/// sets the nucleotide reading frame:
			///		0 	The FASTA file contains amino acid codes. No translation is needed. This is the best and fastest case.
			///		1 	The DNA sequence is scanned left to right (forward direction). The amino acid code starts with the first DNA code.
			///		2 	The DNA sequence is scanned left to right (forward direction). The amino acid code starts with the second DNA code.
			///		3 	The DNA sequence is scanned left to right (forward direction). The amino acid code starts with the third DNA code.
			///		4 	The DNA sequence is scanned right to left (backward direction for the complementary strand). The amino acid code starts with the first DNA code.
			///		5 	The DNA sequence is scanned right to left (backward direction for the complementary strand). The amino acid code starts with the second DNA code.
			///		6 	The DNA sequence is scanned right to left (backward direction for the complementary strand). The amino acid code starts with the third DNA code.
			///		7 	Use each of the DNA translations of the codes 1, 2, 3.
			///		8 	Use each of the DNA translations of the codes 4, 5, 6.
			///		9 	Use each of the DNA translations of the codes 1, 2, 3, 4, 5, 6.
			void setNucleotideReadingFrame(Int nucleotide_reading_frame);

			/// returns the maximum number of internal cleavage sites
			Int getMaxInternalCleavageSites() const;
			/// sets the maximum number of internal cleavage sites
			void setMaxInternalCleavageSites(Int max_internal_cleavage_sites);

			/// returns the number of top abundant peaks to match with theoretical ones
			Int getMatchPeakCount() const;
			/// sets the number of top abundant peaks to with theoretical ones
			void setMatchPeakCount(Int match_peak_count);

			/// returns the number of top abundant peaks that are allowed not to match with a theoretical peak
			Int getMatchPeakAllowedError() const;
			/// sets the number of top abundant peaks that are allowed not to match with a theoretical peak
			void setMatchPeakAllowedError(Int match_peak_allowed_error);


			/// returns whether fragment ions shall be displayed
			bool getShowFragmentIons() const;
			/// sets whether fragment ions shall be displayed
			void setShowFragmentIons(bool show_fragments);

			/// returns whether all proteins containing a found peptide should be displayed
			bool getPrintDuplicateReferences() const;
			/// sets whether all proteins containing a found peptide should be displayed
			void setPrintDuplicateReferences(bool print_duplicate_references);

			/// return whether peaks near (15 amu) the precursor peak are removed
			bool getRemovePrecursorNearPeaks() const;
			/// sets whether peaks near (15 amu) the precursor peak are removed
			void setRemovePrecursorNearPeaks(bool remove_precursor_near_peaks);

			/// return the mass type of the parent (0 - monoisotopic, 1 - average mass)
			bool getMassTypeParent() const;
			/// sets the mass type of the parent (0 - monoisotopic, 1 - average mass)
			void setMassTypeParent(bool mass_type_parent);

			/// return the mass type of the fragments (0 - monoisotopic, 1 - average mass)
			bool getMassTypeFragment() const;
			/// sets the mass type of the fragments (0 - monoisotopic, 1 - average mass)
			void setMassTypeFragment(bool mass_type_fragment);

			/// returns whether normalized xcorr values are displayed
			bool getNormalizeXcorr() const;
			/// sets whether normalized xcorr values are displayed
			void setNormalizeXcorr(bool normalize_xcorr);

			/// returns whether residues are in upper case
			bool getResiduesInUpperCase() const;
			/// sets whether residues are in upper case
			void setResiduesInUpperCase(bool residues_in_upper_case);
			
			/// adds an enzyme to the list and sets is as used
			/// the vector constists of four strings:
			/// name, cut direction: 0 (N to C) / 1, cuts after (list of aa), doesn't cut before (list of aa)
			void addEnzymeInfo(std::vector< String >& enzyme_info);
			
			/// return the modifications (the modification names map to the affected residues, the mass change and the type)
			const std::map< String, std::vector< String > >& getModifications() const;
			
			/// retrieves the name, mass change, affected residues, type and position for all modifications from a string
			void handlePTMs(const String& modification_line, const String& modifications_filename, const bool monoisotopic) throw (Exception::FileNotReadable, Exception::FileNotFound, Exception::ParseError);

		protected:
			/// returns some standard enzymes (used to initialize the enzyme list)
			void setStandardEnzymeInfo_();

			/// the amino acids in one-letter-code
			static const String aas_single_letter_;

			std::map< String, std::vector< String > > enzyme_info_; ///< an endline-delimited list of enzymes; each with cutting direction 0 (N to C) /1; cuts after (list of aa); doesn't cut before (list of aa); the attributes are tab-delimited
			String database_; ///< database used
			String snd_database_; ///< second database used
			String neutral_losses_for_ions_; ///< whether neutral losses are considered for the a-; b- and y-ions (e.g. 011 for b- and y-ions)
			String ion_series_weights_;///< weights for the a-; b-; c-; d-; v-; w-; x-; y- and z-ion series; space delimited
			String partial_sequence_; ///< space-delimited list of sequence parts that have to occur in the theortical spectra
			String sequence_header_filter_;///< space-delimited list of sequences that have to occur or be absend (preceeded by a tilde) in a protein header; to be considered
			String protein_mass_filter_;
			
			Real precursor_mass_tolerance_;///< tolerance for matching a theoretical to an experimental peptide
			Real peak_mass_tolerance_;///< tolerance for matching a theoretical to an experimental peak
			Real match_peak_tolerance_;///< minimum distance between two experimental peaks
			Real ion_cutoff_percentage_;///< cutoff of the ratio matching theoretical peaks/theoretical peaks
			
			UInt peptide_mass_unit_;///< peptide mass unit (0 = amu; 1 = mmu; 2 = ppm)
			UInt output_lines_;///< number of peptides to be displayed
			UInt enzyme_number_;///< number of the enzyme used for cleavage
			UInt max_AA_per_mod_per_peptide_;///< maximum number of amino acids containing the same modification in a peptide
			UInt max_mods_per_peptide_;///< maximum number of modifications per peptide
			UInt nucleotide_reading_frame_;///< nucleotide reading frame:
					/// 0 	The FASTA file contains amino acid codes. No translation is needed. This is the best and fastest case.
					/// 1 	The DNA sequence is scanned left to right (forward direction). The amino acid code starts with the first DNA code.
					/// 2 	The DNA sequence is scanned left to right (forward direction). The amino acid code starts with the second DNA code.
					/// 3 	The DNA sequence is scanned left to right (forward direction). The amino acid code starts with the third DNA code.
					/// 4 	The DNA sequence is scanned right to left (backward direction for the complementary strand). The amino acid code starts with the first DNA code.
					/// 5 	The DNA sequence is scanned right to left (backward direction for the complementary strand). The amino acid code starts with the second DNA code.
					/// 6 	The DNA sequence is scanned right to left (backward direction for the complementary strand). The amino acid code starts with the third DNA code.
					/// 7 	Use each of the DNA translations of the codes 1; 2; 3.
					/// 8 	Use each of the DNA translations of the codes 4; 5; 6.
					/// 9 	Use each of the DNA translations of the codes 1; 2; 3; 4; 5; 6.
			UInt max_internal_cleavage_sites_;///< maximum number of internal cleavage sites
			UInt match_peak_count_;///< number of the top abundant peaks to match with theoretical one
			UInt match_peak_allowed_error_;///< number of peaks that may lack this test
			
			
			bool show_fragment_ions_;///< wether to display fragment ions
			bool print_duplicate_references_;///< whether all proteins containing a found peptide should be displayed
			bool remove_precursor_near_peaks_;///< whether peaks near (15 amu) the precursor peak are removed
			bool mass_type_parent_;///< mass type of the parent peak (0 - monoisotopic; 1 - average)
			bool mass_type_fragment_;///< mass type of fragment peaks (0 - monoisotopic; 1 - average)
			bool normalize_xcorr_;///< whether to display normalized xcorr values
			bool residues_in_upper_case_;///< whether residues are in upper case
			
			std::map< String, std::vector< String > > PTMname_residues_mass_type_;///< the modification names map to the affected residues, the mass change and the type
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_SEQUESTINFILE_H
