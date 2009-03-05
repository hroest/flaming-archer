// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_INSPECTINFILE_H
#define OPENMS_FORMAT_INSPECTINFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
	/**
		@brief Inspect input file adapter.
		
		Creates a file that can be used for Inspect search from a peak list.
  	
  	@ingroup FileIO
	*/
  class OPENMS_DLLAPI InspectInfile
  {
		public:
			
			/// default constructor
			InspectInfile();

			/// copy constructor
			InspectInfile(const InspectInfile& inspect_infile);

			/// destructor
			virtual ~InspectInfile();

			/// assignment operator
			InspectInfile& operator=(const InspectInfile& inspect_infile);

			/// equality operator
			bool operator==(const InspectInfile& inspect_infile) const;

			/** stores the experiment data in an Inspect input file that can be used as input for Inspect shell execution
					@param filename set the given filename
					@throw UnableToCreateFile is thrown if the file could not be created
			*/
			void store(const String& filename);

			/** retrieves the name, mass change, affected residues, type and position for all modifications from a string
					@param modification_line
					@param modifications_filename
					@param monoisotopic if true, masses are considered to be monoisotopic
					@throw FileNotReadable if the modifications_filename could not be read
					@throw FileNotFound if modifications_filename could not be found
					@throw ParseError if modifications_filename could not be parsed
			*/
			void handlePTMs(const String& modification_line, const String& modifications_filename, const bool monoisotopic);
	    
			/**
				@brief Specifies a spectrum file to search.
				
				You can specify the name of a directory to search every file in that directory (non-recursively). Supported spectra file formats are .mzXML, .mzData, .ms2, dta, and .pkl. Multiple spectra in one .dta file are not supported.
			*/
			const String& getSpectra() const;
			void setSpectra(const String& spectra);
			
			/**
				@brief Specifies the name of a database (.trie file) to search.
				
				The .trie file contains one or more protein sequences delimited by asterisks, with no whitespace or other data. Use PrepDB.py (see above) to prepare a .trie file. Most .trie files have a corresponding .index file giving the names of the proteins. You can specify at most one database.
			*/
			const String& getDb() const;
			void setDb(const String& db);
			
			/// Specifies the name of a enzyme. "Trypsin", "None", and "Chymotrypsin" are the available values.
			const String& getEnzyme() const;
			void setEnzyme(const String& enzyme);
			
			/// Number of PTMs permitted in a single peptide.
			Int getModificationsPerPeptide() const;
			void setModificationsPerPeptide(Int modifications_per_peptide);
			
			/**
				@brief run Inspect in a blind mode
				
				If true, use the MS-Alignment algorithm to perform a blind search (allowing arbitrary modification masses). Running a blind search with one mod per peptide is slower than the normal (tag-based) search; running time is approximately 1 second per spectra per megabyte of database. Running a blind search with two mods is significantly slower. We recommend performing "blind" searches against a small database, containing proteins output by an earlier search.
			*/
			UInt getBlind() const;
			void setBlind(UInt blind);
			
			/**
				@brief the maximum modification size (in Da) to consider in a blind search
				
				Defaults to 200. Larger values require more time to search.
			*/
			Real getMaxPTMsize() const;
			void setMaxPTMsize(Real maxptmsize);
			
			/**
				@brief Specifies the parent mass tolerance, in Daltons.
				
				A candidate's flanking mass can differ from the tag's flanking mass by no more than this amount.
			*/
			Real getPrecursorMassTolerance() const;
			void setPrecursorMassTolerance(Real precursor_mass_tolerance);
			
			/**
				@brief How far b and y peaks can be shifted from their expected masses.
				
				Default is 0.5. Higher values produce a more sensitive but much slower search.
			*/
			Real getPeakMassTolerance() const;
			void setPeakMassTolerance(Real peak_mass_tolerance);
			
			/// If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible.
			UInt getMulticharge() const;
			void setMulticharge(UInt multicharge);
			
			/// If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass.
			const String& getInstrument() const;
			void setInstrument(const String& instrument);
			
			/// Number of tags to generate.
			Int getTagCount() const;
			void setTagCount(Int TagCount);

			/// return the modifications (the modification names map to the affected residues, the mass change and the type)
			const Map<String, std::vector<String> >& getModifications() const;
			
		private:
			
			String spectra_; ///< Specifies a spectrum file to search.

			String db_; ///< Specifies the name of a database (.trie file) to search. The .trie file contains one or more protein sequences delimited by asterisks, with no whitespace or other data.

			String enzyme_; ///< Specifies the name of a enzyme. "Trypsin", "None", and "Chymotrypsin" are the available values.

			Int modifications_per_peptide_; ///< allowed number of modifications per peptide

			UInt blind_; ///< If true, use the MS-Alignment algorithm to perform a blind search (allowing arbitrary modification masses). Running a blind search with one mod per peptide is slower than the normal (tag-based) search; running time is approximately 1 second per spectra per megabyte of database. Running a blind search with two mods is significantly slower. We recommend performing "blind" searches against a small database, containing proteins output by an earlier search. (The "Summary.py" script can be used to generate a second-pass database from initial search results)
			/// 0 - false, 1 - true, 2 - not set

			Real maxptmsize_; ///< For blind search, specifies the maximum modification size (in Da) to consider. Defaults to 200. Larger values require more time to search. <0 is not set

			Real precursor_mass_tolerance_; ///< Specifies the parent mass tolerance, in Daltons. A candidate's flanking mass can differ from the tag's flanking mass by no more than ths amount. <0 is not set

			Real peak_mass_tolerance_; ///< How far b and y peaks can be shifted from their expected masses. Default is 0.5. Higher values produce a more sensitive but much slower search. <0 is not set
			
			UInt multicharge_; ///< If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible.
			/// 0 - false, 1 - true, 2 - not set

			String instrument_; ///< If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass.

			Int tag_count_; ///< Number of tags to generate. <0 is not set
			
			Map<String, std::vector<String> > PTMname_residues_mass_type_;///< the modification names map to the affected residues, the mass change and the type
		
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_INSPECTINFILE_H
