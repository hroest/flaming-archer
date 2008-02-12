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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_INSPECTOUTFILE_H
#define OPENMS_FORMAT_INSPECTOUTFILE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

namespace OpenMS
{
	/**
		@brief Representation of an Inspect outfile
		
		This class serves to read in an Inspect outfile and write an IdXML file
		
  	@todo Handle Modifications (Andreas)
		
		@ingroup FileIO
	*/
	class InspectOutfile
	{
		public:
			/// default constructor
			InspectOutfile();

			/// copy constructor
			InspectOutfile(const InspectOutfile& inspect_outfile);

			/// destructor
			virtual ~InspectOutfile();

			/// assignment operator
			InspectOutfile& operator=(const InspectOutfile& inspect_outfile);

			/// equality operator
			bool operator==(const InspectOutfile& inspect_outfile) const;
			
			/// load the results of an Inspect search
			std::vector< UInt > load(const String& result_filename, std::vector< PeptideIdentification >& peptide_identifications, ProteinIdentification& protein_identification, const Real p_value_threshold, const String& database_filename = "") throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument, Exception::FileEmpty);
			
			std::vector< UInt > getWantedRecords(const String& result_filename, Real p_value_threshold) throw (Exception::FileNotFound, Exception::FileEmpty, Exception::IllegalArgument);

			/// generates a trie database from another one, using the wanted records only
			void compressTrieDB(const String& database_filename, const String& index_filename, std::vector< UInt >& wanted_records, const String& snd_database_filename, const String& snd_index_filename, bool append = false) throw (Exception::FileNotFound, Exception::ParseError, Exception::UnableToCreateFile);

			/// generates a trie database from a given one (the type of database is determined by getLabels)
			void generateTrieDB(const String& source_database_filename, const String& database_filename, const String& index_filename, bool append = false, const String species = "") throw (Exception::FileNotFound, Exception::UnableToCreateFile);
			

			/// retrieve the accession type and accession number from a protein description line
			/// (e.g. from FASTA line: >gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus], get ac:AAD44166.1 ac type: GenBank)
			void getACAndACType(String line, String& accession, String& accession_type);

			/// retvrieve the precursor retention time and mz value
			void getPrecursorRTandMZ(const std::vector< std::pair< String, std::vector< std::pair< UInt, UInt > > > >& files_and_peptide_identification_with_scan_number, std::vector< PeptideIdentification >& ids) throw(Exception::ParseError);

			/// retrieve the labes of a given database (at the moment FASTA and Swissprot)
			void getLabels(const String& source_database_filename, String& ac_label, String& sequence_start_label, String& sequence_end_label, String& comment_label, String& species_label) throw (Exception::FileNotFound, Exception::ParseError);

			/// retrieve sequences from a trie database
			std::vector< UInt > getSequences(const String& database_filename, const std::map< UInt, UInt >& wanted_records, std::vector< String >& sequences) throw (Exception::FileNotFound);

			/// get the experiment from a file
			template< typename PeakT > void getExperiment(MSExperiment< PeakT >& exp, String& type, const String& in_filename) throw(Exception::ParseError);

			/// get the search engine and its version from a file with the output of InsPecT without parameters
			void getSearchEngineAndVersion(const String& inspect_output_without_parameters_filename, ProteinIdentification& protein_identification) throw (Exception::FileNotFound);

			/// read the header of an inspect output file and retrieve various informations
			void readOutHeader(const String& filename, const String& header_line, Int& spectrum_file_column, Int& scan_column, Int& peptide_column, Int& protein_column, Int& charge_column, Int& MQ_score_column, Int& p_value_column, Int& record_number_column, Int& DB_file_pos_column, Int& spec_file_pos_column, UInt& number_of_columns) throw (Exception::ParseError);

		protected:
			/// a record in the index file that belongs to a trie database consists of three parts
			/// 1) the protein's position in the original database
			/// 2) the proteins's position in the trie database
			/// 3) the name of the protein (the line with the accession identifier)
			static const UInt db_pos_length_; ///< length of 1)
			static const UInt trie_db_pos_length_; ///< length of 2)
			static const UInt protein_name_length_; ///< length of 3)
			static const UInt record_length_; ///< length of the whole record
			static const char trie_delimiter_; ///< the sequences in the trie database are delimited by this character
			static const String score_type_;///< type of score
	};
	
} //namespace OpenMS

#endif // OPENMS_FORMAT_INSPECTOUTFILE_H
