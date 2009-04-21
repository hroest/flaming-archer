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

#ifndef OPENMS_FORMAT_ANALYSISXMLFILE_H
#define OPENMS_FORMAT_ANALYSISXMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/AnalysisXMLHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/Identification.h>

namespace OpenMS
{
	/**
		@brief File adapter for AnalysisXML files

		If a critical error occurs due to the missing functionality, Exception::NotImplemented is thrown.

		@ingroup FileIO
	*/
	class OPENMS_DLLAPI AnalysisXMLFile
		:	public Internal::XMLFile,
			public ProgressLogger
	{
		public:
			///Default constructor
			AnalysisXMLFile();
			///Destructor
			virtual ~AnalysisXMLFile();

			/**
				@brief Loads a map from a AnalysisXML file.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			void load(const String& filename, Identification& id);

			/**
				@brief Stores a map in a AnalysisXML file.

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
			*/
			void store(const String& filename, const Identification& id) const;

			/**
				@brief Checks if a file is valid with respect to the mapping file and the controlled vocabulary.

				@param filename File name of the file to be checked.
				@param errors Errors during the validation are returned in this output parameter.
				@param warnings Warnings during the validation are returned in this output parameter.

				@exception Exception::FileNotFound is thrown if the file could not be opened
			*/
			bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings);

		private:

	};

} // namespace OpenMS

#endif // OPENMS_FOMAT_ANLYSISXMLFILE_H
