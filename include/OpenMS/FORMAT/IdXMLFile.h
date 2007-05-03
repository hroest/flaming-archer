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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_IdXMLFile_H
#define OPENMS_FORMAT_IdXMLFile_H

#include <OpenMS/METADATA/Identification.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Used to load and store IdXML files
    
    This class is used to load and store documents that implement 
    the schema of IdXML files.
		
		A documented schema for this format can be found at http://open-ms.sourceforge.net/schemas/. 
		
  	@todo Document and upload schema (Nico)
  	
  	@note This format will only be used until the HUPO-PSI AnalysisXML format is finished!
  	
  	@ingroup FileIO
  */
  class IdXMLFile
  {
    public:
      /// Constructor
      IdXMLFile();
      
			/**
      	@brief Loads the identifications of an IdXML file

				The information is read in and the information is stored in the
				corresponding variables
      */
      void load(const String& filename, std::vector<ProteinIdentification>& protein_identifications, std::vector<IdentificationData>& id_data) const throw (Exception::FileNotFound, Exception::ParseError);
      					 
			/**
      	@brief Loads the identifications of an IdXML file

				The information is read in and the information is stored in the
				corresponding variables. This function offers the possibility to load 
				the predicted retention times which can be
				used to filter the peptides via the TOPP tool IDFilter.
      */
      void load(const String& filename, std::vector<ProteinIdentification>& protein_identifications, std::vector<IdentificationData>& id_data, std::map<String, DoubleReal>& predicted_retention_times)  	const throw (Exception::FileNotFound, Exception::ParseError);
      					 
			/**
      	@brief Stores the data in an IdXML file

				The data is read in and stored in the file 'filename'.
      */
      void store(String filename, const std::vector<ProteinIdentification>& protein_identifications, const std::vector<IdentificationData>& id_data) const throw (Exception::UnableToCreateFile); 

			/**
      	@brief Stores the data in an IdXML file

				The data is read in and stored in the file 'filename'.
      */
      void store(String filename, const std::vector<ProteinIdentification>& protein_identifications, const std::vector<IdentificationData>& id_data, const std::map<String, DoubleReal>& predicted_retention_times) 
      	const throw (Exception::UnableToCreateFile); 

  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_IdXMLFile_H
