// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: MascotXMLFile.h,v 1.1 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MASCOTXMLFILE_H
#define OPENMS_FORMAT_MASCOTXMLFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Identification.h>
#include <fstream>

#include <qxml.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Used to load and store analysisXML files
    
    This class is used to load and store documents that implement 
    the schema of analysisXML files.
  
  	@ingroup FileIO
  */
  class MascotXMLFile
  {
    public:
      /// Constructor
      MascotXMLFile();
      /// Copy constructor
      MascotXMLFile(const MascotXMLFile& source);
      /// Destructor
      ~MascotXMLFile();
      
      void load(const String& 								filename,
      					ProteinIdentification*				protein_identification, 
      					std::vector<Identification>* 	identifications, 
      					std::vector<float>* 					precursor_retention_times,
      					std::vector<float>* 					precursor_mz_values)  	
      							const throw (Exception::FileNotFound, 
  							 								 Exception::FileNotReadable, 
  							 								 Exception::FileEmpty,
  							 								 Exception::ParseError);
      					 
  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_MASCOTXMLFILE_H
