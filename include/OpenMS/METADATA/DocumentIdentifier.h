// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_DOCUMENTIDENTIFIER_H
#define OPENMS_METADATA_DOCUMENTIDENTIFIER_H

// OpenMS
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
//~ #include <OpenMS/FORMAT/FileHandler.h>

namespace OpenMS
{
  /**
     @brief Manage source document information.

     This class stored information about the source document.
     Primarily this is the document id e.g. a LSID.

     For source files additional information can be stored:
     - file name
     - file type

     @improvement change file_type_ to FileHandler::Type. Inclusion of FileHandler results in circular inclusions. Forward declarations can't resolve (inheritance, no enum forward declaration with gcc, ...). Possible solution is extraction of FileHandler::Type in extra header file.

     @ingroup Metadata
  */
  class DocumentIdentifier
  {
    public:
      /** @name Constructors and Destructors
      */
      //@{
      /// default constructor
      DocumentIdentifier();

      /// Copy constructor
      DocumentIdentifier(const DocumentIdentifier& source);

      /// Assignment operator
      DocumentIdentifier& operator=(const DocumentIdentifier& source);

			/// Equality operator
			bool operator== (const DocumentIdentifier& rhs) const;

      /// destructor
      virtual ~DocumentIdentifier();
      //@}

      /** @name Acessors
       */
      //@{

			/// retrieve computed zero-charge feature map
      void setIdentifier(const String& id);

      /// retrieve computed zero-charge feature map
      const String& getIdentifier() const;

			/// exchange content with @p from
			void swap(DocumentIdentifier& from);


      /// set the file_name_ according to absolute path of the file loaded from preferrably done whilst loading
      void setLoadedFilePath(const String& file_name);

      /// get the file_name_ which is the absolute path to the file loaded from
      const String& getLoadedFilePath() const;

      /// set the file_type according to the type of the file loaded from (see FileHandler::Type) preferrably done whilst loading
      void setLoadedFileType(const String& file_name);

      /// get the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded from
      const String& getLoadedFileType() const;

      //@}

    protected:
      /// the ID (e.g. LSID)
      String id_;
      /// the path to the loaded file
      String file_path_;
      /// the type name of the loaded file
      String file_type_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_DOCUMENTIDENTIFIER_H
