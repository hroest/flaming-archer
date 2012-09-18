// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer:  Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_TOOLDESCRIPTIONFILE_H
#define OPENMS_FORMAT_TOOLDESCRIPTIONFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/DATASTRUCTURES/ToolDescription.h>

namespace OpenMS
{
	/**
		@brief File adapter for ToolDescriptor files

		If a critical error occurs due to the missing functionality, Exception::NotImplemented is thrown.

		@ingroup FileIO
	*/
	class OPENMS_DLLAPI ToolDescriptionFile
		:	public Internal::XMLFile,
			public ProgressLogger
	{
		public:
			///Default constructor
			ToolDescriptionFile();
			///Destructor
			virtual ~ToolDescriptionFile();

			/**
				@brief Loads a map from a ToolDescriptor file.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
      void load(const String& filename, std::vector <Internal::ToolDescription>& tds);

			/**
				@brief Stores a map in a ToolDescriptor file.

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
			*/
			void store(const String& filename, const std::vector <Internal::ToolDescription>& tds) const;

		private:

	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_TOOLDESCRIPTIONFILE_H
