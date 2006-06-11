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
// $Id: XMLHandler.h,v 1.7 2006/03/28 18:51:06 j-joachim Exp $
// $Author: j-joachim $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_XMLHANDLER_H

#include <OpenMS/CONCEPT/Types.h>

#include <qxml.h>

namespace OpenMS
{
	namespace Internal
	{
	/**
		@brief Base class for XML handlers.
		
		This class extends the QXmlDefaultHandler by some functionality for the handling of errors.
	*/
  class XMLHandler
  	: public QXmlDefaultHandler
  {
    public:
      XMLHandler();

      virtual ~XMLHandler();

      bool error(const QXmlParseException& exception);

      bool fatalError(const QXmlParseException& exception);

      bool warning(const QXmlParseException& exception);

  		QString errorString();

		  virtual bool characters( const QString & chars );

      virtual bool startElement(const QString & uri, const QString & local_name, 
																const QString & qname, const QXmlAttributes & attributes );

      virtual bool endElement( const QString & uri, const QString & local_name,
															 const QString & qname ); 

			void setUseWarnings(bool doUse);
			bool useWarnings();

  	protected:

			QString error_message_;
			QString file_;
			bool no_error_;
	
			/** @brief use QXml-warnings to show unhandled tags or values */
			bool use_warnings_;
	
			inline SignedInt asSignedInt_(const QString& in)
			{
				bool ok = true;
				SignedInt res = in.toInt(&ok);
				if (!ok){
					no_error_ = false;
					error_message_  = QString("SignedInt conversion error of \"%1\" parsed by %2 ").arg(in).arg(file_);
				}
				return res;
			}
	
			inline UnsignedInt asUnsignedInt_(const QString& in)
			{
				bool ok = true;
				UnsignedInt res = in.toUInt(&ok);
				if (!ok){
					no_error_ = false;
					error_message_  = QString("UnsignedInt conversion error of \"%1\" parsed by %2 ").arg(in).arg(file_);
				}
				return res;
			}
	
	 		inline double asDouble_(const QString& in)
			{
				bool ok = true;
				double res = in.toDouble(&ok);
				if (!ok){
					no_error_ = false;
					error_message_  = QString("Double conversion error of \"%1\" parsed by %2 ").arg(in).arg(file_);
				}
				return res;
			}
	
	
	 		inline float asFloat_(const QString& in)
			{
				bool ok = true;
				double res = in.toFloat(&ok);
				if (!ok){
					no_error_ = false;
					error_message_  = QString("Float conversion error of \"%1\" parsed by %2 ").arg(in).arg(file_);
				}
				return res;
			}
	
	 		inline bool asBool_(const QString& in)
			{
				if (in == "true" || in == "1") return true;
				else if (in == "false" || in == "0") return false;
				else {
					no_error_ = false;
					error_message_  = QString("Bool conversion error of \"%1\" parsed by %2 ").arg(in).arg(file_);
				}
				return false;
			}

	};

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
