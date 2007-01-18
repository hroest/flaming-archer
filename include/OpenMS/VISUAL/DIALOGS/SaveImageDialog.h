
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_DIALOGS_SAVEIMAGEDIALOG_H
#define OPENMS_VISUAL_DIALOGS_SAVEIMAGEDIALOG_H

#include <OpenMS/config.h>

#include<qdialog.h>
#include<qcombobox.h>
#include<qlineedit.h>
#include<qcheckbox.h>

namespace OpenMS 
{
	/**
		@brief Dialog for saving an image.
		
		
		
		@ingroup Dialogs
	*/
	class SaveImageDialog : public QDialog
	{
		Q_OBJECT
		
		public:
			///Constructor
			SaveImageDialog( QWidget * parent = 0, const char * name = 0, bool modal = FALSE, WFlags f = 0 );
			///set size and size ratio
			void setSize(int x, int y);
			///accessors for the width
			int getXSize();
			///accessors for the height
			int getYSize();
			///accessors for the format
			QString getFormat();
		
		public slots:
			///changes width keeping proprotions
			void xSizeChanged(const QString& s);
			///changes height keeping proprotions
			void ySizeChanged(const QString& s);	
			///set size ratio when proportions checkbox is activated
			void proportionsActivated(bool state);	
			///checks if the values for width and heigth are ok before accepting the dialog
			void checkSize();	
		
		private:
			//format
			QComboBox* format_;
			//size
			QLineEdit* size_x_;
			QLineEdit* size_y_;
			QCheckBox* size_proportions_;
			//ratio size_x_/size_y_
			float size_ratio_;
	
			//set the size ratio (width/height)
			void setSizeRatio_(float r);
	};
}
#endif


