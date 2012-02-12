// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_DIALOGS_TOPPASVERTEXNAMEDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPASVERTEXNAMEDIALOG_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/UIC/ui_TOPPASVertexNameDialog.h>

namespace OpenMS 
{
	/**
		@brief Dialog which allows to change the name of an input vertex
		
		@ingroup TOPPAS_elements
		@ingroup Dialogs
	*/
	class OPENMS_GUI_DLLAPI TOPPASVertexNameDialog
		: public QDialog,
			public Ui::TOPPASVertexNameDialogTemplate
	{
		Q_OBJECT
				
		public:
			
			/// Constructor
			TOPPASVertexNameDialog(const QString& name);
			
			/// Returns the name
			QString getName();
	};
	
}
#endif // OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILESDIALOG_H
