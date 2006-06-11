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
// $Id: FactoryProductView.h,v 1.3 2006/03/28 10:07:17 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
#ifndef OPENMS_VISUAL_FACTORYPRODUCTVIEW_H
#define OPENMS_VISUAL_FACTORYPRODUCTVIEW_H

#include <qlabel.h>

#include <OpenMS/CONCEPT/FactoryProduct.h>

namespace OpenMS
{

  /**
  displays the params of a FactoryProduct<br>
  */
  class FactoryProductView : public QLabel
  {
    Q_OBJECT
  public:
    FactoryProductView(QWidget* = 0, const char* name = 0);
    ~FactoryProductView();
  public slots:
    void displayFactoryProduct(const FactoryProduct* const conf);
  };

}

#endif //OPENMS_VISUAL_FACTORYPRODUCTVIEW_H
