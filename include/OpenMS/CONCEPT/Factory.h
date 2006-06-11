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
// $Id: Factory.h,v 1.2 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Ole Schulz-Trieglaff  $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FACTORY_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FACTORY_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>
#include <iostream>

namespace OpenMS
{

  /** @brief  Returns FactoryProduct* based on the name of the desired concrete FactoryProduct
    
   		Factory class of the FeatureFinder - creates instances of all modules such as seeders,
   		extenders and fitters. 
   		
   		@ingroup Concept
  */
  template <typename FactoryProduct>
  class Factory
  {
    friend class singletonsNeedNoFriends; //some versions of gcc would warn otherwise

  private:
    /// \typedef Function signature of creator function 
    typedef FactoryProduct* (*FunctionType)();
    typedef std::map<std::string, FunctionType> Map;
    typedef typename Map::const_iterator MapIterator;

    /// destructor 
    virtual ~Factory(){}

    /// create with instance 
    Factory(){}

  public:
    /// singleton access to Factory 
    static Factory* instance()
    {
      if (!instance_ptr_){
				instance_ptr_ = new Factory();
				FactoryProduct::registerChildren();
      }
      return instance_ptr_;
    }

    /// return FactoryProduct according to unique identifier @p name  
    static FactoryProduct* create(const String& name)
    {
    	MapIterator it = instance()->inventory_.find(name);
      if (it != instance()->inventory_.end())
				return ( *(it->second) )();
      else 
      	throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				"This FactoryProduct is not registered!",name.c_str());
    }
    
    /** @brief register new concrete FactoryProduct 
     
       \param name unique name for concrete FactoryProduct
       \param creator default constructor for concrete FactoryProduct 
    */
    static void registerProduct(const String& name, const FunctionType creator)
    {
      instance()->inventory_[name] = creator;
    }

  private:

    Map inventory_;
    static Factory* instance_ptr_;
  };
  
  template<typename FactoryProduct>
    Factory<FactoryProduct>* Factory<FactoryProduct>::instance_ptr_ = 0;
    
}
#endif //OPENMS_TRANSFORMATIONS_FEATUREFINDER_FACTORY_H
