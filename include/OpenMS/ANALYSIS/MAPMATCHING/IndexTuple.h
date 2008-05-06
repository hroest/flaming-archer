// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_INDEXTUPLE_H
#define OPENMS_ANALYSIS_MAPMATCHING_INDEXTUPLE_H

#include <iostream>
#include <vector>

#include <OpenMS/DATASTRUCTURES/DPosition.h>

namespace OpenMS
{
  /**
    @brief This class stores some needful information of an element.
    
    The IndexTuple class is used during map matching. 
    It stores next to an element's index (within a container), a pointer to the element itself,
    an index of the map it is contained as well as the transformed position of the element.
    
    
            
  */
  //TODO Derive from RawDataPoint2D, PeakIndex and add map_index member
  class IndexTuple
  {
    public:
      typedef DPosition<2> PositionType;
      
      /// Default constructor
      IndexTuple()
          : map_index_(0),
          element_index_(0),
          intensity_(0)
      {}

      /// Constructor
      inline IndexTuple(UInt map_index, UInt element_index, DoubleReal intensity, DPosition<2> transformed_pos)
      {
        map_index_ = map_index;
        element_index_ = element_index;
        intensity_ = intensity;
        transformed_position_ = transformed_pos;
      }

      /// Copy constructor
      inline IndexTuple(const IndexTuple& source)
      {
        map_index_ = source.map_index_;
        element_index_ = source.element_index_;
        intensity_ = source.intensity_;
        transformed_position_ = source.transformed_position_;
      }

      /// Assignment operator
      IndexTuple& operator = (const IndexTuple& source)
      {
        if (&source == this)
          return *this;

        map_index_ = source.map_index_;
        element_index_ = source.element_index_;
        intensity_ = source.intensity_;
        transformed_position_ = source.transformed_position_;
        return *this;
      }

      /// Destructor
      virtual ~IndexTuple()
      {}
      
      /// Non-mutable access to the container index
      inline UInt getMapIndex() const
      {
        return map_index_;
      }
      
      /// Set the container index
      inline void setMapIndex(UInt c)
      {
        map_index_ = c;
      }

      /// Non-mutable access to the element index
      inline UInt getElementIndex() const
      {
        return element_index_;
      }
      
      /// Set the element index
      inline void setElementIndex(UInt e)
      {
        element_index_= e;
      }

      /// Non-mutable access to the element
      inline DoubleReal getIntensity() const
      {
        return intensity_;
      }
      
      /// Set the element
      inline void setIntensity(DoubleReal intensity)
      {
        intensity_ = intensity;
      }

      /// Non-mutable access to the transformed position
      inline const PositionType& getTransformedPosition() const
      {
        return transformed_position_;
      }
      
      /// Set the transformed position
      inline void setTransformedPosition(const PositionType& p)
      {
        transformed_position_ = p;
      }

      /// Equality operator
      virtual bool operator == (const IndexTuple& i) const
      {
        return ((map_index_ == i.map_index_) && (element_index_ == i.element_index_) && (intensity_ == i.intensity_));
      }

      /// Equality operator
      virtual bool operator != (const IndexTuple& i) const
      {
        return !((map_index_ == i.map_index_) && (element_index_ == i.element_index_) && (intensity_ == i.intensity_));
      }

      /// Compare by getOverallQuality()
      struct IndexLess
            : std::binary_function < IndexTuple, IndexTuple, bool >
      {
        inline bool operator () ( IndexTuple const & left, IndexTuple const & right ) const
        {
          return ( left.map_index_ < right.map_index_ );
        }
      };

    protected:
      /// Transformed element position
      PositionType transformed_position_;
      /// Int of the element's container
      UInt map_index_;
      /// Int of the element within element's container
      UInt element_index_;
      /// Pointer to the element itself
      DoubleReal intensity_;
  };

  ///Print the contents to a stream.
  template < typename ContainerT >
  std::ostream& operator << (std::ostream& os, const IndexTuple& cons)
  {
    os << "---------- IndexTuple -----------------\n"
    << "Transformed Position: " << cons.getTransformedPosition() << '\n'
    << "Intensity: " << cons.getIntensity() << '\n'
    << "Element Index: " << cons.getElementIndex() << '\n'
    << "Map Index: " << cons.getMapIndex() << std::endl;
    return os;
  }
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMAPPING_INDEXTUPLE_H
