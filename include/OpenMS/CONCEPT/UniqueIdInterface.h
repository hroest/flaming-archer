// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_UNIQUEIDINTERFACE_H
#define OPENMS_CONCEPT_UNIQUEIDINTERFACE_H

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

namespace OpenMS
{

/** @brief A base class defining a common interface for all classes having a unique id.

 Have a look at RichPeak2D for an example how to extend a class to support unique ids.

 @sa UniqueIdGenerator, RichPeak2D

 @ingroup Concept
 */
class OPENMS_DLLAPI UniqueIdInterface
{
  public:

    /**@brief This is the invalid unique id (cast it to a UInt64 if you like)

     It is represented as an \c enum because \c static class members lead to bugs and linker errors all the time...
     */
    enum
    {
      INVALID = 0
    };

    /**@brief Returns true if the unique_id is valid, false otherwise.

     Currently, an invalid unique id is represented by UInt64(0), but please prefer using this method for clarity.
     */
    inline static bool
    isValid(UInt64 unique_id)
    {
      return unique_id != INVALID;
    }

    /// Default constructor - the unique id will be <i>invalid</i>
    UniqueIdInterface()
    {
      clearUniqueId();
    }

    /// Copy constructor - copies the unique id
    UniqueIdInterface(const UniqueIdInterface& rhs) :
      unique_id_(rhs.unique_id_)
    {
    }

    /// Assignment operator - copies the unique id
    UniqueIdInterface &
    operator =(UniqueIdInterface const & rhs)
    {
      unique_id_ = rhs.unique_id_;
      return *this;
    }

    /// Destructor
    ~UniqueIdInterface()
    {
    }

    /// Equality comparison operator - the unique ids must be equal (!)
    bool
    operator ==(UniqueIdInterface const & rhs) const
    {
      return unique_id_ == rhs.unique_id_;
    }

    /// Non-mutable access to unique id - returns the unique id.
    UInt64
    getUniqueId() const
    {
      return unique_id_;
    }

    /// Clear the unique id.  The new unique id will be invalid.  Returns 1 if the unique id was changed, 0 otherwise.
    Size
    clearUniqueId()
    {
      if ( hasValidUniqueId() )
      {
        unique_id_ = 0;
        return 1;
      }
      else
        return 0;
    }

    void
    swap(UniqueIdInterface& from)
    {
      std::swap(unique_id_, from.unique_id_);
      return;
    }

    /// Returns whether the unique id is valid.  Returns 1 if the unique id is valid, 0 otherwise.
    Size
    hasValidUniqueId() const
    {
      return isValid(unique_id_);
    }

    /// Returns whether the unique id is invalid.  Returns 1 if the unique id is invalid, 0 otherwise.
    Size
    hasInvalidUniqueId() const
    {
      return !isValid(unique_id_);
    }

    /// Assigns a new, valid unique id.  Always returns 1.
    Size
    setUniqueId()
    {
      unique_id_ = UniqueIdGenerator::getUniqueId();
      return 1;
    }

    /// Assigns a valid unique id, but only if the present one is invalid.  Returns 1 if the unique id was changed, 0 otherwise.
    Size
    ensureUniqueId()
    {
      if ( !hasValidUniqueId() )
      {
        unique_id_ = UniqueIdGenerator::getUniqueId();
        return 1;
      }
      else
        return 0;
    }

    /// Assigns the given unique id.
    void
    setUniqueId(UInt64 rhs)
    {
      unique_id_ = rhs;
    }

    /**@brief Mutable access to unique id.

     This is designed to work well with \c id attributes in XML, which are prefixed with letters.
     The portion of the String after the last underscore is extracted and parsed as a UInt64.  <i>It must consist of digits only.</i>
     For example, <code>some_feature.setUniqueId("f_12345_00067890")</code> is equivalent to <code>some_feature.setUniqueId(67890)</code>

     */
    void
    setUniqueId(const String & rhs);

  protected:
    /// the unique id
    UInt64 unique_id_;

};

} //namespace OpenMS

#endif // OPENMS_CONCEPT_UNIQUEIDINTERFACE_H
