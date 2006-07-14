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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_PROTEINHIT_H
#define OPENMS_METADATA_PROTEINHIT_H

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /**
    @brief Representation of a protein hit
    
    It contains the fields score, score_type, rank, accession, 
    accession_type and sequence.
		
		@ingroup Metadata
  */
  class ProteinHit
  {
  	public:
		
		/**	@name Constructors and Destructor */
		//@{
		
		/// default constructor
    ProteinHit();
    
		/// values constructor
    ProteinHit(double score, 
    					 std::string score_type, 
    					 uint rank, 
    					 String accession, 
    					 std::string accession_type, 
    					 String sequence);

		/// copy constructor
    ProteinHit(const ProteinHit& source);
				
		/// destructor
    ~ProteinHit();
    //@}
    
		/// assignment operator
    ProteinHit& operator=(const ProteinHit& source);
		/// Equality operator
		bool operator == (const ProteinHit& rhs) const;
		/// Inequality operator
		bool operator != (const ProteinHit& rhs) const;


		/**	@name Accessors */
		//@{
		
    /// returns the score of the protein hit 
    float getScore() const;
    /// returns the type of the score
    const std::string& getScoreType() const;
		/// returns the rank of the protein hit
    UnsignedInt getRank() const;
		/// returns the protein sequence
  	String getSequence() const;
		/// returns the accession of the protein
		String getAccession() const;
		/// returns the type of the accession string of the protein
		std::string getAccessionType() const;   	

    /// sets the score of the protein hit 
    void setScore(const double& score);
    /// sets the type of the score
    void setScoreType(const std::string& score_type);
		/// sets the rank
    void setRank(uint newrank);
		/// sets the protein sequence
		void setSequence(const String& sequence);
		/// sets the accession of the protein
		void setAccession(const String& accession);
		/// sets the type of the accession string of the protein
		void setAccessionType(const std::string& accession_type);  	
				
    //@}

		/// clears all information of the protein hit
    void clear();
  protected:
    float score_;									///< the score of the protein hit
    std::string score_type_;    	///< the score type of the protein hit 
		UnsignedInt rank_;    				///< the position(rank) where the hit 
																	///< appeared in the hit list
    String accession_;						///< the protein identifier
    std::string accession_type_;	///< the type of the accession
    String sequence_;							///< the amino acid sequence of the protein hit

  };

} // namespace OpenMS

#endif // OPENMS_METADATA_PROTEINHIT_H
