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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_PEPTIDEHIT_H
#define OPENMS_METADATA_PEPTIDEHIT_H

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{
  /**
    @brief Representation of a peptide hit
    
    It contains the fields score, score_type, rank, and sequence.
		
		@ingroup Metadata
  */
  class PeptideHit
  	: public MetaInfoInterface
  {
  	public:

			/// @name Comparators for PeptideHit and ProteinHit
			//@{
			/// Greater predicate for scores of hits
			class ScoreMore
			{
			  public:
			  	template<typename Arg>
			    bool operator()(const Arg& a, const Arg& b)
			    {
			      return a.getScore() > b.getScore();
			    }
			};
			
			/// Lesser predicate for scores of hits
			class ScoreLess
			{
			  public:
			  	template<typename Arg>
			    bool operator()(const Arg& a, const Arg& b)
			    {
			      return a.getScore() < b.getScore();
			    }
			};
			//@}

			/**	@name Constructors and Destructor */
			//@{
			/// default constructor
	    PeptideHit();
	    
			/// values constructor
	    PeptideHit(DoubleReal score, 
	    					 UInt rank, 
								 Int charge,
	    					 const String& sequence);
	
			/// copy constructor
	    PeptideHit(const PeptideHit& source);
					
			/// destructor
	    virtual ~PeptideHit();
	    //@}
	    
			/// assignment operator
	    PeptideHit& operator=(const PeptideHit& source);
	
			/// Equality operator
			bool operator == (const PeptideHit& rhs) const;
			
			/// Inequality operator
			bool operator != (const PeptideHit& rhs) const;
	
			/**	@name Accessors 
			*/
			//@{
	    /// returns the score of the peptide hit 
	    Real getScore() const;
	    
			/// returns the rank of the peptide hit
	    UInt getRank() const;
			
			/// returns the peptide sequence without trailing or following spaces
	  	String getSequence() const;
			
			/// returns the carge of the peptide
			Int getCharge() const;
			
			/// returns the corresponding protein accessions
			const std::vector<String>& getProteinAccessions() const;
	
			/// sets the corresponding protein accessions
		  void setProteinAccessions(const std::vector<String>& accessions);
	    
			/// sets the score of the peptide hit 
	    void setScore(DoubleReal score);
	    
			/// sets the rank
	    void setRank(UInt newrank);
	    
			/// sets the peptide sequence
			void setSequence(const String& sequence);
	
			/// sets the charge of the peptide
			void setCharge(Int charge);
			
			/// adds a accession of a protein which contains this peptide hit
			void addProteinAccession(const String& accession); 
	    //@}
	
	  
		protected:
	    Real score_;									///< the score of the peptide hit
			UInt rank_;    				///< the position(rank) where the hit appeared in the hit list
			Int charge_; ///< the charge of the peptide
	    String sequence_;							///< the amino acid sequence of the peptide hit 
	    std::vector<String> corresponding_protein_accessions_; ///< the accessions of the corresponding proteins
	
	};

} // namespace OpenMS

#endif // OPENMS_METADATA_PEPTIDEHIT_H
