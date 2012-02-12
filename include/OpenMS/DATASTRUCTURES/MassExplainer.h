// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_MASSEXPLAINER_H
#define OPENMS_DATASTRUCTURES_MASSEXPLAINER_H

#include <OpenMS/DATASTRUCTURES/Adduct.h>

namespace OpenMS
{
  class Compomer;
  class Adduct;

  /**
    @brief computes empirical formulas for given mass differences using a set of allowed elements
    
  	
  	@ingroup Datastructures
  */
  class OPENMS_DLLAPI MassExplainer
  {

	 public:
			 
		typedef Adduct::AdductsType AdductsType; //vector<Adduct>
		typedef std::vector<Compomer>::const_iterator CompomerIterator;
		
    ///@name Constructors and destructor
    //@{
		/// Default constructor
		MassExplainer();
		
		/// Constructor
		MassExplainer(AdductsType adduct_base);

		/// Constructor
		MassExplainer(Int q_min, Int q_max, Int max_span, DoubleReal thresh_logp);

		/// Constructor
		MassExplainer(AdductsType adduct_base, Int q_min, Int q_max, Int max_span, DoubleReal thresh_logp, Size max_neutrals);

		
	private:
		/// check consistency of input
		/// @param init_thresh_p set default threshold (set to "false" to keep current value)
		void init_(bool init_thresh_p);
	public:
		/// Assignment operator
		MassExplainer& operator = (const MassExplainer& rhs);
		
		/// Destructor
		virtual ~MassExplainer();
    //@}
    
		
		/// fill map with possible mass-differences along with their explanation
		void compute();
		
		
    //@name Accessors
    //@{
		
		/// Sets the set of possible adducts
		void setAdductBase(AdductsType adduct_base);
		/// Returns the set of adducts
		AdductsType getAdductBase() const;
		
		/// return a compomer by its Id (useful after a query() ).
		const Compomer& getCompomerById(Size id) const;
		//@}
				
		
		/// search the mass database for explanations
		/// @param net_charge       net charge of compomer seeked
		/// @param mass_to_explain  mass in Da that needs explanation
		/// @param mass_delta       allowed deviation from exact mass
		/// @param thresh_log_p		  minimal log probability required
		/// @param firstExplanation begin range with candidates according to net_charge and mass
		/// @param lastExplanation  end range
		SignedSize query(const Int net_charge,
										 const float mass_to_explain,
										 const float mass_delta,
										 const float thresh_log_p,
										 std::vector< Compomer >::const_iterator & firstExplanation,
										 std::vector< Compomer >::const_iterator & lastExplanation) const;
	 protected:
    
	  ///check if the generated compomer is valid judged by its probability, charges etc
		bool compomerValid_(const Compomer& cmp);

		/// create a proper adduct from formula and charge and probability
		Adduct createAdduct_(const String & formula, const Int charge, const DoubleReal p) const;
			 
		/// store possible explanations (as formula) for a certain ChargeDifference and MassDifference
		std::vector< Compomer > explanations_;	 
		/// all allowed adducts, whose combination explains the mass difference
		AdductsType adduct_base_;
		/// minimal expected charge
		Int q_min_;
		/// maximal expected charge
		Int q_max_;
		/// maximal span (in terms of charge) for co-features, e.g. a cluster with q={3,6} has span=4
		Int max_span_;
		/// minimum required probability of a compound (all other compounds are discarded)
		DoubleReal thresh_p_;
		/// Maximum number of neutral(q=0) adducts
		Size max_neutrals_;
		
  };

			
} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_MASSEXPLAINER_H

 
