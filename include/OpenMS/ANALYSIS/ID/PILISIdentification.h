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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_ID_PILISIDENTIFICATION_H
#define OPENMS_ANALYSIS_ID_PILISIDENTIFICATION_H

#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/ANALYSIS/ID/PILISSequenceDB.h>

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>

namespace OpenMS
{
	/**
	  @brief This class actually implements a complete identification run with PILIS

		The PILISIdentification class needs a PILISModel and a PILISSequenceDB to generate
		identifications. Simply call getIdentifications with a PeakMap.
	*/
	// forward declarations
	class CompareFunctor;
	
	class PILISIdentification
	{

		public:

			/** @name constructors and destructors
			 */
			//@{
			/// default constructor
			PILISIdentification();
			
			/// copy constructor
			PILISIdentification(const PILISIdentification& source);
			
			/// destructor
			virtual ~PILISIdentification();
			//@}
		
			///
			const PILISIdentification& operator = (const PILISIdentification&);

			/** @name Accessors
			 */
			//@{
			/// sets the sequence DB to be used for the identification runs
			void setSequenceDB(PILISSequenceDB* sequence_db);

			/// sets the model to be used for the identification run
			void setModel(PILISModel* hmm_model);

			/// performs an identification run on a PeakMap
			void getIdentifications(std::vector<Identification>& ids, const PeakMap& exp);

			/// performs an identification run on a PeakSpectrum
			void getIdentification(Identification& id, const PeakSpectrum& spectrum);

			/// mutable access to the parameters
			Param& getParam();

			/// non-mutable access to the parameters
			const Param& getParam() const;
			
			/// sets the parameters
			void setParam(const Param& param);
			//@}

		protected:

			/// fast method to create spectra for pre scoring
			void getSpectrum_(PeakSpectrum& spec, const String& sequence, int charge);
		
			/// 
			void getPreIdentification_(Identification& id, const PeakSpectrum& spec, const std::vector<PILISSequenceDB::PepStruct>& cand_peptides, CompareFunctor* scorer);

			///
			void getFinalIdentification_(Identification& id, const PeakSpectrum& spec, const Identification& pre_id, CompareFunctor* scorer);
			
			Param param_;

			PILISSequenceDB* sequence_db_;

			PILISModel* hmm_model_;

			HashMap<char, double> aa_weight_;

			Peak p_;
	};
}

#endif
