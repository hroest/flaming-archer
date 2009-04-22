// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_RTSIMULATION_H
#define OPENMS_SIMULATION_RTSIMULATION_H

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/SIMULATION/SimTypes.h>

namespace OpenMS {

  /**
   @brief 
   @ingroup Simulation
  */
  class OPENMS_DLLAPI RTSimulation
    : public DefaultParamHandler
  {

  public:
    /** @name Constructors and Destructors
      */
    //@{

    /// Constructor taking a random generator
    RTSimulation(const gsl_rng * random_generator);
    
    /// Copy constructor
    RTSimulation(const RTSimulation& source);

    /// Destructor
    virtual ~RTSimulation();
    //@}

    RTSimulation& operator = (const RTSimulation& source);
    
    /** 
     @brief Predict retention times for given features based on a SVM Model
     */
    void predict_rt(FeatureMap< > &);
 
    /**
     @brief Set retention times randomly for given contaminants
     */
    void predict_contaminants_rt(FeatureMap< > &);
    
    /**
     @brief Returns true if a RT column was simulated
     */
    bool isRTColumnOn();
  private:
    /// Default constructor -> hidden since we need to have a random generator
    RTSimulation();
    
    /// set defaults 
    void setDefaultParams_();

		// Name of the svm model file
		OpenMS::String rtModelFile_;
    
    /// total gradient time
    SimCoordinateType gradientTime_;
  protected:  
		/// Random number generator
		const gsl_rng* rnd_gen_;    
    
    /// Synchronize members with param class
		void updateMembers_();
    
  };

}

#endif
