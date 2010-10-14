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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_RAWMSSIGNALSIMULATION_H
#define OPENMS_SIMULATION_RAWMSSIGNALSIMULATION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/SIMULATION/SimTypes.h>
#include <OpenMS/SIMULATION/EGHModel.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>

namespace OpenMS {

  class IsotopeModel;

  /**
   @brief Simulates MS signales for a given set of peptides

   Simulates MS signales for a given set of peptides, with charge annotation,
   given detectabilities, predicted retention times and charge values.

   @htmlinclude OpenMS_RawMSSignalSimulation.parameters

   @ingroup Simulation
  */
  class OPENMS_DLLAPI RawMSSignalSimulation
    : public DefaultParamHandler,
      public ProgressLogger
  {

  public:
    /** @name Constructors and Destructors
      */
    //@{
    /// Constructor taking a random generator
    RawMSSignalSimulation(const SimRandomNumberGenerator& rng);

    /// Copy constructor
    RawMSSignalSimulation(const RawMSSignalSimulation& source);

    /// Destructor
    virtual ~RawMSSignalSimulation();
    //@}

    RawMSSignalSimulation& operator = (const RawMSSignalSimulation& source);

    /// load the contaminants from contaminants:file param
    /// You do not have to call this function before calling generateRawSignals(), but it might 
    /// be useful to check if the contaminant file is valid
    void loadContaminants();

    /// fill experiment with signals and noise
    void generateRawSignals(FeatureMapSim & features, MSSimExperiment & experiment, FeatureMapSim & contaminants);

  protected:
    /// Default constructor
    RawMSSignalSimulation();

    /// Synchronize members with param class
		void updateMembers_();

    /// Set default parameters
    void setDefaultParams_();

    /**
     @brief Add a 1D signal for a single feature

     @param feature The feature which should be simulated
     @param experiment The experiment to which the simulated signals should be added
     */
    void add1DSignal_(Feature & feature, MSSimExperiment & experiment);

    /**
     @brief Add a 2D signal for a single feature

     @param feature The feature which should be simulated
     @param experiment The experiment to which the simulated signals should be added
     */
    void add2DSignal_(Feature & feature, MSSimExperiment & experiment);

    /**
     @brief Samples signales for the given 1D model

     @param iso The isotope model from which the signales will be sampled
     @param mz_start Start coordinate (in m/z dimension) of the region where the signals will be sampled
     @param mz_end End coordinate (in m/z dimension) of the region where the signals will be sampled
     @param experiment Experiment to which the sampled signales will be added
     @param activeFeature The current feature that is simulated
     */
    void samplePeptideModel1D_(const IsotopeModel & iso,
                               const SimCoordinateType mz_start,  const SimCoordinateType mz_end,
                               MSSimExperiment & experiment, Feature & activeFeature);

    /**
     @brief Samples signales for the given 2D model

     @param pm The product model from which the signales will be sampled
     @param mz_start Start coordinate (in m/z dimension) of the region where the signals will be sampled
     @param mz_end End coordinate (in m/z dimension) of the region where the signals will be sampled
     @param rt_start Start coordinate (in rt dimension) of the region where the signals will be sampled
     @param rt_end End coordinate (in rt dimension) of the region where the signals will be sampled
     @param experiment Experiment to which the sampled signales will be added
     @param activeFeature The current feature that is simulated
     */
    void samplePeptideModel2D_(const ProductModel<2> & pm,
                             const SimCoordinateType mz_start,  const SimCoordinateType mz_end,
                             SimCoordinateType rt_start, SimCoordinateType rt_end,
                             MSSimExperiment & experiment, Feature & activeFeature);

    /**
     @brief Add the correct Elution profile to the passed ProductModel
     */
    void chooseElutionProfile_(EGHModel*& elutionmodel, Feature & feature, const double scale, const DoubleReal rt_sampling_rate, const MSSimExperiment & experiment);

    /**
     @brief build contaminant feature map
    */
    void createContaminants_(FeatureMapSim & contaminants, MSSimExperiment & exp);

    /// Add shot noise to the experimet
    void addShotNoise_(MSSimExperiment & experiment, SimCoordinateType minimal_mz_measurement_limit, SimCoordinateType maximal_mz_measurement_limit);

    /// Add white noise to the experiment
    void addWhiteNoise_(MSSimExperiment & experiment);

    /// Add a base line to the experiment
    void addBaseLine_(MSSimExperiment & experiment, SimCoordinateType minimal_mz_measurement_limit);

    /// Compress signales in a single RT scan given the m/z bin size
    void compressSignals_(MSSimExperiment & experiment);

    /// worker function for compressSignals_
    Size compressSignalsRun_(MSSimExperiment & experiment);

		/// bin size in m/z dimension
		SimCoordinateType mz_sampling_rate_;

		/// Mean of peak m/z error
		SimCoordinateType mz_error_mean_;
		/// Standard deviation of peak m/z error
		SimCoordinateType mz_error_stddev_;

    /**
     * @brief Computes a rescaled feature intensity based on the set parameters for feature intensity scaling and the passed parameter @p natural_scaling_factor.
     *
     * @param feature_intensity Intensity of the current feature.
     * @param natural_scaling_factor Additional scaling factor used by some of the sampling models.
     *
     * @return Rescaled feature intensity.
     */
    SimIntensityType getFeatureScaledIntensity_(const SimIntensityType feature_intensity, const SimIntensityType natural_scaling_factor);

    /// Scaling factor of peak intensities
    SimIntensityType intensity_scale_;
    /// Standard deviation of peak intensity scaling
    SimIntensityType intensity_scale_stddev_;

	  /// Full width at half maximum of simulated peaks
		SimCoordinateType peak_std_;

		/// Random number generator
    SimRandomNumberGenerator const * rnd_gen_;

    enum IONIZATIONMETHOD {IM_ESI=0,IM_MALDI=1,IM_ALL=2};
    enum PROFILESHAPE {RT_RECTANGULAR, RT_GAUSSIAN};

    struct ContaminantInfo
    {
      String name;
      EmpiricalFormula sf;
      DoubleReal rt_start, rt_end, intensity;
      Int q;
      PROFILESHAPE shape;
      IONIZATIONMETHOD im;
    };

    std::vector<ContaminantInfo> contaminants_;

    /**
    @p threaded_random_numbers keeps a set of random numbers for each thread simulating a feature.
      */
    std::vector<std::vector<double> > threaded_random_numbers;

    /**
      Indicates which random numbers each thread has used already and if the random number pool
      should be rebuild.
      */
    std::vector< Size > threaded_random_numbers_index;

    static const Size THREADED_RANDOM_NUMBER_POOL_SIZE = 500;

    bool contaminants_loaded_;
  };

}

#endif
