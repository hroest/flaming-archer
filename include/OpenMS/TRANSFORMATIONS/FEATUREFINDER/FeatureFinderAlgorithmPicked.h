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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>

#include <numeric>

namespace OpenMS
{
	/** 
		@brief FeatureFinderAlgorithm for picked peaks.
		
		@improvment Estimate intensity cutoff from histogram of each bin (Marc)
		
		@ingroup FeatureFinder
	*/
	template<class PeakType, class FeatureType> class FeatureFinderAlgorithmPicked 
		: public FeatureFinderAlgorithm<PeakType, FeatureType>,
			public FeatureFinderDefs
	{
		public:
			///@name Type definitions
			//@{
			typedef typename FeatureFinderAlgorithm<PeakType, FeatureType>::MapType MapType;
			typedef typename MapType::SpectrumType SpectrumType;
			typedef typename PeakType::CoordinateType CoordinateType;
			typedef typename PeakType::IntensityType IntensityType;
			typedef MSExperiment<Peak1D> FilteredMapType;
			//@}
			
			using FeatureFinderAlgorithm<PeakType, FeatureType>::map_;
			using FeatureFinderAlgorithm<PeakType, FeatureType>::param_;
				
		protected:
			///Helper structure for seeds used in FeatureFinderAlgorithmPicked
			struct Seed
			{
				///Spectrum index
				UInt spectrum;
				///Peak index
				UInt peak;
				///Intensity
				Real intensity;

				///Map from isotope pattern charge to isotope pattern score
				std::map<UInt, DoubleReal> charges;
				
				/// Comparison operator
				bool operator<(const Seed& rhs) const
				{
					return intensity<rhs.intensity;
				}
			};
			
			///Helper struct for mass traces used in FeatureFinderAlgorithmPicked
			struct MassTrace
			{
				///Maximum peak pointer
				const FilteredMapType::PeakType* max_peak;
				///RT of maximum peak
				DoubleReal max_rt;
				
				///Contained peaks (pair of RT and pointer to peak)
				std::vector<std::pair<DoubleReal, const FilteredMapType::PeakType*> > peaks;
				
				///determindes the convex hull of the trace
				ConvexHull2D getConvexhull() const
				{
					ConvexHull2D::PointArrayType hull_points(peaks.size());
					for (UInt i=0; i<peaks.size(); ++i)
					{
						hull_points[i][0] = peaks[i].first;
						hull_points[i][1] = peaks[i].second->getMZ();
					}
					return hull_points;
				}
				
				///Sets the maximum to the highest contained peak of the trace
				void updateMaximum()
				{
					if (peaks.size()==0) return;

					max_rt = peaks.begin()->first;					
					max_peak = peaks.begin()->second;
					
					for (UInt i=1; i<peaks.size(); ++i)
					{
						if (peaks[i].second->getIntensity()>max_peak->getIntensity())
						{
							max_rt = peaks[i].first;					
							max_peak = peaks[i].second;
						}
					}
				}
			};

			///Helper structure for a theoretical isotope pattern
			struct TheoreticalIsotopePattern
			{
				///Vector of intensity contributions 
				std::vector<DoubleReal> intensity;
				///Number of optional peaks at the beginning of the pattern
				UInt optional_begin;
				///Number of optional peaks at the end of the pattern
				UInt optional_end;

				/// Returns the size
				UInt size() const
				{
					return intensity.size();
				}
			};

			///Helper structure for a found isotope pattern
			struct IsotopePattern
			{
				///Peak index (-1 if peak was not found, -2 if it was removed to improve the isotope fit)
				std::vector<Int> peak;
				///Spectrum index (undefined if peak index is -1 or -2)
				std::vector<UInt> spectrum;
				///Peak intensity (0 if peak index is -1 or -2)
				std::vector<DoubleReal> intensity;
				///m/z score of peak (0 if peak index is -1 or -2)
				std::vector<DoubleReal> mz_score;
				///Theoretical m/z value of the isotope peak
				std::vector<DoubleReal> theoretical_mz;
				
				/// Constructor that resizes the internal vectors
				IsotopePattern(UInt size)
					: peak(size),
						spectrum(size),
						intensity(size),
						mz_score(size),
						theoretical_mz(size)
				{
				}
			};
			

		public:			
			/// default constructor 
			FeatureFinderAlgorithmPicked() 
				: FeatureFinderAlgorithm<PeakType,FeatureType>(),
					log_("featurefinder.log")
			{
				//debugging
				this->defaults_.setValue("debug",0,"If not 0 debug mode is activated. Then several files with intermediate results are written.");
				//intensity
				this->defaults_.setValue("intensity:bins",10,"Number of bins per dimension (RT and m/z).");
				this->defaults_.setValue("intensity:percentage",35.0,"Percentage of most intense peaks per bin that might be part of a feature.");
				this->defaults_.setSectionDescription("intensity","Settings for the calculation of a score indicating if a peak's intensity is significant (between 0 and 1)");
				//mass trace search parameters
				this->defaults_.setValue("mass_trace:mz_tolerance",0.06,"m/z difference tolerance of peaks belonging to the same mass trace.");
				this->defaults_.setValue("mass_trace:min_spectra",14,"Number of spectra the have to show the same peak mass for a mass trace.");
				this->defaults_.setValue("mass_trace:max_missing",4,"Number of spectra where a high mass deviation or missing peak is acceptable.");
				this->defaults_.setValue("mass_trace:slope_bound",0.1,"The maximum slope of mass trace intensities when extending from the highest peak", true);
				this->defaults_.setSectionDescription("mass_trace","Settings for the calculation of a score indicating if a peak is part of a mass trace (between 0 and 1).");
				//Isotopic pattern search paramters
				this->defaults_.setValue("isotopic_pattern:charge_low",1,"Lowest charge to search for.");
				this->defaults_.setValue("isotopic_pattern:charge_high",4,"Highest charge to search for.");
				this->defaults_.setValue("isotopic_pattern:mz_tolerance",0.06,"Tolerated mass deviation from the theoretical isotopic pattern.");		
				this->defaults_.setValue("isotopic_pattern:intensity_percentage",10.0,"Isotopic peaks that contribute more than this percentage to the overall isotope pattern intensity must be present.", true);
				this->defaults_.setValue("isotopic_pattern:intensity_percentage_optional",0.1,"Isotopic peaks that contribute more than this percentage to the overall isotope pattern intensity can be missing.", true);
				this->defaults_.setValue("isotopic_pattern:optional_fit_improvment",3.0,"Minimal percental improvement of isotope fit to allow leaving out an optional peak.", true);
				this->defaults_.setValue("isotopic_pattern:mass_window_width",100.0,"Window width in Dalton for precalcuation of estimated isotope distribtions.", true);
				this->defaults_.setSectionDescription("isotopic_pattern","Settings for the calculation of a score indicating if a peak is part of a isotoipic pattern (between 0 and 1).");
				//Feature settings
				this->defaults_.setValue("feature:intensity_as_max","true","Determines if feature intensity is reported as the maximum of the feature peaks (true) or the sum of all intensities (false).");
				this->defaults_.setValue("feature:minimum_quality",0.75, "Overall quality threshold for a feature to be reported.");
				this->defaults_.setValue("feature:min_isotope_fit",0.65,"Minimum isotope fit quality.", true);
				this->defaults_.setValue("feature:mass_trace_max_border_intensity",0.7, "Factor how much intensity the border peaks of a mass trace are allowed to have in comarison to the maximum.", true);
				this->defaults_.setSectionDescription("feature","Settings for the features (intensity, quality assessment, ...)");
				
				this->defaultsToParam_();
			}
			
			//TODO try without high_score_map_
			
			/// Main method for actual FeatureFinder
			virtual void run()
			{
				//-------------------------------------------------------------------------
				//Initialization
				//-------------------------------------------------------------------------				
				std::vector< Seed > seeds;
				high_score_map_.resize(map_->size());				
				this->features_->reserve(1000);
				//mass trace search
				UInt min_spectra = std::floor((DoubleReal)param_.getValue("mass_trace:min_spectra")*0.5);
				//isotope pattern search
				UInt charge_low = param_.getValue("isotopic_pattern:charge_low");
				UInt charge_high = param_.getValue("isotopic_pattern:charge_high");
				//quality estimation
				DoubleReal mass_trace_max_border_intensity = param_.getValue("feature:mass_trace_max_border_intensity");
				DoubleReal min_feature_quality = param_.getValue("feature:minimum_quality");
				//feature intensity
				bool max_intensity = String(param_.getValue("feature:intensity_as_max"))=="true";
				
				//---------------------------------------------------------------------------
				//Step 1:
				//For each bin find the intensity threshold of the x% most intense peaks
				//---------------------------------------------------------------------------
				log_ << "Precalculating intensity thresholds ..." << std::endl;
				//new scope to make variables disapear
				{
					DoubleReal percentage = param_.getValue("intensity:percentage");
					DoubleReal rt_start = map_->getMinRT();
					DoubleReal mz_start = map_->getMinMZ();
					intensity_rt_step_ = (map_->getMaxRT() - rt_start ) / (DoubleReal)intensity_bins_;
				 	intensity_mz_step_ = (map_->getMaxMZ() - mz_start ) / (DoubleReal)intensity_bins_;
					intensity_thresholds_.resize(intensity_bins_);
					for (UInt rt=0; rt<intensity_bins_; ++rt)
					{
						intensity_thresholds_[rt].resize(intensity_bins_);
						DoubleReal min_rt = rt_start + rt * intensity_rt_step_;
						DoubleReal max_rt = rt_start + ( rt + 1 ) * intensity_rt_step_;
						std::vector<DoubleReal> tmp;
						for (UInt mz=0; mz<intensity_bins_; ++mz)
						{
							DoubleReal min_mz = mz_start + mz * intensity_mz_step_;
							DoubleReal max_mz = mz_start + ( mz + 1 ) * intensity_mz_step_;
							//std::cout << "rt range: " << min_rt << " - " << max_rt << std::endl;
							//std::cout << "mz range: " << min_mz << " - " << max_mz << std::endl;
							tmp.clear();
							for (typename MapType::ConstAreaIterator it = map_->areaBeginConst(min_rt,max_rt,min_mz,max_mz); it!=map_->areaEndConst(); ++it)
							{
								tmp.push_back(it->getIntensity());
							}
							std::sort(tmp.begin(), tmp.end());
							//std::cout << "number of peaks: " << tmp.size() << std::endl;
							if (tmp.size()==0)
							{
								intensity_thresholds_[rt][mz] = std::make_pair(0.0,0.0);
							}
							else
							{
								//std::cout << "rt:" << rt << " mz:" << mz << " index="<< index << std::endl;
								//std::cout << (min_rt+max_rt)/2.0 << " " << (min_mz+max_mz)/2.0 << " " << tmp.size() << std::endl;
								UInt index = (UInt)std::ceil(tmp.size()*(100.0-percentage)/100.0);
								intensity_thresholds_[rt][mz] = std::make_pair(tmp[index], tmp.back());
							}
						}
					}
				}

				//-------------------------------------------------------------------------
				//debugging
				bool debug = ((UInt)(param_.getValue("debug"))!=0);				
				MapType int_map;
				MapType trace_map;
				MapType pattern_map;
				MapType charge_map;
				FeatureMap<> seed_map;

				if (debug)
				{
					int_map.resize(map_->size());
					trace_map.resize(map_->size());
					pattern_map.resize(map_->size());
					charge_map.resize(map_->size());
					
					CoordinateType rt = (*map_)[0].getRT();
					for (UInt s=0; s<min_spectra; ++s)
					{
						int_map[s].setRT(rt);
						trace_map[s].setRT(rt);
						pattern_map[s].setRT(rt);
						charge_map[s].setRT(rt);
						int_map[map_->size()-1-s].setRT(rt);
						trace_map[map_->size()-1-s].setRT(rt);
						pattern_map[map_->size()-1-s].setRT(rt);
						charge_map[map_->size()-1-s].setRT(rt);
					}
				}
				
				//MAIN LOOP FOR SPECTRA
				this->ff_->startProgress(0, map_->size(), "Finding seeds");
				for (UInt s=0; s<map_->size(); ++s)
				{
					this->ff_->setProgress(s);
					const SpectrumType& spectrum = (*map_)[s];
					high_score_map_[s].setRT(spectrum.getRT());
					
					//do nothing for the first few and last few spectra as
					//the scans required to search for traces are missing
					if (s<min_spectra || s>=map_->size()-min_spectra)
					{
						continue;
					}
					
					//-----------------------------------------------------------
					//Step 2: Calculate IsotopePattern score for all peaks
					// - Test all charges using m/z and intensity fit
					//-----------------------------------------------------------
					//Intermediate storage of all charge scores
					std::vector<std::map<UInt,DoubleReal> > all_pattern_scores(spectrum.size());
					for (UInt p=0; p<spectrum.size(); ++p)
					{
						bool debug_local = false;
						//if (std::fabs(spectrum.getRT()-4348.5)<0.1 && std::fabs(spectrum[p].getMZ()-434.9524)<0.1) debug_local=true;
						if (debug_local) log_ << "Debug isotope score: " << spectrum.getRT() << "/" << spectrum[p].getMZ() << std::endl;
						for (UInt c=charge_low; c<=charge_high; ++c)
						{
							if (debug_local) log_ << "- charge: " << c << std::endl;
							DoubleReal mz = spectrum[p].getMZ();
							//get isotope distribution for this mass
							const TheoreticalIsotopePattern& isotopes = getIsotopeDistribution_(mz*c);
							//determine highest peak in isopope distribution
							UInt max_isotope = std::max_element(isotopes.intensity.begin(), isotopes.intensity.end()) - isotopes.intensity.begin();
							//Look up expected isotopic peaks
							Int peak_index = spectrum.findNearest(mz-((DoubleReal)(isotopes.size()+1)/c));
							IsotopePattern pattern(isotopes.size());
							for (UInt i=0; i<isotopes.size(); ++i)
							{
								if (debug_local) log_ << "  - isotope " << i << std::endl;
								DoubleReal isotope_pos = mz + ((DoubleReal)i-max_isotope)/c;
								if (debug_local) log_ << "    - looking for mass: " << isotope_pos << std::endl;
								peak_index = nearest_(isotope_pos, spectrum, peak_index);
								if (debug_local) log_ << "    - nearest peak mass: " << spectrum[peak_index].getMZ() << std::endl;
								DoubleReal mz_score = positionScore_(isotope_pos, spectrum[peak_index].getMZ(), pattern_tolerance_);
								if (mz_score>0.0) //found
								{
									if (debug_local) log_ << "    - found at " << spectrum[peak_index].getMZ() << std::endl;
									pattern.peak[i] = peak_index;
									pattern.mz_score[i] = mz_score;
									pattern.intensity[i] = spectrum[peak_index].getIntensity();
									pattern.theoretical_mz[i] = isotope_pos;
								}
								else //not found
								{
									if (debug_local) log_ << "    - missing " << std::endl;
									pattern.peak[i] = -1;
									pattern.mz_score[i] = 0.0;
									pattern.intensity[i] = 0.0;
									pattern.theoretical_mz[i] = isotope_pos;
								}
							}
							DoubleReal pattern_score = isotopeScore_(isotopes, pattern, true, false);
							if (debug_local) log_ << "  - final score: " << pattern_score << std::endl;
							
							//update pattern scores and charges of all contained peaks (if necessary)
							if (pattern_score > 0.0)
							{
								for (std::vector<Int>::const_iterator it=pattern.peak.begin(); it!=pattern.peak.end();++it)
								{
									if (*it>=0)
									{
										if (pattern_score>all_pattern_scores[*it][c])
										{
											all_pattern_scores[*it][c] = pattern_score;
										}
									}
								}
							}
						}
					}
					//-----------------------------------------------------------
					//Step 3:
					//Find the mass traces
					//-----------------------------------------------------------					
					if (debug)
					{
						int_map[s].setRT(spectrum.getRT());
						trace_map[s].setRT(spectrum.getRT());
						pattern_map[s].setRT(spectrum.getRT());
						charge_map[s].setRT(spectrum.getRT());
					}
					std::vector<UInt> indices_after(min_spectra+1, 0);
					std::vector<UInt> indices_before(min_spectra+1, 0);
					while( indices_after[0] < spectrum.size() )
					{
						DoubleReal trace_score = 0.0;
						
						//--------------------------------------------------------------
						//Calculate the distances to the nearest peaks in adjacent spectra
						std::vector<DoubleReal> scores;
						
						CoordinateType pos = spectrum[indices_after[0]].getMZ();
						IntensityType inte = spectrum[indices_after[0]].getIntensity();
						
						//log_ << std::endl << "Peak: " << pos << std::endl;
						bool is_max_peak = true; //checking the maximum intensity peaks -> use them later as feature seeds.
						for (UInt i=1; i<=min_spectra; ++i)
						{
							const SpectrumType& spec = (*map_)[s+i];
							indices_after[i] = nearest_(pos, spec, indices_after[i]);
							DoubleReal position_score = positionScore_( pos, spec[indices_after[i]].getMZ(), trace_tolerance_);
							if (position_score >0 && spec[indices_after[i]].getIntensity()>inte) is_max_peak = false;
							scores.push_back(position_score);
						}
						indices_before[0] = indices_after[0];
						for (UInt i=1; i<=min_spectra; ++i)
						{
							const SpectrumType& spec = (*map_)[s-i];
							indices_before[i] = nearest_(pos, spec, indices_before[i]);
							DoubleReal position_score = positionScore_( pos, spec[indices_before[i]].getMZ(), trace_tolerance_);
							if (position_score>0 && spec[indices_before[i]].getIntensity()>inte) is_max_peak = false;
							scores.push_back(position_score);
						}
						
						//--------------------------------------------------------------
						//Calculate a consensus score out of the scores calculated before
						std::sort(scores.begin(), scores.end());
						for(UInt i=max_missing_trace_peaks_; i<scores.size(); ++i)
						{
							trace_score += scores[i];
						}
						trace_score /= (2*min_spectra-max_missing_trace_peaks_);

						//------------------------------------------------------------------
						//Look up best isototope pattern charge and score for this peak
						DoubleReal pattern_charge = 0.0;
						DoubleReal pattern_score = 0.0;
						for (std::map<UInt,DoubleReal>::const_iterator it=all_pattern_scores[indices_after[0]].begin(); it!=all_pattern_scores[indices_after[0]].end(); ++it)
						{
							if (it->second > pattern_score)
							{
								pattern_charge = it->first;
								pattern_score = it->second;
							}
						}
						
						//------------------------------------------------------------------
						//Calculate intensity score
						DoubleReal intensity_score = intensityScore_(inte,s,indices_after[0]);

						//------------------------------------------------------------------
						//Calculate final score (geometric mean). Determine seeds
						DoubleReal final_score = std::pow(intensity_score*trace_score*pattern_score, 1.0/3.0);
						if (final_score!=0.0)
						{
							//peak index
							UInt p = indices_after[0];
							
							FilteredMapType::PeakType feature_peak;
							feature_peak.setIntensity(spectrum[p].getIntensity());
							feature_peak.setMZ(spectrum[p].getMZ());
							high_score_map_[s].push_back(feature_peak);
							//local maximum peaks are considered seeds (if they reach a minimum score and isotope fit)
							if (is_max_peak && final_score>0.1 && pattern_score>=0.5*min_isotope_fit_)
							{
								Seed seed;
								seed.spectrum = s;
								seed.peak = high_score_map_[s].size()-1;
								seed.intensity = inte;
								seed.charges = all_pattern_scores[indices_after[0]];
																
								seeds.push_back(seed);
							}
						}
						
						//-------------------------------------------------------------------------
						//debug output
						if (debug)
						{
							PeakType tmp;
							tmp.setPos(pos);
							
							if (trace_score>0.0)
							{
								tmp.setIntensity(trace_score);
								trace_map[s].push_back(tmp);
							}
							if (pattern_score>0.0)
							{
								tmp.setIntensity(pattern_score);
								pattern_map[s].push_back(tmp);
							}
							if (intensity_score>0.0)
							{
								tmp.setIntensity(intensity_score);
								int_map[s].push_back(tmp);
							}
							tmp.setIntensity(pattern_charge);
							charge_map[s].push_back(tmp);
						}
						++indices_after[0];
					}//for peaks
				}
				this->ff_->endProgress();
				std::cout << "Found " << seeds.size() << " seeds." << std::endl;
				
				//------------------------------------------------------------------
				//Step 4:
				//Extension from a seed
				//------------------------------------------------------------------
				this->ff_->startProgress(0,seeds.size(), "Extending seeds");
				std::sort(seeds.rbegin(),seeds.rend());
				if (debug)
				{
					for (UInt i=0; i<seeds.size(); ++i)
					{
						UInt spectrum = seeds[i].spectrum;
						UInt peak = seeds[i].peak;
						Feature tmp;
						tmp.setIntensity(seeds[i].intensity);
						tmp.setRT(high_score_map_[spectrum].getRT());
						tmp.setMZ(high_score_map_[spectrum][peak].getMZ());
						for (std::map<UInt,DoubleReal>::const_iterator it = seeds[i].charges.begin(); it!=seeds[i].charges.end(); ++it)
						{
							tmp.setMetaValue(String("charge_")+it->first, it->second);
						}
						seed_map.push_back(tmp);
					}
				}
				for (UInt i=0; i<seeds.size(); ++i)
				{
					this->ff_->setProgress(i);
					log_ << std::endl << "Seed " << i+1 << " ~ " << 100.0*i/seeds.size()<< "%" << std::endl;
					//If the intensity is zero this seed is already uses in another feature
					const FilteredMapType::SpectrumType& spectrum = high_score_map_[seeds[i].spectrum];
					const FilteredMapType::PeakType peak = spectrum[seeds[i].peak];
					log_ << " - Int: " << peak.getIntensity() << std::endl;
					log_ << " - RT: " << spectrum.getRT() << std::endl;
					log_ << " - MZ: " << peak.getMZ() << std::endl;
					if (seeds[i].intensity == 0.0)
					{
						abort_("Seed was already used");
						continue;
					}
					//----------------------------------------------------------------
					//determine charge(s) of the seed
					std::vector<UInt> charges;
					UInt max_charge = 0;
					UInt second_charge = 0;
					DoubleReal max_score = 0.0;
					DoubleReal second_score = 0.0;
					
					for (std::map<UInt,DoubleReal>::const_iterator it=seeds[i].charges.begin(); it!=seeds[i].charges.end(); ++it)
					{
						if (it->second>max_score)
						{
							second_charge = max_charge;
							second_score = max_score;
							max_score = it->second;
							max_charge = it->first;
						}
					}
					charges.push_back(max_charge);
					//add second best charge if it is above a score threshold and if it is not too bad compared to the first charge
					if (second_score > 0.5*min_isotope_fit_ && (second_score/max_score)>=0.5)
					{
						charges.push_back(second_charge);
					}
					
					for (UInt j=0;j<charges.size(); ++j)
					{
						UInt charge = charges[j];
						log_ << " - charge: " << charge << " (" << j << "/" << charges.size() << ")" << std::endl;
						//----------------------------------------------------------------
						//Find best fitting isotope pattern for this charge (using averagene)
						IsotopePattern best_pattern(0);
						DoubleReal isotope_fit_quality = findBestIsotopeFit_(seeds[i], charge, best_pattern);
						if (isotope_fit_quality<min_isotope_fit_)
						{
							abort_("Isotope pattern score too low");
							continue;
						}
						
						//extend the convex hull in m/z dimension (starting from the trace peaks)
						//missing traces (index is -1) and removed traces (index is -2) are simply skipped
						log_ << "Collecting mass traces" << std::endl;
						std::vector<MassTrace> traces;
						traces.reserve(best_pattern.peak.size());
						for (UInt p=0; p<best_pattern.peak.size(); ++p)
						{
							log_ << " - Trace " << p << std::endl;
							Seed starting_peak;
							starting_peak.spectrum = best_pattern.spectrum[p];
							starting_peak.peak = best_pattern.peak[p];
							if (best_pattern.peak[p]==-2)
							{
								log_ << "   - removed during isotope fit" << std::endl;
								continue;
							}
							else if (best_pattern.peak[p]==-1)
							{
								log_ << "   - missing" << std::endl;
								continue;
							}
							starting_peak.intensity = spectrum[starting_peak.peak].getIntensity();
							log_ << "   - extending from " << high_score_map_[starting_peak.spectrum].getRT() << " / " << high_score_map_[starting_peak.spectrum][starting_peak.peak].getMZ() << std::endl;
							
							//TODO: Really search for nearby maxima, not only for seeds
							//search for nearby maximum of the mass trace
							//as the extension assumes that it starts at the maximum
							for(UInt s=i+1; s<seeds.size(); ++s)
							{
								//the scan is nearby
								if (std::abs((Int)seeds[s].spectrum-(Int)starting_peak.spectrum)<(Int)min_spectra)
								{
									//the peak mass fits
									if (std::fabs(high_score_map_[starting_peak.spectrum][starting_peak.peak].getMZ()-high_score_map_[seeds[s].spectrum][seeds[s].peak].getMZ())<pattern_tolerance_)
									{
										starting_peak.spectrum = seeds[s].spectrum;
										starting_peak.peak = seeds[s].peak;
										log_ << "   - found nearby seed to extend from at " << high_score_map_[starting_peak.spectrum].getRT() << " / " << high_score_map_[starting_peak.spectrum][starting_peak.peak].getMZ() << std::endl;
										break;
									}
								}
							}
							
							//------------------------------------------------------------------
							//Extend seed to a mass trace
							MassTrace trace;
							const FilteredMapType::PeakType& seed = high_score_map_[starting_peak.spectrum][starting_peak.peak];
							//initialize trace with seed data
							trace.peaks.push_back(std::make_pair(high_score_map_[starting_peak.spectrum].getRT(), &seed));
							trace.max_peak = &seed;
							trace.max_rt = high_score_map_[starting_peak.spectrum].getRT();
							//extend in downstream direction
							extendMassTrace_(trace, starting_peak.spectrum -1, seed.getMZ(), false);
							//invert peak array to bring peaks in the correct cronological order
							std::reverse(trace.peaks.begin(), trace.peaks.end());
							//extend in upstream direction
							extendMassTrace_(trace, starting_peak.spectrum +1, seed.getMZ(), true);
							
							//check if enough peaks were found
							if (trace.peaks.size()<2)
							{
								log_ << "   - could not extend trace " << std::endl;
								continue;
							}
							traces.push_back(trace);
						}
						if (traces.size()<2)
						{
							abort_("Found less than two mass traces");
							continue;
						}
						
						//------------------------------------------------------------------
						//Step 5:
						//Quality estimation
						//------------------------------------------------------------------
						log_ << "Quality estimation" << std::endl;
						
						//------------------------------------------------------------------
						//(1) isotope fit: isotope fit of the collected mass trace maxima
						log_ << " - Isotope fit: " << isotope_fit_quality << std::endl;
						
						//------------------------------------------------------------------
						//(2) overall shape: RT-spread of mass traces decreases with smaller intensities
						std::vector<DoubleReal> rts(traces.size());
						std::vector<DoubleReal> ints(traces.size());
						for (UInt j=0; j<traces.size(); ++j)
						{
							rts[j] = traces[j].peaks.back().first - traces[j].peaks[0].first;
							ints[j] = traces[j].max_peak->getIntensity();
						}
						DoubleReal overall_shape_quality = (Math::BasicStatistics<DoubleReal>::pearsonCorrelationCoefficient(rts.begin(),rts.end(),ints.begin(), ints.end())+1.0)/2.0;
						if (isnan(overall_shape_quality))
						{
							if (traces.size()==2) //for two traces it's ok to have the same width/intensity
							{
								overall_shape_quality = 0.5;
							}
							else //for more than two traces it is not ok to have the same width/intensity
							{
								overall_shape_quality = 0.1;
							}
						}
						log_ << " - overall shape: " << overall_shape_quality << std::endl;
		
						//------------------------------------------------------------------					
						//(3) trace m/z distances
						std::vector<DoubleReal> positions(traces.size());
						for (UInt j=0; j<traces.size(); ++j)
						{
							for (UInt k=0; k<traces[j].peaks.size(); ++k)
							{
								positions[j]+= traces[j].peaks[k].second->getMZ();
							}
							positions[j] /= traces[j].peaks.size();
						}
						DoubleReal mz_distance_quality = 0.0;
						for (UInt j=0; j<positions.size()-1; ++j)
						{
							mz_distance_quality += positionScore_(positions[j+1]-positions[j], 1.0/charge, pattern_tolerance_);
						}
						mz_distance_quality /= positions.size()-1;
						log_ << " - mz distances: " << mz_distance_quality << std::endl;
	
						//------------------------------------------------------------------
						//(4) trace shape: trace intensity goes down towards the border
						UInt error_count = 0;
						for (UInt j=0; j<traces.size(); ++j)
						{
							UInt size = traces[j].peaks.size();
							if (size>=5)
							{
								DoubleReal max = traces[j].max_peak->getIntensity();
								//log_ << "------- max: " << max << std::endl;
								DoubleReal low_int = (traces[j].peaks[0].second->getIntensity() + traces[j].peaks[1].second->getIntensity()) / 2.0;
								//log_ << "------- low_int: " << low_int << std::endl;
								if (low_int / max > mass_trace_max_border_intensity)
								{
									//log_ << "------- error " << std::endl;
									++error_count;
								}
								DoubleReal high_int = (traces[j].peaks[size-2].second->getIntensity() + traces[j].peaks[size-1].second->getIntensity()) / 2.0;
								//log_ << "------- high_int: " << high_int << std::endl;
								if (high_int / max > mass_trace_max_border_intensity)
								{
									//log_ << "------- error " << std::endl;
									++error_count;
								}
							}
							else //increase error count by one as this trace is too small
							{
								++error_count;
							}
						}
						//Score: fraction of mal-formed rt profiles
						DoubleReal rt_shape_quality = 1.0 - (DoubleReal)error_count / (2.0*traces.size());					
						log_ << " - trace shape: " << rt_shape_quality << std::endl;
						
						//------------------------------------------------------------------					
						//(5) quality measure: maxima on one line
						//determine max peak RT and RT spread of that trace in both directions
						DoubleReal max = 0.0;
						DoubleReal max_rt = 0.0;
						DoubleReal spread_low = 0.0;
						DoubleReal spread_high = 0.0;
						for (UInt j=0; j<traces.size(); ++j)
						{
							if (traces[j].max_peak->getIntensity() > max)
							{
								max = traces[j].max_peak->getIntensity();
								max_rt = traces[j].max_rt;
								spread_low = std::max(0.01, max_rt - traces[j].peaks[0].first);
								spread_high = std::max(0.01, traces[j].peaks.back().first - max_rt);
							}
						}
						
						//look at max peak shifts of different scans
						DoubleReal rel_max_deviation = 0.0;
						for (UInt j=0; j<traces.size(); ++j)
						{
							if (traces[j].max_rt>max_rt)
							{
								rel_max_deviation += std::min(1.0,(traces[j].max_rt - max_rt) / spread_high);
							}
							else
							{
								rel_max_deviation += std::min(1.0,(max_rt - traces[j].max_rt) / spread_low);
							}
						}
						rel_max_deviation /= traces.size()-1;
						DoubleReal maxima_quality = 1.0 - rel_max_deviation;
						log_ << " - maxima positions: " << maxima_quality << std::endl;
	
						//----------------------------------------------------------------
						//abort if quality too low
						DoubleReal overall_quality_mean = std::pow(isotope_fit_quality * overall_shape_quality * mz_distance_quality * rt_shape_quality * maxima_quality, 1.0/5.0);
						log_ << " => final score: " << overall_quality_mean << std::endl;
						if (overall_quality_mean<min_feature_quality)
						{
							abort_("Feature quality too low");
							continue;
						}
	
						//------------------------------------------------------------------
						//Step 6:
						//Feature creation
						//------------------------------------------------------------------
						Feature f;
						f.setCharge(charge);
						f.setOverallQuality(overall_quality_mean);
						if (debug)
						{
							f.setMetaValue("rt_shape",rt_shape_quality);
							f.setMetaValue("mz_distance",mz_distance_quality);
							f.setMetaValue("isotope_fit",isotope_fit_quality);
							f.setMetaValue("overall_shape",overall_shape_quality);
							f.setMetaValue("maxima_positions",maxima_quality);
						}
						
						//set feature position and intensity
						//feature RT is the average RT of the mass trace maxima
						//feature intensity and m/z are taken from the highest peak
						DoubleReal rt = 0.0;
						for (UInt j=0; j<traces.size(); ++j)
						{
							if (traces[j].max_peak->getIntensity()>f.getIntensity())
							{
								f.setIntensity(traces[j].max_peak->getIntensity());
								f.setMZ(traces[j].max_peak->getMZ());
							}
							rt += traces[j].max_rt;
						}
						f.setRT(rt/traces.size());
						//feature intensity is the sum of all the peak intensities
						if (!max_intensity)
						{
							DoubleReal int_sum = 0.0;
							for (UInt j=0; j<traces.size(); ++j)
							{
								for (UInt k=0; k<traces[j].peaks.size(); ++k)
								{
									int_sum += traces[j].peaks[k].second->getIntensity();
								}
							}
							f.setIntensity(int_sum);
						} 
						
						//add convex hulls of mass traces
						for (UInt j=0; j<traces.size(); ++j)
						{
							f.getConvexHulls().push_back(traces[j].getConvexhull());
						}
						//add feature to feature list
						this->features_->push_back(f);
						log_ << "Feature number: " << this->features_->size() << std::endl;
	
	
						//----------------------------------------------------------------
						//Remove all seeds that lie inside the convex hull of the new feature
						DBoundingBox<2> bb = f.getConvexHull().getBoundingBox();
						for (UInt j=i+1; j<seeds.size(); ++j)
						{
							DoubleReal rt = high_score_map_[seeds[j].spectrum].getRT();
							DoubleReal mz = high_score_map_[seeds[j].spectrum][seeds[j].peak].getMZ();
							if (bb.encloses(rt,mz) && f.encloses(rt,mz))
							{
								//set intensity to zero => the peak will be skipped!
								seeds[j].intensity = 0.0;
							}
						}
					}
				}
				this->ff_->endProgress();
				std::cout << "Found " << this->features_->size() << " features candidates." << std::endl;
				std::cout << std::endl;
				std::cout << "Abort reasons during feature construction:" << std::endl;
				for (std::map<String,UInt>::const_iterator it=aborts_.begin(); it!=aborts_.end(); ++it)
				{
					std::cout << "- " << it->first << ": " << it->second << std::endl;
				}
				
				//------------------------------------------------------------------
				//Step 7:
				//TODO: Resolve contradicting and overlapping features (of same charge?)
				//------------------------------------------------------------------
				
				
				
				//------------------------------------------------------------------
				//store debug info
				if (debug)
				{
					MzDataFile file;
					file.store("trace_scores.mzData",trace_map);
					file.store("charge_estimates.mzData",charge_map);
					file.store("intensity_scores.mzData",int_map);
					file.store("pattern_scores.mzData",pattern_map);
					file.store("selected_peaks.mzData",high_score_map_);
					
					FeatureXMLFile file2;
					file2.store("seeds.featureXML", seed_map);
				}
			}
			
			static FeatureFinderAlgorithm<PeakType,FeatureType>* create()
			{
				return new FeatureFinderAlgorithmPicked();
			}

			static const String getProductName()
			{
				return "picked_peak";
			}
	
		protected:
			/// The map of possible feature peaks
			FilteredMapType high_score_map_;
			/// Output stream for log/debug info
			std::ofstream log_; 
			/// Array of abort reasons
			std::map<String, UInt> aborts_;
						
			///@name Members for parameters often needed in methods
			//@{
			DoubleReal pattern_tolerance_; ///< Stores mass_trace:mz_tolerance
			DoubleReal trace_tolerance_; ///< Stores isotopic_pattern:mz_tolerance
			UInt max_missing_trace_peaks_; ///< Stores mass_trace:max_missing
			DoubleReal slope_bound_; ///< Max slope of mass trace intensities
			DoubleReal intensity_percentage_; ///< Isotope pattern intensity contribution of required peaks
			DoubleReal intensity_percentage_optional_; ///< Isotope pattern intensity contribution of optional peaks
			DoubleReal optional_fit_improvment_; ///< Minimal imrovment for leaving out optional isotope
			DoubleReal mass_window_width_; ///< Width of the isotope pattern mass bins
			UInt intensity_bins_; ///< Number of bins (in RT and MZ) for intensity significance estimation
			DoubleReal min_isotope_fit_; ///< Mimimum isotope pattern fit for a feature
			//@}

			///@name Members for intensity significance estimation
			//@{			
			///< RT bin width
			DoubleReal intensity_rt_step_;
			///< m/z bin width
			DoubleReal intensity_mz_step_;
			///< Precalculated threshold and maximum stored for each bin (rt,mz)
			std::vector< std::vector< std::pair<IntensityType,IntensityType> > > intensity_thresholds_;
			//@}
			
			///Vector of precalculated isotope distributions for several mass winows
			std::vector< TheoreticalIsotopePattern > isotope_distributions_;

			//Docu in base class
			virtual void updateMembers_()
			{
				pattern_tolerance_ = param_.getValue("mass_trace:mz_tolerance");
				trace_tolerance_ = param_.getValue("isotopic_pattern:mz_tolerance");
				max_missing_trace_peaks_ = param_.getValue("mass_trace:max_missing");
				slope_bound_ = param_.getValue("mass_trace:slope_bound");
				intensity_percentage_ = (DoubleReal)param_.getValue("isotopic_pattern:intensity_percentage")/100.0;
				intensity_percentage_optional_ = (DoubleReal)param_.getValue("isotopic_pattern:intensity_percentage_optional")/100.0;
				optional_fit_improvment_ = (DoubleReal)param_.getValue("isotopic_pattern:optional_fit_improvment")/100.0;
				mass_window_width_ = param_.getValue("isotopic_pattern:mass_window_width");
				intensity_bins_ =  param_.getValue("intensity:bins");
				min_isotope_fit_ = param_.getValue("feature:min_isotope_fit");
			}
			
			///Writes the abort reason to the log file and counts occurences for each reason
			void abort_(const String& reason)
			{
				log_ << "Abort: " << reason << std::endl;
				aborts_[reason]++;
			}
			
			///Returns the isotope distribution for a certain mass window
			const TheoreticalIsotopePattern& getIsotopeDistribution_(DoubleReal mass)
			{
				//calculate index in the vector
				UInt index = (UInt)std::floor(mass/mass_window_width_);
				
				//enlarge vector if necessary
				if (index>=isotope_distributions_.size())
				{
					isotope_distributions_.resize(index+1);
				}
				
				//calculate distribution if necessary
				if (isotope_distributions_[index].intensity.size()==0)
				{
					//log_ << "Calculating iso dist for mass: " << 0.5*mass_window_width_ + index * mass_window_width_ << std::endl;
					IsotopeDistribution d;
					d.setMaxIsotope(10);
					d.estimateFromPeptideWeight(0.5*mass_window_width_ + index * mass_window_width_);
					d.trimLeft(intensity_percentage_optional_);
					d.trimRight(intensity_percentage_optional_);
					for (IsotopeDistribution::Iterator it=d.begin(); it!=d.end(); ++it)
					{
						isotope_distributions_[index].intensity.push_back(it->second);
						//log_ << " - " << it->second << std::endl;
					}
					//determine the number of optional peaks at the beginning/end
					UInt begin = 0;
					UInt end = 0;
					bool is_begin = true;
					bool is_end = false;
					for (UInt i=0; i<isotope_distributions_[index].intensity.size(); ++i)
					{
						if (isotope_distributions_[index].intensity[i]<intensity_percentage_)
						{
							if (!is_end && !is_begin) is_end = true;
							if (is_begin) ++begin;
							else if (is_end) ++end;
						}
						else if (is_begin)
						{
							is_begin = false;
						}
					}
					isotope_distributions_[index].optional_begin = begin;
					isotope_distributions_[index].optional_end = end;
					//log_ << " - optinal begin/end:" << begin << " / " << end << std::endl;
				}
				//Return distribution
				return isotope_distributions_[index];
			}
						
			/**
				@brief Finds the best fitting position of the isotopic pattern estimate defined by @p center
				
				@param center the maximum peak of the isotope distribution (contains charge as well)
				@param charge The charge of the pattern 
				@param best_pattern Returns the indices of the isotopic peaks. If a isopopic peak is missing -1 is returned.
			*/
			DoubleReal findBestIsotopeFit_(const Seed& center, UInt charge, IsotopePattern& best_pattern)
			{
				log_ << "Testing isotope patterns for charge " << charge << ": " << std::endl;			
				const FilteredMapType::SpectrumType& spectrum = high_score_map_[center.spectrum];
				const TheoreticalIsotopePattern& isotopes = getIsotopeDistribution_(spectrum[center.peak].getMZ()*charge);	
				log_ << " - Seed: " << center.peak << " (mz:" << spectrum[center.peak].getMZ()<< ")" << std::endl;
				
				//Find m/z boundaries of search space (linear search as this is local and we have the center already)
				DoubleReal mass_window = (DoubleReal)(isotopes.size()+1) / (DoubleReal)charge;
				log_ << " - Mass window: " << mass_window << std::endl;
				UInt end = center.peak;
				while(end<spectrum.size() && spectrum[end].getMZ()<spectrum[center.peak].getMZ()+mass_window)
				{
					++end;
				}
				--end;
				//search begin
				Int begin = center.peak;
				while(begin>=0 && spectrum[begin].getMZ()>spectrum[center.peak].getMZ()-mass_window)
				{
					--begin;
				}
				++begin;
				log_ << " - Begin: " << begin << " (mz:" << spectrum[begin].getMZ()<< ")" << std::endl;
				log_ << " - End: " << end << " (mz:" << spectrum[end].getMZ()<< ")" << std::endl;

				//fit isotope distribution to peaks
				DoubleReal max_score = 0.0;
				for (UInt start=begin; start<=end; ++start)
				{
					//find isotope peaks for the current start peak
					UInt peak_index = start;
					IsotopePattern pattern(isotopes.size());
					pattern.intensity[0] = spectrum[start].getIntensity();
					pattern.peak[0] = start;
					pattern.spectrum[0] = center.spectrum;
					pattern.mz_score[0] = 1.0;
					pattern.theoretical_mz[0] = spectrum[start].getMZ();
					log_ << " - Fitting at " << start << " (mz:" << spectrum[start].getMZ() << ")" << std::endl;
					log_ << "   - Isotope 0: " << pattern.intensity[0] << " (" << isotopes.intensity[0] << ")" << std::endl;
					for (UInt iso=1; iso<isotopes.size(); ++iso)
					{
						DoubleReal pos = spectrum[start].getMZ() + iso/(DoubleReal)charge;
						peak_index = nearest_(pos, spectrum, peak_index);
						DoubleReal mz_score = positionScore_(pos, spectrum[peak_index].getMZ(), pattern_tolerance_);
						pattern.theoretical_mz[iso] = pos;
						if (mz_score!=0.0) //found
						{
							log_ << "   - Isotope " << iso << ": " << spectrum[peak_index].getIntensity() << " (" << isotopes.intensity[iso] << ")" << std::endl;
							pattern.peak[iso] = peak_index;
							pattern.spectrum[iso] = center.spectrum;
							pattern.mz_score[iso] = mz_score;
							pattern.intensity[iso] = spectrum[peak_index].getIntensity();
						}
						else //missing
						{					
							//try to find the mass in the previous spectrum
							const FilteredMapType::SpectrumType& spectrum_before = high_score_map_[center.spectrum-1];
							Int index_before = spectrum_before.findNearest(pos);
							if (index_before!=-1 && positionScore_(pos, spectrum_before[index_before].getMZ(), pattern_tolerance_)!=0.0)
							{
								log_ << "   - Isotope " << iso << ": " << spectrum_before[index_before].getIntensity() << " (" << isotopes.intensity[iso] << ") - previous spectrum" << std::endl;
								pattern.peak[iso] = index_before;
								pattern.spectrum[iso] = center.spectrum-1;
								pattern.mz_score[iso] = positionScore_(pos, spectrum_before[index_before].getMZ(), pattern_tolerance_);
								pattern.intensity[iso] = spectrum_before[index_before].getIntensity();
							}
							else
							{
								//try to find the mass in the next spectrum
								const FilteredMapType::SpectrumType& spectrum_after = high_score_map_[center.spectrum+1];
								Int index_after = spectrum_after.findNearest(pos);
								if (index_after!=-1 && positionScore_(pos, spectrum_after[index_after].getMZ(), pattern_tolerance_)!=0.0)
								{
									log_ << "   - Isotope " << iso << ": " << spectrum_after[index_after].getIntensity() << " (" << isotopes.intensity[iso] << ") - next spectrum" << std::endl;
									pattern.peak[iso] = index_after;
									pattern.spectrum[iso] = center.spectrum+1;
									pattern.mz_score[iso] = positionScore_(pos, spectrum_after[index_after].getMZ(), pattern_tolerance_);
									pattern.intensity[iso] = spectrum_after[index_after].getIntensity();
								}
								else
								{
									log_ << "   - Isotope " << iso << ": missing" << " (" << isotopes.intensity[iso] << ")" << std::endl;
									pattern.peak[iso] = -1;
									pattern.mz_score[iso] = 0.0;
									pattern.intensity[iso] = 0.0;
								}
							}
						}
					}
					
					//check if the seed is contained, otherwise abort
					bool seed_contained = false;
					for (UInt iso=0; iso<pattern.peak.size(); ++iso)
					{
						if (pattern.peak[iso]==(Int)center.peak && pattern.spectrum[iso]==center.spectrum)
						{
							seed_contained = true;
							break;
						}
					}
					if(!seed_contained)
					{
						log_ << "   - aborting: seed is not contained!" << std::endl;
						continue;
					}

					DoubleReal score = isotopeScore_(isotopes, pattern, false, true);
	
					//check if the seed is still contained, otherwise abort
					seed_contained = false;
					for (UInt iso=0; iso<pattern.peak.size(); ++iso)
					{
						if (pattern.peak[iso]==(Int)center.peak && pattern.spectrum[iso]==center.spectrum)
						{
							seed_contained = true;
							break;
						}
					}
					if(!seed_contained)
					{
						log_ << "   - aborting: seed was removed during isotope fit!" << std::endl;
						continue;
					}
					
					log_ << "   - final score: " << score << std::endl;
					if (score>max_score)
					{
						max_score = score;
						best_pattern = pattern;
					}
				}
				log_ << " - best score: " << max_score << std::endl;
				return max_score;
			}

			//TODO: Find better ways to determine mass trace ends
			//      - Relative change compared to maximum?
			/**
				@brief Extends mass trace from maximum in one RT direction
				
				@note this method assumes that it extends from a local maximum. Otherwise it will not work!
			*/
			void extendMassTrace_(MassTrace& trace, Int spectrum_index, DoubleReal mz, bool inc_rt)
			{
				std::vector<DoubleReal> deltas(3,0.0);
				UInt missing_peaks = 0;
				DoubleReal last_intensity = trace.max_peak->getIntensity();
				UInt added_peaks = 0;
				bool remove_last_peaks = false;
				while((!inc_rt && spectrum_index>=0) || (inc_rt && spectrum_index<(Int)high_score_map_.size()))
				{
					Int peak_index = high_score_map_[spectrum_index].findNearest(mz);
					if (peak_index<0)
					{
						++missing_peaks;
					}
					else
					{
						const FilteredMapType::PeakType& peak = high_score_map_[spectrum_index][peak_index];
						if ( positionScore_( mz, peak.getMZ(), trace_tolerance_)>0 )
						{
							deltas.erase(deltas.begin());
							deltas.push_back((peak.getIntensity() - last_intensity) / last_intensity);
							last_intensity = peak.getIntensity();
							missing_peaks = 0;
							trace.peaks.push_back(std::make_pair(high_score_map_[spectrum_index].getRT(), &peak));
							++added_peaks;
							//update maximum peak of trace
							if (peak.getIntensity()>trace.max_peak->getIntensity())
							{
								trace.max_peak = &peak;
								trace.max_rt = high_score_map_[spectrum_index].getRT();
							}
						}
						else
						{
							++missing_peaks;
						}
					}
					//Abort if too many peaks are missing
					if(missing_peaks>max_missing_trace_peaks_)
					{
						break;
					}
					//Abort if the last three deltas are positive
					if (deltas[0]>=0.0 && deltas[1]>=0.0 && deltas[2]>=0.0)
					{
						remove_last_peaks = true;
						break;
					}
					//Abort if the average delta is too big (as intensity increases then)
					if (std::accumulate(deltas.begin(),deltas.end(),0.0)/3.0>slope_bound_)
					{
						remove_last_peaks = true;
						break;
					}
					//increase/decrease scan index
					if (inc_rt) ++spectrum_index; else --spectrum_index;
				}
				//Go back several peaks if we extended too far
				//Update maximum, if we removed it
				if (remove_last_peaks)
				{
					bool max_removed = false;
					for(UInt i=std::min(added_peaks,3u); i>0; --i)
					{
						if (trace.peaks.back().second==trace.max_peak) max_removed = true;
						trace.peaks.pop_back();
					}
					if (max_removed)
					{
						trace.updateMaximum();
					}
				}
			}

			/// Returns the index of the peak nearest to m/z @p pos in spectrum @p spec (linear search starting from index @p start)
			template <typename SpectrumType>
			UInt nearest_(CoordinateType pos, const SpectrumType& spec, UInt start) const
			{
				UInt index = start;
				CoordinateType dist = std::fabs(pos-spec[index].getMZ());
				++index;
				while (index < spec.size())
				{
					CoordinateType new_dist = std::fabs(pos-spec[index].getMZ());
					if (new_dist<dist)
					{
						dist = new_dist;
						++index;	
					}
					else
					{
						break;
					}
				}
				return --index; 
			}

			/// Calculates a score between 0 and 1 for the m/z deviation of two peaks.
			DoubleReal positionScore_(CoordinateType pos1, CoordinateType pos2, DoubleReal allowed_deviation) const
			{
				DoubleReal diff = fabs(pos1 - pos2);
				if (diff <= 0.5*allowed_deviation)
				{
					return 0.1*(0.5*allowed_deviation-diff)/(0.5*allowed_deviation)+0.9;
				}
				else if (diff <= allowed_deviation)
				{
					return 0.9*(allowed_deviation-diff)/(0.5*allowed_deviation);
				}
				return 0.0;
			}

			/// Calculates a score between 0 and 1 for the correlation between theoretical and found isotope pattern
			DoubleReal isotopeScore_(const TheoreticalIsotopePattern& isotopes, IsotopePattern& pattern, bool consider_mz_distances, bool debug)
			{
				//Abort if a core peak is missing
				for (UInt iso=0+isotopes.optional_begin; iso<pattern.peak.size()-isotopes.optional_end; ++iso)
				{
					if (pattern.peak[iso]==-1)
					{
						if (debug) log_ << "   - aborting: core peak is missing" << std::endl;
						return 0.0;
					}
				}
				//Find best isotope fit
				// - try to leave out optional isotope peaks to improve the fit
				// - do not allow gaps inside the pattern
				DoubleReal best_int_score = 0.1; //Not 0 as this would result in problems when checking for the percental improvement
				UInt best_begin = 0;
				for (UInt i=isotopes.optional_begin; i>0; --i)
				{
					if (pattern.peak[i-1]==-1)
					{
						best_begin = i;
						break;
					}
				}
				UInt best_end = 0;
				for (UInt i=isotopes.optional_end; i>0; --i)
				{
					if (pattern.peak[pattern.peak.size()-i]==-1)
					{
						best_end = i;
						break;
					}
				}
				if (debug) log_ << "   - best_begin/end: " << best_begin << "/" << best_end << std::endl;
				for (UInt b=best_begin; b<=isotopes.optional_begin; ++b)
				{
					for (UInt e=best_end; e<=isotopes.optional_end; ++e)
					{
						//Make sure we have more than 2 peaks (unless in the first loop interation) 
						if (isotopes.size()-b-e>2 || (b==best_begin && e==best_end))
						{
							DoubleReal int_score = Math::BasicStatistics<DoubleReal>::pearsonCorrelationCoefficient(isotopes.intensity.begin()+b, isotopes.intensity.end()-e, pattern.intensity.begin()+b, pattern.intensity.end()-e);	
							if (isnan(int_score)) int_score = 0.0;
							if (isotopes.size()-b-e==2 && int_score>min_isotope_fit_) int_score = min_isotope_fit_; //special case for the first loop iteration (otherwise the score is 1)
							if (debug) log_ << "   - fit (" << b << "/" << e << "): " << int_score;
							if (int_score/best_int_score>=1.0+optional_fit_improvment_)
							{
								if (debug) log_ << " - new best fit ";
								best_int_score = int_score;
								best_begin = b;
								best_end = e;
							}
							if (debug) log_ << std::endl;
						}
					}
				}
				//remove left out peaks from the beginning
				for (UInt i=0; i<best_begin; ++i)
				{
					pattern.peak[i] = -2;
					pattern.intensity[i] = 0.0;
					pattern.mz_score[i] = 0.0;
				}
				//remove left out peaks from the end
				for (UInt i=0; i<best_end; ++i)
				{
					pattern.peak[isotopes.size()-1-i] = -2;
					pattern.intensity[isotopes.size()-1-i] = 0.0;
					pattern.mz_score[isotopes.size()-1-i] = 0.0;
				}
				//calculate m/z score (if required)
				if (consider_mz_distances)
				{
					best_int_score *= std::accumulate(pattern.mz_score.begin()+best_begin, pattern.mz_score.end()-best_end,0.0) / (pattern.mz_score.size()-best_begin-best_end);
				}
				//return final score
				return best_int_score;
			}
			
			DoubleReal intensityScore_(DoubleReal intensity, UInt spectrum, UInt peak)
			{
				UInt rt_bin = std::min(intensity_bins_-1,(UInt)std::floor(((*map_)[spectrum].getRT() - map_->getMinRT()) / intensity_rt_step_));
				UInt mz_bin = std::min(intensity_bins_-1,(UInt)std::floor(((*map_)[spectrum][peak].getMZ() - map_->getMinMZ()) / intensity_mz_step_));
						
				DoubleReal threshold = intensity_thresholds_[rt_bin][mz_bin].first;
				DoubleReal maximum = intensity_thresholds_[rt_bin][mz_bin].second;
				if (intensity>threshold)
				{
					return (intensity-threshold)/(maximum-threshold);
				}
				return 0.0;
			}

		private:
			
			/// Not implemented
			FeatureFinderAlgorithmPicked& operator=(const FeatureFinderAlgorithmPicked&);
			/// Not implemented
			FeatureFinderAlgorithmPicked(const FeatureFinderAlgorithmPicked&);

	};

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_H
