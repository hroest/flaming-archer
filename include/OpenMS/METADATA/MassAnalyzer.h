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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_MASSANALYZER_H
#define OPENMS_METADATA_MASSANALYZER_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/**
		@brief Descripton of a mass analyzer ( Part of a MS Instrument )
		
		
		
		@ingroup Metadata
	*/
  class MassAnalyzer: public MetaInfoInterface
  {
    public:
    	/// analyzer type
	    enum AnalyzerType {ANALYZERNULL,QUADRUPOLE,PAULIONTRAP,RADIALEJECTIONLINEARIONTRAP,AXIALEJECTIONLINEARIONTRAP,TOF,SECTOR,FOURIERTRANSFORM,IONSTORAGE,SIZE_OF_ANALYZERTYPE};
			/// Names of the analyzer types                    
			static const std::string NamesOfAnalyzerType[SIZE_OF_ANALYZERTYPE];                    
			
			/**
				@brief resolution method
				
				Which of the available standard measures is used to define whether two peaks are separate
			*/									 
			enum ResolutionMethod {RESMETHNULL,FWHM,TENPERCENTVALLEY,BASELINE,SIZE_OF_RESOLUTIONMETHOD};
			/// Names of resolustion methods
			static const std::string NamesOfResolutionMethod[SIZE_OF_RESOLUTIONMETHOD];
				
			enum ResolutionType {RESTYPENULL,CONSTANT,PROPORTIONAL,SIZE_OF_RESOLUTIONTYPE};
			/// Names of resulution type
			static const std::string NamesOfResolutionType[SIZE_OF_RESOLUTIONTYPE];
				
			enum ScanFunction {SCANFCTNULL,SELECTEDIONDETECTION,MASSSCAN,SIZE_OF_SCANFUNCTION};
			/// Names of scan functions
			static const std::string NamesOfScanFunction[SIZE_OF_SCANFUNCTION];
				
			/// direction of scanning
			enum ScanDirection {SCANDIRNULL,UP,DOWN,SIZE_OF_SCANDIRECTION};
			/// Names of direction of scanning
			static const std::string NamesOfScanDirection[SIZE_OF_SCANDIRECTION];

			enum ScanLaw	{SCANLAWNULL,EXPONENTIAL,LINEAR,QUADRATIC,SIZE_OF_SCANLAW};
			/// Names of scan laws
			static const std::string NamesOfScanLaw[SIZE_OF_SCANLAW];

			/// MS/MS scan method
			enum TandemScanningMethod {TANDEMNULL,PRODUCTIONSCAN,PRECURSORIONSCAN,CONSTANTNEUTRALLOSS,SINGLEREACTIONMONITORING,MULTIPLEREACTIONMONITORING,SINGLEIONMONITORING,MULTIPLEIONMONITORING,SIZE_OF_TANDEMSCANNINGMETHOD};
			/// Names of MS/MS scan methods
			static const std::string NamesOfTandemScanningMethod[SIZE_OF_TANDEMSCANNINGMETHOD];

			/// reflectron state
			enum ReflectronState {REFLSTATENULL,ON,OFF,NONE,SIZE_OF_REFLECTRONSTATE};
			/// Names of reclectron states
			static const std::string NamesOfReflectronState[SIZE_OF_REFLECTRONSTATE];
			
			/// Constructor
      MassAnalyzer();
      /// Copy constructor
      MassAnalyzer(const MassAnalyzer& source);
      /// Destructor
      ~MassAnalyzer();
			
			/// Assignment operator
      MassAnalyzer& operator= (const MassAnalyzer& source);
 
      /// Equality operator
      bool operator== (const MassAnalyzer& rhs) const;
      /// Equality operator
      bool operator!= (const MassAnalyzer& rhs) const;
			
			/// returns the analyzer type
      AnalyzerType getType() const;
      /// sets the analyzer type
      void setType(AnalyzerType type);
			
			/// returns the method used for determination of the resolution
      ResolutionMethod getResolutionMethod() const;
      /// sets the method used for determination of the resolution
      void setResolutionMethod(ResolutionMethod resolution_method);
			
			/// returns the 
      ResolutionType getResolutionType() const;
      /// sets the 
      void setResolutionType(ResolutionType resolution_type);
			
			/// returns the 
      ScanFunction getScanFunction() const;
      /// sets the 
      void setScanFunction(ScanFunction scan_function);
			
			/// returns the direction of scanning
      ScanDirection getScanDirection() const;
      /// sets the direction of scanning
      void setScanDirection(ScanDirection scan_direction);
			
			/// returns the 
      ScanLaw getScanLaw() const;
      /// sets the 
      void setScanLaw(ScanLaw scan_law);
			
			/// returns the MS/MS scanning method
      TandemScanningMethod getTandemScanMethod() const;
      /// sets the MS/MS scanning method
      void setTandemScanMethod(TandemScanningMethod tandem_scan_method);
			
			/// returns the reflectron state (for TOF)
      ReflectronState getReflectronState() const;
      /// sets the reflectron state (for TOF)
      void setReflectronState(ReflectronState reflecton_state);
			
			/**
				@brief returns the resolution
			
				The maximum m/z value at which two peaks can be resolved, according to one of the standard measures
			*/
      float getResolution() const;
      /// sets the resolution
      void setResolution(float resolution);
			
			/// returns the mass accuracy i.e. how much the theoretical mass differs from the measured mass (in m/z)
      float getAccuracy() const;
      /// sets the accuracy  i.e. how much the theoretical mass differs from the measured mass  (in m/z)
      void setAccuracy(float accuracy);
			
			/// returns the scan rate (in s)
      float getScanRate() const;
      /// sets the scan rate (in s)
      void setScanRate(float scan_rate);
			
			/// returns the scan time for a single scan (in s)
      float getScanTime() const;
      /// sets the scan time for a single scan (in s)
      void setScanTime(float scan_time);
			
			/// returns the path length for a TOF mass analyzer (in mm)
      float getTOFTotalPathLength() const;
      /// sets the path length for a TOF mass analyzer (in mm)
      void setTOFTotalPathLength(float TOF_total_path_length);
			
			/// returns the isolation width i.e. in which m/z range the precursor ion is selected for MS to the n (in m/z)
      float getIsolationWidth() const;
      /// sets the isolation width i.e. in which m/z range the precursor ion is selected for MS to the n (in m/z)
      void setIsolationWidth(float isolation_width);
			
			/// returns the final MS exponent
      Int getFinalMSExponent() const;
      /// sets the final MS exponent
      void setFinalMSExponent(Int final_MS_exponent);
			
			/// returns the strength of the magnetic field (in T)
      float getMagneticFieldStrength() const;
      /// sets the strength of the magnetic field (in T)
      void setMagneticFieldStrength(float magnetic_field_strength);

    protected:
			AnalyzerType type_;
			ResolutionMethod resolution_method_;
			ResolutionType resolution_type_;
			ScanFunction scan_function_;
			ScanDirection scan_direction_;
			ScanLaw scan_law_;
			TandemScanningMethod tandem_scan_method_;
			ReflectronState reflectron_state_;
			float resolution_;
			float accuracy_;
			float scan_rate_;
			float scan_time_;
			float TOF_total_path_length_;
			float isolation_width_;
			Int final_MS_exponent_;
			float magnetic_field_strength_;
	};
} // namespace OpenMS

#endif // OPENMS_METADATA_MASSANALYZER_H
