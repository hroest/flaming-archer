// -*- mode: C++; tab-width: 2; -*-
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


#ifndef OPENMS_VISUAL_SPECTRUM2DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM2DCANVAS_H

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/KERNEL/PeakIndex.h>

// QT
class QPainter;
class QMouseEvent;

namespace OpenMS
{
  /**
  	@brief Canvas for 2D-visualization of peak map, feature map and consensus map data

  	This widget displays a 2D representation of a set of peaks, features or consensus elements.
		
		@image html Spectrum2DCanvas.png
		
		The example image shows %Spectrum2DCanvas displaying a peak layer and a feature layer. 
  	
		@htmlinclude OpenMS_Spectrum2DCanvas.parameters

    @improvement Support different peak icons - cross, star, square, ... (HiWi)
		
  	@ingroup SpectrumWidgets
  */
  class OPENMS_DLLAPI Spectrum2DCanvas 
  	: public SpectrumCanvas
  {
      Q_OBJECT

    public:
      /// Default constructor
      Spectrum2DCanvas(const Param& preferences, QWidget* parent = 0);

      /// Destructor
      ~Spectrum2DCanvas();

			// Docu in base class
			virtual void showCurrentLayerPreferences();

			// Docu in base class
			virtual void saveCurrentLayer(bool visible);
			
			/// Merges the features in @p map into the features layer @p i 
			void mergeIntoLayer(Size i, FeatureMapType& map);

			/// Merges the consensus features in @p map into the features layer @p i 
			void mergeIntoLayer(Size i, ConsensusMapType& map);
			
    signals:
      /// Sets the data for the horizontal projection
      void showProjectionHorizontal(const ExperimentType&, Spectrum1DCanvas::DrawModes);
      /// Sets the data for the vertical projection
      void showProjectionVertical(const ExperimentType&, Spectrum1DCanvas::DrawModes);
      /// Shows the number of peaks and the intensity sum of the projection
      void showProjectionInfo(int, double, double);
			/// Signal emitted when the projections are to be shown/hidden 
			void toggleProjections();
			/// Requests to display the spectrum with index @p index in 1D
			void showSpectrumAs1D(int index);
		
    public slots:

      // Docu in base class
      void activateLayer(int layer_index);
      // Docu in base class
      void removeLayer(int layer_index);
      // Docu in base class
      virtual void horizontalScrollBarChange(int value);
      // Docu in base class
      virtual void verticalScrollBarChange(int value);
      /**
      	@brief Updates the projection data and emits some related signals.
      	
      	Emitted signals are showProjectionHorizontal(const ExperimentType&, Spectrum1DCanvas::DrawModes) and 
      	showProjectionVertical(const ExperimentType&, Spectrum1DCanvas::DrawModes).
      	
      	@see projection_mz_
      	@see projection_rt_
      */
      void updateProjections();
      
    protected:
      // Docu in base class
      bool finishAdding_();
      
      /** @name Reimplemented QT events */
      //@{
      void mousePressEvent(QMouseEvent* e);
      void mouseReleaseEvent(QMouseEvent* e);
      void mouseMoveEvent(QMouseEvent* e);
			void paintEvent(QPaintEvent* e);
			void contextMenuEvent(QContextMenuEvent* e);
			void keyPressEvent(QKeyEvent* e);
      void keyReleaseEvent(QKeyEvent* e);
			void mouseDoubleClickEvent(QMouseEvent* e); 
      //@}

      // Docu in base class
      virtual void updateScrollbars_();

      /**
      	@brief Paints individual peaks.

      	Paints the peaks as small ellipses. The peaks are colored according to the
      	selected dot gradient.
      	
      	@param layer_index The index of the layer.
      	@param p The QPainter to paint on.
      */
      void paintDots_(Size layer_index, QPainter& p);

      /**
      	@brief Paints convex hulls (one for each mass trace) of a features layer.
      	
      	@param layer_index Index of the layer.
      	@param p The QPainter to paint on.
      */
      void paintTraceConvexHulls_(Size layer_index, QPainter& p);

      /**
      	@brief Paints the convex hulls (one for each feature) of a features layer.
      	
      	@param layer_index Index of the layer.
      	@param p The QPainter to paint on.
      */
      void paintFeatureConvexHulls_(Size layer_index, QPainter& p);

      /**
      	@brief Paints the consensus elements of a consensus features layer.
      	
      	@param layer_index Index of the layer.
      	@param p The QPainter to paint on.
      */
      void paintConsensusElements_(Size layer_index, QPainter& p);

      /**
      	@brief Paints one consensus element of a consensus features layer.
      	
      	@param layer_index Index of the layer.
      	@param cf Reference to the feature to be painted.
      	@param p The QPainter to paint on.
      	@param use_buffer Flag to switch between painting on the buffer and screen.
      */		
			void paintConsensusElement_(Size layer_index, const ConsensusFeature& cf, QPainter& p, bool use_buffer);
			
      /**
      	@brief checks if any element of a consensus feature is currently visible.
      	
      	@param layer_index Index of the layer.
      	@param ce The ConsensusFeature that needs checking
      */
			bool isConsensusFeatureVisible_(const ConsensusFeature& ce, Size layer_index);

			/**
      	@brief Paints convex hulls (one for each mass trace) for a single feature.
      	
      	@param hulls Reference to convex hull vector.
      	@param p The QPainter to paint on.
      */
      void paintConvexHulls_(const std::vector<ConvexHull2D>& hulls, QPainter& p);

      // Docu in base class
      virtual void intensityModeChange_();
      // DOcu in base class
      virtual void recalculateSnapFactor_();
			// Docu in base class
			virtual void currentLayerParamtersChanged_();
      /// recalculates the dot gradient of a layer
      void recalculateDotGradient_(Size layer);

      /// m/z projection data
      ExperimentType projection_mz_;
      /// RT projection data
      ExperimentType projection_rt_;

      /**
      	@brief Returns the color associated with @p val for the gradient @p gradient.
      	
      	Takes intensity modes into account.
      */
      inline const QColor& heightColor_(Real val, const MultiGradient& gradient, DoubleReal snap_factor)
			{
				switch (intensity_mode_)
				{
					case IM_NONE:
						return gradient.precalculatedColorAt(val);
						break;
					case IM_PERCENTAGE:
						return gradient.precalculatedColorAt(val*percentage_factor_);
						break;
					case IM_SNAP:
						return gradient.precalculatedColorAt(val*snap_factor);
						break;
					default:
						throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
			}

      /// Highlights a single peak
      void highlightPeak_(QPainter& p, const PeakIndex& peak);

      /// Returns the nearest peak to position @p pos
      PeakIndex findNearestPeak_(const QPoint& pos);

      /// the nearest peak/feature to the mouse cursor
			PeakIndex selected_peak_;
      /// start peak/feature of measuring mode
      PeakIndex measurement_start_;

			//docu in base class
			virtual void updateLayer_(Size i);
      //docu in base class
			virtual void translateLeft_();
			//docu in base class
			virtual void translateRight_();
			//docu in base class
			virtual void translateForward_();
			//docu in base class
			virtual void translateBackward_();
			
			/// Finishes context menu after customization to peaks, features or consensus features
			void finishContextMenu_(QMenu* context_menu, QMenu* settings_menu);
  };
}

#endif
