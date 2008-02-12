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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_STATISTICS_BASICSTATISTICS_H
#define OPENMS_MATH_STATISTICS_BASICSTATISTICS_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>
#include <ostream>
#include <cmath>
#include <numeric>

namespace OpenMS
{
  namespace Math
	{

		/**
			 @brief Calculates some basic statistical parameters of a distribution:
			 sum, mean, variance, and provides the normal approximation.
			 
			 The intended usage is as follows:
			 - <i>create</i> an instance
			 - <i>set</i> the basic statistical parameters by either
			 		- calling one of the update() member functions, or
					- using the set... methods
					.
			 - <i>do</i> something with the basic statistical parameters, e.g.
			 		- using the get... methods, or
					- obtain samples from a normal approximation with these parameters
					- whatever member function you might want to add to this class ;-)
			 .
			 
			 @todo Currently meanSquareError and pearsonCorrelationCoefficient can also be computed via static member functions.  This is against the design idea of this class (it's not using the statistics parameters).  Such functions should be moved to e.g. namespace OpenMS::Math and implemented in another file (how about MathFunctions.h?).  Please contact me if you contradict, I will change this soon.  (Clemens)

			 @ingroup Math
		*/
		template < typename RealT = double >
		class
		BasicStatistics
		{

		 public:

			/// The real type specified as template argument.
			typedef RealT RealType;

			typedef std::vector < RealType > probability_container;
			typedef std::vector < RealType > coordinate_container;

			/// Default constructor.
			BasicStatistics ()
				: mean_(0),
					variance_(0),
					sum_(0)
			{}

			/// Copy constructor.
			BasicStatistics ( BasicStatistics const & arg )
				: mean_(arg.mean_),
					variance_(arg.variance_),
					sum_(arg.sum_)
			{}

			/// Assignment.
			BasicStatistics & operator = ( BasicStatistics const & arg )
			{
				mean_      = arg.mean_;
				variance_  = arg.variance_;
				sum_       = arg.sum_;
				return *this;
			}

			/// Set sum, mean, and variance to zero.
			void clear ()
			{
				mean_ = 0;
				variance_ = 0;
				sum_ = 0;
			}

			/// This does the actual calculation.
			/** You can call this as often as you like, using different input vectors. */
			template < typename ProbabilityIterator >
			void update ( ProbabilityIterator probability_begin,
										ProbabilityIterator const probability_end
									)
			{
				clear();
				unsigned pos = 0;
				ProbabilityIterator iter = probability_begin;

				for ( ; iter != probability_end; ++iter, ++pos )
				{
					sum_  += *iter;
					mean_ += *iter * pos;
				}
				mean_ /= sum_;

				for ( iter = probability_begin, pos = 0; iter != probability_end; ++iter, ++pos )
				{
					RealType diff = RealType(pos) - mean_;
					diff *= diff;
					variance_ += *iter * diff;
				}
				variance_ /= sum_;
			}

			/// This does the actual calculation.
			/** You can call this as often as you like, using different input vectors. */
			template < typename ProbabilityIterator, typename CoordinateIterator >
			void update ( ProbabilityIterator const probability_begin,
										ProbabilityIterator const probability_end,
										CoordinateIterator  const coordinate_begin
									)
			{
				clear();
				ProbabilityIterator prob_iter = probability_begin;
				CoordinateIterator  coord_iter = coordinate_begin;

				for ( ; prob_iter != probability_end; ++prob_iter, ++coord_iter )
				{
					sum_  += *prob_iter;
					mean_ += *prob_iter * *coord_iter;
				}
				mean_ /= sum_;

				for ( prob_iter = probability_begin, coord_iter = coordinate_begin;
							prob_iter != probability_end;
							++prob_iter, ++coord_iter
						)
				{
					RealType diff = *coord_iter - mean_;
					diff *= diff;
					variance_ += *prob_iter * diff;
				}
				variance_ /= sum_;
				return;
			}

			/// Returns the mean.
			RealType mean () const { return mean_; }
			void setMean ( RealType const & mean ) { mean_ = mean; }

			/// Returns the variance.
			RealType variance() const { return variance_; }
			void setVariance ( RealType const & variance ) { variance_ = variance; }

			/// Returns the sum.
			RealType sum() const { return sum_; }
			void setSum ( RealType const & sum ) { sum_ = sum; }


			/**@brief Returns the density of the normal approximation at point,
				 multiplied by sqrt( 2 * pi ).  This saves a division operation compared
				 to normalDensity()
			*/
			RealType normalDensity_sqrt2pi ( RealType coordinate ) const throw()
			{
				coordinate -= mean();
				coordinate *= coordinate;
				return exp ( - coordinate / RealType(2) / variance() );
			}

			/// Returns sqrt( 2 * pi ), which is useful to normalize the result of normalDensity_sqrt2pi().
			static RealType sqrt2pi () { return 2.50662827463100050240; }

			/**@brief See normalDensity_sqrt2pi().  Returns the density of the normal
				 distribution at point.
			*/
			inline RealType normalDensity ( RealType const coordinate ) const throw()
			{
				return normalDensity_sqrt2pi ( coordinate ) / sqrt2pi() ;
			}

			/**@brief The argument \c probability is filled with values according to
				 the normal approximation.  Its \c size() is not changed.  The
				 approximation takes place at coordinate positions 0, 1, ..., size()-1.
			*/
			void normalApproximation ( probability_container & probability )
			{
				normalApproximationHelper_ ( probability, probability.size() );
				return;
			}

			/** The argument \c probability is filled with values according to the
					normal approximation.  Its size() is set to \c size.  The
					approximation takes place at coordinate positions 0, 1, ..., size-1.
			*/
			void normalApproximation ( probability_container & probability,
																 typename probability_container::size_type const size
															 )
			{
				probability.resize ( size );
				normalApproximationHelper_ ( probability, size );
				return;
			}

			/**@brief The argument probability is filled with values according to the
				 normal approximation.  The second argument coordinate contains the
				 positions where the approximation takes place.  probability.size() is
				 set to coordinate.size().
			*/
			void normalApproximation ( probability_container & probability,
																 coordinate_container  const & coordinate
															 )
			{
				probability.resize ( coordinate.size() );
				normalApproximationHelper_ ( probability, coordinate );
				return;
			}

			/**@brief A convenient overload for debugging purposes.

			@relatesalso BasicStatistics
			*/
			friend std::ostream & operator << ( std::ostream & os, BasicStatistics & arg )
			{
				os << "BasicStatistics:  mean=" << arg.mean() << "  variance=" << arg.variance() << "  sum=" << arg.sum();
				return os;
			}

			/**@brief Calculates the mean square error for the values in [begin_a, end_a) and [begin_b, end_b)

			Calculates the mean square error for the data given by the two iterator ranges. If
			one of the ranges contains a smaller number of values the rest of the longer range is omitted.
			*/
			template < typename IteratorType1, typename IteratorType2 >
			static RealType meanSquareError ( IteratorType1 begin_a, const IteratorType1 end_a,
																				IteratorType2 begin_b, const IteratorType2 end_b
																			)
			{
				Int count = 0;
				RealType error = 0;
				IteratorType1 & it_a = begin_a;
				IteratorType2 & it_b = begin_b;

				while(it_a != end_a && it_b != end_b)
				{
					RealType diff = *it_a - *it_b;
					error += diff * diff;
					++count;
					++it_a;
					++it_b;
				}

				return error / count;

			}

			/**@brief Calculates the classification rate for the values in [begin_a, end_a) and [begin_b, end_b)

			Calculates the classification rate for the data given by the two iterator ranges. If
			one of the ranges contains a smaller number of values the rest of the longer range is omitted.
			*/
			template < typename IteratorType1, typename IteratorType2 >
			static RealType classificationRate ( IteratorType1 begin_a, const IteratorType1 end_a,
																					 IteratorType2 begin_b, const IteratorType2 end_b
																				 )
			{
				Int count = 0;
				RealType error = 0;
				IteratorType1 & it_a = begin_a;
				IteratorType2 & it_b = begin_b;

    	  if (it_a == end_a || it_b == end_b)
    	  {
    	  	return 0;
    	  }

				while(it_a != end_a && it_b != end_b)
				{
					if ((*it_a < 0 && *it_b >= 0)
							|| (*it_a >= 0 && *it_b < 0))
					{
						error += 1;
					}
					++count;
					++it_a;
					++it_b;
				}

				return (count - error) / count;

			}

			/**
				@brief calculates the pearson correlation coefficient for the values in [begin_a, end_a) and [begin_b, end_b)

				Calculates the linear correlation coefficient for the data given by the two iterator ranges. 

				If the iterator ranges are not of the same length or empty an exception is thrown.
				
				If one of the ranges contains only the same values 'nan' is returned.	
			*/
			template < typename IteratorType1, typename IteratorType2 >
			static RealType pearsonCorrelationCoefficient ( IteratorType1 begin_a, IteratorType1 end_a, IteratorType2 begin_b, IteratorType2 end_b )
			throw (Exception::InvalidRange)
			{
				UInt count = end_a-begin_a;
				//no data or different lengths
				if (count==0 || end_a-begin_a!=end_b-begin_b)
				{
				  throw Exception::InvalidRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}
				
				//calculate average
				RealType avg_a = std::accumulate(begin_a,end_a,0.0) / count;
				RealType avg_b = std::accumulate(begin_b,end_b,0.0) / count;

				RealType numerator = 0;
				RealType denominator_a = 0;
				RealType denominator_b = 0;
				while(begin_a != end_a)
				{
					RealType temp_a = *begin_a - avg_a;
					RealType temp_b = *begin_b - avg_b;
					numerator += (temp_a * temp_b);
					denominator_a += (temp_a * temp_a);
					denominator_b += (temp_b * temp_b);
					++begin_a;
					++begin_b;
				}
				
				return numerator / sqrt(denominator_a * denominator_b);
			}

		 protected:

			/// @name Protected Members
			//@{

			RealType mean_;
			RealType variance_;
			RealType sum_;

		 private:
			//@}

			/// @name Private Methods
			//@{

			void normalApproximationHelper_ ( probability_container & probability,
																				typename probability_container::size_type const size
																			)
			{
				RealType gaussSum = 0;
				typename coordinate_container::size_type i;

				// precondition size == probability.size() is guaranteed by wrappers.
				for ( i = 0; i < size; ++i ) {
					gaussSum += normalDensity_sqrt2pi ( i );
				}

				for ( i = 0; i < size; ++i ) {
					probability [ i ] = normalDensity_sqrt2pi ( i ) / gaussSum * sum();
				}
				return;
			}

			void normalApproximationHelper_ ( probability_container & probability,
																				coordinate_container const & coordinate
																			)
			{
				RealType gaussSum = 0;
				typename coordinate_container::size_type i;
				typename coordinate_container::size_type const size = coordinate.size();

				for ( i = 0; i < size; ++i ) {
					gaussSum += normalDensity_sqrt2pi ( coordinate[i] );
				}

				for ( i = 0; i < size; ++i ) {
					probability [ i ] = normalDensity_sqrt2pi ( coordinate[i] ) / gaussSum * sum();
				}
				return;
			}

			//@}

		};

	} // namespace Math

} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_BASICSTATISTICS_H
