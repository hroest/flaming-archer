// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Matrix.h>

namespace OpenMS
{
  Matrix<int>    default_matrix_int;
  Matrix<double> default_matrix_double;

#if 0
  template <>
  deprecated_gsl_matrix * Matrix<double>::toGslMatrix()
  {
    deprecated_gsl_matrix * m_ptr = deprecated_gsl_matrix_alloc(rows_, cols_);
    for (size_type i = 0; i < this->rows_; ++i)
    {
      for (size_type j = 0; j < this->cols_; ++j)
      {
        deprecated_gsl_matrix_set(m_ptr, i, j, (*this)(i, j));
      }
    }
    return m_ptr;
  }

  template <>
  deprecated_gsl_matrix * Matrix<float>::toGslMatrix()
  {
    deprecated_gsl_matrix * m_ptr = deprecated_gsl_matrix_alloc(rows_, cols_);
    for (size_type i = 0; i < this->rows_; ++i)
    {
      for (size_type j = 0; j < this->cols_; ++j)
      {
        deprecated_gsl_matrix_set(m_ptr, i, j, (double) (*this)(i, j));
      }
    }
    return m_ptr;
  }

#endif

}
