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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_FEATUREMAP_H
#define OPENMS_KERNEL_FEATUREMAP_H

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/CONCEPT/UniqueIdIndexer.h>

#include <algorithm>
#include <vector>
#include <exception>

namespace OpenMS
{



  /// summary of the peptide identification assigned to each feature of this map.
  /// Each feature contributes one vote (=state)
  struct AnnotationStatistics
  {
    std::vector<Size> states; //< count each state, indexing by BaseFeature::AnnotationState

    AnnotationStatistics()
      : states(BaseFeature::SIZE_OF_ANNOTATIONSTATE, 0) // initialize all with 0
    {
    }

    AnnotationStatistics(const AnnotationStatistics& rhs)
      : states(rhs.states)
    {
    }

    AnnotationStatistics& operator=(const AnnotationStatistics& rhs)
    {
      if (this == &rhs) return *this;

      states = rhs.states;
      return *this;
    }

    bool operator==(const AnnotationStatistics& rhs) const
    {
      return states == rhs.states;
    }

    AnnotationStatistics& operator+=(BaseFeature::AnnotationState state)
    {
      ++states[(Size)state];
      return *this;
    }

  };

  
  /// Print content of an AnnotationStatistics object to a stream
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const AnnotationStatistics& ann);

  /**
    @brief A container for features.

    A map is a container holding 2-dimensional features,
    which in turn represent chemical entities (peptides, proteins, etc.) found
    in a 2-dimensional experiment.

    Maps are implemented as vectors of features and have basically the same interface
    as an STL vector has (model of Random Access Container and Back Insertion Sequence).

    Feature maps are typically created from peak data of 2D runs through the FeatureFinder.

    @ingroup Kernel
  */
  template <typename FeatureT = Feature>
  class FeatureMap :
    private std::vector<FeatureT>,
    public RangeManager<2>,
    public DocumentIdentifier,
    public UniqueIdInterface,
    public UniqueIdIndexer<FeatureMap<FeatureT> >
  {
public:
    /**
      @name Type definitions
    */
    typedef std::vector<FeatureT> privvec;

    // types
    using typename privvec::value_type;
    using typename privvec::iterator;
    using typename privvec::const_iterator;
    using typename privvec::size_type;
    using typename privvec::pointer;          // ConstRefVector
    using typename privvec::reference;        // ConstRefVector
    using typename privvec::const_reference;  // ConstRefVector
    using typename privvec::difference_type;  // ConstRefVector
 
    // functions
    using privvec::begin; 
    using privvec::end; 

    using privvec::size; 
    using privvec::resize;  // ConsensusMap, FeatureXMLFile
    using privvec::empty; 
    using privvec::reserve; 
    using privvec::operator[]; 
    using privvec::at;    // UniqueIdIndexer
    using privvec::back;  // FeatureXMLFile

    using privvec::push_back; 
    using privvec::pop_back;  // FeatureXMLFile
    using privvec::erase;     // source/VISUAL/Spectrum2DCanvas.C 2871, FeatureMap_test 599

    //@{
    typedef FeatureT FeatureType;
    typedef RangeManager<2> RangeManagerType;
    typedef std::vector<FeatureType> Base;
    typedef typename Base::iterator Iterator;
    typedef typename Base::const_iterator ConstIterator;
    typedef typename Base::reverse_iterator ReverseIterator;
    typedef typename Base::const_reverse_iterator ConstReverseIterator;
    typedef FeatureType & Reference;
    typedef const FeatureType & ConstReference;
    //@}

    /**
      @name Constructors and Destructor
    */
    //@{

    /// Default constructor
    FeatureMap() :
      Base(),
      RangeManagerType(),
      DocumentIdentifier(),
      UniqueIdInterface(),
      UniqueIdIndexer<FeatureMap<FeatureT> >(),
      protein_identifications_(),
      unassigned_peptide_identifications_(),
      data_processing_()
    {}

    /// Copy constructor
    FeatureMap(const FeatureMap & source) :
      Base(source),
      RangeManagerType(source),
      DocumentIdentifier(source),
      UniqueIdInterface(source),
      UniqueIdIndexer<FeatureMap<FeatureT> >(source),
      protein_identifications_(source.protein_identifications_),
      unassigned_peptide_identifications_(source.unassigned_peptide_identifications_),
      data_processing_(source.data_processing_)
    {}

    /// Destructor
    virtual ~FeatureMap()
    {}
    //@}

    /// Assignment operator
    FeatureMap & operator=(const FeatureMap & rhs)
    {
      if (&rhs == this) return *this;

      Base::operator=(rhs);
      RangeManagerType::operator=(rhs);
      DocumentIdentifier::operator=(rhs);
      UniqueIdInterface::operator=(rhs);
      protein_identifications_ = rhs.protein_identifications_;
      unassigned_peptide_identifications_ = rhs.unassigned_peptide_identifications_;
      data_processing_ = rhs.data_processing_;

      return *this;
    }

    /// Equality operator
    bool operator==(const FeatureMap & rhs) const
    {
      return std::operator==(*this, rhs) &&
             RangeManagerType::operator==(rhs) &&
             DocumentIdentifier::operator==(rhs) &&
             UniqueIdInterface::operator==(rhs) &&
             protein_identifications_ == rhs.protein_identifications_ &&
             unassigned_peptide_identifications_ == rhs.unassigned_peptide_identifications_ &&
             data_processing_ == rhs.data_processing_;
    }

    /// Equality operator
    bool operator!=(const FeatureMap & rhs) const
    {
      return !(operator==(rhs));
    }

    /**
      @brief Joins two feature maps.

      Features are merged into one container (see operator+= for details).
    */
    FeatureMap operator+(const FeatureMap & rhs) const
    {
      FeatureMap tmp(*this);
      tmp += rhs;
      return tmp;
    }

    /**
      @brief Add one feature map to another.

      Features are merged into one container, simply by appending.
      UnassignedPeptides and ProteinIdentifications are appended.
      Information on DocumentIdentifier, UniqueIdInterface (of container only)
      are reset to default.

      For conflicting UID's, new UID's will be assigned.

      @param rhs The feature to add to this one.
    */
    FeatureMap & operator+=(const FeatureMap & rhs)
    {
      FeatureMap empty_map;
      // reset these:
      RangeManagerType::operator=(empty_map);

      if (!this->getIdentifier().empty() || !rhs.getIdentifier().empty()) LOG_INFO << "DocumentIdentifiers are lost during merge of FeatureMaps\n";
      DocumentIdentifier::operator=(empty_map);

      UniqueIdInterface::operator=(empty_map);

      // merge these:
      protein_identifications_.insert(protein_identifications_.end(), rhs.protein_identifications_.begin(), rhs.protein_identifications_.end());
      unassigned_peptide_identifications_.insert(unassigned_peptide_identifications_.end(), rhs.unassigned_peptide_identifications_.begin(), rhs.unassigned_peptide_identifications_.end());
      data_processing_.insert(data_processing_.end(), rhs.data_processing_.begin(), rhs.data_processing_.end());

      // append features:
      this->insert(this->end(), rhs.begin(), rhs.end());

      // todo: check for double entries
      // features, unassignedpeptides, proteins...

      // consistency
      try
      {
        UniqueIdIndexer<FeatureMap<FeatureT> >::updateUniqueIdToIndex();
      }
      catch (Exception::Postcondition /*&e*/) // assign new UID's for conflicting entries
      {
        Size replaced_uids =  UniqueIdIndexer<FeatureMap<FeatureT> >::resolveUniqueIdConflicts();
        LOG_INFO << "Replaced " << replaced_uids << " invalid uniqueID's\n";
      }

      return *this;
    }

    /**
      @name Sorting.
      These simplified sorting methods are supported in addition to
      the standard sorting methods of std::vector.
    */
    //@{
    /// Sorts the peaks according to ascending intensity.
    void sortByIntensity(bool reverse = false)
    {
      if (reverse)
      {
        std::sort(this->begin(), this->end(), reverseComparator(typename FeatureType::IntensityLess()));
      }
      else
      {
        std::sort(this->begin(), this->end(), typename FeatureType::IntensityLess());
      }
    }

    ///Sort features by position. Lexicographical comparison (first RT then m/z) is done.
    void sortByPosition()
    {
      std::sort(this->begin(), this->end(), typename FeatureType::PositionLess());
    }

    ///Sort features by RT position.
    void sortByRT()
    {
      std::sort(this->begin(), this->end(), typename FeatureType::RTLess());
    }

    ///Sort features by m/z position.
    void sortByMZ()
    {
      std::sort(this->begin(), this->end(), typename FeatureType::MZLess());
    }

    ///Sort features by ascending overall quality.
    void sortByOverallQuality(bool reverse = false)
    {
      if (reverse)
      {
        std::sort(this->begin(), this->end(), reverseComparator(typename FeatureType::OverallQualityLess()));
      }
      else
      {
        std::sort(this->begin(), this->end(), typename FeatureType::OverallQualityLess());
      }
    }

    //@}

    // Docu in base class
    void updateRanges()
    {
      this->clearRanges();
      updateRanges_(this->begin(), this->end());

      //enlarge the range by the convex hull points
      for (Size i = 0; i < this->size(); ++i)
      {
        DBoundingBox<2> box = this->operator[](i).getConvexHull().getBoundingBox();
        if (!box.isEmpty())
        {
          //update RT
          if (box.minPosition()[Peak2D::RT] < this->pos_range_.minPosition()[Peak2D::RT])
          {
            this->pos_range_.setMinX(box.minPosition()[Peak2D::RT]);
          }
          if (box.maxPosition()[Peak2D::RT] > this->pos_range_.maxPosition()[Peak2D::RT])
          {
            this->pos_range_.setMaxX(box.maxPosition()[Peak2D::RT]);
          }
          //update m/z
          if (box.minPosition()[Peak2D::MZ] < this->pos_range_.minPosition()[Peak2D::MZ])
          {
            this->pos_range_.setMinY(box.minPosition()[Peak2D::MZ]);
          }
          if (box.maxPosition()[Peak2D::MZ] > this->pos_range_.maxPosition()[Peak2D::MZ])
          {
            this->pos_range_.setMaxY(box.maxPosition()[Peak2D::MZ]);
          }
        }
      }
    }

    /// Swaps the feature content (plus its range information) of this map with the content of @p from
    void swapFeaturesOnly(FeatureMap& from)
    {
      // TODO used by FeatureFinderAlgorithmPicked -- could it also use regular swap?
      Base::swap(from);
      
      // swap range information (otherwise its false in both maps)
      FeatureMap tmp;
      tmp.RangeManagerType::operator=(* this);
      this->RangeManagerType::operator=(from);
      from.RangeManagerType::operator=(tmp);
    }

    void swap(FeatureMap& from)
    {
      // swap features and ranges
      swapFeaturesOnly(from);

      // swap DocumentIdentifier
      DocumentIdentifier::swap(from);

      // swap unique id
      UniqueIdInterface::swap(from);

      // swap unique id index
      UniqueIdIndexer<FeatureMap<FeatureT> >::swap(from);

      // swap the remaining members
      protein_identifications_.swap(from.protein_identifications_);
      unassigned_peptide_identifications_.swap(from.unassigned_peptide_identifications_);
      data_processing_.swap(from.data_processing_);
    }

    /// non-mutable access to the protein identifications
    const std::vector<ProteinIdentification> & getProteinIdentifications() const
    {
      return protein_identifications_;
    }

    /// mutable access to the protein identifications
    std::vector<ProteinIdentification> & getProteinIdentifications()
    {
      return protein_identifications_;
    }

    /// sets the protein identifications
    void setProteinIdentifications(const std::vector<ProteinIdentification> & protein_identifications)
    {
      protein_identifications_ = protein_identifications;
    }

    /// non-mutable access to the unassigned peptide identifications
    const std::vector<PeptideIdentification> & getUnassignedPeptideIdentifications() const
    {
      return unassigned_peptide_identifications_;
    }

    /// mutable access to the unassigned peptide identifications
    std::vector<PeptideIdentification> & getUnassignedPeptideIdentifications()
    {
      return unassigned_peptide_identifications_;
    }

    /// sets the unassigned peptide identifications
    void setUnassignedPeptideIdentifications(const std::vector<PeptideIdentification> & unassigned_peptide_identifications)
    {
      unassigned_peptide_identifications_ = unassigned_peptide_identifications;
    }

    /// returns a const reference to the description of the applied data processing
    const std::vector<DataProcessing> & getDataProcessing() const
    {
      return data_processing_;
    }

    /// returns a mutable reference to the description of the applied data processing
    std::vector<DataProcessing> & getDataProcessing()
    {
      return data_processing_;
    }

    /// sets the description of the applied data processing
    void setDataProcessing(const std::vector<DataProcessing> & processing_method)
    {
      data_processing_ = processing_method;
    }

    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    void clear(bool clear_meta_data = true)
    {
      Base::clear();

      if (clear_meta_data)
      {
        clearRanges();
        this->DocumentIdentifier::operator=(DocumentIdentifier());             // no "clear" method
        clearUniqueId();
        protein_identifications_.clear();
        unassigned_peptide_identifications_.clear();
        data_processing_.clear();
      }
    }

    /**
      @brief Applies a member function of Type to the container itself and all features (including subordinates).
      The returned values are accumulated.

      <b>Example:</b>  The following will print the number of features with invalid unique ids (plus 1 if the container has an invalid UID as well):
      @code
      FeatureMap<> fm;
      (...)
      std::cout << fm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId) << std::endl;
      @endcode
      See e.g. UniqueIdInterface for what else can be done this way.
    */
    template <typename Type>
    Size applyMemberFunction(Size (Type::* member_function)())
    {
      Size assignments = 0;
      assignments += ((*this).*member_function)();
      for (Iterator iter = this->begin(); iter != this->end(); ++iter)
      {
        assignments += iter->applyMemberFunction(member_function);
      }
      return assignments;
    }

    /// The "const" variant.
    template <typename Type>
    Size applyMemberFunction(Size (Type::* member_function)() const) const
    {
      Size assignments = 0;
      assignments += ((*this).*member_function)();
      for (ConstIterator iter = this->begin(); iter != this->end(); ++iter)
      {
        assignments += iter->applyMemberFunction(member_function);
      }
      return assignments;
    }

    AnnotationStatistics getAnnotationStatistics() const
    {
      AnnotationStatistics result;
      for (ConstIterator iter = this->begin(); iter != this->end(); ++iter)
      {
        result += iter->getAnnotationState();
      }
      return result;
    }

protected:

    /// protein identifications
    std::vector<ProteinIdentification> protein_identifications_;

    /// peptide identifications not matched to a specific feature
    std::vector<PeptideIdentification> unassigned_peptide_identifications_;

    /// applied data processing
    std::vector<DataProcessing> data_processing_;
  };

  /// Print content of a feature map to a stream.
  template <typename FeatureType>
  std::ostream & operator<<(std::ostream & os, const FeatureMap<FeatureType> & map)
  {
    os << "# -- DFEATUREMAP BEGIN --" << "\n";
    os << "# POS \tINTENS\tOVALLQ\tCHARGE\tUniqueID" << "\n";
    for (typename FeatureMap<FeatureType>::const_iterator iter = map.begin(); iter != map.end(); ++iter)
    {
      os << iter->getPosition() << '\t'
      << iter->getIntensity() << '\t'
      << iter->getOverallQuality() << '\t'
      << iter->getCharge() << '\t'
      << iter->getUniqueId() << "\n";
    }
    os << "# -- DFEATUREMAP END --" << std::endl;
    return os;
  }


} // namespace OpenMS

#endif // OPENMS_KERNEL_DFEATUREMAP_H
