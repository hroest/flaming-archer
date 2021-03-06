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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>

#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/util/XMLString.hpp>

namespace OpenMS
{

  OpenMS::Interfaces::SpectrumPtr MzMLSpectrumDecoder::decodeBinaryData(std::vector<BinaryData>& data_)
  {
    Internal::MzMLHandlerHelper::decodeBase64Arrays(data_);
    OpenMS::Interfaces::SpectrumPtr sptr(new OpenMS::Interfaces::Spectrum);

    //look up the precision and the index of the intensity and m/z array
    bool x_precision_64 = true;
    bool int_precision_64 = true;
    SignedSize x_index = -1;
    SignedSize int_index = -1;
    Internal::MzMLHandlerHelper::computeDataProperties_(data_, x_precision_64, x_index, "m/z array");
    Internal::MzMLHandlerHelper::computeDataProperties_(data_, int_precision_64, int_index, "intensity array");

    //Abort if no m/z or intensity array is present
    if (int_index == -1 || x_index == -1)
    {
      // Warning ...
      return sptr;
    }

    // Error if intensity or m/z is encoded as int32|64 - they should be float32|64!
    if ((data_[x_index].ints_32.size() > 0) || (data_[x_index].ints_64.size() > 0))
    {
      //fatalError(LOAD, "Encoding m/z array as integer is not allowed!");
    }
    if ((data_[int_index].ints_32.size() > 0) || (data_[int_index].ints_64.size() > 0))
    {
      //fatalError(LOAD, "Encoding intensity array as integer is not allowed!");
    }

    // Warn if the decoded data has a different size than the defaultArrayLength
    Size mz_size = x_precision_64 ? data_[x_index].floats_64.size() : data_[x_index].floats_32.size();
    Size int_size = int_precision_64 ? data_[int_index].floats_64.size() : data_[int_index].floats_32.size();
    // Check if int-size and mz-size are equal
    if (mz_size != int_size)
    {
      std::cout << "Warning, intensity and m/z array length are unequal" << std::endl;
      return sptr;
    }
    Size default_array_length_ = mz_size;
    // maybe some checks here about internal consistency e.g. with defaultArrayLength ...

    //create meta data arrays and reserve enough space for the content
    if (data_.size() > 2)
    {
      // --> TODO the other arrays ... create and resize those ...
    }

    // Copy meta data from m/z and intensity binary
    // We don't have this as a separate location => store it in spectrum
    // --> maybe TODO

    OpenMS::Interfaces::BinaryDataArrayPtr intensity_array(new OpenMS::Interfaces::BinaryDataArray);
    OpenMS::Interfaces::BinaryDataArrayPtr x_array(new OpenMS::Interfaces::BinaryDataArray);
    for (Size n = 0; n < default_array_length_; n++)
    {
      DoubleReal xcoord = x_precision_64 ? data_[x_index].floats_64[n] : data_[x_index].floats_32[n];
      DoubleReal intensity = int_precision_64 ? data_[int_index].floats_64[n] : data_[int_index].floats_32[n];

      x_array->data.push_back(xcoord);
      intensity_array->data.push_back(intensity);

      // TODO the other arrays
    }
    sptr->setMZArray(x_array);
    sptr->setIntensityArray(intensity_array);
    return sptr;
  }

  OpenMS::Interfaces::ChromatogramPtr MzMLSpectrumDecoder::decodeBinaryDataChrom(std::vector<BinaryData>& data_)
  {
    Internal::MzMLHandlerHelper::decodeBase64Arrays(data_);
    OpenMS::Interfaces::ChromatogramPtr sptr(new OpenMS::Interfaces::Chromatogram);

    //look up the precision and the index of the intensity and m/z array
    bool x_precision_64 = true;
    bool int_precision_64 = true;
    SignedSize x_index = -1;
    SignedSize int_index = -1;
    Internal::MzMLHandlerHelper::computeDataProperties_(data_, x_precision_64, x_index, "time array");
    Internal::MzMLHandlerHelper::computeDataProperties_(data_, int_precision_64, int_index, "intensity array");

    //Abort if no m/z or intensity array is present
    if (int_index == -1 || x_index == -1)
    {
      // Warning ...
      return sptr;
    }

    // Error if intensity or m/z is encoded as int32|64 - they should be float32|64!
    if ((data_[x_index].ints_32.size() > 0) || (data_[x_index].ints_64.size() > 0))
    {
      //fatalError(LOAD, "Encoding m/z array as integer is not allowed!");
    }
    if ((data_[int_index].ints_32.size() > 0) || (data_[int_index].ints_64.size() > 0))
    {
      //fatalError(LOAD, "Encoding intensity array as integer is not allowed!");
    }

    // Warn if the decoded data has a different size than the the defaultArrayLength
    Size mz_size = x_precision_64 ? data_[x_index].floats_64.size() : data_[x_index].floats_32.size();
    Size int_size = int_precision_64 ? data_[int_index].floats_64.size() : data_[int_index].floats_32.size();
    // Check if int-size and mz-size are equal
    if (mz_size != int_size)
    {
      std::cout << "Warning, intensity and m/z array length are unequal" << std::endl;
      return sptr;
    }
    Size default_array_length_ = mz_size;
    // maybe some checks here about internal consistency e.g. with defaultArrayLength ...

    //create meta data arrays and reserve enough space for the content
    if (data_.size() > 2)
    {
      // --> TODO the other arrays ... create and resize those ...
    }

    // Copy meta data from m/z and intensity binary
    // We don't have this as a separate location => store it in spectrum
    // --> maybe TODO

    OpenMS::Interfaces::BinaryDataArrayPtr intensity_array(new OpenMS::Interfaces::BinaryDataArray);
    OpenMS::Interfaces::BinaryDataArrayPtr x_array(new OpenMS::Interfaces::BinaryDataArray);
    for (Size n = 0; n < default_array_length_; n++)
    {
      DoubleReal xcoord = x_precision_64 ? data_[x_index].floats_64[n] : data_[x_index].floats_32[n];
      DoubleReal intensity = int_precision_64 ? data_[int_index].floats_64[n] : data_[int_index].floats_32[n];

      x_array->data.push_back(xcoord);
      intensity_array->data.push_back(intensity);

      // TODO the other arrays
    }
    sptr->setTimeArray(x_array);
    sptr->setIntensityArray(intensity_array);
    return sptr;
  }

  void MzMLSpectrumDecoder::handleBinaryDataArray(xercesc::DOMNode* indexListNode, std::vector<BinaryData>& data_)
  {
    // access result through data_.back()
    data_.push_back(BinaryData());

    XMLCh* TAG_CV = xercesc::XMLString::transcode("cvParam");
    XMLCh* TAG_binary = xercesc::XMLString::transcode("binary");
    XMLCh* TAG_userParam = xercesc::XMLString::transcode("userParam");
    XMLCh* TAG_referenceableParamGroupRef = xercesc::XMLString::transcode("referenceableParamGroupRef");

    // Iterate through binaryDataArray elements
    // only allowed subelements:
    //  - referenceableParamGroupRef (0+)
    //  - cvParam (0+)
    //  - userParam (0+)
    //  - binary (1)
    xercesc::DOMNodeList* index_elems = indexListNode->getChildNodes();
    const  XMLSize_t nodeCount_ = index_elems->getLength();
    for (XMLSize_t j = 0; j < nodeCount_; ++j)
    {
      xercesc::DOMNode* currentNode = index_elems->item(j);
      if (currentNode->getNodeType() &&   // true is not NULL
          currentNode->getNodeType() == xercesc::DOMNode::ELEMENT_NODE)   // is element
      {
        xercesc::DOMElement* currentElement = dynamic_cast<xercesc::DOMElement*>(currentNode);
        if (xercesc::XMLString::equals(currentElement->getTagName(), TAG_binary))
        {
          // TODO do not use expensive xercesc transcoding!
          data_.back().base64 = xercesc::XMLString::transcode(currentNode->getTextContent());
        }
        else if (xercesc::XMLString::equals(currentElement->getTagName(), TAG_CV))
        {
          std::string accession = xercesc::XMLString::transcode(currentElement->getAttribute(xercesc::XMLString::transcode("accession")));
          std::string value = xercesc::XMLString::transcode(currentElement->getAttribute(xercesc::XMLString::transcode("value")));
          std::string name = xercesc::XMLString::transcode(currentElement->getAttribute(xercesc::XMLString::transcode("name")));

          //handleCVParam(data_, accession, value, name);

          Internal::MzMLHandlerHelper::handleBinaryDataArrayCVParam(data_, accession, value, name);
        }
        else if (xercesc::XMLString::equals(currentElement->getTagName(), TAG_userParam))
        {
          std::cout << " unhandled userParam" << std::endl;
        }
        else if (xercesc::XMLString::equals(currentElement->getTagName(), TAG_referenceableParamGroupRef))
        {
          std::cout << " unhandled referenceableParamGroupRef" << std::endl;
        }
        else
        {
          //std::cout << "unhandled" << (string)xercesc::XMLString::transcode(currentNode->getNodeName() << std::endl;
        }
      }
    }

    xercesc::XMLString::release(&TAG_CV);
    xercesc::XMLString::release(&TAG_binary);
    xercesc::XMLString::release(&TAG_userParam);
    xercesc::XMLString::release(&TAG_referenceableParamGroupRef);
  }

  void MzMLSpectrumDecoder::domParseString(const std::string& in, std::vector<BinaryData>& data_)
  {
    //-------------------------------------------------------------
    // Create parser from input string using MemBufInputSource
    //-------------------------------------------------------------
    xercesc::MemBufInputSource myxml_buf(reinterpret_cast<const unsigned char*>(in.c_str()), in.length(), "myxml (in memory)");
    xercesc::XercesDOMParser* parser = new xercesc::XercesDOMParser();
    parser->setDoNamespaces(false);
    parser->setDoSchema(false);
    parser->setLoadExternalDTD(false);
    parser->parse(myxml_buf);

    //-------------------------------------------------------------
    // Start parsing
    // see http://www.yolinux.com/TUTORIALS/XML-Xerces-C.html
    //-------------------------------------------------------------

    // no need to free this pointer - owned by the parent parser object
    xercesc::DOMDocument* doc =  parser->getDocument();
    // Get the top-level element (for example <spectrum> or <chromatogram>)
    xercesc::DOMElement* elementRoot = doc->getDocumentElement();
    if (!elementRoot)
    {
      delete parser;
    }

    // Extract the binaryDataArray tag (there may be multiple)
    XMLCh* tag = xercesc::XMLString::transcode("binaryDataArray");
    xercesc::DOMNodeList* li = elementRoot->getElementsByTagName(tag);
    xercesc::XMLString::release(&tag);

    /// Do the processing of the DOM
    for (Size i = 0; i < li->getLength(); i++)
    {
      handleBinaryDataArray(li->item(i), data_);
    }

    delete parser;
  }

  void MzMLSpectrumDecoder::domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr& sptr)
  {
    std::vector<BinaryData> data_;
    domParseString(in, data_);
    sptr = decodeBinaryData(data_);
  }

  void MzMLSpectrumDecoder::domParseChromatogram(std::string& in, OpenMS::Interfaces::ChromatogramPtr& sptr)
  {
    std::vector<BinaryData> data_;
    domParseString(in, data_);
    sptr = decodeBinaryDataChrom(data_);
  }

}
