// NA6PCCopyright

// Based on:
// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef NA6P_MCTRUTHCONTAINER_H
#define NA6P_MCTRUTHCONTAINER_H

#include "NA6PMCComposedLabel.h"
#include <span>
#include <stdexcept>

class NA6PMCTruthContainer
{
 private:
  std::vector<uint32_t> mHeaderArray; // the header structure array serves as an index into the actual storage
  std::vector<NA6PMCComposedLabel> mTruthArray;   // the buffer containing the actual truth information
  std::vector<char> mStreamerData; // buffer used for streaming a flat raw buffer
  
  size_t getSize(uint32_t dataindex) const
  {
    // calculate size / number of labels from a difference in pointed indices
    const auto size = (dataindex < getIndexedSize() - 1)
      ? getMCTruthHeader(dataindex + 1) - getMCTruthHeader(dataindex)
      : getNElements() - getMCTruthHeader(dataindex);
    return size;
  }

 public:
  // constructor
  NA6PMCTruthContainer() = default;
  // destructor
  ~NA6PMCTruthContainer() = default;
  // copy constructor
  NA6PMCTruthContainer(const NA6PMCTruthContainer& other) = default;
  // move constructor
  NA6PMCTruthContainer(NA6PMCTruthContainer&& other) = default;
  // construct from raw data
  NA6PMCTruthContainer(std::vector<uint32_t>& header, std::vector<NA6PMCComposedLabel>& truthArray)
  {
    setFrom(header, truthArray);
  }
  // assignment operator
  NA6PMCTruthContainer& operator=(const NA6PMCTruthContainer& other) = default;
  // move assignment operator
  NA6PMCTruthContainer& operator=(NA6PMCTruthContainer&& other) = default;

  // Set container directly from header + truthArray.
  void setFrom(std::vector<uint32_t>& header, std::vector<NA6PMCComposedLabel>& truthArray)
  {
    mHeaderArray = std::move(header);
    mTruthArray = std::move(truthArray);
  }

  struct FlatHeader {
    uint8_t version = 1;
    uint8_t sizeofHeaderElement = sizeof(uint32_t);
    uint8_t sizeofTruthElement = sizeof(NA6PMCComposedLabel);
    uint8_t reserved = 0;
    uint32_t nofHeaderElements;
    uint32_t nofTruthElements;
  };

  // access
  uint32_t getMCTruthHeader(uint32_t dataindex) const { return mHeaderArray[dataindex]; }
  // access the element directly (can be encapsulated better away)... needs proper element index
  // which can be obtained from the MCTruthHeader startposition and size
  NA6PMCComposedLabel const& getElement(uint32_t elementindex) const { return mTruthArray[elementindex]; }
  // return the number of original data indexed here
  size_t getIndexedSize() const { return mHeaderArray.size(); }
  // return the number of elements (labels) pointed to in this container
  size_t getNElements() const { return mTruthArray.size(); }
  // return unterlaying vector of elements
  const std::vector<NA6PMCComposedLabel>& getTruthArray() const
  {
    return mTruthArray;
  }

  // get individual "view" container for a given data index
  // the caller can do modifications on this view (such as sorting)
  std::span<NA6PMCComposedLabel> getLabels(uint32_t dataindex)
  {
    if (dataindex >= getIndexedSize()) {
      return std::span<NA6PMCComposedLabel>();
    }
    return std::span<NA6PMCComposedLabel>(&mTruthArray[getMCTruthHeader(dataindex)], getSize(dataindex));
  }
  
  // get individual const "view" container for a given data index
  // the caller can't do modifications on this view
  std::span<const NA6PMCComposedLabel> getLabels(uint32_t dataindex) const
  {
    if (dataindex >= getIndexedSize()) {
      return std::span<const NA6PMCComposedLabel>();
    }
    return std::span<const NA6PMCComposedLabel>(&mTruthArray[getMCTruthHeader(dataindex)], getSize(dataindex));
  }

  void clear()
  {
    mHeaderArray.clear();
    mTruthArray.clear();
  }
  void clear_andfreememory()
  {
    clear();
    // this forces the desctructor being called on existing buffers
    mHeaderArray = std::vector<uint32_t>();
    mTruthArray = std::vector<NA6PMCComposedLabel>();
    mStreamerData = std::vector<char>();
  }
  // add element for a particular dataindex
  // at the moment only strictly consecutive modes are supported
  // noElement indicates that actually no label/element will be added for this dataindex
  void addElement(uint32_t dataindex, NA6PMCComposedLabel const& element, bool noElement = false)
  {
    if (dataindex < mHeaderArray.size()) {
      // look if we have something for this dataindex already
      // must currently be the last one
      // because mTruthArray is a sequential array and
      // labels for a given data-index must occupy a contiguous block.
      if (dataindex != (mHeaderArray.size() - 1)) {
        throw std::runtime_error("NA6PMCTruthContainer: unsupported code path");
      }
    } else {
      // add empty holes = how many data-indices were skipped (had zero labels)
      // between the last one filled and this new one
      int holes = dataindex - mHeaderArray.size();
      assert(holes >= 0);
      for (int i = 0; i < holes; ++i) {
        // For each skipped index, a header entry is created, pointing to the end of mTruthArray 
        mHeaderArray.emplace_back(mTruthArray.size());
      }
      // add the header entry for the new dataindex 
      mHeaderArray.emplace_back(mTruthArray.size());
    }
    if (!noElement) {
      // append the new label to the flat array 
      mTruthArray.emplace_back(element);
    }
  }
  /// adds a data index that has no label
  void addNoLabelIndex(uint32_t dataindex)
  {
    addElement(dataindex, NA6PMCComposedLabel(), true);
  }
  // convenience interface to add multiple labels at once
  void addElements(uint32_t dataindex, const std::vector<NA6PMCComposedLabel>& elements)
  {
    addNoLabelIndex(dataindex);
    for (auto& e : elements) {
      addElement(dataindex, e);
    }
  }
  // Add element at last position or for a previous index
  // (at random access position).
  // This might be a slow process since data has to be moved internally
  // so this function should be used with care.
  void addElementRandomAccess(uint32_t dataindex, NA6PMCComposedLabel const& element)
  {
    if (dataindex > mHeaderArray.size()) {
      throw std::runtime_error("NA6PMCTruthContainer: addElementRandomAccess does not support holes");
    } else if (dataindex == mHeaderArray.size() ) {
      // a new dataindex -> push element at back
      mHeaderArray.resize(dataindex + 1);
      mHeaderArray[dataindex] = mTruthArray.size();
      mTruthArray.emplace_back(element);
    } else if (dataindex == mHeaderArray.size() - 1) {
      // if appending at end use fast function
      addElement(dataindex, element);
    } else {
      // existing dataindex
      // have to:
      // a) move data;
      // b) insert new element;
      // c) adjust indices of all headers right to this
      auto currentindex = mHeaderArray[dataindex];
      auto lastindex = currentindex + getSize(dataindex);
      // resize truth array
      mTruthArray.resize(mTruthArray.size() + 1);
      // move data (have to do from right to left)
      for (size_t i = mTruthArray.size() - 1; i > lastindex; --i) {
        mTruthArray[i] = mTruthArray[i - 1];
      }
      // insert new element
      mTruthArray[lastindex] = element;

      // fix headers
      for (uint32_t i = dataindex + 1; i < mHeaderArray.size(); ++i) {
        auto oldindex = mHeaderArray[i];
        mHeaderArray[i] = (oldindex != (uint32_t)-1) ? oldindex + 1 : oldindex;
      }
    }
  }

  // merge another container to the back of this one
  void mergeAtBack(NA6PMCTruthContainer const& other)
  {
    const auto oldtruthsize = mTruthArray.size();
    const auto oldheadersize = mHeaderArray.size();

    // copy from other
    std::copy(other.mHeaderArray.begin(), other.mHeaderArray.end(), std::back_inserter(mHeaderArray));
    std::copy(other.mTruthArray.begin(), other.mTruthArray.end(), std::back_inserter(mTruthArray));

    // adjust information of newly attached part
    for (uint32_t i = oldheadersize; i < mHeaderArray.size(); ++i) {
      mHeaderArray[i] += oldtruthsize;
    }
  }

  // merge part of another container ("n" entries starting from "from") to the back of this one
  void mergeAtBack(NA6PMCTruthContainer const& other, size_t from, size_t n)
  {
    const auto oldtruthsize = mTruthArray.size();
    const auto oldheadersize = mHeaderArray.size();
    auto endIdx = from + n;
    assert(endIdx <= other.mHeaderArray.size());
    const auto* headBeg = &other.mHeaderArray[from];
    const auto* headEnd = headBeg + n;
    const auto* trtArrBeg = &other.mTruthArray[other.getMCTruthHeader(from)];
    const auto* trtArrEnd = (endIdx == other.mHeaderArray.size()) ? (&other.mTruthArray.back()) + 1 : &other.mTruthArray[other.getMCTruthHeader(endIdx)];

    // copy from other
    std::copy(headBeg, headEnd, std::back_inserter(mHeaderArray));
    std::copy(trtArrBeg, trtArrEnd, std::back_inserter(mTruthArray));
    long offset = long(oldtruthsize) - other.getMCTruthHeader(from);
    // adjust information of newly attached part
    for (uint32_t i = oldheadersize; i < mHeaderArray.size(); ++i) {
      mHeaderArray[i] += offset;
    }
  }
  
  /// Flatten the internal arrays to the provided container
  /// Copies the content of the two vectors of PODs to a contiguous container.
  /// The flattened data starts with a specific header @ref FlatHeader describing
  /// size and content of the two vectors within the raw buffer.
  size_t flatten_to(std::vector<char>& container) const
  {
    size_t bufferSize = sizeof(FlatHeader) + sizeof(uint32_t) * mHeaderArray.size() + sizeof(NA6PMCComposedLabel) * mTruthArray.size();
    container.resize(bufferSize);
    char* target = container.data();
    auto& flatheader = *reinterpret_cast<FlatHeader*>(target);
    target += sizeof(FlatHeader);
    flatheader.version = 1;
    flatheader.sizeofHeaderElement = sizeof(uint32_t);
    flatheader.sizeofTruthElement = sizeof(NA6PMCComposedLabel);
    flatheader.reserved = 0;
    flatheader.nofHeaderElements = mHeaderArray.size();
    flatheader.nofTruthElements = mTruthArray.size();
    size_t copySize = flatheader.sizeofHeaderElement * flatheader.nofHeaderElements;
    memcpy(target, mHeaderArray.data(), copySize);
    target += copySize;
    copySize = flatheader.sizeofTruthElement * flatheader.nofTruthElements;
    memcpy(target, mTruthArray.data(), copySize);
    return bufferSize;
  }
  
  /// Restore internal vectors from a raw buffer
  /// The two vectors are resized according to the information in the \a FlatHeader
  /// struct at the beginning of the buffer. Data is copied to the vectors.
  void restore_from(const char* buffer, size_t bufferSize)
  {
    if (buffer == nullptr || bufferSize < sizeof(FlatHeader)) {
      return;
    }
    auto* source = buffer;
    auto& flatheader = *reinterpret_cast<FlatHeader const*>(source);
    source += sizeof(FlatHeader);
    if (bufferSize < sizeof(FlatHeader) + flatheader.sizeofHeaderElement * flatheader.nofHeaderElements + flatheader.sizeofTruthElement * flatheader.nofTruthElements) {
      throw std::runtime_error("inconsistent buffer size: too small");
    }
    if (flatheader.sizeofHeaderElement != sizeof(uint32_t) || flatheader.sizeofTruthElement != sizeof(NA6PMCComposedLabel)) {
      // not yet handled
      throw std::runtime_error("member element sizes don't match");
    }
    // TODO: with a spectator memory resource the vectors can be built directly
    // over the original buffer, there is the implementation for a memory resource
    // working on a FairMQ message, here we would need two memory resources over
    // the two ranges of the input buffer
    // for now doing a copy
    mHeaderArray.resize(flatheader.nofHeaderElements);
    mTruthArray.resize(flatheader.nofTruthElements);
    size_t copySize = flatheader.sizeofHeaderElement * flatheader.nofHeaderElements;
    memcpy(mHeaderArray.data(), source, copySize);
    source += copySize;
    copySize = flatheader.sizeofTruthElement * flatheader.nofTruthElements;
    memcpy(mTruthArray.data(), source, copySize);
  }
  
  /// Inflate the object from the internal buffer
  /// The class has a specific member to store flattened data. Due to some limitations in ROOT
  /// it is more efficient to first flatten the objects to a raw buffer and empty the two vectors
  /// before serialization. This function restores the vectors from the internal raw buffer.
  /// Called from the custom streamer.
  void inflate()
  {
    if (mHeaderArray.size() > 0) {
      mStreamerData.clear();
      return;
    }
    restore_from(mStreamerData.data(), mStreamerData.size());
    mStreamerData.clear();
  }

  /// Deflate the object to the internal buffer
  /// The class has a specific member to store flattened data. Due to some limitations in ROOT
  /// it is more efficient to first flatten the objects to a raw buffer and empty the two vectors
  /// before serialization. This function stores the vectors to the internal raw buffer.
  /// Called from the custom streamer.
  void deflate()
  {
    if (mStreamerData.size() > 0) {
      clear();
      return;
    }
    mStreamerData.clear();
    flatten_to(mStreamerData);
    clear();
  }
  
  void print() const;
  std::string asString() const;
  
  ClassDefNV(NA6PMCTruthContainer, 1);
};
#endif
