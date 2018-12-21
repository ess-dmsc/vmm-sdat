/// Copyright (C) 2016-2018 European Spallation Source, ERIC. See LICENSE file
//===----------------------------------------------------------------------===//
///
/// \file
///
//===----------------------------------------------------------------------===//

#pragma once

#include <vector>
#include <memory>
#include <iostream>
#include <h5cpp/hdf5.hpp>

#define STRINGIFY(x) #x

#define H5_COMPOUND_DEFINE_TYPE(x) using Type = x; using TypeClass = Compound; static TypeClass create(const Type & = Type())
#define H5_COMPOUND_INIT auto type = datatype::Compound::create(sizeof(Type))
#define H5_COMPOUND_INSERT_MEMBER(member) type.insert(STRINGIFY(member), offsetof(Type, member), datatype::create<decltype(Type::member)>())
#define H5_COMPOUND_RETURN return type


//struct __attribute__ ((packed)) Readout {
struct Readout {
  /// \todo use constexpr string_view when c++17 arrives
  static std::string DatasetName() { return "srs_hits"; }
  static uint16_t FormatVersion() { return 1; }

  /// \todo consider reordering these to optimize
  /// !!! DO NOT MODIFY BELOW - READ HEADER FIRST !!!
  uint8_t fec{0};
  uint8_t chip_id{0};
  uint64_t srs_timestamp{0};
  uint16_t channel{0};
  uint16_t bcid{0};
  uint16_t tdc{0};
  uint16_t adc{0};
  bool over_threshold{false};
  float chiptime{0.0};
  /// !!! DO NOT MODIFY ABOVE -- READ HEADER FIRST !!!

  bool operator==(const Readout &other) const {
    return (
        (fec == other.fec) &&
            (chip_id == other.chip_id) &&
            (srs_timestamp == other.srs_timestamp) &&
            (channel == other.channel) &&
            (bcid == other.bcid) &&
            (tdc == other.tdc) && (adc == other.adc) &&
            (over_threshold == other.over_threshold) &&
            (chiptime == other.chiptime)
    );
  }
};

namespace hdf5 {

namespace datatype {
template<>
class TypeTrait<Readout> {
public:
  H5_COMPOUND_DEFINE_TYPE(Readout) {
    H5_COMPOUND_INIT;
    /// Make sure ALL member variables are inserted
    
    H5_COMPOUND_INSERT_MEMBER(fec);
    
    H5_COMPOUND_INSERT_MEMBER(chip_id);
    H5_COMPOUND_INSERT_MEMBER(srs_timestamp);
    H5_COMPOUND_INSERT_MEMBER(channel);
    H5_COMPOUND_INSERT_MEMBER(bcid);
    H5_COMPOUND_INSERT_MEMBER(tdc);
    H5_COMPOUND_INSERT_MEMBER(adc);
    H5_COMPOUND_INSERT_MEMBER(over_threshold);
    H5_COMPOUND_INSERT_MEMBER(chiptime);
    
    H5_COMPOUND_RETURN;
  }
};
}

}
