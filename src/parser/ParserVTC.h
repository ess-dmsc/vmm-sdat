/* Copyright (C) 2023 European Spallation Source, ERIC. See LICENSE file */
//===----------------------------------------------------------------------===//
///
/// \file
///
/// \brief Class to receive data
/// from VMM3 ASICs via the VTC readout system
///
//===----------------------------------------------------------------------===//

#pragma once

#include <cinttypes>
#include <parser/BitMath.h>
#include <parser/NMXStats.h>
#include <parser/SRSTime.h>
#include <string.h>

static const int MaxHybrids{1}; /// Maximum number of VMMs per VTC card
static const int MaxVTCs{1};    /// Maximum number of VTCs

namespace Gem {

class ParserVTC {
public:
  // bytes
  static const int VTCHeaderSize{12};
  static const int HitSize{8};
  static const int Data1Size{4};

  ///< Do NOT rearrange fields, used for casting to data pointer
  struct VTCHeader {
    uint32_t dataId;       /// data type identifier
    uint32_t frameCounter; /// frame counter packet field
    uint32_t nHits; /// indicates if marker for vmm was sent in last frame
  };

  /// Data related to a single Hit
  struct VMM3Data {
    uint64_t
        clockCounter; /// 47 bits can change within a packet so must be here
    uint16_t bcid;    /// 12 bit - bcid after graydecode
    uint16_t adc;     /// 10 bit - adc value from vmm readout
    uint8_t tdc;      ///  8 bit - tdc value from vmm readout
    uint8_t chno;     ///  6 bit - channel number from readout
    uint8_t
        overThreshold; ///  1 bit - over threshold flag for channel from readout
    uint8_t vmmid;     ///  5 bit - asic identifier - unique id per fec 0 - 15
  };

  /// \brief create a data handler for VMM3 SRS data of fixed size Capacity
  /// \param maxelements The maximum number of readout elements
  ParserVTC(int maxelements, NMXStats &stats, SRSTime time_intepreter)
      : maxHits(maxelements), stats(stats), srsTime(time_intepreter) {
    data = new struct VMM3Data[maxHits];
  }

  /// Delete allocated data, set pointers to nullptr
  ~ParserVTC() {
    delete[] data;
    data = nullptr;
  }

  // \todo use Buffer<char>
  /// \brief reveive readouts from a binary payload buffer, return number of
  /// data elements
  int receive(const char *buffer, int size);

  /// \brief parse the readouts into a data array
  /// \param data1 the raw (unbitreversed) data1 field of a SRS packet
  /// \param data2 the raw (unbitreversed) data2 field of a SRS packet
  /// \param vmd VMM2Data structure holding the parsed data (tdc, bcid, adc,
  /// ...)
  int parse(uint32_t data1, uint32_t data2, struct VMM3Data *vmd);

  /// Holds data common to all readouts in a packet
  struct VTCHeader hdr;

  /// holds all readout data in a packet (up to max_elems)
  struct VMM3Data *data{nullptr};

  int maxHits{0}; /// Maximum capacity of data array

  NMXStats &stats;
  SRSTime srsTime;
  uint64_t last_clock_counter = 0;
  uint32_t last_frame_counter = 0;
};
} // namespace Gem