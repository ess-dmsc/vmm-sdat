/* Copyright (C) 2018 European Spallation Source, ERIC. See LICENSE file */
//===----------------------------------------------------------------------===//
///
/// \file
///
/// \brief Class to receive and generate Gd-GEM detector readout
/// from VMM3 ASICS via the SRS readout system
///
//===----------------------------------------------------------------------===//

#pragma once

#include <cinttypes>
#include <string>
#include <parser/BitMath.h>


static const int MaxVMMs{16}; ///Maximum number of VMMs per FEC card
static const int MaxFECs{16}; ///Maximum number of FECs per EFU

namespace Gem {

  struct ParserStats {
  // Input thread
  int64_t RxPackets{0};
  int64_t RxBytes{0};
  int64_t RxIdle{0};
  int64_t FifoPushErrors{0};
  int64_t PaddingFor64ByteAlignment[4]; // cppcheck-suppress unusedStructMember

  // Processing thread
  int64_t ProcessingIdle {0};
  int64_t FifoSeqErrors{0};

  // Parser stats
  int64_t ParserFrameSeqErrors{0};
  int64_t ParserFrameMissingErrors{0};
  int64_t ParserFramecounterOverflows{0};
  int64_t ParserTimestampLostErrors{0};
  int64_t ParserTimestampSeqErrors{0};
  int64_t ParserTimestampOverflows{0};
  int64_t ParserBadFrames{0};
  int64_t ParserGoodFrames{0};
  int64_t ParserErrorBytes{0};
  int64_t ParserMarkers{0};
  int64_t ParserData{0};
  int64_t ParserReadouts{0};
  int64_t ParserOverThreshold{0};
  };

class ParserSRS {
public:
  // bytes
  static const int SRSHeaderSize{16};
  static const int HitAndMarkerSize{6};
  static const int Data1Size{4};

  ///< Do NOT rearrange fields, used for casting to data pointer
  struct SRSHeader {
    uint32_t frameCounter;   /// frame counter packet field
    uint32_t dataId;         /// data type identifier packet field + ID of the FEC card (0-255)
    uint32_t udpTimeStamp;   /// Transmission time for UDP packet
    uint32_t offsetOverflow; /// indicates if marker for vmm was sent in last frame
  };

  /// 
  struct VMM3Marker {
    uint64_t fecTimeStamp{0};   /// 42 bit or 12 bit if it is relative trigger time
    uint64_t triggerTime{0};   /// 42 bit
    uint64_t triggerCounter{0};   /// 42 bit
  };

  /// Data common to all hits and markers, or other parser related data
  struct ParserData {
    uint8_t fecId{1};
    uint32_t nextFrameCounter{0};
  };

  /// Data related to a single Hit
  struct VMM3Data {
    uint64_t fecTimeStamp; /// 42 bits can change within a packet so must be here
    uint16_t bcid;         /// 12 bit - bcid after graydecode
    uint16_t adc;          /// 10 bit - adc value from vmm readout
    uint8_t tdc;           ///  8 bit - tdc value from vmm readout
    uint8_t chno;          ///  6 bit - channel number from readout
    uint8_t overThreshold; ///  1 bit - over threshold flag for channel from readout
    uint8_t vmmid;         ///  5 bit - asic identifier - unique id per fec 0 - 15
    uint8_t timestampOffset; ///  5 bit
    uint64_t triggerTime; /// Time of the occurence of the FEC trigger signa
    uint64_t triggerCounter; /// Trigger number
  };

  /// \brief create a data handler for VMM3 SRS data of fixed size Capacity
  /// \param maxelements The maximum number of readout elements
  ParserSRS(int maxelements, ParserStats & stats, std::string format) : maxHits(maxelements), 
  stats(stats),  dataFormat(format) {
    markers = new struct VMM3Marker[MaxFECs * MaxVMMs];
    data = new struct VMM3Data[maxHits];
  }

  /// Delete allocated data, set pointers to nullptr
  ~ParserSRS() {
    delete[] data;
    data = nullptr;
    delete[] markers;
    markers = nullptr;
  }

  // \todo use Buffer<char>
  /// \brief reveive readouts from a binary payload buffer, return number of
  /// data elements
  int receive(const char *buffer, int size);

  /// \brief parse the readouts into a data array
  /// \param data1 the raw (unbitreversed) data1 field of a SRS packet
  /// \param data2 the raw (unbitreversed) data2 field of a SRS packet
  /// \param vmm3Data VMM3Data structure holding the parsed data (tdc, bcid, adc, ...)
  int parse(uint32_t data1, uint16_t data2, struct VMM3Data *vmm3Data);

  /// Holds data common to all readouts in a packet
  struct SRSHeader hdr;

  /// holds all readout data in a packet (up to max_elems)
  struct VMM3Data *data{nullptr};

  /// See description above
  struct ParserData pd;

  /// holds time bases for all vmms in a readout
  struct VMM3Marker *markers{nullptr};

  int maxHits{0};       /// Maximum capacity of data array

  ParserStats & stats;
  std::string dataFormat;

 };
}