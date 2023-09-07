/* Copyright (C) 2023 European Spallation Source, ERIC. See LICENSE file */
//===----------------------------------------------------------------------===//
///
/// \file
///
/// \brief Class to receive and generate Gd-GEM detector readout
/// from VMM3 ASICS via the SRS readout system
///
//===----------------------------------------------------------------------===//

#include <arpa/inet.h>
#include <cinttypes>
#include <cstdio>
#include <parser/ParserVTC.h>
#include <parser/Trace.h>
#include <string.h>

//#undef TRC_LEVEL
//#define TRC_LEVEL TRC_L_DEB
namespace Gem {

int ParserVTC::parse(uint32_t data1, uint32_t data2, struct VMM3Data *vd) {
  XTRACE(PROCESS, DEB, "data1: 0x%08x, data2: 0x%04x", data1, data2);

  /// Data
  XTRACE(PROCESS, DEB, "VTC VMM3a Data");
  vd->clockCounter = (data1 >> 6) & 0x3FFFFFF;
  vd->vmmid = (data1 >> 5) & 0x0001;
  vd->overThreshold = (data1 >> 4) & 0x0001;
  vd->bcid = (data1 & 0x000F) * 256 + ((data2 >> 24) & 0x00FF);
  vd->chno = data2 & 0x3f;
  vd->adc = (data2 >> 6) & 0x3FF;
  vd->tdc = (data2 >> 16) & 0x0FF;
  if (vd->clockCounter < last_clock_counter) {
    if (last_clock_counter - vd->clockCounter < 0x0FFFFFFF) {
      stats.ParserTimestampSeqErrors++;
      XTRACE(PROCESS, DEB,
             "ParserTimestampSeqErrors:  ts %llu, last_clock_counter %llu",
             vd->clockCounter, last_clock_counter);
    } else {
      stats.ParserTimestampOverflows++;
    }
  }
  last_clock_counter = vd->clockCounter;
  return 0;
}

int ParserVTC::receive(const char *buffer, int size) {
  int hits = 0;
  if (size < 4) {
    XTRACE(PROCESS, WAR, "Undersize data");
    stats.ParserErrorBytes += size;
    stats.ParserBadFrames++;
    return 0;
  }

  struct VTCHeader *srsHeaderPtr = (struct VTCHeader *)buffer;
  hdr.dataId = ntohl(srsHeaderPtr->dataId);
  hdr.frameCounter = ntohl(srsHeaderPtr->frameCounter);
  hdr.nHits = ntohl(srsHeaderPtr->nHits);

  /// maybe add a protocol error counter here
  if ((hdr.dataId & 0xffffffff) != 0xdadadada) {
    XTRACE(PROCESS, WAR, "Unknown data");
    stats.ParserBadFrames++;
    stats.ParserErrorBytes += size;
    return 0;
  }

  if (hdr.frameCounter > last_frame_counter) {
    if (last_frame_counter + 1 != hdr.frameCounter) {
      if (stats.ParserGoodFrames > 0) {
        stats.ParserFrameMissingErrors +=
            (hdr.frameCounter - last_frame_counter - 1);
        XTRACE(PROCESS, WAR, "ParserFrameMissingErrors: fc %u, last fc %u",
               hdr.frameCounter, last_frame_counter);
      }
    }
  } else {
    if (last_frame_counter - hdr.frameCounter > 0x0FFFFFFF) {
      stats.ParserFramecounterOverflows++;
      XTRACE(PROCESS, DEB, "ParserFramecounterOverflows: fc %d, last fc %d",
             hdr.frameCounter, last_frame_counter);
    } else {
      stats.ParserFrameSeqErrors++;
      XTRACE(PROCESS, WAR, "ParserFrameSeqErrors: fc %d, last fc %d",
             hdr.frameCounter, last_frame_counter);
    }
  }

  last_frame_counter = hdr.frameCounter;

  if (size < VTCHeaderSize + HitSize) {
    XTRACE(PROCESS, WAR, "Undersize data");
    stats.ParserBadFrames++;
    stats.ParserErrorBytes += size;
    return 0;
  }

  auto datalen = size - VTCHeaderSize;
  if ((datalen % 8) != 0) {
    XTRACE(PROCESS, WAR, "Invalid data length: %d", datalen);
    stats.ParserBadFrames++;
    stats.ParserErrorBytes += size;
    return 0;
  }

  int dataIndex = 0;
  int readoutIndex = 0;
  while (datalen >= HitSize) {
    XTRACE(PROCESS, DEB, "readoutIndex: %d, datalen %d", readoutIndex, datalen);
    auto Data1Offset = VTCHeaderSize + HitSize * readoutIndex;
    auto Data2Offset = Data1Offset + Data1Size;
    uint32_t data1 = htonl(*(uint32_t *)&buffer[Data1Offset]);
    uint32_t data2 = htonl(*(uint32_t *)&buffer[Data2Offset]);
    int res = parse(data1, data2, &data[dataIndex]);
    hits++;
    stats.ParserData++;
    dataIndex++;

    stats.ParserReadouts++;
    readoutIndex++;

    datalen -= HitSize;
    if (hits == maxHits && datalen > 0) {
      XTRACE(PROCESS, WAR, "Data overflow, skipping %d bytes", datalen);
      stats.ParserErrorBytes += datalen;
      break;
    }
  }

  if (hits != hdr.nHits) {
    XTRACE(PROCESS, WAR,
           "Number of hits wrong, in header %d hits, in data %d hits",
           hdr.nHits, hits);
    stats.ParserBadFrames++;
  } else {
    stats.ParserGoodFrames++;
  }

  return hits;
}

} // namespace Gem
