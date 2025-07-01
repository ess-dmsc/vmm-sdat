/* Copyright (C) 2018-2021 European Spallation Source, ERIC. See LICENSE file */
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
#include <string.h>
#include <parser/ParserSRS.h>
#include <parser/Trace.h>
#include <iostream>

// #undef TRC_LEVEL
// #define TRC_LEVEL TRC_L_DEB
namespace Gem {

int ParserSRS::parse(uint32_t data1, uint16_t data2, struct VMM3Data *vmm3Data) {
  XTRACE(PROCESS, DEB, "data1: 0x%08x, data2: 0x%04x", data1, data2);
  int dataflag = (data2 >> 15) & 0x1;

  if (dataflag) {
    /// Data
    XTRACE(PROCESS, DEB, "SRS Data");

    vmm3Data->overThreshold = (data2 >> 14) & 0x01;
    stats.ParserOverThreshold += vmm3Data->overThreshold;
    vmm3Data->chno = (data2 >> 8) & 0x3f;
    vmm3Data->tdc = data2 & 0xff;
    vmm3Data->vmmid = (data1 >> 22) & 0x1F;
    vmm3Data->timestampOffset = (data1 >> 27) & 0x1F;
    vmm3Data->adc = (data1 >> 12) & 0x3FF;
    vmm3Data->bcid = BitMath::gray2bin32(data1 & 0xFFF);
    uint16_t idx = (pd.fecId - 1) * MaxVMMs + vmm3Data->vmmid;
    if(markers[idx].fecTimeStamp > 0)  {
      vmm3Data->fecTimeStamp = markers[idx].fecTimeStamp;
      vmm3Data->triggerTime = markers[idx].triggerTime;
      vmm3Data->triggerCounter = markers[idx].triggerCounter + 1;
      //std::cout << "vmm3Data->fecTimeStamp  " <<  idx << " " << vmm3Data->fecTimeStamp  << std::endl;
      //std::cout << "vmm3Data->triggerTime " <<  idx << " " <<  vmm3Data->triggerTime << std::endl;
      //std::cout << "vmm3Data->triggerCounter " <<  idx << " " <<  vmm3Data->triggerCounter << std::endl;
 
    }
   
    XTRACE(PROCESS, DEB, "SRS Data: vmm: %d, channel: %d. adc: %d",
      vmm3Data->vmmid, vmm3Data->chno, vmm3Data->adc);
    return 1;
  } else {
    /// Marker
    uint8_t vmmid = (data2 >> 10) & 0x1F;
    uint16_t idx = (pd.fecId - 1) * MaxVMMs + (vmmid%16);
    
    if(vmmid >= 16) {
      if(dataFormat == "TRG") {
        int triggerFlag = (data1 >> 28) & 0x0F;
        if(triggerFlag == 0xF) {
          uint64_t event_counter_high = data1 & 0x03F;
          uint64_t event_counter_low = data2 & 0x3FF;
          markers[idx].triggerCounter = event_counter_high*1024+event_counter_low;
        } 
        else {
          uint64_t timestamp_lower_10bit = data2 & 0x03FF;
          uint64_t timestamp_upper_32bit = data1;
          uint64_t timestamp_42bit = (timestamp_upper_32bit << 10);
          markers[idx].triggerTime = timestamp_42bit;   
        }
      }
    }
    else {
      uint64_t timestamp_lower_10bit = data2 & 0x03FF;
      uint64_t timestamp_upper_32bit = data1;
      uint64_t timestamp_42bit = (timestamp_upper_32bit << 10)
        + timestamp_lower_10bit;
      //normal data marker  
      if(dataFormat == "SRS") {
         if(markers[idx].fecTimeStamp > timestamp_42bit) {
            if (markers[idx].fecTimeStamp < 0x1FFFFFFF + timestamp_42bit) {
              stats.ParserTimestampSeqErrors++;
              XTRACE(PROCESS, DEB, "ParserTimestampSeqErrors:  ts %llu, marker ts %llu", timestamp_42bit, markers[idx].fecTimeStamp);
            }
            else {
              stats.ParserTimestampOverflows++;
            }
        }
        markers[idx].fecTimeStamp = timestamp_42bit;

   
      }
      // relative trigger time stamp
      else if(timestamp_42bit < 4096) {
        markers[idx].fecTimeStamp = timestamp_42bit;   
      }
    }
    return 0;
  }
}


int ParserSRS::receive(const char *buffer, int size) {
  int hits = 0;
  if (size < 4) {
    XTRACE(PROCESS, WAR, "Undersize data");
    stats.ParserErrorBytes += size;
    stats.ParserBadFrames++;
    return 0;
  }

  struct SRSHeader *srsHeaderPtr = (struct SRSHeader *) buffer;
  hdr.frameCounter = ntohl(srsHeaderPtr->frameCounter);

  if (pd.nextFrameCounter != hdr.frameCounter) {
    if(hdr.frameCounter > pd.nextFrameCounter) {
      if(stats.ParserGoodFrames > 0) {
        stats.ParserFrameMissingErrors +=
         (hdr.frameCounter - pd.nextFrameCounter);
        XTRACE(PROCESS, WAR, "ParserFrameMissingErrors: fc %d, next fc %d",
        hdr.frameCounter, pd.nextFrameCounter);
      }
    }
    else {
      if (pd.nextFrameCounter - hdr.frameCounter > 0x0FFFFFFF) {
        stats.ParserFramecounterOverflows++;
        XTRACE(PROCESS, DEB, "ParserFramecounterOverflows: fc %d, next fc %d",
          hdr.frameCounter, pd.nextFrameCounter);
      }
      else {
        stats.ParserFrameSeqErrors++;
        XTRACE(PROCESS, WAR, "ParserFrameSeqErrors: fc %d, next fc %d",
          hdr.frameCounter, pd.nextFrameCounter);
        /// \todo here we used to clear markers. If this is needed again
        /// some day a new method is required says Dorothea P.
      }
    }
  }
  else {
    if(hdr.frameCounter == 0) {
      stats.ParserFramecounterOverflows++;
    }

  }
  pd.nextFrameCounter = hdr.frameCounter + 1;

  if (size < SRSHeaderSize + HitAndMarkerSize) {
    XTRACE(PROCESS, WAR, "Undersize data");
    stats.ParserBadFrames++;
    stats.ParserErrorBytes += size;
    return 0;
  }

  hdr.dataId = ntohl(srsHeaderPtr->dataId);
  /// maybe add a protocol error counter here
  if ((hdr.dataId & 0xffffff00) != 0x564d3300) {
    XTRACE(PROCESS, WAR, "Unknown data");
    stats.ParserBadFrames++;
    stats.ParserErrorBytes += size;
    return 0;
  }

  pd.fecId = (hdr.dataId >> 4) & 0x0f;

  if (pd.fecId == 0) {
    XTRACE(PROCESS, WAR, "Invalid fecId: %u", pd.fecId);
    stats.ParserBadFrames++;
    stats.ParserErrorBytes += size;
    return 0;
  }
  XTRACE(PROCESS, DEB, "fecId: %u", pd.fecId);
  hdr.udpTimeStamp = ntohl(srsHeaderPtr->udpTimeStamp);
  //This header component will vanish soon
  //and be replaced by a timestamp for each vmm
  hdr.offsetOverflow = ntohl(srsHeaderPtr->offsetOverflow);

  auto datalen = size - SRSHeaderSize;
  if ((datalen % 6) != 0) {
    XTRACE(PROCESS, WAR, "Invalid data length: %d", datalen);
    stats.ParserBadFrames++;
    stats.ParserErrorBytes += size;
    return 0;
  }

  int dataIndex = 0;
  int readoutIndex = 0;
  while (datalen >= HitAndMarkerSize) {
    XTRACE(PROCESS, DEB, "readoutIndex: %d, datalen %d",
    		readoutIndex, datalen);
    auto Data1Offset = SRSHeaderSize + HitAndMarkerSize * readoutIndex;
    auto Data2Offset = Data1Offset + Data1Size;
    uint32_t data1 = htonl(*(uint32_t *) &buffer[Data1Offset]);
    uint16_t data2 = htons(*(uint16_t *) &buffer[Data2Offset]);

    int res = parse(data1, data2, &data[dataIndex]);
    if (res == 1) { // This was data
      hits++;
      stats.ParserData++;
      dataIndex++;
    } else {
      stats.ParserMarkers++;
    }
    stats.ParserReadouts++;
    readoutIndex++;

    datalen -= 6;
    if (hits == maxHits && datalen > 0) {
      XTRACE(PROCESS, WAR, "Data overflow, skipping %d bytes", datalen);
      stats.ParserErrorBytes += datalen;
      break;
    }
  }
  stats.ParserGoodFrames++;

  return hits;
}

}


