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
#include <iostream>
#include "log.h"
#include <parser/ParserSRS.h>

//Maxi-ROC parse function
int ParserSRS::parse(uint64_t data, struct VMM3Data *vmm3Data) {

  uint8_t flag = (data >> 60) & 0x0F;
  if (flag > 0) {    /// Data
   
    vmm3Data->tdc = data & 0xFF;
    vmm3Data->chno = (data >> 8) & 0x3F;
    vmm3Data->overThreshold = (data >> 14) & 0x01;
    vmm3Data->adc = (data >> 15) & 0x3FF;
    vmm3Data->bcid = (data >> 25) & 0xFFF;
    vmm3Data->timestampOffset = (data >> 37) & 0x7FFF;
    vmm3Data->vmmid = (data >> 52) & 0xFF;
    stats.ParserOverThreshold += vmm3Data->overThreshold;
    uint16_t idx = (pd.fecId - 1) * MaxVMMsMaxi + vmm3Data->vmmid%MaxVMMsMaxi;
    if(markers[idx].fecTimeStamp > 0)  {
      vmm3Data->fecTimeStamp = markers[idx].fecTimeStamp;
    }
    return 1;
  } else {
    /// Marker
    uint8_t vmmid = (data >> 52) & 0xFF;
    uint16_t idx = (pd.fecId - 1) * MaxVMMsMaxi + (vmmid%MaxVMMsMaxi);
    uint64_t timestamp_52bit = data & 0xFFFFFFFFFFFFF;
    if(markers[idx].fecTimeStamp > timestamp_52bit) {
      if (markers[idx].fecTimeStamp < 0x1FFFFFFF + timestamp_52bit) {
        stats.ParserTimestampSeqErrors++;
      }
      else {
        stats.ParserTimestampOverflows++;
      }
    }
    markers[idx].fecTimeStamp = timestamp_52bit;

    return 0;
  }
}

//SRS parse function
int ParserSRS::parse(uint32_t data1, uint16_t data2, struct VMM3Data *vmm3Data) {
  int dataflag = (data2 >> 15) & 0x1;

  if (dataflag) {
    /// Data

    vmm3Data->overThreshold = (data2 >> 14) & 0x01;
    stats.ParserOverThreshold += vmm3Data->overThreshold;
    vmm3Data->chno = (data2 >> 8) & 0x3f;
    vmm3Data->tdc = data2 & 0xff;
    vmm3Data->vmmid = (data1 >> 22) & 0x1F;
    vmm3Data->timestampOffset = (data1 >> 27) & 0x1F;
    vmm3Data->adc = (data1 >> 12) & 0x3FF;
    vmm3Data->bcid = BitMath::gray2bin32(data1 & 0xFFF);
    uint16_t idx = (pd.fecId - 1) * MaxVMMsSRS + vmm3Data->vmmid%MaxVMMsSRS;
    if(markers[idx].fecTimeStamp > 0)  {
      vmm3Data->fecTimeStamp = markers[idx].fecTimeStamp;
      vmm3Data->triggerTime = markers[idx].triggerTime;
      vmm3Data->triggerCounter = markers[idx].triggerCounter + 1;
    }
    return 1;
  } else {
    /// Marker
    uint8_t vmmid = (data2 >> 10) & 0x1F;
    uint16_t idx = (pd.fecId - 1) * MaxVMMsSRS + (vmmid%MaxVMMsSRS);

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
      //In normal SRS mode a marker with vmmid 31 contains the NIM trigger timestamp
      else {
        uint64_t timestamp_lower_10bit = data2 & 0x03FF;
        uint64_t timestamp_upper_32bit = data1;
        uint64_t timestamp_42bit = (timestamp_upper_32bit << 10);
        //Since the data does not always come out time ordered, we store the last 5 timestamps
        //Each new timestamp is added in position 0, and the last one at position 4 is delected
        for(int n=4; n>=1;n--) {
          nim->TriggerTime[n] = nim->TriggerTime[n-1];
        }
        nim->TriggerTime[0]  = timestamp_42bit;
      }
    }
    else {
      uint64_t timestamp_lower_10bit = data2 & 0x03FF;
      uint64_t timestamp_upper_32bit = data1;
      uint64_t timestamp_42bit = (timestamp_upper_32bit << 10)
        + timestamp_lower_10bit;
      //normal data marker  
      if(dataFormat == "SRS" || dataFormat == "srs") {
         if(markers[idx].fecTimeStamp > timestamp_42bit) {
            if (markers[idx].fecTimeStamp < 0x1FFFFFFF + timestamp_42bit) {
              stats.ParserTimestampSeqErrors++;
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
  if (dataFormat == "SRS" or dataFormat == "TRG") {
    if(size < 16) {
      stats.ParserErrorBytes += size;
      stats.ParserBadFrames++;
      return 0;
    }

    struct SRSPacketHeader *srsHeaderPtr = (struct SRSPacketHeader *) buffer;   
    hdr.frameCounter = ntohl(srsHeaderPtr->frameCounter);
    if (pd.nextFrameCounter != hdr.frameCounter) {
      if(hdr.frameCounter > pd.nextFrameCounter) {
        if(stats.ParserGoodFrames > 0) {
          stats.ParserFrameMissingErrors +=
          (hdr.frameCounter - pd.nextFrameCounter);
        }
      }
      else {
        if (pd.nextFrameCounter - hdr.frameCounter > 0x0FFFFFFF) {
          stats.ParserFramecounterOverflows++;
        }
        else {
          stats.ParserFrameSeqErrors++;
        }
      }
    }
    else {
      if(hdr.frameCounter == 0) {
        stats.ParserFramecounterOverflows++;
      }

    }
    pd.nextFrameCounter = hdr.frameCounter + 1;

    if (size < SRSHeaderSize + SRSHitAndMarkerSize) {
      stats.ParserBadFrames++;
      stats.ParserErrorBytes += size;
      return 0;
    }

    hdr.dataId = ntohl(srsHeaderPtr->dataId);
    if ((hdr.dataId & 0xffffff00) != 0x564d3300) {
      stats.ParserBadFrames++;
      stats.ParserErrorBytes += size;
      return 0;
    }

    pd.fecId = (hdr.dataId >> 4) & 0x0f;

    if (pd.fecId == 0) {
      stats.ParserBadFrames++;
      stats.ParserErrorBytes += size;
      return 0;
    }
    hdr.udpTimeStamp = ntohl(srsHeaderPtr->udpTimeStamp);

    int dataIndex = 0;
    int readoutIndex = 0;
    auto datalen = size - SRSHeaderSize;
    while (datalen >= SRSHitAndMarkerSize) {
      auto Data1Offset = SRSHeaderSize + SRSHitAndMarkerSize * readoutIndex;
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

      datalen -= SRSHitAndMarkerSize;
      if (hits == maxHits && datalen > 0) {
        stats.ParserErrorBytes += datalen;
        break;
      }
    }
    stats.ParserGoodFrames++;

    return hits;
  }
  else if (dataFormat == "MAX") {
    if(size < 32) {
      stats.ParserErrorBytes += size;
      stats.ParserBadFrames++;
      return 0;
    }
    struct MaxiPacketHeader *maxiHeaderPtr = (struct MaxiPacketHeader *) buffer;    
    uint32_t tmp1 = ntohl(maxiHeaderPtr->frameCounter1);
    uint32_t tmp2 = ntohl(maxiHeaderPtr->frameCounter2);
    hdr.frameCounter = (static_cast<uint64_t>(tmp1) << 32) + static_cast<uint64_t>(tmp2);

    if (pd.nextFrameCounter != hdr.frameCounter) {
      if(hdr.frameCounter > pd.nextFrameCounter) {
        if(stats.ParserGoodFrames > 0) {
          stats.ParserFrameMissingErrors +=
          (hdr.frameCounter - pd.nextFrameCounter);
        }
      }
      else {
        
        if (pd.nextFrameCounter - hdr.frameCounter > 0x0FFFFFFF) {
          stats.ParserFramecounterOverflows++;
        }
        else {
          stats.ParserFrameSeqErrors++;
        }
      }
    }
    else {
      if(hdr.frameCounter == 0) {
        stats.ParserFramecounterOverflows++;
      }
    }

    pd.nextFrameCounter = hdr.frameCounter + 1;

    if (size < MaxiHeaderSize + MaxiHitAndMarkerSize) {
      stats.ParserBadFrames++;
      stats.ParserErrorBytes += size;
      return 0;
    }

    tmp1 = ntohl(maxiHeaderPtr->dataId1);
    tmp2 = ntohl(maxiHeaderPtr->dataId2);
    hdr.dataId = (static_cast<uint64_t>(tmp1) << 32) + static_cast<uint64_t>(tmp2);

    if ((hdr.dataId & 0xffffffffffffff00) != 0x4D415849564D4D00) {
      stats.ParserBadFrames++;
      stats.ParserErrorBytes += size;
      return 0;
    }
    tmp1 = ntohl(maxiHeaderPtr->block_port);
    hdr.block = static_cast<uint16_t>(tmp1>>16);
    hdr.udpPort = static_cast<uint16_t>(tmp1&0xFF);

    tmp2 = ntohl(maxiHeaderPtr->ipAddress);
    hdr.ipAddress = static_cast<uint64_t>(tmp2);
    //MAXI-ROC IP last octect is FEC ID
    pd.fecId = hdr.ipAddress & 0xFF;

    if (pd.fecId == 0) {
      stats.ParserBadFrames++;
      stats.ParserErrorBytes += size;
      return 0;
    }

    tmp1 = ntohl(maxiHeaderPtr->udpTimeStamp1);
    tmp2 = ntohl(maxiHeaderPtr->udpTimeStamp2);
    hdr.udpTimeStamp = (static_cast<uint64_t>(tmp1) << 32) + static_cast<uint64_t>(tmp2);

    auto datalen = size - MaxiHeaderSize;
    int dataIndex = 0;
    int readoutIndex = 0;
    while (datalen >= MaxiHitAndMarkerSize) {
      auto Data1Offset = MaxiHeaderSize + MaxiHitAndMarkerSize * readoutIndex;
      auto Data2Offset = Data1Offset + Data1Size;
      uint32_t data1 = htonl(*(uint32_t *) &buffer[Data1Offset]);
      uint32_t data2 = htonl(*(uint32_t *) &buffer[Data2Offset]);
      uint64_t data64 = (static_cast<uint64_t>(data1) << 32) + static_cast<uint64_t>(data2);
      int res = parse(data64, &data[dataIndex]);
      if (res == 1) { // This was data
        hits++;
        stats.ParserData++;
        dataIndex++;
      } else {
        stats.ParserMarkers++;
      }
      stats.ParserReadouts++;
      readoutIndex++;

      datalen -= MaxiHitAndMarkerSize;
      if (hits == maxHits && datalen > 0) {
        stats.ParserErrorBytes += datalen;
        break;
      }
    }
    stats.ParserGoodFrames++;

    return hits;
  }
  else {
    return 0;
  }
}

