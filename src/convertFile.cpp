/***************************************************************************
**  vmm-sdat
**  Data analysis program for VMM3a data
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see http://www.gnu.org/licenses/.
**
****************************************************************************
**  Contact: dorothea.pfeiffer@cern.ch
**  Date: 12.10.2025
**  Version: 1.0.0
****************************************************************************
**
**  vmm-sdat
**  convertFile.cpp
**
****************************************************************************/

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "log.h"
#include "Clusterer.h"
#include "Configuration.h"
#include <parser/CalibrationFile.h>
#include <parser/ParserSRS.h>
#include <parser/ReaderPcap.h>


int main(int argc, char ** argv) {
  corryvreckan::Log::addStream(std::cout);
  uint64_t total_hits = 0;
  uint64_t triggerSignals = 0;
  std::chrono::time_point < std::chrono::system_clock > timeEnd, timeStart;

  Configuration m_config;
  Statistics m_stats;

  if (!m_config.ParseCommandLine(argc, argv)) {
    return -1;
  }
  if (!m_config.CreateMapping()) {
    return -1;
  }
  if (m_config.pDataFormat != "SRS" && m_config.pDataFormat != "TRG") {
    return -1;
  }

  corryvreckan::LogLevel log_level = corryvreckan::Log::getLevelFromString(m_config.pLogLevel);
  corryvreckan::Log::setReportingLevel(log_level);


  timeStart = std::chrono::system_clock::now();
  uint64_t last_time = 0;

  Clusterer * m_Clusterer = new Clusterer(m_config, m_stats);
  m_stats.CreateFECStats(m_config);
  m_stats.CreateClusterStats(m_config);
  m_stats.CreatePCAPStats(m_config);

  double firstTime = 0;
  char buffer[10000];

  ParserStats stats;
  ParserSRS * parser = new ParserSRS(2000, stats, m_config.pDataFormat);
  CalibrationFile calfile(m_config.pCalFilename);
  ReaderPcap pcap(m_config.pFileName);
  int ret = pcap.open();
  if (ret < 0) {
    corryvreckan::Log::setSection("convertFile");
    LOG(FATAL) << "Error opening file: " << m_config.pFileName <<
      ": return value " << ret;
    return -1;
  }

  uint64_t pcappackets = 0;
  int rdsize;
  bool doContinue = true;
  while (doContinue &&
    (rdsize = pcap.read((char * ) & buffer, sizeof(buffer))) != -1) {
    if (rdsize == 0) {
      continue; // non udp data
    }
    int hits = parser -> receive(buffer, rdsize);
    total_hits += hits;
    for (int i = 0; i < hits; i++) {
      auto & d = parser -> data[i];
      double timestampOffset = -1.0;

      // triggerOffset goes from -1 to 15
      // but presented as uint8_t
      // latency violation
      if (d.timestampOffset <= 15.0) {
        timestampOffset = static_cast < double > (d.timestampOffset);
      } else if (d.timestampOffset == 31) {
        timestampOffset = -1.0;
      } else if (d.timestampOffset == 16) {
        timestampOffset = -99.0;
      }
      if (timestampOffset == -99) {
        continue;
      }
      double srs_timestamp = 0;
      int event_counter = d.triggerCounter;
      if (m_config.pDataFormat == "SRS" || m_config.pDataFormat == "srs") {
        srs_timestamp =
          (static_cast < double > (d.fecTimeStamp) * m_config.pBCTime_ns +
            m_config.pOffsetPeriod * timestampOffset);
      } else {
        srs_timestamp =
          (static_cast < double > (d.triggerTime) * m_config.pBCTime_ns +
            m_config.pOffsetPeriod * timestampOffset -
            static_cast < double > (d.fecTimeStamp) * m_config.pBCTime_ns);
      }
      if (firstTime == 0) {
        firstTime = srs_timestamp;
      }
 
      double t0_correction = 0;
      std::pair < uint8_t, uint8_t > fec_vmm =
        std::make_pair(parser->pd.fecId, d.vmmid);
      auto searchMap = m_config.pFecVMM_time0.find(fec_vmm);
      if (searchMap != m_config.pFecVMM_time0.end()) {
        std::string t0 = m_config.pFecVMM_time0[fec_vmm];
        if (t0 == "run") {
          t0_correction = firstTime;
        } else {
          t0_correction = std::stod(t0);
        }
      }

      srs_timestamp = srs_timestamp - t0_correction;
      auto calib =
        calfile.getCalibration(parser -> pd.fecId, d.vmmid, d.chno);

      double chiptime_corrected =
        static_cast < double > (d.bcid) * m_config.pBCTime_ns +
        (1.5 * m_config.pBCTime_ns -
          static_cast < double > (d.tdc) *
          static_cast < double > (m_config.pTAC) / 255.0 -
          calib.time_offset) *
        calib.time_slope;

      uint16_t corrected_adc = static_cast < uint16_t > (
        (static_cast < double > (d.adc) - calib.adc_offset) *
        calib.adc_slope);

      if (corrected_adc > 1023) {
        corryvreckan::Log::setSection("convertFile");
        LOG(DEBUG) << "After correction, ADC value larger than 1023 "
        "(10bit)!\nUncorrected ADC value " <<
        d.adc << ", corrected adc value: " << corrected_adc;

        corrected_adc = 1023;
      } else if (corrected_adc < 0) {
        corryvreckan::Log::setSection("convertFile");
        LOG(DEBUG) << "After correction, ADC value smaller than 0!"
        "\nUncorrected ADC value " <<
        d.adc << ", corrected adc value: " << corrected_adc;
        corrected_adc = 0;
      }

      double time_without_calib =
        (1.5 * m_config.pBCTime_ns -
          static_cast < double > (d.tdc) *
          static_cast < double > (m_config.pTAC) / 255.0);

      double time_with_calib =
        (1.5 * m_config.pBCTime_ns -
          static_cast < double > (d.tdc) *
          static_cast < double > (m_config.pTAC) / 255.0 -
          calib.time_offset) *
        calib.time_slope;

      m_Clusterer -> FillCalibHistos(parser -> pd.fecId, d.vmmid, d.chno,
        d.adc, corrected_adc,
        time_without_calib, time_with_calib);

      double timewalk_correction =
        calib.timewalk_d +
        (calib.timewalk_a - calib.timewalk_d) /
        (1 +
          pow(corrected_adc / calib.timewalk_c, calib.timewalk_b));

      double corrected_time = chiptime_corrected - timewalk_correction;

      bool result = m_Clusterer -> AnalyzeHits(
        srs_timestamp, parser -> pd.fecId, d.vmmid, d.chno, d.bcid,
        d.tdc, corrected_adc, d.overThreshold != 0, corrected_time, 0, event_counter);
      if (result == false ||
        (total_hits >= m_config.nHits && m_config.nHits > 0)) {
        doContinue = false;
        break;
      }
    }
  }

  pcappackets++;
  uint64_t nextFrameCounter =
    m_stats.GetLastFrameCounter(parser -> pd.fecId) + 1;

  if (nextFrameCounter != parser -> hdr.frameCounter) {
    if (parser -> hdr.frameCounter > nextFrameCounter) {
      if (m_stats.GetCounter("ParserGoodFrames", parser -> pd.fecId) >
        0) {
        m_stats.IncrementCounter(
          "ParserFrameMissingErrors", parser -> pd.fecId,
          parser -> hdr.frameCounter - nextFrameCounter);
      }
    } else {
      if (nextFrameCounter - parser -> hdr.frameCounter > 0x0FFFFFFF) {
        m_stats.IncrementCounter("ParserFramecounterOverflows",
          parser -> pd.fecId, 1);
      } else {
        m_stats.IncrementCounter("ParserFrameSeqErrors",
          parser -> pd.fecId, 1);
      }
    }
  } else {
    if (parser -> hdr.frameCounter == 0) {
      m_stats.IncrementCounter("ParserFramecounterOverflows",
        parser -> pd.fecId, 1);
    }
  }
  m_stats.SetLastFrameCounter(parser -> pd.fecId,
    parser -> hdr.frameCounter);

  // nmxstats.ParserFrameSeqErrors = 0;
  // nmxstats.ParserFrameMissingErrors = 0;
  // nmxstats.ParserFramecounterOverflows = 0;

  m_stats.IncrementCounter("ParserReadouts", parser -> pd.fecId,
    stats.ParserReadouts);
  m_stats.IncrementCounter("ParserMarkers", parser -> pd.fecId,
    stats.ParserMarkers);
  m_stats.IncrementCounter("ParserData", parser -> pd.fecId,
    stats.ParserData);
  stats.ParserReadouts = 0;
  stats.ParserMarkers = 0;
  stats.ParserData = 0;

  m_stats.IncrementCounter("ParserTimestampSeqErrors",
    parser -> pd.fecId,
    stats.ParserTimestampSeqErrors);
  m_stats.IncrementCounter("ParserTimestampOverflows",
    parser -> pd.fecId,
    stats.ParserTimestampOverflows);
  stats.ParserTimestampSeqErrors = 0;
  stats.ParserTimestampOverflows = 0;

  m_stats.IncrementCounter("ParserBadFrames", parser -> pd.fecId,
    stats.ParserBadFrames);
  m_stats.IncrementCounter("ParserGoodFrames", parser -> pd.fecId,
    stats.ParserGoodFrames);

  stats.ParserBadFrames = 0;
  stats.ParserGoodFrames = 0;

  m_Clusterer -> SaveDate(pcap.firstPacketSeconds, pcap.firstPacketDate,
    pcap.lastPacketSeconds, pcap.lastPacketDate, 0);
  delete parser;

  m_Clusterer -> FinishAnalysis();

  delete m_Clusterer;
  timeEnd = std::chrono::system_clock::now();

  int elapsed_seconds = std::chrono::duration_cast < std::chrono::milliseconds > (
      timeEnd - timeStart)
    .count();
  int hit_size = 48;

  corryvreckan::Log::setSection("convertFile");
  LOG(INFO) << "****************************************";
  LOG(INFO) << "Stats (analysis):";
  LOG(INFO) << "****************************************";
  LOG(INFO) << "Analysis time: " << std::setprecision(1) << std::fixed <<
    elapsed_seconds << " ms";
  LOG(INFO) << "Hit rate: " << std::scientific <<
    static_cast < double > (1000 * total_hits / elapsed_seconds) <<
    " hit/s";
  LOG(INFO) << "Data rate: " << std::scientific <<
    static_cast < double > (1000 * total_hits * hit_size /
      elapsed_seconds) <<
    " bit/s";
  LOG(INFO) << "Trigger rate: " << std::scientific <<
    static_cast < double > (triggerSignals / elapsed_seconds) <<
    " trigger/s (total triggers: " << triggerSignals << ")";
  LOG(INFO) << "****************************************";
  pcap.printStats();
  return 0;
}