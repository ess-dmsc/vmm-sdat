#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Clusterer.h"
#include "Configuration.h"
#include <parser/CalibrationFile.h>
#include <parser/ParserSRS.h>
#include <parser/ParserVTC.h>
#include <parser/ReaderPcap.h>
#include <parser/Trace.h>
#include <parser/VMM3Parser.h>

int main(int argc, char **argv) {
  uint64_t total_hits = 0;
  uint64_t triggerSignals = 0;
  std::chrono::time_point<std::chrono::system_clock> timeEnd, timeStart;

  Configuration m_config;
  Statistics m_stats;

  if (m_config.ParseCommandLine(argc, argv)) {
    if (!m_config.CreateMapping()) {
      return -1;
    }
    if (!m_config.CalculateTransform()) {
      return -1;
    }

    timeStart = std::chrono::system_clock::now();
    uint64_t last_time = 0;

    Clusterer *m_Clusterer = new Clusterer(m_config, m_stats);
    if (m_config.pDataFormat == "SRS") {
      m_stats.CreateFECStats(m_config);
    }
    if (m_config.pShowStats) {
      m_stats.CreateClusterStats(m_config);
      if (m_config.pIsPcap) {
        m_stats.CreatePCAPStats(m_config);
      }
    }

    if (m_config.pIsPcap) {
      if (m_config.pDataFormat == "VTC") {
        double firstTime = 0;
        char buffer[10000];
        Gem::SRSTime srs_time;
        srs_time.bc_clock_MHz(m_config.pBC);
        srs_time.tac_slope_ns(m_config.pTAC);

        Gem::NMXStats nmxstats;
        Gem::ParserVTC *parser = new Gem::ParserVTC(2000, nmxstats, srs_time);
        Gem::CalibrationFile calfile(m_config.pCalFilename);
        ReaderPcap pcap(m_config.pFileName);
        int ret = pcap.open();
        if (ret < 0) {
          std::cout << "Error opening file: " << m_config.pFileName
                    << ": return value " << ret << std::endl;
          return -1;
        }
        uint64_t pcappackets = 0;
        int rdsize;
        bool doContinue = true;
        while (doContinue &&
               (rdsize = pcap.read((char *)&buffer, sizeof(buffer))) != -1) {
          if (rdsize == 0) {
            continue; // non udp data
          }
          int hits = parser->receive(buffer, rdsize);

          total_hits += hits;

          for (int i = 0; i < hits; i++) {
            auto &d = parser->data[i];
            double vtc_timestamp =
                static_cast<double>(d.clockCounter) * m_config.pOffsetPeriod;
            if (firstTime == 0) {
              firstTime = vtc_timestamp;
            }
            double t0_correction = 0;
            std::pair<uint8_t, uint8_t> fec_vmm = std::make_pair(0, d.vmmid);
            auto searchMap = m_config.pFecVMM_time0.find(fec_vmm);
            if (searchMap != m_config.pFecVMM_time0.end()) {
              std::string t0 = m_config.pFecVMM_time0[fec_vmm];
              if (t0 == "run") {
                t0_correction = firstTime;
              } else {
                t0_correction = std::stod(t0);
              }
            }
            vtc_timestamp = vtc_timestamp - t0_correction;
            auto calib = calfile.getCalibration(0, d.vmmid, d.chno);
            double chiptime_corrected =
                static_cast<double>(d.bcid) * m_config.pBCTime_ns +
                (1.5 * m_config.pBCTime_ns -
                 static_cast<double>(d.tdc) *
                     static_cast<double>(m_config.pTAC) / 255.0 -
                 calib.time_offset) *
                    calib.time_slope;

            uint16_t corrected_adc = static_cast<uint16_t>(
                (static_cast<double>(d.adc) - calib.adc_offset) *
                calib.adc_slope);

            if (corrected_adc > 1023) {
              DTRACE(DEB,
                     "After correction, ADC value larger than 1023 "
                     "(10bit)!\nUncorrected ADC value %d, uncorrected ADC "
                     "value %d\n",
                     d.adc, corrected_adc);

              corrected_adc = 1023;
            } else if (corrected_adc < 0) {
              DTRACE(DEB,
                     "After correction, ADC value smaller than 0!"
                     "\nUncorrected ADC value %d, uncorrected ADC "
                     "value %d\n",
                     d.adc, corrected_adc);
              corrected_adc = 0;
            }
            double time_without_calib =
                (1.5 * m_config.pBCTime_ns -
                 static_cast<double>(d.tdc) *
                     static_cast<double>(m_config.pTAC) / 255.0);

            double time_with_calib =
                (1.5 * m_config.pBCTime_ns -
                 static_cast<double>(d.tdc) *
                     static_cast<double>(m_config.pTAC) / 255.0 -
                 calib.time_offset) *
                calib.time_slope;

            m_Clusterer->FillCalibHistos(0, d.vmmid, d.chno, d.adc,
                                         corrected_adc, time_without_calib,
                                         time_with_calib);

            double timewalk_correction =
                calib.timewalk_d +
                (calib.timewalk_a - calib.timewalk_d) /
                    (1 +
                     pow(corrected_adc / calib.timewalk_c, calib.timewalk_b));

            double corrected_time = chiptime_corrected - timewalk_correction;

            bool result = m_Clusterer->AnalyzeHits(
                vtc_timestamp, 0, d.vmmid, d.chno, d.bcid, d.tdc, corrected_adc,
                d.overThreshold != 0, corrected_time);
            if (result == false ||
                (total_hits >= m_config.nHits && m_config.nHits > 0)) {
              doContinue = false;
              break;
            }
          }

          if (m_config.pShowStats) {
            pcappackets++;
            uint64_t nextFrameCounter = m_stats.GetLastFrameCounter(0) + 1;

            if (nextFrameCounter != parser->hdr.frameCounter) {
              if (parser->hdr.frameCounter > nextFrameCounter) {
                if (m_stats.GetCounter("ParserGoodFrames", 0) > 0) {
                  m_stats.IncrementCounter("ParserFrameMissingErrors", 0,
                                           parser->hdr.frameCounter -
                                               nextFrameCounter);
                }
              } else {
                if (nextFrameCounter - parser->hdr.frameCounter > 0x0FFFFFFF) {
                  m_stats.IncrementCounter("ParserFramecounterOverflows", 0, 1);
                } else {
                  m_stats.IncrementCounter("ParserFrameSeqErrors", 0, 1);
                }
              }
            } else {
              if (parser->hdr.frameCounter == 0) {
                m_stats.IncrementCounter("ParserFramecounterOverflows", 0, 1);
              }
            }
            m_stats.SetLastFrameCounter(0, parser->hdr.frameCounter);

            // nmxstats.ParserFrameSeqErrors = 0;
            // nmxstats.ParserFrameMissingErrors = 0;
            // nmxstats.ParserFramecounterOverflows = 0;

            m_stats.IncrementCounter("ParserReadouts", 0,
                                     nmxstats.ParserReadouts);
            m_stats.IncrementCounter("ParserMarkers", 0,
                                     nmxstats.ParserMarkers);
            m_stats.IncrementCounter("ParserData", 0, nmxstats.ParserData);
            nmxstats.ParserReadouts = 0;
            nmxstats.ParserMarkers = 0;
            nmxstats.ParserData = 0;

            m_stats.IncrementCounter("ParserTimestampSeqErrors", 0,
                                     nmxstats.ParserTimestampSeqErrors);
            m_stats.IncrementCounter("ParserTimestampOverflows", 0,
                                     nmxstats.ParserTimestampOverflows);
            nmxstats.ParserTimestampSeqErrors = 0;
            nmxstats.ParserTimestampOverflows = 0;

            m_stats.IncrementCounter("ParserBadFrames", 0,
                                     nmxstats.ParserBadFrames);
            m_stats.IncrementCounter("ParserGoodFrames", 0,
                                     nmxstats.ParserGoodFrames);

            nmxstats.ParserBadFrames = 0;
            nmxstats.ParserGoodFrames = 0;
          }
        }
        m_Clusterer->SaveDate(pcap.firstPacketSeconds, pcap.firstPacketDate,
                              pcap.lastPacketSeconds, pcap.lastPacketDate, 0);
        delete parser;
      } else if (m_config.pDataFormat == "SRS") {
        double firstTime = 0;
        char buffer[10000];
        Gem::SRSTime srs_time;
        srs_time.bc_clock_MHz(m_config.pBC);
        srs_time.tac_slope_ns(m_config.pTAC);

        Gem::NMXStats nmxstats;
        Gem::ParserSRS *parser = new Gem::ParserSRS(2000, nmxstats, srs_time);
        Gem::CalibrationFile calfile(m_config.pCalFilename);
        ReaderPcap pcap(m_config.pFileName);
        int ret = pcap.open();
        if (ret < 0) {
          std::cout << "Error opening file: " << m_config.pFileName
                    << ": return value " << ret << std::endl;
          return -1;
        }

        uint64_t pcappackets = 0;
        int rdsize;
        bool doContinue = true;
        while (doContinue &&
               (rdsize = pcap.read((char *)&buffer, sizeof(buffer))) != -1) {
          if (rdsize == 0) {
            continue; // non udp data
          }
          int hits = parser->receive(buffer, rdsize);
          total_hits += hits;

          for (int i = 0; i < hits; i++) {
            auto &d = parser->data[i];

            double triggerOffset = -1.0;

            // triggerOffset goes from -1 to 15
            // but presented as uint8_t
            // latency violation
            if (d.triggerOffset <= 15.0) {
              triggerOffset = static_cast<double>(d.triggerOffset);
            } else if (d.triggerOffset == 31) {
              triggerOffset = -1.0;
            } else if (d.triggerOffset == 16) {
              triggerOffset = -99.0;
            }
            if (d.hasDataMarker && d.fecTimeStamp > 0 && triggerOffset != -99) {

              double srs_timestamp =
                  (static_cast<double>(d.fecTimeStamp) * m_config.pBCTime_ns +
                   m_config.pOffsetPeriod * triggerOffset);
              if (firstTime == 0) {
                firstTime = srs_timestamp;
              }
              double t0_correction = 0;
              std::pair<uint8_t, uint8_t> fec_vmm =
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
                  calfile.getCalibration(parser->pd.fecId, d.vmmid, d.chno);

              double chiptime_corrected =
                  static_cast<double>(d.bcid) * m_config.pBCTime_ns +
                  (1.5 * m_config.pBCTime_ns -
                   static_cast<double>(d.tdc) *
                       static_cast<double>(m_config.pTAC) / 255.0 -
                   calib.time_offset) *
                      calib.time_slope;

              uint16_t corrected_adc = static_cast<uint16_t>(
                  (static_cast<double>(d.adc) - calib.adc_offset) *
                  calib.adc_slope);

              if (corrected_adc > 1023) {
                DTRACE(DEB,
                       "After correction, ADC value larger than 1023 "
                       "(10bit)!\nUncorrected ADC value %d, uncorrected ADC "
                       "value %d\n",
                       d.adc, corrected_adc);

                corrected_adc = 1023;
              } else if (corrected_adc < 0) {
                DTRACE(DEB,
                       "After correction, ADC value smaller than 0!"
                       "\nUncorrected ADC value %d, uncorrected ADC "
                       "value %d\n",
                       d.adc, corrected_adc);
                corrected_adc = 0;
              }

              double time_without_calib =
                  (1.5 * m_config.pBCTime_ns -
                   static_cast<double>(d.tdc) *
                       static_cast<double>(m_config.pTAC) / 255.0);

              double time_with_calib =
                  (1.5 * m_config.pBCTime_ns -
                   static_cast<double>(d.tdc) *
                       static_cast<double>(m_config.pTAC) / 255.0 -
                   calib.time_offset) *
                  calib.time_slope;

              m_Clusterer->FillCalibHistos(parser->pd.fecId, d.vmmid, d.chno,
                                           d.adc, corrected_adc,
                                           time_without_calib, time_with_calib);

              double timewalk_correction =
                  calib.timewalk_d +
                  (calib.timewalk_a - calib.timewalk_d) /
                      (1 +
                       pow(corrected_adc / calib.timewalk_c, calib.timewalk_b));

              double corrected_time = chiptime_corrected - timewalk_correction;

              bool result = m_Clusterer->AnalyzeHits(
                  srs_timestamp, parser->pd.fecId, d.vmmid, d.chno, d.bcid,
                  d.tdc, corrected_adc, d.overThreshold != 0, corrected_time);
              if (result == false ||
                  (total_hits >= m_config.nHits && m_config.nHits > 0)) {
                doContinue = false;
                break;
              }
            }
          }

          if (m_config.pShowStats) {
            pcappackets++;
            uint64_t nextFrameCounter =
                m_stats.GetLastFrameCounter(parser->pd.fecId) + 1;

            if (nextFrameCounter != parser->hdr.frameCounter) {
              if (parser->hdr.frameCounter > nextFrameCounter) {
                if (m_stats.GetCounter("ParserGoodFrames", parser->pd.fecId) >
                    0) {
                  m_stats.IncrementCounter(
                      "ParserFrameMissingErrors", parser->pd.fecId,
                      parser->hdr.frameCounter - nextFrameCounter);
                }
              } else {
                if (nextFrameCounter - parser->hdr.frameCounter > 0x0FFFFFFF) {
                  m_stats.IncrementCounter("ParserFramecounterOverflows",
                                           parser->pd.fecId, 1);
                } else {
                  m_stats.IncrementCounter("ParserFrameSeqErrors",
                                           parser->pd.fecId, 1);
                }
              }
            } else {
              if (parser->hdr.frameCounter == 0) {
                m_stats.IncrementCounter("ParserFramecounterOverflows",
                                         parser->pd.fecId, 1);
              }
            }
            m_stats.SetLastFrameCounter(parser->pd.fecId,
                                        parser->hdr.frameCounter);

            // nmxstats.ParserFrameSeqErrors = 0;
            // nmxstats.ParserFrameMissingErrors = 0;
            // nmxstats.ParserFramecounterOverflows = 0;

            m_stats.IncrementCounter("ParserReadouts", parser->pd.fecId,
                                     nmxstats.ParserReadouts);
            m_stats.IncrementCounter("ParserMarkers", parser->pd.fecId,
                                     nmxstats.ParserMarkers);
            m_stats.IncrementCounter("ParserData", parser->pd.fecId,
                                     nmxstats.ParserData);
            nmxstats.ParserReadouts = 0;
            nmxstats.ParserMarkers = 0;
            nmxstats.ParserData = 0;

            m_stats.IncrementCounter("ParserTimestampSeqErrors",
                                     parser->pd.fecId,
                                     nmxstats.ParserTimestampSeqErrors);
            m_stats.IncrementCounter("ParserTimestampOverflows",
                                     parser->pd.fecId,
                                     nmxstats.ParserTimestampOverflows);
            nmxstats.ParserTimestampSeqErrors = 0;
            nmxstats.ParserTimestampOverflows = 0;

            m_stats.IncrementCounter("ParserBadFrames", parser->pd.fecId,
                                     nmxstats.ParserBadFrames);
            m_stats.IncrementCounter("ParserGoodFrames", parser->pd.fecId,
                                     nmxstats.ParserGoodFrames);

            nmxstats.ParserBadFrames = 0;
            nmxstats.ParserGoodFrames = 0;
          }
        }
        m_Clusterer->SaveDate(pcap.firstPacketSeconds, pcap.firstPacketDate,
                              pcap.lastPacketSeconds, pcap.lastPacketDate, 0);
        delete parser;
      } else if (m_config.pDataFormat == "ESS") {
        double firstTime = 0;
        char buffer[10000];
        Gem::SRSTime srs_time;
        srs_time.bc_clock_MHz(m_config.pBC);
        srs_time.tac_slope_ns(m_config.pTAC);

        VMM3Parser *parser = new VMM3Parser();
        ReadoutParser readoutParser;
        Gem::CalibrationFile calfile(m_config.pCalFilename);
        ReaderPcap pcap(m_config.pFileName);
        int ret = pcap.open();
        if (ret < 0) {
          std::cout << "Error opening file: " << m_config.pFileName
                    << ": return value " << ret << std::endl;
          return -1;
        }
        uint64_t pcappackets = 0;
        uint64_t goodFrames = 0;
        uint64_t badFrames = 0;

        int rdsize;
        bool doContinue = true;
        long seqNumError = 0;
        double pulseTime = 0;

        double t0_correction = 0;
        std::pair<uint8_t, uint8_t> fec_vmm = std::make_pair(0, 0);
        auto searchMap = m_config.pFecVMM_time0.find(fec_vmm);
        if (searchMap != m_config.pFecVMM_time0.end()) {
          std::string t0 = m_config.pFecVMM_time0[fec_vmm];
          t0_correction = std::stod(t0);
        }

        while (doContinue &&
               (rdsize = pcap.read((char *)&buffer, sizeof(buffer))) != -1) {
          if (rdsize == 0) {
            continue; // non udp data
          }
          // Freia 0x48
          // NMX 0x44
          // Loki4Amp 0x30
          int ret = readoutParser.validate((char *)&buffer, rdsize,
                                           ReadoutParser::Loki4Amp);

          if (seqNumError != readoutParser.Stats.ErrorSeqNum) {
            printf("Sequence number error at packet %" PRIu64 "\n",
                   pcappackets);
          }
          seqNumError = readoutParser.Stats.ErrorSeqNum;

          if (m_config.pShowStats) {
            pcappackets++;
            if (ret != ReadoutParser::OK) {
              badFrames++;
              continue;
            } else {
              goodFrames++;
            }
          }

          double temp_pulseTime = 0;
          if (t0_correction > 0) {
            temp_pulseTime =
                (readoutParser.Packet.HeaderPtr->PulseHigh - t0_correction) *
                    1.0E+09 +
                readoutParser.Packet.HeaderPtr->PulseLow * m_config.pBCTime_ns *
                    0.5;

            // Filter out pulse times that come directly after a valid pulse
            // time due to jitter
            if (temp_pulseTime - pulseTime > 1.0E+06) {
              pulseTime = temp_pulseTime;
              triggerSignals++;
            }
          }

          int hits = parser->parse(readoutParser.Packet.DataPtr,
                                   readoutParser.Packet.DataLength);
          total_hits += hits;
          for (int i = 0; i < hits; i++) {
            auto &hit = parser->Result[i];
            uint16_t assisterId =
                static_cast<uint8_t>(hit.RingId / 2) * 32 + hit.FENId;
            if (firstTime == 0) {
              firstTime = hit.TimeHigh * 1.0E+09;
            }

            double complete_timestamp = 0;
            if (t0_correction == 0) {
              //  ESS time format use 64bit timestamp in nanoseconds
              //  To be able to use double as type for timestamp calculation,
              //  the timestamp has to be truncated to 52 bit
              //  The easiest to have a relative timestamp with respect to the
              //  start of the run
              complete_timestamp = (0 + hit.TimeHigh) * 1.0E+09 - firstTime +
                                   hit.TimeLow * m_config.pBCTime_ns * 0.5;
            } else {
              // If several files will be joined later, it is recommended to
              // just remove the most significant bits of the 64 bit timestamp
              complete_timestamp = (hit.TimeHigh - t0_correction) * 1.0E+09 +
                                   hit.TimeLow * m_config.pBCTime_ns * 0.5;
            }

            uint16_t adc = hit.OTADC & 0x03ff;
            bool overThreshold = hit.OTADC & 0x8000;

            m_stats.IncrementCounter("ParserDataReadouts", assisterId, 1);
            auto calib =
                calfile.getCalibration(assisterId, hit.VMM, hit.Channel);
            double chiptime_corrected =
                (1.5 * m_config.pBCTime_ns -
                 static_cast<double>(hit.TDC) *
                     static_cast<double>(m_config.pTAC) / 255 -
                 calib.time_offset) *
                calib.time_slope;

            if (calib.adc_slope == 0) {
              std::cout << "Error in calibration file: adc_slope correction "
                           "for assister "
                        << (int)assisterId << ", chip " << (int)hit.VMM
                        << ", channel " << (int)hit.Channel
                        << " is 0!\nIs that intentional?" << std::endl;
            }
            int corrected_adc = (adc - calib.adc_offset) * calib.adc_slope;

            if (corrected_adc > 1023) {
              DTRACE(DEB,
                     "After correction, ADC value larger than 1023 "
                     "(10bit)!\nUncorrected ADC value %d, uncorrected ADC "
                     "value %d\n",
                     adc, corrected_adc);
              corrected_adc = 1023;
            } else if (corrected_adc < 0) {
              DTRACE(DEB,
                     "After correction, ADC value smaller than 0!"
                     "\nUncorrected ADC value %d, uncorrected ADC "
                     "value %d\n",
                     adc, corrected_adc);
              corrected_adc = 0;
            }
            double timewalk_correction =
                calib.timewalk_d +
                (calib.timewalk_a - calib.timewalk_d) /
                    (1 +
                     pow(corrected_adc / calib.timewalk_c, calib.timewalk_b));

            double corrected_time = chiptime_corrected - timewalk_correction;

            double time_without_calib =
                (1.5 * m_config.pBCTime_ns -
                 static_cast<double>(hit.TDC) *
                     static_cast<double>(m_config.pTAC) / 255.0);

            double time_with_calib =
                (1.5 * m_config.pBCTime_ns -
                 static_cast<double>(hit.TDC) *
                     static_cast<double>(m_config.pTAC) / 255.0 -
                 calib.time_offset) *
                calib.time_slope;

            m_Clusterer->FillCalibHistos(assisterId, hit.VMM, hit.Channel, adc,
                                         corrected_adc, time_without_calib,
                                         time_with_calib);

            bool result = m_Clusterer->AnalyzeHits(
                complete_timestamp, assisterId, hit.VMM, hit.Channel, hit.BC,
                hit.TDC, static_cast<uint16_t>(corrected_adc),
                overThreshold != 0, corrected_time, hit.GEO, pulseTime);
            if (result == false ||
                (total_hits >= m_config.nHits && m_config.nHits > 0)) {
              doContinue = false;
              break;
            }
          }
        }
        m_Clusterer->SaveDate(pcap.firstPacketSeconds, pcap.firstPacketDate,
                              pcap.lastPacketSeconds, pcap.lastPacketDate,
                              triggerSignals);
        m_stats.IncrementCounter("ErrorBuffer", 384,
                                 readoutParser.Stats.ErrorBuffer);
        m_stats.IncrementCounter("ErrorSize", 384,
                                 readoutParser.Stats.ErrorSize);
        m_stats.IncrementCounter("ErrorVersion", 384,
                                 readoutParser.Stats.ErrorVersion);
        m_stats.IncrementCounter("ErrorCookie", 384,
                                 readoutParser.Stats.ErrorCookie);
        m_stats.IncrementCounter("ErrorPad", 384, readoutParser.Stats.ErrorPad);
        m_stats.IncrementCounter("ErrorOutputQueue", 384,
                                 readoutParser.Stats.ErrorOutputQueue);
        m_stats.IncrementCounter("ErrorTypeSubType", 384,
                                 readoutParser.Stats.ErrorTypeSubType);
        m_stats.IncrementCounter("ErrorSeqNum", 384,
                                 readoutParser.Stats.ErrorSeqNum);
        m_stats.IncrementCounter("ErrorTimeHigh", 384,
                                 readoutParser.Stats.ErrorTimeHigh);
        m_stats.IncrementCounter("ErrorTimeFrac", 384,
                                 readoutParser.Stats.ErrorTimeFrac);
        m_stats.IncrementCounter("HeartBeats", 384,
                                 readoutParser.Stats.HeartBeats);
        m_stats.IncrementCounter("GoodFrames", 384, goodFrames);
        m_stats.IncrementCounter("BadFrames", 384, badFrames);
        m_stats.IncrementCounter("TotalFrames", 384, pcappackets);

        m_stats.IncrementCounter("ParserErrorSize", 384,
                                 parser->Stats.ErrorSize);
        m_stats.IncrementCounter("ParserErrorRing", 384,
                                 parser->Stats.ErrorRing);
        m_stats.IncrementCounter("ParserErrorFEN", 384, parser->Stats.ErrorFEN);
        m_stats.IncrementCounter("ParserErrorDataLength", 384,
                                 parser->Stats.ErrorDataLength);
        m_stats.IncrementCounter("ParserErrorTimeFrac", 384,
                                 parser->Stats.ErrorTimeFrac);
        m_stats.IncrementCounter("ParserErrorBC", 384, parser->Stats.ErrorBC);
        m_stats.IncrementCounter("ParserErrorADC", 384, parser->Stats.ErrorADC);
        m_stats.IncrementCounter("ParserErrorVMM", 384, parser->Stats.ErrorVMM);
        m_stats.IncrementCounter("ParserErrorChannel", 384,
                                 parser->Stats.ErrorChannel);
        m_stats.IncrementCounter("ParserReadouts", 384, parser->Stats.Readouts);
        m_stats.IncrementCounter("ParserCalibReadouts", 384,
                                 parser->Stats.CalibReadouts);
        m_stats.IncrementCounter("ParserDataReadouts", 384,
                                 parser->Stats.DataReadouts);
        m_stats.IncrementCounter("ParserOverThreshold", 384,
                                 parser->Stats.OverThreshold);
      }
    }
    m_Clusterer->FinishAnalysis();

    delete m_Clusterer;
    timeEnd = std::chrono::system_clock::now();

    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(
                              timeEnd - timeStart)
                              .count();
    int hit_size = 160;
    if (m_config.pDataFormat == "VTC") {
      hit_size = 64;
    } else if (m_config.pDataFormat == "SRS") {
      hit_size = 48;
    }
    std::cout << "\n****************************************" << std::endl;
    std::cout << "Stats (analysis):" << std::endl;
    std::cout << "****************************************" << std::endl;
    std::cout << "Analysis time: " << std::setprecision(1) << std::fixed
              << elapsed_seconds << " ms" << std::endl;
    std::cout << "Hit rate: " << std::scientific
              << static_cast<double>(1000 * total_hits / elapsed_seconds)
              << " hit/s" << std::endl;
    std::cout << "Data rate: " << std::scientific
              << static_cast<double>(1000 * total_hits * hit_size /
                                     elapsed_seconds)
              << " bit/s" << std::endl;
    std::cout << "Trigger rate: " << std::scientific
              << static_cast<double>(triggerSignals / elapsed_seconds)
              << " trigger/s (total triggers: " << triggerSignals << ")"
              << std::endl;
    std::cout << "****************************************" << std::endl;
  }
  return 0;
}
