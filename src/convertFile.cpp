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
#include <parser/ReaderPcap.h>
#include <parser/ReadoutSRS.h>
#include <parser/VMM3Parser.h>

using namespace hdf5;

int main(int argc, char **argv) {
  uint64_t total_hits = 0;
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
      if (m_config.pDataFormat == "SRS") {
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
        int lastFecID = 0;
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

            int64_t triggerOffset = -1;

            // triggerOffset goes from -1 to 15, a value of -16 indicates a
            // latency violation
            if (d.triggerOffset <= 15) {
              triggerOffset = d.triggerOffset;
            }
            if (d.hasDataMarker && d.fecTimeStamp > 0) {
              double srs_timestamp =
                  (static_cast<uint64_t>(d.fecTimeStamp) * m_config.pBCTime_ns +
                   m_config.pOffsetPeriod * triggerOffset);
              auto calib =
                  calfile.getCalibration(parser->pd.fecId, d.vmmid, d.chno);
              float chiptime =
                  static_cast<double>(d.bcid) * m_config.pBCTime_ns +
                  (m_config.pBCTime_ns -
                   static_cast<double>(d.tdc) *
                       static_cast<double>(m_config.pTAC) / 255 -
                   calib.time_offset) *
                      calib.time_slope;
              if (calib.adc_slope == 0) {
                std::cout << "Error in calibration file: adc_slope correction "
                             "for fec "
                          << parser->pd.fecId << ", chip " << d.vmmid
                          << ", channel " << d.chno
                          << " is 0!\nIs that intentional?" << std::endl;
              }
              int corrected_adc = (d.adc - calib.adc_offset) * calib.adc_slope;

              if (corrected_adc > 2047) {
                std::cout << "After correction, ADC value much larger than "
                             "10bit!  Uncorrected ADC value:"
                          << d.adc << ", corrected value " << corrected_adc
                          << std::endl;
                corrected_adc = 2047;
              } else if (corrected_adc < 0) {
                std::cout << "After correction, ADC value smaller than 0!  "
                             "Uncorrected ADC value:"
                          << d.adc << ", corrected value " << corrected_adc
                          << std::endl;
                corrected_adc = 0;
              }
              uint16_t adc = static_cast<uint16_t>(corrected_adc);
              bool result = m_Clusterer->AnalyzeHits(
                  srs_timestamp, parser->pd.fecId, d.vmmid, d.chno, d.bcid,
                  d.tdc, adc, d.overThreshold != 0, chiptime);
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
        while (doContinue &&
               (rdsize = pcap.read((char *)&buffer, sizeof(buffer))) != -1) {
          if (rdsize == 0) {
            continue; // non udp data
          }
          // Freia 0x48
          // NMX 0x44
          int ret =
              readoutParser.validate((char *)&buffer, rdsize, ReadoutParser::FREIA);
          if (m_config.pShowStats) {
            pcappackets++;
            if (ret != ReadoutParser::OK) {
              badFrames++;
              continue;
            }
            else {
              goodFrames++;
            } 
          }
          int hits = parser->parse(readoutParser.Packet.DataPtr, readoutParser.Packet.DataLength);
          total_hits += hits;
          
          for (int i = 0; i < hits; i++) {
            auto &hit = parser->Result[i];
            if(firstTime == 0) {
              firstTime = hit.TimeHigh * 1.0E+09;
            }
            double complete_timestamp = hit.TimeHigh * 1.0E+09 - firstTime + hit.TimeLow * m_config.pBCTime_ns * 0.5;
            uint16_t adc = hit.OTADC & 0x03ff;
            bool overThreshold = hit.OTADC & 0x8000;
            uint16_t assisterId =
                static_cast<uint8_t>(hit.RingId / 2) * 32 + hit.FENId;
            auto calib =
                calfile.getCalibration(assisterId, hit.VMM, hit.Channel);
            float chiptime_correction =
                (m_config.pBCTime_ns -
                 static_cast<double>(hit.TDC) *
                     static_cast<double>(m_config.pTAC) / 255 -
                 calib.time_offset) *
                calib.time_slope;
            if (calib.adc_slope == 0) {
              std::cout << "Error in calibration file: adc_slope correction "
                           "for assister "
                        << assisterId << ", chip " << hit.VMM << ", channel "
                        << hit.Channel << " is 0!\nIs that intentional?"
                        << std::endl;
            }
            int corrected_adc = (adc - calib.adc_offset) * calib.adc_slope;

            if (corrected_adc > 2047) {
              std::cout << "After correction, ADC value much larger than "
                           "10bit!  Uncorrected ADC value:"
                        << adc << ", corrected value " << corrected_adc
                        << std::endl;
              corrected_adc = 2047;
            } else if (corrected_adc < 0) {
              std::cout << "After correction, ADC value smaller than 0!  "
                           "Uncorrected ADC value:"
                        << adc << ", corrected value " << corrected_adc
                        << std::endl;
              corrected_adc = 0;
            }
  
            bool result = m_Clusterer->AnalyzeHits(
                complete_timestamp, assisterId, hit.VMM, hit.Channel, hit.BC,
                hit.TDC, static_cast<uint16_t>(corrected_adc),
                overThreshold != 0, chiptime_correction);
            if (result == false ||
                (total_hits >= m_config.nHits && m_config.nHits > 0)) {
              doContinue = false;
              break;
            }
          }
        }
        m_stats.IncrementCounter("ErrorBuffer", 384,readoutParser.Stats.ErrorBuffer);
        m_stats.IncrementCounter("ErrorSize", 384,readoutParser.Stats.ErrorSize);
        m_stats.IncrementCounter("ErrorVersion", 384,readoutParser.Stats.ErrorVersion);
        m_stats.IncrementCounter("ErrorCookie", 384,readoutParser.Stats.ErrorCookie);
        m_stats.IncrementCounter("ErrorPad", 384,readoutParser.Stats.ErrorPad);
        m_stats.IncrementCounter("ErrorOutputQueue", 384,readoutParser.Stats.ErrorOutputQueue);
        m_stats.IncrementCounter("ErrorTypeSubType", 384,readoutParser.Stats.ErrorTypeSubType);
        m_stats.IncrementCounter("ErrorSeqNum", 384,readoutParser.Stats.ErrorSeqNum);
        m_stats.IncrementCounter("ErrorTimeHigh", 384,readoutParser.Stats.ErrorTimeHigh);
        m_stats.IncrementCounter("ErrorTimeFrac", 384,readoutParser.Stats.ErrorTimeFrac);
        m_stats.IncrementCounter("HeartBeats", 384,readoutParser.Stats.HeartBeats);
        m_stats.IncrementCounter("GoodFrames", 384,goodFrames);
        m_stats.IncrementCounter("BadFrames", 384,badFrames);
        m_stats.IncrementCounter("TotalFrames", 384,pcappackets);
        
        m_stats.IncrementCounter("ParserErrorSize", 384,parser->Stats.ErrorSize);
        m_stats.IncrementCounter("ParserErrorRing", 384,parser->Stats.ErrorRing);        
        m_stats.IncrementCounter("ParserErrorFEN", 384,parser->Stats.ErrorFEN);
        m_stats.IncrementCounter("ParserErrorDataLength", 384,parser->Stats.ErrorDataLength);        
        m_stats.IncrementCounter("ParserErrorTimeFrac", 384,parser->Stats.ErrorTimeFrac);
        m_stats.IncrementCounter("ParserErrorBC", 384,parser->Stats.ErrorBC);        
        m_stats.IncrementCounter("ParserErrorADC", 384,parser->Stats.ErrorADC);
        m_stats.IncrementCounter("ParserErrorVMM", 384,parser->Stats.ErrorVMM);        
        m_stats.IncrementCounter("ParserErrorChannel", 384,parser->Stats.ErrorChannel);
        m_stats.IncrementCounter("ParserReadouts", 384,parser->Stats.Readouts);        
        m_stats.IncrementCounter("ParserCalibReadouts", 384,parser->Stats.CalibReadouts);
        m_stats.IncrementCounter("ParserDataReadouts", 384,parser->Stats.DataReadouts);        
        m_stats.IncrementCounter("ParserOverThreshold", 384,parser->Stats.OverThreshold);        
     
      }
    } else {
      auto DataFile = file::open(m_config.pFileName);
      auto RootGroup = DataFile.root();
      auto Dataset = RootGroup.get_dataset("srs_hits");
      dataspace::Simple Dataspace(Dataset.dataspace());
      std::vector<Gem::ReadoutSRS> AllElements(Dataspace.size());

      Dataset.read(AllElements);
      /*
      std::sort(AllElements.begin(), AllElements.end(), [](const Readout& lhs,
      const Readout& rhs) { return lhs.srs_timestamp < rhs.srs_timestamp;
      });
      */
      Gem::CalibrationFile calfile(m_config.pCalFilename);
      for (auto RowData : AllElements) {
        /*
        if(RowData.chip_id == 4 || RowData.chip_id == 5) {
                if(RowData.channel%2 == 1) {
                        RowData.channel = RowData.channel - 1;
                }
                else {
                        RowData.channel = RowData.channel + 1;
                }
        }
        */
        auto calib = calfile.getCalibration(RowData.fec, RowData.chip_id,
                                            RowData.channel);

        // User calibration added to analysis
        // The chiptime has to be recalculated in this case from bcid and tdc
        if (calib.time_offset != 0 || calib.time_slope != 1.0) {
          double bcTime = m_config.pBCTime_ns * RowData.bcid;
          double tdcTime = RowData.tdc * m_config.pTAC / 255;
          RowData.chiptime =
              bcTime + (m_config.pBCTime_ns - tdcTime - calib.time_offset) *
                           calib.time_slope;
        }
        if (calib.adc_slope == 0) {
          std::cout
              << "Error in calibration file: adc_slope correction for fec "
              << RowData.fec << ", chip " << RowData.chip_id << ", channel "
              << RowData.channel << " is 0!\nIs that intentional?" << std::endl;
        }
        int corrected_adc = (RowData.adc - calib.adc_offset) * calib.adc_slope;

        if (corrected_adc > 2047) {
          std::cout << "After correction, ADC value much than 10bit!  "
                       "Uncorrected ADC value:"
                    << RowData.adc << ", corrected value " << corrected_adc
                    << std::endl;
          corrected_adc = 2047;
        } else if (corrected_adc < 0) {
          std::cout << "After correction, ADC value smaller than 0!  "
                       "Uncorrected ADC value:"
                    << RowData.adc << ", corrected value " << corrected_adc
                    << std::endl;
          corrected_adc = 0;
        }
        RowData.adc = static_cast<uint16_t>(corrected_adc);
        bool result = m_Clusterer->AnalyzeHits(
            static_cast<double>(RowData.srs_timestamp), RowData.fec,
            RowData.chip_id, RowData.channel, RowData.bcid, RowData.tdc,
            RowData.adc, RowData.over_threshold, RowData.chiptime);
        total_hits++;
        m_stats.IncrementCounter("ParserData", RowData.fec);
        if (result == false ||
            (total_hits >= m_config.nHits && m_config.nHits > 0))
          break;
      }
    }

    m_Clusterer->FinishAnalysis();

    delete m_Clusterer;
    timeEnd = std::chrono::system_clock::now();

    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(
                              timeEnd - timeStart)
                              .count();

    std::cout << "\n****************************************" << std::endl;
    std::cout << "Stats (analysis):" << std::endl;
    std::cout << "****************************************" << std::endl;
    std::cout << "Analysis time: " << std::setprecision(1) << std::fixed
              << elapsed_seconds << " ms" << std::endl;
    std::cout << "Hit rate: " << std::scientific
              << static_cast<double>(1000 * total_hits / elapsed_seconds)
              << " hit/s" << std::endl;
    std::cout << "Data rate: " << std::scientific
              << static_cast<double>(1000 * total_hits * 48 / elapsed_seconds)
              << " bit/s" << std::endl;
    std::cout << "****************************************" << std::endl;
  }
  return 0;
}
