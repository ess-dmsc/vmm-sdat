#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <chrono>

#include "Clusterer.h"
#include "Configuration.h"
#include <gdgem/srs/CalibrationFile.h>
#include <gdgem/srs/ParserVMM3.h>
#include <udpgenpcap/ReaderPcap.h>
#include <gdgem/nmx/Readout.h>

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
    m_stats.CreateFECStats(m_config);
    if (m_config.pShowStats) {
      m_stats.CreateClusterStats(m_config);
      if (m_config.pIsPcap) {
        m_stats.CreatePCAPStats(m_config);
      }
    }
    if (m_config.pIsPcap) {
      char buffer[10000];
      Gem::SRSTime srs_time;
      srs_time.bc_clock_MHz(m_config.pBC);
      srs_time.tac_slope_ns(m_config.pTAC);

      Gem::NMXStats nmxstats;
      Gem::ParserVMM3 *parser = new Gem::ParserVMM3(2000, nmxstats, srs_time);
      Gem::CalibrationFile calfile(m_config.pCalFilename);
      ReaderPcap pcap(m_config.pFileName);
      int ret = pcap.open();
      if (ret < 0) {
        std::cout << "Error opening file: " << m_config.pFileName 
        << ": return value " << ret<<std::endl;
        return -1;
      }
      uint64_t pcappackets = 0;

      int rdsize;
      bool doContinue = true;
      while (doContinue && (rdsize = pcap.read((char *)&buffer, sizeof(buffer))) != -1) {
        if (rdsize == 0) {
          continue; // non udp data
        }
        int hits = parser->receive(buffer, rdsize);
        total_hits += hits;

        for (int i = 0; i < hits; i++) {
          auto &d = parser->data[i];
          if (d.hasDataMarker && d.fecTimeStamp > 0) {

            double srs_timestamp =
                (static_cast<uint64_t>(d.fecTimeStamp) *
                     Gem::SRSTime::internal_SRS_clock_period_ns +
                 static_cast<uint64_t>(d.triggerOffset) *
                     srs_time.trigger_period_ns());
          
            auto calib = calfile.getCalibration(
                parser->pd.fecId, d.vmmid, d.chno);
            float chiptime = static_cast<float>(srs_time.chip_time_ns(
                d.bcid, d.tdc, calib.time_offset, calib.time_slope));
            uint16_t adc = (d.adc - calib.adc_offset) * calib.adc_slope;

            bool result = m_Clusterer->AnalyzeHits(
                srs_timestamp, parser->pd.fecId, d.vmmid, d.chno, d.bcid, d.tdc,
                adc, d.overThreshold != 0, chiptime);
            if (result == false ||
                (total_hits >= m_config.nHits && m_config.nHits > 0)) {
                doContinue = false;
              break;
            }
          }
        }
        if (m_config.pShowStats) {
          pcappackets++;
          m_stats.IncrementCounter("ParserFrameSeqErrors", parser->pd.fecId,
                                   nmxstats.ParserFrameSeqErrors);
          m_stats.IncrementCounter("ParserFrameMissingErrors", parser->pd.fecId,
                                   nmxstats.ParserFrameMissingErrors);
          m_stats.IncrementCounter("ParserFramecounterOverflows", parser->pd.fecId,
                                   nmxstats.ParserFramecounterOverflows);
          nmxstats.ParserFrameSeqErrors = 0;
          nmxstats.ParserFrameMissingErrors = 0;
          nmxstats.ParserFramecounterOverflows = 0;

          m_stats.IncrementCounter("ParserReadouts", parser->pd.fecId,
                                   nmxstats.ParserReadouts);
          m_stats.IncrementCounter("ParserMarkers", parser->pd.fecId,
                                   nmxstats.ParserMarkers);
          m_stats.IncrementCounter("ParserData", parser->pd.fecId,
                                   nmxstats.ParserData);
          nmxstats.ParserReadouts = 0;
          nmxstats.ParserMarkers = 0;
          nmxstats.ParserData = 0;

          m_stats.IncrementCounter("ParserTimestampSeqErrors", parser->pd.fecId,
                                   nmxstats.ParserTimestampSeqErrors);
          m_stats.IncrementCounter("ParserTimestampOverflows", parser->pd.fecId,
                                   nmxstats.ParserTimestampOverflows);
          nmxstats.ParserTimestampSeqErrors = 0;
          nmxstats.ParserTimestampOverflows = 0;
        }
      }
      if (m_config.pShowStats) {
        m_stats.IncrementCounter("ParserBadFrames", parser->pd.fecId,
                                 nmxstats.ParserBadFrames);
        m_stats.IncrementCounter("ParserGoodFrames", parser->pd.fecId,
                                 nmxstats.ParserGoodFrames);
      }
      delete parser;
    } else {
      auto DataFile = file::open(m_config.pFileName);
      auto RootGroup = DataFile.root();
      auto Dataset = RootGroup.get_dataset("srs_hits");
      dataspace::Simple Dataspace(Dataset.dataspace());
      std::vector<Gem::Readout> AllElements(Dataspace.size());

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
        RowData.adc = static_cast<uint16_t>((RowData.adc - calib.adc_offset) *
                                            calib.adc_slope);

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
    std::cout << "Analysis time: " << elapsed_seconds << " ms" << std::endl;
    std::cout << "Hit rate: " << std::scientific << static_cast<double>(1000*total_hits/elapsed_seconds) << " hit/s" << std::endl;
    std::cout << "Data rate: " << std::scientific  << static_cast<double>(1000*total_hits*48/elapsed_seconds) << " bit/s" << std::endl;
    std::cout << "****************************************" << std::endl;
          
  }
  return 0;
}
