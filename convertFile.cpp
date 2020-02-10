#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <chrono>

#include "Clusterer.h"
#include "Configuration.h"
#include <gdgem/nmx/Readout.h>
#include <gdgem/srs/CalibrationFile.h>
#include <gdgem/srs/ParserVMM3.h>
#include <tools/ReaderPcap.h>

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
    if(m_config.pShowStats) {
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
      Gem::Readout readout;
      Gem::CalibrationFile calfile(m_config.pCalFilename);
      ReaderPcap pcap(m_config.pFileName);
      if (pcap.open() < 0) {
        std::cout << "Error opening file: " << m_config.pFileName
                  << std::endl;
        return -1;
      }
      uint64_t pcappackets = 0;

      int rdsize;
      while ((rdsize = pcap.read((char *)&buffer, sizeof(buffer))) != -1) {
        if (rdsize == 0) {
          continue; // non udp data
        }
        int hits = parser->receive(buffer, rdsize);
        total_hits += hits;
        // field fec id starts at 1
        readout.fec = parser->pd.fecId;
        /*
        Readsout
        uint8_t fec{0};
        uint8_t chip_id{0};
        uint64_t srs_timestamp{0};
        uint16_t channel{0};
        uint16_t bcid{0};
        uint16_t tdc{0};
        uint16_t adc{0};
        bool over_threshold{false};
        float chiptime{0.0};
        */

	    for (int i = 0; i < hits; i++) {
          auto &d = parser->data[i];
          if (d.hasDataMarker && d.fecTimeStamp > 0) {

            readout.srs_timestamp =
                (static_cast<uint64_t>(d.fecTimeStamp) *
                     Gem::SRSTime::internal_SRS_clock_period_ns +
                 static_cast<uint64_t>(d.triggerOffset) *
                     srs_time.trigger_period_ns());

            readout.chip_id = d.vmmid;
            readout.channel = d.chno;
            readout.bcid = d.bcid;
            readout.tdc = d.tdc;
            readout.over_threshold = (d.overThreshold != 0);

            auto calib = calfile.getCalibration(readout.fec, readout.chip_id,
                                                readout.channel);
            readout.chiptime = static_cast<float>(srs_time.chip_time_ns(
                d.bcid, d.tdc, calib.time_offset, calib.time_slope));
            readout.adc = (d.adc - calib.adc_offset) * calib.adc_slope;

            bool result = m_Clusterer->AnalyzeHits(
                static_cast<double>(readout.srs_timestamp), readout.fec,
                readout.chip_id, readout.channel, readout.bcid, readout.tdc,
                readout.adc, readout.over_threshold, readout.chiptime);
            if (result == false || (total_hits >= m_config.nHits &&
                                    m_config.nHits > 0))
              break;
          }
        }
        if(m_config.pShowStats) {
          pcappackets++;
          m_stats.IncrementCounter("parser_frame_seq_errors", readout.fec,
                                  nmxstats.parser_frame_seq_errors);
          m_stats.IncrementCounter("parser_frame_missing_errors", readout.fec,
                                  nmxstats.parser_frame_missing_errors);
          m_stats.IncrementCounter("parser_framecounter_overflows", readout.fec,
                                  nmxstats.parser_framecounter_overflows);
          nmxstats.parser_frame_seq_errors = 0;
          nmxstats.parser_frame_missing_errors = 0;
          nmxstats.parser_framecounter_overflows = 0;

          m_stats.IncrementCounter("parser_readouts",readout.fec,nmxstats.parser_readouts);
          m_stats.IncrementCounter("parser_markers",readout.fec,nmxstats.parser_markers);
          m_stats.IncrementCounter("parser_data",readout.fec,nmxstats.parser_data);
          nmxstats.parser_readouts=0;
          nmxstats.parser_markers=0;
          nmxstats.parser_data=0;

          m_stats.IncrementCounter("parser_timestamp_seq_errors",readout.fec,nmxstats.parser_timestamp_seq_errors);
          m_stats.IncrementCounter("parser_timestamp_overflows",readout.fec,nmxstats.parser_timestamp_overflows);
          nmxstats.parser_timestamp_seq_errors=0;
          nmxstats.parser_timestamp_overflows=0;
        }
      
      }
      if(m_config.pShowStats) {
        m_stats.IncrementCounter("parser_bad_frames", readout.fec,
                               nmxstats.parser_bad_frames);
        m_stats.IncrementCounter("parser_good_frames", readout.fec,
                               nmxstats.parser_good_frames);
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
          RowData.chiptime = bcTime + (m_config.pBCTime_ns - tdcTime -
                                       calib.time_offset) *
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
		    m_stats.IncrementCounter("parser_data",RowData.fec);
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

    std::cout << "Analyzed " << total_hits << " hits. Finished computation in "
              << elapsed_seconds << " ms\n";
  }
  return 0;
}
