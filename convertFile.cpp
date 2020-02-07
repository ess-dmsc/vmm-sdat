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

  Configuration m_configuration;
  Statistics m_stats;

  if (m_configuration.ParseCommandLine(argc, argv)) {
    if (!m_configuration.CreateMapping()) {
      return -1;
    }
    if (!m_configuration.CalculateTransform()) {
      return -1;
    }

    timeStart = std::chrono::system_clock::now();
	uint64_t last_time = 0;

    Clusterer *m_Clusterer = new Clusterer(m_configuration, m_stats);
    m_stats.CreateFECStats(m_configuration);
    m_stats.CreateClusterStats(m_configuration);
    if (m_configuration.isPcap) {
      m_stats.CreatePCAPStats(m_configuration);
      char buffer[10000];
      Gem::SRSTime srs_time;
	  srs_time.bc_clock_MHz(m_configuration.pBC);
  	  srs_time.tac_slope_ns(m_configuration.pTAC);
 
      Gem::NMXStats nmxstats;
      Gem::ParserVMM3 *parser = new Gem::ParserVMM3(2000, nmxstats, srs_time);
      Gem::Readout readout;
      Gem::CalibrationFile calfile(m_configuration.pCalFilename);
      ReaderPcap pcap(m_configuration.pFileName);
      if (pcap.open() < 0) {
        std::cout << "Error opening file: " << m_configuration.pFileName
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
            // std::cout << (int)readout.channel  << " " << (int) readout.adc <<
            // std::endl;
            readout.bcid = d.bcid;
            readout.tdc = d.tdc;
            readout.over_threshold = (d.overThreshold != 0);

            auto calib = calfile.getCalibration(readout.fec, readout.chip_id,
                                                readout.channel);
            readout.chiptime = static_cast<float>(srs_time.chip_time_ns(
                d.bcid, d.tdc, calib.time_offset, calib.time_slope));
            readout.adc = (d.adc - calib.adc_offset) * calib.adc_slope;

            // std::cout << "pcappackets: " << pcappackets << ": " << (uint32_t)
            // readout.fec << ", " << (uint32_t) readout.chip_id << ", " <<
            // readout.srs_timestamp << ", "
            //<< (uint32_t)readout.channel << ", " << (uint32_t)readout.bcid <<
            //", " << (uint32_t)readout.tdc << ", " << (uint32_t)readout.adc <<
            //", "
            //<< (uint32_t)readout.over_threshold << ", " << readout.chiptime <<
            // std::endl;

            bool result = m_Clusterer->AnalyzeHits(
                static_cast<double>(readout.srs_timestamp), readout.fec,
                readout.chip_id, readout.channel, readout.bcid, readout.tdc,
                readout.adc, readout.over_threshold, readout.chiptime);
            if (result == false || (total_hits >= m_configuration.nHits &&
                                    m_configuration.nHits > 0))
              break;
          }
        }
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
      m_stats.IncrementCounter("parser_bad_frames", readout.fec,
                               nmxstats.parser_bad_frames);
      m_stats.IncrementCounter("parser_good_frames", readout.fec,
                               nmxstats.parser_good_frames);

      delete parser;
    } else {
      auto DataFile = file::open(m_configuration.pFileName);
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
      Gem::CalibrationFile calfile(m_configuration.pCalFilename);
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
          double bcTime = m_configuration.pBCTime_ns * RowData.bcid;
          double tdcTime = RowData.tdc * m_configuration.pTAC / 255;
          RowData.chiptime = bcTime + (m_configuration.pBCTime_ns - tdcTime -
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
        // std::cout << "hit: " << lines << ": " << (uint32_t) RowData.fec << ",
        // " << (uint32_t) RowData.chip_id << ", " << RowData.srs_timestamp <<
        // ", " << RowData.channel << ", " << RowData.bcid << ", " <<
        // RowData.tdc <<
        // ", " << RowData.adc << ", " << RowData.over_threshold << ", " <<
        // RowData.chiptime << std::endl;
        total_hits++;
		m_stats.IncrementCounter("parser_data",RowData.fec);
        if (result == false ||
            (total_hits >= m_configuration.nHits && m_configuration.nHits > 0))
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
