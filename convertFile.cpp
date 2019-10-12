#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include <chrono>

#include "Configuration.h"
#include "Clusterer.h"
#include <gdgem/nmx/Readout.h>
#include <gdgem/srs/CalibrationFile.h>

using namespace hdf5;

int main(int argc, char**argv) {
    int lines = 0;
    std::chrono::time_point<std::chrono::system_clock> timeEnd, timeStart;

    Configuration m_configuration;
    Statistics m_stats;

    if (m_configuration.ParseCommandLine(argc, argv)) {
        if(!m_configuration.CreateMapping())
        {
            return -1;
        }
        if(!m_configuration.CalculateTransform())
        {
            return -1;
        }
        m_stats.CreateStats(m_configuration);
        timeStart = std::chrono::system_clock::now();
       
        Clusterer* m_Clusterer = new Clusterer(m_configuration, m_stats);


        auto DataFile = file::open(m_configuration.pFileName);
        auto RootGroup = DataFile.root();
        auto Dataset = RootGroup.get_dataset("srs_hits");
        dataspace::Simple Dataspace(Dataset.dataspace());
        std::vector<Gem::Readout> AllElements(Dataspace.size());

        Dataset.read(AllElements);
        /*
        std::sort(AllElements.begin(), AllElements.end(), [](const Readout& lhs, const Readout& rhs) {
            return lhs.srs_timestamp < rhs.srs_timestamp;
        });
        */
        Gem::CalibrationFile calfile(m_configuration.pCalFilename);
        for (auto RowData : AllElements ) {
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
            auto calib = calfile.getCalibration(RowData.fec, RowData.chip_id, RowData.channel);
            
            //User calibration added to analysis
            //The chiptime has to be recalculated in this case from bcid and tdc
            if(calib.time_offset != 0 || calib.time_slope != 1.0) {
                double bcTime = m_configuration.pBCTime_ns * RowData.bcid;
                double tdcTime = RowData.tdc * m_configuration.pTAC / 255;
                RowData.chiptime = bcTime + (m_configuration.pBCTime_ns - tdcTime - calib.time_offset)*calib.time_slope;
            }
            if(calib.adc_slope == 0) {
                std::cout << "Error in calibration file: adc_slope correction for fec " << RowData.fec << ", chip " 
                << RowData.chip_id << ", channel " << RowData.channel << " is 0!\nIs that intentional?" << std::endl;
            }
            RowData.adc = static_cast<uint16_t>((RowData.adc- calib.adc_offset)*calib.adc_slope);
                
            bool result = m_Clusterer->AnalyzeHits(
                static_cast<double>(RowData.srs_timestamp), 
                RowData.fec, RowData.chip_id, 
                RowData.channel, RowData.bcid, 
                RowData.tdc, RowData.adc, 
                RowData.over_threshold, RowData.chiptime);
            //std::cout << "hit: " << lines << ": " << (uint32_t) RowData.fec << ", " << (uint32_t) RowData.chip_id << ", " << RowData.srs_timestamp << ", " << RowData.channel << ", " << RowData.bcid << ", " << RowData.tdc << ", " << RowData.adc << ", " << RowData.over_threshold << ", " << RowData.chiptime << std::endl;
            lines++;
            if (result == false || (lines >= m_configuration.nHits && m_configuration.nHits > 0))
                break;
        }


        m_Clusterer->FinishAnalysis();

        delete m_Clusterer;
        timeEnd = std::chrono::system_clock::now();

        int elapsed_seconds = std::chrono::duration_cast < std::chrono::milliseconds > (timeEnd - timeStart).count();

        std::cout << "Analyzed " << lines << " lines. Finished computation in " << elapsed_seconds << " ms\n";
    }
    return 0;

}

