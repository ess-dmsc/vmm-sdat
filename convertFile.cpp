#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include <chrono>

#include "Configuration.h"
#include "Clusterer.h"
#include "Readout.h"

using namespace hdf5;

int main(int argc, char**argv) {
    int lines = 0;
    std::chrono::time_point<std::chrono::system_clock> timeEnd, timeStart;

    Configuration m_configuration;
    Statistics m_stats;

    if (m_configuration.ParseCommandLine(argc, argv)) {
        m_configuration.CreateMapping();
        m_stats.CreateStats(m_configuration);
        timeStart = std::chrono::system_clock::now();
       
        Clusterer* m_Clusterer = new Clusterer(m_configuration, m_stats);


        auto DataFile = file::open(m_configuration.pFileName);
        auto RootGroup = DataFile.root();
        auto Dataset = RootGroup.get_dataset("srs_hits");
        dataspace::Simple Dataspace(Dataset.dataspace());
        std::vector<Readout> AllElements(Dataspace.size());

        Dataset.read(AllElements);
        /*
        std::sort(AllElements.begin(), AllElements.end(), [](const Readout& lhs, const Readout& rhs) {
            return lhs.srs_timestamp < rhs.srs_timestamp;
        });
        */
        for (auto RowData : AllElements ) {

            bool result = m_Clusterer->AnalyzeHits(static_cast<double>(RowData.srs_timestamp), RowData.fec, RowData.chip_id, RowData.channel, RowData.bcid, RowData.tdc, RowData.adc, RowData.over_threshold, RowData.chiptime);
            //std::cout << "hit: " << lines << ": " << (uint32_t) RowData.fec << ", " << (uint32_t) RowData.chip_id << ", " << RowData.srs_timestamp << ", " << RowData.channel << ", " << RowData.bcid << ", " << RowData.tdc << ", " << RowData.adc << ", " << RowData.over_threshold << ", " << RowData.chiptime << std::endl;
            lines++;
            if (result == false || (lines >= m_configuration.nHits && m_configuration.nHits > 0))
                break;
        }


        m_Clusterer->PrintStats();

        delete m_Clusterer;
        timeEnd = std::chrono::system_clock::now();

        int elapsed_seconds = std::chrono::duration_cast < std::chrono::milliseconds > (timeEnd - timeStart).count();

        std::cout << "Analyzed " << lines << " lines. Finished computation in " << elapsed_seconds << " ms\n";
    }
    return 0;

}

