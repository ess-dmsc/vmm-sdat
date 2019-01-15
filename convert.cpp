#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include <chrono>

#include "Configuration.h"
#include "NMXClusterer.h"
#include "Readout.h"

using namespace hdf5;

int main(int argc, char**argv) {
    int lines = 0;
    std::chrono::time_point<std::chrono::system_clock> timeEnd, timeStart;

    Configuration m_configuration;

    if (m_configuration.parseCommandLine(argc, argv)) {

        timeStart = std::chrono::system_clock::now();

        NMXClusterer* m_NMXClusterer = new NMXClusterer(m_configuration.rootFilename, m_configuration.pBC, m_configuration.pTAC, m_configuration.pXChips, m_configuration.pYChips,
                                                        m_configuration.pADCThreshold, m_configuration.pMinClusterSize, m_configuration.pXYClusterSize, m_configuration.pDeltaTimeHits,
                                                        m_configuration.pMissingStripsCluster, m_configuration.pSpanClusterTime, m_configuration.pDeltaTimePlanes, m_configuration.analyzeChannels,
                                                        m_configuration.useUTPC, m_configuration.useHits);


        auto AnotherFile = file::open(m_configuration.fileName);
        auto RootGroup = AnotherFile.root();
        auto Dataset = RootGroup.get_dataset("srs_hits");
        dataspace::Simple Dataspace(Dataset.dataspace());
        std::vector<Readout> AllElements(Dataspace.size());

        Dataset.read(AllElements);
        /*
        std::sort(AllElements.begin(), AllElements.end(), [](const Readout& lhs, const Readout& rhs) {
            return lhs.srs_timestamp < rhs.srs_timestamp;
        });
        */
        for (auto RowData : AllElements) {

            int result = m_NMXClusterer->AnalyzeHits(static_cast<double>(RowData.srs_timestamp), RowData.fec, RowData.chip_id, RowData.channel, RowData.bcid, RowData.tdc, RowData.adc, RowData.over_threshold, RowData.chiptime);
            //std::cout << "hit: " << lines << ": " << (uint32_t) RowData.fec << ", " << (uint32_t) RowData.chip_id << ", " << RowData.srs_timestamp << ", " << RowData.channel << ", " << RowData.bcid << ", " << RowData.tdc << ", " << RowData.adc << ", " << RowData.over_threshold << ", " << RowData.chiptime << std::endl;
            lines++;
            if (result == -1)
                break;
        }


        m_NMXClusterer->PrintStats();

        delete m_NMXClusterer;
        timeEnd = std::chrono::system_clock::now();

        int elapsed_seconds = std::chrono::duration_cast < std::chrono::milliseconds > (timeEnd - timeStart).count();

        std::cout << "Analyzed " << lines << " lines. Finished computation in " << elapsed_seconds << " ms\n";
    }
    return 0;

}

