#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>


class Configuration  {
public:
        Configuration() = default;
        ~Configuration() = default;
        bool parseCommandLine(int argc, char**argv);

        bool printUsage(const std::string & errorMessage, char* argv);
        int channels_x = 0;
        int channels_y = 0;
        std::string fileName = "";
        std::string rootFilename = "";
        bool isConfigured = false;
        bool fFound = false;
        bool tacFound = false;
        bool bcFound = false;
        bool xFound = false;
        bool yFound = false;
        bool thresholdFound = false;
        bool clusterSizeFound = false;
        bool clusterSizeXYFound = false;
        bool deltaTimeHitsFound = false;
        bool missingStripsClusterFound = false;
        bool spanClusterTimeFound = false;
        bool deltaTimePlanesFound = false;
        bool analyzeChannelsFound = false;

        bool useUTPCFound = false;
        bool useHitsFound = false;

        std::vector<std::tuple<uint8_t, uint8_t, uint8_t>> pXChips { {1, 1, 6 }, { 1, 1, 7 }, {1, 1, 0 }, {1, 1, 1 } };
        std::vector<std::tuple<uint8_t, uint8_t, uint8_t>> pYChips { { 1, 1, 14 }, { 1, 1, 15 }, { 1, 1, 4 }, {1, 1, 5 } };

        uint16_t pTAC = 100;
        uint16_t pBC = 20;
        uint16_t pADCThreshold = 0;
        uint16_t pMinClusterSize = 3;
        uint16_t pXYClusterSize = 6;
        //Maximum time difference between strips in time sorted cluster (x or y)
        uint16_t pDeltaTimeHits = 200;
        //Number of missing strips in strip sorted cluster (x or y)
        uint16_t pMissingStripsCluster = 2;
        //Maximum time span for total cluster (x or y)
        uint16_t pSpanClusterTime = 500;
        //Maximum cluster time difference between matching clusters in x and y
        //Cluster time is either calculated with center-of-mass or uTPC method
        uint16_t pDeltaTimePlanes = 200;

        bool analyzeChannels = false;

        bool useHits = true;
        bool useUTPC = true;
        int lines = 0;

};






