#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include "Configuration.h"


bool Configuration::printUsage(const std::string & errorMessage, char* argv)
{
    if (argv != nullptr) {
        std::cout << "\nERROR: " << errorMessage << ": " << argv << std::endl;
    } else {
        std::cout << "\nERROR: " << errorMessage << std::endl;
    }
    printf("\nUsages:\n");
    printf("analyse raw data:\n\t./convert -f ./data/a00010.h5 -x 1,0,1,1 -y 2,14,2,15 "
           "-bc 20 -tac 100 -th 0 -cs 3 -cxys 6 -dt 200 -mst 2 -spc 500 -dp 200 -cha 0 -utpc 1 -hits 1 ");

    printf("\nFlags:\n");
    printf("-f: h5 data file with the extension .h5\n\tThe data file was created by ESS tool.\n");
    printf("-x: mapping of chips, list of fecs and chips in x direction separated by comma (fec,chip, fec,chip etc) \n\n");
    printf("-y: mapping of chips, list of fecs and chips in y direction separated by comma (fec,chip, fec,chip etc) \n\n");
    printf("-bc: bunch crossing clock. Optional argument (default 20 MHz)\n\n");
    printf("-tac: tac slope. Optional argument (default 100 ns)\n\n");
    printf("-th: threshold value in ADC counts. Optional argument (default 0)\n\n");
    printf("-cs: minimum cluster size. Optional argument (default 3)\n\n");
    printf("-cxys: minimum cluster size in x and y together. Optional argument (default 6)\n\n");
    printf("-dt: maximum time difference between strips in time sorted vector. Optional argument (default 200)\n\n");
    printf("-mst: maximum missing strips in strip sorted vector. Optional argument (default 2)\n\n");
    printf("-spc: maximum time span of cluster in one dimension (determined by drift size and speed). Optional argument (default 500)\n\n");
    printf("-dp: maximum time between matched clusters in x and y. Optional argument (default 200)\n\n");
    printf("-cha: analyze TDC, BCID and ADC of all channels. Takes a long time. Optional argument (default 0)\n\n");
    printf("-utpc: use uTPC method to determine timestamp of cluster. Optional argument (default 1)\n\n");
    printf("-hits: store not only clusters but all hits (a hit is a VMM3 channel over threshold). Creates large files. Optional argument (default 1)\n\n");

    return false;



}

bool Configuration::parseCommandLine(int argc, char**argv)
{
    if (argc == 1 || argc % 2 == 0) {
        return printUsage("Wrong number of arguments!", argv[argc - 1]);
    }
    for (int i = 1; i < argc; i += 2) {
        if (strncmp(argv[i], "-f", 2) == 0) {
            fFound = true;
            fileName = argv[i + 1];
            std::cout << fileName << std::endl;
        } else if (strncmp(argv[i], "-bc", 3) == 0) {
            bcFound = true;
            pBC = atoi(argv[i + 1]);
        } else if (strncmp(argv[i], "-x", 2) == 0) {
            xFound = true;
            std::string xString = argv[i + 1];
            std::string delims = ",";
            size_t lastOffset = 0;
            pXChips.clear();
            int n = 0;
            int det = 0;
            int fec = 0;
            int vmm = 0;
            while (true) {
                size_t offset = xString.find_first_of(delims, lastOffset);
                if (n % 3 == 0) {
                    det = atoi(xString.substr(lastOffset, offset - lastOffset).c_str());
                } else if(n % 3 == 1) {
                    fec = atoi(xString.substr(lastOffset, offset - lastOffset).c_str());
                }
                else {
                   vmm = atoi(xString.substr(lastOffset, offset - lastOffset).c_str());
                   pXChips.emplace_back(std::make_tuple(det, fec, vmm));
                }
                n++;
                if (offset == std::string::npos) {
                    break;
                } else {
                    lastOffset = offset + 1; // add one to skip the delimiter
                }

            }
            if (pXChips.size() != (int) (n / 3)) {
                return printUsage("Wrong number of Detectors, FECs and VMMs in x!", argv[i]);
            }

        } else if (strncmp(argv[i], "-y", 2) == 0) {
            yFound = true;
            std::string yString = argv[i + 1];
            std::string delims = ",";
            size_t lastOffset = 0;
            pYChips.clear();
            int n = 0;
            int det = 0;
            int fec = 0;
            int vmm = 0;
            while (true) {
                size_t offset = yString.find_first_of(delims, lastOffset);
                if (n % 3 == 0) {
                    det = atoi(yString.substr(lastOffset, offset - lastOffset).c_str());
                } else if(n % 3 == 1) {
                    fec = atoi(yString.substr(lastOffset, offset - lastOffset).c_str());
                }
                else {
                   vmm = atoi(yString.substr(lastOffset, offset - lastOffset).c_str());
                   pYChips.emplace_back(std::make_tuple(det, fec, vmm));
                }
                n++;
                if (offset == std::string::npos) {
                    break;
                } else {
                    lastOffset = offset + 1; // add one to skip the delimiter
                }
            }
            if (pYChips.size() != (int) (n / 3)) {
                return printUsage("Wrong number of Detectors, FECs and VMMs in y!", argv[i]);
            }

        } else if (strncmp(argv[i], "-tac", 4) == 0) {
            tacFound = true;
            pTAC = atoi(argv[i + 1]);
        } else if (strncmp(argv[i], "-th", 3) == 0) {
            thresholdFound = true;
            pADCThreshold = atoi(argv[i + 1]);
        } else if (strncmp(argv[i], "-cs", 3) == 0) {
            clusterSizeFound = true;
            pMinClusterSize = atoi(argv[i + 1]);
        } else if (strncmp(argv[i], "-cxys", 5) == 0) {
            clusterSizeXYFound = true;
            pXYClusterSize = atoi(argv[i + 1]);
        } else if (strncmp(argv[i], "-dt", 3) == 0) {
            deltaTimeHitsFound = true;
            pDeltaTimeHits = atoi(argv[i + 1]);
        } else if (strncmp(argv[i], "-mst", 4) == 0) {
            missingStripsClusterFound = true;
            pMissingStripsCluster = atoi(argv[i + 1]);
        } else if (strncmp(argv[i], "-spc", 4) == 0) {
            spanClusterTimeFound = true;
            pSpanClusterTime = atoi(argv[i + 1]);
        } else if (strncmp(argv[i], "-dp", 3) == 0) {
            deltaTimePlanesFound = true;
            pDeltaTimePlanes = atoi(argv[i + 1]);
        } else if (strncmp(argv[i], "-cha", 4) == 0) {
            analyzeChannelsFound = true;
            analyzeChannels = atoi(argv[i + 1]);
        } else if (strncmp(argv[i], "-utpc", 5) == 0) {
            useUTPCFound = true;
            useUTPC = atoi(argv[i + 1]);
        } else if (strncmp(argv[i], "-hits", 5) == 0) {
            useHitsFound = true;
            useHits = atoi(argv[i + 1]);
        } else {
            return printUsage("Wrong type of argument!", argv[i]);
        }
    }

    if (!fFound) {
        return printUsage("Data file has to be loaded with -f data.h5!", nullptr);

    }

    if (fFound && fileName.find(".h5") == std::string::npos) {
        return printUsage("Wrong extension: .h5 file required!", nullptr);
    }

    channels_x = 64 * pXChips.size();
    channels_y = 64 * pYChips.size();
    rootFilename = fileName;
    if (rootFilename.find(".h5") != std::string::npos) {
        rootFilename.replace(rootFilename.size() - 3, rootFilename.size(), "");
    }
    std::string strParams;

    if (useHits) {
        strParams += "_HITS";
    }
    strParams += ".root";
    rootFilename = rootFilename + "_bc_" + std::to_string(pBC) + "_tac_" + std::to_string((int) pTAC) + "_cxys" + std::to_string((int) pXYClusterSize) + "_cs" + std::to_string((int) pMinClusterSize) + "_dt" + std::to_string((int) pDeltaTimeHits) + "_mst" + std::to_string((int) pMissingStripsCluster) + "_spc" + std::to_string((int) pSpanClusterTime) + "_dp" + std::to_string((int) pDeltaTimePlanes) + strParams;

    isConfigured = true;
    return true;

}
