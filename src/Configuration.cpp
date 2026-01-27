/***************************************************************************
**  vmm-sdat
**  Data analysis program for VMM3a data
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see http://www.gnu.org/licenses/.
**
****************************************************************************
**  Contact: dorothea.pfeiffer@cern.ch
**  Date: 12.10.2025
**  Version: 1.0.0
****************************************************************************
**
**  vmm-sdat
**  Configuration.cpp
**
****************************************************************************/

#include <cmath>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <regex>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>

#include "log.h"
#include "Configuration.h"


using json = nlohmann::json;

Configuration::Configuration()
{
  corryvreckan::Log::setSection("Configuration");
}

bool Configuration::PrintUsage(const std::string & errorMessage, char * argv) {
  std::cout << "\nUsage:" << std::endl;
  std::cout <<
    "./convertFile -f ../../data.pcapng " <<
    "-log INFO -geo example_geometry_pad.json " <<
  "-bc 40 -tac 60 -th [0,20] -cs [1,2] -ccs [3,4] -dt [200,300] -mst "
  "[1,0] -spc [500,400] " <<
  "-dp [200,250] -coin center-of-mass -crl [0.75,0.5] -cru [3.0,5] "
  "-save [[1],[1,2],[1,2]] -df SRS" <<
  std::endl;

  std::cout << "\n\nFlags:\n" << std::endl;
  std::cout << "-f:     PCAPNG file created with wireshark or tcdump "
  "(*.pcapng), coming from the SRS or the ESS readout.\n" <<
  std::endl;
  std::cout << "-log:   Corryvreckan style log, set to TRACE,DEBUG, INFO, STATUS, WARNING, ERROR, FATAL." <<
  std::endl;
  std::cout << "-geo:   The detector "
  "geometry has to be defined in a JSON file.\n" <<
  "        Examples of geometry files (for strips and pads) "
  "are in the run folder.\n" <<
  std::endl;

  std::cout <<
    "-bc:    bunch crossing clock. Optional argument (default 40 MHz).\n" <<
    std::endl;
  std::cout << "-tac:   tac slope. Optional argument (default 60 ns).\n" <<
    std::endl;
  std::cout <<
    "-t0:    Time correction in ns for each vmm in the format "
  "[[fec0, vmm0, correction0],[fec1, vmm1, correction1]]. "
  "The correction value is subtracted from all timestamps.\n"
  "        If instead of a number the word 'run' is put as correction,"
  " the first timestamp of the run is used as correction. "
  "Optional argument for SRS data format (default 0)\n " <<
  std::endl;
  std::cout <<
    "-th:    threshold value in ADC counts. Optional argument (default 0, "
  "if -1, only hits with over threshold flag 1 are expected, one value "
  "per detector).\n" <<
  std::endl;
  std::cout << "-cs:    minimum cluster size per plane. Optional argument "
  "(default 1), one value per detector.\n" <<
  std::endl;
  std::cout << "-ccs:   minimum cluster size in plane 0 and plane 1 together. "
  "Optional argument (default 2), one value per detector.\n" <<
  std::endl;
  std::cout <<
    "-dt:    maximum time difference between strips in time sorted "
  "vector. Optional argument (default 200), one value per detector.\n" <<
  std::endl;
  std::cout << "-mst:   maximum missing strips in strip sorted vector. "
  "Optional argument (default 0), one value per detector.\n" <<
  std::endl;
  std::cout <<
    "-mp0:   maximum missing pads in dimension 0 in pad sorted "
  "vector. Optional argument (default 0), one value per detector.\n" <<
  std::endl;
  std::cout <<
    "-mp1:   maximum missing pads in dimension 1 in pad sorted "
  "vector. Optional argument (default 0), one value per detector.\n" <<
  std::endl;
  std::cout <<
    "-spc:   maximum time span of cluster in one dimension (determined by "
  "drift size and speed). Optional argument (default 500), one value "
  "per detector.\n" <<
  std::endl;
  std::cout << "-dp:    maximum time between matched clusters in x and y. "
  "Optional argument (default 200), one value per detector.\n" <<
  std::endl;
  std::cout << "-coin:  Valid clusters normally occur at the same time in "
  "plane 0 and plane 1 of a detctor. The parameter -dp determines "
  "the permitted time difference between the planes.\n" <<
  "        The time can be calculated with the center-of-mass "
  "algorithm (center-of-mass), the uTPC method (utpc) or the "
  "center-of-mass squared method (charge2).\n" <<
  "        Optional argument (default center-of-mass).\n" <<
  std::endl;
  std::cout << "-algo:  Select with algorithm is used in pos_algo and "
  "time_algo field in clusters" <<
  std::endl;
  std::cout << "        0: no algorithm" << std::endl;
  std::cout << "        1: Gaussian fit" << std::endl;
  std::cout << "        2: COG including only over Threshold hits" << std::endl;
  std::cout << "        3: COG2 including only over Threshold hits" <<
    std::endl;
  std::cout << "        4: position and time of largest ADC" << std::endl;
  std::cout << "        5: trigger pattern (NIP box), the trigger pattern is "
  "stored as integer in time_algo2" <<
  std::endl;
  std::cout << "           The vmm that is connected to the NIP box has to be "
  "defined as plane 2 of the detector." <<
  std::endl;
  std::cout << "           The channels of the VMM have to be mapped to the "
  "strips in the form:" <<
  std::endl;
  std::cout << "           channel representing bit 0 = strip 0, channel for "
  "bit 1  = strip 1 and so on." <<
  std::endl;
  std::cout << "        6: utpc with COG" << std::endl;
  std::cout << "        7: utpc with COG2" << std::endl;

  std::cout << "-crl:   Valid clusters normally have the same amount of charge "
  "in both detector planes (ratio of charge plane 0/charge plane "
  "1 is 100\% or 1, one value per detector.\n" <<
  "        Depending on the readout, the charge sharing can be "
  "different, e.g. in a standard GEM strip readout the total "
  "charge is divided 60/40 between plane 0/ plane 1\n" <<
  "        With -crl one sets the lower threshold for the "
  "plane0/plane1 charge ratio. Optional argument (default 0.5)" <<
  std::endl;
  std::cout << "-cru:   With -cru one sets the upper threshold for the "
  "plane0/plane1 charge ratio. Optional argument (default 2), one "
  "value per detector.\n" <<
  std::endl;
  std::cout <<
    "-hm:    High-multiplicity matching mode (values 0 or 1).\n"
  "        During the normal matching vmm-sdat searches for each "
  "cluster in one "
  "plane the best match in the other plane, based on the minimum time "
  "difference between the clusters.\n"
  "        Before considering a cluster in the other plane as match "
  "candidate, "
  "conditions like cluster size, and charge sharing are checked.\n"
  "        In high-multiplicity mode, each combination of clusters in "
  "the two "
  "planes that fulfils the conditions is stored as detector_cluster.\n"
  "        That means each plane cluster can appear several times as "
  "part of a "
  "detector_cluster.\n" <<
  std::endl;
  std::cout << "-save:  select which data to store in root file. Input is a "
  "list of lists of detectors, e.g. [[1,2],[1,2],[1,2,3]]." <<
  std::endl;
  std::cout << "        first list : detectors for which to write the hits "
  "(hit is a VMM3a channel over threshold)" <<
  std::endl;
  std::cout << "        second list : clusters plane" << std::endl;
  std::cout << "        third list : clusters detector" << std::endl;
  std::cout << "        Examples:" << std::endl;
  std::cout << "            [[1,2],[],[]]: hits for detectors 1 and 2 only" <<
    std::endl;
  std::cout << "            [[],[],[1,2]]: clusters detector for detector 1 "
  "and 2 only" <<
  std::endl;
  std::cout << "            [[2],[1],[1]]: hits for detector 2, clusters "
  "plane, clusters detector for detector 1 \n" <<
  std::endl;
  std::cout << "-json:  create a json file of the detector images. Optional "
  "argument (default 1).\n" <<
  std::endl;
  std::cout << "-n:     number of hits to analyze. Optional argument (default "
  "0, i.e. all hits).\n" <<
  std::endl;
  std::cout << "-stats: Show statistics of the run (default 0, do not show any "
  "stats).\n" <<
  std::endl;
  std::cout << "-cal:   Name of the calibration file. A calibration file is a "
  "JSON file containing an ADC and/or time correction in the form "
  "of a slope and an offset value. Optional parameter.\n" <<
  std::endl;
  std::cout <<
    "-df:    Data format: The pcap files can have different data formats, " <<
    "depending on the firmware on the FEC.\n" <<
    "        SRS (default): FEC card, and time stamps and offsets are "
  "sent in markers, offset is interpreted as signed number and valid "
  "offsets goes from -1 to 15\n" <<
  "        TRIG: Triggered firmware, event_counter in hit tree" <<
  std::endl;
  std::cout <<
    "-cahi:    Calibration histograms: If a calibration file is used, " <<
    "histograms of the calibrated and uncalibrated adc and time can be "
  "produced.\n" <<
  "        With the help of these histograms the effect of the "
  "calibration can "
  "be checked." <<
  "        default: 0" << std::endl;
  std::cout << "-info:  Additional info the user wants to be added to the end "
  "of the newly created file name.\n" <<
  std::endl;
  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
  "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" <<
  std::endl;
  if (argv != nullptr) {
    std::cout << "ERROR: " << errorMessage << ": " << argv << std::endl;
  } else {
    std::cout << "ERROR: " << errorMessage << std::endl;
  }
  std::cout << "\nFor meaning of the flags and the correct usage of "
  "convertFile, please see above!" <<
  std::endl;
  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
  "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" <<
  std::endl;

  return false;
}

bool Configuration::ParseCommandLine(int argc, char ** argv) {

  if (argc == 1 || argc % 2 == 0) {
    return PrintUsage("Wrong number of arguments!", argv[argc - 1]);
  }
  for (int i = 1; i < argc; i += 2) {
    if (strncmp(argv[i], "-f", 2) == 0) {
      fFound = true;
      pFileName = argv[i + 1];
    } else if (strncmp(argv[i], "-log", 4) == 0) {
      pLogLevel = argv[i + 1];
      auto it = find(pLogLevels.begin(), pLogLevels.end(), pLogLevel);
	  if (it == pLogLevels.end()) {
	  	pLogLevel = "INFO";
	  }
    } else if (strncmp(argv[i], "-info", 5) == 0) {
      pInfo = argv[i + 1];
    } else if (strncmp(argv[i], "-t0", 3) == 0) {
      std::string t0String = argv[i + 1];
      char removeChars[] = "[]";
      for (unsigned int i = 0; i < strlen(removeChars); ++i) {
        t0String.erase(
          std::remove(t0String.begin(), t0String.end(), removeChars[i]),
          t0String.end());
      }
      std::string delims = ",";
      size_t lastOffset = 0;
      pFecVMM_time0.clear();

      int n = 0;
      std::string correction = "";
      int fec = 0;
      int vmm = 0;
      while (true) {
        size_t offset = t0String.find_first_of(delims, lastOffset);
        if (n % 3 == 0) {
          fec = atoi(t0String.substr(lastOffset, offset - lastOffset).c_str());
        } else if (n % 3 == 1) {
          vmm = atoi(t0String.substr(lastOffset, offset - lastOffset).c_str());
        } else {
          correction = t0String.substr(lastOffset, offset - lastOffset);
          pFecVMM_time0[std::make_pair(fec, vmm)] = correction;
        }
        n++;
        if (offset == std::string::npos) {
          break;
        } else {
          lastOffset = offset + 1; // add one to skip the delimiter
        }
      }
      if (pFecVMM_time0.size() != (int)(n / 3)) {
        return PrintUsage("Wrong number of FECs, VMMs and correction!",
          argv[i]);
      }
    } else if (strncmp(argv[i], "-bc", 3) == 0) {
      pBC = atof(argv[i + 1]);
      // ESS SRS firmware has 44.444444 MHz clock
      if (pBC >= 44.4 && pBC <= 44.5) {
        pBCTime_ns = 22.5;
        pOffsetPeriod = 4096.0 * pBCTime_ns;
      }
      // 88.0525 MHz ESS clock, half of that is BC clock
      else if (pBC >= 44.0 && pBC <= 44.1) {
        pBCTime_ns = 22.7137219;
        pOffsetPeriod = 4096.0 * pBCTime_ns;
      } else {
        pBCTime_ns = (1000.0 / pBC);
        pOffsetPeriod = 4096.0 * 25.0;
      }
    } else if (strncmp(argv[i], "-geo", 4) == 0) {
      pGeometryFile = argv[i + 1];
      if (pGeometryFile.find(".json") == std::string::npos) {
        return PrintUsage("The geometry parameter -geo requires a JSON file!",
          nullptr);
      } 
    } else if (strncmp(argv[i], "-tac", 4) == 0) {
      pTAC = atof(argv[i + 1]);
    } else if (strncmp(argv[i], "-th", 3) == 0) {
      GetDetectorParameters(argv[i + 1], pADCThreshold);
    } else if (strncmp(argv[i], "-cs", 3) == 0) {
      GetDetectorParameters(argv[i + 1], pMinClusterSize);
    } else if (strncmp(argv[i], "-ccs", 4) == 0) {
      GetDetectorParameters(argv[i + 1], pCoincidentClusterSize);
    } else if (strncmp(argv[i], "-dt", 3) == 0) {
      GetDetectorParameters(argv[i + 1], pDeltaTimeHits);
    } else if (strncmp(argv[i], "-mst", 4) == 0 ||
      strncmp(argv[i], "-mp0", 4) == 0) {
      GetDetectorParameters(argv[i + 1], pMissingStripsClusterX);
    } else if (strncmp(argv[i], "-mp1", 4) == 0) {
      GetDetectorParameters(argv[i + 1], pMissingStripsClusterY);
    } else if (strncmp(argv[i], "-spc", 4) == 0) {
      GetDetectorParameters(argv[i + 1], pSpanClusterTime);
    } else if (strncmp(argv[i], "-dp", 3) == 0) {
      GetDetectorParameters(argv[i + 1], pDeltaTimePlanes);
    } else if (strncmp(argv[i], "-crl", 4) == 0) {
      GetDetectorParameters(argv[i + 1], pChargeRatioLower);
    } else if (strncmp(argv[i], "-cru", 4) == 0) {
      GetDetectorParameters(argv[i + 1], pChargeRatioUpper);
    } else if (strncmp(argv[i], "-cahi", 5) == 0) {
      if (atoi(argv[i + 1]) == 1) {
        calibrationHistogram = true;
      } else {
        calibrationHistogram = false;
      }
    } else if (strncmp(argv[i], "-save", 5) == 0) {
      std::string parameterString = argv[i + 1];
      char removeChars[] = " ";
      for (unsigned int i = 0; i < strlen(removeChars); ++i) {
        parameterString.erase(std::remove(parameterString.begin(),
            parameterString.end(),
            removeChars[i]),
          parameterString.end());
      }

      parameterString =
        std::regex_replace(parameterString, std::regex("\\[\\["), "; ");
      parameterString =
        std::regex_replace(parameterString, std::regex("\\]\\]"), ";");
      parameterString =
        std::regex_replace(parameterString, std::regex("\\],\\["), "| ");
      parameterString.erase(
        std::remove(parameterString.begin(), parameterString.end(), ';'),
        parameterString.end());
      std::vector < std::string > vTokens;
      std::string token;
      std::istringstream tokenStream(parameterString);
      int n = 0;
      pSaveWhat = 0;
      while (std::getline(tokenStream, token, '|')) {
        vTokens.push_back(token);
        if (token != " ") {
          if (n == 0) {
            pSaveWhat = 1;
          } else if (n == 1) {
            pSaveWhat += 10;
          } else if (n == 2) {
            pSaveWhat += 100;
          }
        }
        n++;
      }
      if (vTokens.size() != 3) {
        return PrintUsage(
          "The -save parameter accepts only a list of lists of detectors "
          "(numbers) in the format [[],[1,2],[1,2,3]]!",
          nullptr);
      }
      pSaveHits.clear();
      pSaveClustersPlane.clear();
      pSaveClustersDetector.clear();
      n = 0;
      pHighMultiplicity = false;

      for (auto & s: vTokens) {
        std::vector < std::string > v;
        std::string token;
        for (unsigned int i = 0; i < strlen(removeChars); ++i) {
          s.erase(std::remove(s.begin(), s.end(), removeChars[i]), s.end());
        }
        std::istringstream tokenStream(s);
        while (std::getline(tokenStream, token, ',')) {
          if (token != " ") {
            if (token.find_first_not_of("0123456789") != std::string::npos) {
              return PrintUsage(
                "The -save parameter accepts only a list of lists of "
                "detectors (numbers) in the format [[],[1,2],[1,2,3]]!",
                nullptr);
            }
          }
          if (n == 0 && token != " ") {
            pSaveHits.push_back(std::stoi(token));
          } else if (n == 1 && token != " ") {
            pSaveClustersPlane.push_back(std::stoi(token));
          } else if (n == 2 && token != " ") {
            pSaveClustersDetector.push_back(std::stoi(token));
          }
        }
        n++;
      }
    } else if (strncmp(argv[i], "-n", 2) == 0) {
      nHits = atoi(argv[i + 1]);
    } else if (strncmp(argv[i], "-cal", 4) == 0) {
      pCalFilename = argv[i + 1];
      useCalibration = true;
    } else if (strncmp(argv[i], "-json", 5) == 0) {
      createJSON = atoi(argv[i + 1]);
    } else if (strncmp(argv[i], "-coin", 5) == 0) {
      pConditionCoincidence = "center-of-mass";
      if (strncmp(argv[i + 1], "utpc", 4) == 0 ||
        strncmp(argv[i + 1], "charge2", 7) == 0) {
        pConditionCoincidence = argv[i + 1];
      }
    } else if (strncmp(argv[i], "-algo", 5) == 0) {
      pAlgo = atoi(argv[i + 1]);
    } else if (strncmp(argv[i], "-hm", 3) == 0) {
      if (atoi(argv[i + 1]) == 1) {
        pHighMultiplicity = true;
      }
    } else if (strncmp(argv[i], "-df", 3) == 0) {
      pDataFormat = argv[i + 1];
      std::vector < std::string > v_valid_values = {
        "SRS",
        "TRG",
        "srs",
        "trg"
      };
      auto searchValid =
        std::find(v_valid_values.begin(), v_valid_values.end(), pDataFormat);
      if (searchValid == v_valid_values.end()) {
        return PrintUsage("The data format parameter -df accepts only the "
          "values SRS, TRG!",
          nullptr);
      }
    } else {
      return PrintUsage("Wrong type of argument!", argv[i]);
    }
  }
  if (!fFound) {
    return PrintUsage("Data file has to be loaded with -f data.h5!", nullptr);
  }

  if (pFileName.find(".pcapng") == std::string::npos) {
    return PrintUsage(
      "Wrong extension: .pcapng file required for data files!",
      nullptr);
  }

  if (useCalibration && pCalFilename.find(".json") == std::string::npos) {
    return PrintUsage("Wrong extension: .json file required for calibration!",
      nullptr);
  }
  if (!vmmsFound && pGeometryFile.find(".json") == std::string::npos) {
    return PrintUsage("Detectors, planes, fecs and VMMs have to be defined, or "
      "a geometry file loaded!",
      nullptr);
  }
  LOG(STATUS) << "Analyzing " << pFileName << " ...";
  pRootFilename = pFileName;
  if (pRootFilename.find(".pcapng") != std::string::npos) {
    pRootFilename.replace(pRootFilename.size() - 7, pRootFilename.size(), "");
  }

  time_t ttime = time(0);
  tm * local_time = localtime( & ttime);
  auto t = std::time(nullptr);
  auto tm = * std::localtime( & t);
  std::stringstream sTime;
  sTime << std::put_time( & tm, "%Y%m%d%H%M%S");

  std::string strParams = "_";
  strParams += sTime.str();

  if (pInfo.length() > 0) {
    strParams += "_";
    strParams += pInfo;
  }
  strParams += ".root";
  pRootFilename = pRootFilename + strParams;
  return true;
}



bool Configuration::GetDetectorPlane(std::pair < uint8_t, uint8_t > dp) {
  auto searchDetPlane = p_DetPlane_idx.find(dp);
  if (searchDetPlane == p_DetPlane_idx.end()) {
    return false;
  }
  return true;
}

bool Configuration::CreateMapping() {
  if (pGeometryFile.find(".json") == std::string::npos) {
    return PrintUsage("Geometry definiton missing! Define geometry by -geo!",
      nullptr);
  }
  for (int f = 0; f < NUMFECS; f++) {
    for (int v = 0; v < 16; v++) {
      pDetectors[f][v] = -1;
      pPlanes[f][v] = -1;
      if (f == 0) {
        pIsPads[v] = false;
      }
      for (int ch = 0; ch < 64; ch++) {
        pPositions0[f][v][ch] = -1;
        pPositions1[f][v][ch] = -1;
      }
    }
  }
  
    std::ifstream t(pGeometryFile);
    std::string jsonstring((std::istreambuf_iterator < char > (t)),
      std::istreambuf_iterator < char > ());
    if (!t.good()) {
      return PrintUsage("Invalid JSON file format!", nullptr);
    }
    nlohmann::json Root;
    try {
      Root = nlohmann::json::parse(jsonstring);
    } catch (...) {
      throw std::runtime_error("Invalid Json in geometry file.");
    }
    pVMMs.clear();
    try {
      auto vmm_geos = Root["vmm_geometry"];
      for (auto & geo: vmm_geos) {
        auto fec = geo["fec"].get < uint16_t > ();
        auto vmm = geo["vmm"].get < uint8_t > ();
        auto detector = geo["detector"].get < uint8_t > ();
        auto strips0 = geo["id0"];
        auto strips1 = geo["id1"];
        uint8_t plane = 0;
        pIsPads[detector] = false;
        if (strips0.size() != 64 ||
          (strips1.size() > 0 && strips1.size() < 64)) {
          throw std::runtime_error(
            "Wrong lengths of id0 or id1 arrays in geometry file.");
        } else if (strips1.size() == 64) {
          pIsPads[detector] = true;
        } 
        plane = geo["plane"].get < uint8_t > ();
        auto searchDetPlane =
          p_DetPlane_idx.find(std::make_pair(detector, plane));
        if (searchDetPlane == p_DetPlane_idx.end()) {
          // Add det/plane pair to the list and set index
          p_DetPlane_idx.emplace(
            std::make_pair(std::make_pair(detector, plane), 0));
        }
        auto searchTuple =
          std::find(std::begin(pVMMs), std::end(pVMMs),
            std::make_tuple(detector, plane, fec, vmm));
        if (searchTuple == pVMMs.end()) {
          pVMMs.emplace_back(std::make_tuple(detector, plane, fec, vmm));
          auto searchTuple = pChannels.find(std::make_pair(detector, plane));
          int strips = 0;
          for (size_t ch = 0; ch < strips0.size(); ch++) {
            int s0 = strips0[ch].get < int > ();
            if (s0 > -1) {
              strips++;
            }
          }
          if (searchTuple == pChannels.end()) {
            pChannels[std::make_tuple(detector, plane)] = strips;
          } else {
            pChannels[std::make_tuple(detector, plane)] += strips;
          }
        }
        
        auto searchFec = std::find(std::begin(pFecs), std::end(pFecs), fec);
        if (searchFec == pFecs.end()) {
          pFecs.push_back(fec);
        }

        auto searchDet = pDets.find(detector);
        if (searchDet == pDets.end()) {
          pDets.emplace(detector, pDets.size());
          pChannels0[detector] = 0;
          pChannels1[detector] = 0;
        }

        bool found = false;
        // Search whether there is a new det/plane/fec combination
        for (auto
          const & searchDetPlaneFec: pDetectorPlane_Fec) {
          if (searchDetPlaneFec.first == std::make_pair(detector, plane) &&
            searchDetPlaneFec.second == fec) {
            found = true;
            break;
          }
        }
        if (found == false) {
          pDetectorPlane_Fec.emplace(
            std::make_pair(std::make_pair(detector, plane), fec));
        }

        // Search whether there is a new fec/chip combination
        auto searchFecChip =
          pFecChip_DetectorPlane.find(std::make_pair(fec, vmm));
        if (searchFecChip == pFecChip_DetectorPlane.end()) {
          pDetectors[fec][vmm] = (int) detector;
          pPlanes[fec][vmm] = (int) plane;
          for (size_t ch = 0; ch < strips0.size(); ch++) {
            int s0 = strips0[ch].get < int > ();
            pPositions0[fec][vmm][ch] = s0;
            if (pIsPads[detector]) {
              int s1 = strips1[ch].get < int > ();
              pPositions1[fec][vmm][ch] = s1;
              if (s0 > pChannels0[detector]) {
                pChannels0[detector] = s0;
              }
              if (s1 > pChannels1[detector]) {
                pChannels1[detector] = s1;
              }
            }
          }

          // Add the new fec/chip pair to the list
          pFecChip_DetectorPlane.emplace(std::make_pair(
            std::make_pair(fec, vmm), std::make_pair(detector, plane)));
        }
      }
    } catch (const std::exception & exc) {
      throw std::runtime_error("Invalid json in geometry file.");
    }
 

  bool ret = CheckDetectorParameters("pMinClusterSize", pMinClusterSize);
  if (ret == false) {
    return false;
  }

  ret =
    CheckDetectorParameters("pCoincidentClusterSize", pCoincidentClusterSize);
  if (ret == false) {
    return false;
  }

  ret = CheckDetectorParameters("pDeltaTimeHits", pDeltaTimeHits);
  if (ret == false) {
    return false;
  }

  ret =
    CheckDetectorParameters("pMissingStripsClusterX", pMissingStripsClusterX);
  if (ret == false) {
    return false;
  }

  ret =
    CheckDetectorParameters("pMissingStripsClusterY", pMissingStripsClusterY);
  if (ret == false) {
    return false;
  }

  ret = CheckDetectorParameters("pSpanClusterTime", pSpanClusterTime);
  if (ret == false) {
    return false;
  }

  ret = CheckDetectorParameters("pDeltaTimePlanes", pDeltaTimePlanes);
  if (ret == false) {
    return false;
  }

  ret = CheckDetectorParameters("pChargeRatioLower", pChargeRatioLower);
  if (ret == false) {
    return false;
  }

  ret = CheckDetectorParameters("pChargeRatioUpper", pChargeRatioUpper);
  if (ret == false) {
    return false;
  }

  ret = CheckDetectorParameters("pADCThreshold", pADCThreshold);
  if (ret == false) {
    return false;
  }

  return true;
}

void Configuration::GetDetectorParameters(std::string input,
  std::vector < double > & v) {
  v.clear();
  std::string parameterString = input;
  char removeChars[] = "[ ]";
  for (unsigned int i = 0; i < strlen(removeChars); ++i) {
    parameterString.erase(std::remove(parameterString.begin(),
        parameterString.end(), removeChars[i]),
      parameterString.end());
  }
  std::string token;
  std::istringstream tokenStream(parameterString);
  while (std::getline(tokenStream, token, ',')) {
    v.push_back(std::stof(token));
  }
}

bool Configuration::CheckDetectorParameters(std::string name,
  std::vector < double > & v) {
  if (v.size() == 1) {
    float val = v[0];
    for (int n = 1; n < pDets.size(); n++) {
      v.push_back(val);
    }
  } else {
    if (v.size() != pDets.size()) {
      return PrintUsage(
        "Wrong number of parameters, one per detector or one for all!",
        (char * ) name.c_str());
    }
  }
  return true;
}