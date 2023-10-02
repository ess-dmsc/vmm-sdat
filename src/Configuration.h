#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#define NUMFECS 385

class Configuration {
public:
  Configuration() = default;
  ~Configuration() = default;
  bool ParseCommandLine(int argc, char **argv);

  bool PrintUsage(const std::string &errorMessage, char *argv);
  bool CreateMapping();
  bool CalculateTransform();
  // bool GetAxes(std::pair<uint8_t, uint8_t> dp);
  bool GetDetectorPlane(std::pair<uint8_t, uint8_t> dp);
  //**************************************************************
  // BEGIN INPUT PARAMETERS
  //**************************************************************

  // The tuples for the VMMs are defined as follows:
  // detector (choose a number between 0 and 255)
  // plane (0 or 1)
  // fec (fecID set in firmware)
  // vmm (depends on connection of hybrid to FEC, FEC channel 1 equals VMMs 0
  // and 1, FEC channel 2 VMMs 2 and 3, FEC channel 8 VMMs 14 and 15) When
  // looking at the detector, the following conventions are used:
  //   - top side of the hybrids is visible (if the hybrids are mounted in the
  //   readout plane)
  //   - OR: Side of the Hirose connector is visible (if hybrids are mounted on
  //   the side of the detector)
  //   - plane 0 is at the bottom (HDMI cables go downwards)
  //   - plane 1 is at the right side (HDMI cables go to the right)
  // If one looks at a VMM3a hybrid (connector to detector readout is on the
  // bottom side), the channel 0 of the VMM 0 is always where the HDMI cable is
  // connected If the planes are correctly used as described above, the VMM IDs
  // are always in icreasing order PER HYBRID (e.g. 14, 15 or e.g. 0, 1)
  std::vector<std::tuple<uint8_t, uint8_t, uint16_t, uint8_t>> pVMMs{
      {1, 0, 1, 0}, {1, 0, 1, 1}, {1, 0, 1, 2}, {1, 0, 1, 3},
      {1, 1, 1, 6}, {1, 1, 1, 7}, {1, 1, 1, 8}, {1, 1, 1, 9}};

  // The tuples for the axes are defined as follows:
  // detector (choose a number between 0 and 255)
  // plane (0 or 1)
  // flip axis flag (0 or 1)
  // Using the convention described above, if the axis is NOT FLIPPED:
  //   - plane 0 is at the bottom and goes from left (0) to right (255)
  //   - plane 1 is at the right and goes from bottom (0) to top (255)
  // If the axis is FLIPPED:
  //   - plane 0 is at the bottom and goes from right (255) to left (0)
  //   - plane 1 is at the right and goes from top (0) to bottom (255)
  std::map<std::pair<uint8_t, uint8_t>, uint8_t> pAxes{{{1, 0}, 0},
                                                       {{1, 1}, 0}};

  std::vector<std::tuple<double, double, double>> pTranslation;

  std::vector<std::tuple<double, double, double>> pScale;

  std::vector<std::tuple<double, double, double>> pRotation;

  std::vector<std::vector<std::string>> pTransform;

  std::vector<std::tuple<double, double, double, double>> pTransformX;
  std::vector<std::tuple<double, double, double, double>> pTransformY;
  std::vector<std::tuple<double, double, double, double>> pTransformZ;

  std::string pFileName = "";

  double pTAC = 60.0;
  double pBC = 40.0;
  double pADCThreshold = 0;
  uint16_t pMinClusterSize = 1;
  uint16_t pCoincidentClusterSize = 1;
  // Maximum time difference between strips in time sorted cluster (x or y)
  double pDeltaTimeHits = 200.0;
  // Number of missing strips in strip sorted cluster (x or y)
  uint16_t pMissingStripsClusterX = 1;
  uint16_t pMissingStripsClusterY = 1;
  // Maximum time span for total cluster (x or y)
  double pSpanClusterTime = 500.0;
  // Maximum cluster time difference between matching clusters in two planes
  double pDeltaTimePlanes = 200.0;
  std::string pChannelMapping = "gem";
  std::string pGeometryFile = "";
  bool createJSON = false;
  bool useCalibration = false;
  bool calibrationHistogram = false;
  int pSaveWhat = 111;
  std::string pConditionCoincidence = "center-of-mass";
  double pChargeRatioLower = 0.5;
  double pChargeRatioUpper = 2;
  uint64_t nHits = 0;
  //**************************************************************
  // END PARAMETERS
  //**************************************************************

  std::string pRootFilename = "";
  std::string pCalFilename = "";
  std::string pInfo = "";

  double pBCTime_ns = 1000.0 / (double)pBC;
  double pOffsetPeriod = 1000.0 * 4096.9 / (double)pBC;

  std::map<std::tuple<uint8_t, uint8_t>, int> pChannels;
  std::map<uint8_t, int> pChannels0;
  std::map<uint8_t, int> pChannels1;
  std::map<std::pair<uint8_t, uint8_t>, std::pair<uint8_t, uint8_t>>
      pFecChip_DetectorPlane;
  std::multimap<std::pair<uint8_t, uint8_t>, uint8_t> pDetectorPlane_Fec;
  std::map<std::pair<uint8_t, uint8_t>, uint32_t> pOffsets;
  std::map<std::pair<uint8_t, uint8_t>, uint32_t> p_DetPlane_idx;
  std::map<uint8_t, uint8_t> pDets;
  std::vector<uint16_t> pFecs;
  std::map<std::pair<uint8_t, uint8_t>, std::string> pFecVMM_time0;

  std::vector<uint8_t> pSaveHits;
  std::vector<uint8_t> pSaveClustersPlane;
  std::vector<uint8_t> pSaveClustersDetector;

  bool fFound = false;
  bool vmmsFound = false;
  int pAlgo = 0;
  bool pIsPcap = false;
  bool pShowStats = true;
  bool pTimeZero = false;
  bool pIsPads[16];
  std::string pDataFormat = "ESS";
  int pPositions0[NUMFECS][16][64];
  int pPositions1[NUMFECS][16][64];
  int pDetectors[NUMFECS][16];
  int pPlanes[NUMFECS][16];
};