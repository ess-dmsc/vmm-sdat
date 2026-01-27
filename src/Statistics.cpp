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
**  Statistics.cpp
**
****************************************************************************/

#include <iomanip>
#include <iostream>
#include "log.h"
#include "Statistics.h"

Statistics::Statistics() {

}
void Statistics::CreatePCAPStats(Configuration &config) {
  m_counter_names.push_back("ParserFrameSeqErrors");
  m_counter_names.push_back("ParserFrameMissingErrors");
  m_counter_names.push_back("ParserFramecounterOverflows");
  m_counter_names.push_back("ParserTimestampSeqErrors");
  m_counter_names.push_back("ParserTimestampOverflows");
  m_counter_names.push_back("ParserBadFrames");
  m_counter_names.push_back("ParserGoodFrames");
  m_counter_names.push_back("ParserReadouts");
  m_counter_names.push_back("ParserMarkers");
  m_counter_names.push_back("ParserData");
  for (auto const &fec : config.pFecs) {
    m_counters.emplace(std::make_pair(std::make_pair(fec, "ParserFrameSeqErrors"), 0));
    m_counters.emplace(std::make_pair(std::make_pair(fec, "ParserFrameMissingErrors"), 0));
    m_counters.emplace(std::make_pair(std::make_pair(fec, "ParserFramecounterOverflows"), 0));
    m_counters.emplace(std::make_pair(std::make_pair(fec, "ParserTimestampSeqErrors"), 0));
    m_counters.emplace(std::make_pair(std::make_pair(fec, "ParserTimestampOverflows"), 0));
    m_counters.emplace(std::make_pair(std::make_pair(fec, "ParserBadFrames"), 0));
    m_counters.emplace(std::make_pair(std::make_pair(fec, "ParserGoodFrames"), 0));
    m_counters.emplace(std::make_pair(std::make_pair(fec, "ParserReadouts"), 0));
    m_counters.emplace(std::make_pair(std::make_pair(fec, "ParserMarkers"), 0));
    m_counters.emplace(std::make_pair(std::make_pair(fec, "ParserData"), 0));
  }
}

void Statistics::CreateFECStats(Configuration &config) {
  m_counter_names.push_back("TimestampTooLarge");
  m_counter_names.push_back("TimestampOverflow");

  for (auto const &fec : config.pFecs) {
    m_deltaTriggerTimestamp.emplace(std::make_pair(fec, 0));
    m_oldTriggerTimestamp.emplace(std::make_pair(fec, 0));
    m_lastFrameCounter.emplace(std::make_pair(fec, 0));
    m_counters.emplace(std::make_pair(std::make_pair(fec, "TimestampTooLarge"), 0));
    m_counters.emplace(std::make_pair(std::make_pair(fec, "TimestampOverflow"), 0));
  }
}


void Statistics::CreateClusterStats(Configuration &config) {
  int size = 0;
  for (auto const &det : config.pDets) {
    int the_det = config.pDets[det.first];
    auto plane0 = std::make_pair(det.first, 0);
    auto plane1 = std::make_pair(det.first, 1);
    auto plane2 = std::make_pair(det.first, 2);

    // initialize timestamps
    m_lowestCommonTriggerTimestamp_det[det.first] = 0;
    m_lowestCommonTriggerTimestamp_plane[plane0] = 0;
    m_lowestCommonTriggerTimestamp_plane[plane1] = 0;
    m_lowestCommonTriggerTimestamp_plane[plane2] = 0;
    //det.second == 0: first detector in config.pDets
    //Names ans limits only have to be created once
    if (det.second == 0) {
      m_units.emplace(std::make_pair("DeltaTimeHits", "ns"));
      m_factors.emplace(std::make_pair("DeltaTimeHits", 0.02));
      m_limits.emplace(std::make_pair("DeltaTimeHits", 1 + config.pDeltaTimeHits[the_det] * m_factors["DeltaTimeHits"]));
      m_stats_plane_names.push_back("DeltaTimeHits");

      m_units.emplace(std::make_pair("MissingStripsCluster", "strips"));
      m_factors.emplace(std::make_pair("MissingStripsCluster", 1));
      m_limits.emplace(std::make_pair("MissingStripsCluster", 1 + config.pMissingStripsClusterX[the_det] * m_factors["MissingStripsCluster"]));
      m_stats_plane_names.push_back("MissingStripsCluster");

      m_units.emplace(std::make_pair("SpanClusterTime", "ns"));
      m_factors.emplace(std::make_pair("SpanClusterTime", 0.02));
      m_limits.emplace(std::make_pair("SpanClusterTime", 1 + config.pSpanClusterTime[the_det] * m_factors["SpanClusterTime"]));
      m_stats_plane_names.push_back("SpanClusterTime");

      m_units.emplace(std::make_pair("ClusterSize", "strips"));
      m_factors.emplace(std::make_pair("ClusterSize", 1));
      m_limits.emplace(std::make_pair("ClusterSize", 1 + 64 * m_factors["ClusterSize"]));
      m_stats_plane_names.push_back("ClusterSize");

      m_units.emplace(std::make_pair("ClusterCntPlane", ""));
      m_factors.emplace(std::make_pair("ClusterCntPlane", 1));
      m_limits.emplace(std::make_pair("ClusterCntPlane", 1));
      m_stats_plane_names.push_back("ClusterCntPlane");

      m_units.emplace(std::make_pair("DeltaTimePlanes_0_1", "ns"));
      m_factors.emplace(std::make_pair("DeltaTimePlanes_0_1", 0.02));
      m_limits.emplace(std::make_pair("DeltaTimePlanes_0_1", 1 + config.pDeltaTimePlanes[the_det] * m_factors["DeltaTimePlanes_0_1"]));
      m_stats_detector_names.push_back("DeltaTimePlanes_0_1");

      m_units.emplace(std::make_pair("DeltaTimePlanes_1_2", "ns"));
      m_factors.emplace(std::make_pair("DeltaTimePlanes_1_2", 0.02));
      m_limits.emplace(std::make_pair("DeltaTimePlanes_1_2", 1 + config.pDeltaTimePlanes[the_det] * m_factors["DeltaTimePlanes_1_2"]));
      m_stats_detector_names.push_back("DeltaTimePlanes_1_2");

      m_units.emplace(std::make_pair("DeltaTimePlanes_0_2", "ns"));
      m_factors.emplace(std::make_pair("DeltaTimePlanes_0_2", 0.02));
      m_limits.emplace(std::make_pair("DeltaTimePlanes_0_2", 1 + config.pDeltaTimePlanes[the_det] * m_factors["DeltaTimePlanes_0_2"]));
      m_stats_detector_names.push_back("DeltaTimePlanes_0_2");

      m_units.emplace(std::make_pair("ChargeRatio_0_1", "%"));
      m_factors.emplace(std::make_pair("ChargeRatio_0_1", 0.1));
      m_limits.emplace(std::make_pair("ChargeRatio_0_1", 11));
      m_stats_detector_names.push_back("ChargeRatio_0_1");

      m_units.emplace(std::make_pair("ChargeRatio_1_0", "%"));
      m_factors.emplace(std::make_pair("ChargeRatio_1_0", 0.1));
      m_limits.emplace(std::make_pair("ChargeRatio_1_0", 10));
      m_stats_detector_names.push_back("ChargeRatio_1_0");

      m_units.emplace(std::make_pair("ChargeRatio_0_2", "%"));
      m_factors.emplace(std::make_pair("ChargeRatio_0_2", 0.1));
      m_limits.emplace(std::make_pair("ChargeRatio_0_2", 11));
      m_stats_detector_names.push_back("ChargeRatio_0_2");

      m_units.emplace(std::make_pair("ChargeRatio_2_0", "%"));
      m_factors.emplace(std::make_pair("ChargeRatio_2_0", 0.1));
      m_limits.emplace(std::make_pair("ChargeRatio_2_0", 10));
      m_stats_detector_names.push_back("ChargeRatio_2_0");

      m_units.emplace(std::make_pair("ChargeRatio_1_2", "%"));
      m_factors.emplace(std::make_pair("ChargeRatio_1_2", 0.1));
      m_limits.emplace(std::make_pair("ChargeRatio_1_2", 11));
      m_stats_detector_names.push_back("ChargeRatio_1_2");

      m_units.emplace(std::make_pair("ChargeRatio_2_1", "%"));
      m_factors.emplace(std::make_pair("ChargeRatio_2_1", 0.1));
      m_limits.emplace(std::make_pair("ChargeRatio_2_1", 10));
      m_stats_detector_names.push_back("ChargeRatio_2_1");

      m_units.emplace(std::make_pair("ClusterCntDetector", ""));
      m_factors.emplace(std::make_pair("ClusterCntDetector", 1));
      m_limits.emplace(std::make_pair("ClusterCntDetector", 1));
      m_stats_detector_names.push_back("ClusterCntDetector");

      m_units.emplace(std::make_pair("ClusterSizePads", "pads"));
      m_factors.emplace(std::make_pair("ClusterSizePads", 1));
      m_limits.emplace(std::make_pair("ClusterSizePads", 1 + 16 * m_factors["ClusterSizePads"]));
      m_stats_pads_names.push_back("ClusterSizePads");

      m_units.emplace(std::make_pair("ClusterExtensionPads_0", "pads"));
      m_factors.emplace(std::make_pair("ClusterExtensionPads_0", 1));
      m_limits.emplace(std::make_pair("ClusterExtensionPads_0", 1 + 8 * m_factors["ClusterExtensionPads_0"]));
      m_stats_pads_names.push_back("ClusterExtensionPads_0");

      m_units.emplace(std::make_pair("ClusterExtensionPads_1", "pads"));
      m_factors.emplace(std::make_pair("ClusterExtensionPads_1", 1));
      m_limits.emplace(std::make_pair("ClusterExtensionPads_1", 1 + 8 * m_factors["ClusterExtensionPads_1"]));
      m_stats_pads_names.push_back("ClusterExtensionPads_1");
     
      m_units.emplace(std::make_pair("ClusterCntPadDetector", ""));
      m_factors.emplace(std::make_pair("ClusterCntPadDetector", 1));
      m_limits.emplace(std::make_pair("ClusterCntPadDetector", 1));
      m_stats_pads_names.push_back("ClusterCntPadDetector");
    }
    size = static_cast<int>(m_limits["DeltaTimeHits"]);
    std::vector<long> v(size, 0);
    std::fill(v.begin(), v.end(), 0);
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane0, "DeltaTimeHits"), v));
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane1, "DeltaTimeHits"), v));
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane2, "DeltaTimeHits"), v));

    size = static_cast<int>(m_limits["MissingStripsCluster"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane0, "MissingStripsCluster"), v));
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane1, "MissingStripsCluster"), v));
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane2, "MissingStripsCluster"), v));

    size = static_cast<int>(m_limits["SpanClusterTime"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane0, "SpanClusterTime"), v));
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane1, "SpanClusterTime"), v));
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane2, "SpanClusterTime"), v));

    size = static_cast<int>(m_limits["ClusterSize"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane0, "ClusterSize"), v));
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane1, "ClusterSize"), v));
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane2, "ClusterSize"), v));

    size = static_cast<int>(m_limits["ClusterCntPlane"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane0, "ClusterCntPlane"), v));
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane1, "ClusterCntPlane"), v));
    m_stats_plane.emplace(std::make_pair(std::make_pair(plane2, "ClusterCntPlane"), v));

    size = static_cast<int>(m_limits["DeltaTimePlanes_0_1"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector.emplace(std::make_pair(std::make_pair(det.first, "DeltaTimePlanes_0_1"), v));

    size = static_cast<int>(m_limits["DeltaTimePlanes_1_2"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector.emplace(std::make_pair(std::make_pair(det.first, "DeltaTimePlanes_1_2"), v));

    size = static_cast<int>(m_limits["DeltaTimePlanes_0_2"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector.emplace(std::make_pair(std::make_pair(det.first, "DeltaTimePlanes_0_2"), v));

    size = static_cast<int>(m_limits["ChargeRatio_0_1"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector.emplace(std::make_pair(std::make_pair(det.first, "ChargeRatio_0_1"), v));

    size = static_cast<int>(m_limits["ChargeRatio_1_0"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector.emplace(std::make_pair(std::make_pair(det.first, "ChargeRatio_1_0"), v));

    size = static_cast<int>(m_limits["ChargeRatio_0_2"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector.emplace(std::make_pair(std::make_pair(det.first, "ChargeRatio_0_2"), v));

    size = static_cast<int>(m_limits["ChargeRatio_2_0"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector.emplace(std::make_pair(std::make_pair(det.first, "ChargeRatio_2_0"), v));

    size = static_cast<int>(m_limits["ChargeRatio_1_2"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector.emplace(std::make_pair(std::make_pair(det.first, "ChargeRatio_1_2"), v));

    size = static_cast<int>(m_limits["ChargeRatio_2_1"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector.emplace(std::make_pair(std::make_pair(det.first, "ChargeRatio_2_1"), v));

    size = static_cast<int>(m_limits["ClusterCntDetector"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector.emplace(std::make_pair(std::make_pair(det.first, "ClusterCntDetector"), v));

    size = static_cast<int>(m_limits["ClusterSizePads"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_pads.emplace(std::make_pair(std::make_pair(det.first, "ClusterSizePads"), v));

     size = static_cast<int>(m_limits["ClusterExtensionPads_0"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_pads.emplace(std::make_pair(std::make_pair(det.first, "ClusterExtensionPads_0"), v));

     size = static_cast<int>(m_limits["ClusterExtensionPads_1"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_pads.emplace(std::make_pair(std::make_pair(det.first, "ClusterExtensionPads_1"), v));

     size = static_cast<int>(m_limits["ClusterCntPadDetector"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_pads.emplace(std::make_pair(std::make_pair(det.first, "ClusterCntPadDetector"), v));

  }
}

long Statistics::GetStatsPads(std::string stats, uint8_t det, int n) {
  if (n < m_limits[stats.c_str()]) {
    return m_stats_pads[std::make_pair(det, stats.c_str())][n];
  }
  return -1;
}

void Statistics::SetStatsPads(std::string stats, uint8_t det, double value) {
  if (value * m_factors[stats] < m_limits[stats]) {
    m_stats_pads[std::make_pair(det, stats)][static_cast<unsigned int>(value * m_factors[stats])]++;
  } else {
    m_stats_pads[std::make_pair(det, stats)][m_limits[stats] - 1]++;
  }
}

long Statistics::GetStatsDetector(std::string stats, uint8_t det, int n) {
  if (n < m_limits[stats.c_str()]) {
    return m_stats_detector[std::make_pair(det, stats.c_str())][n];
  }
  return -1;
}

void Statistics::SetStatsDetector(std::string stats, uint8_t det, double value) {
  if (value * m_factors[stats] < m_limits[stats]) {
    m_stats_detector[std::make_pair(det, stats)][static_cast<unsigned int>(value * m_factors[stats])]++;
  } else {
    m_stats_detector[std::make_pair(det, stats)][m_limits[stats] - 1]++;
  }
}

long Statistics::GetStatsPlane(std::string stats, std::pair<uint8_t, uint8_t> dp, int n) {
  if (n < m_limits[stats]) {
    return m_stats_plane[std::make_pair(dp, stats)][n];
  }
  return -1;
}

void Statistics::SetStatsPlane(std::string stats, std::pair<uint8_t, uint8_t> dp, double value) {
  if (value * m_factors[stats] < m_limits[stats]) {
    m_stats_plane[std::make_pair(dp, stats)][static_cast<unsigned int>(value * m_factors[stats])]++;
  } else {
    m_stats_plane[std::make_pair(dp, stats)][m_limits[stats] - 1]++;
  }
}

void Statistics::IncrementCounter(std::string error, uint16_t fecId, uint64_t increment) { m_counters[std::make_pair(fecId, error)] += increment; }

long Statistics::GetCounter(std::string error, uint16_t fecId) { return m_counters[std::make_pair(fecId, error)]; }

double Statistics::GetOldTriggerTimestamp(uint16_t fecId) { return m_oldTriggerTimestamp[fecId]; }

void Statistics::SetOldTriggerTimestamp(uint16_t fecId, double srsTimestamp) { m_oldTriggerTimestamp[fecId] = srsTimestamp; }

double Statistics::GetFirstTriggerTimestamp(uint16_t fecId) { return m_firstTriggerTimestamp[fecId]; }

void Statistics::SetFirstTriggerTimestamp(uint16_t fecId, double srsTimestamp) { m_firstTriggerTimestamp[fecId] = srsTimestamp; }

double Statistics::GetMaxTriggerTimestamp(uint16_t fecId) { return m_maxTriggerTimestamp[fecId]; }

void Statistics::SetMaxTriggerTimestamp(uint16_t fecId, double srsTimestamp) { m_maxTriggerTimestamp[fecId] = srsTimestamp; }

uint64_t Statistics::GetLastFrameCounter(uint16_t fecId) { return m_lastFrameCounter[fecId]; }

void Statistics::SetLastFrameCounter(uint16_t fecId, uint64_t frameCounter) { m_lastFrameCounter[fecId] = frameCounter; }

double Statistics::GetLowestCommonTriggerTimestampDet(uint8_t det) { return m_lowestCommonTriggerTimestamp_det[det]; }

void Statistics::SetLowestCommonTriggerTimestampDet(uint8_t det, double val) { m_lowestCommonTriggerTimestamp_det[det] = val; }

double Statistics::GetLowestCommonTriggerTimestampPlane(std::pair<uint8_t, uint8_t> dp) { return m_lowestCommonTriggerTimestamp_plane[dp]; }

void Statistics::SetLowestCommonTriggerTimestampPlane(std::pair<uint8_t, uint8_t> dp, double val) { m_lowestCommonTriggerTimestamp_plane[dp] = val; }

void Statistics::PrintClusterStats(Configuration &config) {
  long totalPlane0 = 0;
  long totalPlane1 = 0;
  long totalDetector = 0;
  bool bothPlanes = false;
  for (auto const &det : config.pDets) {
    auto dp0 = std::make_pair(det.first, 0);
    auto dp1 = std::make_pair(det.first, 1);
    auto dp2 = std::make_pair(det.first, 2);
    corryvreckan::Log::setSection("Statistics");
    LOG(INFO) << "****************************************";
    LOG(INFO) << "Stats detector " << (int)det.first;
    LOG(INFO) << "****************************************";
    if(config.pIsPads[det.first]) {
      long cntPads = m_stats_pads[std::make_pair(det.first, "ClusterCntPadDetector")][0];
  
      for (auto const &stat : m_stats_pads_names) {
        LOG(INFO) << "****************************************";
        LOG(INFO) << stat;
        LOG(INFO) << "****************************************";
        std::vector<long> v = m_stats_pads[std::make_pair(det.first, stat)];
        for (unsigned int n = 0; n < static_cast<unsigned int>(m_limits[stat]); n++) {
          StatsOutput(n, v[n], stat, cntPads);
        }
        LOG(INFO) << "****************************************";
      }
  
    }
    else {
      long cnt = m_stats_detector[std::make_pair(det.first, "ClusterCntDetector")][0];
      long cnt0 = 1;
      if (m_stats_plane[std::make_pair(dp0, "ClusterCntPlane")][0] > 0) {
      cnt0 = m_stats_plane[std::make_pair(dp0, "ClusterCntPlane")][0];
      }
      long cnt1 = 1;
      if (m_stats_plane[std::make_pair(dp1, "ClusterCntPlane")][0] > 0) {
      cnt1 = m_stats_plane[std::make_pair(dp1, "ClusterCntPlane")][0];
      }
      long cnt2 = 1;
      if (m_stats_plane[std::make_pair(dp2, "ClusterCntPlane")][0] > 0) {
      cnt2 = m_stats_plane[std::make_pair(dp2, "ClusterCntPlane")][0];
      }
      for (auto const &stat : m_stats_plane_names) {
        if (config.GetDetectorPlane(dp0) == true) {
          LOG(INFO) << "****************************************";
          LOG(INFO) << "Plane 0: " << stat;
          LOG(INFO) << "****************************************";
          std::vector<long> v = m_stats_plane[std::make_pair(dp0, stat)];
          for (unsigned int n = 0; n < static_cast<unsigned int>(m_limits[stat]); n++) {
            StatsOutput(n, v[n], stat, cnt0);
          }
        }
        if (config.GetDetectorPlane(dp1) == true) {
          LOG(INFO) << "****************************************";
          LOG(INFO) << "****************************************";
          LOG(INFO) << "Plane 1: " << stat;
          LOG(INFO) << "****************************************";
          std::vector<long> v = m_stats_plane[std::make_pair(dp1, stat)];
          for (unsigned int n = 0; n < static_cast<unsigned int>(m_limits[stat]); n++) {
            StatsOutput(n, v[n], stat, cnt1);
          }
        }
        if (config.GetDetectorPlane(dp2) == true) {
          LOG(INFO) << "****************************************";
          LOG(INFO) << "****************************************";
          LOG(INFO) << "Plane 2: " << stat;
          LOG(INFO) << "****************************************";
          std::vector<long> v = m_stats_plane[std::make_pair(dp2, stat)];
          for (unsigned int n = 0; n < static_cast<unsigned int>(m_limits[stat]); n++) {
            StatsOutput(n, v[n], stat, cnt2);
          }
        }
        LOG(INFO) << "****************************************";
      }

      if ((config.GetDetectorPlane(dp0) == true && config.GetDetectorPlane(dp1)) == true) {
        for (auto const &stat : m_stats_detector_names) {
          if (config.GetDetectorPlane(dp2) == true || (stat.find("_0_2") == std::string::npos && stat.find("_1_2") == std::string::npos && stat.find("_2_0") == std::string::npos && stat.find("_2_1") == std::string::npos)) {
            LOG(INFO) << "****************************************";
            LOG(INFO) << stat;
            LOG(INFO) << "****************************************";
            std::vector<long> v = m_stats_detector[std::make_pair(det.first, stat)];
            for (unsigned int n = 0; n < static_cast<unsigned int>(m_limits[stat]); n++) {
              if (config.GetDetectorPlane(dp2) == true) {
                StatsOutput(n, v[n], stat, cnt, cnt0, cnt1, cnt2);
              } else {
                StatsOutput(n, v[n], stat, cnt, cnt0, cnt1);
              }
            }
            LOG(INFO) << "****************************************";
          }
        }
      }
    }
  }
}

void Statistics::PrintFECStats(Configuration &config) {
  corryvreckan::Log::setSection("Statistics");
  for (auto const &fec : config.pFecs) {
    LOG(INFO) << "****************************************";
    LOG(INFO) << "FEC " << (int)fec;
    LOG(INFO) << "****************************************";

    for (unsigned int n = 0; n < m_counter_names.size(); n++) {
      LOG(INFO) << m_counter_names[n] << ": " << GetCounter(m_counter_names[n], fec);
    }
    uint64_t first = GetFirstTriggerTimestamp(fec);
    uint64_t max = GetMaxTriggerTimestamp(fec);
    uint64_t last = GetOldTriggerTimestamp(fec);

    int overflow = GetCounter("TimestampOverflow", fec);
    m_acq_time = 0;

    if (overflow >= 1) {
      if (max <= 4294967295) {
        max = 4294967295;
      }
      if (max > 4294967295 && max <= 109951162777575) {
        max = 109951162777575;
      }
      m_acq_time = ((max - first) + (overflow - 1) * max + last) / 1000000.0;
    } else {
      m_acq_time = (max - first) / 1000000.0;
    }
    LOG(INFO) << "****************************************";
    LOG(INFO) << "Stats (acquisition):";
    LOG(INFO) << "****************************************";
    LOG(INFO) << "Acq time: " << std::setprecision(1) << std::fixed << m_acq_time << " ms";
    LOG(INFO) << "Hit rate FEC: " << std::scientific << 1000 * GetCounter("ParserData", fec) / m_acq_time << " hit/s";
    LOG(INFO) << "Data rate FEC: " << std::scientific << 1000 * GetCounter("ParserData", fec) * 48 / m_acq_time << " bit/s";
    LOG(INFO) << "****************************************";
  }
  LOG(INFO) << "****************************************";
  long cnt = 0;
  for (auto const &det : config.pDets) {
    cnt += GetStatsDetector("ClusterCntDetector", det.first, 0);
  }
  LOG(INFO) << "Total Cluster rate: " << std::scientific << 1000 * cnt / m_acq_time << " particles/s";
  LOG(INFO) << "****************************************";
}

void Statistics::StatsOutput(int n, long val, std::string stat, long cnt, long cnt0, long cnt1, long cnt2) {
  corryvreckan::Log::setSection("Statistics");
  if (m_limits[stat] > 1 && cnt > 0) {
    // if (cnt == 0)
    // cnt = 1;
    if (m_factors[stat] != 1) {
      LOG(INFO) << static_cast<unsigned int>(n / m_factors[stat]) << "-" << static_cast<unsigned int>(n / m_factors[stat] + 1 / m_factors[stat] - 1) << " " << m_units[stat] << ":  " << val << " (" << std::setprecision(1) << std::fixed
                << (100 * (double)val / (double)cnt) << " %)";
    } else {
      LOG(INFO) << static_cast<unsigned int>(n / m_factors[stat]) << " " << m_units[stat] << ":  " << val << " (" << std::setprecision(1) << std::fixed << (100 * (double)val / (double)cnt) << " %)";
    }
  } else {
    if (cnt0 > 0 && cnt1 > 0 && cnt2 > 0) {
      LOG(INFO) << val << " (common cluster in detector, " << std::setprecision(1) << std::fixed << (100 * (double)val / (double)cnt0) << " % plane 0, " << std::setprecision(1) << std::fixed << (100 * (double)val / (double)cnt1) << " % plane 1, "
                << std::setprecision(1) << std::fixed << (100 * (double)val / (double)cnt2) << " % plane 2)";
    } else if (cnt0 > 0 && cnt1 > 0) {
      LOG(INFO) << val << " (common cluster in detector, " << std::setprecision(1) << std::fixed << (100 * (double)val / (double)cnt0) << " % plane 0, " << std::setprecision(1) << std::fixed << (100 * (double)val / (double)cnt1) << " % plane 1)";
    } else {
      LOG(INFO) << val << " (" << std::setprecision(1) << std::fixed << (100 * (double)val / (double)cnt) << " %)";
    }
  }
}