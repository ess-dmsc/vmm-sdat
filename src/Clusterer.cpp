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
**  Clusterer.cpp
**
****************************************************************************/

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <future>
#include <iomanip>
#include <limits>
#include <thread>
#include <log.h>

#include <Math/MinimizerOptions.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "Clusterer.h"

#define UNUSED __attribute__((unused))

auto now = std::chrono::steady_clock::now;

auto timethis(std::function<void()> thunk) -> decltype((now() - now()).count()) {
  auto start = now();
  thunk();
  auto stop = now();
  return (stop - start).count();
}

Clusterer::Clusterer(Configuration& config, Statistics& stats) : m_config(config), m_stats(stats) {
  m_rootFile = RootFile::GetInstance(config);
}

Clusterer::~Clusterer() {
  RootFile::Dispose();
}

//====================================================================================================================
bool Clusterer::AnalyzeHits(double srsTimestamp, uint8_t fecId, uint8_t vmmId, uint16_t chNo, uint16_t bcid, uint16_t tdc, uint16_t adc, bool overThresholdFlag, double chipTime, uint8_t geoId, int event_counter) {
  int pos0 = m_config.pPositions0[fecId][vmmId][chNo];
  int pos1 = m_config.pPositions1[fecId][vmmId][chNo];

  if (pos0 == -1) {
    corryvreckan::Log::setSection("Clusterer");
    LOG(TRACE) << "Detector or Plane not defined for fecId " << (int)fecId << " and vmmId " << (int)vmmId;
    return true;
  }
  if ((tdc == 0 && adc == 0) || (overThresholdFlag == 1 && adc < 16)) {
    corryvreckan::Log::setSection("Clusterer");
    LOG(TRACE) << "Invalid data fecId " << (int)fecId << ", vmmId " << (int)vmmId << ", adc " << adc << ", tdc " << tdc << ", ovThr " << overThresholdFlag << ", time " << srsTimestamp;
    return true;
  }

  if (m_config.pDataFormat == "SRS") {
    // Biggest possible time should be:
    // from FEC: 2^42-1 = 0x3FFFFFFFFFF clock cycles
    // converted to ns: 0x3FFFFFFFFFF * bc period = 109951162777575 ns
    // added 31 offsets of 4096 * bc period ns: 3174400
    // total: 109951165951975 ns
    if (srsTimestamp > 109951165951975) {
      m_stats.IncrementCounter("TimestampTooLarge", fecId);
      corryvreckan::Log::setSection("Clusterer");
      LOG(WARNING) << "Timestamp " << static_cast<unsigned long long>(srsTimestamp) << " larger than 42 bit and 31 times trigger period for fecId " << (int)fecId << " and vmmId " << (int)vmmId;
    }

    if (srsTimestamp < m_stats.GetMaxTriggerTimestamp(fecId)) {
      // 42 bit: 0x1FFFFFFFFFF
      // 32 bit: 0xFFFFFFFF
      if (m_stats.GetMaxTriggerTimestamp(fecId) > 0x1FFFFFFFFFF + srsTimestamp) {
        m_stats.IncrementCounter("TimestampOverflow", fecId);
        corryvreckan::Log::setSection("Clusterer");
        LOG(WARNING) << "Overflow: fecId" << (int)fecId << ",  hitNr " << m_hitNr << ",  srsTimestamp " << static_cast<unsigned long long>(srsTimestamp) << ",  old srsTimestamp " << static_cast<unsigned long long>(m_stats.GetMaxTriggerTimestamp(fecId));

      }
    }
  }

  bool newData = false;
  double factor = 16.0;
  if (srsTimestamp >= m_stats.GetOldTriggerTimestamp(fecId) + factor * m_config.pOffsetPeriod) {
    newData = true;
  }

  if (newData) {
    if (m_config.pSaveWhat % 2 == 1) {
      m_rootFile->SaveHits();
    }

    if (m_config.pSaveWhat >= 10) {
      uint64_t ts = 0;
      for (auto const& fec : m_config.pFecs) {
        if (fec != 384) {
          if (ts == 0 || ts > m_stats.GetOldTriggerTimestamp(fec)) {
            ts = m_stats.GetOldTriggerTimestamp(fec);
          }
        }
      }

      for (auto const& det : m_config.pDets) {
        auto dp0 = std::make_pair(det.first, 0);
        auto dp1 = std::make_pair(det.first, 1);
        auto dp2 = std::make_pair(det.first, 2);
        if (m_stats.GetLowestCommonTriggerTimestampDet(det.first) < ts) {
          m_stats.SetLowestCommonTriggerTimestampPlane(dp0, ts);
          m_stats.SetLowestCommonTriggerTimestampPlane(dp1, ts);
          m_stats.SetLowestCommonTriggerTimestampPlane(dp2, ts);
          m_stats.SetLowestCommonTriggerTimestampDet(det.first, ts);
          AnalyzeClustersPlane(dp0);
          AnalyzeClustersPlane(dp1);
          AnalyzeClustersPlane(dp2);
          AnalyzeClustersDetector(det.first);
        }
      }
    }
    m_stats.SetOldTriggerTimestamp(fecId, srsTimestamp);
  }

  m_hitNr++;
  double totalTime = srsTimestamp + chipTime;
  auto det = m_config.pDetectors[fecId][vmmId];
  auto plane = m_config.pPlanes[fecId][vmmId];

  if (std::find(m_config.pSaveHits.begin(), m_config.pSaveHits.end(), det) != m_config.pSaveHits.end()) {
    if (m_config.pIsPads[det]) {
      Hit theHit0;
      theHit0.id = m_hitNr;
      theHit0.det = det;
      theHit0.plane = 0;
      theHit0.fec = fecId;
      theHit0.vmm = vmmId;
      theHit0.geo_id = geoId;
      theHit0.readout_time = srsTimestamp;
      theHit0.ch = chNo;
      theHit0.pos = (uint16_t)pos0;
      theHit0.bcid = bcid;
      theHit0.tdc = tdc;
      theHit0.adc = adc;
      theHit0.over_threshold = overThresholdFlag;
      theHit0.chip_time = chipTime;
      theHit0.time = totalTime;
      theHit0.event_counter = event_counter;
      m_rootFile->AddHits(std::move(theHit0));

      Hit theHit1;
      theHit1.id = m_hitNr;
      theHit1.det = det;
      theHit1.plane = 1;
      theHit1.fec = fecId;
      theHit1.vmm = vmmId;
      theHit1.geo_id = geoId;
      theHit1.readout_time = srsTimestamp;
      theHit1.ch = chNo;
      theHit1.pos = (uint16_t)pos1;
      theHit1.bcid = bcid;
      theHit1.tdc = tdc;
      theHit1.adc = adc;
      theHit1.over_threshold = overThresholdFlag;
      theHit1.chip_time = chipTime;
      theHit1.time = totalTime;
      theHit1.event_counter = event_counter;
      m_rootFile->AddHits(std::move(theHit1));
    } else {
      Hit theHit;
      theHit.id = m_hitNr;
      theHit.det = det;
      theHit.plane = plane;
      theHit.fec = fecId;
      theHit.vmm = vmmId;
      theHit.geo_id = geoId;
      theHit.readout_time = srsTimestamp;
      theHit.ch = chNo;
      theHit.pos = (uint16_t)pos0;
      theHit.bcid = bcid;
      theHit.tdc = tdc;
      theHit.adc = adc;
      theHit.over_threshold = overThresholdFlag;
      theHit.chip_time = chipTime;
      theHit.time = totalTime;
      theHit.event_counter = event_counter;
      m_rootFile->AddHits(std::move(theHit));
    }
  }

  if (overThresholdFlag) {
    // keep the overThresholdFlag as bit 15 of the ADC
    adc = adc + 32768;
  }
  int the_det = m_config.pDets[det];
  if (m_config.pADCThreshold[the_det] < 0) {
    if (overThresholdFlag) {
      if (m_config.pIsPads[det]) {
        m_hits_new[std::make_pair(det, 0)].emplace_back(totalTime, (uint16_t)pos0, adc, (uint16_t)pos1, event_counter);
        m_hits_new[std::make_pair(det, 1)].emplace_back(totalTime, (uint16_t)pos0, adc, (uint16_t)pos1, event_counter);
      } else {
        m_hits_new[std::make_pair(det, plane)].emplace_back(totalTime, (uint16_t)pos0, adc, (uint16_t)pos1, event_counter);
      }
    }
  } else {
    if ((adc >= m_config.pADCThreshold[the_det])) {
      if (m_config.pIsPads[det]) {
        m_hits_new[std::make_pair(det, 0)].emplace_back(totalTime, (uint16_t)pos0, adc, (uint16_t)pos1, event_counter);
        m_hits_new[std::make_pair(det, 1)].emplace_back(totalTime, (uint16_t)pos0, adc, (uint16_t)pos1, event_counter);
      } else {
        m_hits_new[std::make_pair(det, plane)].emplace_back(totalTime, (uint16_t)pos0, adc, (uint16_t)pos1, event_counter);
      }
    }
  }
  corryvreckan::Log::setSection("Clusterer");
  if (m_oldFecId != fecId || newData) {
    LOG(TRACE) << "fecId" << (int)fecId;
  }
  if (m_oldVmmId != vmmId || newData) {
    LOG(TRACE) << "detector " << (int)det << ", plane " << (int)plane << ", vmmId " << (int)vmmId;
  }
  LOG(TRACE) << "id0 " << pos0 << ", id1 " << pos1 << ", chNo " << chNo << ", overThreshold " << (int)overThresholdFlag;
  LOG(TRACE) << "bcid " << bcid << ", tdc " << tdc << ", adc " << (adc & 0x3FF);
  LOG(TRACE) << "total time " << totalTime << ", chip time " << chipTime;

  if (m_stats.GetFirstTriggerTimestamp(fecId) == 0) {
    m_stats.SetFirstTriggerTimestamp(fecId, srsTimestamp);
  }

  if (m_stats.GetMaxTriggerTimestamp(fecId) < srsTimestamp) {
    m_stats.SetMaxTriggerTimestamp(fecId, srsTimestamp);
  }
  m_oldVmmId = vmmId;
  m_oldFecId = fecId;

  return true;
}

//====================================================================================================================
int Clusterer::ClusterByTime(std::pair<uint8_t, uint8_t> dp) {

  ClusterContainer cluster;
  double maxDeltaTime = 0;
  int clusterCount = 0;
  double time1 = 0, time2 = 0;
  uint32_t adc1 = 0;
  uint16_t strip1 = 0;
  uint16_t strip2 = 0;
  uint16_t event_counter;

  for (auto& itHits : m_hits[dp]) {
    time2 = time1;

    time1 = (double)std::get<0>(itHits);
    strip1 = std::get<1>(itHits);
    adc1 = std::get<2>(itHits);
    strip2 = std::get<3>(itHits);
    event_counter = std::get<4>(itHits);
    if (!cluster.empty()) {
      if (std::fabs(time1 - time2) > m_config.pDeltaTimeHits[m_config.pDets[dp.first]]) {
        if (m_config.pIsPads[dp.first]) {
          clusterCount += ClusterByPad(dp, cluster, maxDeltaTime);
        } else {
          clusterCount += ClusterByStrip(dp, cluster, maxDeltaTime);
        }
        cluster.clear();
        maxDeltaTime = 0.0;
      } else {
        if (maxDeltaTime < std::fabs(time1 - time2)) {
          maxDeltaTime = (time1 - time2);
        }
      }
    }
    cluster.emplace_back(strip1, time1, adc1, strip2, event_counter);
  }

  if (!cluster.empty()) {
    if (m_config.pIsPads[dp.first]) {
      clusterCount += ClusterByPad(dp, cluster, maxDeltaTime);
    } else {
      clusterCount += ClusterByStrip(dp, cluster, maxDeltaTime);
    }
  }
  return clusterCount;
}

//====================================================================================================================
int Clusterer::ClusterByPad(std::pair<uint8_t, uint8_t> dp, ClusterContainer& cluster, double maxDeltaTime) {
  corryvreckan::Log::setSection("Clusterer");
  int stripCount = 0;
  uint16_t spanCluster = 0;
  int adcTotal = 0;
  int clusterCount = 0;
  double largestTime = 0;
  double smallestTime = 0;
  int delta_x = -1;
  int delta_y = -1;
  double time1 = 0, time2 = 0;
  double x_adc = 0, y_adc = 0, xy_time = 0;
  int adc1 = 0, adc2 = 0;
  int idX1 = 0, idX2 = 0, idY1 = 0, idY2 = 0;
  
  uint16_t event_counter1 = 0, event_counter2 = 0;

  std::vector<double> vADC, vIdx, vIdy, vTimes;

  const auto det = std::get<0>(dp);
  const auto plane = std::get<1>(dp);
  const bool useMaxADC = (m_config.pAlgo == 10);

  int np = 0;
  while (!cluster.empty()) {
    // --- choose seed ---------------------------------------------------------
    ClusterContainer::const_iterator seedIt = cluster.end();

    if (useMaxADC) {
      // pick pad with maximum ADC anywhere
      int maxAdc = -1;
      for (auto it = cluster.begin(); it != cluster.end(); ++it) {
        const int adc = (std::get<2>(*it) & 0x3FF);
        if (adc > maxAdc) {
          maxAdc = adc;
          seedIt = it;
        }
      }

    } else {
      seedIt = cluster.end();
      --seedIt;
    }

    // store properties of seed pad
    const int padX = std::get<0>(*seedIt);
    const double padTime = std::get<1>(*seedIt);
    const int padAdc = (std::get<2>(*seedIt) & 0x3FF);
    const int padY = std::get<3>(*seedIt);
    const uint16_t padEvtCtr = std::get<4>(*seedIt);

    np = 1;
    LOG(DEBUG) << "New cluster";
    LOG(DEBUG) << "\tpad " << np << ", (ADC,x,y)=(" << padAdc << ", " << padX << ", " << padY << ")";

    // start new cluster with seed
    ClusterContainer clusterFound;
    clusterFound.emplace_back(padX, padTime, padAdc, padY, padEvtCtr);

    vIdx.clear();
    vIdy.clear();
    vTimes.clear();
    vADC.clear();
    vIdx.emplace_back(padX);
    vIdy.emplace_back(padY);
    vTimes.emplace_back(padTime);
    vADC.emplace_back(padAdc);

    smallestTime = padTime;
    largestTime = padTime;
    adcTotal = padAdc;
    x_adc = padAdc * padX;
    y_adc = padAdc * padY;
    xy_time = padAdc * padTime;
    stripCount = 1;

    // remove seed from pool
    cluster.erase(seedIt);

    // --- grow cluster --------------------------------------------------------
    for (int n = 0; n < static_cast<int>(clusterFound.size()); ++n) {
      idX1 = std::get<0>(clusterFound[n]);
      time1 = std::get<1>(clusterFound[n]);
      adc1 = (std::get<2>(clusterFound[n]) & 0x3FF);
      idY1 = std::get<3>(clusterFound[n]);
      event_counter1 = std::get<4>(clusterFound[n]);

      auto itCluster = cluster.begin();
      while (itCluster != cluster.end()) {
        idX2 = std::get<0>(*itCluster);
        time2 = std::get<1>(*itCluster);
        adc2 = (std::get<2>(*itCluster) & 0x3FF);
        idY2 = std::get<3>(*itCluster);
        event_counter2 = std::get<4>(*itCluster);

        const bool timeOK = (event_counter1 == event_counter2) && (time2 - smallestTime <= m_config.pSpanClusterTime[m_config.pDets[dp.first]]) && (largestTime - time2 <= m_config.pSpanClusterTime[m_config.pDets[dp.first]]);

        bool match = false;
        if (useMaxADC) {
          match = (std::abs(padX - idX2) <= 1) && (std::abs(padY - idY2) <= 1) && timeOK;
        } else {
          match = (std::abs(idX1 - idX2) - 1 <= m_config.pMissingStripsClusterX[m_config.pDets[dp.first]]) && (std::abs(idY1 - idY2) - 1 <= m_config.pMissingStripsClusterY[m_config.pDets[dp.first]]) && timeOK;
        }

        if (match) {
          np++;
          LOG(DEBUG) << "\tpad " << np << ", (ADC,x,y)=(" << adc2 << ", " << idX2 << ", " << idY2 << ")";
          clusterFound.emplace_back(idX2, time2, adc2, idY2, event_counter2);

          ++stripCount;

          if (time1 < smallestTime)
            smallestTime = time1;
          if (time1 > largestTime)
            largestTime = time1;

          if (std::abs(idX1 - idX2) > delta_x)
            delta_x = std::abs(idX1 - idX2) - 1;
          if (std::abs(idY1 - idY2) > delta_y)
            delta_y = std::abs(idY1 - idY2) - 1;

          vIdx.emplace_back(idX2);
          vIdy.emplace_back(idY2);
          vTimes.emplace_back(time2);
          vADC.emplace_back(adc2);

          adcTotal += adc2;
          x_adc += adc2 * idX2;
          y_adc += adc2 * idY2;
          xy_time += adc2 * time2;

          itCluster = cluster.erase(itCluster);
        } else {
          ++itCluster;
        }
      }

      if (useMaxADC)
        break;
    }

    if (clusterFound.size() >= m_config.pMinClusterSize[m_config.pDets[dp.first]]) {
      clusterCount++;
      spanCluster = (largestTime - smallestTime);

      ClusterDetector clusterDetector;
      clusterDetector.id = m_cluster_detector_id;
      clusterDetector.det = det;
      clusterDetector.id0 = 0;
      clusterDetector.id1 = 0;
      clusterDetector.id2 = 0;
      clusterDetector.det = det;
      clusterDetector.event_counter = event_counter2;
      int max0 = *std::max_element( vIdx.begin(), vIdx.end() );
      int max1 = *std::max_element( vIdy.begin(), vIdy.end() );
      int min0 = *std::min_element( vIdx.begin(), vIdx.end() );
      int min1 = *std::min_element( vIdy.begin(), vIdy.end() );
            
      clusterDetector.size0 = max0 - min0+1;
      clusterDetector.size1 = max1 - min1+1;
      clusterDetector.size2 = stripCount;
      clusterDetector.adc0 = adcTotal;
      clusterDetector.adc1 = adcTotal;
      clusterDetector.adc2 = 0;
      clusterDetector.pos0 = x_adc / adcTotal;
      clusterDetector.pos1 = y_adc / adcTotal;

      clusterDetector.time0 = xy_time / adcTotal;
      clusterDetector.time1 = xy_time / adcTotal;
      
      clusterDetector.strips0 = std::move(vIdx);
      clusterDetector.strips1 = std::move(vIdy);
      clusterDetector.times0 = vTimes;
      clusterDetector.times1 = std::move(vTimes);
      clusterDetector.adcs0 = vADC;
      clusterDetector.adcs1 = std::move(vADC);
      clusterDetector.max_delta_time0 = maxDeltaTime;
      clusterDetector.max_delta_time1 = maxDeltaTime;
      clusterDetector.max_missing_strip0 = delta_x;
      clusterDetector.max_missing_strip1 = delta_y;
      clusterDetector.span_cluster0 = spanCluster;
      clusterDetector.span_cluster1 = spanCluster;
      
      m_clusters_detector[det].emplace_back(std::move(clusterDetector));
   
      LOG(TRACE) << "cluster: det " << (int)det << ", id: " << clusterDetector.id;
      LOG(TRACE) << "pos x/pos y: " << x_adc / adcTotal << "/" << y_adc / adcTotal;

      LOG(TRACE) << "time: " << clusterDetector.time0;
      LOG(TRACE) << "adc: " << clusterDetector.adc0;
      LOG(TRACE) << "extension plane0: " << clusterDetector.size0;
      LOG(TRACE) << "extension plane1: " << clusterDetector.size1;
      LOG(TRACE) << "number of pads: " << clusterDetector.size2;

      
      m_stats.SetStatsDetector("ClusterCntDetector", det, 0);
      m_cluster_detector_id++;
    }
  }
  if (m_config.pSaveWhat >= 100) {
    m_rootFile->SaveClustersDetector(std::move(m_clusters_detector[det]));
  }
  m_clusters_detector[det].clear();

  return clusterCount;
}

//====================================================================================================================
int Clusterer::ClusterByStrip(std::pair<uint8_t, uint8_t> dp, ClusterContainer& cluster, double maxDeltaTime) {
  corryvreckan::Log::setSection("Clusterer");
  int maxMissingStrip = 0;
  double spanCluster = 0;

  double startTime = 0;
  double largestTime = 0;
  double largestADCTime = 0;
  double largestADCPos = 0;
  double centerOfGravity = 0;
  double centerOfTime = 0;
  double centerOfGravity2 = 0;
  double centerOfTime2 = 0;
  double centerOfGravity_ovTh = 0;
  double centerOfTime_ovTh = 0;
  double centerOfGravity2_ovTh = 0;
  double centerOfTime2_ovTh = 0;
  long int totalADC = 0;
  long int totalADC2 = 0;
  long int totalADC_ovTh = 0;
  long int totalADC2_ovTh = 0;

  double time1 = 0;
  int idx_left = 0;
  int idx_right = 0;
  int adc1 = 0;
  int adc2 = 0;
  bool ovTh = false;
  int strip1 = 0;
  int strip2 = 0;
  int stripCount = 0;
  int clusterCount = 0;
  uint16_t event_counter1 = 0;
  std::vector<double> vADC;
  std::vector<double> vStrips;
  std::vector<double> vTimes;
  auto det = std::get<0>(dp);
  auto plane = std::get<1>(dp);

  std::sort(begin(cluster), end(cluster), [](const ClusterTuple& t1, const ClusterTuple& t2) { return std::get<0>(t1) < std::get<0>(t2) || (std::get<0>(t1) == std::get<0>(t2) && std::get<1>(t1) > std::get<1>(t2)); });
  for (auto& itCluster : cluster) {
    adc2 = adc1;
    strip2 = strip1;
    strip1 = std::get<0>(itCluster);
    time1 = std::get<1>(itCluster);
    adc1 = std::get<2>(itCluster);
    event_counter1 = std::get<4>(itCluster);
    if (adc1 > 1024) {
      ovTh = true;
    } else {
      ovTh = false;
    }
    adc1 = adc1 & 0x3FF;

    // At beginning of cluster, set start time of cluster
    if (stripCount == 0) {
      maxMissingStrip = 0;
      idx_left = 0;
      idx_right = 0;
      startTime = time1;
      largestTime = time1;
      largestADCTime = time1;
      largestADCPos = strip1;
      LOG(TRACE) << "detector: " << (int)det << " plane: " << (int)plane;
    }

    // Add members of a cluster, if it is either the beginning of a cluster,
    // or if strip gap and time span is correct
    if (stripCount == 0 || (((plane == 2 && m_config.pAlgo == 5) || (std::fabs(strip1 - strip2) - 1 <= m_config.pMissingStripsClusterX[m_config.pDets[dp.first]])) && time1 - startTime <= m_config.pSpanClusterTime[m_config.pDets[dp.first]] && largestTime - time1 <= m_config.pSpanClusterTime[m_config.pDets[dp.first]])) {
      LOG(TRACE) << "strip: " << strip1 << " time: " << static_cast<unsigned long long>(time1) << " adc: " << adc1;

      if (adc1 > adc2) {
        largestADCTime = time1;
        largestADCPos = strip1;
      }

      if (time1 == largestTime) {
        idx_right = stripCount;
      }
      if (time1 > largestTime) {
        idx_left = stripCount;
        idx_right = stripCount;
        largestTime = time1;
      }
      if (time1 < startTime) {
        startTime = time1;
      }
      if (stripCount > 0 && maxMissingStrip < std::fabs(strip1 - strip2) - 1) {
        maxMissingStrip = std::fabs(strip1 - strip2) - 1;
      }
      spanCluster = (largestTime - startTime);
      totalADC += adc1;
      totalADC2 += adc1 * adc1;
      centerOfGravity += strip1 * adc1;
      centerOfTime += time1 * adc1;
      centerOfGravity2 += strip1 * adc1 * adc1;
      centerOfTime2 += time1 * adc1 * adc1;

      if (ovTh) {
        totalADC_ovTh += adc1;
        totalADC2_ovTh += adc1 * adc1;
        centerOfGravity_ovTh += strip1 * adc1;
        centerOfTime_ovTh += time1 * adc1;
        centerOfGravity2_ovTh += strip1 * adc1 * adc1;
        centerOfTime2_ovTh += time1 * adc1 * adc1;
      }

      vStrips.emplace_back(strip1);
      vTimes.emplace_back(time1);
      vADC.emplace_back(adc1);
      stripCount++;
    }
    // Stop clustering if gap between strips is too large or time span too
    // long
    else if ((!(plane == 2 && m_config.pAlgo == 5) && (std::fabs(strip1 - strip2) - 1 > m_config.pMissingStripsClusterX[m_config.pDets[dp.first]])) || time1 - startTime > m_config.pSpanClusterTime[m_config.pDets[dp.first]] || largestTime - time1 > m_config.pSpanClusterTime[m_config.pDets[dp.first]]) {
      // Valid cluster
      if (stripCount < m_config.pMinClusterSize[m_config.pDets[dp.first]] || totalADC == 0) {
        LOG(TRACE) << "Invalid clusters size: " << stripCount;

      } else {
        spanCluster = (largestTime - startTime);
        centerOfGravity = (centerOfGravity / totalADC);
        centerOfTime = (centerOfTime / totalADC);
        centerOfGravity2 = (centerOfGravity2 / totalADC2);
        centerOfTime2 = (centerOfTime2 / totalADC2);

        if (totalADC_ovTh > 0) {
          centerOfGravity_ovTh = (centerOfGravity_ovTh / totalADC_ovTh);
          centerOfTime_ovTh = (centerOfTime_ovTh / totalADC_ovTh);
        }
        if (totalADC2_ovTh > 0) {
          centerOfGravity2_ovTh = (centerOfGravity2_ovTh / totalADC2_ovTh);
          centerOfTime2_ovTh = (centerOfTime2_ovTh / totalADC2_ovTh);
        }

        m_stats.SetStatsPlane("DeltaTimeHits", dp, maxDeltaTime);
        m_stats.SetStatsPlane("MissingStripsCluster", dp, maxMissingStrip);
        m_stats.SetStatsPlane("SpanClusterTime", dp, spanCluster);
        m_stats.SetStatsPlane("ClusterSize", dp, stripCount);

        ClusterPlane clusterPlane;
        clusterPlane.event_counter = event_counter1;
        clusterPlane.size = stripCount;
        clusterPlane.adc = totalADC;
        clusterPlane.time = centerOfTime;
        clusterPlane.pos = centerOfGravity;

        clusterPlane.time_charge2 = centerOfTime2;
        clusterPlane.pos_charge2 = centerOfGravity2;

        double time_utpc = 0;
        double pos_utpc = 0;
        double time_algo = 0;
        double pos_algo = 0;
        AlgorithmUTPC(idx_left, idx_right, vADC, vStrips, vTimes, pos_utpc, time_utpc, pos_algo, time_algo);

        if (m_config.pAlgo == 0) {
          pos_algo = 0;
          time_algo = 0;
        } else if (m_config.pAlgo == 1) {
          pos_algo = FitGaussianWeights(vStrips, vADC, centerOfGravity);
          time_algo = centerOfTime;
        }
        // COG over Threshold only
        else if (m_config.pAlgo == 2) {
          pos_algo = centerOfGravity_ovTh;
          time_algo = centerOfTime_ovTh;
        }
        // COG2 over Threshold only
        else if (m_config.pAlgo == 3) {
          pos_algo = centerOfGravity2_ovTh;
          time_algo = centerOfTime2_ovTh;
        }
        // time of highest ADC
        else if (m_config.pAlgo == 4) {
          pos_algo = largestADCPos;
          time_algo = largestADCTime;
        }
        // trigger pattern
        else if (m_config.pAlgo == 5) {
          pos_algo = 0;
          time_algo = 0;
          if (plane == 2) {
            for (int n = 0; n < vStrips.size(); n++) {
              time_algo = time_algo + pow(2.0, vStrips[n]);
            }
          }
        }

        clusterPlane.time_utpc = time_utpc;
        clusterPlane.pos_utpc = pos_utpc;
        clusterPlane.time_algo = time_algo;
        clusterPlane.pos_algo = pos_algo;

        clusterPlane.plane_coincidence = false;
        clusterPlane.max_delta_time = maxDeltaTime;
        clusterPlane.max_missing_strip = maxMissingStrip;
        clusterPlane.span_cluster = spanCluster;
        clusterPlane.strips = std::move(vStrips);
        clusterPlane.times = std::move(vTimes);
        clusterPlane.adcs = std::move(vADC);

        m_cluster_id++;
        LOG(TRACE) << "Cluster id: " << m_cluster_id;
        clusterPlane.id = static_cast<uint32_t>(m_cluster_id);
        clusterPlane.det = det;
        clusterPlane.plane = plane;
        /*
        if (plane == 0) {
          auto det_plane = std::make_pair(det, 1);
          if (m_config.GetDetectorPlane(det_plane) == false) {
            ClusterPlane clusterOtherPlane;
            clusterOtherPlane = clusterPlane;
            clusterOtherPlane.plane = 1;
            m_clusters_new[det_plane].emplace_back(std::move(clusterOtherPlane));
          }
        } else {
          auto det_plane = std::make_pair(det, 0);
          if (m_config.GetDetectorPlane(det_plane) == false) {
            ClusterPlane clusterOtherPlane;
            clusterOtherPlane = clusterPlane;
            clusterOtherPlane.plane = 0;
            m_clusters_new[det_plane].emplace_back(std::move(clusterOtherPlane));
          }
        }
        */
        m_clusters_new[dp].emplace_back(std::move(clusterPlane));
        m_stats.SetStatsPlane("ClusterCntPlane", dp, 0);

        clusterCount++;
      }

      // Clear vectors
      vADC.clear();
      vStrips.clear();
      vTimes.clear();
      // Strip that caused gap in cluster is added as first strip of new
      // cluster
      vStrips.emplace_back(strip1);
      vTimes.emplace_back(time1);
      vADC.emplace_back(adc1);
      stripCount = 1;
      largestADCTime = time1;
      largestADCPos = strip1;
      idx_right = 0;
      idx_left = 0;
      largestTime = time1;
      startTime = time1;
      maxMissingStrip = 0;
      spanCluster = 0;
      totalADC = adc1;
      totalADC2 = adc1 * adc1;
      centerOfGravity = strip1 * adc1;
      centerOfTime = time1 * adc1;
      centerOfGravity2 = strip1 * adc1 * adc1;
      centerOfTime2 = time1 * adc1 * adc1;

      if (ovTh) {
        totalADC_ovTh = adc1;
        totalADC2_ovTh = adc1 * adc1;
        centerOfGravity_ovTh = strip1 * adc1;
        centerOfTime_ovTh = time1 * adc1;
        centerOfGravity2_ovTh = strip1 * adc1 * adc1;
        centerOfTime2_ovTh = time1 * adc1 * adc1;
      }
    }
  }

  // At the end of the clustering, check again if there is a last valid
  // cluster
  if (stripCount >= m_config.pMinClusterSize[m_config.pDets[dp.first]] && totalADC > 0) {
    spanCluster = (largestTime - startTime);
    centerOfGravity = (centerOfGravity / totalADC);
    centerOfTime = (centerOfTime / totalADC);
    centerOfGravity2 = (centerOfGravity2 / totalADC2);
    centerOfTime2 = (centerOfTime2 / totalADC2);

    m_stats.SetStatsPlane("DeltaTimeHits", dp, maxDeltaTime);
    m_stats.SetStatsPlane("MissingStripsCluster", dp, maxMissingStrip);
    m_stats.SetStatsPlane("SpanClusterTime", dp, spanCluster);
    m_stats.SetStatsPlane("ClusterSize", dp, stripCount);

    ClusterPlane clusterPlane;
    clusterPlane.event_counter = event_counter1;
    clusterPlane.size = stripCount;
    clusterPlane.adc = totalADC;
    clusterPlane.time = centerOfTime;
    clusterPlane.pos = centerOfGravity;
    clusterPlane.time_charge2 = centerOfTime2;
    clusterPlane.pos_charge2 = centerOfGravity2;

    double time_utpc = 0;
    double pos_utpc = 0;
    double time_algo = 0;
    double pos_algo = 0;
    // Functions calculates uTPC position and time independent of the algo
    // for algo 6 and 7 it calculates uTPC with COG and COG2
    AlgorithmUTPC(idx_left, idx_right, vADC, vStrips, vTimes, pos_utpc, time_utpc, pos_algo, time_algo);

    if (m_config.pAlgo == 0) {
      pos_algo = 0;
      time_algo = 0;
    } else if (m_config.pAlgo == 1) {
      pos_algo = FitGaussianWeights(vStrips, vADC, centerOfGravity);
      time_algo = centerOfTime;
    }
    // COG over Threshold only
    else if (m_config.pAlgo == 2) {
      pos_algo = centerOfGravity_ovTh;
      time_algo = centerOfTime_ovTh;
    }
    // COG2 over Threshold only
    else if (m_config.pAlgo == 3) {
      pos_algo = centerOfGravity2_ovTh;
      time_algo = centerOfTime2_ovTh;
    }
    // time of highest ADC
    else if (m_config.pAlgo == 4) {
      pos_algo = largestADCPos;
      time_algo = largestADCTime;
    }
    // trigger pattern
    else if (m_config.pAlgo == 5) {
      pos_algo = 0;
      time_algo = 0;
      if (plane == 2) {
        for (int n = 0; n < vStrips.size(); n++) {
          time_algo = time_algo + pow(2.0, vStrips[n]);
        }
      }
    }

    clusterPlane.time_utpc = time_utpc;
    clusterPlane.pos_utpc = pos_utpc;
    clusterPlane.time_algo = time_algo;
    clusterPlane.pos_algo = pos_algo;

    clusterPlane.plane_coincidence = false;
    clusterPlane.max_delta_time = maxDeltaTime;
    clusterPlane.max_missing_strip = maxMissingStrip;
    clusterPlane.span_cluster = spanCluster;
    clusterPlane.strips = std::move(vStrips);
    clusterPlane.times = std::move(vTimes);
    clusterPlane.adcs = std::move(vADC);
    m_cluster_id++;
    LOG(TRACE) << "Cluster id: " << m_cluster_id;

    clusterPlane.id = static_cast<uint32_t>(m_cluster_id);
    clusterPlane.det = det;
    clusterPlane.plane = plane;
    /*
    if (plane == 0) {
      auto det_plane = std::make_pair(det, 1);
      if (m_config.GetDetectorPlane(det_plane) == false) {
        ClusterPlane clusterOtherPlane;
        clusterOtherPlane = clusterPlane;
        clusterOtherPlane.plane = 1;
        m_clusters_new[det_plane].emplace_back(std::move(clusterOtherPlane));
      }
    } else {
      auto det_plane = std::make_pair(det, 0);
      if (m_config.GetDetectorPlane(det_plane) == false) {
        ClusterPlane clusterOtherPlane;
        clusterOtherPlane = clusterPlane;
        clusterOtherPlane.plane = 0;
        m_clusters_new[det_plane].emplace_back(std::move(clusterOtherPlane));
      }
    }
    */
    m_clusters_new[dp].emplace_back(std::move(clusterPlane));
    m_stats.SetStatsPlane("ClusterCntPlane", dp, 0);

    clusterCount++;
  }
  return clusterCount;
}

void Clusterer::AlgorithmUTPC(int idx_min_largest_time, int idx_max_largest_time, std::vector<double>& vADC, std::vector<double>& vStrips, std::vector<double>& vTimes, double& positionUTPC, double& timeUTPC, double& positionAlgo, double& timeAlgo) {
  double a1 = 0, a2 = 0, a3 = 0, p1 = 0, p2 = 0, p3 = 0, t1 = 0, t2 = 0, t3 = 0;
  int idx_largest_time = 0;
  // One largest time exists
  if (idx_max_largest_time == idx_min_largest_time) {
    idx_largest_time = idx_max_largest_time;
    positionUTPC = vStrips[idx_largest_time];
    timeUTPC = vTimes[idx_largest_time];
  } else {
    // More than one largest time, the right most largest time strip
    // is closer to the end of the track than the lest most
    if (vStrips.size() - 1 - idx_max_largest_time < idx_min_largest_time) {
      idx_largest_time = idx_max_largest_time;
    }
    // More than one largest time, the left most largest time strip
    // is closer to the end of the track than the right most
    else if (vStrips.size() - 1 - idx_max_largest_time > idx_min_largest_time) {
      idx_largest_time = idx_min_largest_time;
    }
    // More than one largest time, the left and right most largest time strips
    // have an identical distance to the start/end of the track
    // Take the strip with the larges ADC
    else {
      if (vADC[idx_min_largest_time] > vADC[idx_max_largest_time]) {
        idx_largest_time = idx_min_largest_time;
      } else if (vADC[idx_min_largest_time] <= vADC[idx_max_largest_time]) {
        idx_largest_time = idx_max_largest_time;
      }
    }
  }

  positionUTPC = vStrips[idx_largest_time];
  timeUTPC = vTimes[idx_largest_time];

  p2 = vStrips[idx_largest_time];
  a2 = vADC[idx_largest_time];
  t2 = vTimes[idx_largest_time];
  if (idx_largest_time > 0) {
    p1 = vStrips[idx_largest_time - 1];
    a1 = vADC[idx_largest_time - 1];
    t1 = vTimes[idx_largest_time - 1];
  }
  if (idx_largest_time < vStrips.size() - 1) {
    p3 = vStrips[idx_largest_time + 1];
    a3 = vADC[idx_largest_time + 1];
    t3 = vTimes[idx_largest_time + 1];
  }
  if (m_config.pAlgo == 6) {
    positionAlgo = (p1 * a1 * a1 + p2 * a2 * a2 + p3 * a3 * a3) / (a1 * a1 + a2 * a2 + a3 * a3);
    timeAlgo = (t1 * a1 * a1 + t2 * a2 * a2 + t3 * a3 * a3) / (a1 * a1 + a2 * a2 + a3 * a3);
  } else if (m_config.pAlgo == 7) {
    positionAlgo = (p1 * a1 + p2 * a2 + p3 * a3) / (a1 + a2 + a3);
    timeAlgo = (t1 * a1 + t2 * a2 + t3 * a3) / (a1 + a2 + a3);
  }
}

double Clusterer::FitGaussianWeights(const std::vector<double>& x, const std::vector<double>& y, double cog) {
  const size_t n = x.size();
  if (n <= 2 || y.size() != n)
    return cog;

  auto [xmin_it, xmax_it] = std::minmax_element(x.begin(), x.end());
  const double xmin = *xmin_it, xmax = *xmax_it;
  const double spanX = std::max(1e-9, xmax - xmin);

  auto itmaxY = std::max_element(y.begin(), y.end());
  const int iMax = int(std::distance(y.begin(), itmaxY));
  const double A0 = std::min(1023.0, std::max(1e-12, *itmaxY));
  const double mu0 = cog;

  const double sigma0 = std::clamp(spanX / 6.0, 0.3, 10.0);
  const double sigma_lo = std::max(0.2, spanX / 20.0);
  const double sigma_hi = std::max(1.0, 1.0 * spanX);

  TGraphErrors gr;
  gr.Set((int)n);

  int k = 0;
  double ymax = 0.0;
  for (size_t i = 0; i < n; ++i) {
    gr.SetPoint(k, x[i], y[i]);
    double ey = std::sqrt(std::max(1.0, y[i]));
    gr.SetPointError(k, 0.0, ey);
    ++k;
    if (y[i] > ymax)
      ymax = y[i];
  }

  // ----- model & bounds -----
  TF1 f("f_mu_err", "gaus(0)", xmin, xmax);
  f.SetParNames("A", "mu", "sigma");
  f.SetParameters(A0, mu0, sigma0);

  f.SetParLimits(0, 0.0, 5.0 * 1023.0);
  f.SetParLimits(1, xmin, xmax);
  f.SetParLimits(2, 0.3, 10.0);

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-6);

  TFitResultPtr r = gr.Fit(&f, "RSWQ");
  if (r.Get() == nullptr)
    return cog;
  const bool ok = (r->Status() == 0) && (f.GetNDF() >= 0 && (r->CovMatrixStatus() >= 2));
  if (!ok) {
    return cog;
  }
  return f.GetParameter(1);
}

//====================================================================================================================
int Clusterer::MatchClustersDetector(uint8_t det) {
  corryvreckan::Log::setSection("Clusterer");
  int clusterCount = 0;
  auto dp0 = std::make_pair(det, 0);
  auto dp1 = std::make_pair(det, 1);
  auto dp2 = std::make_pair(det, 2);

  bool found1 = false;
  bool found2 = true;
  ClusterVectorPlane::iterator itStartPlane1 = begin(m_clusters[dp1]);
  ClusterVectorPlane::iterator itStartPlane2 = begin(m_clusters[dp2]);

  for (auto& c0 : m_clusters[dp0]) {
    if (m_config.GetDetectorPlane(dp2) == true) {
      found2 = false;
    }
    double minDelta1 = std::numeric_limits<double>::max();
    double lastDelta_t1 = std::numeric_limits<double>::max();
    double delta_t1 = std::numeric_limits<double>::max();
    bool isFirstMatch1 = true;
    ClusterVectorPlane::iterator bestMatchPlane1 = end(m_clusters[dp1]);
    ClusterVectorPlane::iterator bestMatchPlane2 = end(m_clusters[dp2]);
    for (ClusterVectorPlane::iterator c1 = itStartPlane1; c1 != end(m_clusters[dp1]); ++c1) {
      if ((*c1).plane_coincidence == false) {
        double chargeRatio1 = (double)(c0).adc / (double)(*c1).adc;
        lastDelta_t1 = delta_t1;
        delta_t1 = (*c1).time - c0.time;
        if (m_config.pConditionCoincidence == "utpc") {
          delta_t1 = (*c1).time_utpc - c0.time_utpc;
        } else if (m_config.pConditionCoincidence == "charge2") {
          delta_t1 = (*c1).time_charge2 - c0.time_charge2;
        }
        if (std::fabs(delta_t1) <= m_config.pDeltaTimePlanes[m_config.pDets[det]]) {
          if (isFirstMatch1) {
            itStartPlane1 = c1;
            isFirstMatch1 = false;
          }
        }

        if (chargeRatio1 >= m_config.pChargeRatioLower[m_config.pDets[det]] && chargeRatio1 <= m_config.pChargeRatioUpper[m_config.pDets[det]] && std::fabs(delta_t1) < minDelta1 && std::fabs(delta_t1) <= m_config.pDeltaTimePlanes[m_config.pDets[det]] && (c0.size + (*c1).size >= m_config.pCoincidentClusterSize[m_config.pDets[det]])) {
          minDelta1 = std::fabs(delta_t1);
          bestMatchPlane1 = c1;
        }
        if (std::fabs(delta_t1) > std::fabs(lastDelta_t1) && std::fabs(delta_t1) > m_config.pDeltaTimePlanes[m_config.pDets[det]]) {
          break;
        }
      }
    }

    if (bestMatchPlane1 != end(m_clusters[dp1])) {
      found1 = true;
      if (m_config.GetDetectorPlane(dp2) == true) {
        double minDelta2 = std::numeric_limits<double>::max();
        double lastDelta_t2 = std::numeric_limits<double>::max();
        double delta_t2 = std::numeric_limits<double>::max();
        bool isFirstMatch2 = true;

        for (ClusterVectorPlane::iterator c2 = itStartPlane2; c2 != end(m_clusters[dp2]); ++c2) {
          if ((*c2).plane_coincidence == false) {
            double chargeRatio2 = (double)(c0).adc / (double)(*c2).adc;
            lastDelta_t2 = delta_t2;
            delta_t2 = (*c2).time - c0.time;
            if (m_config.pConditionCoincidence == "utpc") {
              delta_t2 = (*c2).time_utpc - c0.time_utpc;
            } else if (m_config.pConditionCoincidence == "charge2") {
              delta_t2 = (*c2).time_charge2 - c0.time_charge2;
            }
            if (std::fabs(delta_t2) <= m_config.pDeltaTimePlanes[m_config.pDets[det]]) {
              if (isFirstMatch2) {
                itStartPlane2 = c2;
                isFirstMatch2 = false;
              }
            }
            if (chargeRatio2 >= m_config.pChargeRatioLower[m_config.pDets[det]] && chargeRatio2 <= m_config.pChargeRatioUpper[m_config.pDets[det]] && std::fabs(delta_t2) < minDelta2 && std::fabs(delta_t2) <= m_config.pDeltaTimePlanes[m_config.pDets[det]] && (c0.size + (*bestMatchPlane1).size + (*c2).size >= m_config.pCoincidentClusterSize[m_config.pDets[det]])) {
              minDelta2 = std::fabs(delta_t2);
              bestMatchPlane2 = c2;
            }
            if (std::fabs(delta_t2) > std::fabs(lastDelta_t2) && std::fabs(delta_t2) > m_config.pDeltaTimePlanes[m_config.pDets[det]]) {
              break;
            }
          }
        }
        if (bestMatchPlane2 != end(m_clusters[dp2])) {
          found2 = true;
        }
      }
      if (found1 == true && found2 == true) {
        m_cluster_detector_id++;
        ClusterDetector clusterDetector;
        c0.plane_coincidence = true;
        (*bestMatchPlane1).plane_coincidence = true;
        clusterDetector.id = m_cluster_detector_id;
        clusterDetector.det = det;
        clusterDetector.event_counter = (*bestMatchPlane1).event_counter;
        clusterDetector.id0 = c0.id;
        clusterDetector.id1 = (*bestMatchPlane1).id;
        clusterDetector.id2 = 0;
        clusterDetector.size0 = c0.size;
        clusterDetector.size1 = (*bestMatchPlane1).size;
        clusterDetector.size2 = 0;
        clusterDetector.adc0 = c0.adc;
        clusterDetector.adc1 = (*bestMatchPlane1).adc;
        clusterDetector.adc2 = 0;
        clusterDetector.pos0 = c0.pos;
        clusterDetector.pos1 = (*bestMatchPlane1).pos;
        clusterDetector.pos2 = 0;
        clusterDetector.pos0_utpc = c0.pos_utpc;
        clusterDetector.pos1_utpc = (*bestMatchPlane1).pos_utpc;
        clusterDetector.pos2_utpc = 0;
        clusterDetector.pos0_charge2 = c0.pos_charge2;
        clusterDetector.pos1_charge2 = (*bestMatchPlane1).pos_charge2;
        clusterDetector.pos2_charge2 = 0;
        clusterDetector.pos0_algo = c0.pos_algo;
        clusterDetector.pos1_algo = (*bestMatchPlane1).pos_algo;
        clusterDetector.pos2_algo = 0;
        clusterDetector.time0 = c0.time;
        clusterDetector.time1 = (*bestMatchPlane1).time;
        clusterDetector.time2 = 0;
        clusterDetector.time0_utpc = c0.time_utpc;
        clusterDetector.time1_utpc = (*bestMatchPlane1).time_utpc;
        clusterDetector.time2_utpc = 0;
        clusterDetector.time0_charge2 = c0.time_charge2;
        clusterDetector.time1_charge2 = (*bestMatchPlane1).time_charge2;
        clusterDetector.time2_charge2 = 0;
        clusterDetector.time0_algo = c0.time_algo;
        clusterDetector.time1_algo = (*bestMatchPlane1).time_algo;
        clusterDetector.time2_algo = 0;
        clusterDetector.dt0 = clusterDetector.time0 - last_time0;
        clusterDetector.dt1 = clusterDetector.time1 - last_time1;
        clusterDetector.dt2 = 0;
        clusterDetector.max_delta_time0 = c0.max_delta_time;
        clusterDetector.max_delta_time1 = (*bestMatchPlane1).max_delta_time;
        clusterDetector.max_delta_time2 = 0;
        clusterDetector.max_missing_strip0 = c0.max_missing_strip;
        clusterDetector.max_missing_strip1 = (*bestMatchPlane1).max_missing_strip;
        clusterDetector.max_missing_strip2 = 0;
        clusterDetector.span_cluster0 = c0.span_cluster;
        clusterDetector.span_cluster1 = (*bestMatchPlane1).span_cluster;
        clusterDetector.span_cluster2 = 0;
        clusterDetector.strips0 = c0.strips;
        clusterDetector.times0 = c0.times;
        clusterDetector.adcs0 = c0.adcs;
        clusterDetector.strips1 = (*bestMatchPlane1).strips;
        clusterDetector.times1 = (*bestMatchPlane1).times;
        clusterDetector.adcs1 = (*bestMatchPlane1).adcs;
        if (m_config.GetDetectorPlane(dp2) == true) {
          (*bestMatchPlane2).plane_coincidence = true;
          clusterDetector.id2 = (*bestMatchPlane2).id;
          clusterDetector.size2 = (*bestMatchPlane2).size;
          clusterDetector.adc2 = (*bestMatchPlane2).adc;
          clusterDetector.pos2 = (*bestMatchPlane2).pos;
          clusterDetector.pos2_utpc = (*bestMatchPlane2).pos_utpc;
          clusterDetector.pos2_charge2 = (*bestMatchPlane2).pos_charge2;
          clusterDetector.pos2_algo = (*bestMatchPlane2).pos_algo;
          clusterDetector.time2 = (*bestMatchPlane2).time;
          clusterDetector.time2_utpc = (*bestMatchPlane2).time_utpc;
          clusterDetector.time2_charge2 = (*bestMatchPlane2).time_charge2;
          clusterDetector.time2_algo = (*bestMatchPlane2).time_algo;
          clusterDetector.dt2 = clusterDetector.time2 - last_time2;
          clusterDetector.max_delta_time2 = (*bestMatchPlane2).max_delta_time;
          clusterDetector.max_missing_strip2 = (*bestMatchPlane2).max_missing_strip;
          clusterDetector.span_cluster2 = (*bestMatchPlane2).span_cluster;
          clusterDetector.strips2 = (*bestMatchPlane2).strips;
          clusterDetector.times2 = (*bestMatchPlane2).times;
          clusterDetector.adcs2 = (*bestMatchPlane2).adcs;
        }

        last_time0 = clusterDetector.time0;
        last_time1 = clusterDetector.time1;
        last_time2 = clusterDetector.time2;

        clusterDetector.delta_plane_0_1 = clusterDetector.time1 - clusterDetector.time0;
        clusterDetector.delta_plane_1_2 = clusterDetector.time2 - clusterDetector.time1;
        clusterDetector.delta_plane_0_2 = clusterDetector.time2 - clusterDetector.time0;

        if (m_config.pConditionCoincidence == "utpc") {
          clusterDetector.delta_plane_0_1 = clusterDetector.time1_utpc - clusterDetector.time0_utpc;
          clusterDetector.delta_plane_1_2 = clusterDetector.time2_utpc - clusterDetector.time1_utpc;
          clusterDetector.delta_plane_0_2 = clusterDetector.time2_utpc - clusterDetector.time0_utpc;
        } else if (m_config.pConditionCoincidence == "charge2") {
          clusterDetector.delta_plane_0_1 = clusterDetector.time1_charge2 - clusterDetector.time0_charge2;
          clusterDetector.delta_plane_1_2 = clusterDetector.time2_charge2 - clusterDetector.time1_charge2;
          clusterDetector.delta_plane_0_2 = clusterDetector.time2_charge2 - clusterDetector.time0_charge2;
        } else if (m_config.pConditionCoincidence == "algo") {
          clusterDetector.delta_plane_0_1 = clusterDetector.time1_algo - clusterDetector.time0_algo;
          clusterDetector.delta_plane_1_2 = clusterDetector.time2_algo - clusterDetector.time1_algo;
          clusterDetector.delta_plane_0_2 = clusterDetector.time2_algo - clusterDetector.time0_algo;
        }

        m_stats.SetStatsDetector("DeltaTimePlanes_0_1", det, std::fabs(clusterDetector.delta_plane_0_1));
        m_stats.SetStatsDetector("DeltaTimePlanes_1_2", det, std::fabs(clusterDetector.delta_plane_1_2));
        m_stats.SetStatsDetector("DeltaTimePlanes_0_2", det, std::fabs(clusterDetector.delta_plane_0_2));
        double ratio = 100 * (double)clusterDetector.adc0 / (double)clusterDetector.adc1;
        if (ratio > 100.0) {
          ratio = 100 * (double)clusterDetector.adc1 / (double)clusterDetector.adc0;
          m_stats.SetStatsDetector("ChargeRatio_1_0", det, ratio);
        } else {
          m_stats.SetStatsDetector("ChargeRatio_0_1", det, ratio);
        }
        if (m_config.GetDetectorPlane(dp2) == true) {
          double ratio = 100 * (double)clusterDetector.adc1 / (double)clusterDetector.adc2;
          if (ratio > 100.0) {
            ratio = 100 * (double)clusterDetector.adc2 / (double)clusterDetector.adc1;
            m_stats.SetStatsDetector("ChargeRatio_2_1", det, ratio);
          } else {
            m_stats.SetStatsDetector("ChargeRatio_1_2", det, ratio);
          }

          ratio = 100 * (double)clusterDetector.adc0 / (double)clusterDetector.adc2;
          if (ratio > 100.0) {
            ratio = 100 * (double)clusterDetector.adc2 / (double)clusterDetector.adc0;
            m_stats.SetStatsDetector("ChargeRatio_2_0", det, ratio);
          } else {
            m_stats.SetStatsDetector("ChargeRatio_0_2", det, ratio);
          }
        }
        m_stats.SetStatsDetector("ClusterCntDetector", det, 0);
        clusterCount++;

        if (m_config.GetDetectorPlane(dp2) == true) {
          LOG(TRACE) << "Common cluster " << (int)det << " id x/id y/id u: " << clusterDetector.id0 << "/" << clusterDetector.id1 << "/" << clusterDetector.id2;
          LOG(TRACE) << "pos x/pos y/pos u: " << clusterDetector.pos0 << "/" << clusterDetector.pos0 << "/" << clusterDetector.pos2;
          LOG(TRACE) << "time x/time y/time u: " << static_cast<unsigned long long>(clusterDetector.time0) << "/" << static_cast<unsigned long long>(clusterDetector.time1) << "/" << static_cast<unsigned long long>(clusterDetector.time2);
          LOG(TRACE) << "adc x/adc y/adc u: " << clusterDetector.adc0 << "/" << clusterDetector.adc1 << "/" << clusterDetector.adc2;
          LOG(TRACE) << "size x/size y/size u: " << clusterDetector.size0 << "/" << clusterDetector.size1 << "/" << clusterDetector.size2;
          LOG(TRACE) << "delta time planes 0 - 1: " << (int)clusterDetector.delta_plane_0_1;
          LOG(TRACE) << "delta time planes 1 - 2: " << (int)clusterDetector.delta_plane_1_2;
          LOG(TRACE) << "delta time planes 0 - 2: " << (int)clusterDetector.delta_plane_0_2;

        } else {
          LOG(TRACE) << "Common cluster " << (int)det << " id x/id y: " << clusterDetector.id0 << "/" << clusterDetector.id1;
          LOG(TRACE) << "pos x/pos y: " << clusterDetector.pos0 << "/" << clusterDetector.pos1;
          LOG(TRACE) << "time x/time y: " << static_cast<unsigned long long>(clusterDetector.time0) << "/" << static_cast<unsigned long long>(clusterDetector.time1);
          LOG(TRACE) << "adc x/adc y: " << clusterDetector.adc0 << "/" << clusterDetector.adc1;
          LOG(TRACE) << "size x/size y: " << clusterDetector.size0 << "/" << clusterDetector.size1;
          LOG(TRACE) << "delta time planes: " << (int)clusterDetector.delta_plane_0_1;
        }
        m_clusters_detector[det].emplace_back(std::move(clusterDetector));
      }
    }
  }

  return clusterCount;
}

//====================================================================================================================
int Clusterer::MatchClustersDetector_HighMultiplicity(uint8_t det) {
  int clusterCount = 0;
  auto dp0 = std::make_pair(det, 0);
  auto dp1 = std::make_pair(det, 1);

  bool found1 = false;
  ClusterVectorPlane::iterator itStartPlane1 = begin(m_clusters[dp1]);
  for (auto& c0 : m_clusters[dp0]) {
    double delta_t1 = std::numeric_limits<double>::max();
    double lastDelta_t1 = std::numeric_limits<double>::max();
    double minDelta1 = std::numeric_limits<double>::max();
    bool isFirstMatch1 = true;
    ClusterVectorPlane::iterator bestMatchPlane1 = end(m_clusters[dp1]);
    int bestMatchIdx = -1;
    for (ClusterVectorPlane::iterator c1 = itStartPlane1; c1 != end(m_clusters[dp1]); ++c1) {
      double chargeRatio1 = (double)(c0).adc / (double)(*c1).adc;
      lastDelta_t1 = delta_t1;
      delta_t1 = (*c1).time - c0.time;
      if (m_config.pConditionCoincidence == "utpc") {
        delta_t1 = (*c1).time_utpc - c0.time_utpc;
      } else if (m_config.pConditionCoincidence == "charge2") {
        delta_t1 = (*c1).time_charge2 - c0.time_charge2;
      }
      if (std::fabs(delta_t1) <= m_config.pDeltaTimePlanes[m_config.pDets[det]]) {
        if (isFirstMatch1) {
          itStartPlane1 = c1;
          isFirstMatch1 = false;
        }
      }
      if (std::fabs(delta_t1) > m_config.pDeltaTimePlanes[m_config.pDets[det]] && std::fabs(delta_t1) > std::fabs(lastDelta_t1)) {
        break;
      }
      if (chargeRatio1 >= m_config.pChargeRatioLower[m_config.pDets[det]] && chargeRatio1 <= m_config.pChargeRatioUpper[m_config.pDets[det]] && std::fabs(delta_t1) <= m_config.pDeltaTimePlanes[m_config.pDets[det]] && (c0.size + (*c1).size >= m_config.pCoincidentClusterSize[m_config.pDets[det]])) {
        if (std::fabs(delta_t1) < minDelta1 && std::fabs(delta_t1) <= m_config.pDeltaTimePlanes[m_config.pDets[det]]) {
          if ((*c1).plane_coincidence == false) {
            minDelta1 = std::fabs(delta_t1);
            bestMatchIdx = m_clusters_detector[det].size();
            bestMatchPlane1 = c1;
          }
        }

        m_cluster_detector_id++;
        ClusterDetector clusterDetector;
        clusterDetector.id = m_cluster_detector_id;
        clusterDetector.det = det;
        clusterDetector.event_counter = (*c1).event_counter;
        clusterDetector.id0 = c0.id;
        clusterDetector.id1 = (*c1).id;
        clusterDetector.id2 = 0;
        clusterDetector.size0 = c0.size;
        clusterDetector.size1 = (*c1).size;
        clusterDetector.size2 = 0;
        clusterDetector.adc0 = c0.adc;
        clusterDetector.adc1 = (*c1).adc;
        clusterDetector.adc2 = 0;
        clusterDetector.pos0 = c0.pos;
        clusterDetector.pos1 = (*c1).pos;
        clusterDetector.pos2 = 0;
        clusterDetector.pos0_utpc = c0.pos_utpc;
        clusterDetector.pos1_utpc = (*c1).pos_utpc;
        clusterDetector.pos2_utpc = 0;
        clusterDetector.pos0_charge2 = c0.pos_charge2;
        clusterDetector.pos1_charge2 = (*c1).pos_charge2;
        clusterDetector.pos2_charge2 = 0;
        clusterDetector.pos0_algo = c0.pos_algo;
        clusterDetector.pos1_algo = (*c1).pos_algo;
        clusterDetector.pos2_algo = 0;
        clusterDetector.time0 = c0.time;
        clusterDetector.time1 = (*c1).time;
        clusterDetector.time2 = 0;
        clusterDetector.time0_utpc = c0.time_utpc;
        clusterDetector.time1_utpc = (*c1).time_utpc;
        clusterDetector.time2_utpc = 0;
        clusterDetector.time0_charge2 = c0.time_charge2;
        clusterDetector.time1_charge2 = (*c1).time_charge2;
        clusterDetector.time2_charge2 = 0;
        clusterDetector.time0_algo = c0.time_algo;
        clusterDetector.time1_algo = (*c1).time_algo;
        clusterDetector.time2_algo = 0;
        clusterDetector.dt0 = clusterDetector.time0 - last_time0;
        clusterDetector.dt1 = clusterDetector.time1 - last_time1;
        clusterDetector.dt2 = 0;
        clusterDetector.max_delta_time0 = c0.max_delta_time;
        clusterDetector.max_delta_time1 = (*c1).max_delta_time;
        clusterDetector.max_delta_time2 = 0;
        clusterDetector.max_missing_strip0 = c0.max_missing_strip;
        clusterDetector.max_missing_strip1 = (*c1).max_missing_strip;
        clusterDetector.max_missing_strip2 = 0;
        clusterDetector.span_cluster0 = c0.span_cluster;
        clusterDetector.span_cluster1 = (*c1).span_cluster;
        clusterDetector.span_cluster2 = 0;
        clusterDetector.strips0 = c0.strips;
        clusterDetector.times0 = c0.times;
        clusterDetector.adcs0 = c0.adcs;
        clusterDetector.strips1 = (*c1).strips;
        clusterDetector.times1 = (*c1).times;
        clusterDetector.adcs1 = (*c1).adcs;

        last_time0 = clusterDetector.time0;
        last_time1 = clusterDetector.time1;

        clusterDetector.delta_plane_0_1 = clusterDetector.time1 - clusterDetector.time0;
        if (m_config.pConditionCoincidence == "utpc") {
          clusterDetector.delta_plane_0_1 = clusterDetector.time1_utpc - clusterDetector.time0_utpc;
        } else if (m_config.pConditionCoincidence == "charge2") {
          clusterDetector.delta_plane_0_1 = clusterDetector.time1_charge2 - clusterDetector.time0_charge2;
        }

        m_stats.SetStatsDetector("DeltaTimePlanes_0_1", det, std::fabs(clusterDetector.delta_plane_0_1));
        double ratio = 100 * (double)clusterDetector.adc0 / (double)clusterDetector.adc1;
        if (ratio > 100.0) {
          ratio = 100 * (double)clusterDetector.adc1 / (double)clusterDetector.adc0;
          m_stats.SetStatsDetector("ChargeRatio_1_0", det, ratio);
        } else {
          m_stats.SetStatsDetector("ChargeRatio_0_1", det, ratio);
        }
        m_stats.SetStatsDetector("ClusterCntDetector", det, 0);
        LOG(TRACE) << "Common cluster " << (int)det << " id x/id y: " << clusterDetector.id0 << "/" << clusterDetector.id1;
        LOG(TRACE) << "pos x/pos y: " << clusterDetector.pos0 << "/" << clusterDetector.pos0;
        LOG(TRACE) << "time x/time y: " << static_cast<unsigned long long>(clusterDetector.time0) << "/" << static_cast<unsigned long long>(clusterDetector.time1);
        LOG(TRACE) << "adc x/adc y: " << clusterDetector.adc0 << "/" << clusterDetector.adc1;
        LOG(TRACE) << "size x/size y: " << clusterDetector.size0 << "/" << clusterDetector.size1;
        LOG(TRACE) << "delta time planes: " << (int)clusterDetector.delta_plane_0_1;

        m_clusters_detector[det].emplace_back(std::move(clusterDetector));
        clusterCount++;
      }
    }
    if (bestMatchIdx != -1) {
      m_clusters_detector[det][bestMatchIdx].id2 = 1;
      (*bestMatchPlane1).plane_coincidence = true;
    }
  }
  return clusterCount;
}

//====================================================================================================================
int Clusterer::FillClustersDetector(uint8_t det, uint8_t plane) {
  corryvreckan::Log::setSection("Clusterer");
  int clusterCount = 0;
  auto dp = std::make_pair(det, plane);

  bool found = true;
  ClusterVectorPlane::iterator itStartPlane = begin(m_clusters[dp]);

  for (auto& c : m_clusters[dp]) {
    m_cluster_detector_id++;
    ClusterDetector clusterDetector;
    clusterDetector.id = m_cluster_detector_id;
    clusterDetector.det = det;
    clusterDetector.event_counter = c.event_counter;
    clusterDetector.id0 = 0;
    clusterDetector.id1 = 0;
    clusterDetector.id2 = 0;
    clusterDetector.size0 = 0;
    clusterDetector.size1 = 0;
    clusterDetector.size2 = 0;
    clusterDetector.adc0 = 0;
    clusterDetector.adc1 = 0;
    clusterDetector.adc2 = 0;
    clusterDetector.pos0 = 0;
    clusterDetector.pos1 = 0;
    clusterDetector.pos2 = 0;
    clusterDetector.pos0_utpc = 0;
    clusterDetector.pos1_utpc = 0;
    clusterDetector.pos2_utpc = 0;
    clusterDetector.pos0_charge2 = 0;
    clusterDetector.pos1_charge2 = 0;
    clusterDetector.pos2_charge2 = 0;
    clusterDetector.pos0_algo = 0;
    clusterDetector.pos1_algo = 0;
    clusterDetector.pos2_algo = 0;
    clusterDetector.time0 = 0;
    clusterDetector.time1 = 0;
    clusterDetector.time2 = 0;
    clusterDetector.time0_utpc = 0;
    clusterDetector.time1_utpc = 0;
    clusterDetector.time2_utpc = 0;
    clusterDetector.time0_charge2 = 0;
    clusterDetector.time1_charge2 = 0;
    clusterDetector.time2_charge2 = 0;
    clusterDetector.time0_algo = 0;
    clusterDetector.time1_algo = 0;
    clusterDetector.time2_algo = 0;
    clusterDetector.dt0 = 0;
    clusterDetector.dt1 = 0;
    clusterDetector.dt2 = 0;
    clusterDetector.max_delta_time0 = 0;
    clusterDetector.max_delta_time1 = 0;
    clusterDetector.max_delta_time2 = 0;
    clusterDetector.max_missing_strip0 = 0;
    clusterDetector.max_missing_strip1 = 0;
    clusterDetector.max_missing_strip2 = 0;
    clusterDetector.span_cluster0 = 0;
    clusterDetector.span_cluster1 = 0;
    clusterDetector.span_cluster2 = 0; 

    if (plane == 0) {
      clusterDetector.id0 = c.id;
      clusterDetector.size0 = c.size;
      clusterDetector.adc0 = c.adc;
      clusterDetector.pos0 = c.pos;
      clusterDetector.pos0_utpc = c.pos_utpc;
      clusterDetector.pos0_charge2 = c.pos_charge2;
      clusterDetector.pos0_algo = c.pos_algo;
      clusterDetector.time0 = c.time;
      clusterDetector.time0_utpc = c.time_utpc;
      clusterDetector.time0_charge2 = c.time_charge2;
      clusterDetector.time0_algo = c.time_algo;
      clusterDetector.dt0 = clusterDetector.time0 - last_time0;
      clusterDetector.max_delta_time0 = c.max_delta_time;
      clusterDetector.max_missing_strip0 = c.max_missing_strip;
      clusterDetector.span_cluster0 = c.span_cluster;
      clusterDetector.strips0 = c.strips;
      clusterDetector.times0 = c.times;
      clusterDetector.adcs0 = c.adcs;
      last_time0 = clusterDetector.time0;
    } else if (plane == 1) {
      clusterDetector.id1 = c.id;
      clusterDetector.size1 = c.size;
      clusterDetector.adc1 = c.adc;
      clusterDetector.pos1 = c.pos;
      clusterDetector.pos1_utpc = c.pos_utpc;
      clusterDetector.pos1_charge2 = c.pos_charge2;
      clusterDetector.pos1_algo = c.pos_algo;
      clusterDetector.time1 = c.time;
      clusterDetector.time1_utpc = c.time_utpc;
      clusterDetector.time1_charge2 = c.time_charge2;
      clusterDetector.time1_algo = c.time_algo;
      clusterDetector.dt1 = clusterDetector.time1 - last_time1;
      clusterDetector.max_delta_time1 = c.max_delta_time;
      clusterDetector.max_missing_strip1 = c.max_missing_strip;
      clusterDetector.span_cluster1 = c.span_cluster;
      clusterDetector.strips1 = c.strips;
      clusterDetector.times1 = c.times;
      clusterDetector.adcs1 = c.adcs;
      last_time1 = clusterDetector.time1;
    } else {
      clusterDetector.id2 = c.id;
      clusterDetector.size2 = c.size;
      clusterDetector.adc2 = c.adc;
      clusterDetector.pos2 = c.pos;
      clusterDetector.pos2_utpc = c.pos_utpc;
      clusterDetector.pos2_charge2 = c.pos_charge2;
      clusterDetector.pos2_algo = c.pos_algo;
      clusterDetector.time2 = c.time;
      clusterDetector.time2_utpc = c.time_utpc;
      clusterDetector.time2_charge2 = c.time_charge2;
      clusterDetector.time2_algo = c.time_algo;
      clusterDetector.dt2 = clusterDetector.time2 - last_time2;
      clusterDetector.max_delta_time2 = c.max_delta_time;
      clusterDetector.max_missing_strip2 = c.max_missing_strip;
      clusterDetector.span_cluster2 = c.span_cluster;
      clusterDetector.strips2 = c.strips;
      clusterDetector.times2 = c.times;
      clusterDetector.adcs2 = c.adcs;
      last_time2 = clusterDetector.time2;
    }

    m_stats.SetStatsDetector("ClusterCntDetector", det, 0);
    clusterCount++;
    m_clusters_detector[det].emplace_back(std::move(clusterDetector));
  }
  return clusterCount;
}

void Clusterer::AnalyzeClustersPlane(std::pair<uint8_t, uint8_t> dp) {
  if (ChooseHitsToBeClustered(dp) == false && m_hits[dp].empty()) {
    return;
  }
  int cnt = ClusterByTime(dp);
  corryvreckan::Log::setSection("Clusterer");
  LOG(TRACE) << "Number of clusters in detector " << (int)std::get<0>(dp) << " plane " << (int)std::get<1>(dp) << ": " << cnt;
  m_hits[dp].clear();
}

void Clusterer::AnalyzeClustersDetector(uint8_t det) {
  if(m_config.pIsPads[det]) {
    return;
  }
  int cnt = 0;
  auto dp0 = std::make_pair(det, 0);
  auto dp1 = std::make_pair(det, 1);
  auto dp2 = std::make_pair(det, 2);

  int num_planes = 0;
  if (m_config.GetDetectorPlane(dp0) == true) {
    num_planes++;
  }
  if (m_config.GetDetectorPlane(dp1) == true) {
    num_planes++;
  }
  if (m_config.GetDetectorPlane(dp2) == true) {
    num_planes++;
  }
  if (num_planes == 0) {
    return;
  }

  if (m_config.GetDetectorPlane(dp0) == true) {
    if (ChooseClustersToBeMatched(dp0) == false && m_clusters[dp0].empty()) {
      return;
    }
  }
  if (m_config.GetDetectorPlane(dp1) == true) {
    if (ChooseClustersToBeMatched(dp1) == false && m_clusters[dp1].empty()) {
      return;
    }
  }
  if (m_config.GetDetectorPlane(dp2) == true) {
    if (ChooseClustersToBeMatched(dp2) == false && m_clusters[dp2].empty()) {
      return;
    }
  }

  if (num_planes > 1) {
    if (m_config.pHighMultiplicity == true) {
      cnt = MatchClustersDetector_HighMultiplicity(det);
    } else {
      cnt = MatchClustersDetector(det);
    }
  } else {
    if (m_config.GetDetectorPlane(dp0) == true) {
      cnt = FillClustersDetector(det, 0);
    } else if (m_config.GetDetectorPlane(dp1) == true) {
      cnt = FillClustersDetector(det, 1);
    } else if (m_config.GetDetectorPlane(dp2) == true) {
      cnt = FillClustersDetector(det, 2);
    }
  }

  if (m_config.pSaveWhat == 10 || m_config.pSaveWhat == 11 || m_config.pSaveWhat == 110 || m_config.pSaveWhat == 111) {
    m_rootFile->SaveClustersPlane(std::move(m_clusters[dp0]));
    m_rootFile->SaveClustersPlane(std::move(m_clusters[dp1]));
    m_rootFile->SaveClustersPlane(std::move(m_clusters[dp2]));
  }

  if (m_config.pSaveWhat >= 100) {
    m_rootFile->SaveClustersDetector(std::move(m_clusters_detector[det]));
  }
  m_clusters[std::make_pair(det, 0)].clear();
  m_clusters[std::make_pair(det, 1)].clear();
  m_clusters[std::make_pair(det, 2)].clear();
  m_clusters_detector[det].clear();
}

//====================================================================================================================
bool Clusterer::ChooseHitsToBeClustered(std::pair<uint8_t, uint8_t> dp) {
  // std::pair<uint8_t, uint8_t> dp = std::make_pair(det, plane);
  double timeReadyToCluster = m_stats.GetLowestCommonTriggerTimestampPlane(dp);
  // Nothing to cluster, newHits vector empty
  if (m_hits_new[dp].empty()) {
    return false;
  }

  auto theMin = std::min_element(m_hits_new[dp].begin(), m_hits_new[dp].end(), [](const HitTuple& t1, const HitTuple& t2) { return std::get<0>(t1) < std::get<0>(t2); });

  // Nothing to cluster, tuples in newHits vector too recent
  if (std::get<0>(*theMin) > timeReadyToCluster) {
    //(smallest timestamp larger than
    // m_stats.GetLowestCommonTriggerTimestampPlane(dp)) Will be clustered
    // later
    return false;
  }

  // Sort vector newHits
  std::sort(begin(m_hits_new[dp]), end(m_hits_new[dp]), [](const HitTuple& t1, const HitTuple& t2) { return std::get<0>(t1) < std::get<0>(t2); });

  // First tuple with timestamp larger than
  // m_stats.GetLowestCommonTriggerTimestampPlane(dp)
  auto it = std::upper_bound(m_hits_new[dp].begin(), m_hits_new[dp].end(), std::make_tuple(m_stats.GetLowestCommonTriggerTimestampPlane(dp), 0, 0, 0, 0), [](const HitTuple& t1, const HitTuple& t2) { return std::get<0>(t1) < std::get<0>(t2); });

  // Find elements in vector that could still be part of a cluster,
  // since they are close in time to
  // m_stats.GetLowestCommonTriggerTimestampPlane(dp)
  while (it != m_hits_new[dp].end()) {
    if (std::get<0>(*it) - timeReadyToCluster > m_config.pDeltaTimeHits[m_config.pDets[dp.first]]) {
      break;
    }
    timeReadyToCluster = std::get<0>(*it);
    ++it;
  }
  int index = std::distance(m_hits_new[dp].begin(), it);
  // Insert the data that is ready to be clustered from newHits into hits
  m_hits[dp].insert(m_hits[dp].end(), std::make_move_iterator(m_hits_new[dp].begin()), std::make_move_iterator(m_hits_new[dp].begin() + index));
  // Delete the data from newHits
  m_hits_new[dp].erase(m_hits_new[dp].begin(), m_hits_new[dp].begin() + index);

  return true;
}

bool Clusterer::ChooseClustersToBeMatched(std::pair<uint8_t, uint8_t> dp) {
  int index = 0;
  // std::pair<uint8_t, uint8_t> dp = std::make_pair(det, plane);
  double timeReadyToMatch = m_stats.GetLowestCommonTriggerTimestampPlane(dp);

  // Nothing to match, newClusters vector empty
  if (m_clusters_new[dp].empty()) {
    return false;
  }
  if (m_config.pConditionCoincidence == "utpc") {
    auto theMin = std::min_element(m_clusters_new[dp].begin(), m_clusters_new[dp].end(), [](const ClusterPlane& t1, const ClusterPlane& t2) { return t1.time_utpc < t2.time_utpc; });

    // Nothing to cluster, clusters in newClusters vector too recent
    if ((*theMin).time_utpc > timeReadyToMatch) {
      //(smallest time larger than timeReadyToMatch)
      // Will be matched later
      return false;
    }

    // Sort vector newClusters based on time
    std::sort(begin(m_clusters_new[dp]), end(m_clusters_new[dp]), [](const ClusterPlane& t1, const ClusterPlane& t2) { return t1.time_utpc < t2.time_utpc; });

    ClusterPlane theCluster;
    theCluster.time_utpc = timeReadyToMatch;

    // First ClusterPlane with time bigger than timeReadyToMatch
    auto it = std::upper_bound(m_clusters_new[dp].begin(), m_clusters_new[dp].end(), theCluster, [](const ClusterPlane& t1, const ClusterPlane& t2) { return t1.time_utpc < t2.time_utpc; });

    // Find elements in vector that could still be matched with another
    // cluster since they are close in time to timeReadyToMatch
    while (it != m_clusters_new[dp].end()) {
      if ((*it).time_utpc - timeReadyToMatch > m_config.pDeltaTimeHits[m_config.pDets[dp.first]]) {
        break;
      }
      timeReadyToMatch = (*it).time_utpc;
      ++it;
    }
    index = std::distance(m_clusters_new[dp].begin(), it);
  } else if (m_config.pConditionCoincidence == "charge2") {
    auto theMin = std::min_element(m_clusters_new[dp].begin(), m_clusters_new[dp].end(), [](const ClusterPlane& t1, const ClusterPlane& t2) { return t1.time_charge2 < t2.time_charge2; });

    // Nothing to cluster, clusters in newClusters vector too recent
    if ((*theMin).time_charge2 > timeReadyToMatch) {
      //(smallest time larger than timeReadyToMatch)
      // Will be matched later
      return false;
    }

    // Sort vector newClusters based on time
    std::sort(begin(m_clusters_new[dp]), end(m_clusters_new[dp]), [](const ClusterPlane& t1, const ClusterPlane& t2) { return t1.time_charge2 < t2.time_charge2; });

    ClusterPlane theCluster;
    theCluster.time_charge2 = timeReadyToMatch;

    // First ClusterPlane with time that bigger than timeReadyToMatch
    auto it = std::upper_bound(m_clusters_new[dp].begin(), m_clusters_new[dp].end(), theCluster, [](const ClusterPlane& t1, const ClusterPlane& t2) { return t1.time_charge2 < t2.time_charge2; });

    // Find elements in vector that could still be matched with another
    // cluster since they are close in time to timeReadyToMatch
    while (it != m_clusters_new[dp].end()) {
      if ((*it).time_charge2 - timeReadyToMatch > m_config.pDeltaTimeHits[m_config.pDets[dp.first]]) {
        break;
      }
      timeReadyToMatch = (*it).time_charge2;
      ++it;
    }
    index = std::distance(m_clusters_new[dp].begin(), it);
  } else {
    auto theMin = std::min_element(m_clusters_new[dp].begin(), m_clusters_new[dp].end(), [](const ClusterPlane& t1, const ClusterPlane& t2) { return t1.time < t2.time; });

    // Nothing to cluster, clusters in newClusters vector too recent
    if ((*theMin).time > timeReadyToMatch) {
      //(smallest time larger than timeReadyToMatch)
      // Will be matched later
      return false;
    }

    // Sort vector newClusters based on time
    std::sort(begin(m_clusters_new[dp]), end(m_clusters_new[dp]), [](const ClusterPlane& t1, const ClusterPlane& t2) { return t1.time < t2.time; });

    ClusterPlane theCluster;
    theCluster.time = timeReadyToMatch;

    // First ClusterPlane with time that bigger than timeReadyToMatch
    auto it = std::upper_bound(m_clusters_new[dp].begin(), m_clusters_new[dp].end(), theCluster, [](const ClusterPlane& t1, const ClusterPlane& t2) { return t1.time < t2.time; });

    // Find elements in vector that could still be matched with another
    // cluster since they are close in time to timeReadyToMatch
    while (it != m_clusters_new[dp].end()) {
      if ((*it).time - timeReadyToMatch > m_config.pDeltaTimeHits[m_config.pDets[dp.first]]) {
        break;
      }
      timeReadyToMatch = (*it).time;
      ++it;
    }
    index = std::distance(m_clusters_new[dp].begin(), it);
  }
  // Insert the clusters that are ready to be matched from newClusters into
  // clusters
  m_clusters[dp].insert(m_clusters[dp].end(), std::make_move_iterator(m_clusters_new[dp].begin()), std::make_move_iterator(m_clusters_new[dp].begin() + index));
  // Delete the clusters from newClusters
  m_clusters_new[dp].erase(m_clusters_new[dp].begin(), m_clusters_new[dp].begin() + index);

  return true;
}

void Clusterer::FinishAnalysis() {
  double ts = 0;
  for (auto const& fec : m_config.pFecs) {
    if (ts < m_stats.GetMaxTriggerTimestamp(fec)) {
      ts = m_stats.GetMaxTriggerTimestamp(fec);
    }
  }
  for (auto const& det : m_config.pDets) {
    auto dp0 = std::make_pair(det.first, 0);
    auto dp1 = std::make_pair(det.first, 1);
    auto dp2 = std::make_pair(det.first, 2);

    // Set the largest timestamp of plane to detector
    // cluster all remaining data in plane
    m_stats.SetLowestCommonTriggerTimestampPlane(dp0, ts);
    m_stats.SetLowestCommonTriggerTimestampPlane(dp1, ts);
    m_stats.SetLowestCommonTriggerTimestampPlane(dp2, ts);
    m_stats.SetLowestCommonTriggerTimestampDet(det.first, ts);

    AnalyzeClustersPlane(dp0);
    AnalyzeClustersPlane(dp1);
    AnalyzeClustersPlane(dp2);

    AnalyzeClustersDetector(det.first);

    if (m_config.pSaveWhat % 2 == 1) {
      m_rootFile->SaveHits();
    }
  }

  if (m_config.pSaveWhat >= 10) {
    m_stats.PrintClusterStats(m_config);
  }
  m_stats.PrintFECStats(m_config);
}

void Clusterer::SaveDate(double the_seconds_start, std::string the_date_start, double the_seconds_end, std::string the_date_end, uint64_t triggers) {
  m_rootFile->SaveDate(the_seconds_start, the_date_start, the_seconds_end, the_date_end, triggers);
}

void Clusterer::FillCalibHistos(uint16_t fec, uint8_t vmm, uint8_t ch, float adc, float adc_corrected, float chip_time, float chip_time_corrected) {
  if (m_config.useCalibration && m_config.calibrationHistogram) {
    m_rootFile->FillCalibHistos(fec, vmm, ch, adc, adc_corrected, chip_time, chip_time_corrected);
  }
}