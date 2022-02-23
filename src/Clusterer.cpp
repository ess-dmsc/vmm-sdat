#include "Clusterer.h"
#include <algorithm>
#include <cmath>
#include <parser/Trace.h>

#include <chrono>
#include <functional>
#include <iomanip>
#define UNUSED __attribute__((unused))

//#undef TRC_LEVEL
//#define TRC_LEVEL TRC_L_DEB

auto now = std::chrono::steady_clock::now;

auto timethis(std::function<void()> thunk)
    -> decltype((now() - now()).count()) {
  auto start = now();
  thunk();
  auto stop = now();
  return (stop - start).count();
}

Clusterer::Clusterer(Configuration &config, Statistics &stats)
    : m_config(config), m_stats(stats) {
  m_rootFile = RootFile::GetInstance(config);
}

Clusterer::~Clusterer() { RootFile::Dispose(); }

//====================================================================================================================
bool Clusterer::AnalyzeHits(double srsTimestamp, uint8_t fecId, uint8_t vmmId,
                            uint16_t chNo, uint16_t bcid, uint16_t tdc,
                            uint16_t adc, bool overThresholdFlag,
                            double chipTime) {

  int pos0 = m_config.pPositions0[fecId][vmmId][chNo];
  int pos1 = m_config.pPositions1[fecId][vmmId][chNo];
  // std::cout << "pos0/pos1 " << pos0 << " " << pos1 << std::endl;
  if (pos0 == -1) {
    DTRACE(DEB, "\t\tDetector or Plane not defined for FEC %d and vmmId %d!\n",
           (int)fecId, (int)vmmId);
    return true;
  }

  if (m_config.pShowStats && m_config.pDataFormat == "SRS") {
    // Biggest possible time should be:
    // from FEC: 2^42-1 = 0x3FFFFFFFFFF clock cycles
    // converted to ns: 0x3FFFFFFFFFF * bc period = 109951162777575 ns
    // added 31 offsets of 4096 * bc period ns: 3174400
    // total: 109951165951975 ns
    if (srsTimestamp > 109951165951975) {
      m_stats.IncrementCounter("TimestampTooLarge", fecId);
      DTRACE(
          DEB,
          "\t\tTimestamp %llu larger than 42 bit and 31 times trigger periodd "
          "for FEC %d and vmmId %d!\n",
          static_cast<uint64_t>(srsTimestamp), (int)fecId, (int)vmmId);
    }

    if (srsTimestamp < m_stats.GetMaxTriggerTimestamp(fecId)) {
      // 42 bit: 0x1FFFFFFFFFF
      // 32 bit: 0xFFFFFFFF
      if (m_stats.GetMaxTriggerTimestamp(fecId) >
          0x1FFFFFFFFFF + srsTimestamp) {
        m_stats.IncrementCounter("TimestampOverflow", fecId);
        DTRACE(DEB,
               "\n*********************************** OVERFLOW  fecId %d, "
               "m_hitNr %d, "
               "srsTimestamp %llu, old srsTimestamp %llu\n",
               fecId, m_hitNr, static_cast<uint64_t>(srsTimestamp),
               static_cast<uint64_t>(m_stats.GetMaxTriggerTimestamp(fecId)));

      } else {
        m_stats.IncrementCounter("TimestampOrderError", fecId);
        DTRACE(DEB,
               "\n*********************************** TIME ERROR  fecId %d, "
               "m_hitNr %d, "
               "srsTimestamp %llu, old srsTimestamp %llu\n",
               fecId, m_hitNr, static_cast<uint64_t>(srsTimestamp),
               static_cast<uint64_t>(m_stats.GetMaxTriggerTimestamp(fecId)));
      }
    }

    double remainder = std::fmod(m_stats.GetDeltaTriggerTimestamp(fecId),
                                 m_config.pOffsetPeriod);
    if (remainder > 0) {
      m_stats.IncrementCounter("TriggerPeriodError", fecId);
      uint64_t offset =
          m_stats.GetDeltaTriggerTimestamp(fecId) / m_config.pOffsetPeriod;
      DTRACE(
          DEB,
          "\n******* ERROR: SRS timestamp wrong increment: fec "
          "%d,vmmId %d, chNo %d, hit %d, "
          "offset period %f, offset %llu, remainder %llu, new time %llu, old "
          "time %llu\n",
          fecId, vmmId, chNo, m_hitNr, m_config.pOffsetPeriod, offset,
          static_cast<uint64_t>(remainder), static_cast<uint64_t>(srsTimestamp),
          static_cast<uint64_t>(m_stats.GetOldTriggerTimestamp(fecId)));
    }
  }

  bool newData = false;

  if (m_config.pDataFormat == "SRS") {
    double factor = 16.0;
    if (srsTimestamp >= m_stats.GetOldTriggerTimestamp(fecId) +
                            factor * m_config.pOffsetPeriod) {
      if (m_config.pShowStats) {
        m_stats.SetDeltaTriggerTimestamp(
            fecId, srsTimestamp - m_stats.GetOldTriggerTimestamp(fecId));
      }
      newData = true;
    }
  } else if (m_config.pDataFormat == "ESS") {
    double buffer_interval_ns = 1000000.0;
    if (srsTimestamp >=
        m_stats.GetOldTriggerTimestamp(fecId) + buffer_interval_ns) {
      newData = true;
    }
  }

  if (newData) {
    if (m_config.pSaveWhat % 2 == 1) {
      m_rootFile->SaveHits();
    }

    if (m_config.pSaveWhat >= 10) {
      uint64_t ts = 0;
      for (auto const &fec : m_config.pFecs) {
        if (ts == 0 || ts > m_stats.GetOldTriggerTimestamp(fec)) {
          ts = m_stats.GetOldTriggerTimestamp(fec);
        }
      }

      for (auto const &det : m_config.pDets) {
        auto dp0 = std::make_pair(det.first, 0);
        auto dp1 = std::make_pair(det.first, 1);
        if (m_stats.GetLowestCommonTriggerTimestampDet(det.first) < ts) {
          m_stats.SetLowestCommonTriggerTimestampPlane(dp0, ts);
          m_stats.SetLowestCommonTriggerTimestampPlane(dp1, ts);
          m_stats.SetLowestCommonTriggerTimestampDet(det.first, ts);
          AnalyzeClustersPlane(dp0);
          AnalyzeClustersPlane(dp1);
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
  if (m_config.pSaveWhat % 2 == 1) {
    if (std::find(m_config.pSaveHits.begin(), m_config.pSaveHits.end(), det) !=
        m_config.pSaveHits.end()) {

      if (m_config.pIsPads[det]) {
        Hit theHit0;
        theHit0.id = m_hitNr;
        theHit0.det = det;
        theHit0.plane = 0;
        theHit0.fec = fecId;
        theHit0.vmm = vmmId;
        theHit0.readout_time = srsTimestamp;
        theHit0.ch = chNo;
        theHit0.pos = (uint16_t)pos0;
        theHit0.bcid = bcid;
        theHit0.tdc = tdc;
        theHit0.adc = adc;
        theHit0.over_threshold = overThresholdFlag;
        theHit0.chip_time = chipTime;
        theHit0.time = totalTime;
        m_rootFile->AddHits(std::move(theHit0));

        Hit theHit1;
        theHit1.id = m_hitNr;
        theHit1.det = det;
        theHit1.plane = 1;
        theHit1.fec = fecId;
        theHit1.vmm = vmmId;
        theHit1.readout_time = srsTimestamp;
        theHit1.ch = chNo;
        theHit1.pos = (uint16_t)pos1;
        theHit1.bcid = bcid;
        theHit1.tdc = tdc;
        theHit1.adc = adc;
        theHit1.over_threshold = overThresholdFlag;
        theHit1.chip_time = chipTime;
        theHit1.time = totalTime;
        m_rootFile->AddHits(std::move(theHit1));
      } else {
        Hit theHit;
        theHit.id = m_hitNr;
        theHit.det = det;
        theHit.plane = plane;
        theHit.fec = fecId;
        theHit.vmm = vmmId;
        theHit.readout_time = srsTimestamp;
        theHit.ch = chNo;
        theHit.pos = (uint16_t)pos0;
        theHit.bcid = bcid;
        theHit.tdc = tdc;
        theHit.adc = adc;
        theHit.over_threshold = overThresholdFlag;
        theHit.chip_time = chipTime;
        theHit.time = totalTime;
        m_rootFile->AddHits(std::move(theHit));
      }
    }
  }
  if (overThresholdFlag) {
    // keep the overThresholdFlag as bit 15 of the ADC
    adc = adc + 32768;
  }
  if (m_config.pADCThreshold < 0) {
    if (overThresholdFlag) {
      if (m_config.pIsPads[det]) {
        m_hits_new[std::make_pair(det, 0)].emplace_back(
            totalTime, (uint16_t)pos0, adc, (uint16_t)pos1);
        m_hits_new[std::make_pair(det, 1)].emplace_back(
            totalTime, (uint16_t)pos0, adc, (uint16_t)pos1);
      } else {
        m_hits_new[std::make_pair(det, plane)].emplace_back(
            totalTime, (uint16_t)pos0, adc, (uint16_t)pos1);
      }
    }
  } else {
    if ((adc >= m_config.pADCThreshold)) {
      if (m_config.pIsPads[det]) {
        m_hits_new[std::make_pair(det, 0)].emplace_back(
            totalTime, (uint16_t)pos0, adc, (uint16_t)pos1);
        m_hits_new[std::make_pair(det, 1)].emplace_back(
            totalTime, (uint16_t)pos0, adc, (uint16_t)pos1);
      } else {
        m_hits_new[std::make_pair(det, plane)].emplace_back(
            totalTime, (uint16_t)pos0, adc, (uint16_t)pos1);
      }
    }
  }

  if (m_stats.GetDeltaTriggerTimestamp(fecId) > 0) {
    DTRACE(DEB, "\n\tTriggerTimestamp %llu [ns]\n",
           static_cast<uint64_t>(srsTimestamp));
    DTRACE(DEB, "\tTime since last trigger %f us (%.4f kHz)\n",
           m_stats.GetDeltaTriggerTimestamp(fecId) * 0.001,
           (double)(1000000 / m_stats.GetDeltaTriggerTimestamp(fecId)));
  }

  if (m_oldFecId != fecId || newData) {
    DTRACE(DEB, "\tfecId  %d\n", fecId);
  }
  if (m_oldVmmId != vmmId || newData) {
    DTRACE(DEB, "\tDetector %d, plane %d, vmmId  %d\n", (int)det, (int)plane,
           vmmId);
  }
  DTRACE(DEB, "\t\tId0 %d, id1 %d (chNo  %d) - overThresholdFlag %d\n", pos0,
         pos1, chNo, (int)overThresholdFlag);
  DTRACE(DEB, "\t\t\tbcid %d, tdc %d, adc %d\n", bcid, tdc, adc & 0x3FF);
  DTRACE(DEB, "\t\t\ttotal time %f, chip time %f ns\n", totalTime, chipTime);

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

  // std::pair<uint8_t, uint8_t> dp = std::make_pair(det, plane);

  ClusterContainer cluster;
  double maxDeltaTime = 0;
  int clusterCount = 0;
  double time1 = 0, time2 = 0;
  uint32_t adc1 = 0;
  uint16_t strip1 = 0;
  uint16_t strip2 = 0;

  for (auto &itHits : m_hits[dp]) {
    time2 = time1;

    time1 = (double)std::get<0>(itHits);
    strip1 = std::get<1>(itHits);
    adc1 = std::get<2>(itHits);
    strip2 = std::get<3>(itHits);
    if (!cluster.empty()) {
      if (std::fabs(time1 - time2) > m_config.pDeltaTimeHits) {
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
    cluster.emplace_back(strip1, time1, adc1, strip2);
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
int Clusterer::ClusterByPad(std::pair<uint8_t, uint8_t> dp,
                            ClusterContainer &cluster, double maxDeltaTime) {

  int stripCount = 0;
  uint16_t spanCluster = 0;
  int adcTotal = 0;
  int clusterCount = 0;
  double largestTime = 0;
  double smallestTime = 0;
  int delta_x = -1;
  int delta_y = -1;
  double time1 = 0;
  double time2 = 0;
  double x_adc = 0;
  double y_adc = 0;
  double xy_time = 0;

  int adc1 = 0;
  int adc2 = 0;
  int idX1 = 0;
  int idX2 = 0;
  int idY1 = 0;
  int idY2 = 0;

  std::vector<double> vADC;
  std::vector<double> vIdx;
  std::vector<double> vIdy;
  std::vector<double> vTimes;
  auto det = std::get<0>(dp);
  auto plane = std::get<1>(dp);
  /*std::sort(begin(cluster), end(cluster),
            [](const ClusterTuple &t1, const ClusterTuple &t2) {
              return std::get<0>(t1) < std::get<0>(t2) ||
                     (std::get<0>(t1) == std::get<0>(t2) &&
                      std::get<3>(t1) < std::get<3>(t2)) ||
                     (std::get<0>(t1) == std::get<0>(t2) &&
                      std::get<3>(t1) == std::get<3>(t2) &&
                      std::get<1>(t1) > std::get<1>(t2));
            });
  */
  while (cluster.size() > 0) {
    ClusterContainer::const_iterator it = cluster.end();
    it--;
    idX1 = std::get<0>(*it);
    time1 = std::get<1>(*it);
    adc1 = std::get<2>(*it);
    adc1 = adc1 & 0x3FF;
    idY1 = std::get<3>(*it);
    // Start of new cluster
    ClusterContainer clusterFound;
    clusterFound.emplace_back(idX1, time1, adc1, idY1);
    vIdx.clear();
    vIdy.clear();
    vTimes.clear();
    vADC.clear();
    vIdx.emplace_back(idX1);
    vIdy.emplace_back(idY1);
    vTimes.emplace_back(time1);
    vADC.emplace_back(adc1);
    smallestTime = time1;
    largestTime = time1;
    adcTotal = adc1;
    x_adc = adc1 * idX1;
    y_adc = adc1 * idY1;
    xy_time = adc1 * time1;
    stripCount = 1;
    // Remove from list of pads
    cluster.pop_back();

    // Loop over pads in cluster
    for (int n = 0; n < clusterFound.size(); n++) {
      idX1 = std::get<0>(clusterFound[n]);
      time1 = std::get<1>(clusterFound[n]);
      adc1 = std::get<2>(clusterFound[n]);
      adc1 = adc1 & 0x3FF;
      idY1 = std::get<3>(clusterFound[n]);

      ClusterContainer::iterator itCluster = cluster.begin();
      // Loop over all pads to check whether they belong to cluster
      while (itCluster != cluster.end()) {
        idX2 = std::get<0>(*itCluster);
        time2 = std::get<1>(*itCluster);
        adc2 = std::get<2>(*itCluster);
        adc2 = adc2 & 0x3FF;
        idY2 = std::get<3>(*itCluster);
        // Pad in vector belongs to cluster
        if (std::fabs(idX1 - idX2) - 1 <= m_config.pMissingStripsClusterX &&
            std::fabs(idY1 - idY2) - 1 <= m_config.pMissingStripsClusterY &&
            time2 - smallestTime <= m_config.pSpanClusterTime &&
            largestTime - time2 <= m_config.pSpanClusterTime) {
          clusterFound.emplace_back(idX2, time2, adc2, idY2);
          stripCount++;
          if (time1 < smallestTime) {
            smallestTime = time1;
          }
          if (time1 > largestTime) {
            largestTime = time1;
          }
          if (std::fabs(idX1 - idX2) > delta_x) {
            delta_x = std::fabs(idX1 - idX2) - 1;
          }
          if (std::fabs(idY1 - idY2) > delta_y) {
            delta_y = std::fabs(idY1 - idY2) - 1;
          }
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
          // Pad in vector does not belong to cluster yet
          ++itCluster;
        }
      }
    }

    if (clusterFound.size() >= m_config.pMinClusterSize) {
      clusterCount++;
      spanCluster = (largestTime - smallestTime);
      if (m_config.pShowStats) {
        m_stats.SetStatsPlane("DeltaTimeHits", dp, maxDeltaTime);
        m_stats.SetStatsPlane("SpanClusterTime", dp, spanCluster);
        m_stats.SetStatsPlane("ClusterSize", dp, stripCount);
      }

      ClusterPlane clusterPlane;
      clusterPlane.id = m_cluster_detector_id;
      clusterPlane.det = det;
      clusterPlane.plane = plane;
      clusterPlane.size = stripCount;
      clusterPlane.adc = adcTotal;
      if (plane == 0) {
        clusterPlane.pos = x_adc / adcTotal;
      } else {
        clusterPlane.pos = y_adc / adcTotal;
      }
      clusterPlane.time = xy_time / adcTotal;
      clusterPlane.time_utpc = xy_time / adcTotal;
      clusterPlane.pos_utpc = clusterPlane.pos;
      clusterPlane.time_algo = xy_time / adcTotal;
      clusterPlane.pos_algo = clusterPlane.pos;
      clusterPlane.strips = std::move(vIdx);
      clusterPlane.times = std::move(vTimes);
      clusterPlane.adcs = std::move(vADC);
      clusterPlane.max_delta_time = maxDeltaTime;
      if (delta_x > delta_y) {
        if (delta_x == -1)
          delta_x = 0;
        clusterPlane.max_missing_strip = delta_x;
      } else {
        if (delta_y == -1)
          delta_y = 0;
        clusterPlane.max_missing_strip = delta_y;
      }
      clusterPlane.span_cluster = spanCluster;
      clusterPlane.plane_coincidence = false;
      m_clusters_new[dp].emplace_back(std::move(clusterPlane));

      DTRACE(DEB, "\ncluster det %d, id %d", (int)det, clusterPlane.id);
      DTRACE(DEB, "\tpos x/pos y: %f/%f", x_adc / adcTotal, y_adc / adcTotal);
      DTRACE(DEB, "\ttime: %f", clusterPlane.time);
      DTRACE(DEB, "\tadc: %u", clusterPlane.adc);
      DTRACE(DEB, "\tsize: %u\n", clusterPlane.size);
      if (m_config.pShowStats) {
        m_stats.SetStatsPlane("ClusterCntPlane", dp, 0);
      }
      m_cluster_detector_id++;
    }
  }
  return clusterCount;
}

//====================================================================================================================
int Clusterer::ClusterByStrip(std::pair<uint8_t, uint8_t> dp,
                              ClusterContainer &cluster, double maxDeltaTime) {
  int maxMissingStrip = 0;
  uint16_t spanCluster = 0;

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
  std::vector<double> vADC;
  std::vector<double> vStrips;
  std::vector<double> vTimes;
  auto det = std::get<0>(dp);
  auto plane = std::get<1>(dp);

  std::sort(begin(cluster), end(cluster),
            [](const ClusterTuple &t1, const ClusterTuple &t2) {
              return std::get<0>(t1) < std::get<0>(t2) ||
                     (std::get<0>(t1) == std::get<0>(t2) &&
                      std::get<1>(t1) > std::get<1>(t2));
            });
  for (auto &itCluster : cluster) {
    adc2 = adc1;
    strip2 = strip1;
    strip1 = std::get<0>(itCluster);
    time1 = std::get<1>(itCluster);
    adc1 = std::get<2>(itCluster);
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
      // position_utpc = (double)strip1;
      DTRACE(DEB, "\nDetector %d, plane %d cluster:\n", (int)det, (int)plane);
    }

    // Add members of a cluster, if it is either the beginning of a cluster,
    // or if strip gap and time span is correct
    if (stripCount == 0 ||
        (std::fabs(strip1 - strip2) > 0 &&
         std::fabs(strip1 - strip2) - 1 <= m_config.pMissingStripsClusterX &&
         time1 - startTime <= m_config.pSpanClusterTime &&
         largestTime - time1 <= m_config.pSpanClusterTime)) {
      DTRACE(DEB, "\tstrip %d, time %llu, adc %d:\n", strip1, (uint64_t)time1,
             adc1);
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
        // position_utpc = (double)strip1;
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
    // Stop clustering if gap between strips is too large or time span too long
    else if (std::fabs(strip1 - strip2) - 1 > m_config.pMissingStripsClusterX ||
             time1 - startTime > m_config.pSpanClusterTime ||
             largestTime - time1 > m_config.pSpanClusterTime) {
      // Valid cluster
      if (stripCount < m_config.pMinClusterSize || totalADC == 0) {
        DTRACE(DEB, "******** INVALID ********\n\n");
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
        if (m_config.pShowStats) {
          m_stats.SetStatsPlane("DeltaTimeHits", dp, maxDeltaTime);
          m_stats.SetStatsPlane("MissingStripsCluster", dp, maxMissingStrip);
          m_stats.SetStatsPlane("SpanClusterTime", dp, spanCluster);
          m_stats.SetStatsPlane("ClusterSize", dp, stripCount);
        }

        ClusterPlane clusterPlane;
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
        AlgorithmUTPC(idx_left, idx_right, vADC, vStrips, vTimes, pos_utpc,
                      time_utpc, pos_algo, time_algo);
        // COT only over Threshold
        if (m_config.pAlgo == 2) {
          pos_algo = centerOfGravity_ovTh;
          time_algo = centerOfTime_ovTh;
        }
        // COT2 only over Threshold
        else if (m_config.pAlgo == 3) {
          pos_algo = centerOfGravity2_ovTh;
          time_algo = centerOfTime2_ovTh;
        }
        // time of highest ADC
        else if (m_config.pAlgo == 4) {
          pos_algo = largestADCPos;
          time_algo = largestADCTime;
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

        DTRACE(DEB, "Cluster id %d\n", m_cluster_id);
        clusterPlane.id = static_cast<uint32_t>(m_cluster_id);
        clusterPlane.det = det;
        clusterPlane.plane = plane;
        if (plane == 0) {
          auto det_plane = std::make_pair(det, 1);
          if (m_config.GetDetectorPlane(det_plane) == false) {
            ClusterPlane clusterOtherPlane;
            clusterOtherPlane = clusterPlane;
            clusterOtherPlane.plane = 1;
            m_clusters_new[det_plane].emplace_back(
                std::move(clusterOtherPlane));
          }
        } else {
          auto det_plane = std::make_pair(det, 0);
          if (m_config.GetDetectorPlane(det_plane) == false) {
            ClusterPlane clusterOtherPlane;
            clusterOtherPlane = clusterPlane;
            clusterOtherPlane.plane = 0;
            m_clusters_new[det_plane].emplace_back(
                std::move(clusterOtherPlane));
          }
        }
        m_clusters_new[dp].emplace_back(std::move(clusterPlane));
        if (m_config.pShowStats) {
          m_stats.SetStatsPlane("ClusterCntPlane", dp, 0);
        }
        clusterCount++;
      }
      // Reset all parameters
      startTime = 0;
      largestTime = 0;
      largestADCTime = 0;
      largestADCPos = 0;
      stripCount = 0;
      centerOfGravity = 0;
      centerOfTime = 0;
      totalADC = 0;
      centerOfGravity2 = 0;
      centerOfTime2 = 0;
      totalADC2 = 0;
      centerOfGravity_ovTh = 0;
      centerOfTime_ovTh = 0;
      totalADC_ovTh = 0;
      centerOfGravity2_ovTh = 0;
      centerOfTime2_ovTh = 0;
      totalADC2_ovTh = 0;
      strip1 = 0;
      maxMissingStrip = 0;
      vADC.clear();
      vStrips.clear();
      vTimes.clear();
    }
  }

  // At the end of the clustering, check again if there is a last valid cluster
  if (stripCount < m_config.pMinClusterSize || totalADC == 0) {
    DTRACE(DEB, "******** INVALID ********\n\n");
  } else {
    spanCluster = (largestTime - startTime);
    centerOfGravity = (centerOfGravity / totalADC);
    centerOfTime = (centerOfTime / totalADC);
    centerOfGravity2 = (centerOfGravity2 / totalADC2);
    centerOfTime2 = (centerOfTime2 / totalADC2);
    if (m_config.pShowStats) {
      m_stats.SetStatsPlane("DeltaTimeHits", dp, maxDeltaTime);
      m_stats.SetStatsPlane("MissingStripsCluster", dp, maxMissingStrip);
      m_stats.SetStatsPlane("SpanClusterTime", dp, spanCluster);
      m_stats.SetStatsPlane("ClusterSize", dp, stripCount);
    }
    ClusterPlane clusterPlane;
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
    AlgorithmUTPC(idx_left, idx_right, vADC, vStrips, vTimes, pos_utpc,
                  time_utpc, pos_algo, time_algo);

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

    DTRACE(DEB, "Cluster id %d\n", m_cluster_id);
    clusterPlane.id = static_cast<uint32_t>(m_cluster_id);
    clusterPlane.det = det;
    clusterPlane.plane = plane;

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
    m_clusters_new[dp].emplace_back(std::move(clusterPlane));
    if (m_config.pShowStats) {
      m_stats.SetStatsPlane("ClusterCntPlane", dp, 0);
    }
    clusterCount++;
  }
  return clusterCount;
}

void Clusterer::AlgorithmUTPC(int idx_min_largest_time,
                              int idx_max_largest_time,
                              std::vector<double> &vADC,
                              std::vector<double> &vStrips,
                              std::vector<double> &vTimes, double &positionUTPC,
                              double &timeUTPC, double &positionAlgo,
                              double &timeAlgo) {
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
  if (m_config.pAlgo == 1) {
    positionAlgo = (p1 * a1 * a1 + p2 * a2 * a2 + p3 * a3 * a3) /
                   (a1 * a1 + a2 * a2 + a3 * a3);
    timeAlgo = (t1 * a1 * a1 + t2 * a2 * a2 + t3 * a3 * a3) /
               (a1 * a1 + a2 * a2 + a3 * a3);
  } else if (m_config.pAlgo == 0) {
    positionAlgo = (p1 * a1 + p2 * a2 + p3 * a3) / (a1 + a2 + a3);
    timeAlgo = (t1 * a1 + t2 * a2 + t3 * a3) / (a1 + a2 + a3);
  }
}

//====================================================================================================================
int Clusterer::MatchClustersDetector(uint8_t det) {
  int clusterCount = 0;
  auto dp0 = std::make_pair(det, 0);
  auto dp1 = std::make_pair(det, 1);
  ClusterVectorPlane::iterator itStartPlane1 = begin(m_clusters[dp1]);
  for (auto &c0 : m_clusters[dp0]) {
    double minDelta = 99999999.0;
    double lastDelta_t = 99999999.0;
    double delta_t = 99999999.0;
    bool isFirstMatch = true;
    ClusterVectorPlane::iterator bestMatchPlane1 = end(m_clusters[dp1]);

    for (ClusterVectorPlane::iterator c1 = itStartPlane1;
         c1 != end(m_clusters[dp1]); ++c1) {
      if ((*c1).plane_coincidence == false) {

        double chargeRatio = (double)(c0).adc / (double)(*c1).adc;
        lastDelta_t = delta_t;
        delta_t = (*c1).time - c0.time;
        if (m_config.pConditionCoincidence == "utpc") {
          delta_t = (*c1).time_utpc - c0.time_utpc;
        } else if (m_config.pConditionCoincidence == "charge2") {
          delta_t = (*c1).time_charge2 - c0.time_charge2;
        }
        if (chargeRatio >= m_config.pChargeRatioLower &&
            chargeRatio <= m_config.pChargeRatioUpper &&
            std::fabs(delta_t) < minDelta &&
            std::fabs(delta_t) <= m_config.pDeltaTimePlanes &&
            (c0.size + (*c1).size >= m_config.pCoincidentClusterSize)) {
          minDelta = std::fabs(delta_t);
          bestMatchPlane1 = c1;
          if (isFirstMatch) {
            itStartPlane1 = c1;
            isFirstMatch = false;
          }
        }
        if (std::fabs(delta_t) > std::fabs(lastDelta_t)) {
          break;
        }
      }
    }

    if (bestMatchPlane1 != end(m_clusters[dp1])) {
      c0.plane_coincidence = true;
      (*bestMatchPlane1).plane_coincidence = true;
      ClusterDetector clusterDetector;
      m_cluster_detector_id++;
      clusterDetector.id = m_cluster_detector_id;
      clusterDetector.det = det;
      clusterDetector.id0 = c0.id;
      clusterDetector.id1 = (*bestMatchPlane1).id;
      clusterDetector.size0 = c0.size;
      clusterDetector.size1 = (*bestMatchPlane1).size;
      clusterDetector.adc0 = c0.adc;
      clusterDetector.adc1 = (*bestMatchPlane1).adc;

      if (m_config.pTransform.size() == m_config.pDets.size()) {
        auto tx = m_config.pTransformX[m_config.pDets[det]];
        auto ty = m_config.pTransformY[m_config.pDets[det]];
        auto tz = m_config.pTransformZ[m_config.pDets[det]];

        clusterDetector.pos0 = c0.pos * std::get<0>(tx) +
                               (*bestMatchPlane1).pos * std::get<1>(tx) +
                               std::get<3>(tx);
        clusterDetector.pos1 = c0.pos * std::get<0>(ty) +
                               (*bestMatchPlane1).pos * std::get<1>(ty) +
                               std::get<3>(ty);
        clusterDetector.pos2 = c0.pos * std::get<0>(tz) +
                               (*bestMatchPlane1).pos * std::get<1>(tz) +
                               std::get<3>(tz);

        clusterDetector.pos0_utpc =
            c0.pos_utpc * std::get<0>(tx) +
            (*bestMatchPlane1).pos_utpc * std::get<1>(tx) + std::get<3>(tx);
        clusterDetector.pos1_utpc =
            c0.pos_utpc * std::get<0>(ty) +
            (*bestMatchPlane1).pos_utpc * std::get<1>(ty) + std::get<3>(ty);
        clusterDetector.pos2_utpc =
            c0.pos_utpc * std::get<0>(tz) +
            (*bestMatchPlane1).pos_utpc * std::get<1>(tz) + std::get<3>(tz);

        clusterDetector.pos0_charge2 =
            c0.pos_charge2 * std::get<0>(tx) +
            (*bestMatchPlane1).pos_charge2 * std::get<1>(tx) + std::get<3>(tx);
        clusterDetector.pos1_charge2 =
            c0.pos_charge2 * std::get<0>(ty) +
            (*bestMatchPlane1).pos_charge2 * std::get<1>(ty) + std::get<3>(ty);
        clusterDetector.pos2_charge2 =
            c0.pos_charge2 * std::get<0>(tz) +
            (*bestMatchPlane1).pos_charge2 * std::get<1>(tz) + std::get<3>(tz);

        clusterDetector.pos0_algo =
            c0.pos_algo * std::get<0>(tx) +
            (*bestMatchPlane1).pos_algo * std::get<1>(tx) + std::get<3>(tx);
        clusterDetector.pos1_algo =
            c0.pos_algo * std::get<0>(ty) +
            (*bestMatchPlane1).pos_algo * std::get<1>(ty) + std::get<3>(ty);
        clusterDetector.pos2_algo =
            c0.pos_algo * std::get<0>(tz) +
            (*bestMatchPlane1).pos_algo * std::get<1>(tz) + std::get<3>(tz);

      } else {
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
      }

      clusterDetector.time0 = c0.time;
      clusterDetector.time1 = (*bestMatchPlane1).time;
      clusterDetector.time0_utpc = c0.time_utpc;
      clusterDetector.time1_utpc = (*bestMatchPlane1).time_utpc;
      clusterDetector.time0_charge2 = c0.time_charge2;
      clusterDetector.time1_charge2 = (*bestMatchPlane1).time_charge2;
      clusterDetector.time0_algo = c0.time_algo;
      clusterDetector.time1_algo = (*bestMatchPlane1).time_algo;
      clusterDetector.dt0 = clusterDetector.time0 - last_time0;
      clusterDetector.dt1 = clusterDetector.time1 - last_time1;

      // std::cout << clusterDetector.dt1 << " " << clusterDetector.time1 << " "
      // << last_time1 << std::endl;
      /*
      clusterDetector.dt0_utpc = clusterDetector.time0_utpc - last_time0_utpc;
      clusterDetector.dt1_utpc = clusterDetector.time1_utpc - last_time1_utpc;
      clusterDetector.dt0_charge2 = clusterDetector.time0_charge2 -
      last_time0_charge2; clusterDetector.dt1_charge2 =
      clusterDetector.time1_charge2 - last_time1_charge2;
      */
      last_time0 = clusterDetector.time0;
      last_time1 = clusterDetector.time1;
      /*
      last_time0_utpc = clusterDetector.time0_utpc;
      last_time1_utpc = clusterDetector.time1_utpc;
      last_time0_charge2 = clusterDetector.time0_charge2;
      last_time1_charge2 = clusterDetector.time1_charge2;
      */
      clusterDetector.delta_plane =
          clusterDetector.time1 - clusterDetector.time0;

      if (m_config.pConditionCoincidence == "utpc") {
        clusterDetector.delta_plane =
            clusterDetector.time1_utpc - clusterDetector.time0_utpc;
      } else if (m_config.pConditionCoincidence == "charge2") {
        clusterDetector.delta_plane =
            clusterDetector.time1_charge2 - clusterDetector.time0_charge2;
      }

      if (m_config.pShowStats) {
        m_stats.SetStatsDetector("DeltaTimePlanes", det,
                                 std::fabs(clusterDetector.delta_plane));
        double ratio =
            100 * (double)clusterDetector.adc0 / (double)clusterDetector.adc1;
        if (ratio > 100.0) {
          ratio =
              100 * (double)clusterDetector.adc1 / (double)clusterDetector.adc0;
          m_stats.SetStatsDetector("ChargeRatio_1_0", det, ratio);
        } else {
          m_stats.SetStatsDetector("ChargeRatio_0_1", det, ratio);
        }
        m_stats.SetStatsDetector("ClusterCntDetector", det, 0);
        clusterCount++;
      }

      clusterDetector.max_delta_time0 = c0.max_delta_time;
      clusterDetector.max_delta_time1 = (*bestMatchPlane1).max_delta_time;
      clusterDetector.max_missing_strip0 = c0.max_missing_strip;
      clusterDetector.max_missing_strip1 = (*bestMatchPlane1).max_missing_strip;
      clusterDetector.span_cluster0 = c0.span_cluster;
      clusterDetector.span_cluster1 = (*bestMatchPlane1).span_cluster;
      clusterDetector.strips0 = c0.strips;
      clusterDetector.times0 = c0.times;
      clusterDetector.adcs0 = c0.adcs;
      clusterDetector.strips1 = (*bestMatchPlane1).strips;
      clusterDetector.times1 = (*bestMatchPlane1).times;
      clusterDetector.adcs1 = (*bestMatchPlane1).adcs;
      DTRACE(DEB, "\ncommon cluster det %d x/y: %d/%d", (int)det,
             clusterDetector.id0, clusterDetector.id1);
      DTRACE(DEB, "\tpos x/pos y: %f/%f", clusterDetector.pos0,
             clusterDetector.pos1);
      DTRACE(DEB, "\ttime x/time y: %llu/%llu", (uint64_t)clusterDetector.time0,
             (uint64_t)clusterDetector.time1);
      DTRACE(DEB, "\tadc x/adc y: %u/%u", clusterDetector.adc0,
             clusterDetector.adc1);
      DTRACE(DEB, "\tsize x/size y: %u/%u", clusterDetector.size0,
             clusterDetector.size1);
      DTRACE(DEB, "\tdelta time planes: %d\n",
             (int)clusterDetector.delta_plane);
      m_clusters_detector[det].emplace_back(std::move(clusterDetector));
    }
  }

  return clusterCount;
}

// void Clusterer::AnalyzeClustersPlane(uint8_t det, uint8_t plane,
// HitContainer& hits, HitContainer& newHits, double timeReadyToCluster,
// uint16_t correctionTime)
void Clusterer::AnalyzeClustersPlane(std::pair<uint8_t, uint8_t> dp) {

  if (ChooseHitsToBeClustered(dp) == false && m_hits[dp].empty()) {
    return;
  }
  int cnt = ClusterByTime(dp);

  DTRACE(DEB, "%d cluster in detector %d plane %d\n", cnt, (int)std::get<0>(dp),
         (int)std::get<1>(dp));

  m_hits[dp].clear();
}

void Clusterer::AnalyzeClustersDetector(uint8_t det) {
  int cnt = 0;
  auto dp0 = std::make_pair(det, 0);
  auto dp1 = std::make_pair(det, 1);
  if (ChooseClustersToBeMatched(dp0) == false && m_clusters[dp0].empty()) {
    return;
  }
  if (ChooseClustersToBeMatched(dp1) == false && m_clusters[dp1].empty()) {
    return;
  }
  cnt = MatchClustersDetector(det);

  if (m_config.pSaveWhat == 10 || m_config.pSaveWhat == 11 ||
      m_config.pSaveWhat == 110 || m_config.pSaveWhat == 111) {
    m_rootFile->SaveClustersPlane(std::move(m_clusters[dp0]));
    m_rootFile->SaveClustersPlane(std::move(m_clusters[dp1]));
  }

  if (m_config.pSaveWhat >= 100) {
    m_rootFile->SaveClustersDetector(std::move(m_clusters_detector[det]));
  }
  m_clusters[std::make_pair(det, 0)].clear();
  m_clusters[std::make_pair(det, 1)].clear();
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

  auto theMin = std::min_element(m_hits_new[dp].begin(), m_hits_new[dp].end(),
                                 [](const HitTuple &t1, const HitTuple &t2) {
                                   return std::get<0>(t1) < std::get<0>(t2);
                                 });

  // Nothing to cluster, tuples in newHits vector too recent
  if (std::get<0>(*theMin) > timeReadyToCluster) {

    //(smallest timestamp larger than
    // m_stats.GetLowestCommonTriggerTimestampPlane(dp)) Will be clustered
    // later
    return false;
  }

  // Sort vector newHits
  std::sort(begin(m_hits_new[dp]), end(m_hits_new[dp]),
            [](const HitTuple &t1, const HitTuple &t2) {
              return std::get<0>(t1) < std::get<0>(t2);
            });

  // First tuple with timestamp larger than
  // m_stats.GetLowestCommonTriggerTimestampPlane(dp)
  auto it = std::upper_bound(
      m_hits_new[dp].begin(), m_hits_new[dp].end(),
      std::make_tuple(m_stats.GetLowestCommonTriggerTimestampPlane(dp), 0, 0,
                      0),
      [](const HitTuple &t1, const HitTuple &t2) {
        return std::get<0>(t1) < std::get<0>(t2);
      });

  // Find elements in vector that could still be part of a cluster,
  // since they are close in time to
  // m_stats.GetLowestCommonTriggerTimestampPlane(dp)
  while (it != m_hits_new[dp].end()) {
    if (std::get<0>(*it) - timeReadyToCluster > m_config.pDeltaTimeHits) {
      break;
    }
    timeReadyToCluster = std::get<0>(*it);
    ++it;
  }

  // std::cout << "still in cluster " << std::get < 0 > (*it) << "\n" <<
  // std::endl;

  int index = std::distance(m_hits_new[dp].begin(), it);
  // Insert the data that is ready to be clustered from newHits into hits
  m_hits[dp].insert(m_hits[dp].end(),
                    std::make_move_iterator(m_hits_new[dp].begin()),
                    std::make_move_iterator(m_hits_new[dp].begin() + index));
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
    auto theMin =
        std::min_element(m_clusters_new[dp].begin(), m_clusters_new[dp].end(),
                         [](const ClusterPlane &t1, const ClusterPlane &t2) {
                           return t1.time_utpc < t2.time_utpc;
                         });

    // Nothing to cluster, clusters in newClusters vector too recent
    if ((*theMin).time_utpc > timeReadyToMatch) {

      //(smallest time larger than timeReadyToMatch)
      // Will be matched later
      return false;
    }

    // Sort vector newClusters based on time
    std::sort(begin(m_clusters_new[dp]), end(m_clusters_new[dp]),
              [](const ClusterPlane &t1, const ClusterPlane &t2) {
                return t1.time_utpc < t2.time_utpc;
              });

    ClusterPlane theCluster;
    theCluster.time_utpc = timeReadyToMatch;

    // First ClusterPlane with time bigger than timeReadyToMatch
    auto it = std::upper_bound(
        m_clusters_new[dp].begin(), m_clusters_new[dp].end(), theCluster,
        [](const ClusterPlane &t1, const ClusterPlane &t2) {
          return t1.time_utpc < t2.time_utpc;
        });

    // Find elements in vector that could still be matched with another
    // cluster since they are close in time to timeReadyToMatch
    while (it != m_clusters_new[dp].end()) {
      if ((*it).time_utpc - timeReadyToMatch > m_config.pDeltaTimeHits) {
        break;
      }
      timeReadyToMatch = (*it).time_utpc;
      ++it;
    }
    index = std::distance(m_clusters_new[dp].begin(), it);
  } else if (m_config.pConditionCoincidence == "charge2") {
    auto theMin =
        std::min_element(m_clusters_new[dp].begin(), m_clusters_new[dp].end(),
                         [](const ClusterPlane &t1, const ClusterPlane &t2) {
                           return t1.time_charge2 < t2.time_charge2;
                         });

    // Nothing to cluster, clusters in newClusters vector too recent
    if ((*theMin).time_charge2 > timeReadyToMatch) {

      //(smallest time larger than timeReadyToMatch)
      // Will be matched later
      return false;
    }

    // Sort vector newClusters based on time
    std::sort(begin(m_clusters_new[dp]), end(m_clusters_new[dp]),
              [](const ClusterPlane &t1, const ClusterPlane &t2) {
                return t1.time_charge2 < t2.time_charge2;
              });

    ClusterPlane theCluster;
    theCluster.time_charge2 = timeReadyToMatch;

    // First ClusterPlane with time that bigger than timeReadyToMatch
    auto it = std::upper_bound(
        m_clusters_new[dp].begin(), m_clusters_new[dp].end(), theCluster,
        [](const ClusterPlane &t1, const ClusterPlane &t2) {
          return t1.time_charge2 < t2.time_charge2;
        });

    // Find elements in vector that could still be matched with another
    // cluster since they are close in time to timeReadyToMatch
    while (it != m_clusters_new[dp].end()) {
      if ((*it).time_charge2 - timeReadyToMatch > m_config.pDeltaTimeHits) {
        break;
      }
      timeReadyToMatch = (*it).time_charge2;
      ++it;
    }
    index = std::distance(m_clusters_new[dp].begin(), it);
  } else {
    auto theMin =
        std::min_element(m_clusters_new[dp].begin(), m_clusters_new[dp].end(),
                         [](const ClusterPlane &t1, const ClusterPlane &t2) {
                           return t1.time < t2.time;
                         });

    // Nothing to cluster, clusters in newClusters vector too recent
    if ((*theMin).time > timeReadyToMatch) {
      //(smallest time larger than timeReadyToMatch)
      // Will be matched later
      return false;
    }

    // Sort vector newClusters based on time
    std::sort(begin(m_clusters_new[dp]), end(m_clusters_new[dp]),
              [](const ClusterPlane &t1, const ClusterPlane &t2) {
                return t1.time < t2.time;
              });

    ClusterPlane theCluster;
    theCluster.time = timeReadyToMatch;

    // First ClusterPlane with time that bigger than timeReadyToMatch
    auto it = std::upper_bound(
        m_clusters_new[dp].begin(), m_clusters_new[dp].end(), theCluster,
        [](const ClusterPlane &t1, const ClusterPlane &t2) {
          return t1.time < t2.time;
        });

    // Find elements in vector that could still be matched with another
    // cluster since they are close in time to timeReadyToMatch
    while (it != m_clusters_new[dp].end()) {
      if ((*it).time - timeReadyToMatch > m_config.pDeltaTimeHits) {
        break;
      }
      timeReadyToMatch = (*it).time;
      ++it;
    }
    index = std::distance(m_clusters_new[dp].begin(), it);
  }
  // Insert the clusters that are ready to be matched from newClusters into
  // clusters
  m_clusters[dp].insert(
      m_clusters[dp].end(), std::make_move_iterator(m_clusters_new[dp].begin()),
      std::make_move_iterator(m_clusters_new[dp].begin() + index));
  // Delete the clusters from newClusters
  m_clusters_new[dp].erase(m_clusters_new[dp].begin(),
                           m_clusters_new[dp].begin() + index);

  return true;
}

void Clusterer::FinishAnalysis() {
  double ts = 0;
  for (auto const &fec : m_config.pFecs) {
    if (ts < m_stats.GetMaxTriggerTimestamp(fec)) {
      ts = m_stats.GetMaxTriggerTimestamp(fec);
    }
  }
  for (auto const &det : m_config.pDets) {
    auto dp0 = std::make_pair(det.first, 0);
    auto dp1 = std::make_pair(det.first, 1);

    // Set the largest timestamp of plane to detector
    // cluster all remaining data in plane
    m_stats.SetLowestCommonTriggerTimestampPlane(dp0, ts);
    m_stats.SetLowestCommonTriggerTimestampPlane(dp1, ts);
    m_stats.SetLowestCommonTriggerTimestampDet(det.first, ts);

    AnalyzeClustersPlane(dp0);
    AnalyzeClustersPlane(dp1);
    AnalyzeClustersDetector(det.first);
    if (m_config.pSaveWhat % 2 == 1) {
      m_rootFile->SaveHits();
    }
  }
  if (m_config.pShowStats) {
    if (m_config.pSaveWhat >= 10) {
      m_stats.PrintClusterStats(m_config);
    }
    m_stats.PrintFECStats(m_config);
  }
}
