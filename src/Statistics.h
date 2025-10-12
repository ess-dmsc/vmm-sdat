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
**  Statistics.h
**
****************************************************************************/
#pragma once

#include "Configuration.h"

class Statistics {
public:
  Statistics();
  ~Statistics() = default;

  void CreateClusterStats(Configuration &config);
  void CreateFECStats(Configuration &config);
  void CreatePCAPStats(Configuration &config);

  long GetStatsDetector(std::string stats, uint8_t det, int n);
  void SetStatsDetector(std::string stats, uint8_t det, double value);
  long GetStatsPlane(std::string stats, std::pair<uint8_t, uint8_t> dp, int n);
  void SetStatsPlane(std::string stats, std::pair<uint8_t, uint8_t> dp,
                     double value);

  void IncrementCounter(std::string error, uint16_t fecId,
                        uint64_t increment = 1);
  long GetCounter(std::string error, uint16_t fecId);

  double GetOldTriggerTimestamp(uint16_t fecId);
  void SetOldTriggerTimestamp(uint16_t fecId, double srsTimestamp);
  double GetFirstTriggerTimestamp(uint16_t fecId);
  void SetFirstTriggerTimestamp(uint16_t fecId, double srsTimestamp);
  double GetMaxTriggerTimestamp(uint16_t fecId);
  void SetMaxTriggerTimestamp(uint16_t fecId, double srsTimestamp);
  uint64_t GetLastFrameCounter(uint16_t fecId);
  void SetLastFrameCounter(uint16_t fecId, uint64_t frameCounter);

  double GetLowestCommonTriggerTimestampDet(uint8_t det);
  void SetLowestCommonTriggerTimestampDet(uint8_t det, double val);
  double GetLowestCommonTriggerTimestampPlane(std::pair<uint8_t, uint8_t> dp);
  void SetLowestCommonTriggerTimestampPlane(std::pair<uint8_t, uint8_t> dp,
                                            double val);

  void PrintClusterStats(Configuration &config);
  void PrintFECStats(Configuration &config);

  void StatsOutput(int n, long val, std::string stat, long cnt, long cnt0 = 0,
                   long cnt1 = 0, long cnt2 = 0);

private:
  double m_acq_time;
  std::map<std::pair<std::pair<uint8_t, uint8_t>, std::string>,
           std::vector<long>>
      m_stats_plane;
  std::map<std::pair<uint8_t, std::string>, std::vector<long>> m_stats_detector;
  std::vector<std::string> m_stats_plane_names;
  std::vector<std::string> m_stats_detector_names;
  std::map<std::string, double> m_factors;
  std::map<std::string, double> m_limits;
  std::map<std::string, std::string> m_units;
  std::map<std::pair<uint8_t, std::string>, long> m_counters;
  std::vector<std::string> m_counter_names;

  // per plane
  std::map<std::pair<uint8_t, uint8_t>, double>
      m_lowestCommonTriggerTimestamp_plane;
  // per detector
  std::map<uint8_t, double> m_lowestCommonTriggerTimestamp_det;
  // per FEC
  std::map<uint8_t, double> m_deltaTriggerTimestamp;
  std::map<uint8_t, double> m_oldTriggerTimestamp;
  std::map<uint8_t, double> m_maxTriggerTimestamp;
  std::map<uint8_t, double> m_firstTriggerTimestamp;
  std::map<uint8_t, double> m_lastTriggerTimestamp;
  std::map<uint8_t, double> m_lastFrameCounter;
};

