#include <iostream>
//#include <sstream>
//#include <cstring>

#include "Statistics.h"



void Statistics::CreatePCAPStats(Configuration &config) {
  for (auto const &fec : config.pFecs) {
    m_counter_names.push_back("parser_frame_seq_errors");
    m_counters.emplace(
        std::make_pair(std::make_pair(fec, "parser_frame_seq_errors"), 0));

    m_counter_names.push_back("parser_frame_missing_errors");
    m_counters.emplace(
        std::make_pair(std::make_pair(fec, "parser_frame_missing_errors"), 0));

    m_counter_names.push_back("parser_framecounter_overflows");
    m_counters.emplace(
        std::make_pair(std::make_pair(fec, "parser_framecounter_overflows"), 0));

    m_counter_names.push_back("parser_timestamp_seq_errors");
    m_counters.emplace(
        std::make_pair(std::make_pair(fec, "parser_timestamp_seq_errors"), 0));

    m_counter_names.push_back("parser_timestamp_overflows");
    m_counters.emplace(
        std::make_pair(std::make_pair(fec, "parser_timestamp_overflows"), 0));

    m_counter_names.push_back("parser_bad_frames");
    m_counters.emplace(std::make_pair(std::make_pair(fec, "parser_bad_frames"), 0));

    m_counter_names.push_back("parser_good_frames");
    m_counters.emplace(std::make_pair(std::make_pair(fec, "parser_good_frames"), 0));

    m_counter_names.push_back("parser_readouts");
    m_counters.emplace(std::make_pair(std::make_pair(fec, "parser_readouts"), 0));

    m_counter_names.push_back("parser_markers");
    m_counters.emplace(std::make_pair(std::make_pair(fec, "parser_markers"), 0));

    m_counter_names.push_back("parser_data");
    m_counters.emplace(std::make_pair(std::make_pair(fec, "parser_data"), 0));
  }
}

void Statistics::CreateFECStats(Configuration &config) {
  for (auto const &fec : config.pFecs) {
    m_deltaTriggerTimestamp.emplace(std::make_pair(fec, 0));
    m_oldTriggerTimestamp.emplace(std::make_pair(fec, 0));
    m_counter_names.push_back("time_stamp_too_large");
    m_counters.emplace(
        std::make_pair(std::make_pair(fec, "time_stamp_too_large"), 0));
    m_counter_names.push_back("time_stamp_order_error");
    m_counters.emplace(
        std::make_pair(std::make_pair(fec, "time_stamp_order_error"), 0));
    m_counter_names.push_back("time_stamp_overflow");
    m_counters.emplace(
        std::make_pair(std::make_pair(fec, "time_stamp_overflow"), 0));
    m_counter_names.push_back("trigger_period_error");
    m_counters.emplace(
        std::make_pair(std::make_pair(fec, "trigger_period_error"), 0));

    if(!config.pIsPcap) {
      m_counter_names.push_back("parser_data");
      m_counters.emplace(std::make_pair(std::make_pair(fec, "parser_data"), 0));
    }  

  }
}

void Statistics::CreateClusterStats(Configuration &config) {
  int size = 0;
  for (auto const &det : config.pDets) {
    auto plane0 = std::make_pair(det.first, 0);
    auto plane1 = std::make_pair(det.first, 1);

    // initialize timestamps
    m_lowestCommonTriggerTimestamp_det[det.first] = 0;
    m_lowestCommonTriggerTimestamp_plane[plane0] = 0;
    m_lowestCommonTriggerTimestamp_plane[plane1] = 0;

    m_units.emplace(std::make_pair("delta_time_hits", "ns"));
    m_factors.emplace(std::make_pair("delta_time_hits", 0.02));
    m_limits.emplace(std::make_pair("delta_time_hits",
                                    1 + config.pDeltaTimeHits *
                                            m_factors["delta_time_hits"]));
    size = static_cast<int>(m_limits["delta_time_hits"]);
    std::vector<long> v(size, 0);
    m_stats_plane_names.push_back("delta_time_hits");
    m_stats_plane.emplace(
        std::make_pair(std::make_pair(plane0, "delta_time_hits"), v));
    m_stats_plane.emplace(
        std::make_pair(std::make_pair(plane1, "delta_time_hits"), v));

    m_units.emplace(std::make_pair("missing_strips_cluster", "strips"));
    m_factors.emplace(std::make_pair("missing_strips_cluster", 1));
    m_limits.emplace(std::make_pair(
        "missing_strips_cluster", 1 + config.pMissingStripsCluster *
                                          m_factors["missing_strips_cluster"]));
    size = static_cast<int>(m_limits["missing_strips_cluster"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_plane_names.push_back("missing_strips_cluster");
    m_stats_plane.emplace(
        std::make_pair(std::make_pair(plane0, "missing_strips_cluster"), v));
    m_stats_plane.emplace(
        std::make_pair(std::make_pair(plane1, "missing_strips_cluster"), v));

    m_units.emplace(std::make_pair("span_cluster_time", "ns"));
    m_factors.emplace(std::make_pair("span_cluster_time", 0.02));
    m_limits.emplace(std::make_pair("span_cluster_time",
                                    1 + config.pSpanClusterTime *
                                            m_factors["span_cluster_time"]));
    size = static_cast<int>(m_limits["span_cluster_time"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_plane_names.push_back("span_cluster_time");
    m_stats_plane.emplace(
        std::make_pair(std::make_pair(plane0, "span_cluster_time"), v));
    m_stats_plane.emplace(
        std::make_pair(std::make_pair(plane1, "span_cluster_time"), v));

    m_units.emplace(std::make_pair("cluster_size", "strips"));
    m_factors.emplace(std::make_pair("cluster_size", 1));
    m_limits.emplace(
        std::make_pair("cluster_size", 1 + 64 * m_factors["cluster_size"]));
    size = static_cast<int>(m_limits["cluster_size"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_plane_names.push_back("cluster_size");
    m_stats_plane.emplace(
        std::make_pair(std::make_pair(plane0, "cluster_size"), v));
    m_stats_plane.emplace(
        std::make_pair(std::make_pair(plane1, "cluster_size"), v));

    m_units.emplace(std::make_pair("cluster_cnt_plane", ""));
    m_factors.emplace(std::make_pair("cluster_cnt_plane", 1));
    m_limits.emplace(std::make_pair("cluster_cnt_plane", 1));
    size = static_cast<int>(m_limits["cluster_cnt_plane"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_plane_names.push_back("cluster_cnt_plane");
    m_stats_plane.emplace(
        std::make_pair(std::make_pair(plane0, "cluster_cnt_plane"), v));
    m_stats_plane.emplace(
        std::make_pair(std::make_pair(plane1, "cluster_cnt_plane"), v));

    m_units.emplace(std::make_pair("delta_time_planes", "ns"));
    m_factors.emplace(std::make_pair("delta_time_planes", 0.02));
    m_limits.emplace(std::make_pair("delta_time_planes",
                                    1 + config.pDeltaTimePlanes *
                                            m_factors["delta_time_planes"]));
    size = static_cast<int>(m_limits["delta_time_planes"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector_names.push_back("delta_time_planes");
    m_stats_detector.emplace(
        std::make_pair(std::make_pair(det.first, "delta_time_planes"), v));

    m_units.emplace(std::make_pair("charge_ratio_0_1", "%"));
    m_factors.emplace(std::make_pair("charge_ratio_0_1", 0.1));
    m_limits.emplace(std::make_pair("charge_ratio_0_1", 11));
    size = static_cast<int>(m_limits["charge_ratio_0_1"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector_names.push_back("charge_ratio_0_1");
    m_stats_detector.emplace(
        std::make_pair(std::make_pair(det.first, "charge_ratio_0_1"), v));

    m_units.emplace(std::make_pair("charge_ratio_1_0", "%"));
    m_factors.emplace(std::make_pair("charge_ratio_1_0", 0.1));
    m_limits.emplace(std::make_pair("charge_ratio_1_0", 10));
    size = static_cast<int>(m_limits["charge_ratio_1_0"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector_names.push_back("charge_ratio_1_0");
    m_stats_detector.emplace(
        std::make_pair(std::make_pair(det.first, "charge_ratio_1_0"), v));

    m_units.emplace(std::make_pair("cluster_cnt_detector", ""));
    m_factors.emplace(std::make_pair("cluster_cnt_detector", 1));
    m_limits.emplace(std::make_pair("cluster_cnt_detector", 1));
    size = static_cast<int>(m_limits["cluster_cnt_detector"]);
    v.resize(size);
    std::fill(v.begin(), v.end(), 0);
    m_stats_detector_names.push_back("cluster_cnt_detector");
    m_stats_detector.emplace(
        std::make_pair(std::make_pair(det.first, "cluster_cnt_detector"), v));
  }
}

long Statistics::GetStatsDetector(std::string stats, uint8_t det, int n) {
  if (n < m_limits[stats.c_str()]) {
    return m_stats_detector[std::make_pair(det, stats.c_str())][n];
  }
  return -1;
}

void Statistics::SetStatsDetector(std::string stats, uint8_t det,
                                  double value) {
  if (value * m_factors[stats] < m_limits[stats]) {
    m_stats_detector[std::make_pair(det, stats)]
                    [static_cast<unsigned int>(value * m_factors[stats])]++;
  } else {
    m_stats_detector[std::make_pair(det, stats)][m_limits[stats] - 1]++;
  }
}

long Statistics::GetStatsPlane(std::string stats,
                               std::pair<uint8_t, uint8_t> dp, int n) {
  if (n < m_limits[stats]) {
    return m_stats_plane[std::make_pair(dp, stats)][n];
  }
  return -1;
}

void Statistics::SetStatsPlane(std::string stats,
                               std::pair<uint8_t, uint8_t> dp, double value) {
  if (value * m_factors[stats] < m_limits[stats]) {
    m_stats_plane[std::make_pair(dp, stats)]
                 [static_cast<unsigned int>(value * m_factors[stats])]++;
  } else {
    m_stats_plane[std::make_pair(dp, stats)][m_limits[stats] - 1]++;
  }
}

void Statistics::IncrementCounter(std::string error, uint8_t fecId, uint64_t increment) {
  m_counters[std::make_pair(fecId, error)] += increment;
}

long Statistics::GetCounter(std::string error, uint8_t fecId) {
  return m_counters[std::make_pair(fecId, error)];
}

double Statistics::GetDeltaTriggerTimestamp(uint8_t fecId) {
  return m_deltaTriggerTimestamp[fecId];
}

void Statistics::SetDeltaTriggerTimestamp(uint8_t fecId, double val) {
  m_deltaTriggerTimestamp[fecId] = val;
}

double Statistics::GetOldTriggerTimestamp(uint8_t fecId) {
  return m_oldTriggerTimestamp[fecId];
}

void Statistics::SetOldTriggerTimestamp(uint8_t fecId, double srsTimestamp) {
  m_oldTriggerTimestamp[fecId] = srsTimestamp;
}

double Statistics::GetFirstTriggerTimestamp(uint8_t fecId) {
  return m_firstTriggerTimestamp[fecId];
}

void Statistics::SetFirstTriggerTimestamp(uint8_t fecId, double srsTimestamp) {
  m_firstTriggerTimestamp[fecId] = srsTimestamp;
}

double Statistics::GetMaxTriggerTimestamp(uint8_t fecId) {
  return m_maxTriggerTimestamp[fecId];
}

void Statistics::SetMaxTriggerTimestamp(uint8_t fecId, double srsTimestamp) {
  m_maxTriggerTimestamp[fecId] = srsTimestamp;
}

double Statistics::GetLowestCommonTriggerTimestampDet(uint8_t det) {
  return m_lowestCommonTriggerTimestamp_det[det];
}

void Statistics::SetLowestCommonTriggerTimestampDet(uint8_t det, double val) {
  m_lowestCommonTriggerTimestamp_det[det] = val;
}

double Statistics::GetLowestCommonTriggerTimestampPlane(
    std::pair<uint8_t, uint8_t> dp) {
  return m_lowestCommonTriggerTimestamp_plane[dp];
}

void Statistics::SetLowestCommonTriggerTimestampPlane(
    std::pair<uint8_t, uint8_t> dp, double val) {
  m_lowestCommonTriggerTimestamp_plane[dp] = val;
}

void Statistics::PrintClusterStats(Configuration &config) {
  long totalPlane0 = 0;
  long totalPlane1 = 0;
  long totalDetector = 0;
  bool bothPlanes = false;

  for (auto const &det : config.pDets) {
    auto dp0 = std::make_pair(det.first, 0);
    auto dp1 = std::make_pair(det.first, 1);
    long cnt =
        m_stats_detector[std::make_pair(det.first, "cluster_cnt_detector")][0];
    long cnt0 = 1;
    if (m_stats_plane[std::make_pair(dp0, "cluster_cnt_plane")][0] > 0) {
      cnt0 = m_stats_plane[std::make_pair(dp0, "cluster_cnt_plane")][0];
    }
    long cnt1 = 1;
    if (m_stats_plane[std::make_pair(dp1, "cluster_cnt_plane")][0] > 0) {
      cnt1 = m_stats_plane[std::make_pair(dp1, "cluster_cnt_plane")][0];
    }
    if (config.GetAxes(dp0) && config.GetAxes(dp1)) {
      bothPlanes = true;
    }

    std::cout << "\n\n****************************************" << std::endl;
    std::cout << "Stats detector " << (int)det.first << std::endl;
    std::cout << "****************************************" << std::endl;

    for (auto const &stat : m_stats_plane_names) {
      std::cout << "\n****************************************" << std::endl;
      std::cout << "Plane 0: " << stat << std::endl;
      std::cout << "****************************************" << std::endl;
      std::vector<long> v = m_stats_plane[std::make_pair(dp0, stat)];
      for (unsigned int n = 0; n < static_cast<unsigned int>(m_limits[stat]);
           n++) {
        StatsOutput(n, v[n], stat, cnt0);
      }
      std::cout << "****************************************" << std::endl;
      if (bothPlanes) {
        std::cout << "\n****************************************" << std::endl;
        std::cout << "Plane 1: " << stat << std::endl;
        std::cout << "****************************************" << std::endl;
        std::vector<long> v = m_stats_plane[std::make_pair(dp1, stat)];
        for (unsigned int n = 0; n < static_cast<unsigned int>(m_limits[stat]);
             n++) {
          StatsOutput(n, v[n], stat, cnt1);
        }
        std::cout << "****************************************" << std::endl;
      }
    }
    for (auto const &stat : m_stats_detector_names) {
      std::cout << "\n****************************************" << std::endl;
      std::cout << stat << std::endl;
      std::cout << "****************************************" << std::endl;
      std::vector<long> v = m_stats_detector[std::make_pair(det.first, stat)];
      for (unsigned int n = 0; n < static_cast<unsigned int>(m_limits[stat]);
           n++) {
        StatsOutput(n, v[n], stat, cnt, cnt0, cnt1);
      }
      std::cout << "****************************************" << std::endl;
    }
  }
}

void Statistics::PrintFECStats(Configuration &config) {
  for (auto const &fec : config.pFecs) {
    std::cout << "\n****************************************" << std::endl;
    std::cout << "FEC " << (int)fec << std::endl;
    std::cout << "****************************************" << std::endl;
    for (unsigned int n = 0; n < m_counter_names.size(); n++) {
      std::cout << m_counter_names[n] << ": "
                << GetCounter(m_counter_names[n], fec) << std::endl;
    }
    double first = GetFirstTriggerTimestamp(fec);
    double max = GetMaxTriggerTimestamp(fec);
    double last = GetOldTriggerTimestamp(fec);
    int overflow = GetCounter("time_stamp_overflow", fec);
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
    std::cout << "\n****************************************" << std::endl;
    std::cout << "Stats (acquisition):" << std::endl;
    std::cout << "****************************************" << std::endl;
    std::cout << "Acq time: " << m_acq_time << " ms" << std::endl;
    std::cout << "Hit rate FEC: " << std::scientific << 1000*GetCounter("parser_data", fec)/m_acq_time << " hit/s" << std::endl;
    std::cout << "Data rate FEC: " << std::scientific << 1000*GetCounter("parser_data", fec)*48/m_acq_time << " bit/s" << std::endl;
    std::cout << "****************************************" << std::endl;

  }
  std::cout << "\n****************************************" << std::endl;
  long cnt = 0;
    for (auto const &det : config.pDets) {
      cnt += GetStatsDetector("cluster_cnt_detector", det.first, 0); 
    }
    std::cout << "Total Cluster rate: " << std::scientific << 1000*cnt/m_acq_time << " particles/s" << std::endl;
  std::cout << "****************************************" << std::endl;
}

void Statistics::StatsOutput(int n, long val, std::string stat, long cnt,
                             long cnt0, long cnt1) {
  if (cnt == 0)
    cnt = 1;
  if (m_limits[stat] > 1) {
    if (m_factors[stat] != 1) {
      std::cout << static_cast<unsigned int>(n / m_factors[stat]) << "-"
                << static_cast<unsigned int>(n / m_factors[stat] +
                                             1 / m_factors[stat] - 1)
                << " " << m_units[stat] << ":  " << val << " ("
                << (100 * val / cnt) << " %)" << std::endl;
    } else {
      std::cout << static_cast<unsigned int>(n / m_factors[stat]) << " "
                << m_units[stat] << ":  " << val << " (" << (100 * val / cnt)
                << " %)" << std::endl;
    }
  } else {
    if (cnt0 > 0 && cnt1 > 0) {
      std::cout << val << " (common cluster in detector, " << (100 * val / cnt0)
                << " % plane 0, " << (100 * val / cnt1) << " % plane 1)"
                << std::endl;
    } else {
      std::cout << val << " (" << (100 * val / cnt) << " %)" << std::endl;
    }
  }
}