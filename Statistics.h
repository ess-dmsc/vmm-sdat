#ifndef STATS_H
#define STATS_H

#include "Configuration.h"

class Statistics {
public:
    Statistics() = default;
    ~Statistics() = default;

    void CreateStats(Configuration & config); 
    int GetStatsDetector(std::string stats, uint8_t det, int n);
    void SetStatsDetector(std::string stats, uint8_t det, double value);
    int GetStatsPlane(std::string stats, std::pair<uint8_t, uint8_t> dp, int n);
    void SetStatsPlane(std::string stats, std::pair<uint8_t, uint8_t> dp, double value);

    void IncrementErrorCount(std::string error, uint8_t fecId);
    int GetErrorCount(std::string error, uint8_t fecId);
   
    double GetDeltaTriggerTimestamp(uint8_t fecId);
    void SetDeltaTriggerTimestamp(uint8_t fecId, double val);
    double GetOldTriggerTimestamp(uint8_t fecId);
    void SetOldTriggerTimestamp(uint8_t fecId, double srsTimestamp);
    
    double GetLowestCommonTriggerTimestampDet(uint8_t det);
    void SetLowestCommonTriggerTimestampDet(uint8_t det, double val);
    double GetLowestCommonTriggerTimestampPlane(std::pair<uint8_t, uint8_t> dp);
    void SetLowestCommonTriggerTimestampPlane(std::pair<uint8_t, uint8_t> dp, double val);
   
    void PrintStats(Configuration& config);
    void StatsOutput(int n, int val, std::string stat, int cnt,int cnt0=0,int cnt1=0);
private:
    std::map<std::pair<std::pair<uint8_t, uint8_t>, std::string>, std::vector<int>> m_stats_plane;
    std::map<std::pair<uint8_t, std::string>, std::vector<int>> m_stats_detector;
    std::vector<std::string> m_stats_plane_names;
    std::vector<std::string> m_stats_detector_names;
    std::map<std::string, double> m_factors;
    std::map<std::string, double> m_limits;
    std::map<std::string, std::string> m_units;
    std::map<std::pair<uint8_t, std::string>, int> m_errors;
    std::vector<std::string> m_error_names;


    // per plane
    std::map<std::pair<uint8_t, uint8_t>, double> m_lowestCommonTriggerTimestamp_plane;
    // per detector
    std::map<uint8_t, double> m_lowestCommonTriggerTimestamp_det;
    // per FEC
    std::map<uint8_t, double> m_deltaTriggerTimestamp;
    std::map<uint8_t, double> m_oldTriggerTimestamp;
};



#endif // STATS_H
