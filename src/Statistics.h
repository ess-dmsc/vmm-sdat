#ifndef STATS_H
#define STATS_H

#include "Configuration.h"

class Statistics {
public:
    Statistics() = default;
    ~Statistics() = default;

    void CreateClusterStats(Configuration & config); 
    void CreateFECStats(Configuration & config); 
    void CreatePCAPStats(Configuration & config); 
    
    long GetStatsDetector(std::string stats, uint8_t det, int n);
    void SetStatsDetector(std::string stats, uint8_t det, double value);
    long GetStatsPlane(std::string stats, std::pair<uint8_t, uint8_t> dp, int n);
    void SetStatsPlane(std::string stats, std::pair<uint8_t, uint8_t> dp, double value);

    void IncrementCounter(std::string error, uint8_t fecId, uint64_t increment=1);
    long GetCounter(std::string error, uint8_t fecId);
   
    double GetDeltaTriggerTimestamp(uint8_t fecId);
    void SetDeltaTriggerTimestamp(uint8_t fecId, double val);
    double GetOldTriggerTimestamp(uint8_t fecId);
    void SetOldTriggerTimestamp(uint8_t fecId, double srsTimestamp);
    double GetFirstTriggerTimestamp(uint8_t fecId);
    void SetFirstTriggerTimestamp(uint8_t fecId, double srsTimestamp);
    double GetMaxTriggerTimestamp(uint8_t fecId);
    void SetMaxTriggerTimestamp(uint8_t fecId, double srsTimestamp);
    uint64_t GetLastFrameCounter(uint8_t fecId);
    void SetLastFrameCounter(uint8_t fecId, uint64_t frameCounter);
    double GetLowestCommonTriggerTimestampDet(uint8_t det);
    void SetLowestCommonTriggerTimestampDet(uint8_t det, double val);
    double GetLowestCommonTriggerTimestampPlane(std::pair<uint8_t, uint8_t> dp);
    void SetLowestCommonTriggerTimestampPlane(std::pair<uint8_t, uint8_t> dp, double val);
   
    void PrintClusterStats(Configuration& config);
    void PrintFECStats(Configuration& config);
    
    void StatsOutput(int n, long val, std::string stat, long cnt,long cnt0=0,long cnt1=0);
private:
    double m_acq_time;
    std::map<std::pair<std::pair<uint8_t, uint8_t>, std::string>, std::vector<long>> m_stats_plane;
    std::map<std::pair<uint8_t, std::string>, std::vector<long>> m_stats_detector;
    std::vector<std::string> m_stats_plane_names;
    std::vector<std::string> m_stats_detector_names;
    std::map<std::string, double> m_factors;
    std::map<std::string, double> m_limits;
    std::map<std::string, std::string> m_units;
    std::map<std::pair<uint8_t, std::string>, long> m_counters;
    std::vector<std::string> m_counter_names;
   

    // per plane
    std::map<std::pair<uint8_t, uint8_t>, double> m_lowestCommonTriggerTimestamp_plane;
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



#endif // STATS_H
