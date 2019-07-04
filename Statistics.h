#ifndef STATS_H
#define STATS_H

#include "Configuration.h"

class Statistics {
public:
    Statistics() = default;
    ~Statistics() = default;

    void CreateStats(Configuration & config); 
    void SetClusterStatsPlane(std::pair<uint8_t, uint8_t> dp, uint16_t maxDeltaTime, uint16_t maxMissingStrip, uint16_t deltaSpan);
    void SetStatsDetector(uint8_t det, double deltaPlane, double ratio, int plane);
    int maxDeltaTime(std::pair<uint8_t, uint8_t> dp, int n);
    int maxMissingStrip(std::pair<uint8_t, uint8_t> dp, int n);
    int deltaSpan(std::pair<uint8_t, uint8_t> dp, int n);
    int deltaPlane(uint8_t det, int n);
    int chargeRatio(uint8_t det, int n, int plane);

    void SetOldTriggerTimestamp(uint8_t fecId, double srsTimestamp);
    void SetDeltaTriggerTimestamp(uint8_t fecId, double val);
    void IncrementOverflow(uint8_t fecId);
    void IncrementTimeError(uint8_t fecId);
    void IncrementTimestampTooLarge(uint8_t fecId);
    double oldTriggerTimestamp(uint8_t fecId);
    double deltaTriggerTimestamp(uint8_t fecId);
    double lowestCommonTriggerTimestampDet(uint8_t det);
    void SetLowestCommonTriggerTimestampDet(uint8_t det, double val);
    double lowestCommonTriggerTimestampPlane(std::pair<uint8_t, uint8_t> dp);
    void SetLowestCommonTriggerTimestampPlane(std::pair<uint8_t, uint8_t> dp, double val);
    int timeError(uint8_t fecId);
    int overflow(uint8_t fecId);
    int timestampTooLargeError(uint8_t fecId);
private:
    std::map<std::pair<uint8_t, uint8_t>, std::vector<int>> m_maxDeltaTime;
    std::map<std::pair<uint8_t, uint8_t>, std::vector<int>> m_maxMissingStrip;
    std::map<std::pair<uint8_t, uint8_t>, std::vector<int>> m_deltaSpan;
    std::map<uint8_t, std::vector<int>> m_deltaPlane;
    std::map<uint8_t, std::vector<int>> m_chargeRatio_0;
    std::map<uint8_t, std::vector<int>> m_chargeRatio_1;
    std::map<uint8_t, int> m_overflow;
    std::map<uint8_t, int> m_timeError;
    std::map<uint8_t, int> m_timeStampTooLargeError;

    std::map<std::pair<uint8_t, uint8_t>, double> m_lowestCommonTriggerTimestamp_plane;
    std::map<uint8_t, double> m_lowestCommonTriggerTimestamp_det;
    std::map<uint8_t, double> m_deltaTriggerTimestamp;
    std::map<uint8_t, double> m_oldTriggerTimestamp;
};



#endif // STATS_H
