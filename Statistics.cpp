#include <iostream>
#include <sstream>

#include <cstring>

#include "Statistics.h"

void Statistics::CreateStats(Configuration &config)
{

    for (auto const &fec : config.pFecs)
    {
        m_deltaTriggerTimestamp.emplace(std::make_pair(fec, 0));
        m_oldTriggerTimestamp.emplace(std::make_pair(fec, 0));
        m_overflow.emplace(std::make_pair(fec, 0));
        m_timeError.emplace(std::make_pair(fec, 0));
        m_timeStampTooLargeError.emplace(std::make_pair(fec, 0));
    }

    for (auto const &det : config.pDets)
    {
        auto plane0 = std::make_pair(det.first, 0);
        auto plane1 = std::make_pair(det.first, 1);

        m_lowestCommonTriggerTimestamp_det[det.first] = 0;
        m_lowestCommonTriggerTimestamp_plane[plane0] = 0;
        m_lowestCommonTriggerTimestamp_plane[plane1] = 0;
        for (unsigned int n = 0; n <= static_cast<unsigned int>(config.pDeltaTimeHits / 50); n++)
        {
            m_maxDeltaTime[plane0].push_back(0);
            m_maxDeltaTime[plane1].push_back(0);
        }
        for (unsigned int n = 0; n <= config.pMissingStripsCluster; n++)
        {
            m_maxMissingStrip[plane0].push_back(0);
            m_maxMissingStrip[plane1].push_back(0);
        }
        for (unsigned int n = 0; n <= static_cast<unsigned int>(config.pSpanClusterTime / 50); n++)
        {
            m_deltaSpan[plane0].push_back(0);
            m_deltaSpan[plane1].push_back(0);
        }
        for (unsigned int n = 0; n <= static_cast<unsigned int>(config.pDeltaTimePlanes / 50); n++)
        {

            m_deltaPlane[det.first].push_back(0);
        }
        for (unsigned int n = 0; n <= 10; n++)
        {
            m_chargeRatio_0[det.first].push_back(0);
            m_chargeRatio_1[det.first].push_back(0);
        }
    }
}

void Statistics::SetClusterStatsPlane(std::pair<uint8_t, uint8_t> dp, uint16_t maxDeltaTime, uint16_t maxMissingStrip, uint16_t deltaSpan)
{
    m_maxDeltaTime[dp][(unsigned int)(maxDeltaTime / 50)]++;
    m_maxMissingStrip[dp][(unsigned int)(maxMissingStrip)]++;
    m_deltaSpan[dp][(unsigned int)(deltaSpan / 50)]++;
}

void Statistics::SetStatsDetector(uint8_t det, double deltaPlane, double ratio, int plane)
{
    if(plane == 0)
    {
        m_chargeRatio_0[det][static_cast<unsigned int>(ratio * 10)]++;
    }
    else
    {
        m_chargeRatio_1[det][static_cast<unsigned int>(ratio * 10)]++;
    }
    
    m_deltaPlane[det][std::abs((int)deltaPlane) / 50]++;
    
}

int Statistics::maxDeltaTime(std::pair<uint8_t, uint8_t> dp, int n)
{
    return m_maxDeltaTime[dp][n];
}

int Statistics::maxMissingStrip(std::pair<uint8_t, uint8_t> dp, int n)
{
    return m_maxMissingStrip[dp][n];
}

int Statistics::deltaSpan(std::pair<uint8_t, uint8_t> dp, int n)
{
    return m_deltaSpan[dp][n];
}

int Statistics::deltaPlane(uint8_t det, int n)
{
    return m_deltaPlane[det][n];
}

int Statistics::chargeRatio(uint8_t det, int n, int plane)
{
    if(plane == 0)
    {
        return m_chargeRatio_0[det][n];
    }
    return m_chargeRatio_1[det][n];
}

void Statistics::SetOldTriggerTimestamp(uint8_t fecId, double srsTimestamp)
{
    m_oldTriggerTimestamp[fecId] = srsTimestamp;
}

void Statistics::SetDeltaTriggerTimestamp(uint8_t fecId, double val)
{
    m_deltaTriggerTimestamp[fecId] = val;
}

void Statistics::IncrementOverflow(uint8_t fecId)
{
    m_overflow[fecId]++;
}

void Statistics::IncrementTimeError(uint8_t fecId)
{
    m_timeError[fecId]++;
}

void Statistics::IncrementTimestampTooLarge(uint8_t fecId)
{
    m_timeStampTooLargeError[fecId]++;
}



int Statistics::timeError(uint8_t fecId)
{
    return m_timeError[fecId];
}

int Statistics::overflow(uint8_t fecId)
{
    return m_overflow[fecId];
}

int Statistics::timestampTooLargeError(uint8_t fecId)
{
    return m_timeStampTooLargeError[fecId];
}

double Statistics::oldTriggerTimestamp(uint8_t fecId)
{
    return m_oldTriggerTimestamp[fecId];
}

double Statistics::deltaTriggerTimestamp(uint8_t fecId)
{
    return m_deltaTriggerTimestamp[fecId];
}

double Statistics::lowestCommonTriggerTimestampDet(uint8_t det)
{
    return m_lowestCommonTriggerTimestamp_det[det];
}

void Statistics::SetLowestCommonTriggerTimestampDet(uint8_t det, double val)
{
    m_lowestCommonTriggerTimestamp_det[det] = val;
}

double Statistics::lowestCommonTriggerTimestampPlane(std::pair<uint8_t, uint8_t> dp)
{
    return m_lowestCommonTriggerTimestamp_plane[dp];
}

void Statistics::SetLowestCommonTriggerTimestampPlane(std::pair<uint8_t, uint8_t> dp, double val)
{
    m_lowestCommonTriggerTimestamp_plane[dp] = val;
}
