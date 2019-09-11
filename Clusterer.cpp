#include <algorithm>
#include <cmath>
#include "Clusterer.h"
#include "Trace.h"

#include <chrono>
#include <functional>

#define UNUSED __attribute__((unused))

//#undef TRC_LEVEL
//#define TRC_LEVEL TRC_L_DEB

auto now = std::chrono::steady_clock::now;

auto timethis(std::function<void()> thunk) -> decltype((now() - now()).count())
{
    auto start = now();
    thunk();
    auto stop = now();
    return (stop - start).count();
}

Clusterer::Clusterer(Configuration &config, Statistics &stats) : m_config(config), m_stats(stats)
{
    m_rootFile = RootFile::GetInstance(config);
}

Clusterer::~Clusterer()
{
    RootFile::Dispose();
}

//====================================================================================================================
bool Clusterer::AnalyzeHits(double srsTimestamp, uint8_t fecId, uint8_t vmmId, uint16_t chNo, uint16_t bcid, uint16_t tdc, uint16_t adc, bool overThresholdFlag, float chipTime)
{
    auto searchFecChip = m_config.pFecChip_DetectorPlane.find(std::make_pair(fecId, vmmId));
    if (searchFecChip == m_config.pFecChip_DetectorPlane.end())
    {
        DTRACE(DEB, "\t\tDetector or Plane not defined for FEC %d and vmmId %d not defined!\n", (int)fecId, (int)vmmId);
        return true;
    }
    auto det_plane = GetDetectorPlane(std::make_pair(fecId, vmmId));
    auto det = std::get<0>(det_plane);
    auto plane = std::get<1>(det_plane);

    // Plane 0: x
    // plane 1: y
    int pos = GetChannel(std::make_pair(fecId, vmmId), chNo);
    if (pos == -1)
    {
        DTRACE(DEB, "\t\tChannel not defined for FEC %d and vmmId %d not defined!\n", (int)fecId, (int)vmmId);
        return true;
    }
    m_stats.SetDeltaTriggerTimestamp(fecId, 0);
    if (srsTimestamp > 0xFFFFFFFFFFF)
    {
        m_stats.IncrementTimestampTooLarge(fecId);
        DTRACE(DEB, "\t\tTimestamp %llu larger than 42 bit and 31 times trigger periodd for FEC %d and vmmId %d!\n", static_cast<uint64_t>(srsTimestamp), (int)fecId, (int)vmmId);
        return true;
    }
    if (srsTimestamp < m_stats.oldTriggerTimestamp(fecId))
    {
        if (m_stats.oldTriggerTimestamp(fecId) - srsTimestamp > 0x7FFFFFFFFFF)
        {
            m_stats.IncrementOverflow(fecId);
            DTRACE(DEB, "\n*********************************** OVERFLOW  fecId %d, m_lineNr %d, eventNr  %d, "
                        "srsTimestamp %llu, old srsTimestamp %llu\n",
                   fecId, m_lineNr, m_eventNr, static_cast<uint64_t>(srsTimestamp), static_cast<uint64_t>(m_stats.oldTriggerTimestamp(fecId)));
        }
        else
        {
            m_stats.IncrementTimeError(fecId);
            DTRACE(DEB, "\n*********************************** TIME ERROR  fecId %d, m_lineNr %d, eventNr  %d, "
                        "srsTimestamp %llu, old srsTimestamp %llu\n",
                   fecId, m_lineNr, m_eventNr, static_cast<uint64_t>(srsTimestamp), static_cast<uint64_t>(m_stats.oldTriggerTimestamp(fecId)));
        }
    }

    double remainder = std::fmod(m_stats.deltaTriggerTimestamp(fecId), m_config.pTriggerPeriod);
    if (remainder > 0)
    {
        uint64_t offset = m_stats.deltaTriggerTimestamp(fecId) / m_config.pTriggerPeriod;
        printf("\n******* ABORT ANALYSIS: SRS timestamp wrong increment: fec %d,vmmId %d, chNo %d, line %d, "
               "trigger period %d, offset %llu, remainder %llu, new time %llu, old time %llu\n",
               fecId, vmmId, chNo, m_lineNr, m_config.pTriggerPeriod, offset, static_cast<uint64_t>(remainder),
               static_cast<uint64_t>(srsTimestamp), static_cast<uint64_t>(m_stats.oldTriggerTimestamp(fecId)));
        return false;
    }

    bool newEvent = false;
    if (srsTimestamp > m_stats.oldTriggerTimestamp(fecId))
    {
        m_stats.SetDeltaTriggerTimestamp(fecId, srsTimestamp - m_stats.oldTriggerTimestamp(fecId));
        newEvent = true;
    }

    if (newEvent)
    {
        m_eventNr++;

        int factor = 10;
        for (auto const &searchDetPlaneFec : m_config.pDetectorPlane_Fec)
        {
            auto fec = searchDetPlaneFec.second;
            auto det_plane = searchDetPlaneFec.first;
            auto det = std::get<0>(det_plane);
            auto plane = std::get<1>(det_plane);
            auto dp0 = std::make_pair(det, 0);
            auto dp1 = std::make_pair(det, 1);
            if (m_stats.oldTriggerTimestamp(fecId) > factor * m_config.pTriggerPeriod && m_stats.lowestCommonTriggerTimestampPlane(det_plane) < m_stats.oldTriggerTimestamp(fecId))
            {
                m_stats.SetLowestCommonTriggerTimestampPlane(det_plane, m_stats.oldTriggerTimestamp(fecId) - factor * m_config.pTriggerPeriod);
            }
            m_stats.SetLowestCommonTriggerTimestampDet(det, std::min(m_stats.lowestCommonTriggerTimestampPlane(dp0), m_stats.lowestCommonTriggerTimestampPlane(dp1)));
        }
        for (auto const &searchDet : m_config.pDets)
        {
            AnalyzeClustersPlane(searchDet.first, 0);
            AnalyzeClustersPlane(searchDet.first, 1);
            AnalyzeClustersDetector(searchDet.first);
        }
    }

    // TDC has reduced resolution due to most significant bit problem of current
    // sources (like ADC)
    //int theTDC =  8 * (tdc / 8) + 4;
    //theTDC = tdc;

    //uint64_t bcTime = pBCTime_ns * (bcid + 1);
    //uint64_t tdcTime = tdc * pTAC / 256;
    m_lineNr++;
    //uint64_t theChiptime = (bcTime - tdcTime);
    uint64_t totalTime = srsTimestamp + (int)chipTime;
    if (m_config.pCreateHits)
    {
        Hit theHit;
        theHit.id = m_lineNr;
        theHit.event = m_eventNr;
        theHit.det = det;
        theHit.plane = plane;
        theHit.fec = fecId;
        theHit.vmm = vmmId;
        theHit.readout_time = srsTimestamp;
        theHit.ch = chNo;
        theHit.pos = (uint16_t)pos;

        theHit.bcid = bcid;
        theHit.tdc = tdc;
        theHit.adc = adc;
        theHit.over_threshold = overThresholdFlag;
        theHit.chip_time = chipTime;
        theHit.time = totalTime;
        m_rootFile->AddHits(std::move(theHit));
    }

    StoreHits(det, plane, pos, adc, totalTime, overThresholdFlag);

    if (newEvent)
    {
        DTRACE(DEB, "\neventNr  %d\n", m_eventNr);
        //DTRACE(DEB, "fecId  %d\n", fecId);
    }
    if (m_stats.deltaTriggerTimestamp(fecId) > 0)
    {
        DTRACE(DEB, "\tTriggerTimestamp %llu [ns]\n", static_cast<uint64_t>(srsTimestamp));
        DTRACE(DEB, "\tTime since last trigger %f us (%.4f kHz)\n", m_stats.deltaTriggerTimestamp(fecId) * 0.001, (double)(1000000 / m_stats.deltaTriggerTimestamp(fecId)));
    }

    if (m_oldFecId != fecId || newEvent)
    {
        DTRACE(DEB, "\tfecId  %d\n", fecId);
    }
    if (m_oldVmmId != vmmId || newEvent)
    {
        DTRACE(DEB, "\tDetector %d, plane %d, vmmId  %d\n", (int)det, (int)plane, vmmId);
    }
    DTRACE(DEB, "\t\tChannel %d (chNo  %d) - overThresholdFlag %d\n", pos, chNo, (int)overThresholdFlag);
    DTRACE(DEB, "\t\t\tbcid %d, tdc %d, adc %d\n", bcid, tdc, adc);
    DTRACE(DEB, "\t\t\ttotal time %llu, chip time %f ns\n", totalTime, chipTime);

    m_stats.SetOldTriggerTimestamp(fecId, srsTimestamp);
    m_oldBcId = bcid;
    m_oldVmmId = vmmId;
    m_oldFecId = fecId;

    return true;
}

//====================================================================================================================
void Clusterer::StoreHits(uint8_t det, uint8_t plane, int pos, uint16_t adc, double totalTime, bool overThresholdFlag)
{
    if ((adc >= m_config.pADCThreshold || overThresholdFlag))
    {
        m_hits_new[std::make_pair(det, plane)].emplace_back(totalTime, (uint16_t)pos, adc);
    }
}

//====================================================================================================================
int Clusterer::ClusterByTime(uint8_t det, uint8_t plane)
{

    std::pair<uint8_t, uint8_t> dp = std::make_pair(det, plane);

    ClusterContainer cluster;
    uint16_t maxDeltaTime = 0;
    int clusterCount = 0;
    double time1 = 0, time2 = 0;
    uint32_t adc1 = 0;
    uint16_t strip1 = 0;

    for (auto &itHits : m_hits[dp])
    {
        time2 = time1;

        time1 = (double)std::get<0>(itHits);
        strip1 = std::get<1>(itHits);
        adc1 = std::get<2>(itHits);
        if (!cluster.empty())
        {
            if (abs(time1 - time2) > m_config.pDeltaTimeHits)
            {
                clusterCount += ClusterByStrip(det, plane, cluster, maxDeltaTime);
                cluster.clear();
                maxDeltaTime = 0;
            }
            else
            {
                if (maxDeltaTime < abs(time1 - time2))
                {
                    maxDeltaTime = (time1 - time2);
                }
            }
        }
        cluster.emplace_back(strip1, time1, adc1);
    }

    if (!cluster.empty())
    {
        clusterCount += ClusterByStrip(det, plane, cluster, maxDeltaTime);
    }
    return clusterCount;
}

//====================================================================================================================
int Clusterer::ClusterByStrip(uint8_t det, uint8_t plane, ClusterContainer &cluster, uint16_t maxDeltaTime)
{
    int maxMissingStrip = 0;
    uint16_t spanCluster = 0;

    double startTime = 0;
    double largestTime = 0;
    float position_utpc = -1;
    double centerOfGravity = 0;
    double centerOfTime = 0;
    double centerOfGravity2 = 0;
    double centerOfTime2 = 0;
    long int totalADC = 0;
    long int totalADC2 = 0;

    double time1 = 0;
    int adc1 = 0;
    int strip1 = 0;
    int strip2 = 0;
    int stripCount = 0;
    int clusterCount = 0;
    std::vector<double> vStrips;
    std::vector<double> vTimes;
    std::sort(begin(cluster), end(cluster), [](const ClusterTuple &t1, const ClusterTuple &t2) {
        return std::get<0>(t1) < std::get<0>(t2);
    });
    for (auto &itCluster : cluster)
    {
        strip2 = strip1;
        strip1 = std::get<0>(itCluster);
        time1 = std::get<1>(itCluster);
        adc1 = std::get<2>(itCluster);

        // At beginning of cluster, set start time of cluster
        if (stripCount == 0)
        {
            maxMissingStrip = 0;
            startTime = time1;
            largestTime = time1;
            DTRACE(DEB, "\nDetector %d, plane %d cluster:\n", (int)det, (int)plane);
        }
        // Add members of a cluster, if it is either the beginning of a cluster,
        // or if strip gap and time span is correct
        if (stripCount == 0 || (std::abs(strip1 - strip2) > 0 && std::abs(strip1 - strip2) - 1 <= m_config.pMissingStripsCluster && time1 - startTime <= m_config.pSpanClusterTime && largestTime - time1 <= m_config.pSpanClusterTime))
        {
            DTRACE(DEB, "\tstrip %d, time %llu, adc %d:\n", strip1, (uint64_t)time1, adc1);

            if (time1 >= largestTime)
            {
                largestTime = time1;
                position_utpc = strip1;
            }
            if (time1 <= startTime)
            {
                startTime = time1;
            }
            if (stripCount > 0 && maxMissingStrip < std::abs(strip1 - strip2) - 1)
            {
                maxMissingStrip = std::abs(strip1 - strip2) - 1;
            }
            spanCluster = (largestTime - startTime);
            totalADC += adc1;
            totalADC2 += adc1 * adc1;
            centerOfGravity += strip1 * adc1;
            centerOfTime += time1 * adc1;
            centerOfGravity2 += strip1 * adc1 * adc1;
            centerOfTime2 += time1 * adc1 * adc1;
            vStrips.emplace_back(strip1);
            vTimes.emplace_back(time1);
            stripCount++;
        }
        // Stop clustering if gap between strips is too large or time span too long
        else if (std::abs(strip1 - strip2) - 1 > m_config.pMissingStripsCluster || time1 - startTime > m_config.pSpanClusterTime || largestTime - time1 > m_config.pSpanClusterTime)
        {
            // Valid cluster
            if (stripCount < m_config.pMinClusterSize || totalADC == 0)
            {
                DTRACE(DEB, "******** INVALID ********\n\n");
            }
            else
            {
                spanCluster = (largestTime - startTime);
                centerOfGravity = (centerOfGravity / totalADC);
                centerOfTime = (centerOfTime / totalADC);
                centerOfGravity2 = (centerOfGravity2 / totalADC2);
                centerOfTime2 = (centerOfTime2 / totalADC2);
                StoreClusters(det, plane, vStrips, vTimes, position_utpc, largestTime, centerOfGravity, centerOfTime, centerOfGravity2, centerOfTime2, stripCount, totalADC, maxDeltaTime, maxMissingStrip, spanCluster);
                clusterCount++;
                maxMissingStrip = 0;
            }
            // Reset all parameters
            startTime = 0;
            largestTime = 0;
            stripCount = 0;
            centerOfGravity = 0;
            centerOfTime = 0;
            totalADC = 0;
            centerOfGravity2 = 0;
            centerOfTime2 = 0;
            totalADC2 = 0;
            strip1 = 0;
        }
    }
    // At the end of the clustering, check again if there is a last valid cluster
    if (stripCount < m_config.pMinClusterSize || totalADC == 0)
    {
        DTRACE(DEB, "******** INVALID ********\n\n");
    }
    else
    {
        spanCluster = (largestTime - startTime);
        centerOfGravity = (centerOfGravity / (float)totalADC);
        centerOfTime = (centerOfTime / (float)totalADC);
        centerOfGravity2 = (centerOfGravity2 / (float)totalADC2);
        centerOfTime2 = (centerOfTime2 / (float)totalADC2);
        StoreClusters(det, plane, vStrips, vTimes, position_utpc, largestTime, centerOfGravity, centerOfTime, centerOfGravity2, centerOfTime2, stripCount, totalADC, maxDeltaTime, maxMissingStrip, spanCluster);
        clusterCount++;
    }
    return clusterCount;
}

//====================================================================================================================
void Clusterer::StoreClusters(uint8_t det, uint8_t plane, std::vector<double> &strips, std::vector<double> &times,
                              double pos_utpc, double time_utpc,
                              double pos_center_charge, double time_center_charge,
                              double pos_center_charge2, double time_center_charge2,
                              uint16_t clusterSize, uint32_t clusterADC, uint16_t maxDeltaTime, uint16_t maxMissingStrip, uint16_t deltaSpan)
{
    ClusterPlane clusterPlane;
    m_cluster_id++;

    DTRACE(DEB, "Cluster id %d\n", m_cluster_id);
    clusterPlane.id = static_cast<uint32_t>(m_cluster_id);
    clusterPlane.det = det;
    clusterPlane.plane = plane;
    clusterPlane.size = clusterSize;
    clusterPlane.adc = clusterADC;
    clusterPlane.time_utpc = time_utpc;
    clusterPlane.pos_utpc = pos_utpc;
    clusterPlane.time = time_center_charge;
    clusterPlane.pos = pos_center_charge;
    clusterPlane.time_charge2 = time_center_charge2;
    clusterPlane.pos_charge2 = pos_center_charge2;
    clusterPlane.plane_coincidence = false;
    clusterPlane.max_delta_time = maxDeltaTime;
    clusterPlane.max_missing_strip = maxMissingStrip;
    clusterPlane.span_cluster = deltaSpan;
    clusterPlane.strips = std::move(strips);
    clusterPlane.times = std::move(times);
    if (clusterPlane.pos > -1.0)
    {
        auto dp = std::make_pair(det, plane);
        m_stats.SetClusterStatsPlane(dp, maxDeltaTime, maxMissingStrip, deltaSpan);
        m_clusters_new[dp].emplace_back(std::move(clusterPlane));
    }
    else
    {
        std::cout << (int)clusterPlane.id << " " << (int)clusterPlane.det << " " << (int)clusterPlane.plane << " " << clusterPlane.pos << std::endl;
    }
}

//====================================================================================================================
void Clusterer::MatchClustersDetector(uint8_t det)
{

    ClusterVectorPlane::iterator itStartPlane1 = begin(m_clusters[std::make_pair(det, 1)]);

    for (auto &c0 : m_clusters[std::make_pair(det, 0)])
    {
        double minDelta = 99999999;
        double lastDelta_t = 99999999;
        double delta_t = 99999999;
        bool isFirstMatch = true;
        ClusterVectorPlane::iterator bestMatchPlane1 = end(m_clusters[std::make_pair(det, 1)]);

        for (ClusterVectorPlane::iterator c1 = itStartPlane1;
             c1 != end(m_clusters[std::make_pair(det, 1)]); ++c1)
        {
            if ((*c1).plane_coincidence == false)
            {

                double chargeRatio = (double)(*c1).adc / (double)c0.adc;
                lastDelta_t = delta_t;
                delta_t = (*c1).time - c0.time;
                if (m_config.pConditionCoincidence == "utpc")
                {
                    delta_t = (*c1).time_utpc - c0.time_utpc;
                }
                else if (m_config.pConditionCoincidence == "charge2")
                {
                    delta_t = (*c1).time_charge2 - c0.time_charge2;
                }
                if (chargeRatio >= 1 / m_config.pChargeRatio && chargeRatio <= m_config.pChargeRatio && std::abs(delta_t) < minDelta && std::abs(delta_t) <= m_config.pDeltaTimePlanes && (c0.size + (*c1).size >= m_config.pCoincidentClusterSize))
                {
                    minDelta = std::abs(delta_t);
                    bestMatchPlane1 = c1;
                    if (isFirstMatch)
                    {
                        itStartPlane1 = c1;
                        isFirstMatch = false;
                    }
                }
                if (std::abs(delta_t) > std::abs(lastDelta_t))
                {
                    break;
                }
            }
        }

        if (bestMatchPlane1 != end(m_clusters[std::make_pair(det, 1)]))
        {
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
 
            if (m_config.pTransform.size() == m_config.pDets.size())
            {
                auto tx = m_config.pTransformX[m_config.pDets[det]];
                auto ty = m_config.pTransformY[m_config.pDets[det]];
                auto tz = m_config.pTransformZ[m_config.pDets[det]];

                clusterDetector.pos0 = c0.pos * std::get<0>(tx) + (*bestMatchPlane1).pos * std::get<1>(tx) + std::get<3>(tx);
                clusterDetector.pos1 = c0.pos * std::get<0>(ty) + (*bestMatchPlane1).pos * std::get<1>(ty) + std::get<3>(ty);
                clusterDetector.pos2 = c0.pos * std::get<0>(tz) + (*bestMatchPlane1).pos * std::get<1>(tz) + std::get<3>(tz);

                clusterDetector.pos0_utpc = c0.pos_utpc * std::get<0>(tx) + (*bestMatchPlane1).pos_utpc * std::get<1>(tx) + std::get<3>(tx);
                clusterDetector.pos1_utpc = c0.pos_utpc * std::get<0>(ty) + (*bestMatchPlane1).pos_utpc * std::get<1>(ty) + std::get<3>(ty);
                clusterDetector.pos2_utpc = c0.pos_utpc * std::get<0>(tz) + (*bestMatchPlane1).pos_utpc * std::get<1>(tz) + std::get<3>(tz);

                clusterDetector.pos0_charge2 = c0.pos_charge2 * std::get<0>(tx) + (*bestMatchPlane1).pos_charge2 * std::get<1>(tx) + std::get<3>(tx);
                clusterDetector.pos1_charge2 = c0.pos_charge2 * std::get<0>(ty) + (*bestMatchPlane1).pos_charge2 * std::get<1>(ty) + std::get<3>(ty);
                clusterDetector.pos2_charge2 = c0.pos_charge2 * std::get<0>(tz) + (*bestMatchPlane1).pos_charge2 * std::get<1>(tz) + std::get<3>(tz);
            }
            else
            {
                clusterDetector.pos0 = c0.pos;
                clusterDetector.pos1 = (*bestMatchPlane1).pos;
                clusterDetector.pos2 = 0;
                clusterDetector.pos0_utpc = c0.pos_utpc;
                clusterDetector.pos1_utpc = (*bestMatchPlane1).pos_utpc;
                clusterDetector.pos2_utpc = 0;
                clusterDetector.pos0_charge2 = c0.pos_charge2;
                clusterDetector.pos1_charge2 = (*bestMatchPlane1).pos_charge2;
                clusterDetector.pos2_charge2 = 0;
            }

            clusterDetector.time0 = c0.time;
            clusterDetector.time1 = (*bestMatchPlane1).time;
            clusterDetector.time0_utpc = c0.time_utpc;
            clusterDetector.time1_utpc = (*bestMatchPlane1).time_utpc;
            clusterDetector.time0_charge2 = c0.time_charge2;
            clusterDetector.time1_charge2 = (*bestMatchPlane1).time_charge2;
            clusterDetector.dt0 = clusterDetector.time0 - last_time0;
            clusterDetector.dt1 = clusterDetector.time1 - last_time1;
            /*
            clusterDetector.dt0_utpc = clusterDetector.time0_utpc - last_time0_utpc;
            clusterDetector.dt1_utpc = clusterDetector.time1_utpc - last_time1_utpc;
            clusterDetector.dt0_charge2 = clusterDetector.time0_charge2 - last_time0_charge2;
            clusterDetector.dt1_charge2 = clusterDetector.time1_charge2 - last_time1_charge2;
            */
            last_time0 = clusterDetector.time0;
            last_time1 = clusterDetector.time1;
            /*
            last_time0_utpc = clusterDetector.time0_utpc;
            last_time1_utpc = clusterDetector.time1_utpc;
            last_time0_charge2 = clusterDetector.time0_charge2;
            last_time1_charge2 = clusterDetector.time1_charge2;
            */
            clusterDetector.delta_plane = clusterDetector.time1 - clusterDetector.time0;

            if (m_config.pConditionCoincidence == "utpc")
            {
                clusterDetector.delta_plane = clusterDetector.time1_utpc - clusterDetector.time0_utpc;
            }
            else if (m_config.pConditionCoincidence == "charge2")
            {
                clusterDetector.delta_plane = clusterDetector.time1_charge2 - clusterDetector.time0_charge2;
            }

            double ratio = (double)clusterDetector.adc0 / (double)clusterDetector.adc1;
            if (ratio > 1.0)
            {
                ratio = (double)clusterDetector.adc1 / (double)clusterDetector.adc0;
                m_stats.SetStatsDetector(det, clusterDetector.delta_plane, ratio, 1);
            }
            else
            {
                m_stats.SetStatsDetector(det, clusterDetector.delta_plane, ratio, 0);
            }
            clusterDetector.max_delta_time0 = c0.max_delta_time;
            clusterDetector.max_delta_time1 = (*bestMatchPlane1).max_delta_time;
            clusterDetector.max_missing_strip0 = c0.max_missing_strip;
            clusterDetector.max_missing_strip1 = (*bestMatchPlane1).max_missing_strip;
            clusterDetector.span_cluster0 = c0.span_cluster;
            clusterDetector.span_cluster1 = (*bestMatchPlane1).span_cluster;
            clusterDetector.strips0 = c0.strips;
            clusterDetector.times0 = c0.times;
            clusterDetector.strips1 = (*bestMatchPlane1).strips;
            clusterDetector.times1 = (*bestMatchPlane1).times;

            DTRACE(DEB, "\ncommon cluster det %d x/y: %d/%d", (int)det, clusterDetector.id0, clusterDetector.id1);
            DTRACE(DEB, "\tpos x/pos y: %f/%f", clusterDetector.pos0, clusterDetector.pos1);
            DTRACE(DEB, "\ttime x/time y: : %llu/%llu", (uint64_t)clusterDetector.time0, (uint64_t)clusterDetector.time1);
            DTRACE(DEB, "\tadc x/adc y: %u/%u", clusterDetector.adc0, clusterDetector.adc1);
            DTRACE(DEB, "\tsize x/size y: %u/%u", clusterDetector.size0, clusterDetector.size1);
            DTRACE(DEB, "\tdelta time planes: %d", (int)clusterDetector.delta_plane);
            m_clusters_detector[det].emplace_back(std::move(clusterDetector));
        }
    }
}

//void Clusterer::AnalyzeClustersPlane(uint8_t det, uint8_t plane, HitContainer& hits, HitContainer& newHits, double timeReadyToCluster, uint16_t correctionTime)
void Clusterer::AnalyzeClustersPlane(uint8_t det, uint8_t plane)
{

    if (ChooseHitsToBeClustered(det, plane) == false)
    {
        return;
    }

    int cnt = ClusterByTime(det, plane);

    DTRACE(DEB, "%d cluster in detector %d plane %d\n", cnt, (int)det, (int)plane);

    std::pair<uint8_t, uint8_t> dp = std::make_pair(det, plane);
    if (!m_hits[dp].empty())
    {
        m_hits[dp].clear();
    }
}

void Clusterer::AnalyzeClustersDetector(uint8_t det)
{

    if (ChooseClustersToBeMatched(det, 0) == false)
    {
        return;
    }
    if (m_config.GetAxes(det, 0) && m_config.GetAxes(det, 1))
    {  
        if (ChooseClustersToBeMatched(det, 1) == false)
        {
            return;
        }

        MatchClustersDetector(det);
    }
    m_clusters_detector_cnt[det] += m_clusters_detector[det].size();
    m_clusters_cnt[std::make_pair(det, 0)] += m_clusters[std::make_pair(det, 0)].size();
    m_clusters_cnt[std::make_pair(det, 1)] += m_clusters[std::make_pair(det, 1)].size();

    m_rootFile->SaveHits();
    m_rootFile->SaveClustersPlane(std::move(m_clusters[std::make_pair(det, 0)]));
    m_rootFile->SaveClustersPlane(std::move(m_clusters[std::make_pair(det, 1)]));
    m_rootFile->SaveClustersDetector(std::move(m_clusters_detector[det]));

    m_clusters[std::make_pair(det, 0)].clear();
    m_clusters[std::make_pair(det, 1)].clear();
    m_clusters_detector[det].clear();
}

//====================================================================================================================
bool Clusterer::ChooseHitsToBeClustered(uint8_t det, uint8_t plane)
{

    std::pair<uint8_t, uint8_t> dp = std::make_pair(det, plane);
    double timeReadyToCluster = m_stats.lowestCommonTriggerTimestampPlane(dp);
    //Nothing to cluster, newHits vector empty
    if (m_hits_new[dp].empty())
    {
        return false;
    }

    auto theMin = std::min_element(m_hits_new[dp].begin(), m_hits_new[dp].end(), [](const HitTuple &t1, const HitTuple &t2) {
        return std::get<0>(t1) < std::get<0>(t2);
    });

    //Nothing to cluster, tuples in newHits vector too recent
    if (std::get<0>(*theMin) > timeReadyToCluster)
    {

        //(smallest timestamp larger than m_stats.lowestCommonTriggerTimestampPlane(dp))
        //Will be clustered later
        return false;
    }

    //Sort vector newHits
    std::sort(begin(m_hits_new[dp]), end(m_hits_new[dp]), [](const HitTuple &t1, const HitTuple &t2) {
        return std::get<0>(t1) < std::get<0>(t2);
    });

    //First tuple with timestamp larger than m_stats.lowestCommonTriggerTimestampPlane(dp)
    auto it = std::upper_bound(m_hits_new[dp].begin(), m_hits_new[dp].end(), std::make_tuple(m_stats.lowestCommonTriggerTimestampPlane(dp), 0, 0), [](const HitTuple &t1, const HitTuple &t2) {
        return std::get<0>(t1) < std::get<0>(t2);
    });

    //Find elements in vector that could still be part of a cluster,
    //since they are close in time to m_stats.lowestCommonTriggerTimestampPlane(dp)
    while (it != m_hits_new[dp].end())
    {
        if (std::get<0>(*it) - timeReadyToCluster > m_config.pDeltaTimeHits)
        {
            break;
        }
        timeReadyToCluster = std::get<0>(*it);
        ++it;
    }

    //std::cout << "still in cluster " << std::get < 0 > (*it) << "\n" << std::endl;

    int index = std::distance(m_hits_new[dp].begin(), it);
    //Insert the data that is ready to be clustered from newHits into hits
    m_hits[dp].insert(m_hits[dp].end(), std::make_move_iterator(m_hits_new[dp].begin()), std::make_move_iterator(m_hits_new[dp].begin() + index));
    //Delete the data from newHits
    m_hits_new[dp].erase(m_hits_new[dp].begin(), m_hits_new[dp].begin() + index);

    return true;
}

bool Clusterer::ChooseClustersToBeMatched(uint8_t det, uint8_t plane)
{
    int index = 0;
    std::pair<uint8_t, uint8_t> dp = std::make_pair(det, plane);
    double timeReadyToMatch = m_stats.lowestCommonTriggerTimestampPlane(dp);

    //Nothing to match, newClusters vector empty
    if (m_clusters_new[dp].empty())
    {
        return false;
    }
    if (m_config.pConditionCoincidence == "utpc")
    {
        auto theMin = std::min_element(m_clusters_new[dp].begin(), m_clusters_new[dp].end(), [](const ClusterPlane &t1, const ClusterPlane &t2) {
            return t1.time_utpc < t2.time_utpc;
        });

        //Nothing to cluster, clusters in newClusters vector too recent
        if ((*theMin).time_utpc > timeReadyToMatch)
        {

            //(smallest time larger than timeReadyToMatch)
            //Will be matched later
            return false;
        }

        //Sort vector newClusters based on time
        std::sort(begin(m_clusters_new[dp]), end(m_clusters_new[dp]), [](const ClusterPlane &t1, const ClusterPlane &t2) {
            return t1.time_utpc < t2.time_utpc;
        });

        ClusterPlane theCluster;
        theCluster.time_utpc = timeReadyToMatch;

        //First ClusterPlane with time bigger than timeReadyToMatch
        auto it = std::upper_bound(m_clusters_new[dp].begin(), m_clusters_new[dp].end(), theCluster, [](const ClusterPlane &t1, const ClusterPlane &t2) {
            return t1.time_utpc < t2.time_utpc;
        });

        //Find elements in vector that could still be matched with another cluster
        //since they are close in time to timeReadyToMatch
        while (it != m_clusters_new[dp].end())
        {
            if ((*it).time_utpc - timeReadyToMatch > m_config.pDeltaTimeHits)
            {
                break;
            }
            timeReadyToMatch = (*it).time_utpc;
            ++it;
        }
        index = std::distance(m_clusters_new[dp].begin(), it);
    }
    else if (m_config.pConditionCoincidence == "charge2")
    {
        auto theMin = std::min_element(m_clusters_new[dp].begin(), m_clusters_new[dp].end(), [](const ClusterPlane &t1, const ClusterPlane &t2) {
            return t1.time_charge2 < t2.time_charge2;
        });

        //Nothing to cluster, clusters in newClusters vector too recent
        if ((*theMin).time_charge2 > timeReadyToMatch)
        {

            //(smallest time larger than timeReadyToMatch)
            //Will be matched later
            return false;
        }

        //Sort vector newClusters based on time
        std::sort(begin(m_clusters_new[dp]), end(m_clusters_new[dp]), [](const ClusterPlane &t1, const ClusterPlane &t2) {
            return t1.time_charge2 < t2.time_charge2;
        });

        ClusterPlane theCluster;
        theCluster.time_charge2 = timeReadyToMatch;

        //First ClusterPlane with time that bigger than timeReadyToMatch
        auto it = std::upper_bound(m_clusters_new[dp].begin(), m_clusters_new[dp].end(), theCluster, [](const ClusterPlane &t1, const ClusterPlane &t2) {
            return t1.time_charge2 < t2.time_charge2;
        });

        //Find elements in vector that could still be matched with another cluster
        //since they are close in time to timeReadyToMatch
        while (it != m_clusters_new[dp].end())
        {
            if ((*it).time_charge2 - timeReadyToMatch > m_config.pDeltaTimeHits)
            {
                break;
            }
            timeReadyToMatch = (*it).time_charge2;
            ++it;
        }
        index = std::distance(m_clusters_new[dp].begin(), it);
    }
    else
    {
        auto theMin = std::min_element(m_clusters_new[dp].begin(), m_clusters_new[dp].end(), [](const ClusterPlane &t1, const ClusterPlane &t2) {
            return t1.time < t2.time;
        });

        //Nothing to cluster, clusters in newClusters vector too recent
        if ((*theMin).time > timeReadyToMatch)
        {

            //(smallest time larger than timeReadyToMatch)
            //Will be matched later
            return false;
        }

        //Sort vector newClusters based on time
        std::sort(begin(m_clusters_new[dp]), end(m_clusters_new[dp]), [](const ClusterPlane &t1, const ClusterPlane &t2) {
            return t1.time < t2.time;
        });

        ClusterPlane theCluster;
        theCluster.time = timeReadyToMatch;

        //First ClusterPlane with time that bigger than timeReadyToMatch
        auto it = std::upper_bound(m_clusters_new[dp].begin(), m_clusters_new[dp].end(), theCluster, [](const ClusterPlane &t1, const ClusterPlane &t2) {
            return t1.time < t2.time;
        });

        //Find elements in vector that could still be matched with another cluster
        //since they are close in time to timeReadyToMatch
        while (it != m_clusters_new[dp].end())
        {
            if ((*it).time - timeReadyToMatch > m_config.pDeltaTimeHits)
            {
                break;
            }
            timeReadyToMatch = (*it).time;
            ++it;
        }
        index = std::distance(m_clusters_new[dp].begin(), it);
    }

    //Insert the clusters that are ready to be matched from newClusters into clusters
    m_clusters[dp].insert(m_clusters[dp].end(), std::make_move_iterator(m_clusters_new[dp].begin()), std::make_move_iterator(m_clusters_new[dp].begin() + index));
    //Delete the clusters from newClusters
    m_clusters_new[dp].erase(m_clusters_new[dp].begin(), m_clusters_new[dp].begin() + index);

    return true;
}

void Clusterer::PrintStats()
{
    int totalPlane0 = 0;
    int totalPlane1 = 0;
    int totalDetector = 0;
    bool bothPlanes = false;
    for (auto const &searchDet : m_config.pDets)
    {
        if (m_config.GetAxes(searchDet.first, 0) && m_config.GetAxes(searchDet.first, 1))
        {
            bothPlanes = true;
        }
        std::cout << "\n\n*******************************************************************" << std::endl;
        std::cout << "*********************        Stats detector " << (int)searchDet.first << " ******************" << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        auto dp0 = std::make_pair(searchDet.first, 0);
        auto dp1 = std::make_pair(searchDet.first, 1);
        m_stats.SetLowestCommonTriggerTimestampPlane(dp1, m_stats.lowestCommonTriggerTimestampPlane(dp0));
        m_stats.SetLowestCommonTriggerTimestampDet(searchDet.first, m_stats.lowestCommonTriggerTimestampPlane(dp0));
        m_stats.SetLowestCommonTriggerTimestampPlane(dp0, 100 * m_config.pTriggerPeriod + std::max(m_stats.lowestCommonTriggerTimestampPlane(dp0), m_stats.lowestCommonTriggerTimestampPlane(dp1)));
        AnalyzeClustersPlane(searchDet.first, 0);
        AnalyzeClustersPlane(searchDet.first, 1);
        AnalyzeClustersDetector(searchDet.first);

        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "Plane 0: Max delta time" << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        for (unsigned int n = 0; n <= static_cast<unsigned int>(m_config.pDeltaTimeHits / 50); n++)
        {
            std::cout << n * 50 << "-" << n * 50 + 49 << " ns:  " << m_stats.maxDeltaTime(dp0, n) << std::endl;
        }
        std::cout << "*******************************************************************" << std::endl;

        if (bothPlanes)
        {
            std::cout << "\n*******************************************************************" << std::endl;
            std::cout << "Plane 1: Max delta time" << std::endl;
            std::cout << "*******************************************************************" << std::endl;
            for (unsigned int n = 0; n <= static_cast<int>(m_config.pDeltaTimeHits / 50); n++)
            {
                std::cout << n * 50 << "-" << n * 50 + 49 << " ns:  " << m_stats.maxDeltaTime(dp1, n) << std::endl;
            }
            std::cout << "*******************************************************************" << std::endl;
        }
        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "Plane 0: Max missing strip" << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        for (unsigned int n = 0; n <= m_config.pMissingStripsCluster; n++)
        {
            std::cout << n << ": " << m_stats.maxMissingStrip(dp0, n) << std::endl;
        }
        std::cout << "*******************************************************************" << std::endl;

        if (bothPlanes)
        {
            std::cout << "\n*******************************************************************" << std::endl;
            std::cout << "Plane 1: Max missing strip" << std::endl;
            std::cout << "*******************************************************************" << std::endl;
            for (unsigned int n = 0; n <= m_config.pMissingStripsCluster; n++)
            {
                std::cout << n << ": " << m_stats.maxMissingStrip(dp1, n) << std::endl;
            }
            std::cout << "*******************************************************************" << std::endl;
        }
        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "Plane 0: Span cluster time" << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        for (unsigned int n = 0; n <= static_cast<unsigned int>(m_config.pSpanClusterTime / 50); n++)
        {
            std::cout << n * 50 << "-" << n * 50 + 49 << " ns:  " << m_stats.deltaSpan(dp0, n) << std::endl;
        }
        std::cout << "*******************************************************************" << std::endl;

        if (bothPlanes)
        {
            std::cout << "\n*******************************************************************" << std::endl;
            std::cout << "Plane 1: Span cluster time" << std::endl;
            std::cout << "*******************************************************************" << std::endl;
            for (unsigned int n = 0; n <= static_cast<unsigned int>(m_config.pSpanClusterTime / 50); n++)
            {
                std::cout << n * 50 << "-" << n * 50 + 49 << " ns:  " << m_stats.deltaSpan(dp1, n) << std::endl;
            }
            std::cout << "*******************************************************************" << std::endl;

            std::cout << "\n*******************************************************************" << std::endl;
            std::cout << "Plane 0/plane 1: Max delta plane" << std::endl;
            std::cout << "*******************************************************************" << std::endl;
            for (unsigned int n = 0; n <= static_cast<unsigned int>(m_config.pDeltaTimePlanes / 50); n++)
            {
                std::cout << n * 50 << "-" << n * 50 + 49 << " ns:  " << m_stats.deltaPlane(searchDet.first, n) << std::endl;
            }
            std::cout << "*******************************************************************" << std::endl;

            std::cout << "\n*******************************************************************" << std::endl;
            std::cout << "Plane 0/plane 1: charge ratio" << std::endl;
            std::cout << "*******************************************************************" << std::endl;
            for (unsigned int n = 0; n <= 9; n++)
            {
                std::cout << n * 10 << "-" << n * 10 + 9 << " %:  " << m_stats.chargeRatio(searchDet.first, n, 0) << std::endl;
            }
            std::cout << "100 %   :  " << m_stats.chargeRatio(searchDet.first, 10, 0) << std::endl;

            std::cout << "*******************************************************************" << std::endl;

            std::cout << "\n*******************************************************************" << std::endl;
            std::cout << "Plane 1/plane 0: charge ratio" << std::endl;
            std::cout << "*******************************************************************" << std::endl;
            for (unsigned int n = 0; n <= 9; n++)
            {
                std::cout << n * 10 << "-" << n * 10 + 9 << " %:  " << m_stats.chargeRatio(searchDet.first, n, 1) << std::endl;
            }
            std::cout << "*******************************************************************" << std::endl;
        }

        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "Clusters in plane 0: " << m_clusters_cnt[dp0] << std::endl;
        if (bothPlanes)
        {
            std::cout << "Clusters in plane 1: " << m_clusters_cnt[dp1] << std::endl;
            std::cout << "Common clusters in plane 0 / plane 1: " << m_clusters_detector_cnt[searchDet.first] << std::endl;
        }
        std::cout << "*******************************************************************" << std::endl;
        totalPlane0 += m_clusters_cnt[dp0];
        totalPlane1 += m_clusters_cnt[dp1];
        totalDetector += m_clusters_detector_cnt[searchDet.first];
        bothPlanes = false;
    }

    for (auto const &fec : m_config.pFecs)
    {
        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "Time stamp too large errors fec " << (int)fec << ": " << m_stats.timestampTooLargeError(fec) << std::endl;
        std::cout << "Overflows fec " << (int)fec << ": " << m_stats.overflow(fec) << std::endl;
        std::cout << "Time errors fec " << (int)fec << ": " << m_stats.timeError(fec) << std::endl;
        std::cout << "*******************************************************************" << std::endl;
    }
    std::cout << "\n*******************************************************************" << std::endl;
    std::cout << "Total clusters in plane 0: " << totalPlane0 << std::endl;
    std::cout << "Total clusters in plane 1: " << totalPlane1 << std::endl;
    std::cout << "Total common clusters in plane 0 / plane 1: " << totalDetector << std::endl;
    std::cout << "*******************************************************************" << std::endl;
}

//====================================================================================================================
std::pair<int, int> Clusterer::GetDetectorPlane(std::pair<uint8_t, uint8_t> fecChip)
{
    auto pair = m_config.pFecChip_DetectorPlane.find(fecChip);
    if (pair != end(m_config.pFecChip_DetectorPlane))
    {
        return pair->second;
    }
    return std::make_pair(-1, -1);
    ;
}

//====================================================================================================================
int Clusterer::GetChannel(std::pair<uint8_t, uint8_t> fecChip, int chNo)
{
    auto search = m_config.pOffsets.find(fecChip);
    auto det_plane = m_config.pFecChip_DetectorPlane[fecChip];
    auto flag = m_config.pAxes[det_plane];
    if (search != end(m_config.pOffsets))
    {
        uint32_t ch = chNo + search->second;
        if (flag == 1)
        {
            ch = search->second - chNo;
        }
        return ch;
    }
    else
    {
        return -1;
    }
}
