#include <algorithm>
#include <cmath>
#include "NMXClusterer.h"
#include "Trace.h"

#include <chrono>
#include <functional>

#define UNUSED __attribute__((unused))

//#undef TRC_LEVEL
//#define TRC_LEVEL TRC_L_DEB


auto now = std::chrono::steady_clock::now;

auto timethis(std::function<void()> thunk) -> decltype((now()-now()).count()) {
    auto start = now();
    thunk();
    auto stop = now();
    return (stop - start).count();
}

NMXClusterer::NMXClusterer(std::string rootFileName, int bc, int tac, std::vector<
        std::tuple<uint8_t,uint8_t, uint8_t>> xChips, std::vector<
        std::tuple<uint8_t,uint8_t, uint8_t>> yChips, uint16_t adcThreshold, uint16_t minClusterSize,
        uint16_t xyClusterSize, uint16_t deltaTimeHits, uint16_t missingStripCluster, uint16_t spanClusterTime, uint16_t deltaTimePlanes,
                           bool analyzeChannels, bool useUTPC, bool createHits) :
                                                            pBC(bc), pTAC(tac), pXChipIDs(xChips), pYChipIDs(yChips), pADCThreshold(adcThreshold), pMinClusterSize(minClusterSize), pXYClusterSize(xyClusterSize), pDeltaTimeHits(deltaTimeHits), pMissingStripsCluster(missingStripCluster), pSpanClusterTime(spanClusterTime), pDeltaTimePlanes(deltaTimePlanes), pAnalyzeChannels(analyzeChannels), pUseUTPC(useUTPC), pCreateHits(createHits), m_eventNr(0) {

    pBCTime_ns = 1000 / (int) pBC;
    pTriggerPeriod = 1000 * 4096 / (int) pBC;
    //pTriggerPeriod = 204800;
    m_errors = 0;
    uint8_t idx = 0;

    for (int i = 0; i < pXChipIDs.size(); i++) {
        uint8_t plane = 0;
        auto tuple = pXChipIDs[i];
        auto det = std::get < 0 > (tuple);
        auto fec = std::get < 1 > (tuple);
        auto chip = std::get < 2 > (tuple);
        m_lowestCommonTriggerTimestamp_plane[std::make_pair(det, plane)] = 0;
        m_lowestCommonTriggerTimestamp_det[det] = 0;
        auto searchDet = pDets.find(det);
        if (searchDet == pDets.end()) {
            pDets.emplace(det,pDets.size());
        }
        auto searchFec = m_deltaTriggerTimestamp.find(fec);
        if (searchFec == m_deltaTriggerTimestamp.end()) {
            m_deltaTriggerTimestamp.emplace(std::make_pair(fec, 0));
            m_oldTriggerTimestamp.emplace(std::make_pair(fec, 0));
            m_stats_overflow.emplace(std::make_pair(fec, 0));
            m_stats_timeError.emplace(std::make_pair(fec, 0));
        }
        bool found = false;
        //Search whether there is a new det/plane/fec combination
        for(auto const &searchDetPlaneFec: pDetectorPlane_Fec)
        {
            if(searchDetPlaneFec.first == std::make_pair(det, plane) && searchDetPlaneFec.second == fec)
            {
                found = true;
                break;
            }
        }
        if(found == false)
        {
            pDetectorPlane_Fec.emplace(std::make_pair(std::make_pair(det, plane),fec));
        }

        //Search whether there is a new fec/chip combination
        auto searchFecChip = pFecChip_DetectorPlane.find(std::make_pair(fec, chip));
        if (searchFecChip == pFecChip_DetectorPlane.end())
        {
            //Add the new fec/chip pair to the list
            pFecChip_DetectorPlane.emplace(std::make_pair(std::make_pair(fec, chip), std::make_pair(det, plane)));
            //Search for det/plane pairs
            auto searchDetPlane = p_DetPlane_idx.find(std::make_pair(det, plane));
            uint32_t offset = 0;
            if (searchDetPlane == p_DetPlane_idx.end()) {
                //Add det/plane pair to the list and set index
                p_DetPlane_idx.emplace(std::make_pair(std::make_pair(det, plane), offset));
            }
            else
            {
                //Increment det/plane index
                offset = searchDetPlane->second + 1;
                p_DetPlane_idx[std::make_pair(det, 0)] = offset;
            }
            //Set offset for the new fec/chip combination
            pOffsets.emplace(std::make_pair(std::make_pair(fec, chip), offset*64));
        }
    }

    for (int i = 0; i < pYChipIDs.size(); i++) {
        uint8_t plane = 1;
        auto tuple = pYChipIDs[i];
        auto det = std::get < 0 > (tuple);
        auto fec = std::get < 1 > (tuple);
        auto chip = std::get < 2 > (tuple);
        m_lowestCommonTriggerTimestamp_plane[std::make_pair(det, plane)] = 0;
        m_lowestCommonTriggerTimestamp_det[det] = 0;

        auto searchDet = pDets.find(det);
        if (searchDet == pDets.end()) {
            pDets.emplace(det,pDets.size());
        }
        auto searchFec = m_deltaTriggerTimestamp.find(fec);
        if (searchFec == m_deltaTriggerTimestamp.end()) {
            m_deltaTriggerTimestamp.emplace(std::make_pair(fec, 0));
            m_oldTriggerTimestamp.emplace(std::make_pair(fec, 0));
            m_stats_overflow.emplace(std::make_pair(fec, 0));
            m_stats_timeError.emplace(std::make_pair(fec, 0));
        }
        bool found = false;
        //Search whether there is a new det/plane/fec combination
        for(auto const &searchDetPlaneFec: pDetectorPlane_Fec)
        {
            if(searchDetPlaneFec.first == std::make_pair(det, plane) && searchDetPlaneFec.second == fec)
            {
                found = true;
                break;
            }
        }
        if(found == false)
        {
            pDetectorPlane_Fec.emplace(std::make_pair(std::make_pair(det, plane),fec));
        }
        //Search whether there is a new fec/chip combination
        auto searchFecChip = pFecChip_DetectorPlane.find(std::make_pair(fec, chip));
        if (searchFecChip == pFecChip_DetectorPlane.end())
        {
            //Add the new fec/chip pair to the list
            pFecChip_DetectorPlane.emplace(std::make_pair(std::make_pair(fec, chip), std::make_pair(det, plane)));
            //Search for det/plane pairs
            auto searchDetPlane = p_DetPlane_idx.find(std::make_pair(det, plane));
            uint32_t offset = 0;
            if (searchDetPlane == p_DetPlane_idx.end()) {
                //Add det/plane pair to the list and set index
                p_DetPlane_idx.emplace(std::make_pair(std::make_pair(det, plane), offset));
            }
            else
            {
                //Increment det/plane index
                offset = searchDetPlane->second + 1;
                p_DetPlane_idx[std::make_pair(det, 1)] = offset;
            }
            //Set offset for the new fec/chip combination
            pOffsets.emplace(std::make_pair(std::make_pair(fec, chip), offset*64));
        }
    }


    for(auto const & det: pDets)
    {
        auto x = std::make_pair(det.first,0);
        auto y = std::make_pair(det.first,1);
        for (unsigned int n = 0; n <= static_cast<unsigned int>(pDeltaTimeHits / 50); n++) {
            m_stats_maxDeltaTime[x].push_back(0);
            m_stats_maxDeltaTime[y].push_back(0);
        }
        for (unsigned int n = 0; n <= pMissingStripsCluster; n++) {
            m_stats_maxMissingStrip[x].push_back(0);
            m_stats_maxMissingStrip[y].push_back(0);
        }
        for (unsigned int n = 0; n <= static_cast<unsigned int>(pSpanClusterTime / 50); n++) {
            m_stats_deltaSpan[x].push_back(0);
            m_stats_deltaSpan[y].push_back(0);
        }
        for (unsigned int n = 0; n <=static_cast<unsigned int>(pDeltaTimePlanes / 50); n++) {

            m_stats_deltaPlane[det.first].push_back(0);
        }
    }


#ifdef USE_ROOT
    m_rootFile = RootFile::GetInstance(rootFileName, pDets, p_DetPlane_idx, analyzeChannels, createHits);
#endif

}

NMXClusterer::~NMXClusterer() {
    RootFile::Dispose();
}

//====================================================================================================================
bool NMXClusterer::AnalyzeHits(double srsTimestamp, uint8_t fecId, uint8_t vmmId, uint16_t chNo, uint16_t bcid, uint16_t tdc, uint16_t adc, bool overThresholdFlag, float chipTime) {

    auto searchFecChip = pFecChip_DetectorPlane.find(std::make_pair(fecId, vmmId));
    if (searchFecChip == pFecChip_DetectorPlane.end()) {
        DTRACE(DEB, "\t\tDetector or Plane not defined for FEC %d and vmmId %d not defined!\n", (int)fecId, (int)vmmId);
        return true;
    }
    auto det_plane = GetDetectorPlane(std::make_pair(fecId, vmmId));
    auto det = std::get < 0 > (det_plane);
    auto plane = std::get < 1 > (det_plane);

    // Plane 0: x
    // plane 1: y
    int pos = GetChannel(std::make_pair(fecId, vmmId), chNo);
    if (pos == -1) {
        DTRACE(DEB, "\t\tChannel not defined for FEC %d and vmmId %d not defined!\n", (int)fecId, (int)vmmId);
        return true;
    }
    m_deltaTriggerTimestamp[fecId] = 0;
    if (srsTimestamp < m_oldTriggerTimestamp[fecId]) {
        if (m_oldTriggerTimestamp[fecId] - srsTimestamp > 0x1FFFFFFFFFF) {
            m_stats_overflow[fecId]++;
                        DTRACE(DEB,"\n*********************************** OVERFLOW  fecIndex %d , fecId %d, m_lineNr %d, eventNr  %d, "
                        "srsTimestamp %llu, old srsTimestamp %llu\n", fecId, fecId, m_lineNr, m_eventNr, static_cast<uint64_t>(srsTimestamp), static_cast<uint64_t>(m_oldTriggerTimestamp[fecId]));
                        srsTimestamp += 0x3FFFFFFFFFF;
        } else {
            m_stats_timeError[fecId]++;
                        DTRACE(DEB,"\n*********************************** TIME ERROR  fecIndex %d , fecId %d, m_lineNr %d, eventNr  %d, "
                        "srsTimestamp %llu, old srsTimestamp %llu\n", fecId, fecId, m_lineNr, m_eventNr, static_cast<uint64_t>(srsTimestamp), static_cast<uint64_t>(m_oldTriggerTimestamp[fecId]));
        }
    }

    double remainder = std::fmod(m_deltaTriggerTimestamp[fecId], pTriggerPeriod);
    if (remainder > 0) {
        uint64_t offset = m_deltaTriggerTimestamp[fecId] / pTriggerPeriod;
        m_deltaTriggerTimestamp[fecId] = m_deltaTriggerTimestamp[fecId] - offset * pTriggerPeriod - remainder + (offset + 1) * pTriggerPeriod;
        srsTimestamp = m_oldTriggerTimestamp[fecId] + m_deltaTriggerTimestamp[fecId];
                DTRACE(DEB,"\n******* SRS timestamp wrong increment: fec %d,vmmId %d, chNo %d, line %d, "
				"trigger period %d, offset %llu, remainder %llu,  delta time %llu, "
                                "new time %llu, old time %llu\n", fecId, vmmId, chNo, m_lineNr, pTriggerPeriod, offset, static_cast<uint64_t>(remainder), static_cast<uint64_t>(m_deltaTriggerTimestamp[fecId]), static_cast<uint64_t>(srsTimestamp), static_cast<uint64_t>(m_oldTriggerTimestamp[fecId]));
    }

    bool newEvent = false;
    if (srsTimestamp > m_oldTriggerTimestamp[fecId]) {
        m_deltaTriggerTimestamp[fecId] = srsTimestamp - m_oldTriggerTimestamp[fecId];
        newEvent = true;
    }



    if (newEvent) {
        m_eventNr++;

        int factor = 4;
        for(auto const &searchDetPlaneFec: pDetectorPlane_Fec)
        {
            auto fec = searchDetPlaneFec.second;
            auto det_plane = searchDetPlaneFec.first;
            auto det = std::get < 0 > (det_plane);
            auto plane = std::get < 1 > (det_plane);
            if (m_oldTriggerTimestamp[fec] > factor * pTriggerPeriod && m_lowestCommonTriggerTimestamp_plane[det_plane] < m_oldTriggerTimestamp[fec]) {
                m_lowestCommonTriggerTimestamp_plane[det_plane] = m_oldTriggerTimestamp[fec] - factor * pTriggerPeriod;
            }
            m_lowestCommonTriggerTimestamp_det[det] = std::min(m_lowestCommonTriggerTimestamp_plane[std::make_pair(det,0)], m_lowestCommonTriggerTimestamp_plane[std::make_pair(det,1)]);
        }
        for(auto const &searchDet: pDets)
        {
            AnalyzeClusters(searchDet.first, 0, m_hits[std::make_pair(searchDet.first,0)], m_hits_new[std::make_pair(searchDet.first,0)], m_lowestCommonTriggerTimestamp_plane[std::make_pair(searchDet.first,0)], pDeltaTimeHits);
            AnalyzeClusters(searchDet.first, 1, m_hits[std::make_pair(searchDet.first,1)], m_hits_new[std::make_pair(searchDet.first,1)], m_lowestCommonTriggerTimestamp_plane[std::make_pair(searchDet.first,1)], pDeltaTimeHits);
            AnalyzeClustersXY(searchDet.first, m_lowestCommonTriggerTimestamp_det[searchDet.first]);
        }

    }



    // TDC has reduced resolution due to most significant bit problem of current
    // sources (like ADC)
    //int theTDC =  8 * (tdc / 8) + 4;
    //theTDC = tdc;

    //uint64_t bcTime = pBCTime_ns * (bcid + 1);
    //uint64_t tdcTime = tdc * pTAC / 256;
    m_lineNr++;
    if (pCreateHits) {

        HitNMX theHit;
        theHit.id = m_lineNr;
        theHit.eventNr = m_eventNr;
        theHit.detId = det;
        theHit.planeId = plane;
        theHit.fecId = fecId;
        theHit.vmmId = vmmId;
        theHit.triggerTimestamp = srsTimestamp;
        theHit.chNo = chNo;
        theHit.position = pos;

        theHit.bcid = bcid;
        theHit.tdc = tdc;
        theHit.adc = adc;
        theHit.overThresholdFlag = overThresholdFlag;
        theHit.chipTime = chipTime;
        theHit.totalTime = srsTimestamp + chipTime;
        theHit.position = (uint16_t) pos;
        m_rootFile->AddHits(std::move(theHit));
    }

    //uint64_t theChiptime = (bcTime - tdcTime);
    uint64_t totalTime = srsTimestamp + (int) chipTime;
    StoreHits(det, plane, pos, adc, bcid, totalTime, overThresholdFlag);

    if (newEvent) {
        DTRACE(DEB, "\neventNr  %d\n", m_eventNr);
        //DTRACE(DEB, "fecId  %d\n", fecId);
    }
    if (m_deltaTriggerTimestamp[fecId] > 0) {
        DTRACE(DEB, "\tTriggerTimestamp %llu [ns]\n", static_cast<uint64_t>(srsTimestamp));
        DTRACE(DEB, "\tTime since last trigger %f us (%.4f kHz)\n", m_deltaTriggerTimestamp[fecId] * 0.001, (double )(1000000 / m_deltaTriggerTimestamp[fecId]));
    }

    if (m_oldFecId != fecId || newEvent) {
        DTRACE(DEB, "\tfecId  %d\n", fecId);

    }
    if (m_oldVmmId != vmmId || newEvent) {
        DTRACE(DEB, "\tDetector %d, plane %d, vmmId  %d\n",(int)det, (int)plane, vmmId);
    }
    DTRACE(DEB, "\t\tChannel %d (chNo  %d) - overThresholdFlag %d\n",  pos, chNo, (int)overThresholdFlag);
    DTRACE(DEB, "\t\t\tbcid %d, tdc %d, adc %d\n", bcid, tdc, adc);
    DTRACE(DEB, "\t\t\ttotal time %llu, chip time %f ns\n", totalTime, chipTime);

    m_oldTriggerTimestamp[fecId] = srsTimestamp;

    m_oldBcId = bcid;
    m_oldVmmId = vmmId;
    m_oldFecId = fecId;

    return true;
}

//====================================================================================================================
void NMXClusterer::StoreHits(uint8_t det, uint8_t plane, int pos, uint16_t adc, uint16_t bcid, double chipTime, bool overThresholdFlag) {
    if ((adc >= pADCThreshold || overThresholdFlag)) {
        m_hits_new[std::make_pair(det, plane)].emplace_back(chipTime, (uint16_t) pos, adc);
    }
}

//====================================================================================================================
int NMXClusterer::ClusterByTime(uint8_t det, uint8_t plane, HitContainer &hits, uint16_t dTime, uint16_t missingStrips, uint16_t spanTime) {

    ClusterContainer cluster;
    uint16_t maxDeltaTime = 0;
    int clusterCount = 0;
    int stripCount = 0;
    double time1 = 0, time2 = 0;
    uint32_t adc1 = 0;
    uint16_t strip1 = 0;

    for (auto& itHits : hits) {
        time2 = time1;

        time1 = (double) std::get < 0 > (itHits);
        strip1 = std::get < 1 > (itHits);
        adc1 = std::get < 2 > (itHits);

        if (abs(time1 - time2) <= dTime && stripCount > 0 && maxDeltaTime < abs(time1 - time2)) {
            maxDeltaTime = (time1 - time2);
        }

        if (abs(time1 - time2) > dTime && stripCount > 0) {
            clusterCount += ClusterByStrip(det, plane, cluster, missingStrips, spanTime, maxDeltaTime);
            cluster.clear();
            maxDeltaTime = 0;
        }
        cluster.emplace_back(strip1, time1, adc1);
        stripCount++;
    }

    if (stripCount > 0) {
        clusterCount += ClusterByStrip(det, plane, cluster, missingStrips, spanTime, maxDeltaTime);
    }
    return clusterCount;
}

//====================================================================================================================
int NMXClusterer::ClusterByStrip(uint8_t det, uint8_t plane, ClusterContainer &cluster, uint16_t missingStrips, uint16_t spanTime, uint16_t maxDeltaTime) {
    int maxMissingStrip = 0;
    uint16_t spanCluster = 0;

    double startTime = 0;
    double largestTime = 0;
    float clusterPosition = -1;
    float centerOfGravity = 0;
    double centerOfTime = 0;
    long int totalADC = 0;
    //long int totalADCSquared = 0;

    double time1 = 0;
    int adc1 = 0;
    int strip1 = 0;
    int strip2 = 0;
    int stripCount = 0;
    int clusterCount = 0;
    std::vector<float> vStrips;
    std::vector<double> vTimes;
    std::sort(begin(cluster), end(cluster), [](const ClusterTuple &t1, const ClusterTuple &t2)
              {
                  return std::get<0>(t1) < std::get<0>(t2);
              });
    for (auto& itCluster : cluster) {
        strip2 = strip1;
        strip1 = std::get < 0 > (itCluster);
        time1 = std::get < 1 > (itCluster);
        adc1 = std::get < 2 > (itCluster);

        // At beginning of cluster, set start time of cluster
        if (stripCount == 0) {
            maxMissingStrip = 0;
            startTime = time1;
            largestTime = time1;
            DTRACE(DEB, "\nDetector %d, plane %d cluster:\n", (int)det, (int)plane);
        }
        // Add members of a cluster, if it is either the beginning of a cluster,
        // or if strip gap and time span is correct
        if (stripCount == 0 || (std::abs(strip1 - strip2) > 0 && std::abs(strip1 - strip2) - 1 <= missingStrips && time1 - startTime <= spanTime && largestTime - time1 <= spanTime)) {
            DTRACE(DEB, "\tstrip %d, time %llu, adc %d:\n", strip1, (uint64_t )time1, adc1);

            if (time1 >= largestTime) {
                largestTime = time1;
                clusterPosition = strip1;
            }
            if (time1 <= startTime) {
                startTime = time1;
            }
            if (stripCount > 0 && maxMissingStrip < std::abs(strip1 - strip2) - 1) {
                maxMissingStrip = std::abs(strip1 - strip2) - 1;
            }
            spanCluster = (largestTime - startTime);
            totalADC += adc1;
            centerOfGravity += strip1 * adc1;
            centerOfTime += time1 * adc1;
            vStrips.emplace_back(strip1);
            vTimes.emplace_back(time1);
            stripCount++;
        }
        // Stop clustering if gap between strips is too large or time span too long
        else if (std::abs(strip1 - strip2) - 1 > missingStrips || time1 - startTime > spanTime || largestTime - time1 > spanTime) {
            // Valid cluster
            if (stripCount < pMinClusterSize || totalADC == 0) {
                DTRACE(DEB, "******** INVALID ********\n\n");
            } else {
                spanCluster = (largestTime - startTime);
                centerOfGravity = (centerOfGravity / (float) totalADC);
                centerOfTime = (centerOfTime / (float) totalADC);
                StoreClusters(det, plane, vStrips, vTimes, clusterPosition, largestTime, centerOfGravity, centerOfTime, stripCount, totalADC, maxDeltaTime, maxMissingStrip, spanCluster);
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
            strip1 = 0;
        }
    }
    // At the end of the clustering, check again if there is a last valid cluster
    if (stripCount < pMinClusterSize || totalADC == 0) {
        DTRACE(DEB, "******** INVALID ********\n\n");
    } else {
        spanCluster = (largestTime - startTime);
        centerOfGravity = (centerOfGravity / (float) totalADC);
        centerOfTime = (centerOfTime / (float) totalADC);
        StoreClusters(det, plane, vStrips, vTimes, clusterPosition, largestTime, centerOfGravity, centerOfTime,
                      stripCount, totalADC, maxDeltaTime, maxMissingStrip, spanCluster);

        clusterCount++;
    }
    return clusterCount;
}

//====================================================================================================================
void NMXClusterer::StoreClusters(uint8_t det, uint8_t plane, std::vector<float>& strips, std::vector<double>& times, float clusterPosition, double clusterTime, float centerOfCharge, double centerOfTime, uint16_t clusterSize, uint32_t clusterADC, uint16_t maxDeltaTime, uint16_t maxMissingStrip, uint16_t deltaSpan) {
    ClusterNMX theCluster;
    m_cluster_id++;

    DTRACE(DEB, "Cluster id %d\n", m_cluster_id);
    theCluster.id = static_cast<uint32_t>(m_cluster_id);
    theCluster.detId = det;
    theCluster.planeId = plane;
    theCluster.size = clusterSize;
    theCluster.adc = clusterADC;
    theCluster.time = clusterTime;
    theCluster.position = clusterPosition;
    if (pUseUTPC) {
        theCluster.time = clusterTime;
        theCluster.position = clusterPosition;

    } else {
        theCluster.time = centerOfTime;
        theCluster.position = centerOfCharge;
    }
    theCluster.centerOfTime = centerOfTime;
    theCluster.centerOfCharge = centerOfCharge;
    theCluster.utpcTime = clusterTime;
    theCluster.utpcPosition = clusterPosition;
    theCluster.clusterXAndY = false;
    theCluster.maxDeltaTime = maxDeltaTime;
    theCluster.maxMissingStrip = maxMissingStrip;
    theCluster.spanCluster = deltaSpan;
    theCluster.strips = std::move(strips);
    theCluster.times = std::move(times);
    if (theCluster.position > -1.0) {

        m_stats_maxDeltaTime[std::make_pair(det, plane)][(unsigned int)(maxDeltaTime / 50)]++;
        m_stats_maxMissingStrip[std::make_pair(det, plane)][(unsigned int) (maxMissingStrip)]++;
        m_stats_deltaSpan[std::make_pair(det, plane)][(unsigned int) (deltaSpan / 50)]++;
        m_clusters_new[std::make_pair(det, plane)].emplace_back(std::move(theCluster));
    }
    else
    {
        std::cout << (int)theCluster.id << " " << (int) theCluster.detId << " " <<  (int)theCluster.planeId << " " << theCluster.position << std::endl;
    }
}

//====================================================================================================================
void NMXClusterer::MatchClustersXY(uint8_t det, uint16_t dPlane) {
    for (auto & nx : m_clusters[std::make_pair(det, 0)]) {
        double ctx = nx.centerOfTime;
        double tx = nx.time;

        double minDelta = 99999999;
        double deltaT = 0;
        double deltaCT = 0;
        ClusterVector::iterator it = end(m_clusters[std::make_pair(det, 1)]);

        double ty = 0;
        double cty = 0;

        for (ClusterVector::iterator ny = begin(m_clusters[std::make_pair(det, 1)]);
             ny != end(m_clusters[std::make_pair(det, 1)]); ++ny) {
            if ((*ny).clusterXAndY == false) {
                cty = (*ny).centerOfTime;
                ty = (*ny).time;
                deltaT = std::abs(ty - tx);
                deltaCT = std::abs(cty - ctx);
                uint64_t chargeRatio = (*ny).adc / nx.adc;
                //if (deltaCT < minDelta && deltaCT <= dPlane && chargeRatio >= 0.8
                //&& chargeRatio < 1.25 && (nx.size + (*ny).size >= pXYClusterSize)) {
                if (deltaCT < minDelta && deltaCT <= dPlane && deltaT <= dPlane && (nx.size + (*ny).size >= pXYClusterSize)) {
                    minDelta = deltaCT;
                    it = ny;
                }
            }
        }

        if (it != end(m_clusters[std::make_pair(det, 1)])) {
            nx.clusterXAndY = true;
            (*it).clusterXAndY = true;
            CommonClusterNMX theCommonCluster;
            m_clusterXY_id++;
            theCommonCluster.id = m_clusterXY_id;
            theCommonCluster.detId = det;
            theCommonCluster.idX = nx.id;
            theCommonCluster.idY = (*it).id;
            theCommonCluster.sizeX = nx.size;
            theCommonCluster.sizeY = (*it).size;
            theCommonCluster.adcX = nx.adc;
            theCommonCluster.adcY = (*it).adc;
            if (pUseUTPC) {
                theCommonCluster.positionX = nx.position;
                theCommonCluster.positionY = (*it).position;
                theCommonCluster.timeX = nx.time;
                theCommonCluster.timeY = (*it).time;
            } else {
                theCommonCluster.positionX = nx.centerOfCharge;
                theCommonCluster.positionY = (*it).centerOfCharge;
                theCommonCluster.timeX = nx.centerOfTime;
                theCommonCluster.timeY = (*it).centerOfTime;
            }
            theCommonCluster.deltaPlane = theCommonCluster.timeX - theCommonCluster.timeY;

            m_stats_deltaPlane[det][std::abs(theCommonCluster.deltaPlane/50)]++;

            theCommonCluster.maxDeltaTimeX = nx.maxDeltaTime;
            theCommonCluster.maxDeltaTimeY = (*it).maxDeltaTime;
            theCommonCluster.maxMissingStripX = nx.maxMissingStrip;
            theCommonCluster.maxMissingStripY = (*it).maxMissingStrip;
            theCommonCluster.spanClusterX = nx.spanCluster;
            theCommonCluster.spanClusterY = (*it).spanCluster;
            theCommonCluster.stripsX = nx.strips;
            theCommonCluster.timesX = nx.times;
            theCommonCluster.stripsY = (*it).strips;
            theCommonCluster.timesY = (*it).times;

            DTRACE(DEB, "\ncommon cluster det %d x/y: %d/%d", (int)det, theCommonCluster.idX, theCommonCluster.idY);
            DTRACE(DEB, "\tpos x/pos y: %f/%f", theCommonCluster.positionX, theCommonCluster.positionY);
            DTRACE(DEB, "\ttime x/time y: : %llu/%llu", (uint64_t )theCommonCluster.timeX, (uint64_t )theCommonCluster.timeY);
            DTRACE(DEB, "\tadc x/adc y: %u/%u", theCommonCluster.adcX, theCommonCluster.adcY);
            DTRACE(DEB, "\tsize x/size y: %u/%u", theCommonCluster.sizeX, theCommonCluster.sizeY);
            DTRACE(DEB, "\tdelta time planes: %d", theCommonCluster.deltaPlane);
            m_clustersXY[det].emplace_back(std::move(theCommonCluster));

        }

    }

}

void NMXClusterer::AnalyzeClusters(uint8_t det, uint8_t plane, HitContainer& hits, HitContainer& newHits, double timeReadyToCluster, uint16_t correctionTime) {

    if (ChooseHitsToBeClustered(hits, newHits, timeReadyToCluster, pDeltaTimeHits) == false) {
        return;
    }

    int cnt = ClusterByTime(det, plane, hits, pDeltaTimeHits, pMissingStripsCluster, pSpanClusterTime);

    DTRACE(DEB, "%d cluster in detector %d plane %d\n", cnt, (int)det, (int)plane);

    if (!hits.empty()) {
        hits.clear();
    }
}

void NMXClusterer::AnalyzeClustersXY(uint8_t det, double timeReadyToCluster) {

    if (ChooseClustersToBeMatched(m_clusters[std::make_pair(det, 0)], m_clusters_new[std::make_pair(det, 0)], timeReadyToCluster, pDeltaTimePlanes) == false) {
        return;
    }

    if (ChooseClustersToBeMatched(m_clusters[std::make_pair(det, 1)], m_clusters_new[std::make_pair(det, 1)], timeReadyToCluster, pDeltaTimePlanes) == false) {
        return;
    }

    MatchClustersXY(det, pDeltaTimePlanes);

    m_clusterCnt_XY[det] += m_clustersXY[det].size();
    m_clusterCnt[std::make_pair(det, 0)] += m_clusters[std::make_pair(det, 0)].size();
    m_clusterCnt[std::make_pair(det, 1)] += m_clusters[std::make_pair(det, 1)].size();

    m_rootFile->SaveHits();
    m_rootFile->SaveClusters(std::move(m_clusters[std::make_pair(det, 0)]));
    m_rootFile->SaveClusters(std::move(m_clusters[std::make_pair(det, 1)]));
    m_rootFile->SaveClustersXY(std::move(m_clustersXY[det]));

    m_clusters[std::make_pair(det, 0)].clear();
    m_clusters[std::make_pair(det, 1)].clear();
    m_clustersXY[det].clear();

}

//====================================================================================================================
bool NMXClusterer::ChooseHitsToBeClustered(HitContainer &data, HitContainer &newData, double timeReady, uint16_t correctionTime) {
    //Nothing to cluster, newHits vector empty
    if (newData.empty()) {
        return false;
    }

    auto theMin = std::min_element(newData.begin(), newData.end(), [](const HitTuple &t1, const HitTuple &t2)
                                   {
                                       return std::get<0>(t1) < std::get<0>(t2);
                                   });

    //Nothing to cluster, tuples in newHits vector too recent
    if (std::get < 0 > (*theMin) > timeReady) {

        //(smallest timestamp larger than timeReady)
        //Will be clustered later
        return false;
    }

    //Sort vector newHits
    std::sort(begin(newData), end(newData), [](const HitTuple &t1, const HitTuple &t2)
              {
                  return std::get<0>(t1) < std::get<0>(t2);
              });

    //First tuple with timestamp larger than timeReady
    auto it = std::upper_bound(newData.begin(), newData.end(), std::make_tuple(timeReady, 0, 0), [](const HitTuple &t1, const HitTuple &t2)
                               {
                                   return std::get<0>(t1) < std::get<0>(t2);
                               });

    //Find elements in vector that could still be part of a cluster,
    //since they are close in time to timeReadyToCluster
    while (it != newData.end()) {
        if (std::get < 0 > (*it) - timeReady > correctionTime) {
            break;
        }
        timeReady = std::get < 0 > (*it);
        //std::cout << "XOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXO " << std::get < 1 > (*it) << " " << (uint64_t) std::get < 0 > (*it) << std::endl;
        ++it;
    }

    //std::cout << "still in cluster " << std::get < 0 > (*it) << "\n" << std::endl;

    int index = std::distance(newData.begin(), it);
    //Insert the data that is ready to be clustered from newHits into hits
    data.insert(data.end(), std::make_move_iterator(newData.begin()), std::make_move_iterator(newData.begin() + index));
    //Delete the data from newHits
    newData.erase(newData.begin(), newData.begin() + index);


    return true;
}

bool NMXClusterer::ChooseClustersToBeMatched(ClusterVector &data, ClusterVector &newData, double timeReady, uint16_t correctionTime) {
    //Nothing to match, newClusters vector empty
    if (newData.empty()) {
        return false;
    }

    auto theMin = std::min_element(newData.begin(), newData.end(), [](const ClusterNMX &t1, const ClusterNMX &t2)
                                   {
                                       return t1.time < t2.time;
                                   });

    //Nothing to cluster, clusters in newClusters vector too recent
    if ((*theMin).time > timeReady) {

        //(smallest time larger than timeReadyToMatch)
        //Will be matched later
        return false;
    }

    //Sort vector newClusters based on time
    std::sort(begin(newData), end(newData), [](const ClusterNMX &t1, const ClusterNMX &t2)
              {
                  return t1.time < t2.time;
              });

    ClusterNMX theCluster;
    theCluster.time = timeReady;

    //First ClusterNMX with time that bigger than timeReadyToMatch
    auto it = std::upper_bound(newData.begin(), newData.end(), theCluster, [](const ClusterNMX &t1, const ClusterNMX &t2)
                               {
                                   return t1.time < t2.time;
                               });

    //Find elements in vector that could still be matched with another cluster
    //since they are close in time to timeReadyToCluster
    while (it != newData.end()) {
        if ((*it).time - timeReady > correctionTime) {
            break;
        }
        timeReady = (*it).time;
        //std::cout << "XOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXOXO " << (*it).position << " " << (*it).time << std::endl;
        ++it;
    }

    int index = std::distance(newData.begin(), it);

    //Insert the clusters that are ready to be matched from newClusters into clusters
    data.insert(data.end(), std::make_move_iterator(newData.begin()), std::make_move_iterator(newData.begin() + index));
    //Delete the clusters from newClusters
    newData.erase(newData.begin(), newData.begin() + index);


    return true;

}

void NMXClusterer::PrintStats() {
    int totalX = 0;
    int totalY = 0;
    int totalXY = 0;
    for(auto const &searchDet: pDets)
    {
        std::cout << "\n\n*******************************************************************" << std::endl;
        std::cout << "*********************        Stats detector " << (int)searchDet.first << "******************" << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        auto x = std::make_pair(searchDet.first,0);
        auto y = std::make_pair(searchDet.first,1);
        m_lowestCommonTriggerTimestamp_plane[x] = std::max(m_lowestCommonTriggerTimestamp_plane[x], m_lowestCommonTriggerTimestamp_plane[y]);
        m_lowestCommonTriggerTimestamp_plane[x] += 100 * pTriggerPeriod;
        m_lowestCommonTriggerTimestamp_plane[y] = m_lowestCommonTriggerTimestamp_plane[x];
        m_lowestCommonTriggerTimestamp_det[searchDet.first] = m_lowestCommonTriggerTimestamp_plane[x];

        AnalyzeClusters(searchDet.first, 0, m_hits[x], m_hits_new[x], m_lowestCommonTriggerTimestamp_plane[x], pDeltaTimeHits);
        AnalyzeClusters(searchDet.first, 1, m_hits[y], m_hits_new[y], m_lowestCommonTriggerTimestamp_plane[y], pDeltaTimeHits);

        AnalyzeClustersXY(searchDet.first, m_lowestCommonTriggerTimestamp_det[searchDet.first]);



        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "x: Max delta time" << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        for (unsigned int n = 0; n <= static_cast<unsigned int>(pDeltaTimeHits / 50); n++) {
            std::cout << "[" << n * 50 << "-" << n*50+49  << " ns:  " << m_stats_maxDeltaTime[x][n] << std::endl;
        }
        std::cout << "*******************************************************************" << std::endl;

        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "y: Max delta time" << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        for (unsigned int n = 0; n <= static_cast<int>(pDeltaTimeHits / 50); n++) {
            std::cout << "[" << n * 50 << "-" << n* 50+49 << " ns:  " << m_stats_maxDeltaTime[y][n] << std::endl;
        }
        std::cout << "*******************************************************************" << std::endl;

        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "x: Max missing strip" << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        for (unsigned int n = 0; n <= pMissingStripsCluster; n++) {
            std::cout << n << ": " << m_stats_maxMissingStrip[x][n] << std::endl;
        }
        std::cout << "*******************************************************************" << std::endl;

        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "y: Max missing strip" << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        for (unsigned int n = 0; n <= pMissingStripsCluster; n++) {
            std::cout << n << ": " << m_stats_maxMissingStrip[y][n] << std::endl;
        }
        std::cout << "*******************************************************************" << std::endl;

        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "x: Span cluster time" << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        for (unsigned int n = 0; n <= static_cast<unsigned int>(pSpanClusterTime / 50); n++) {
            std::cout << n * 50 << "-" << n* 50+49 << " ns:  " << m_stats_deltaSpan[x][n] << std::endl;
        }
        std::cout << "*******************************************************************" << std::endl;

        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "y: Span cluster time" << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        for (unsigned int n = 0; n <= static_cast<unsigned int>(pSpanClusterTime / 50); n++) {
            std::cout << n * 50 << "-" << n* 50+49 << " ns:  " << m_stats_deltaSpan[y][n] << std::endl;
        }
        std::cout << "*******************************************************************" << std::endl;

        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "x/y: Max delta place" << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        for (unsigned int n = 0; n <= static_cast<unsigned int>(pDeltaTimePlanes / 50); n++) {
            std::cout << n * 50 << "-" << n* 50+49 << " ns:  " << m_stats_deltaPlane[searchDet.first][n] << std::endl;
        }
        std::cout << "*******************************************************************" << std::endl;


        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "Clusters in x: " << m_clusterCnt[x] << std::endl;
        std::cout << "Clusters in y: " << m_clusterCnt[y] << std::endl;
        std::cout << "Common clusters in x/y: " << m_clusterCnt_XY[searchDet.first] << std::endl;
        std::cout << "*******************************************************************" << std::endl;
        totalX+=m_clusterCnt[x] ;
        totalY+=m_clusterCnt[y] ;
        totalXY+=m_clusterCnt_XY[searchDet.first] ;
    }

    for (auto const &fec: m_stats_overflow) {
        std::cout << "\n*******************************************************************" << std::endl;
        std::cout << "Overflows fec " << (int)fec.first << ": " << m_stats_overflow[fec.first] << std::endl;
        std::cout << "Time errors fec " << (int)fec.first << ": " << m_stats_timeError[fec.first] << std::endl;
        std::cout << "*******************************************************************" << std::endl;
    }
    std::cout << "\n*******************************************************************" << std::endl;
    std::cout << "Total clusters in x: " << totalX << std::endl;
    std::cout << "Toltal clusters in y: " << totalY << std::endl;
    std::cout << "Total common clusters in x/y: " << totalXY << std::endl;
    std::cout << "*******************************************************************" << std::endl;
}


//====================================================================================================================
std::pair<int, int> NMXClusterer::GetDetectorPlane(std::pair<uint8_t, uint8_t> fecChip) {
    auto pair = pFecChip_DetectorPlane.find(fecChip);
    if (pair != end(pFecChip_DetectorPlane)) {
        return pair->second;
    }
    return std::make_pair(-1,-1);;
}


//====================================================================================================================
int NMXClusterer::GetChannel(std::pair<uint8_t, uint8_t> fecChip, int chNo) {
    auto search = pOffsets.find(fecChip);
    if (search != end(pOffsets)) {
        uint32_t ch = chNo + search->second;
        if(fecChip.second == 6 || fecChip.second ==7 || (fecChip.first == 2 && (fecChip.second == 4 || fecChip.second == 5)))
        {
            return ch;
        }
        //Chips for adapter


        if (ch % 2 == 0) {
            ch += 1;
        } else {
            ch -= 1;
        }

        return ch;
    } else {
        return -1;
    }
}

