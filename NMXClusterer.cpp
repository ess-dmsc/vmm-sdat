#include <algorithm>
#include <cmath>
#include "NMXClusterer.h"
#include "Trace.h"

#include <chrono>
#include <functional>

#define UNUSED __attribute__((unused))



auto now = std::chrono::steady_clock::now;

auto timethis(std::function<void()> thunk) -> decltype((now()-now()).count()) {
	auto start = now();
	thunk();
	auto stop = now();
	return (stop - start).count();
}

NMXClusterer::NMXClusterer(int bc, int tac, std::vector<
		std::pair<uint8_t, uint8_t>> xChips, std::vector<
		std::pair<uint8_t, uint8_t>> yChips, std::vector<
		std::pair<uint8_t, uint8_t>> ignoreChips, uint16_t adcThreshold, uint16_t minClusterSize, uint16_t xyClusterSize, uint16_t deltaTimeHits, uint16_t missingStripCluster, uint16_t spanClusterTime, uint16_t deltaTimePlanes, bool analyzeChannels, bool useUTPC, bool createHits) :
		pBC(bc), pTAC(tac), pXChipIDs(xChips), pYChipIDs(yChips), pIgnoreChipIDs(ignoreChips), pADCThreshold(adcThreshold), pMinClusterSize(minClusterSize), pXYClusterSize(xyClusterSize), pDeltaTimeHits(deltaTimeHits), pMissingStripsCluster(missingStripCluster), pSpanClusterTime(spanClusterTime), pDeltaTimePlanes(deltaTimePlanes), pAnalyzeChannels(analyzeChannels), pUseUTPC(useUTPC), pCreateHits(createHits), m_eventNr(0) {
#ifdef USE_ROOT
	m_rootFile = RootFile::GetInstance();
#endif
	pBCTime_ns = 1000 / (int) pBC;
	pTriggerPeriod = 1000 * 4096 / (int) pBC;
	//pTriggerPeriod = 204800;
	m_errors = 0;
	uint8_t idx = 0;
	for (int i = 0; i < pXChipIDs.size(); i++) {
		auto pair = pXChipIDs[i];
		auto fec = std::get < 0 > (pair);
		if (!IsFecInPlane(fec, 0)) {
			p_Fec_Index.emplace(std::make_pair(fec, p_Fec_Index.size()));
			p_Index_Fec.emplace(std::make_pair(p_Index_Fec.size(), fec));
			p_Fec_Plane.push_back(std::make_pair(fec, 0));
		}
	}
	for (int i = 0; i < pYChipIDs.size(); i++) {
		auto pair = pYChipIDs[i];
		auto fec = std::get < 0 > (pair);
		if (!IsFecInPlane(fec, 1)) {
			p_Fec_Index.emplace(std::make_pair(fec, p_Fec_Index.size()));
			p_Index_Fec.emplace(std::make_pair(p_Index_Fec.size(), fec));
			p_Fec_Plane.push_back(std::make_pair(fec, 1));
		}
	}

	for (int i = 0; i < p_Fec_Index.size(); i++) {
		m_oldTriggerTimestamp.push_back(0);
		m_deltaTriggerTimestamp.push_back(0);
		m_stats_overflow.push_back(0);
		m_stats_timeError.push_back(0);
	}
	for (int n = 0; n <= static_cast<int>(pDeltaTimeHits / 50); n++) {
		m_statsX_maxDeltaTime.push_back(0);
		m_statsY_maxDeltaTime.push_back(0);
	}
	for (int n = 0; n <= pMissingStripsCluster; n++) {
		m_statsX_maxMissingStrip.push_back(0);
		m_statsY_maxMissingStrip.push_back(0);
	}
	for (int n = 0; n <= static_cast<int>(pSpanClusterTime / 50); n++) {
		m_statsX_deltaSpan.push_back(0);
		m_statsY_deltaSpan.push_back(0);
	}
	for (int n = 0; n <= static_cast<int>(pDeltaTimePlanes / 50); n++) {
		m_statsXY_deltaPlane.push_back(0);
	}
}

NMXClusterer::~NMXClusterer() {

}

//====================================================================================================================
int NMXClusterer::AnalyzeHits(double srsTimestamp, uint8_t fecID, uint8_t vmmID, uint16_t chNo, uint16_t bcid, uint16_t tdc, uint16_t adc, bool overThresholdFlag, float chipTime) {
	//if(m_eventNr > 5) return -1;

	auto search = std::find(pIgnoreChipIDs.begin(), pIgnoreChipIDs.end(), std::make_pair(fecID, vmmID));
	if (search != pIgnoreChipIDs.end()) {
		return true;
	}

	m_deltaTriggerTimestamp[p_Fec_Index[fecID]] = 0;
	if (srsTimestamp < m_oldTriggerTimestamp[p_Fec_Index[fecID]]) {
		if (m_oldTriggerTimestamp[p_Fec_Index[fecID]] - srsTimestamp > 0x1FFFFFFFFFF) {
			m_stats_overflow[p_Fec_Index[fecID]]++;
			DTRACE(DEB,"\n*********************************** OVERFLOW  fecIndex %d , fecID %d, m_lineNr %d, eventNr  %d, "
				"srsTimestamp %llu, old srsTimestamp %llu\n", p_Fec_Index[fecID], fecID, m_lineNr, m_eventNr, static_cast<uint64_t>(srsTimestamp), static_cast<uint64_t>(m_oldTriggerTimestamp[p_Fec_Index[fecID]]));
			srsTimestamp += 0x3FFFFFFFFFF;
		} else {
			m_stats_timeError[p_Fec_Index[fecID]]++;
			DTRACE(DEB,"\n*********************************** TIME ERROR  fecIndex %d , fecID %d, m_lineNr %d, eventNr  %d, "
			"srsTimestamp %llu, old srsTimestamp %llu\n", p_Fec_Index[fecID], fecID, m_lineNr, m_eventNr, static_cast<uint64_t>(srsTimestamp), static_cast<uint64_t>(m_oldTriggerTimestamp[p_Fec_Index[fecID]]));
		}
	}

	double remainder = std::fmod(m_deltaTriggerTimestamp[p_Fec_Index[fecID]], pTriggerPeriod);
	if (remainder > 0) {
		uint64_t offset = m_deltaTriggerTimestamp[p_Fec_Index[fecID]] / pTriggerPeriod;
		m_deltaTriggerTimestamp[p_Fec_Index[fecID]] = m_deltaTriggerTimestamp[p_Fec_Index[fecID]] - offset * pTriggerPeriod - remainder + (offset + 1) * pTriggerPeriod;
		srsTimestamp = m_oldTriggerTimestamp[p_Fec_Index[fecID]] + m_deltaTriggerTimestamp[p_Fec_Index[fecID]];
		DTRACE(DEB,"\n******* SRS timestamp wrong increment: fec %d,vmmId %d, chNo %d, line %d, "
				"trigger period %d, offset %llu, remainder %llu,  delta time %llu, "
				"new time %llu, old time %llu\n", fecID, vmmID, chNo, m_lineNr, pTriggerPeriod, offset, static_cast<uint64_t>(remainder), static_cast<uint64_t>(m_deltaTriggerTimestamp[p_Fec_Index[fecID]]), static_cast<uint64_t>(srsTimestamp), static_cast<uint64_t>(m_oldTriggerTimestamp[p_Fec_Index[fecID]]));
	}

	bool newEvent = false;
	if (srsTimestamp > m_oldTriggerTimestamp[p_Fec_Index[fecID]]) {
		m_deltaTriggerTimestamp[p_Fec_Index[fecID]] = srsTimestamp - m_oldTriggerTimestamp[p_Fec_Index[fecID]];
		newEvent = true;
	}

	if (newEvent) {
		m_eventNr++;

		for (int n = 0; n < p_Index_Fec.size(); n++) {
			int factor = 10;
			if (IsFecInPlane(p_Index_Fec[n], 0)) {
				if (m_oldTriggerTimestamp[n] > factor * pTriggerPeriod && m_lowestCommonTriggerTimestamp_x < m_oldTriggerTimestamp[n]) {
					m_lowestCommonTriggerTimestamp_x = m_oldTriggerTimestamp[n] - factor * pTriggerPeriod;
				}
			}
			if (IsFecInPlane(p_Index_Fec[n], 1)) {
				if (m_oldTriggerTimestamp[n] > factor * pTriggerPeriod && m_lowestCommonTriggerTimestamp_y < m_oldTriggerTimestamp[n]) {
					m_lowestCommonTriggerTimestamp_y = m_oldTriggerTimestamp[n] - factor * pTriggerPeriod;
				}
			}
		}
		m_lowestCommonTriggerTimestamp_xy = std::min(m_lowestCommonTriggerTimestamp_x, m_lowestCommonTriggerTimestamp_y);

		AnalyzeClusters(m_hitsX, m_hitsX_new, m_lowestCommonTriggerTimestamp_x, pDeltaTimeHits, "x");

		AnalyzeClusters(m_hitsY, m_hitsY_new, m_lowestCommonTriggerTimestamp_y, pDeltaTimeHits, "y");

		AnalyzeClustersXY(m_lowestCommonTriggerTimestamp_xy);

	}

	// Plane 0: x
	// plane 1: y
	int x = -1;
	int y = -1;
	int planeID = GetPlaneID(std::make_pair(fecID, vmmID));
	if (planeID == 0) {
		x = GetChannel(pXChipIDs, std::make_pair(fecID, vmmID), chNo);

	} else if (planeID == 1) {
		y = GetChannel(pYChipIDs, std::make_pair(fecID, vmmID), chNo);
	}

	// TDC has reduced resolution due to most significant bit problem of current
	// sources (like ADC)
	//int theTDC =  8 * (tdc / 8) + 4;
	//theTDC = tdc;

	//uint64_t bcTime = pBCTime_ns * (bcid + 1);
	//uint64_t tdcTime = tdc * pTAC / 256;

	if (pCreateHits) {
		HitNMX theHit;
		theHit.eventNr = m_eventNr;
		theHit.fecID = fecID;
		theHit.vmmID = vmmID;
		theHit.triggerTimestamp = srsTimestamp;
		theHit.chNo = chNo;
		theHit.position = 0;

		theHit.bcid = bcid;
		theHit.tdc = tdc;
		theHit.adc = adc;
		theHit.overThresholdFlag = overThresholdFlag;
		theHit.chipTime = chipTime;
		theHit.totalTime = srsTimestamp + chipTime;
		if (x >= 0) {
			theHit.position = (uint16_t) x;
			m_rootFile->SaveHitsX(std::move(theHit));
		} else if (y >= 0) {
			theHit.position = (uint16_t) y;
			m_rootFile->SaveHitsY(std::move(theHit));
		}
	}

	//uint64_t theChiptime = (bcTime - tdcTime);
	uint64_t totalTime = srsTimestamp + (int) chipTime;
	StoreHits(x, y, adc, bcid, totalTime, overThresholdFlag);

	if (newEvent) {
		DTRACE(DEB, "\neventNr  %d\n", m_eventNr);
		//DTRACE(DEB, "fecID  %d\n", fecID);
	}
	if (m_deltaTriggerTimestamp[p_Fec_Index[fecID]] > 0) {
		DTRACE(DEB, "\tTriggerTimestamp %llu [ns]\n", static_cast<uint64_t>(srsTimestamp));
		DTRACE(DEB, "\tTime since last trigger %f us (%.4f kHz)\n", m_deltaTriggerTimestamp[p_Fec_Index[fecID]] * 0.001, (double )(1000000 / m_deltaTriggerTimestamp[p_Fec_Index[fecID]]));
	}

	if (m_oldFecID != fecID || newEvent) {
		DTRACE(DEB, "\tfecID  %d\n", fecID);

	}
	if (m_oldVmmID != vmmID || newEvent) {
		DTRACE(DEB, "\tvmmID  %d\n", vmmID);
	}
	if (planeID == 0) {
		DTRACE(DEB, "\t\tx-channel %d (chNo  %d) - overThresholdFlag %d\n", x, chNo, overThresholdFlag);
	} else if (planeID == 1) {
		DTRACE(DEB, "\t\ty-channel %d (chNo  %d) - overThresholdFlag %d\n", y, chNo, overThresholdFlag);
	} else {
		DTRACE(DEB, "\t\tPlane for vmmID %d not defined!\n", vmmID);
	}
	DTRACE(DEB, "\t\t\tbcid %d, tdc %d, adc %d\n", bcid, tdc, adc);
	DTRACE(DEB, "\t\t\ttotal time %llu, chip time %f ns\n", totalTime, chipTime);

	m_oldTriggerTimestamp[p_Fec_Index[fecID]] = srsTimestamp;

	m_oldBcID = bcid;
	m_oldVmmID = vmmID;
	m_oldFecID = fecID;
	m_lineNr++;
	return 0;
}

//====================================================================================================================
void NMXClusterer::StoreHits(int x, int y, uint16_t adc, uint16_t bcid, double chipTime, bool overThresholdFlag) {
	if ((adc >= pADCThreshold || overThresholdFlag)) {
		if (x > -1) {
			m_hitsX_new.emplace_back(chipTime, (uint16_t) x, adc);
		} else if (y > -1) {
			m_hitsY_new.emplace_back(chipTime, (uint16_t) y, adc);
		}

	}
}

//====================================================================================================================
int NMXClusterer::ClusterByTime(HitContainer &hits, uint16_t dTime, uint16_t missingStrips, uint16_t spanTime, string coordinate) {

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
			clusterCount += ClusterByStrip(cluster, missingStrips, spanTime, coordinate, maxDeltaTime);
			cluster.clear();
			maxDeltaTime = 0;
		}
		cluster.emplace_back(strip1, time1, adc1);
		stripCount++;
	}

	if (stripCount > 0) {
		clusterCount += ClusterByStrip(cluster, missingStrips, spanTime, coordinate, maxDeltaTime);
	}
	return clusterCount;
}

//====================================================================================================================
int NMXClusterer::ClusterByStrip(ClusterContainer &cluster, uint16_t missingStrips, uint16_t spanTime, string coordinate, uint16_t maxDeltaTime) {
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
	std::map<double, float, std::greater<float>> mClusters;
	ClusterContainer stripClusters;
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
			DTRACE(DEB, "\n%s cluster:\n", coordinate.c_str());
		}

// Add members of a cluster, if it is either the beginning of a cluster,
// or if strip gap and time span is correct
		if (stripCount == 0 || (std::abs(strip1 - strip2) > 0 && std::abs(strip1 - strip2) - 1 <= missingStrips && time1 - startTime <= spanTime)) {
			DTRACE(DEB, "\tstrip %d, time %llu, adc %d:\n", strip1, (uint64_t )time1, adc1);

			if (time1 > largestTime) {
				largestTime = time1;
				clusterPosition = strip1;
			}
			if (time1 < startTime) {
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
			mClusters.emplace(std::make_pair(time1, strip1));
			stripCount++;
		}
// Stop clustering if gap between strips is too large or time span too long
		else if (std::abs(strip1 - strip2) - 1 > missingStrips || largestTime - startTime > spanTime) {
// Valid cluster
			if (stripCount < pMinClusterSize || totalADC == 0) {
				DTRACE(DEB, "******** INVALID ********\n\n");
			} else {
				centerOfGravity = (centerOfGravity / (float) totalADC);
				centerOfTime = (centerOfTime / (float) totalADC);
				StoreClusters(vStrips, vTimes, clusterPosition, largestTime, centerOfGravity, centerOfTime, stripCount, totalADC, coordinate, maxDeltaTime, maxMissingStrip, spanCluster);
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
		StoreClusters(vStrips, vTimes, clusterPosition, largestTime, centerOfGravity, centerOfTime, stripCount, totalADC, coordinate, maxDeltaTime, maxMissingStrip, spanCluster);
		clusterCount++;
	}
	return clusterCount;
}

//====================================================================================================================
void NMXClusterer::StoreClusters(std::vector<float>& strips, std::vector<double>& times, float clusterPosition, double clusterTime, float centerOfCharge, double centerOfTime, uint16_t clusterSize, uint32_t clusterADC, string coordinate, uint16_t maxDeltaTime, uint16_t maxMissingStrip, uint16_t deltaSpan) {
	ClusterNMX theCluster;
	m_cluster_id++;
	DTRACE(DEB, "id %d\n", m_cluster_id);
	theCluster.id = m_cluster_id;
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
	theCluster.deltaSpan = deltaSpan;
	theCluster.strips = std::move(strips);
	theCluster.times = std::move(times);
	if (coordinate == "x" && clusterPosition > -1.0) {
		m_clusterX_new.emplace_back(std::move(theCluster));
		m_statsX_maxDeltaTime[((int) maxDeltaTime / 50)]++;
		m_statsX_maxMissingStrip[(int) (maxMissingStrip)]++;
		m_statsX_deltaSpan[(int) (deltaSpan / 50)]++;
	}
	if (coordinate == "y" && clusterPosition > -1.0) {
		m_clusterY_new.emplace_back(std::move(theCluster));
		m_statsY_maxDeltaTime[(int) (maxDeltaTime / 50)]++;
		m_statsY_maxMissingStrip[(int) (maxMissingStrip)]++;
		m_statsY_deltaSpan[(int) (deltaSpan / 50)]++;
	}

}

//====================================================================================================================
void NMXClusterer::MatchClustersXY(uint16_t dPlane) {

	for (auto & nx : m_clusterX) {
		double ctx = nx.centerOfTime;
		double tx = nx.time;

		double minDelta = 99999999;
		double deltaT = 0;
		double deltaCT = 0;
		ClusterVector::iterator it = end(m_clusterY);

		double ty = 0;
		double cty = 0;

		for (ClusterVector::iterator ny = begin(m_clusterY);
				ny != end(m_clusterY); ++ny) {
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

		if (it != end(m_clusterY)) {
			nx.clusterXAndY = true;
			(*it).clusterXAndY = true;

			CommonClusterNMX theCommonCluster;
			m_clusterXY_id++;
			theCommonCluster.id = m_clusterXY_id;
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
			m_statsXY_deltaPlane[std::abs(theCommonCluster.deltaPlane / 50)]++;
			theCommonCluster.maxDeltaTimeX = nx.maxDeltaTime;
			theCommonCluster.maxDeltaTimeY = (*it).maxDeltaTime;
			theCommonCluster.maxMissingStripX = nx.maxMissingStrip;
			theCommonCluster.maxMissingStripY = (*it).maxMissingStrip;

			theCommonCluster.stripsX = nx.strips;
			theCommonCluster.timesX = nx.times;
			theCommonCluster.stripsY = (*it).strips;
			theCommonCluster.timesY = (*it).times;
			DTRACE(DEB, "\ncommon cluster x/y: %d/%d", theCommonCluster.idX, theCommonCluster.idY);
			DTRACE(DEB, "\tpos x/pos y: %f/%f", theCommonCluster.positionX, theCommonCluster.positionY);
			DTRACE(DEB, "\ttime x/time y: : %llu/%llu", (uint64_t )theCommonCluster.timeX, (uint64_t )theCommonCluster.timeY);
			DTRACE(DEB, "\tadc x/adc y: %u/%u", theCommonCluster.adcX, theCommonCluster.adcY);
			DTRACE(DEB, "\tsize x/size y: %u/%u", theCommonCluster.sizeX, theCommonCluster.sizeY);
			DTRACE(DEB, "\tdelta time planes: %d", theCommonCluster.deltaPlane);
			m_clusterXY.emplace_back(std::move(theCommonCluster));

		}

	}

}

void NMXClusterer::AnalyzeClusters(HitContainer& hits, HitContainer& newHits, double timeReadyToCluster, uint16_t correctionTime, std::string coordinate) {

	if (ChooseHitsToBeClustered(hits, newHits, timeReadyToCluster, pDeltaTimeHits, coordinate) == false) {
		return;
	}

	int cntX = ClusterByTime(hits, pDeltaTimeHits, pMissingStripsCluster, pSpanClusterTime, coordinate);

	DTRACE(DEB, "%d cluster in %s\n", cntX, coordinate.c_str());

	if (!hits.empty()) {
		hits.clear();
	}
}

void NMXClusterer::AnalyzeClustersXY(double timeReadyToCluster) {

	if (ChooseClustersToBeMatched(m_clusterX, m_clusterX_new, timeReadyToCluster, pDeltaTimePlanes, "x") == false) {
		return;
	}

	if (ChooseClustersToBeMatched(m_clusterY, m_clusterY_new, timeReadyToCluster, pDeltaTimePlanes, "y") == false) {
		return;
	}

	MatchClustersXY(pDeltaTimePlanes);

	m_clusterCnt_XY += m_clusterXY.size();
	m_clusterCnt_X += m_clusterX.size();
	m_clusterCnt_Y += m_clusterY.size();

	m_rootFile->SaveClustersX(std::move(m_clusterX));
	m_rootFile->SaveClustersY(std::move(m_clusterY));
	m_rootFile->SaveClustersXY(std::move(m_clusterXY));
	m_rootFile->FillTree();
	int m_clusterCnt_X = 0;
	int m_clusterCnt_Y = 0;
	int m_clusterCnt_XY = 0;

	m_clusterX.clear();
	m_clusterY.clear();
	m_clusterXY.clear();

}

//====================================================================================================================
bool NMXClusterer::ChooseHitsToBeClustered(HitContainer &data, HitContainer &newData, double timeReady, uint16_t correctionTime, std::string coordinate) {
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

	/*
	 std::cout << " " << std::endl;
	 for (auto x : data) {
	 std::cout << "\t\t\t\t\t\thits " << coordinate << " " << (uint64_t) timeReady << " " << (uint64_t)(std::get < 1 > (x)) << " " << (uint64_t)(std::get < 0 > (x)) << std::endl;
	 }
	 std::cout << " " << std::endl;
	 for (auto x : newData) {
	 std::cout << "\t\t\t\t\t\tnewHits " << coordinate << " " << (uint64_t) timeReady << " " << (uint64_t)(std::get < 1 > (x)) << " " << (uint64_t)(std::get < 0 > (x)) << std::endl;
	 }
	 */
	return true;
}

bool NMXClusterer::ChooseClustersToBeMatched(ClusterVector &data, ClusterVector &newData, double timeReady, uint16_t correctionTime, std::string coordinate) {
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

	/*
	 #ifdef TRC_LEVEL
	 std::cout << " " << std::endl;
	 for (auto x : data) {
	 std::cout << "\t\t\t\t\t\tclusters " << coordinate << " " << x.id << " " << (uint64_t) timeReady << " " << x.position << " " << (uint64_t) x.time << std::endl;
	 }
	 std::cout << " " << std::endl;
	 for (auto x : newData) {
	 std::cout << "\t\t\t\t\t\tnewClusters " << coordinate << " " << x.id << " " << (uint64_t) timeReady << " " << x.position << " " << (uint64_t) x.time << std::endl;
	 }
	 #endif
	 */
	return true;

}

void NMXClusterer::PrintStats() {

	m_lowestCommonTriggerTimestamp_x = std::max(m_lowestCommonTriggerTimestamp_x, m_lowestCommonTriggerTimestamp_y);
	m_lowestCommonTriggerTimestamp_x += 100 * pTriggerPeriod;
	m_lowestCommonTriggerTimestamp_y = m_lowestCommonTriggerTimestamp_x;
	m_lowestCommonTriggerTimestamp_xy = m_lowestCommonTriggerTimestamp_x;
	AnalyzeClusters(m_hitsX, m_hitsX_new, m_lowestCommonTriggerTimestamp_x, pDeltaTimeHits, "x");
	AnalyzeClusters(m_hitsY, m_hitsY_new, m_lowestCommonTriggerTimestamp_y, pDeltaTimeHits, "y");
	AnalyzeClustersXY(m_lowestCommonTriggerTimestamp_xy);

	std::cout << "*******************************************************************" << std::endl;
	std::cout << "*********************        Stats        *************************" << std::endl;
	std::cout << "*******************************************************************" << std::endl;

	std::cout << "\n*******************************************************************" << std::endl;
	std::cout << "x: Max delta time" << std::endl;
	std::cout << "*******************************************************************" << std::endl;
	for (int n = 0; n <= static_cast<int>(pDeltaTimeHits / 50); n++) {
		std::cout << "[" << n * 50 << "-" << (n + 1) * 50 << "[ ns:  " << m_statsX_maxDeltaTime[n] << std::endl;
	}
	std::cout << "*******************************************************************" << std::endl;

	std::cout << "\n*******************************************************************" << std::endl;
	std::cout << "x: Max missing strip" << std::endl;
	std::cout << "*******************************************************************" << std::endl;
	for (int n = 0; n <= pMissingStripsCluster; n++) {
		std::cout << n << ": " << m_statsX_maxMissingStrip[n] << std::endl;
	}
	std::cout << "*******************************************************************" << std::endl;

	std::cout << "\n*******************************************************************" << std::endl;
	std::cout << "x: Span cluster time" << std::endl;
	for (int n = 0; n <= static_cast<int>(pSpanClusterTime / 50); n++) {
		std::cout << n * 50 << "-" << (n + 1) * 50 << "[ ns:  " << m_statsX_deltaSpan[n] << std::endl;
	}
	std::cout << "*******************************************************************" << std::endl;

	std::cout << "\n*******************************************************************" << std::endl;
	std::cout << "y: Max delta time" << std::endl;
	for (int n = 0; n <= static_cast<int>(pDeltaTimeHits / 50); n++) {
		std::cout << n * 50 << "-" << (n + 1) * 50 << "[ ns:  " << m_statsY_maxDeltaTime[n] << std::endl;
	}
	std::cout << "*******************************************************************" << std::endl;

	std::cout << "\n*******************************************************************" << std::endl;
	std::cout << "y: Max missing strips" << std::endl;
	std::cout << "*******************************************************************" << std::endl;
	for (int n = 0; n <= pMissingStripsCluster; n++) {
		std::cout << n << ": " << m_statsY_maxMissingStrip[n] << std::endl;
	}
	std::cout << "*******************************************************************" << std::endl;

	std::cout << "\n*******************************************************************" << std::endl;
	std::cout << "y: Span cluster time" << std::endl;
	std::cout << "*******************************************************************" << std::endl;
	for (int n = 0; n < static_cast<int>(pSpanClusterTime / 50); n++) {
		std::cout << n * 50 << "-" << (n + 1) * 50 << "[ ns:  " << m_statsY_deltaSpan[n] << std::endl;
	}
	std::cout << "*******************************************************************" << std::endl;

	std::cout << "\n*******************************************************************" << std::endl;
	std::cout << "x/y: Max delta place" << std::endl;
	std::cout << "*******************************************************************" << std::endl;
	for (int n = 0; n < static_cast<int>(pDeltaTimePlanes / 50); n++) {
		std::cout << n * 50 << "-" << (n + 1) * 50 << "[ ns:  " << m_statsXY_deltaPlane[n] << std::endl;
	}
	std::cout << "*******************************************************************" << std::endl;

	std::cout << "\n*******************************************************************" << std::endl;
	std::cout << "Clusters in x: " << m_clusterCnt_X << std::endl;
	std::cout << "Clusters in y: " << m_clusterCnt_Y << std::endl;
	std::cout << "Common clusters in x/y: " << m_clusterCnt_XY << std::endl;
	std::cout << "*******************************************************************" << std::endl;

	std::cout << "\n*******************************************************************" << std::endl;
	for (int n = 0; n < p_Fec_Index.size(); n++) {
		std::cout << "Overflows fec " << (int)p_Index_Fec[n] << ": " << m_stats_overflow[n] << std::endl;
	}
	for (int n = 0; n < p_Fec_Index.size(); n++) {
			std::cout << "Time errors fec " << (int)p_Index_Fec[n] << ": " << m_stats_timeError[n] << std::endl;
	}
	std::cout << "*******************************************************************" << std::endl;
}

bool NMXClusterer::IsFecInPlane(uint8_t fecId, uint8_t planeId) {

	auto search = std::find(p_Fec_Plane.begin(), p_Fec_Plane.end(), std::make_pair(fecId, planeId));
	if (search != p_Fec_Plane.end()) {
		return true;
	}

	return false;
}

//====================================================================================================================
int NMXClusterer::GetPlaneID(std::pair<uint8_t, uint8_t> chipID) {
	auto chip = std::find(begin(pXChipIDs), end(pXChipIDs), chipID);
	if (chip != end(pXChipIDs)) {
		return 0;
	} else {
		auto chip = std::find(begin(pYChipIDs), end(pYChipIDs), chipID);
		if (chip != end(pYChipIDs)) {
			return 1;
		} else {
			return -1;
		}
	}
}

//====================================================================================================================
uint32_t NMXClusterer::GetChannel(std::vector<std::pair<uint8_t, uint8_t>>& chipIDs, std::pair<
		uint8_t, uint8_t> chipID, int channelID) {
	auto chip = std::find(begin(chipIDs), end(chipIDs), chipID);
	if (chip != end(chipIDs)) {
		uint32_t ch = channelID + (chip - begin(chipIDs)) * 64;
		/*
		 if (ch % 2 == 0) {
		 ch += 1;
		 } else {
		 ch -= 1;
		 }
		 */
		return ch;
	} else {
		return -1;
	}
}

