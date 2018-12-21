#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>

#include <thread>   // sleep_for
#include <future>   // async
#include <atomic>   // async
#include <queue>   // async
#include "TH1F.h"
#include "TF1.h"
#include "TFitResult.h"

#define USE_ROOT

#ifdef USE_ROOT
#include "RootFile.h"
#else
#include "DataStructures.h"
#endif

class NMXClusterer {
public:
	NMXClusterer(int bc, int tac, std::vector<std::pair<uint8_t, uint8_t>> xChips, std::vector<
			std::pair<uint8_t, uint8_t>> yChips, std::vector<
			std::pair<uint8_t, uint8_t>> ignoreChips, uint16_t adcThreshold, uint16_t minClusterSize,
			uint16_t xyClusterSize, uint16_t deltaTimeHits, uint16_t missingStripCluster,
			uint16_t spanClusterTime, uint16_t deltaTimePlanes, bool analyzeChannels, bool useUTPC, bool createHits);

	~NMXClusterer();

	// Analyzing and storing the hits
	int AnalyzeHits(double srsTimestamp, uint8_t fecID, uint8_t vmmID, uint16_t chNo, uint16_t bcid, uint16_t tdc, uint16_t adc, bool overThresholdFlag, float chipTime);
	void StoreHits(int x, int y, uint16_t adc, uint16_t bcid, double chipTime, bool overThresholdFlag);

	// Analyzing and storing the clusters in x or y
	void AnalyzeClusters(HitContainer& hits, HitContainer& newHits, double timeReadyToCluster, uint16_t correctionTime, std::string coordinate);


	//Select hits that are ready to be clustered in time
	bool ChooseHitsToBeClustered(HitContainer &data, HitContainer &newData, double timeReady, uint16_t correctionTime,std::string coordinate);

	bool ChooseClustersToBeMatched(ClusterVector &data, ClusterVector &newData, double timeReady, uint16_t correctionTime, std::string coordinate);

	// Matching the clusters that are common between x and y
	void AnalyzeClustersXY(double timeStamp);

	int ClusterByTime(HitContainer& hits, uint16_t dTime, uint16_t missingStrips, uint16_t spanTime, string coordinate);
	int ClusterByStrip(ClusterContainer &cluster, uint16_t missingStrips, uint16_t spanTime, string coordinate, uint16_t maxDeltaTime);

	void StoreClusters(std::vector<float>& strips, std::vector<double>& times,
			float clusterPosition, double clusterTime, float centerOfCharge,
			double centerOfTime, uint16_t clusterSize, uint32_t clusterADC, string coordinate,
			uint16_t maxDeltaTime, uint16_t maxDeltaStrip, uint16_t deltaSpan);

	void MatchClustersXY(uint16_t dPlane);

	void PrintStats();

	bool IsFecInPlane(uint8_t fecId, uint8_t planeId);

// Helper methods that map channels to strips
	int GetPlaneID(std::pair<uint8_t, uint8_t> chipID);
	uint32_t GetChannel(std::vector<std::pair<uint8_t, uint8_t>>& chipIDs, std::pair<
			uint8_t, uint8_t> chipID, int channelID);

	int getNumClustersX() {
		return m_clusterX.size();
	}
	;
	int getNumClustersY() {
		return m_clusterY.size();
	}
	;
	int getNumClustersXY() {
		return m_clusterXY.size();
	}
	;

#ifdef USE_ROOT
	void createRootFile(string fileName);
#endif
private:
	uint16_t pBC;
	uint16_t pTAC;
	uint32_t pBCTime_ns = 0;
	uint32_t pTriggerPeriod = 0;
	std::vector<std::pair<uint8_t, uint8_t>> pXChipIDs;
	std::vector<std::pair<uint8_t, uint8_t>> pYChipIDs;
	std::vector<std::pair<uint8_t, uint8_t>> pIgnoreChipIDs;

	std::map<uint8_t, uint8_t> p_Fec_Index;
	std::map<uint8_t, uint8_t> p_Index_Fec;
	std::vector<std::pair<uint8_t, uint8_t>> p_Fec_Plane;

	uint16_t pADCThreshold;
	uint16_t pMinClusterSize;
	uint16_t pXYClusterSize;
	uint16_t pDeltaTimeHits;
	uint16_t pMissingStripsCluster;
	uint16_t pSpanClusterTime;
	uint16_t pDeltaTimePlanes;
	bool pAnalyzeChannels;
	bool pUseUTPC = true;
	bool pCreateHits;

	int m_lineNr = 1;
	int m_eventNr = 0;
	double m_lowestCommonTriggerTimestamp_x = 0;
	double m_lowestCommonTriggerTimestamp_y = 0;
	double m_lowestCommonTriggerTimestamp_xy = 0;

	std::vector<double> m_deltaTriggerTimestamp;
	std::vector<double> m_oldTriggerTimestamp;

	//double m_timeStamp_ms = 0;
	uint16_t m_oldBcID = 0;
	uint8_t m_oldVmmID = 0;
	uint8_t m_oldFecID = 0;

	HitContainer m_hitsX;
	HitContainer m_hitsY;
	HitContainer m_hitsX_new;
	HitContainer m_hitsY_new;

	CommonClusterVector m_clusterXY;
	ClusterVector m_clusterX_new;
	ClusterVector m_clusterY_new;
	ClusterVector m_clusterX;
	ClusterVector m_clusterY;
	int m_cluster_id = 0;
	int m_clusterXY_id = 0;
	int m_clusterCnt_X = 0;
	int m_clusterCnt_Y = 0;
	int m_clusterCnt_XY = 0;
	RootFile* m_rootFile;

	int m_errors = 0;

	std::vector<int> m_statsX_maxDeltaTime;
	std::vector<int> m_statsX_maxMissingStrip;
	std::vector<int> m_statsX_deltaSpan;
	std::vector<int> m_statsY_maxDeltaTime;
	std::vector<int> m_statsY_maxMissingStrip;
	std::vector<int> m_statsY_deltaSpan;
	std::vector<int> m_statsXY_deltaPlane;

	std::vector<int> m_stats_overflow;
	std::vector<int> m_stats_timeError;

};

