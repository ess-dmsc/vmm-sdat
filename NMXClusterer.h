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



#define USE_ROOT

#ifdef USE_ROOT
#include "RootFile.h"
#else
#include "DataStructures.h"
#endif

class NMXClusterer {
public:
        NMXClusterer(std::string rootFileName, int bc, int tac, std::vector<std::tuple<uint8_t,uint8_t, uint8_t>> xChips, std::vector<std::tuple<uint8_t,uint8_t, uint8_t>> yChips, uint16_t adcThreshold, uint16_t minClusterSize,
			uint16_t xyClusterSize, uint16_t deltaTimeHits, uint16_t missingStripCluster,
			uint16_t spanClusterTime, uint16_t deltaTimePlanes, bool analyzeChannels, bool useUTPC, bool createHits);

	~NMXClusterer();

	// Analyzing and storing the hits
        bool AnalyzeHits(double srsTimestamp, uint8_t fecId, uint8_t vmmId, uint16_t chNo, uint16_t bcid, uint16_t tdc, uint16_t adc, bool overThresholdFlag, float chipTime);
        void StoreHits(uint8_t det,uint8_t plane, int pos, uint16_t adc, uint16_t bcid, double chipTime, bool overThresholdFlag);

        // Analyzing and storing the clusters in x or y
        void AnalyzeClusters(uint8_t det, uint8_t plane, HitContainer& hits, HitContainer& newHits, double timeReadyToCluster, uint16_t correctionTime);


	//Select hits that are ready to be clustered in time
        bool ChooseHitsToBeClustered(HitContainer &data, HitContainer &newData, double timeReady, uint16_t correctionTime);

        bool ChooseClustersToBeMatched(ClusterVector &data, ClusterVector &newData, double timeReady, uint16_t correctionTime);

        // Matching the clusters that are common between x and y
        void AnalyzeClustersXY(uint8_t det, double timeStamp);

        int ClusterByTime(uint8_t det,uint8_t plane,HitContainer& hits, uint16_t dTime, uint16_t missingStrips, uint16_t spanTime);
        int ClusterByStrip(uint8_t det,uint8_t plane,ClusterContainer &cluster, uint16_t missingStrips, uint16_t spanTime, uint16_t maxDeltaTime);

        void StoreClusters(uint8_t det, uint8_t plane, std::vector<float>& strips, std::vector<double>& times,
			float clusterPosition, double clusterTime, float centerOfCharge,
                        double centerOfTime, uint16_t clusterSize, uint32_t clusterADC, uint16_t maxDeltaTime, uint16_t maxDeltaStrip, uint16_t deltaSpan);

        void MatchClustersXY(uint8_t det,uint16_t dPlane);

	void PrintStats();

        //bool IsFecChipInDetectorPlane(uint8_t fecId, uint8_t det, uint8_t chip, uint8_t planeId);

        std::pair<int, int> GetDetectorPlane(std::pair<uint8_t, uint8_t> fecChip);
        // Helper methods that map channels to strips
        int GetChannel(std::pair<uint8_t, uint8_t> fecChip, int chNo);



#ifdef USE_ROOT
	void createRootFile(string fileName);
#endif
private:
	uint16_t pBC;
	uint16_t pTAC;
	uint32_t pBCTime_ns = 0;
        uint32_t pTriggerPeriod = 0;
        std::vector<std::tuple<uint8_t,uint8_t, uint8_t>> pXChipIDs;
        std::vector<std::tuple<uint8_t,uint8_t, uint8_t>> pYChipIDs;

        std::map<std::pair<uint8_t,uint8_t>, std::pair<uint8_t,uint8_t>> pFecChip_DetectorPlane;
        std::multimap<std::pair<uint8_t,uint8_t>, uint8_t> pDetectorPlane_Fec;
        std::map<std::pair<uint8_t,uint8_t>, uint32_t> pOffsets;
        std::map<std::pair<uint8_t, uint8_t>, uint32_t>  p_DetPlane_idx;
        std::map<uint8_t, uint8_t> pDets;
        std::vector<uint8_t> pFecs;

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

        int m_lineNr = 0;
        int m_eventNr = 0;
        std::map<std::pair<uint8_t, uint8_t>, double> m_lowestCommonTriggerTimestamp_plane;
        std::map<uint8_t, double> m_lowestCommonTriggerTimestamp_det;

        std::map<uint8_t, double> m_deltaTriggerTimestamp;
        std::map<uint8_t, double> m_oldTriggerTimestamp;

        uint16_t m_oldBcId = 0;
        uint8_t m_oldVmmId = 0;
        uint8_t m_oldFecId = 0;

        std::map<std::pair<uint8_t,uint8_t>,HitContainer> m_hits;
        std::map<std::pair<uint8_t,uint8_t>,HitContainer> m_hits_new;
        std::map<std::pair<uint8_t,uint8_t>,ClusterVector> m_clusters;
        std::map<std::pair<uint8_t,uint8_t>,ClusterVector> m_clusters_new;
        std::map<uint8_t,CommonClusterVector> m_clustersXY;

	int m_cluster_id = 0;
	int m_clusterXY_id = 0;
        std::map<std::pair<uint8_t,uint8_t>,int> m_clusterCnt;
        std::map<uint8_t,int> m_clusterCnt_XY;


        RootFile* m_rootFile;

	int m_errors = 0;

        std::map<std::pair<uint8_t, uint8_t>, std::vector<int>> m_stats_maxDeltaTime;
        std::map<std::pair<uint8_t, uint8_t>, std::vector<int>> m_stats_maxMissingStrip;
        std::map<std::pair<uint8_t, uint8_t>, std::vector<int>> m_stats_deltaSpan;
        std::map<uint8_t, std::vector<int>> m_stats_deltaPlane;

        std::map<uint8_t, int> m_stats_overflow;
        std::map<uint8_t, int> m_stats_timeError;

};

