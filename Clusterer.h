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

#include "RootFile.h"
#include "Statistics.h"

class Clusterer {
public:
        Clusterer(Configuration& config, Statistics &stats);


	~Clusterer();

	// Analyzing and storing the hits
        bool AnalyzeHits(double srsTimestamp, uint8_t fecId, uint8_t vmmId, uint16_t chNo, uint16_t bcid, uint16_t tdc, uint16_t adc, bool overThresholdFlag, float chipTime);
        void StoreHits(uint8_t det, uint8_t plane, int pos, uint16_t adc, double chipTime, bool overThresholdFlag);

        // Analyzing and storing the clusters in plane 0 and 1
        void AnalyzeClustersPlane(uint8_t det, uint8_t plane);

	//Select hits that are ready to be clustered in time
        bool ChooseHitsToBeClustered(uint8_t det, uint8_t plane);

        bool ChooseClustersToBeMatched(uint8_t det, uint8_t plane);

        // Matching the clusters that are common between detector planes
        void AnalyzeClustersDetector(uint8_t det);

        int ClusterByTime(uint8_t det,uint8_t plane);
        int ClusterByStrip(uint8_t det,uint8_t plane,ClusterContainer &cluster, uint16_t maxDeltaTime);

        void StoreClusters(uint8_t det, uint8_t plane, std::vector<double>& strips, std::vector<double>& times,
			double pos_utpc, double time_utpc, 
                        double pos_center_charge, double time_center_charge, 
                        double pos_center_charge2, double time_center_charge2, 
                        uint16_t clusterSize, uint32_t clusterADC, uint16_t maxDeltaTime, uint16_t maxMissingStrip, uint16_t deltaSpan);

        void MatchClustersDetector(uint8_t det);

	void PrintStats();

        //bool IsFecChipInDetectorPlane(uint8_t fecId, uint8_t det, uint8_t chip, uint8_t planeId);

        std::pair<int, int> GetDetectorPlane(std::pair<uint8_t, uint8_t> fecChip);
        // Helper methods that map channels to strips
        int GetChannel(std::pair<uint8_t, uint8_t> fecChip, int chNo);



#ifdef USE_ROOT
	void createRootFile(string fileName);
#endif
private:
        Configuration& m_config;
        Statistics &m_stats;
      
        int m_lineNr = 0;
        int m_eventNr = 0;
        
        double last_time0 = 0;
        double last_time1 = 0;
        double last_time0_utpc = 0;
        double last_time1_utpc = 0;
        double last_time0_charge2 = 0;
        double last_time1_charge2 = 0;

        uint16_t m_oldBcId = 0;
        uint8_t m_oldVmmId = 0;
        uint8_t m_oldFecId = 0;

        std::map<std::pair<uint8_t,uint8_t>,HitContainer> m_hits;
        std::map<std::pair<uint8_t,uint8_t>,HitContainer> m_hits_new;
        std::map<std::pair<uint8_t,uint8_t>,ClusterVectorPlane> m_clusters;
        std::map<std::pair<uint8_t,uint8_t>,ClusterVectorPlane> m_clusters_new;
        std::map<uint8_t,ClusterVectorDetector> m_clusters_detector;

	int m_cluster_id = 0;
	int m_cluster_detector_id = 0;
        std::map<std::pair<uint8_t,uint8_t>,int> m_clusters_cnt;
        std::map<uint8_t,int> m_clusters_detector_cnt;


        RootFile* m_rootFile;

	

};

