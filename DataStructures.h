#pragma once
#include <vector>
struct HitNMX {
        uint32_t id;
	uint32_t eventNr;
        uint8_t detId;
        uint8_t planeId;
        uint8_t fecId;
        uint8_t vmmId;
	double triggerTimestamp;
	double totalTime;
	uint8_t chNo;
	uint16_t position;
	uint16_t bcid;
	uint16_t tdc;
	uint16_t adc;
	bool overThresholdFlag;
	float chipTime;
};

struct ClusterNMX {
        uint32_t id;
        uint8_t detId;
        uint8_t planeId;
        uint16_t size;
        uint32_t adc;
	float position;
	double time;
	double utpcTime;
        double utpcPosition;
	double centerOfTime;
        double centerOfCharge;
	bool clusterXAndY;
	uint16_t maxDeltaTime;
	uint16_t maxMissingStrip;
        uint16_t spanCluster;
	std::vector<float> strips;
	std::vector<double> times;

};

struct CommonClusterNMX {
        uint32_t id;
        uint8_t detId;
	uint32_t idX;
	uint32_t idY;
	uint16_t sizeX;
	uint16_t sizeY;
	uint32_t adcX;
	uint32_t adcY;
	float positionX;
	float positionY;
	double timeX;
	double timeY;
	int16_t deltaPlane;
        uint16_t spanClusterX;
        uint16_t spanClusterY;
	uint16_t maxDeltaTimeX;
	uint16_t maxDeltaTimeY;
	uint16_t maxMissingStripX;
	uint16_t maxMissingStripY;
	std::vector<float> stripsX;
	std::vector<double> timesX;
	std::vector<float> stripsY;
	std::vector<double> timesY;
};

using std::string;


using HitTuple = std::tuple<double, uint16_t, uint16_t>;
using ClusterTuple = std::tuple<uint16_t, double, uint16_t>;
using HitContainer = std::vector<HitTuple>;
using ClusterContainer = std::vector<ClusterTuple>;

using ClusterVector = std::vector<ClusterNMX>;
using CommonClusterVector = std::vector<CommonClusterNMX>;

using RootHitVector = std::vector<HitNMX>;

