#pragma once

#define NUMDET 10

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <TObject.h>
#include <TString.h>
#include "DataStructures.h"

class RootFile: public TObject {

public:

        static RootFile* GetInstance();
        static RootFile* GetInstance(TString fileName, std::map<uint8_t, uint8_t> dets, std::map<std::pair<uint8_t,uint8_t>, uint32_t> channels, bool analyzeChannels,  bool createHits);
        static void Dispose();

	void WriteRootFile();
	void FillTree();
	void SaveHistograms();
        void AddHits(HitNMX&& hits);
        void SaveHits();
        void SaveClusters(ClusterVector&& clusters);
        void SaveClustersXY(CommonClusterVector&& clustersXY);

private:

    RootFile(TString fileName, std::map<uint8_t, uint8_t> dets, std::map<std::pair<uint8_t,uint8_t>, uint32_t> channels, bool analyzeChannels,  bool createHits);

	~RootFile();
	void CreateChannelHistograms();
	void AnalyzeChannels();
	static RootFile* m_rootFile;

        TFile * m_file;
        TTree * m_tree;
	TString m_fileName;
        std::map<uint8_t, uint8_t> m_dets;
        std::map<std::pair<uint8_t, uint8_t>,uint32_t> m_channels;

	bool m_analyzeChannels = false;
	bool m_createHits = false;
        int m_eventNr;

        RootHitVector m_hits;
        CommonClusterVector m_clustersXY;
        ClusterVector m_clusters;


        TH2D * m_TH2D_clusterXY[NUMDET];
        TH2D * m_TH2D_imageXY[NUMDET];
        TH1D * m_TH1D_clusterXYDeltaTime[NUMDET];
        TH2D * m_TH2D_chargeX[NUMDET];
        TH2D * m_TH2D_chargeY[NUMDET];
        TH2D * m_TH2D_chargeXY[NUMDET];
        TH2D * m_TH2D_sizeX[NUMDET];
        TH2D * m_TH2D_sizeY[NUMDET];
        TH2D * m_TH2D_sizeXY[NUMDET];



        TH1D * m_tdc_mean_x[NUMDET];
        TH1D * m_tdc_min_x[NUMDET];
        TH1D * m_tdc_max_x[NUMDET];
        TH1D * m_tdc_stddev_x[NUMDET];
        TH1D * m_tdc_totalrange_x[NUMDET];
        TH1D * m_tdc_range_x[NUMDET];
        TH1D * m_tdc_mean_y[NUMDET];
        TH1D * m_tdc_min_y[NUMDET];
        TH1D * m_tdc_max_y[NUMDET];
        TH1D * m_tdc_stddev_y[NUMDET];
        TH1D * m_tdc_totalrange_y[NUMDET];
        TH1D * m_tdc_range_y[NUMDET];

        TH1D * m_adc_mean_x[NUMDET];
        TH1D * m_adc_min_x[NUMDET];
        TH1D * m_adc_max_x[NUMDET];
        TH1D * m_adc_stddev_x[NUMDET];
        TH1D * m_adc_totalrange_x[NUMDET];
        TH1D * m_adc_range_x[NUMDET];
        TH1D * m_adc_mean_y[NUMDET];
        TH1D * m_adc_min_y[NUMDET];
        TH1D * m_adc_max_y[NUMDET];
        TH1D * m_adc_stddev_y[NUMDET];
        TH1D * m_adc_totalrange_y[NUMDET];
        TH1D * m_adc_range_y[NUMDET];

        TH1D * m_bcid_mean_x[NUMDET];
        TH1D * m_bcid_min_x[NUMDET];
        TH1D * m_bcid_max_x[NUMDET];
        TH1D * m_bcid_stddev_x[NUMDET];
        TH1D * m_bcid_totalrange_x[NUMDET];
        TH1D * m_bcid_range_x[NUMDET];
        TH1D * m_bcid_mean_y[NUMDET];
        TH1D * m_bcid_min_y[NUMDET];
        TH1D * m_bcid_max_y[NUMDET];
        TH1D * m_bcid_stddev_y[NUMDET];
        TH1D * m_bcid_totalrange_y[NUMDET];
        TH1D * m_bcid_range_y[NUMDET];

	ClassDef(RootFile,1)
};

