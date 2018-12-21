#pragma once

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
	static RootFile* GetInstance(std::string fileName, int channels_x, int channels_y, bool analyzeChannels,
			bool createHits);
	static void Dispose();

	void WriteRootFile();
	void FillTree();
	void SaveHistograms();

	void SaveHitsX(HitNMX&& hit);
	void SaveHitsY(HitNMX&& hit);
	void SaveClustersX(ClusterVector&& clusters);
	void SaveClustersY(ClusterVector&& clusters);
	void SaveClustersXY(CommonClusterVector&& clustersXY);

private:

	RootFile(TString fileName, int channels_x, int channels_y, bool analyzeChannels,  bool createHits);

	~RootFile();
	void CreateChannelHistograms();
	void AnalyzeChannels();
	static RootFile* m_rootFile;

	TFile * m_file = 0;
	TTree * m_tree = 0;
	TString m_fileName;



	int m_channels_x = 0;
	int m_channels_y = 0;
	bool m_analyzeChannels = false;
	bool m_createHits = false;
	int m_eventNr = 0;
	int m_nch = 0;
	int m_nchX = 0;
	int m_nchY = 0;

	int m_nclX = 0;
	int m_nclY = 0;
	int m_nclXY = 0;
	int m_nclXY_uTPC = 0;

	RootHitVector m_hitsX;
	RootHitVector m_hitsY;

	CommonClusterVector m_clusterXY;
	ClusterVector m_clusterX;
	ClusterVector m_clusterY;

	TH2D * m_TH2D_clusterXY;
	TH2D * m_TH2D_imageXY;
	TH1D * m_TH1D_clusterXYDeltaTime;
	TH2D * m_TH2D_chargeX;
	TH2D * m_TH2D_chargeY;
	TH2D * m_TH2D_chargeXY;
	TH2D * m_TH2D_sizeX;
	TH2D * m_TH2D_sizeY;
	TH2D * m_TH2D_sizeXY;



	TH1D * m_tdc_mean_x = 0;
	TH1D * m_tdc_min_x = 0;
	TH1D * m_tdc_max_x = 0;
	TH1D * m_tdc_stddev_x = 0;
	TH1D * m_tdc_totalrange_x = 0;
	TH1D * m_tdc_range_x = 0;
	TH1D * m_tdc_mean_y = 0;
	TH1D * m_tdc_min_y = 0;
	TH1D * m_tdc_max_y = 0;
	TH1D * m_tdc_stddev_y = 0;
	TH1D * m_tdc_totalrange_y = 0;
	TH1D * m_tdc_range_y = 0;

	TH1D * m_adc_mean_x = 0;
	TH1D * m_adc_min_x = 0;
	TH1D * m_adc_max_x = 0;
	TH1D * m_adc_stddev_x = 0;
	TH1D * m_adc_totalrange_x = 0;
	TH1D * m_adc_range_x = 0;
	TH1D * m_adc_mean_y = 0;
	TH1D * m_adc_min_y = 0;
	TH1D * m_adc_max_y = 0;
	TH1D * m_adc_stddev_y = 0;
	TH1D * m_adc_totalrange_y = 0;
	TH1D * m_adc_range_y = 0;

	TH1D * m_bcid_mean_x = 0;
	TH1D * m_bcid_min_x = 0;
	TH1D * m_bcid_max_x = 0;
	TH1D * m_bcid_stddev_x = 0;
	TH1D * m_bcid_totalrange_x = 0;
	TH1D * m_bcid_range_x = 0;
	TH1D * m_bcid_mean_y = 0;
	TH1D * m_bcid_min_y = 0;
	TH1D * m_bcid_max_y = 0;
	TH1D * m_bcid_stddev_y = 0;
	TH1D * m_bcid_totalrange_y = 0;
	TH1D * m_bcid_range_y = 0;

	ClassDef(RootFile,1)
};

