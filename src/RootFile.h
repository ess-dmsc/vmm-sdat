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
#include "Configuration.h"

class RootFile: public TObject {

public:

    static RootFile* GetInstance();
    static RootFile* GetInstance(Configuration& config);
    static void Dispose();

    void WriteRootFile();
    void SaveHistograms();
    void AddHits(Hit&& the_hit);
    void SaveHits();
    void SaveClustersPlane(ClusterVectorPlane&& clusters_plane);
    void SaveClustersDetector(ClusterVectorDetector&& clusters_detector);

private:
    Configuration& m_config;
    RootFile(Configuration& config);

    ~RootFile();
    static RootFile* m_rootFile;
    TString m_fileName;
    TFile * m_file;
    TTree * m_tree;
    int m_eventNr;

    HitVector m_hits;
    ClusterVectorDetector m_clusters_detector;
    ClusterVectorPlane m_clusters_plane;
    TrackVector m_tracks;
	
    std::map<std::pair<uint8_t, std::string>, int> m_map_TH2D;
    std::map<std::pair<uint8_t, std::string>, int> m_map_TH1D;

    std::vector<TH1D*> m_TH1D;
    std::vector<TH2D*> m_TH2D;
 

    ClassDef(RootFile,1)
};

