#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "Configuration.h"
#include "DataStructures.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include <TObject.h>
#include <TString.h>

class RootFile : public TObject {

public:
  static RootFile *GetInstance();
  static RootFile *GetInstance(Configuration &config);
  static void Dispose();

  void WriteRootFile();
  void SaveHistograms();
  void AddHits(Hit &&the_hit);
  void SaveHits();
  void SaveClustersPlane(ClusterVectorPlane &&clusters_plane);
  void SaveClustersDetector(ClusterVectorDetector &&clusters_detector);
  void SaveDate(double the_seconds, std::string the_date);

private:
  Configuration &m_config;
  RootFile(Configuration &config);

  ~RootFile();
  static RootFile *m_rootFile;
  TString m_fileName;
  TFile *m_file;
  TTree *m_tree_hits;
  TTree *m_tree_clusters_plane;
  TTree *m_tree_clusters_detector;

  int m_eventNr;

  Hit m_hit;
  ClusterPlane m_cluster_plane;
  ClusterDetector m_cluster_detector;

  HitVector m_hits;

  std::map<std::pair<uint8_t, std::string>, int> m_map_TH2D;
  std::map<std::pair<uint8_t, std::string>, int> m_map_TH1D;

  std::vector<TH1D *> m_TH1D;
  std::vector<TH2D *> m_TH2D;

  ClassDef(RootFile, 1)
};
