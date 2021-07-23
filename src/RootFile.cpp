#include "RootFile.h"
#include "TBufferJSON.h"
#include "TMath.h"
#include <TStyle.h>
#include <time.h>

RootFile *RootFile::m_rootFile = nullptr;

RootFile *RootFile::GetInstance() { return m_rootFile; }

RootFile *RootFile::GetInstance(Configuration &config) {
  if (!m_rootFile) {
    m_rootFile = new RootFile(config);
  }

  return m_rootFile;
}

void RootFile::Dispose() {
  if (m_rootFile) {
    m_rootFile->WriteRootFile();
    delete m_rootFile;
    m_rootFile = 0;
  }
}

void RootFile::WriteRootFile() {
  if(m_config.pSaveWhat >= 100) {
    SaveHistograms();
  }
  m_file->Write("", TObject::kOverwrite);
  m_file->Close();
}

RootFile::RootFile(Configuration &config) : m_config(config) {
  m_fileName = m_config.pRootFilename.c_str();
  m_file = TFile::Open(m_fileName, "RECREATE");
  m_eventNr = 0;

  m_tree_hits = new TTree("hits", "hits");
  m_tree_hits->SetDirectory(m_file);
  m_tree_clusters_plane = new TTree("clusters_plane", "clusters plane");
  m_tree_clusters_plane->SetDirectory(m_file);
  m_tree_clusters_detector = new TTree("clusters_detector", "clusters detector");
  m_tree_clusters_detector->SetDirectory(m_file);

  switch (m_config.pSaveWhat) {
    case 1: m_tree_hits->Branch("hits", &m_hits);
            break;
    case 10: m_tree_clusters_plane->Branch("clusters_plane", &m_clusters_plane);
            break;
    case 11: m_tree_hits->Branch("hits", &m_hits);
            m_tree_clusters_plane->Branch("clusters_plane", &m_clusters_plane);
            break;
    case 100: m_tree_clusters_detector->Branch("clusters_detector", &m_clusters_detector);
            break;
    case 101: m_tree_hits->Branch("hits", &m_hits);
              m_tree_clusters_detector->Branch("clusters_detector", &m_clusters_detector);
            break;
    case 110: m_tree_clusters_plane->Branch("clusters_plane", &m_clusters_plane);
              m_tree_clusters_detector->Branch("clusters_detector", &m_clusters_detector);
            break;
    case 111: m_tree_hits->Branch("hits", &m_hits);
              m_tree_clusters_plane->Branch("clusters_plane", &m_clusters_plane);
              m_tree_clusters_detector->Branch("clusters_detector", &m_clusters_detector);
            break;

  }

  if(m_config.pSaveWhat >= 100) {
    TH2D *h2;
    TH1D *h1;
    std::string name = "";
    int cnt1D = 0;
    int cnt2D = 0;
    for (auto const &det : m_config.pDets) {
      auto dp0 = std::make_pair(det.first, 0);
      auto dp1 = std::make_pair(det.first, 1);
      if (m_config.GetAxes(dp0) && m_config.GetAxes(dp1)) {
        name = std::to_string(det.first) + "_delta_time_planes";
        h1 = new TH1D(name.c_str(), name.c_str(), 1000, -500, 500);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(std::make_pair(
            std::make_pair(det.first, "delta_time_planes"), cnt1D));
        cnt1D++;

        name = std::to_string(det.first) + "_delta_time_utpc_planes";
        h1 = new TH1D(name.c_str(), name.c_str(), 1000, -500, 500);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(std::make_pair(
            std::make_pair(det.first, "delta_time_utpc_planes"), cnt1D));
        cnt1D++;

        name = std::to_string(det.first) + "_delta_time_charge2_planes";
        h1 = new TH1D(name.c_str(), name.c_str(), 1000, -500, 500);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(std::make_pair(
            std::make_pair(det.first, "delta_time_charge2_planes"), cnt1D));
        cnt1D++;

        name = std::to_string(det.first) + "_dt0";
        h1 = new TH1D(name.c_str(), name.c_str(), 11000, -500, 50000);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(
            std::make_pair(std::make_pair(det.first, "dt0"), cnt1D));
        cnt1D++;

        name = std::to_string(det.first) + "_dt1";
        h1 = new TH1D(name.c_str(), name.c_str(), 11000, -500, 50000);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(
            std::make_pair(std::make_pair(det.first, "dt1"), cnt1D));
        cnt1D++;

        int n0 = m_config.pChannels[std::make_pair(det.first, 0)];
        int n1 = m_config.pChannels[std::make_pair(det.first, 1)];
        int r0 = 0;
        int r1 = 0;

        if (m_config.pTransform.size() == m_config.pDets.size()) {
          auto tx = m_config.pTransformX[m_config.pDets[det.first]];
          auto ty = m_config.pTransformY[m_config.pDets[det.first]];
          auto tz = m_config.pTransformZ[m_config.pDets[det.first]];

          double t0 =
              n0 * std::get<0>(tx) + 0 * std::get<1>(tx) + std::get<3>(tx);
          double t1 =
              n0 * std::get<0>(ty) + 0 * std::get<1>(ty) + std::get<3>(ty);
          double t2 =
              0 * std::get<0>(tx) + n1 * std::get<1>(tx) + std::get<3>(tx);
          double t3 =
              0 * std::get<0>(ty) + n1 * std::get<1>(ty) + std::get<3>(ty);
          double max1 = std::max(t0, t1);
          double max2 = std::max(t2, t3);
          n0 = std::max(max1, max2) + 1;
          n1 = n0;
          double min1 = std::min(t0, t1);
          double min2 = std::min(t2, t3);
          r0 = std::min(min1, min2) - 1;
          r1 = r0;
          // std::cout << t0 << " " << t1 << " " << t2 << " " << t3 << " " << n0
          // << " " << r0 << std::endl;
        }

        name = std::to_string(det.first) + "_cluster";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, r0, n0, n1 * 4, r1, n1);
        h2->GetXaxis()->SetTitle("x");
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "cluster"), cnt2D));

        cnt2D++;

        name = std::to_string(det.first) + "_cluster_utpc";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, r0, n0, n1 * 4, r1, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "cluster_utpc"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_cluster_charge2";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, r0, n0, n1 * 4, r1, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "cluster_charge2"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_cluster_algo";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, r0, n0, n1 * 4, r1, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "cluster_algo"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_size_plane0";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, r0, n0, n1 * 4, r1, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "size_plane0"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_size_plane1";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, r0, n0, n1 * 4, r1, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "size_plane1"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_size_plane01";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, r0, n0, n1 * 4, r1, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "size_plane01"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_charge_plane0";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, r0, n0, n1 * 4, r1, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "charge_plane0"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_charge_plane1";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, r0, n0, n1 * 4, r1, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "charge_plane1"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_charge_plane01";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, r0, n0, n1 * 4, r1, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "charge_plane01"), cnt2D));
        cnt2D++;
      }
    }
  }
  std::cout << "ROOT file " << m_fileName << " created!" << std::endl;
}

RootFile::~RootFile() {}

void RootFile::AddHits(Hit &&the_hit) { m_hits.emplace_back(the_hit); }

void RootFile::SaveHits() {
  if (m_hits.size() > 0) {
    m_tree_hits->Fill();
    m_hits.clear();
  }
}

void RootFile::SaveClustersPlane(ClusterVectorPlane &&clusters_plane) {
  if (clusters_plane.size() > 0) {
    for (auto &it : clusters_plane) {
      if(std::find (m_config.pSaveClustersPlane.begin(), m_config.pSaveClustersPlane.end(), it.det) != m_config.pSaveClustersPlane.end()) {
        m_clusters_plane.push_back(it);
      }
    }
    m_tree_clusters_plane->Fill();
    m_clusters_plane.clear();
  }
}

void RootFile::SaveClustersDetector(ClusterVectorDetector &&clusters_detector) {
  
  for (auto &it : clusters_detector) {
    if(std::find (m_config.pSaveClustersDetector.begin(), m_config.pSaveClustersDetector.end(), it.det) != m_config.pSaveClustersDetector.end()) {
      m_clusters_detector.push_back(it);
      
      int idx = m_map_TH1D[std::make_pair(it.det, "delta_time_planes")];
      m_TH1D[idx]->Fill(it.time0 - it.time1);

      idx = m_map_TH1D[std::make_pair(it.det, "delta_time_utpc_planes")];
      m_TH1D[idx]->Fill(it.time0_utpc - it.time1_utpc);

      idx = m_map_TH1D[std::make_pair(it.det, "delta_time_charge2_planes")];
      m_TH1D[idx]->Fill(it.time0_charge2 - it.time1_charge2);

      idx = m_map_TH1D[std::make_pair(it.det, "dt0")];
      m_TH1D[idx]->Fill(it.dt0);

      idx = m_map_TH1D[std::make_pair(it.det, "dt1")];
      m_TH1D[idx]->Fill(it.dt1);

      idx = m_map_TH2D[std::make_pair(it.det, "cluster")];
      m_TH2D[idx]->Fill(it.pos0, it.pos1);

      idx = m_map_TH2D[std::make_pair(it.det, "cluster_utpc")];
      m_TH2D[idx]->Fill(it.pos0_utpc, it.pos1_utpc);

      idx = m_map_TH2D[std::make_pair(it.det, "cluster_charge2")];
      m_TH2D[idx]->Fill(it.pos0_charge2, it.pos1_charge2);

      idx = m_map_TH2D[std::make_pair(it.det, "cluster_algo")];
      m_TH2D[idx]->Fill(it.pos0_algo, it.pos1_algo);

      idx = m_map_TH2D[std::make_pair(it.det, "size_plane0")];
      m_TH2D[idx]->Fill(it.pos0, it.pos1, it.size0);

      idx = m_map_TH2D[std::make_pair(it.det, "size_plane1")];
      m_TH2D[idx]->Fill(it.pos0, it.pos1, it.size1);

      idx = m_map_TH2D[std::make_pair(it.det, "size_plane01")];
      m_TH2D[idx]->Fill(it.pos0, it.pos1, it.size0 + it.size1);

      idx = m_map_TH2D[std::make_pair(it.det, "charge_plane0")];
      m_TH2D[idx]->Fill(it.pos0, it.pos1, it.adc0);

      idx = m_map_TH2D[std::make_pair(it.det, "charge_plane1")];
      m_TH2D[idx]->Fill(it.pos0, it.pos1, it.adc1);

      idx = m_map_TH2D[std::make_pair(it.det, "charge_plane01")];
      m_TH2D[idx]->Fill(it.pos0, it.pos1, it.adc0 + it.adc1);
    }
  }
  if (m_clusters_detector.size() > 0) {
    m_tree_clusters_detector->Fill();
    m_clusters_detector.clear();
  }
}

void RootFile::SaveHistograms() {
  for (auto const &h1 : m_TH1D) {
    h1->Write("", TObject::kOverwrite);
  }
  for (auto const &det : m_config.pDets) {
    auto dp0 = std::make_pair(det.first, 0);
    auto dp1 = std::make_pair(det.first, 1);
    if (m_config.GetAxes(dp0) && m_config.GetAxes(dp1)) {
      if (m_config.createJSON) {
        int id = m_map_TH2D[std::make_pair(det.first, "cluster")];

        TString jsonFilename = m_fileName;
        jsonFilename.ReplaceAll(".root", "");

        TString json = TBufferJSON::ToJSON(m_TH2D[id], 3);
        std::ofstream f1;
        f1.open(jsonFilename + "_detector" + std::to_string(det.first) +
                    "_cluster.json",
                std::ios::out);
        f1 << json;
        f1.close();

        id = m_map_TH2D[std::make_pair(det.first, "cluster_utpc")];
        json = TBufferJSON::ToJSON(m_TH2D[id], 3);
        std::ofstream f2;
        f2.open(jsonFilename + "_detector" + std::to_string(det.first) +
                    "_cluster_utpc.json",
                std::ios::out);
        f2 << json;
        f2.close();

        id = m_map_TH2D[std::make_pair(det.first, "cluster_charge2")];
        json = TBufferJSON::ToJSON(m_TH2D[id], 3);
        std::ofstream f3;
        f3.open(jsonFilename + "_detector" + std::to_string(det.first) +
                    "_cluster_charge2.json",
                std::ios::out);
        f3 << json;
        f3.close();

        id = m_map_TH2D[std::make_pair(det.first, "cluster_algo")];
        json = TBufferJSON::ToJSON(m_TH2D[id], 3);
        std::ofstream f4;
        f4.open(jsonFilename + "_detector" + std::to_string(det.first) +
                    "_cluster_algo.json",
                std::ios::out);
        f4 << json;
        f4.close();
      }
      int n0 = m_config.pChannels[std::make_tuple(det.first, 0)] * 4;
      int n1 = m_config.pChannels[std::make_tuple(det.first, 1)] * 4;
      int n = 0;
      for (int p0 = 1; p0 <= n0; p0++) {
        for (int p1 = 1; p1 <= n1; p1++) {
          int idx = m_map_TH2D[std::make_pair(det.first, "cluster")];
          int cnt = m_TH2D[idx]->GetBinContent(p0, p1);
          double val = 0;
          if (cnt > 0) {
            n++;
            idx = m_map_TH2D[std::make_pair(det.first, "size_plane0")];
            val = m_TH2D[idx]->GetBinContent(p0, p1) / cnt;
            m_TH2D[idx]->SetBinContent(p0, p1, val);

            idx = m_map_TH2D[std::make_pair(det.first, "size_plane1")];
            val = m_TH2D[idx]->GetBinContent(p0, p1) / cnt;
            m_TH2D[idx]->SetBinContent(p0, p1, val);

            idx = m_map_TH2D[std::make_pair(det.first, "size_plane01")];
            val = m_TH2D[idx]->GetBinContent(p0, p1) / cnt;
            m_TH2D[idx]->SetBinContent(p0, p1, val);

            idx = m_map_TH2D[std::make_pair(det.first, "charge_plane0")];
            val = m_TH2D[idx]->GetBinContent(p0, p1) / cnt;
            m_TH2D[idx]->SetBinContent(p0, p1, val);

            idx = m_map_TH2D[std::make_pair(det.first, "charge_plane1")];
            val = m_TH2D[idx]->GetBinContent(p0, p1) / cnt;
            m_TH2D[idx]->SetBinContent(p0, p1, val);

            idx = m_map_TH2D[std::make_pair(det.first, "charge_plane01")];
            val = m_TH2D[idx]->GetBinContent(p0, p1) / cnt;
            m_TH2D[idx]->SetBinContent(p0, p1, val);
          }
        }
      }

      int idx = m_map_TH2D[std::make_pair(det.first, "size_plane0")];
      m_TH2D[idx]->SetEntries(n);
      idx = m_map_TH2D[std::make_pair(det.first, "size_plane1")];
      m_TH2D[idx]->SetEntries(n);
      idx = m_map_TH2D[std::make_pair(det.first, "size_plane01")];
      m_TH2D[idx]->SetEntries(n);
      idx = m_map_TH2D[std::make_pair(det.first, "charge_plane0")];
      m_TH2D[idx]->SetEntries(n);
      idx = m_map_TH2D[std::make_pair(det.first, "charge_plane1")];
      m_TH2D[idx]->SetEntries(n);
      idx = m_map_TH2D[std::make_pair(det.first, "charge_plane01")];
      m_TH2D[idx]->SetEntries(n);
    }
  }
}
