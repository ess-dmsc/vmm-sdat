#include "RootFile.h"
#include "TBufferJSON.h"
#include "TMath.h"
#include <RooDouble.h>
#include <TStyle.h>
#include <time.h>

#define BINNING_FACTOR 1

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
  if (m_config.pSaveWhat >= 100) {
    SaveHistograms();
  }
  m_file->Write("", TObject::kOverwrite);
  m_file->Close();
}

void RootFile::SaveDate(double the_seconds, std::string the_date) {
  TString str_date = the_date;
  TString str_time = Form("%f", the_seconds);
  TNamed unixtime;
  unixtime.SetName("unixtime");
  unixtime.SetTitle(str_time);
  unixtime.Write();

  TNamed datetime;
  datetime.SetName("date");
  datetime.SetTitle(str_date);
  datetime.Write();

  TNamed channelMapping;
  channelMapping.SetName("channel mapping");
  channelMapping.SetTitle(m_config.pChannelMapping.c_str());
  channelMapping.Write();

  TNamed chargeRatioLower;
  chargeRatioLower.SetName("charge_ratio_lower");
  std::string sChargeRatioLower = std::to_string(m_config.pChargeRatioLower);
  chargeRatioLower.SetTitle(sChargeRatioLower.c_str());
  chargeRatioLower.Write();

  TNamed chargeRatioUpper;
  chargeRatioUpper.SetName("charge_ratio_lower");
  std::string sChargeRatioUpper = std::to_string(m_config.pChargeRatioUpper);
  chargeRatioUpper.SetTitle(sChargeRatioUpper.c_str());
  chargeRatioUpper.Write();

  TNamed bc;
  bc.SetName("bc");
  std::string sBC = std::to_string(m_config.pBC);
  bc.SetTitle(sBC.c_str());
  bc.Write();

  TNamed tac;
  tac.SetName("tac");
  std::string sTAC = std::to_string(m_config.pTAC);
  tac.SetTitle(sTAC.c_str());
  tac.Write();

  TNamed minClusterSize;
  minClusterSize.SetName("min_cluster_size");
  std::string sMinClusterSize = std::to_string(m_config.pMinClusterSize);
  minClusterSize.SetTitle(sMinClusterSize.c_str());
  minClusterSize.Write();

  TNamed coincidentClusterSize;
  coincidentClusterSize.SetName("coincident_cluster_size");
  std::string sCoincidentClusterSize =
      std::to_string(m_config.pCoincidentClusterSize);
  coincidentClusterSize.SetTitle(sCoincidentClusterSize.c_str());
  coincidentClusterSize.Write();

  TNamed deltaTimeHits;
  deltaTimeHits.SetName("delta_time_hits");
  std::string sDeltaTimeHits = std::to_string(m_config.pDeltaTimeHits);
  deltaTimeHits.SetTitle(sDeltaTimeHits.c_str());
  deltaTimeHits.Write();

  TNamed missingStripsClusterX;
  missingStripsClusterX.SetName("missing_strips");
  std::string sMissingStripsClusterX =
      std::to_string(m_config.pMissingStripsClusterX);
  missingStripsClusterX.SetTitle(sMissingStripsClusterX.c_str());
  missingStripsClusterX.Write();

  TNamed spanClusterTime;
  spanClusterTime.SetName("span_cluster_time");
  std::string sSpanClusterTime = std::to_string(m_config.pSpanClusterTime);
  spanClusterTime.SetTitle(sSpanClusterTime.c_str());
  spanClusterTime.Write();

  TNamed deltaTimePlanes;
  deltaTimePlanes.SetName("delta_time_planes");
  std::string sDeltaTimePlanes = std::to_string(m_config.pDeltaTimePlanes);
  deltaTimePlanes.SetTitle(sDeltaTimePlanes.c_str());
  deltaTimePlanes.Write();
}

void RootFile::FillCalibHistos(uint16_t fec, uint8_t vmm, uint8_t ch, float adc,
                               float adc_corrected, float chip_time,
                               float chip_time_corrected) {

  int idx = m_map_calib_TH2D[std::make_tuple(fec, vmm, "adc_without_calib")];
  m_calib_TH2D[idx]->Fill(ch, adc);
  idx = m_map_calib_TH2D[std::make_tuple(fec, vmm, "adc_with_calib")];
  m_calib_TH2D[idx]->Fill(ch, adc_corrected);
  idx = m_map_calib_TH2D[std::make_tuple(fec, vmm, "time_without_calib")];
  m_calib_TH2D[idx]->Fill(ch, chip_time);
  idx = m_map_calib_TH2D[std::make_tuple(fec, vmm, "time_with_calib")];
  m_calib_TH2D[idx]->Fill(ch, chip_time_corrected);
}

RootFile::RootFile(Configuration &config) : m_config(config) {
  m_fileName = m_config.pRootFilename.c_str();
  m_file = TFile::Open(m_fileName, "RECREATE");
  m_eventNr = 0;

  if (m_config.useCalibration && m_config.calibrationHistogram) {
    TH2D *h2;
    std::string name = "";
    int cntCal = 0;
    for (int i = 0; i < m_config.pVMMs.size(); i++) {
      auto tuple = m_config.pVMMs[i];
      auto det = std::get<0>(tuple);
      auto plane = std::get<1>(tuple);
      auto fec = std::get<2>(tuple);
      auto vmm = std::get<3>(tuple);

      name = "fec" + std::to_string(fec) + "_vmm" + std::to_string(vmm) +
             "_adc_with_calib";
      h2 = new TH2D(name.c_str(), name.c_str(), 64, 0, 64, 1023, 0, 1023);
      m_calib_TH2D.push_back(h2);
      m_map_calib_TH2D.emplace(
          std::make_pair(std::make_tuple(fec, vmm, "adc_with_calib"), cntCal));
      cntCal++;
      name = "fec" + std::to_string(fec) + "_vmm" + std::to_string(vmm) +
             "_adc_without_calib";
      h2 = new TH2D(name.c_str(), name.c_str(), 64, 0, 64, 1023, 0, 1023);
      m_calib_TH2D.push_back(h2);
      m_map_calib_TH2D.emplace(std::make_pair(
          std::make_tuple(fec, vmm, "adc_without_calib"), cntCal));
      cntCal++;
      name = "fec" + std::to_string(fec) + "_vmm" + std::to_string(vmm) +
             "_time_with_calib";
      h2 = new TH2D(name.c_str(), name.c_str(), 64, 0, 64, 350, 0, 35);
      m_calib_TH2D.push_back(h2);
      m_map_calib_TH2D.emplace(
          std::make_pair(std::make_tuple(fec, vmm, "time_with_calib"), cntCal));
      cntCal++;
      name = "fec" + std::to_string(fec) + "_vmm" + std::to_string(vmm) +
             "_time_without_calib";
      h2 = new TH2D(name.c_str(), name.c_str(), 64, 0, 64, 350, 0, 35);
      m_calib_TH2D.push_back(h2);
      m_map_calib_TH2D.emplace(std::make_pair(
          std::make_tuple(fec, vmm, "time_without_calib"), cntCal));
      cntCal++;
    }
  }

  switch (m_config.pSaveWhat) {
  case 1:
    m_tree_hits = new TTree("hits", "hits");
    m_tree_hits->SetDirectory(m_file);
    m_tree_hits->Branch("hits", &m_hit);
    break;
  case 10:
    m_tree_clusters_plane = new TTree("clusters_plane", "clusters plane");
    m_tree_clusters_plane->SetDirectory(m_file);
    m_tree_clusters_plane->Branch("clusters_plane", &m_cluster_plane);
    break;
  case 11:
    m_tree_hits = new TTree("hits", "hits");
    m_tree_hits->SetDirectory(m_file);
    m_tree_clusters_plane = new TTree("clusters_plane", "clusters plane");
    m_tree_clusters_plane->SetDirectory(m_file);
    m_tree_hits->Branch("hits", &m_hit);
    m_tree_clusters_plane->Branch("clusters_plane", &m_cluster_plane);
    break;
  case 100:
    m_tree_clusters_detector =
        new TTree("clusters_detector", "clusters detector");
    m_tree_clusters_detector->SetDirectory(m_file);
    m_tree_clusters_detector->Branch("clusters_detector", &m_cluster_detector);
    break;
  case 101:
    m_tree_hits = new TTree("hits", "hits");
    m_tree_hits->SetDirectory(m_file);
    m_tree_clusters_detector =
        new TTree("clusters_detector", "clusters detector");
    m_tree_clusters_detector->SetDirectory(m_file);
    m_tree_hits->Branch("hits", &m_hit);
    m_tree_clusters_detector->Branch("clusters_detector", &m_cluster_detector);
    break;
  case 110:
    m_tree_clusters_plane = new TTree("clusters_plane", "clusters plane");
    m_tree_clusters_plane->SetDirectory(m_file);
    m_tree_clusters_detector =
        new TTree("clusters_detector", "clusters detector");
    m_tree_clusters_detector->SetDirectory(m_file);
    m_tree_clusters_plane->Branch("clusters_plane", &m_cluster_plane);
    m_tree_clusters_detector->Branch("clusters_detector", &m_cluster_detector);
    break;
  case 111:
    m_tree_hits = new TTree("hits", "hits");
    m_tree_hits->SetDirectory(m_file);
    m_tree_clusters_plane = new TTree("clusters_plane", "clusters plane");
    m_tree_clusters_plane->SetDirectory(m_file);
    m_tree_clusters_detector =
        new TTree("clusters_detector", "clusters detector");
    m_tree_clusters_detector->SetDirectory(m_file);
    m_tree_hits->Branch("hits", &m_hit);
    m_tree_clusters_plane->Branch("clusters_plane", &m_cluster_plane);
    m_tree_clusters_detector->Branch("clusters_detector", &m_cluster_detector);

    break;
  }

  if (m_config.pSaveWhat >= 100) {
    TH2D *h2;
    TH1D *h1;
    std::string name = "";
    int cnt1D = 0;
    int cnt2D = 0;

    m_max0 = -999999999;
    m_max1 = -999999999;
    m_min0 = 9999999;
    m_min1 = 9999999;
    m_bins0 = 0;
    m_bins1 = 0;

    for (auto const &det : m_config.pDets) {
      auto dp0 = std::make_pair(det.first, 0);
      auto dp1 = std::make_pair(det.first, 1);
      // 2D detectors
      if (!m_config.pIsPads[det.first] &&
          m_config.GetDetectorPlane(dp0) == true &&
          m_config.GetDetectorPlane(dp1) == true) {
        int n0 = m_config.pChannels[dp0];
        int n1 = m_config.pChannels[dp1];

        if (m_config.pTransform.size() != m_config.pDets.size()) {
          m_min0 = 0;
          m_min1 = 0;
          m_max0 = n0;
          m_max1 = n1;
          m_bins0 = n0 * BINNING_FACTOR;
          m_bins1 = n1 * BINNING_FACTOR;
        }
        // Transformation is applied
        else {
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

          if (t0 > m_max0) {
            m_max0 = t0;
          }
          if (t1 > m_max0) {
            m_max0 = t1;
          }
          if (t2 > m_max0) {
            m_max0 = t2;
          }
          if (t3 > m_max0) {
            m_max0 = t3;
          }
          if (t0 < m_min0) {
            m_min0 = t0;
          }
          if (t1 < m_min0) {
            m_min0 = t1;
          }
          if (t2 < m_min0) {
            m_min0 = t2;
          }
          if (t3 < m_min0) {
            m_min0 = t3;
          }
          m_max1 = m_max0;
          m_min1 = m_min0;
          m_bins0 = (m_max0 - m_min0) * BINNING_FACTOR;
          m_bins1 = (m_max1 - m_min1) * BINNING_FACTOR;
        }

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
        h1 = new TH1D(name.c_str(), name.c_str(), 11000, -1000, 109000);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(
            std::make_pair(std::make_pair(det.first, "dt0"), cnt1D));
        cnt1D++;

        name = std::to_string(det.first) + "_dt1";
        h1 = new TH1D(name.c_str(), name.c_str(), 11000, -1000, 109000);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(
            std::make_pair(std::make_pair(det.first, "dt1"), cnt1D));
        cnt1D++;

        name = std::to_string(det.first) + "_cluster";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);

        h2->GetXaxis()->SetTitle("x");
        h2->GetYaxis()->SetTitle("y");
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "cluster"), cnt2D));

        cnt2D++;

        name = std::to_string(det.first) + "_cluster_utpc";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "cluster_utpc"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_cluster_charge2";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(std::make_pair(
            std::make_pair(det.first, "cluster_charge2"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_cluster_algo";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "cluster_algo"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_size_plane0";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "size_plane0"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_size_plane1";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "size_plane1"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_size_plane01";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "size_plane01"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_charge_plane0";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "charge_plane0"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_charge_plane1";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "charge_plane1"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_charge_plane01";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "charge_plane01"), cnt2D));
        cnt2D++;
      }
      // Pad detectors
      else if (m_config.pIsPads[det.first]) {
        m_max0 = m_config.pChannels0[det.first] + 1;
        m_max1 = m_config.pChannels1[det.first] + 1;
        m_min0 = 0;
        m_min1 = 0;
        m_bins0 = m_max0 * BINNING_FACTOR;
        m_bins1 = m_max1 * BINNING_FACTOR;

        name = std::to_string(det.first) + "_cluster";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "cluster"), cnt2D));

        cnt2D++;

        name = std::to_string(det.first) + "_size";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "size"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_charge";
        h2 = new TH2D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0,
                      m_bins1, m_min1, m_max1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(
            std::make_pair(std::make_pair(det.first, "charge"), cnt2D));
        cnt2D++;
      }
      // 1D detector
      else if (m_config.GetDetectorPlane(dp0) == true) {
        if (m_config.GetDetectorPlane(dp0)) {
          m_max0 = m_config.pChannels[std::make_tuple(det.first, 0)];
        } else {
          m_max0 = m_config.pChannels[std::make_tuple(det.first, 1)];
        }
        m_min0 = 0;
        m_bins0 = m_max0 * BINNING_FACTOR;
        name = std::to_string(det.first) + "_cluster";
        h1 = new TH1D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(
            std::make_pair(std::make_pair(det.first, "cluster"), cnt1D));

        cnt1D++;

        name = std::to_string(det.first) + "_size";
        h1 = new TH1D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(
            std::make_pair(std::make_pair(det.first, "size"), cnt1D));
        cnt1D++;

        name = std::to_string(det.first) + "_charge";
        h1 = new TH1D(name.c_str(), name.c_str(), m_bins0, m_min0, m_max0);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(
            std::make_pair(std::make_pair(det.first, "charge"), cnt1D));
        cnt1D++;
      }
    }
  }
  if (m_config.useCalibration) {
  }
  std::cout << "ROOT file " << m_fileName << " created!" << std::endl;
}

RootFile::~RootFile() {}

void RootFile::AddHits(Hit &&the_hit) { m_hits.emplace_back(the_hit); }

void RootFile::SaveHits() {
  if (m_hits.size() > 0) {
    for (int n = 0; n < m_hits.size(); n++) {
      m_hit = m_hits[n];
      m_tree_hits->Fill();
    }
    m_hits.clear();
  }
}

void RootFile::SaveClustersPlane(ClusterVectorPlane &&clusters_plane) {
  if (clusters_plane.size() > 0) {
    for (auto &it : clusters_plane) {
      if (std::find(m_config.pSaveClustersPlane.begin(),
                    m_config.pSaveClustersPlane.end(),
                    it.det) != m_config.pSaveClustersPlane.end()) {
        auto detector_plane = std::make_pair(it.det, it.plane);
        if (m_config.GetDetectorPlane(detector_plane) ||
            m_config.pIsPads[it.det]) {
          m_cluster_plane = it;
          m_tree_clusters_plane->Fill();
        }
      }
    }
  }
}

void RootFile::SaveClustersDetector(ClusterVectorDetector &&clusters_detector) {
  for (auto &it : clusters_detector) {

    if (std::find(m_config.pSaveClustersDetector.begin(),
                  m_config.pSaveClustersDetector.end(),
                  it.det) != m_config.pSaveClustersDetector.end()) {
      auto dp0 = std::make_pair(it.det, 0);
      auto dp1 = std::make_pair(it.det, 1);
      // 2D detector
      if (!m_config.pIsPads[it.det] && m_config.GetDetectorPlane(dp0) == true &&
          m_config.GetDetectorPlane(dp1) == true) {
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
      // Pad detector
      else if (m_config.pIsPads[it.det]) {
        int idx = m_map_TH2D[std::make_pair(it.det, "cluster")];
        m_TH2D[idx]->Fill(it.pos0, it.pos1);

        idx = m_map_TH2D[std::make_pair(it.det, "size")];
        m_TH2D[idx]->Fill(it.pos0, it.pos1, it.size0);

        idx = m_map_TH2D[std::make_pair(it.det, "charge")];
        m_TH2D[idx]->Fill(it.pos0, it.pos1, it.adc0);
      }
      // 1D detector
      else if (m_config.GetDetectorPlane(dp0) == true) {

        int idx = m_map_TH1D[std::make_pair(it.det, "cluster")];
        m_TH1D[idx]->Fill(it.pos0);

        idx = m_map_TH1D[std::make_pair(it.det, "size")];
        m_TH1D[idx]->Fill(it.pos0, it.size0);

        idx = m_map_TH1D[std::make_pair(it.det, "charge")];
        m_TH1D[idx]->Fill(it.pos0, it.adc0);
      }
      m_cluster_detector = it;

      m_tree_clusters_detector->Fill();
    }
  }
}

void RootFile::SaveHistograms() {
  for (auto const &h1 : m_TH1D) {
    h1->Write("", TObject::kOverwrite);
  }
  for (auto const &det : m_config.pDets) {
    auto dp0 = std::make_pair(det.first, 0);
    auto dp1 = std::make_pair(det.first, 1);

    if (!m_config.pIsPads[det.first] && m_config.GetDetectorPlane(dp0) &&
        m_config.GetDetectorPlane(dp1)) {
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
      int n = 0;
      for (int b0 = 1; b0 <= m_bins0; b0++) {
        for (int b1 = 1; b1 <= m_bins1; b1++) {
          int idx = m_map_TH2D[std::make_pair(det.first, "cluster")];
          int cnt = m_TH2D[idx]->GetBinContent(b0, b1);
          double val = 0;
          if (cnt > 0) {
            n++;
            idx = m_map_TH2D[std::make_pair(det.first, "size_plane0")];
            val = m_TH2D[idx]->GetBinContent(b0, b1) / cnt;
            m_TH2D[idx]->SetBinContent(b0, b1, val);

            idx = m_map_TH2D[std::make_pair(det.first, "size_plane1")];
            val = m_TH2D[idx]->GetBinContent(b0, b1) / cnt;
            m_TH2D[idx]->SetBinContent(b0, b1, val);

            idx = m_map_TH2D[std::make_pair(det.first, "size_plane01")];
            val = m_TH2D[idx]->GetBinContent(b0, b1) / cnt;
            m_TH2D[idx]->SetBinContent(b0, b1, val);

            idx = m_map_TH2D[std::make_pair(det.first, "charge_plane0")];
            val = m_TH2D[idx]->GetBinContent(b0, b1) / cnt;
            m_TH2D[idx]->SetBinContent(b0, b1, val);

            idx = m_map_TH2D[std::make_pair(det.first, "charge_plane1")];
            val = m_TH2D[idx]->GetBinContent(b0, b1) / cnt;
            m_TH2D[idx]->SetBinContent(b0, b1, val);

            idx = m_map_TH2D[std::make_pair(det.first, "charge_plane01")];
            val = m_TH2D[idx]->GetBinContent(b0, b1) / cnt;
            m_TH2D[idx]->SetBinContent(b0, b1, val);
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
    } else if (m_config.pIsPads[det.first]) {
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
      }
      int n = 0;

      for (int b0 = 1; b0 <= m_bins0; b0++) {
        for (int b1 = 1; b1 <= m_bins1; b1++) {
          int idx = m_map_TH2D[std::make_pair(det.first, "cluster")];
          int cnt = m_TH2D[idx]->GetBinContent(b0, b1);
          double val = 0;
          if (cnt > 0) {
            n++;
            idx = m_map_TH2D[std::make_pair(det.first, "size")];
            val = m_TH2D[idx]->GetBinContent(b0, b1) / cnt;
            m_TH2D[idx]->SetBinContent(b0, b1, val);

            idx = m_map_TH2D[std::make_pair(det.first, "charge")];
            val = m_TH2D[idx]->GetBinContent(b0, b1) / cnt;
            m_TH2D[idx]->SetBinContent(b0, b1, val);
          }
        }
      }
      int idx = m_map_TH2D[std::make_pair(det.first, "size")];
      m_TH2D[idx]->SetEntries(n);
      idx = m_map_TH2D[std::make_pair(det.first, "charge")];
      m_TH2D[idx]->SetEntries(n);
    } else if (m_config.GetDetectorPlane(dp0) == true) {
      if (m_config.createJSON) {
        int id = m_map_TH1D[std::make_pair(det.first, "cluster")];

        TString jsonFilename = m_fileName;
        jsonFilename.ReplaceAll(".root", "");

        TString json = TBufferJSON::ToJSON(m_TH1D[id], 3);
        std::ofstream f1;
        f1.open(jsonFilename + "_detector" + std::to_string(det.first) +
                    "_cluster.json",
                std::ios::out);
        f1 << json;
        f1.close();
      }
      int n = 0;
      for (int b = 1; b <= m_bins0; b++) {
        int idx = m_map_TH1D[std::make_pair(det.first, "cluster")];
        // Number of clusters in bin
        int cnt = m_TH1D[idx]->GetBinContent(b);
        double val = 0;
        if (cnt > 0) {
          n++;
          idx = m_map_TH1D[std::make_pair(det.first, "size")];
          // Size divided by cnt
          val = m_TH1D[idx]->GetBinContent(b) / cnt;
          m_TH1D[idx]->SetBinContent(b, val);

          idx = m_map_TH1D[std::make_pair(det.first, "charge")];
          val = m_TH1D[idx]->GetBinContent(b) / cnt;
          m_TH1D[idx]->SetBinContent(b, val);
        }
      }

      int idx = m_map_TH1D[std::make_pair(det.first, "size")];
      m_TH1D[idx]->SetEntries(n);
      m_TH1D[idx]->Sumw2(false);
      idx = m_map_TH1D[std::make_pair(det.first, "charge")];
      m_TH1D[idx]->SetEntries(n);
      m_TH1D[idx]->Sumw2(false);
    }
  }
}
