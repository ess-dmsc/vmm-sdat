#include "TBufferJSON.h"
#include <TStyle.h>
#include "TMath.h"
#include <time.h>
#include "RootFile.h"

RootFile *RootFile::m_rootFile = nullptr;

RootFile *RootFile::GetInstance()
{
    return m_rootFile;
}

RootFile *RootFile::GetInstance(Configuration &config)
{
    if (!m_rootFile)
    {
        m_rootFile = new RootFile(config);
    }

    return m_rootFile;
}

void RootFile::Dispose()
{
    if (m_rootFile)
    {
        m_rootFile->WriteRootFile();
        delete m_rootFile;
        m_rootFile = 0;
    }
}

void RootFile::WriteRootFile()
{
    SaveHistograms();
    m_file->Write("", TObject::kOverwrite);
    m_file->Close();
}

RootFile::RootFile(Configuration &config) : m_config(config)
{
    m_fileName = m_config.pRootFilename.c_str();
    m_file = TFile::Open(m_fileName, "RECREATE");
    m_eventNr = 0;

    m_tree = new TTree("events", "vmm3 events");
    m_tree->SetDirectory(m_file);

    if (m_config.pCreateHits)
    {
        m_tree->Branch("hits", &m_hits);
    }
    m_tree->Branch("clusters_plane", &m_clusters_plane);
    m_tree->Branch("clusters_detector", &m_clusters_detector);
    m_tree->Branch("tracks", &m_tracks);

    TH2D *h2;
    TH1D *h1;
    std::string name = "";
    int cnt1D = 0;
    int cnt2D = 0;
    for (auto const &det : m_config.pDets)
    {
        name = std::to_string(det.first) + "_delta_time_planes";
        h1 = new TH1D(name.c_str(), name.c_str(), 1000, -500, 500);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(std::make_pair(std::make_pair(det.first, "delta_time_planes"), cnt1D));
        cnt1D++;
        
        name = std::to_string(det.first) + "_delta_time_utpc_planes";
        h1 = new TH1D(name.c_str(), name.c_str(), 1000, -500, 500);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(std::make_pair(std::make_pair(det.first, "delta_time_utpc_planes"), cnt1D));
        cnt1D++;

        name = std::to_string(det.first) + "_delta_time_charge2_planes";
        h1 = new TH1D(name.c_str(), name.c_str(), 1000, -500, 500);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(std::make_pair(std::make_pair(det.first, "delta_time_charge2_planes"), cnt1D));
        cnt1D++;

        name = std::to_string(det.first) + "_dt0";
        h1 = new TH1D(name.c_str(), name.c_str(), 11000, -500, 50000);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(std::make_pair(std::make_pair(det.first, "dt0"), cnt1D));
        cnt1D++;

        name = std::to_string(det.first) + "_dt1";
        h1 = new TH1D(name.c_str(), name.c_str(), 11000, -500, 50000);
        m_TH1D.push_back(h1);
        m_map_TH1D.emplace(std::make_pair(std::make_pair(det.first, "dt1"), cnt1D));
        cnt1D++;

        int n0 = m_config.pChannels[std::make_pair(det.first, 0)];
        int n1 = m_config.pChannels[std::make_pair(det.first, 1)];
        name = std::to_string(det.first) + "_cluster";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, 0, n0, n1 * 4, 0, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(std::make_pair(std::make_pair(det.first, "cluster"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_cluster_utpc";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, 0, n0, n1 * 4, 0, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(std::make_pair(std::make_pair(det.first, "cluster_utpc"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_cluster_charge2";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, 0, n0, n1 * 4, 0, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(std::make_pair(std::make_pair(det.first, "cluster_charge2"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_size_plane0";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, 0, n0, n1 * 4, 0, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(std::make_pair(std::make_pair(det.first, "size_plane0"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_size_plane1";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, 0, n0, n1 * 4, 0, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(std::make_pair(std::make_pair(det.first, "size_plane1"), cnt2D));
        cnt2D++;
        
        name = std::to_string(det.first) + "_size_plane01";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, 0, n0, n1 * 4, 0, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(std::make_pair(std::make_pair(det.first, "size_plane01"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_charge_plane0";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, 0, n0, n1 * 4, 0, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(std::make_pair(std::make_pair(det.first, "charge_plane0"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_charge_plane1";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, 0, n0, n1 * 4, 0, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(std::make_pair(std::make_pair(det.first, "charge_plane1"), cnt2D));
        cnt2D++;

        name = std::to_string(det.first) + "_charge_plane01";
        h2 = new TH2D(name.c_str(), name.c_str(), n0 * 4, 0, n0, n1 * 4, 0, n1);
        m_TH2D.push_back(h2);
        m_map_TH2D.emplace(std::make_pair(std::make_pair(det.first, "charge_plane01"), cnt2D));
        cnt2D++;

    }

    std::cout << "ROOT file " << m_fileName << " created!" << std::endl;
}

RootFile::~RootFile()
{
}


void RootFile::AddHits(Hit &&the_hit)
{
    m_hits.emplace_back(the_hit);
}

void RootFile::SaveHits()
{
    if (m_hits.size() > 0)
    {
        m_tree->Fill();
        m_hits.clear();
    }
}

void RootFile::SaveClustersPlane(ClusterVectorPlane &&clusters_plane)
{
    m_clusters_plane = clusters_plane;
    if (m_clusters_plane.size() > 0)
    {
        m_tree->Fill();
        m_clusters_plane.clear();
    }
}

void RootFile::SaveClustersDetector(ClusterVectorDetector &&clusters_detector)
{
    m_clusters_detector = clusters_detector;
    for (auto &it : m_clusters_detector)
    {
        int idx = m_map_TH1D[std::make_pair(it.det,"delta_time_planes")];
        m_TH1D[idx]->Fill(it.time0 - it.time1);

        idx = m_map_TH1D[std::make_pair(it.det,"delta_time_utpc_planes")];
        m_TH1D[idx]->Fill(it.time0_utpc - it.time1_utpc);

        idx = m_map_TH1D[std::make_pair(it.det,"delta_time_charge2_planes")];
        m_TH1D[idx]->Fill(it.time0_charge2 - it.time1_charge2);

        idx = m_map_TH1D[std::make_pair(it.det,"dt0")];
        m_TH1D[idx]->Fill(it.dt0);

        idx = m_map_TH1D[std::make_pair(it.det,"dt1")];
        m_TH1D[idx]->Fill(it.dt1);

        idx = m_map_TH2D[std::make_pair(it.det,"cluster")];
        m_TH2D[idx]->Fill(it.pos0, it.pos1);
        
        idx = m_map_TH2D[std::make_pair(it.det,"cluster_utpc")];
        m_TH2D[idx]->Fill(it.pos0_utpc, it.pos1_utpc);

        idx = m_map_TH2D[std::make_pair(it.det,"cluster_charge2")];
        m_TH2D[idx]->Fill(it.pos0_charge2, it.pos1_charge2);
        
        idx = m_map_TH2D[std::make_pair(it.det,"size_plane0")];
        m_TH2D[idx]->Fill(it.pos0, it.pos1, it.size0);

        idx = m_map_TH2D[std::make_pair(it.det,"size_plane1")];
        m_TH2D[idx]->Fill(it.pos0, it.pos1, it.size1);

        idx = m_map_TH2D[std::make_pair(it.det,"size_plane01")];
        m_TH2D[idx]->Fill(it.pos0, it.pos1, it.size0 + it.size1);

        idx = m_map_TH2D[std::make_pair(it.det,"charge_plane0")];
        m_TH2D[idx]->Fill(it.pos0, it.pos1, it.adc0);

        idx = m_map_TH2D[std::make_pair(it.det,"charge_plane1")];
        m_TH2D[idx]->Fill(it.pos0, it.pos1, it.adc1);

        idx = m_map_TH2D[std::make_pair(it.det,"charge_plane01")];
        m_TH2D[idx]->Fill(it.pos0, it.pos1, it.adc0 + it.adc1);

    }
    if (m_clusters_detector.size() > 0)
    {
        m_tree->Fill();
        m_clusters_detector.clear();
    }
}


void RootFile::SaveHistograms()
{
    for (auto const &h1 : m_TH1D)
    {
        h1->Write("", TObject::kOverwrite);
    }
    for (auto const &det : m_config.pDets)
    {
        if(m_config.createJSON)
        {
            int id = m_map_TH2D[std::make_pair(det.first,"cluster")];
        
            TString jsonFilename = m_fileName;
            jsonFilename.ReplaceAll(".root", "");

            TString json = TBufferJSON::ToJSON(m_TH2D[id], 3);
            std::ofstream f1;
            f1.open(jsonFilename  + "_detector" + std::to_string(det.first) + "_cluster.json", std::ios::out);
            f1 << json;
            f1.close();
        }

        int n0 = m_config.pChannels[std::make_tuple(det.first, 0)]*4;
        int n1 = m_config.pChannels[std::make_tuple(det.first, 1)]*4;
        int n = 0;
        for (int p0 = 1; p0 <= n0; p0++)
        {
            for (int p1 = 1; p1 <= n1; p1++)
            {
                int idx = m_map_TH2D[std::make_pair(det.first,"cluster")];
                int cnt = m_TH2D[idx]->GetBinContent(p0, p1);
                double val = 0;
                if (cnt > 0)
                {
                    n++;
                    idx = m_map_TH2D[std::make_pair(det.first,"size_plane0")];
                    val = m_TH2D[idx]->GetBinContent(p0, p1) / cnt;
                    m_TH2D[idx]->SetBinContent(p0, p1,val);

                    idx = m_map_TH2D[std::make_pair(det.first,"size_plane0")];
                    val = m_TH2D[idx]->GetBinContent(p0, p1) / cnt;
                    m_TH2D[idx]->SetBinContent(p0, p1,val);

                    idx = m_map_TH2D[std::make_pair(det.first,"size_plane01")];
                    val = m_TH2D[idx]->GetBinContent(p0, p1) / cnt;
                    m_TH2D[idx]->SetBinContent(p0, p1,val);

                    idx = m_map_TH2D[std::make_pair(det.first,"charge_plane0")];
                    val = m_TH2D[idx]->GetBinContent(p0, p1) / cnt;
                    m_TH2D[idx]->SetBinContent(p0, p1,val);

                    idx = m_map_TH2D[std::make_pair(det.first,"charge_plane1")];
                    val = m_TH2D[idx]->GetBinContent(p0, p1) / cnt;
                    m_TH2D[idx]->SetBinContent(p0, p1,val);

                    idx = m_map_TH2D[std::make_pair(det.first,"charge_plane01")];
                    val = m_TH2D[idx]->GetBinContent(p0, p1) / cnt;
                    m_TH2D[idx]->SetBinContent(p0, p1,val);

                }
            }
        }
        int idx = m_map_TH2D[std::make_pair(det.first,"size_plane0")];
        m_TH2D[idx]->SetEntries(n);
        idx = m_map_TH2D[std::make_pair(det.first,"size_plane1")];
        m_TH2D[idx]->SetEntries(n);
        idx = m_map_TH2D[std::make_pair(det.first,"size_plane01")];
        m_TH2D[idx]->SetEntries(n);
        idx = m_map_TH2D[std::make_pair(det.first,"charge_plane0")];
        m_TH2D[idx]->SetEntries(n);
        idx = m_map_TH2D[std::make_pair(det.first,"charge_plane1")];
        m_TH2D[idx]->SetEntries(n);
        idx = m_map_TH2D[std::make_pair(det.first,"charge_plane01")];
        m_TH2D[idx]->SetEntries(n);
    }
}

