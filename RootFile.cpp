#include "TBufferJSON.h"
#include <TStyle.h>
#include "TMath.h"
#include <time.h>
#include "RootFile.h"

RootFile* RootFile::m_rootFile = nullptr;

RootFile* RootFile::GetInstance() {
    return m_rootFile;
}

RootFile* RootFile::GetInstance(TString fileName, std::map<uint8_t, uint8_t> dets, std::map<std::pair<uint8_t,uint8_t>, uint32_t> channels, bool analyzeChannels, bool createHits) {
    if (!m_rootFile) {
        m_rootFile = new RootFile(fileName, dets, channels, analyzeChannels, createHits);
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
    SaveHistograms();
    m_file->Write("", TObject::kOverwrite);
    m_file->Close();
}

RootFile::RootFile(TString fileName, std::map<uint8_t, uint8_t> dets, std::map<std::pair<uint8_t,uint8_t>, uint32_t> channels,  bool analyzeChannels,  bool createHits) :
m_fileName(fileName), m_dets(dets), m_channels(channels), m_analyzeChannels(analyzeChannels), m_createHits(createHits) {

    m_file = TFile::Open(m_fileName, "RECREATE");
    m_eventNr = 0;

    m_tree = new TTree("events", "vmm3 events");
    m_tree->SetDirectory(m_file);
    if (analyzeChannels) {
        CreateChannelHistograms();
    }


    if (m_createHits) {
        m_tree->Branch("hits", &m_hits);
    }
    m_tree->Branch("cluster", &m_clusters);
    m_tree->Branch("clusterXY", &m_clustersXY);


    for(auto const & det: m_dets)
    {
        std::string name = std::to_string(det.first) + "_clusterXY";
        int x = (m_channels[std::make_pair(det.first, 0)]+1) *64;
        int y = (m_channels[std::make_pair(det.first, 1)]+1) *64;
        m_TH2D_clusterXY[det.second] = new TH2D(name.c_str(), name.c_str(), x* 8, 0, x,y*8,0,y);
        name = std::to_string(det.first) + "_clusterXYDeltaTime";
        m_TH1D_clusterXYDeltaTime[det.second] = new TH1D(name.c_str(), name.c_str(), 2000, -1000, 1000);
        name = std::to_string(det.first) + "_imageXY";
        m_TH2D_imageXY[det.second] = new TH2D(name.c_str(), name.c_str(), x, 0, x,y,0,y);
        name = std::to_string(det.first) + "_chargeX";
        m_TH2D_chargeX[det.second] = new TH2D(name.c_str(), name.c_str(), x, 0, x,y,0,y);
        name = std::to_string(det.first) + "_chargeY";
        m_TH2D_chargeY[det.second] =new TH2D(name.c_str(), name.c_str(), x, 0, x,y,0,y);
        name = std::to_string(det.first) + "_chargeXY";
        m_TH2D_chargeXY[det.second] = new TH2D(name.c_str(), name.c_str(), x, 0, x,y,0,y);
        name = std::to_string(det.first) + "_sizeX";
        m_TH2D_sizeX[det.second] = new TH2D(name.c_str(), name.c_str(), x, 0, x,y,0,y);
        name = std::to_string(det.first) + "_sizeY";
        m_TH2D_sizeY[det.second] = new TH2D(name.c_str(), name.c_str(), x, 0, x,y,0,y);
        name = std::to_string(det.first) + "_sizeXY";
        m_TH2D_sizeXY[det.second] = new TH2D(name.c_str(), name.c_str(), x, 0, x,y,0,y);

    }



    std::cout << "ROOT file " << m_fileName << " created!" << std::endl;

}

RootFile::~RootFile() {

}

void RootFile::FillTree() {
    if (m_hits.size() > 0 || m_clusters.size() > 0) {
        m_tree->Fill();
        m_hits.clear();
        m_clusters.clear();
        m_clustersXY.clear();
    }
}


void RootFile::SaveHits(HitNMX&& theHit) {
    m_hits.emplace_back(theHit);
}



void RootFile::SaveClusters(ClusterVector&& clusters) {
    m_clusters = clusters;
}



void RootFile::SaveClustersXY(CommonClusterVector&& clustersXY) {
    m_clustersXY = clustersXY;
    for (auto & it : m_clustersXY) {
        m_TH2D_clusterXY[m_dets[it.detID]]->Fill(it.positionX, it.positionY);
        m_TH2D_imageXY[m_dets[it.detID]]->Fill(it.positionX, it.positionY);
        m_TH1D_clusterXYDeltaTime[m_dets[it.detID]]->Fill(it.timeX - it.timeY);
        m_TH2D_chargeX[m_dets[it.detID]]->Fill(it.positionX, it.positionY, it.adcX);
        m_TH2D_chargeY[m_dets[it.detID]]->Fill(it.positionX, it.positionY, it.adcY);
        m_TH2D_chargeXY[m_dets[it.detID]]->Fill(it.positionX, it.positionY, it.adcX + it.adcY);
        m_TH2D_sizeX[m_dets[it.detID]]->Fill(it.positionX, it.positionY, it.sizeX);
        m_TH2D_sizeY[m_dets[it.detID]]->Fill(it.positionX, it.positionY, it.sizeY);
        m_TH2D_sizeXY[m_dets[it.detID]]->Fill(it.positionX, it.positionY, it.sizeX + it.sizeY);
    }
}

void RootFile::SaveHistograms() {
    if (m_analyzeChannels) {
        AnalyzeChannels();
    }
    for(auto const & det: m_dets)
    {
        m_TH2D_clusterXY[det.second]->Write("", TObject::kOverwrite);
        m_TH2D_imageXY[det.second]->Write("", TObject::kOverwrite);
        m_TH1D_clusterXYDeltaTime[det.second]->Write("", TObject::kOverwrite);

        TString jsonFilename = m_fileName;
        jsonFilename.ReplaceAll(".root", "");

        TString json = TBufferJSON::ToJSON(m_TH2D_imageXY[det.second], 3);
        std::ofstream f1;
        f1.open(jsonFilename + "_imageXY.json", std::ios::out);
        f1 << json;
        f1.close();

        json = TBufferJSON::ToJSON(m_TH2D_clusterXY[det.second], 3);
        std::ofstream f2;
        f2.open(jsonFilename + "_clusterXY.json", std::ios::out);
        f2 << json;
        f2.close();

        int xm = (m_channels[std::make_pair(det.first, 0)]+1) *64;
        int ym = (m_channels[std::make_pair(det.first, 1)]+1) *64;

        for (int x = 1; x <= xm; x++) {
            for (int y = 1; y <= ym; y++) {
                int cnt = m_TH2D_imageXY[det.second]->GetBinContent(x, y);
                if (cnt > 0) {
                    double val = m_TH2D_chargeX[det.second]->GetBinContent(x, y) / cnt;
                    m_TH2D_chargeX[det.second]->SetBinContent(x, y, val);

                    val = m_TH2D_chargeY[det.second]->GetBinContent(x, y) / cnt;
                    m_TH2D_chargeY[det.second]->SetBinContent(x, y, val);

                    val = m_TH2D_chargeXY[det.second]->GetBinContent(x, y) / cnt;
                    m_TH2D_chargeXY[det.second]->SetBinContent(x, y, val);

                    val = m_TH2D_sizeX[det.second]->GetBinContent(x, y) / cnt;
                    m_TH2D_sizeX[det.second]->SetBinContent(x, y, val);

                    val = m_TH2D_sizeY[det.second]->GetBinContent(x, y) / cnt;
                    m_TH2D_sizeY[det.second]->SetBinContent(x, y, val);

                    val = m_TH2D_sizeXY[det.second]->GetBinContent(x, y) / cnt;
                    m_TH2D_sizeXY[det.second]->SetBinContent(x, y, val);
                }

            }
        }
        m_TH2D_chargeX[det.second]->Write("", TObject::kOverwrite);
        m_TH2D_chargeY[det.second]->Write("", TObject::kOverwrite);
        m_TH2D_chargeXY[det.second]->Write("", TObject::kOverwrite);

        m_TH2D_sizeX[det.second]->Write("", TObject::kOverwrite);
        m_TH2D_sizeY[det.second]->Write("", TObject::kOverwrite);
        m_TH2D_sizeXY[det.second]->Write("", TObject::kOverwrite);
    }
}

void RootFile::CreateChannelHistograms() {
    for(auto const & det: m_dets)
    {
        int x = (m_channels[std::make_pair(det.first, 0)]+1) *64;
        int y = (m_channels[std::make_pair(det.first, 1)]+1) *64;

        std::string name = std::to_string(det.first) + "_tdc_min_x";
        m_tdc_min_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);

        name = std::to_string(det.first) + "_tdc_max_x";
        m_tdc_max_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);

        name = std::to_string(det.first) + "_tdc_mean_x";
        m_tdc_mean_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);

        name = std::to_string(det.first) + "_tdc_stddev_x";
        m_tdc_stddev_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);

        name = std::to_string(det.first) + "_tdc_range_x";
        m_tdc_range_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);

        name = std::to_string(det.first) + "_tdc_totalrange_x";
        m_tdc_totalrange_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);

        name = std::to_string(det.first) + "_tdc_mean_y";
        m_tdc_mean_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);

        name = std::to_string(det.first) + "_tdc_max_y";
        m_tdc_max_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);

        name = std::to_string(det.first) + "_tdc_min_y";
        m_tdc_min_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);

        name = std::to_string(det.first) + "_tdc_stddev_y";
        m_tdc_stddev_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);

        name = std::to_string(det.first) + "_tdc_range_y";
        m_tdc_range_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);

        name = std::to_string(det.first) + "_tdc_totalrange_y";
        m_tdc_totalrange_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);

        name = std::to_string(det.first) + "_adc_mean_x";
        m_adc_mean_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);

        name = std::to_string(det.first) + "_adc_max_x";
        m_adc_max_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);

        name = std::to_string(det.first) + "_adc_min_x";
        m_adc_min_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);

        name = std::to_string(det.first) + "_adc_stddev_x";
        m_adc_stddev_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);

        name = std::to_string(det.first) + "_adc_range_x";
        m_adc_range_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);

        name = std::to_string(det.first) + "_adc_totalrange_x";
        m_adc_totalrange_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);

        name = std::to_string(det.first) + "_adc_mean_y";
        m_adc_mean_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);

        name = std::to_string(det.first) + "_adc_max_y";
        m_adc_max_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);

        name = std::to_string(det.first) + "_adc_min_y";
        m_adc_min_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);

        name = std::to_string(det.first) + "_adc_stddev_y";
        m_adc_stddev_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);

        name = std::to_string(det.first) + "_adc_range_y";
        m_adc_range_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);
        name = std::to_string(det.first) + "_adc_totalrange_y";
        m_adc_totalrange_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);
        name = std::to_string(det.first) + "_bcid_mean_x";
        m_bcid_mean_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);
        name = std::to_string(det.first) + "_bcid_max_x";
        m_bcid_max_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);
        name = std::to_string(det.first) + "_bcid_min_x";
        m_bcid_min_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);
        name = std::to_string(det.first) + "_bcid_stddev_x";
        m_bcid_stddev_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);
        name = std::to_string(det.first) + "_bcid_range_x";
        m_bcid_range_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);
        name = std::to_string(det.first) + "_bcid_totalrange_x";
        m_bcid_totalrange_x[det.second] = new TH1D(name.c_str(), name.c_str(), x,0,x);
        name = std::to_string(det.first) + "_bcid_mean_y";
        m_bcid_mean_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);
        name = std::to_string(det.first) + "_bcid_max_y";
        m_bcid_max_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);
        name = std::to_string(det.first) + "_bcid_min_y";
        m_bcid_min_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);
        name = std::to_string(det.first) + "_bcid_stddev_y";
        m_bcid_stddev_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);
        name = std::to_string(det.first) + "_bcid_range_y";
        m_bcid_range_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);
        name = std::to_string(det.first) + "_bcid_totalrange_y";
        m_bcid_totalrange_y[det.second] = new TH1D(name.c_str(), name.c_str(), y,0,y);
    }
}

void RootFile::AnalyzeChannels() {
    for(auto const & det: m_dets)
    {
        int x = (m_channels[std::make_pair(det.first, 0)]+1) *64;
        int y = (m_channels[std::make_pair(det.first, 1)]+1) *64;
        std::cout << "Detector " << (int)det.first << std::endl;
        std::cout << "Creating TDC histograms X!" << std::endl;
        for (int i = 0; i < x; i++) {
                TH1D *h1 = new TH1D("h1", "h1", 256, 0, 256);
                TString cut = "hits.position==" + std::to_string(i) + " && hits.detID ==" + std::to_string(det.first) + " && hits.planeID == " + std::to_string(1);
                m_tree->Draw("hits.tdc>>h1", cut, "goff");
                double mean = h1->GetMean();

                double stddev = h1->GetStdDev();
                double max = h1->FindLastBinAbove(0);

                double min = h1->FindFirstBinAbove(0);

                double totalrange = max - min;
                m_tdc_mean_x[det.second]->Fill(i, mean);
                m_tdc_min_x[det.second]->Fill(i, min);
                m_tdc_max_x[det.second]->Fill(i, max);
                m_tdc_stddev_x[det.second]->Fill(i, stddev);
                m_tdc_totalrange_x[det.second]->Fill(i, totalrange);
                double start = mean - stddev;
                if (start < 0) {
                        start = 0;
                }
                double end = mean + stddev;
                if (end > 255) {
                        end = 255;
                }
                m_tdc_range_x[det.second]->Fill(i, end - start);
                delete h1;
        }
        std::cout << "Creating TDC histograms Y!" << std::endl;
        for (int i = 0; i < y; i++) {
                TH1D *h1 = new TH1D("h1", "h1", 256, 0, 256);
                TString cut = "hits.position==" + std::to_string(i) + " && hits.detID ==" + std::to_string(det.first) + " && hits.planeID == " + std::to_string(1);
                m_tree->Draw("hits.tdc>>h1", cut, "goff");
                double mean = h1->GetMean();
                double stddev = h1->GetStdDev();
                double max = h1->FindLastBinAbove(0);
                double min = h1->FindFirstBinAbove(0);
                double totalrange = max - min;
                m_tdc_mean_y[det.second]->Fill(i, mean);
                m_tdc_min_y[det.second]->Fill(i, min);
                m_tdc_max_y[det.second]->Fill(i, max);
                m_tdc_stddev_y[det.second]->Fill(i, stddev);
                m_tdc_totalrange_y[det.second]->Fill(i, totalrange);

                double start = mean - stddev;
                if (start < 0) {
                        start = 0;
                }
                double end = mean + stddev;
                if (end > 255) {
                        end = 255;
                }
                m_tdc_range_y[det.second]->Fill(i, end - start);

                delete h1;
        }

        std::cout << "Creating ADC histograms X!" << std::endl;
        for (int i = 0; i < x; i++) {
                TH1D *h1 = new TH1D("h1", "h1", 1024, 0, 1024);
                TString cut = "hits.position==" + std::to_string(i) + " && hits.detID ==" + std::to_string(det.first) + " && hits.planeID == " + std::to_string(0);
                m_tree->Draw("hits.adc>>h1", cut, "goff");
                double mean = h1->GetMean();
                double stddev = h1->GetStdDev();
                double max = h1->FindLastBinAbove(0);
                double min = h1->FindFirstBinAbove(0);
                double totalrange = max - min;
                m_adc_mean_x[det.second]->Fill(i, mean);
                m_adc_min_x[det.second]->Fill(i, min);
                m_adc_max_x[det.second]->Fill(i, max);
                m_adc_stddev_x[det.second]->Fill(i, stddev);
                m_adc_totalrange_x[det.second]->Fill(i, totalrange);
                double start = mean - stddev;
                if (start < 0) {
                        start = 0;
                }
                double end = mean + stddev;
                if (end > 1023) {
                        end = 1023;
                }
                m_adc_range_x[det.second]->Fill(i, end - start);
                delete h1;

        }
        std::cout << "Creating ADC histograms Y!" << std::endl;
        for (int i = 0; i < y; i++) {
                TH1D *h1 = new TH1D("h1", "h1", 1024, 0, 1024);
                TString cut = "hits.position==" + std::to_string(i) + " && hits.detID ==" + std::to_string(det.first) + " && hits.planeID == " + std::to_string(1);

                m_tree->Draw("hits.adc>>h1", cut, "goff");
                double mean = h1->GetMean();
                double stddev = h1->GetStdDev();
                double max = h1->FindLastBinAbove(0);
                double min = h1->FindFirstBinAbove(0);
                double totalrange = max - min;
                m_adc_mean_y[det.second]->Fill(i, mean);
                m_adc_min_y[det.second]->Fill(i, min);
                m_adc_max_y[det.second]->Fill(i, max);
                m_adc_stddev_y[det.second]->Fill(i, stddev);
                m_adc_totalrange_y[det.second]->Fill(i, totalrange);

                double start = mean - stddev;
                if (start < 0) {
                        start = 0;
                }
                double end = mean + stddev;
                if (end > 1023) {
                        end = 1023;
                }
                m_adc_range_y[det.second]->Fill(i, end - start);

                delete h1;

        }
        std::cout << "Creating BCID histograms X!" << std::endl;
        for (int i = 0; i < x; i++) {
                TH1D *h1 = new TH1D("h1", "h1", 4096, 0, 4096);
                TString cut = "hits.position==" + std::to_string(i) + " && hits.detID ==" + std::to_string(det.first) + " && hits.planeID == " + std::to_string(0);

                m_tree->Draw("hits.bcid>>h1", cut, "goff");
                double mean = h1->GetMean();
                double stddev = h1->GetStdDev();
                double max = h1->FindLastBinAbove(0);
                double min = h1->FindFirstBinAbove(0);
                double totalrange = max - min;
                m_bcid_mean_x[det.second]->Fill(i, mean);
                m_bcid_min_x[det.second]->Fill(i, min);
                m_bcid_max_x[det.second]->Fill(i, max);
                m_bcid_stddev_x[det.second]->Fill(i, stddev);
                m_bcid_totalrange_x[det.second]->Fill(i, totalrange);
                double start = mean - stddev;
                if (start < 0) {
                        start = 0;
                }
                double end = mean + stddev;
                if (end > 4095) {
                        end = 4095;
                }
                m_bcid_range_x[det.second]->Fill(i, end - start);
                delete h1;
        }
        std::cout << "Creating BCID histograms Y!" << std::endl;
        for (int i = 0; i < y; i++) {
                TH1D *h1 = new TH1D("h1", "h1", 4096, 0, 4096);
                TString cut = "hits.position==" + std::to_string(i) + " && hits.detID ==" + std::to_string(det.first) + " && hits.planeID == " + std::to_string(1);
                m_tree->Draw("hits.bcid>>h1", cut, "goff");
                double mean = h1->GetMean();
                double stddev = h1->GetStdDev();
                double max = h1->FindLastBinAbove(0);
                double min = h1->FindFirstBinAbove(0);
                double totalrange = max - min;
                m_bcid_mean_y[det.second]->Fill(i, mean);
                m_bcid_min_y[det.second]->Fill(i, min);
                m_bcid_max_y[det.second]->Fill(i, max);
                m_bcid_stddev_y[det.second]->Fill(i, stddev);
                m_bcid_totalrange_y[det.second]->Fill(i, totalrange);

                double start = mean - stddev;
                if (start < 0) {
                        start = 0;
                }
                double end = mean + stddev;
                if (end > 4095) {
                        end = 4095;
                }
                m_bcid_range_y[det.second]->Fill(i, end - start);

                delete h1;
        }

        m_bcid_mean_x[det.second]->Write("", TObject::kOverwrite);
        m_bcid_min_x[det.second]->Write("", TObject::kOverwrite);
        m_bcid_max_x[det.second]->Write("", TObject::kOverwrite);
        m_bcid_stddev_x[det.second]->Write("", TObject::kOverwrite);
        m_bcid_totalrange_x[det.second]->Write("", TObject::kOverwrite);
        m_bcid_range_x[det.second]->Write("", TObject::kOverwrite);
        m_bcid_mean_y[det.second]->Write("", TObject::kOverwrite);
        m_bcid_min_y[det.second]->Write("", TObject::kOverwrite);
        m_bcid_max_y[det.second]->Write("", TObject::kOverwrite);
        m_bcid_stddev_y[det.second]->Write("", TObject::kOverwrite);
        m_bcid_totalrange_y[det.second]->Write("", TObject::kOverwrite);
        m_bcid_range_y[det.second]->Write("", TObject::kOverwrite);
    }
}
