#include "TBufferJSON.h"
#include <TStyle.h>
#include "TMath.h"
#include <time.h>
#include "RootFile.h"

RootFile* RootFile::m_rootFile = 0;

RootFile* RootFile::GetInstance() {
	return m_rootFile;
}

RootFile* RootFile::GetInstance(std::string fileName, int channels_x, int channels_y, bool analyzeChannels, bool createHits) {
	if (!m_rootFile) {
		m_rootFile = new RootFile(fileName, channels_x, channels_y, analyzeChannels, createHits);
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

RootFile::RootFile(TString fileName, int channels_x, int channels_y, bool analyzeChannels, bool createHits) :
		m_fileName(fileName), m_channels_x(channels_x), m_channels_y(channels_y), m_analyzeChannels(analyzeChannels), m_createHits(createHits) {

	m_file = TFile::Open(m_fileName, "RECREATE");
	m_eventNr = 0;
	m_nch = 0;
	m_nchX = 0;
	m_nchY = 0;

	m_nclX = 0;
	m_nclY = 0;
	m_nclXY = 0;
	m_nclXY_uTPC = 0;

	m_tree = new TTree("events", "vmm3 events");
	m_tree->SetDirectory(m_file);
	if (analyzeChannels) {
		CreateChannelHistograms();
	}

	if (m_createHits) {
		m_tree->Branch("hitsX", &m_hitsX);
		m_tree->Branch("hitsY", &m_hitsY);
	}
	m_tree->Branch("clusterX", &m_clusterX);
	m_tree->Branch("clusterY", &m_clusterY);
	m_tree->Branch("clusterXY", &m_clusterXY);
	m_tree->Branch("clusterXY_uTPC");

	m_TH2D_clusterXY = new TH2D("clusterXY", "clusterXY", m_channels_x * 8, 0, m_channels_x, m_channels_y * 8, 0, m_channels_y);
	m_TH1D_clusterXYDeltaTime = new TH1D("clusterXYDeltaTime", "clusterXYDeltaTime", 2000, -1000, 1000);

	m_TH2D_imageXY = new TH2D("imageXY", "imageXY", m_channels_x, 0, m_channels_x, m_channels_y, 0, m_channels_y);

	m_TH2D_chargeX = new TH2D("chargeX", "chargeX", m_channels_x, 0, m_channels_x, m_channels_y, 0, m_channels_y);
	m_TH2D_chargeY = new TH2D("chargeY", "chargeY", m_channels_x, 0, m_channels_x, m_channels_y, 0, m_channels_y);
	m_TH2D_chargeXY = new TH2D("chargeXY", "chargeXY", m_channels_x, 0, m_channels_x, m_channels_y, 0, m_channels_y);

	m_TH2D_sizeX = new TH2D("sizeX", "sizeX", m_channels_x, 0, m_channels_x, m_channels_y, 0, m_channels_y);
	m_TH2D_sizeY = new TH2D("sizeY", "sizeY", m_channels_x, 0, m_channels_x, m_channels_y, 0, m_channels_y);
	m_TH2D_sizeXY = new TH2D("sizeXY", "sizeXY", m_channels_x, 0, m_channels_x, m_channels_y, 0, m_channels_y);

	std::cout << "ROOT file " << m_fileName << " created!" << std::endl;

}

RootFile::~RootFile() {

}

void RootFile::FillTree() {
	if (m_hitsX.size() > 0 || m_hitsY.size() || m_clusterX.size() > 0 || m_clusterY.size() > 0) {
		m_tree->Fill();
		m_hitsX.clear();
		m_hitsY.clear();
		m_clusterX.clear();
		m_clusterY.clear();
		m_clusterXY.clear();
	}
}

void RootFile::SaveHitsX(HitNMX&& theHit) {
	m_hitsX.emplace_back(theHit);
}

void RootFile::SaveHitsY(HitNMX&& theHit) {
	m_hitsY.emplace_back(theHit);
}

void RootFile::SaveClustersX(ClusterVector&& clusters) {
	m_clusterX = clusters;
}

void RootFile::SaveClustersY(ClusterVector&& clusters) {
	m_clusterY = clusters;
}

void RootFile::SaveClustersXY(CommonClusterVector&& clustersXY) {
	m_clusterXY = clustersXY;
	for (auto & it : m_clusterXY) {
		m_TH2D_clusterXY->Fill(it.positionX, it.positionY);
		m_TH2D_imageXY->Fill(it.positionX, it.positionY);
		m_TH1D_clusterXYDeltaTime->Fill(it.timeX - it.timeY);
		m_TH2D_chargeX->Fill(it.positionX, it.positionY, it.adcX);
		m_TH2D_chargeY->Fill(it.positionX, it.positionY, it.adcY);
		m_TH2D_chargeXY->Fill(it.positionX, it.positionY, it.adcX + it.adcY);
		m_TH2D_sizeX->Fill(it.positionX, it.positionY, it.sizeX);
		m_TH2D_sizeY->Fill(it.positionX, it.positionY, it.sizeY);
		m_TH2D_sizeXY->Fill(it.positionX, it.positionY, it.sizeX + it.sizeY);
	}
}

void RootFile::SaveHistograms() {

	m_TH2D_clusterXY->Write("", TObject::kOverwrite);
	m_TH2D_imageXY->Write("", TObject::kOverwrite);
	m_TH1D_clusterXYDeltaTime->Write("", TObject::kOverwrite);
	if (m_analyzeChannels) {
		AnalyzeChannels();
	}
	TString jsonFilename = m_fileName;
	jsonFilename.ReplaceAll(".root", "");

	TString json = TBufferJSON::ToJSON(m_TH2D_imageXY, 3);
	std::ofstream f1;
	f1.open(jsonFilename + "_imageXY.json", std::ios::out);
	f1 << json;
	f1.close();

	json = TBufferJSON::ToJSON(m_TH2D_clusterXY, 3);
	std::ofstream f2;
	f2.open(jsonFilename + "_clusterXY.json", std::ios::out);
	f2 << json;
	f2.close();

	for (int x = 1; x <= 256; x++) {
		for (int y = 1; y <= 256; y++) {
			int cnt = m_TH2D_imageXY->GetBinContent(x, y);
			if (cnt > 0) {
				double val = m_TH2D_chargeX->GetBinContent(x, y) / cnt;
				m_TH2D_chargeX->SetBinContent(x, y, val);

				val = m_TH2D_chargeY->GetBinContent(x, y) / cnt;
				m_TH2D_chargeY->SetBinContent(x, y, val);

				val = m_TH2D_chargeXY->GetBinContent(x, y) / cnt;
				m_TH2D_chargeXY->SetBinContent(x, y, val);

				val = m_TH2D_sizeX->GetBinContent(x, y) / cnt;
				m_TH2D_sizeX->SetBinContent(x, y, val);

				val = m_TH2D_sizeY->GetBinContent(x, y) / cnt;
				m_TH2D_sizeY->SetBinContent(x, y, val);

				val = m_TH2D_sizeXY->GetBinContent(x, y) / cnt;
				m_TH2D_sizeXY->SetBinContent(x, y, val);
			}

		}
	}
	m_TH2D_chargeX->Write("", TObject::kOverwrite);
	m_TH2D_chargeY->Write("", TObject::kOverwrite);
	m_TH2D_chargeXY->Write("", TObject::kOverwrite);

	m_TH2D_sizeX->Write("", TObject::kOverwrite);
	m_TH2D_sizeY->Write("", TObject::kOverwrite);
	m_TH2D_sizeXY->Write("", TObject::kOverwrite);

}

void RootFile::CreateChannelHistograms() {
	m_tdc_mean_x = new TH1D("tdc_mean_x", "tdc_mean_x", m_channels_x, 0, m_channels_x);
	m_tdc_max_x = new TH1D("tdc_max_x", "tdc_max_x", m_channels_x, 0, m_channels_x);
	m_tdc_min_x = new TH1D("tdc_min_x", "tdc_min_x", m_channels_x, 0, m_channels_x);
	m_tdc_stddev_x = new TH1D("tdc_stddev_x", "tdc_stddev_x", m_channels_x, 0, m_channels_x);
	m_tdc_range_x = new TH1D("tdc_range_x", "tdc_range_x", m_channels_x, 0, m_channels_x);
	m_tdc_totalrange_x = new TH1D("tdc_totalrange_x", "tdc_totalrange_x", m_channels_x, 0, m_channels_x);

	m_tdc_mean_y = new TH1D("tdc_mean_y", "tdc_mean_y", m_channels_y, 0, m_channels_y);
	m_tdc_max_y = new TH1D("tdc_max_y", "tdc_max_y", m_channels_y, 0, m_channels_y);
	m_tdc_min_y = new TH1D("tdc_min_y", "tdc_min_y", m_channels_y, 0, m_channels_y);
	m_tdc_stddev_y = new TH1D("tdc_stddev_y", "tdc_stddev_y", m_channels_y, 0, m_channels_y);
	m_tdc_range_y = new TH1D("tdc_range_y", "tdc_range_y", m_channels_y, 0, m_channels_y);
	m_tdc_totalrange_y = new TH1D("tdc_totalrange_y", "tdc_totalrange_y", m_channels_y, 0, m_channels_y);

	m_adc_mean_x = new TH1D("adc_mean_x", "adc_mean_x", m_channels_x, 0, m_channels_x);
	m_adc_max_x = new TH1D("adc_max_x", "adc_max_x", m_channels_x, 0, m_channels_x);
	m_adc_min_x = new TH1D("adc_min_x", "adc_min_x", m_channels_x, 0, m_channels_x);
	m_adc_stddev_x = new TH1D("adc_stddev_x", "adc_stddev_x", m_channels_x, 0, m_channels_x);
	m_adc_range_x = new TH1D("adc_range_x", "adc_range_x", m_channels_x, 0, m_channels_x);
	m_adc_totalrange_x = new TH1D("adc_totalrange_x", "adc_totalrange_x", m_channels_x, 0, m_channels_x);

	m_adc_mean_y = new TH1D("adc_mean_y", "adc_mean_y", m_channels_y, 0, m_channels_y);
	m_adc_max_y = new TH1D("adc_max_y", "adc_max_y", m_channels_y, 0, m_channels_y);
	m_adc_min_y = new TH1D("adc_min_y", "adc_min_y", m_channels_y, 0, m_channels_y);
	m_adc_stddev_y = new TH1D("adc_stddev_y", "adc_stddev_y", m_channels_y, 0, m_channels_y);
	m_adc_range_y = new TH1D("adc_range_y", "adc_range_y", m_channels_y, 0, m_channels_y);
	m_adc_totalrange_y = new TH1D("adc_totalrange_y", "adc_totalrange_y", m_channels_y, 0, m_channels_y);

	m_bcid_mean_x = new TH1D("bcid_mean_x", "bcid_mean_x", m_channels_x, 0, m_channels_x);
	m_bcid_max_x = new TH1D("bcid_max_x", "bcid_max_x", m_channels_x, 0, m_channels_x);
	m_bcid_min_x = new TH1D("bcid_min_x", "bcid_min_x", m_channels_x, 0, m_channels_x);
	m_bcid_stddev_x = new TH1D("bcid_stddev_x", "bcid_stddev_x", m_channels_x, 0, m_channels_x);
	m_bcid_range_x = new TH1D("bcid_range_x", "bcid_range_x", m_channels_x, 0, m_channels_x);
	m_bcid_totalrange_x = new TH1D("bcid_totalrange_x", "bcid_totalrange_x", m_channels_x, 0, m_channels_x);

	m_bcid_mean_y = new TH1D("bcid_mean_y", "bcid_mean_y", m_channels_y, 0, m_channels_y);
	m_bcid_max_y = new TH1D("bcid_max_y", "bcid_max_y", m_channels_y, 0, m_channels_y);
	m_bcid_min_y = new TH1D("bcid_min_y", "bcid_min_y", m_channels_y, 0, m_channels_y);
	m_bcid_stddev_y = new TH1D("bcid_stddev_y", "bcid_stddev_y", m_channels_y, 0, m_channels_y);
	m_bcid_range_y = new TH1D("bcid_range_y", "bcid_range_y", m_channels_y, 0, m_channels_y);
	m_bcid_totalrange_y = new TH1D("bcid_totalrange_y", "bcid_totalrange_y", m_channels_y, 0, m_channels_y);
}

void RootFile::AnalyzeChannels() {
	std::cout << "Creating TDC histograms X!" << std::endl;
	for (int i = 0; i < m_channels_x; i++) {
		TH1D *h1 = new TH1D("h1", "h1", 256, 0, 256);
		TString cut = "hitsX.position==" + std::to_string(i);
		m_tree->Draw("hitsX.tdc>>h1", cut, "goff");
		double mean = h1->GetMean();

		double stddev = h1->GetStdDev();
		double max = h1->FindLastBinAbove(0);

		double min = h1->FindFirstBinAbove(0);

		double totalrange = max - min;
		m_tdc_mean_x->Fill(i, mean);
		m_tdc_min_x->Fill(i, min);
		m_tdc_max_x->Fill(i, max);
		m_tdc_stddev_x->Fill(i, stddev);
		m_tdc_totalrange_x->Fill(i, totalrange);
		double start = mean - stddev;
		if (start < 0) {
			start = 0;
		}
		double end = mean + stddev;
		if (end > 255) {
			end = 255;
		}
		m_tdc_range_x->Fill(i, end - start);
		delete h1;
	}
	std::cout << "Creating TDC histograms Y!" << std::endl;
	for (int i = 0; i < m_channels_y; i++) {
		TH1D *h1 = new TH1D("h1", "h1", 256, 0, 256);
		TString cut = "hitsY.position==" + std::to_string(i);
		m_tree->Draw("hitsY.tdc>>h1", cut, "goff");
		double mean = h1->GetMean();
		double stddev = h1->GetStdDev();
		double max = h1->FindLastBinAbove(0);
		double min = h1->FindFirstBinAbove(0);
		double totalrange = max - min;
		m_tdc_mean_y->Fill(i, mean);
		m_tdc_min_y->Fill(i, min);
		m_tdc_max_y->Fill(i, max);
		m_tdc_stddev_y->Fill(i, stddev);
		m_tdc_totalrange_y->Fill(i, totalrange);

		double start = mean - stddev;
		if (start < 0) {
			start = 0;
		}
		double end = mean + stddev;
		if (end > 255) {
			end = 255;
		}
		m_tdc_range_y->Fill(i, end - start);

		delete h1;
	}

	std::cout << "Creating ADC histograms!" << std::endl;
	for (int i = 0; i < m_channels_x; i++) {
		TH1D *h1 = new TH1D("h1", "h1", 1024, 0, 1024);
		TString cut = "hitsX.position==" + std::to_string(i);
		m_tree->Draw("hitsX.bcid>>h1", cut, "goff");
		double mean = h1->GetMean();
		double stddev = h1->GetStdDev();
		double max = h1->FindLastBinAbove(0);
		double min = h1->FindFirstBinAbove(0);
		double totalrange = max - min;
		m_adc_mean_x->Fill(i, mean);
		m_adc_min_x->Fill(i, min);
		m_adc_max_x->Fill(i, max);
		m_adc_stddev_x->Fill(i, stddev);
		m_adc_totalrange_x->Fill(i, totalrange);
		double start = mean - stddev;
		if (start < 0) {
			start = 0;
		}
		double end = mean + stddev;
		if (end > 1023) {
			end = 1023;
		}
		m_adc_range_x->Fill(i, end - start);
		delete h1;

	}

	for (int i = 0; i < m_channels_y; i++) {
		TH1D *h1 = new TH1D("h1", "h1", 1024, 0, 1024);
		TString cut = "hitsY.position==" + std::to_string(i);
		m_tree->Draw("hitsY.bcid>>h1", cut, "goff");
		double mean = h1->GetMean();
		double stddev = h1->GetStdDev();
		double max = h1->FindLastBinAbove(0);
		double min = h1->FindFirstBinAbove(0);
		double totalrange = max - min;
		m_adc_mean_y->Fill(i, mean);
		m_adc_min_y->Fill(i, min);
		m_adc_max_y->Fill(i, max);
		m_adc_stddev_y->Fill(i, stddev);
		m_adc_totalrange_y->Fill(i, totalrange);

		double start = mean - stddev;
		if (start < 0) {
			start = 0;
		}
		double end = mean + stddev;
		if (end > 1023) {
			end = 1023;
		}
		m_adc_range_y->Fill(i, end - start);

		delete h1;

	}
	std::cout << "Creating BCID histograms!" << std::endl;
	for (int i = 0; i < m_channels_x; i++) {
		TH1D *h1 = new TH1D("h1", "h1", 4096, 0, 4096);
		TString cut = "hitsX.position==" + std::to_string(i);
		m_tree->Draw("hitsX.bcid>>h1", cut, "goff");
		double mean = h1->GetMean();
		double stddev = h1->GetStdDev();
		double max = h1->FindLastBinAbove(0);
		double min = h1->FindFirstBinAbove(0);
		double totalrange = max - min;
		m_bcid_mean_x->Fill(i, mean);
		m_bcid_min_x->Fill(i, min);
		m_bcid_max_x->Fill(i, max);
		m_bcid_stddev_x->Fill(i, stddev);
		m_bcid_totalrange_x->Fill(i, totalrange);
		double start = mean - stddev;
		if (start < 0) {
			start = 0;
		}
		double end = mean + stddev;
		if (end > 4095) {
			end = 4095;
		}
		m_bcid_range_x->Fill(i, end - start);
		delete h1;
	}

	for (int i = 0; i < m_channels_y; i++) {
		TH1D *h1 = new TH1D("h1", "h1", 4096, 0, 4096);
		TString cut = "hitsY.position==" + std::to_string(i);
		m_tree->Draw("hitsY.bcid>>h1", cut, "goff");
		double mean = h1->GetMean();
		double stddev = h1->GetStdDev();
		double max = h1->FindLastBinAbove(0);
		double min = h1->FindFirstBinAbove(0);
		double totalrange = max - min;
		m_bcid_mean_y->Fill(i, mean);
		m_bcid_min_y->Fill(i, min);
		m_bcid_max_y->Fill(i, max);
		m_bcid_stddev_y->Fill(i, stddev);
		m_bcid_totalrange_y->Fill(i, totalrange);

		double start = mean - stddev;
		if (start < 0) {
			start = 0;
		}
		double end = mean + stddev;
		if (end > 4095) {
			end = 4095;
		}
		m_bcid_range_y->Fill(i, end - start);

		delete h1;
	}

	m_bcid_mean_x->Write("", TObject::kOverwrite);
	m_bcid_min_x->Write("", TObject::kOverwrite);
	m_bcid_max_x->Write("", TObject::kOverwrite);
	m_bcid_stddev_x->Write("", TObject::kOverwrite);
	m_bcid_totalrange_x->Write("", TObject::kOverwrite);
	m_bcid_range_x->Write("", TObject::kOverwrite);
	m_bcid_mean_y->Write("", TObject::kOverwrite);
	m_bcid_min_y->Write("", TObject::kOverwrite);
	m_bcid_max_y->Write("", TObject::kOverwrite);
	m_bcid_stddev_y->Write("", TObject::kOverwrite);
	m_bcid_totalrange_y->Write("", TObject::kOverwrite);
	m_bcid_range_y->Write("", TObject::kOverwrite);

}
