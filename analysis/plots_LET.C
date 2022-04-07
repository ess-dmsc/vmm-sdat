#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"

#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TMath.h"
#include "TLegend.h"
#include <sstream>
#include <iostream>
#include <TStyle.h>

struct Hit {
  uint32_t id;
  uint8_t det;
  uint8_t plane;
  uint8_t fec;
  uint8_t vmm;
  double readout_time;
  double time;
  uint8_t ch;
  uint16_t pos;
  uint16_t bcid;
  uint16_t tdc;
  uint16_t adc;
  bool over_threshold;
  double chip_time;
};

struct ClusterDetector {
  uint32_t id;
  uint32_t id0;
  uint32_t id1;
  uint8_t det;
  uint16_t size0;
  uint16_t size1;
  uint16_t adc0;
  uint16_t adc1;
  double pos0;
  double pos1;
  double pos2;
  double time0;
  double time1;
  double pos0_utpc;
  double pos1_utpc;
  double pos2_utpc;
  double time0_utpc;
  double time1_utpc;
  double pos0_charge2;
  double pos1_charge2;
  double pos2_charge2;
  double time0_charge2;
  double time1_charge2;
  double pos0_algo;
  double pos1_algo;
  double pos2_algo;
  double time0_algo;
  double time1_algo;
  double dt0;
  double dt1;
  double delta_plane;
  uint16_t span_cluster0;
  uint16_t span_cluster1;
  uint16_t max_delta_time0;
  uint16_t max_delta_time1;
  uint16_t max_missing_strip0;
  uint16_t max_missing_strip1;
  std::vector<double> strips0;
  std::vector<double> times0;
  std::vector<double> strips1;
  std::vector<double> times1;
  std::vector<double> adcs0;
  std::vector<double> adcs1;
};

const int DET_START = 0;
const int DET_END = 2;
const int GRIDS = 51;
const int WIRES = 96;
const int TIME_MEASUREMENT=7;
const double TIME_BIN_SECONDS=1;


TH1D * h_hits_adc_grids[DET_END-DET_START];
TH1D * h_hits_adc_wires[DET_END-DET_START];
TH1D * h_hits_pos_grids[DET_END-DET_START];
TH1D * h_hits_pos_wires[DET_END-DET_START];
TH1D * h_hits_time_grids[DET_END-DET_START];
TH1D * h_hits_time_wires[DET_END-DET_START];

TH2D * h_hits_pos_adc_wires[DET_END-DET_START];
TH2D * h_hits_pos_adc_grids[DET_END-DET_START];

TH1D * h_clusters_time[DET_END-DET_START];
TH1D * h_clusters_pos_wires[DET_END-DET_START];
TH1D * h_clusters_pos_grids[DET_END-DET_START];
TH1D * h_clusters_adc_wires[DET_END-DET_START];
TH1D * h_clusters_adc_grids[DET_END-DET_START];
TH1D * h_clusters_size_wires[DET_END-DET_START];
TH1D * h_clusters_size_grids[DET_END-DET_START];
TH1D * h_clusters_delta_time[DET_END-DET_START];

TH2D * h_clusters_pos_adc_wires[DET_END-DET_START];
TH2D * h_clusters_pos_adc_grids[DET_END-DET_START];
	
TH2D * h_clusters_adc[DET_END-DET_START];
TH2D * h_clusters_pos[DET_END-DET_START];
TH2D * h_clusters_charge_wires[DET_END-DET_START];
TH2D * h_clusters_charge_grids[DET_END-DET_START];
TH2D * h_clusters_charge[DET_END-DET_START];
TH2D * h_clusters_average_size_wires[DET_END-DET_START];
TH2D * h_clusters_average_size_grids[DET_END-DET_START];
TH2D * h_clusters_average_size[DET_END-DET_START];

TMultiGraph *mg; 
Hit hit;
ClusterDetector cluster_detector;
std::string extension;
TTree * th = 0;
TTree * tc = 0;
TCanvas * c = 0;
				
void fillChargeHistogram(bool isLog);				
void setup_histograms();
void delete_histograms();


int plots_LET(TString nameFilter, TString path=".", bool isLog=false) {
	gStyle->SetPadRightMargin(0.12);
	gStyle->SetPadLeftMargin(0.11);
	//gStyle->SetPadBottomMargin(0.10);
	//gStyle->SetPadTopMargin(0.03);
	gStyle->SetErrorX(0);
	
	TSystemDirectory dir(path, path);
	TList *files = dir.GetListOfFiles();
	if (files) {
		TSystemFile *file;
		TString fname;
		TIter next(files);

		while ((file = (TSystemFile*) next())) {
			fname = file->GetName();
	
			if (!file->IsDirectory()
					&& fname.BeginsWith(nameFilter)
					&& fname.EndsWith(".root")) {

				std::cout << "Analysing " << fname << std::endl;
				
				
				TFile* f = new TFile(path + "/" + fname);
				setup_histograms();
				

				th = (TTree*) f->Get("hits");
			    tc = (TTree*) f->Get("clusters_detector");

				tc->SetBranchAddress("det", &cluster_detector.det);
				tc->SetBranchAddress("size0", &cluster_detector.size0);
				tc->SetBranchAddress("size1", &cluster_detector.size1);
				tc->SetBranchAddress("adc0", &cluster_detector.adc0);
				tc->SetBranchAddress("adc1", &cluster_detector.adc1);
				tc->SetBranchAddress("pos0", &cluster_detector.pos0);
				tc->SetBranchAddress("pos1", &cluster_detector.pos1);
				tc->SetBranchAddress("time0", &cluster_detector.time0);
				tc->SetBranchAddress("time1", &cluster_detector.time1);
			
				c = new TCanvas("c", "c", 1000, 750);
				
				
				if(isLog) {
					c->SetLogy(true);
					extension = "_log.png";
				}
				else {
					c->SetLogy(false);
					extension = ".png";
				}
				
				
				for(int n=DET_START; n<DET_END; n++) {
					
					gStyle->SetOptStat(1111);
					//1D hits
					std::string variable = "adc>>hits_adc_grids_" + std::to_string(n);
					std::string filter = "plane==0&&det=="+ std::to_string(n);
					std::string name = "hits_adc_grids_" + std::to_string(n) + extension;
					th->Draw(variable.c_str(), filter.c_str());
					c->SaveAs(name.c_str());
				
					variable = "adc>>hits_adc_wires_" + std::to_string(n);
					filter = "plane==1&det=="+ std::to_string(n);
					name = "hits_adc_wires_" + std::to_string(n) + extension;
					th->Draw(variable.c_str(), filter.c_str());
					c->SaveAs(name.c_str());
					
					gStyle->SetOptStat(11);
					
					variable = "pos>>hits_pos_grids_" + std::to_string(n);
					filter = "plane==0&&det=="+ std::to_string(n);
					name = "hits_pos_grids_" + std::to_string(n) + extension;
					th->Draw(variable.c_str(), filter.c_str());
					c->SaveAs(name.c_str());
					
					variable = "pos>>hits_pos_wires_" + std::to_string(n);
					filter = "plane==1&det=="+ std::to_string(n);
					name = "hits_pos_wires_" + std::to_string(n) + extension;
					th->Draw(variable.c_str(), filter.c_str());
					c->SaveAs(name.c_str());
					
					variable = "time*1e-09>>hits_time_grids_" + std::to_string(n);
					filter = "plane==0&&det=="+ std::to_string(n);	
					name = "hits_time_grids_" + std::to_string(n) + extension;
					th->Draw(variable.c_str(), filter.c_str());
					h_hits_time_grids[n]->Scale(1.0/TIME_BIN_SECONDS);
					h_hits_time_grids[n]->SetMarkerStyle(kFullTriangleDown);
					h_hits_time_grids[n]->SetMarkerColor(kRed);
					c->SaveAs(name.c_str());
					
					variable = "time*1e-09>>hits_time_wires_" + std::to_string(n);
					filter = "plane==1&&det=="+ std::to_string(n);
					name = "hits_time_wires_" + std::to_string(n) + extension;
					th->Draw(variable.c_str(), filter.c_str());
					h_hits_time_wires[n]->Scale(1.0/TIME_BIN_SECONDS);
					h_hits_time_wires[n]->SetMarkerStyle(kFullSquare);
					h_hits_time_wires[n]->SetMarkerColor(kRed);
					c->SaveAs(name.c_str());	
					
					
					gStyle->SetOptStat(1111);
			
					//1D clusters
					variable = "adc0>>clusters_adc_grids_" + std::to_string(n);
					filter = "det=="+ std::to_string(n);
					name = "clusters_adc_grids_" + std::to_string(n) + extension;
					tc->Draw(variable.c_str(), filter.c_str());
					c->SaveAs(name.c_str());
					
					variable = "adc1>>clusters_adc_wires_" + std::to_string(n);
					filter = "det=="+ std::to_string(n);
					name = "clusters_adc_wires_" + std::to_string(n) + extension;
					tc->Draw(variable.c_str(), filter.c_str());
					c->SaveAs(name.c_str());
					
					gStyle->SetOptStat(11);
					
					variable = "pos0>>clusters_pos_grids_" + std::to_string(n);
					filter = "det=="+ std::to_string(n);
					name = "clusters_pos_grids_" + std::to_string(n) + extension;
					tc->Draw(variable.c_str(), filter.c_str());
					c->SaveAs(name.c_str());
					
					variable = "pos1>>clusters_pos_wires_" + std::to_string(n);
					filter = "det=="+ std::to_string(n);
					name = "clusters_pos_wires_" + std::to_string(n) + extension;
					tc->Draw(variable.c_str(), filter.c_str());
					c->SaveAs(name.c_str());
					
					gStyle->SetOptStat(1111);
					
					variable = "size0>>clusters_size_grids_" + std::to_string(n);
					filter = "det=="+ std::to_string(n);
					name = "clusters_size_grids_" + std::to_string(n) + extension;
					tc->Draw(variable.c_str(), filter.c_str());
					c->SaveAs(name.c_str());
					
					variable = "size1>>clusters_size_wires_" + std::to_string(n);
					filter = "det=="+ std::to_string(n);
					name = "clusters_size_wires_" + std::to_string(n) + extension;
					tc->Draw(variable.c_str(), filter.c_str());
					c->SaveAs(name.c_str());
					
					variable = "time1-time0>>clusters_delta_time_" + std::to_string(n);
					filter = "det=="+ std::to_string(n);
					name = "clusters_delta_time_" + std::to_string(n) + extension;
					tc->Draw(variable.c_str(), filter.c_str());
					c->SaveAs(name.c_str());
					
					gStyle->SetOptStat(11);
					
					variable = "adc0:pos0>>clusters_pos_adc_grids_" + std::to_string(n);
					filter = "det=="+ std::to_string(n);
					name = "clusters_pos_adc_grids_" + std::to_string(n) + extension;
					tc->Draw(variable.c_str(), filter.c_str());
					h_clusters_pos_adc_grids[n]->Draw("COLZ");
					gPad->SetLogz(isLog);
					h_clusters_pos_adc_grids[n]->GetZaxis()->SetTitle("counts");
					h_clusters_pos_adc_grids[n]->GetZaxis()->CenterTitle(true);
  					c->SaveAs(name.c_str());
			
					variable = "adc1:pos1>>clusters_pos_adc_wires_" + std::to_string(n);
					filter = "det=="+ std::to_string(n);
					name = "clusters_pos_adc_wires_" + std::to_string(n) + extension;
					tc->Draw(variable.c_str(), filter.c_str());
					h_clusters_pos_adc_wires[n]->Draw("COLZ");
					gPad->SetLogz(isLog);
					h_clusters_pos_adc_wires[n]->GetZaxis()->SetTitle("counts");
					h_clusters_pos_adc_wires[n]->GetZaxis()->CenterTitle(true);
					c->SaveAs(name.c_str());
					
					variable = "adc1:adc0>>clusters_adc_" + std::to_string(n);
					filter = "det=="+ std::to_string(n);
					name = "clusters_adc_" + std::to_string(n) + extension;
					tc->Draw(variable.c_str(), filter.c_str());
					h_clusters_adc[n]->Draw("COLZ");
					gPad->SetLogz(isLog);
					h_clusters_adc[n]->GetZaxis()->SetTitle("counts");
					h_clusters_adc[n]->GetZaxis()->CenterTitle(true);
					c->SaveAs(name.c_str());
					
					variable = "pos0:pos1>>clusters_pos_" + std::to_string(n);
					filter = "det=="+ std::to_string(n);
					name = "clusters_pos_" + std::to_string(n) + extension;
					tc->Draw(variable.c_str(), filter.c_str());
					h_clusters_pos[n]->Draw("COLZ");
					gPad->SetLogz(isLog);
					h_clusters_pos[n]->GetZaxis()->SetTitle("counts");
					h_clusters_pos[n]->GetZaxis()->CenterTitle(true);
					c->SaveAs(name.c_str());
	
				}
				fillChargeHistogram(isLog);
				c->Close();
				f->Close();
							
			}
		}
	}
	return 0;
}


void fillChargeHistogram(bool isLog) {
	gStyle->SetOptStat(11);
	for(int n=DET_START; n<DET_END;n++) {	
		std::string name = "clusters_charge_wires_" + std::to_string(n);
		h_clusters_charge_wires[n] = new TH2D(name.c_str(), name.c_str(), WIRES, 0, WIRES, GRIDS, 0, GRIDS);
		h_clusters_charge_wires[n]->GetYaxis()->SetTitle("position [grid]");
		h_clusters_charge_wires[n]->GetXaxis()->CenterTitle(true);
		h_clusters_charge_wires[n]->GetXaxis()->SetTitle("position [wire]");
		h_clusters_charge_wires[n]->GetYaxis()->CenterTitle(true);
		std::string title = "Column " + std::to_string(n) + " - cluster charge wire vs. position";
		h_clusters_charge_wires[n]->SetTitle(title.c_str());
	
		name = "clusters_charge_grids_" + std::to_string(n);
		h_clusters_charge_grids[n] = new TH2D(name.c_str(), name.c_str(),WIRES, 0, WIRES, GRIDS, 0, GRIDS);
		h_clusters_charge_grids[n]->GetYaxis()->SetTitle("position [grid]");
		h_clusters_charge_grids[n]->GetXaxis()->CenterTitle(true);
		h_clusters_charge_grids[n]->GetXaxis()->SetTitle("position [wire]");
		h_clusters_charge_grids[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - cluster charge grid vs. position";
		h_clusters_charge_grids[n]->SetTitle(title.c_str());
	
		name = "clusters_charge_" + std::to_string(n);
		h_clusters_charge[n] = new TH2D(name.c_str(), name.c_str(),WIRES, 0, WIRES, GRIDS, 0, GRIDS);
		h_clusters_charge[n]->GetYaxis()->SetTitle("position [grid]");
		h_clusters_charge[n]->GetXaxis()->CenterTitle(true);
		h_clusters_charge[n]->GetXaxis()->SetTitle("position [wire]");
		h_clusters_charge[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - cluster charge wire plus grid vs. position";
		h_clusters_charge_wires[n]->SetTitle(title.c_str());
	
		name = "clusters_average_size_wires_" + std::to_string(n);
		h_clusters_average_size_wires[n] = new TH2D(name.c_str(), name.c_str(), WIRES, 0, WIRES, GRIDS, 0, GRIDS);
		h_clusters_average_size_wires[n]->GetYaxis()->SetTitle("position [grid]");
		h_clusters_average_size_wires[n]->GetXaxis()->CenterTitle(true);
		h_clusters_average_size_wires[n]->GetXaxis()->SetTitle("position [wire]");
		h_clusters_average_size_wires[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - cluster size wire vs. position";
		h_clusters_average_size_wires[n]->SetTitle(title.c_str());
	
		name = "clusters_average_size_grids_" + std::to_string(n);
		h_clusters_average_size_grids[n] = new TH2D(name.c_str(), name.c_str(),WIRES, 0, WIRES, GRIDS, 0, GRIDS);
		h_clusters_average_size_grids[n]->GetYaxis()->SetTitle("position [grid]");
		h_clusters_average_size_grids[n]->GetXaxis()->CenterTitle(true);
		h_clusters_average_size_grids[n]->GetXaxis()->SetTitle("position [wire]");
		h_clusters_average_size_grids[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - cluster size grid vs. position";
		h_clusters_average_size_grids[n]->SetTitle(title.c_str());
	
		name = "clusters_average_size_" + std::to_string(n);
		h_clusters_average_size[n] = new TH2D(name.c_str(), name.c_str(),WIRES, 0, WIRES, GRIDS, 0, GRIDS);
		h_clusters_average_size[n]->GetYaxis()->SetTitle("position [grid]");
		h_clusters_average_size[n]->GetXaxis()->CenterTitle(true);
		h_clusters_average_size[n]->GetXaxis()->SetTitle("position [wire]");
		h_clusters_average_size[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - cluster size wire plus grid vs. position";
		h_clusters_average_size[n]->SetTitle(title.c_str());
		
		name = "clusters_time_" + std::to_string(n);
		h_clusters_time[n] = new TH1D(name.c_str(), name.c_str(), TIME_MEASUREMENT/TIME_BIN_SECONDS, 0, TIME_MEASUREMENT);
		h_clusters_time[n]->GetYaxis()->SetTitle("rate [Hz]");
		h_clusters_time[n]->GetXaxis()->CenterTitle(true);
		h_clusters_time[n]->GetXaxis()->SetTitle("measurement time [s]");
		h_clusters_time[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - cluster rate";
		h_clusters_time[n]->SetTitle(title.c_str());		
	}
		
	int numEvents = tc->GetEntries();

	for (int i = 0; i < numEvents; i++) {
		tc->GetEntry(i);
		if(i%1000000==0) {
			std::cout << "Entry " << i << std::endl;
		}
		if(cluster_detector.det<DET_END && cluster_detector.det>= DET_START) {
			h_clusters_charge_grids[cluster_detector.det]->Fill(cluster_detector.pos1, cluster_detector.pos0, cluster_detector.adc0);  
			h_clusters_charge_wires[cluster_detector.det]->Fill(cluster_detector.pos1, cluster_detector.pos0, cluster_detector.adc1); 
			h_clusters_charge[cluster_detector.det]->Fill(cluster_detector.pos1, cluster_detector.pos0, cluster_detector.adc0+cluster_detector.adc1);
			h_clusters_average_size_grids[cluster_detector.det]->Fill(cluster_detector.pos1, cluster_detector.pos0, cluster_detector.size0);
			h_clusters_average_size_wires[cluster_detector.det]->Fill(cluster_detector.pos1, cluster_detector.pos0, cluster_detector.size1);
			h_clusters_average_size[cluster_detector.det]->Fill(cluster_detector.pos1, cluster_detector.pos0, cluster_detector.size0+cluster_detector.size1);
			h_clusters_time[cluster_detector.det]->Fill(cluster_detector.time0*1e-09,1.0/TIME_BIN_SECONDS);
			
		} 
		    
    }
    
    for(int n=DET_START; n<DET_END;n++) {
    	int events=0;
    	for (int b0 = 1; b0 <= GRIDS; b0++) {
        	for (int b1 = 1; b1 <= WIRES; b1++) {
          		int cnt = h_clusters_pos[n]->GetBinContent(b1, b0);
          		double val = 0;
          		if (cnt > 0) {
            		events++;
            		val = h_clusters_average_size_grids[n]->GetBinContent(b1, b0) / cnt;
            		h_clusters_average_size_grids[n]->SetBinContent(b1, b0, val);
            		
            		val = h_clusters_average_size_wires[n]->GetBinContent(b1, b0) / cnt;
            		h_clusters_average_size_wires[n]->SetBinContent(b1, b0, val);

            		val = h_clusters_average_size[n]->GetBinContent(b1, b0) / cnt;
            		h_clusters_average_size[n]->SetBinContent(b1, b0, val);
            		
            		val = h_clusters_charge_grids[n]->GetBinContent(b1, b0) / cnt;
            		h_clusters_charge_grids[n]->SetBinContent(b1, b0, val);
            		
            		val = h_clusters_charge_wires[n]->GetBinContent(b1, b0) / cnt;
            		h_clusters_charge_wires[n]->SetBinContent(b1, b0, val);
            		
            		val = h_clusters_charge[n]->GetBinContent(b1, b0) / cnt;
            		h_clusters_charge[n]->SetBinContent(b1, b0, val);
          		}
        	}
      	}
      	h_clusters_average_size_grids[n]->SetEntries(events);
		h_clusters_average_size_wires[n]->SetEntries(events);
		h_clusters_average_size[n]->SetEntries(events);
		h_clusters_charge_grids[n]->SetEntries(events);
		h_clusters_charge_wires[n]->SetEntries(events);
		h_clusters_charge[n]->SetEntries(events);
		
		gPad->SetLogz(isLog);
		
		std::string name = "clusters_average_size_grids_" + std::to_string(n) + extension;
		h_clusters_average_size_grids[n]->Draw("COLZ");
		h_clusters_average_size_grids[n]->GetZaxis()->SetTitle("mean size");
		h_clusters_average_size_grids[n]->GetZaxis()->CenterTitle(true);
		c->SaveAs(name.c_str());

		name = "clusters_average_size_wires_" + std::to_string(n) + extension;
		h_clusters_average_size_wires[n]->Draw("COLZ");
		h_clusters_average_size_wires[n]->GetZaxis()->SetTitle("mean size");
		h_clusters_average_size_wires[n]->GetZaxis()->CenterTitle(true);
		c->SaveAs(name.c_str());
		
		name = "clusters_average_size_" + std::to_string(n) + extension;
		h_clusters_average_size[n]->Draw("COLZ");
		h_clusters_average_size[n]->GetZaxis()->SetTitle("mean size");
		h_clusters_average_size[n]->GetZaxis()->CenterTitle(true);
		c->SaveAs(name.c_str());
		
		name = "clusters_charge_grids_" + std::to_string(n) + extension;
		h_clusters_charge_grids[n]->Draw("COLZ");
		h_clusters_charge_grids[n]->GetZaxis()->SetTitle("mean charge [ADC]");
		h_clusters_charge_grids[n]->GetZaxis()->CenterTitle(true);
		c->SaveAs(name.c_str());
		
		name = "clusters_charge_wires_" + std::to_string(n) + extension;
		h_clusters_charge_wires[n]->Draw("COLZ");
		h_clusters_charge_wires[n]->GetZaxis()->SetTitle("mean charge [ADC]");
		h_clusters_charge_wires[n]->GetZaxis()->CenterTitle(true);
		c->SaveAs(name.c_str());
		
		name = "clusters_charge_" + std::to_string(n) + extension;
		h_clusters_charge[n]->Draw("COLZ");
		h_clusters_charge[n]->GetZaxis()->SetTitle("mean charge [ADC]");
		h_clusters_charge[n]->GetZaxis()->CenterTitle(true);
		c->SaveAs(name.c_str());			
		
		name = "clusters_time_" + std::to_string(n) + extension;
		h_clusters_time[n]->Draw("hist p");
		h_clusters_time[n]->SetMarkerStyle(kFullSquare);
		h_clusters_time[n]->SetMarkerColor(kRed);
		c->SaveAs(name.c_str());
		
		
    }
   
}



void setup_histograms() {
	for(int n=DET_START; n<DET_END; n++) {
		//1D hits
		std::string name = "hits_adc_grids_" + std::to_string(n);
		h_hits_adc_grids[n] = new TH1D(name.c_str(), name.c_str(), 1024, 0, 1024);
		h_hits_adc_grids[n]->GetYaxis()->SetTitle("counts");
		h_hits_adc_grids[n]->GetXaxis()->CenterTitle(true);
		h_hits_adc_grids[n]->GetXaxis()->SetTitle("ADC");
		h_hits_adc_grids[n]->GetYaxis()->CenterTitle(true);
		std::string title = "Grids - column " + std::to_string(n) + " - ADC";
		h_hits_adc_grids[n]->SetTitle(title.c_str());

		name = "hits_adc_wires_" + std::to_string(n);
		h_hits_adc_wires[n] = new TH1D(name.c_str(), name.c_str(), 1024, 0, 1024);
		h_hits_adc_wires[n]->GetYaxis()->SetTitle("counts");
		h_hits_adc_wires[n]->GetXaxis()->CenterTitle(true);
		h_hits_adc_wires[n]->GetXaxis()->SetTitle("ADC");
		h_hits_adc_wires[n]->GetYaxis()->CenterTitle(true);
		title = "Wires - column " + std::to_string(n) + " - ADC";
		h_hits_adc_wires[n]->SetTitle(title.c_str());
		
		name = "hits_pos_grids_" + std::to_string(n);
		h_hits_pos_grids[n] = new TH1D(name.c_str(), name.c_str(), GRIDS, 0, GRIDS);
		h_hits_pos_grids[n]->GetYaxis()->SetTitle("counts");
		h_hits_pos_grids[n]->GetXaxis()->CenterTitle(true);
		h_hits_pos_grids[n]->GetXaxis()->SetTitle("position [grid]");
		h_hits_pos_grids[n]->GetYaxis()->CenterTitle(true);
		title = "Grids - column " + std::to_string(n) + " - Position";
		h_hits_pos_grids[n]->SetTitle(title.c_str());
		
		name = "hits_pos_wires_" + std::to_string(n);
		h_hits_pos_wires[n] = new TH1D(name.c_str(), name.c_str(), WIRES, 0, WIRES);
		h_hits_pos_wires[n]->GetYaxis()->SetTitle("counts");
		h_hits_pos_wires[n]->GetXaxis()->CenterTitle(true);
		h_hits_pos_wires[n]->GetXaxis()->SetTitle("position [wire]");
		h_hits_pos_wires[n]->GetYaxis()->CenterTitle(true);
		title = "Wires - column " + std::to_string(n) + " - Position";
		h_hits_pos_wires[n]->SetTitle(title.c_str());	
		
				
		name = "hits_time_grids_" + std::to_string(n);
		h_hits_time_grids[n] = new TH1D(name.c_str(), name.c_str(), TIME_MEASUREMENT/TIME_BIN_SECONDS, 0, TIME_MEASUREMENT);
		h_hits_time_grids[n]->GetYaxis()->SetTitle("rate [Hz]");
		h_hits_time_grids[n]->GetXaxis()->CenterTitle(true);
		h_hits_time_grids[n]->GetXaxis()->SetTitle("measurement time [s]");
		h_hits_time_grids[n]->GetYaxis()->CenterTitle(true);
		title = "Grids - column " + std::to_string(n) + " - Time";
		h_hits_time_grids[n]->SetTitle(title.c_str());
		
		name = "hits_time_wires_" + std::to_string(n);
		h_hits_time_wires[n] = new TH1D(name.c_str(), name.c_str(), TIME_MEASUREMENT/TIME_BIN_SECONDS, 0, TIME_MEASUREMENT);
		h_hits_time_wires[n]->GetYaxis()->SetTitle("rate [Hz]");
		h_hits_time_wires[n]->GetXaxis()->CenterTitle(true);
		h_hits_time_wires[n]->GetXaxis()->SetTitle("measurement time [s]");
		h_hits_time_wires[n]->GetYaxis()->CenterTitle(true);
		title = "Wires - column " + std::to_string(n) + " - Time";
		h_hits_time_wires[n]->SetTitle(title.c_str());
		
		//1D clusters
		name = "clusters_adc_grids_" + std::to_string(n);
		h_clusters_adc_grids[n] = new TH1D(name.c_str(), name.c_str(), 3072, 0, 3072);
		h_clusters_adc_grids[n]->GetYaxis()->SetTitle("counts");
		h_clusters_adc_grids[n]->GetXaxis()->CenterTitle(true);
		h_clusters_adc_grids[n]->GetXaxis()->SetTitle("ADC");
		h_clusters_adc_grids[n]->GetYaxis()->CenterTitle(true);
		title = "Grids - column " + std::to_string(n) + " - cluster ADC";
		h_clusters_adc_grids[n]->SetTitle(title.c_str());
		
		
		name = "clusters_adc_wires_" + std::to_string(n);
		h_clusters_adc_wires[n] = new TH1D(name.c_str(), name.c_str(), 1024, 0, 1024);
		h_clusters_adc_wires[n]->GetYaxis()->SetTitle("counts");
		h_clusters_adc_wires[n]->GetXaxis()->CenterTitle(true);
		h_clusters_adc_wires[n]->GetXaxis()->SetTitle("ADC");
		h_clusters_adc_wires[n]->GetYaxis()->CenterTitle(true);
		title = "Wires - column " + std::to_string(n) + " - cluster ADC";
		h_clusters_adc_wires[n]->SetTitle(title.c_str());
		
		name = "clusters_pos_grids_" + std::to_string(n);
		h_clusters_pos_grids[n] = new TH1D(name.c_str(), name.c_str(), GRIDS, 0, GRIDS);
		h_clusters_pos_grids[n]->GetYaxis()->SetTitle("counts");
		h_clusters_pos_grids[n]->GetXaxis()->CenterTitle(true);
		h_clusters_pos_grids[n]->GetXaxis()->SetTitle("position [grid]");
		h_clusters_pos_grids[n]->GetYaxis()->CenterTitle(true);
		title = "Grids - column " + std::to_string(n) + " - cluster position";
		h_clusters_pos_grids[n]->SetTitle(title.c_str());
		
		name = "clusters_pos_wires_" + std::to_string(n);
		h_clusters_pos_wires[n] = new TH1D(name.c_str(), name.c_str(), WIRES, 0, WIRES); 
		h_clusters_pos_wires[n]->GetYaxis()->SetTitle("counts");
		h_clusters_pos_wires[n]->GetXaxis()->CenterTitle(true);
		h_clusters_pos_wires[n]->GetXaxis()->SetTitle("position [wire]");
		h_clusters_pos_wires[n]->GetYaxis()->CenterTitle(true);
		title = "Wires - column " + std::to_string(n) + " - cluster position";
		h_clusters_pos_wires[n]->SetTitle(title.c_str());

		name = "clusters_size_grids_" + std::to_string(n);
		h_clusters_size_grids[n] = new TH1D(name.c_str(), name.c_str(), 10, 0, 10);
		h_clusters_size_grids[n]->GetYaxis()->SetTitle("counts");
		h_clusters_size_grids[n]->GetXaxis()->CenterTitle(true);
		h_clusters_size_grids[n]->GetXaxis()->SetTitle("size");
		h_clusters_size_grids[n]->GetYaxis()->CenterTitle(true);
		title = "Grids - column " + std::to_string(n) + " - cluster size";
		h_clusters_size_grids[n]->SetTitle(title.c_str());
				
		name = "clusters_size_wires_" + std::to_string(n);
		h_clusters_size_wires[n] = new TH1D(name.c_str(), name.c_str(), 10, 0, 10);
		h_clusters_size_wires[n]->GetYaxis()->SetTitle("counts");
		h_clusters_size_wires[n]->GetXaxis()->CenterTitle(true);
		h_clusters_size_wires[n]->GetXaxis()->SetTitle("size");
		h_clusters_size_wires[n]->GetYaxis()->CenterTitle(true);
		title = "Wires - column " + std::to_string(n) + " - cluster size";
		h_clusters_size_wires[n]->SetTitle(title.c_str());
		
		name = "clusters_delta_time_" + std::to_string(n);
		h_clusters_delta_time[n] = new TH1D(name.c_str(), name.c_str(), 1000, -500, 500);
		h_clusters_delta_time[n]->GetYaxis()->SetTitle("counts");
		h_clusters_delta_time[n]->GetXaxis()->CenterTitle(true);
		h_clusters_delta_time[n]->GetXaxis()->SetTitle("delta time wire minus grid [ns]");
		h_clusters_delta_time[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - Delta time wire minus grid";
		h_clusters_delta_time[n]->SetTitle(title.c_str());
		
		//2D hits
		name = "hits_pos_adc_grids_" + std::to_string(n);
		h_hits_pos_adc_grids[n] = new TH2D(name.c_str(), name.c_str(),  GRIDS, 0, GRIDS, 1024,0,1024);
		h_hits_pos_adc_grids[n]->GetYaxis()->SetTitle("ADC");
		h_hits_pos_adc_grids[n]->GetXaxis()->CenterTitle(true);
		h_hits_pos_adc_grids[n]->GetXaxis()->SetTitle("position [grid]");
		h_hits_pos_adc_grids[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - hit ADC vs. grid position";
		h_hits_pos_adc_grids[n]->SetTitle(title.c_str());
		
		name = "hits_pos_adc_wires_" + std::to_string(n);
		h_hits_pos_adc_wires[n] = new TH2D(name.c_str(), name.c_str(),  WIRES, 0, WIRES, 1024,0,1024);
		h_hits_pos_adc_wires[n]->GetYaxis()->SetTitle("ADC");
		h_hits_pos_adc_wires[n]->GetXaxis()->CenterTitle(true);
		h_hits_pos_adc_wires[n]->GetXaxis()->SetTitle("position [wire]");
		h_hits_pos_adc_wires[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - hit ADC vs. wire position";
		h_hits_pos_adc_wires[n]->SetTitle(title.c_str());
				

		//2D clusters
		name = "clusters_pos_adc_grids_" + std::to_string(n);
		h_clusters_pos_adc_grids[n] = new TH2D(name.c_str(), name.c_str(),  GRIDS, 0, GRIDS, 1024,0,3*1024);
		h_clusters_pos_adc_grids[n]->GetYaxis()->SetTitle("ADC");
		h_clusters_pos_adc_grids[n]->GetXaxis()->CenterTitle(true);
		h_clusters_pos_adc_grids[n]->GetXaxis()->SetTitle("position [grid]");
		h_clusters_pos_adc_grids[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - cluster ADC vs. grid position";
		h_clusters_pos_adc_grids[n]->SetTitle(title.c_str());

		name = "clusters_pos_adc_wires_" + std::to_string(n);
		h_clusters_pos_adc_wires[n] = new TH2D(name.c_str(), name.c_str(),  WIRES, 0, WIRES, 1024,0,1024);
		h_clusters_pos_adc_wires[n]->GetYaxis()->SetTitle("ADC");
		h_clusters_pos_adc_wires[n]->GetXaxis()->CenterTitle(true);
		h_clusters_pos_adc_wires[n]->GetXaxis()->SetTitle("position [wire]");
		h_clusters_pos_adc_wires[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - cluster ADC vs. wire position";
		h_clusters_pos_adc_wires[n]->SetTitle(title.c_str());
	
	
		name = "clusters_adc_" + std::to_string(n);
		h_clusters_adc[n] = new TH2D(name.c_str(), name.c_str(), 1024,0,3072,1024,0,1024);
		h_clusters_adc[n]->GetYaxis()->SetTitle("grid ADC");
		h_clusters_adc[n]->GetXaxis()->CenterTitle(true);
		h_clusters_adc[n]->GetXaxis()->SetTitle("wire ADC");
		h_clusters_adc[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - cluster grid ADC vs. wire ADC";
		h_clusters_adc[n]->SetTitle(title.c_str());
	
	
		name = "clusters_pos_" + std::to_string(n);
		h_clusters_pos[n] = new TH2D(name.c_str(), name.c_str(), WIRES, 0, WIRES, GRIDS, 0, GRIDS);
		h_clusters_pos[n]->GetYaxis()->SetTitle("position [grid]");
		h_clusters_pos[n]->GetXaxis()->CenterTitle(true);
		h_clusters_pos[n]->GetXaxis()->SetTitle("position [wire]");
		h_clusters_pos[n]->GetYaxis()->CenterTitle(true);
		title = "Column " + std::to_string(n) + " - cluster grid position vs. wire position";
		h_clusters_pos[n]->SetTitle(title.c_str());
	
	}
}

