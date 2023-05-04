#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include <vector>
#include <map>
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <TObject.h>
#include <TString.h>
#include "TBufferJSON.h"

#include <chrono>

TString fileInName = "";
TString fileOutName = "";
TString fileNormName = "out";
std::string pAction = "";
std::string pPlot = "";

bool fileInFound = false;
bool fileOutFound = false;
bool fileNormFound = false;

bool plotLog = false;
bool plotPdf = false;
bool plotJSON = false;
int pMap = 56;
int pRebinFactor = 1;
double pZmin = 0;
double pZmax = 0;
double pTimeSignal = 0;
double pTimeBackground = 0;
double pXbins = 640;
double pXmin = 0;
double pXmax=640;
double pYbins = 640;
double pYmin = 0;
double pYmax=640;
std::string pXvar = "pos0";
std::string pYvar = "pos1";
TString cutIn = "";
TString cutNorm = "";

int printUsage(std::string errorMessage, char* argv);

int main(int argc, char**argv) {
	std::chrono::time_point<std::chrono::system_clock> timeEnd, timeStart;

	if (argc == 1 || argc % 2 == 0) {
		return printUsage("Wrong number of arguments!", argv[argc - 1]);
	}
	for (int i = 1; i < argc; i += 2) {
		if (strncmp(argv[i], "-fin", 3) == 0) {
			fileInFound = true;
			fileInName = argv[i + 1];
		} else if (strncmp(argv[i], "-fnorm", 6) == 0) {
			fileNormFound = true;
			fileNormName = argv[i + 1];
		} else if (strncmp(argv[i], "-fout", 5) == 0) {
			fileOutFound = true;
			fileOutName = argv[i + 1];
		} else if (strncmp(argv[i], "-ns", 3) == 0) {
			pTimeSignal = atof(argv[i + 1]);
		} else if (strncmp(argv[i], "-nb", 3) == 0) {
			pTimeBackground = atof(argv[i + 1]);
		} else if (strncmp(argv[i], "-act", 4) == 0) {
			pAction = argv[i + 1];
		} else if (strncmp(argv[i], "-plot", 5) == 0) {
			pPlot = argv[i + 1];
			if (pPlot.find("j") != std::string::npos) {
				plotJSON = true;
			}
			if (pPlot.find("l") != std::string::npos) {
				plotLog = true;
			}
			if (pPlot.find("p") != std::string::npos) {
				plotPdf = true;
			}

		} else if (strncmp(argv[i], "-map", 4) == 0) {
			pMap = atoi(argv[i + 1]);

		}
		else if (strncmp(argv[i], "-xbin", 5) == 0) {
			pXbins = atof(argv[i + 1]);
		} 
		else if (strncmp(argv[i], "-xmin", 5) == 0) {
			pXmin = atof(argv[i + 1]);
		}
		else if (strncmp(argv[i], "-xmax", 5) == 0) {
			pXmax = atof(argv[i + 1]);
		}
		else if (strncmp(argv[i], "-ybin", 5) == 0) {
			pYbins = atof(argv[i + 1]);
		} 
		else if (strncmp(argv[i], "-ymin", 5) == 0) {
			pYmin = atof(argv[i + 1]);
		}
		else if (strncmp(argv[i], "-ymax", 5) == 0) {
			pYmax = atof(argv[i + 1]);
		}
		else if (strncmp(argv[i], "-zmin", 5) == 0) {
			pZmin = atof(argv[i + 1]);
		}
		else if (strncmp(argv[i], "-zmax", 5) == 0) {
			pZmax = atof(argv[i + 1]);
		}
		else if (strncmp(argv[i], "-xvar", 5) == 0) {
			pXvar = argv[i + 1];
		}
		else if (strncmp(argv[i], "-yvar", 5) == 0) {
			pYvar = argv[i + 1];
		}
		else if (strncmp(argv[i], "-cuti", 5) == 0) {
			cutIn = argv[i + 1];
		}
		else if (strncmp(argv[i], "-cutn", 5) == 0) {
			cutNorm = argv[i + 1];
		}
		else {
			return printUsage("Wrong type of argument!", argv[i]);
		}
	}

	if (!fileInFound ||  !fileInName.EndsWith(".root")) {
		return printUsage("Data file has to be loaded with -fin in.root!", nullptr);

	}
	if (!fileNormFound ||  !fileNormName.EndsWith(".root")) {
		return printUsage("Data file has to be loaded with -fnorm norm.root!", nullptr);

	}
	if (fileOutName.EndsWith(".root")) {
		fileOutName.ReplaceAll(".root", "");
	}	

	
	timeStart = std::chrono::system_clock::now();
	TFile *fIn = new TFile(fileInName, "read");
	TH2D *hIn = new TH2D("hIn", "hIn",pXbins,pXmin,pXmax,pYbins, pYmin,pYmax);
	TTree *tIn = (TTree*)fIn->Get("clusters_detector");
	TString drawChoice = pYvar + ":" + pXvar + ">>hIn";

	tIn->Draw(drawChoice, cutIn);
	
	
	TFile *fNorm = new TFile(fileNormName, "read");
	TH2D *hNorm = new TH2D("hNorm", "hNorm",pXbins,pXmin,pXmax,pYbins, pYmin,pYmax);
	TTree *tNorm = (TTree*)fNorm->Get("clusters_detector");
	drawChoice = pYvar + ":" + pXvar + ">>hNorm";
	tNorm->Draw(drawChoice, cutNorm);	
	
	
	TFile *fOut = new TFile(fileOutName+".root", "recreate");
	std::cout << fileOutName+".root" << " opened .." << std::endl;
	
	TH2D *hOut = new TH2D("hOut", "hOut",pXbins,pXmin,pXmax,pYbins, pYmin,pYmax);	
	for (Int_t x = 0; x < pXbins; x++) {
		for (Int_t y = 0; y < pYbins; y++) {
			if (pTimeSignal * pTimeBackground > 0) {
				double valSignal =  pTimeBackground * hIn->GetBinContent(x + 1, y + 1);
				double valBackground = pTimeSignal * hNorm->GetBinContent(x + 1, y + 1);
				if (pAction == "n") {
					double val = 0.0;
					if (valBackground > 0) {
						val = valSignal / valBackground;
					}
					hOut->SetBinContent(x + 1, y + 1, val);
	
				} else {
					double val = 0.0;
					if (valSignal - valBackground > 0) {
						val = valSignal - valBackground;
					}
					hOut->SetBinContent(x + 1, y + 1, val);
				}
			}
			else {
				double val =  hIn->GetBinContent(x + 1, y + 1);
				hOut->SetBinContent(x + 1, y + 1, val);	
			}
		}
	}
	gStyle->SetPalette(pMap);
	gStyle->SetPadRightMargin(0.13);

	TCanvas *c1 = new TCanvas("c1", "c1", 630, 600);

	c1->SetFillColor(0);
	c1->SetBorderMode(0);
	c1->SetBorderSize(2);
	c1->SetFrameBorderMode(0);
	hOut->GetXaxis()->SetTitle("x coordinate [mm]");
	hOut->GetYaxis()->SetTitle("y coordinate [mm]");
	gPad->SetRightMargin(0.18);
	hOut->GetXaxis()->SetTitleOffset(1.4);	
	hOut->GetYaxis()->SetTitleOffset(1.4);
	hOut->GetZaxis()->SetTitleOffset(1.5);						
	hOut->GetZaxis()->SetTitle("counts");
	hOut->GetXaxis()->CenterTitle();
	hOut->GetYaxis()->CenterTitle();
	hOut->GetZaxis()->CenterTitle();
	hOut->SetTitle("");
	hOut->SetStats(0);
	if(pZmax > 0) {
		hOut->GetZaxis()->SetRangeUser(pZmin, pZmax);
	}
	hOut->GetXaxis()->SetRangeUser(pXmin,pXmax);
	hOut->GetYaxis()->SetRangeUser(pYmin,pYmax);
	hOut->Draw("COLZ");
	
				
	if (plotJSON) {
		TString jsonOut = TBufferJSON::ToJSON(hOut, 3);
		std::ofstream fJSON;
		fJSON.open(fileOutName+".json", std::ios::out);
		fJSON << jsonOut;
		fJSON.close();
		std::cout << "JSON" << std::endl;
	}

	if (plotPdf) {
		TString pdfName = "";
		if (plotLog) {
			c1->SetLogz(true);
			pdfName = fileOutName+"_log.pdf";
		} else {
			pdfName = fileOutName+".pdf";
		}
		c1->SaveAs(pdfName);
	}
	fOut->cd();
	c1->Write();
	hOut->Write();
	fOut->Close();
	fIn->Close();
	fNorm->Close();
	
	timeEnd = std::chrono::system_clock::now();

	int elapsed_seconds = std::chrono::duration_cast < std::chrono::milliseconds > (timeEnd - timeStart).count();

	std::cout << "finished computation in " << elapsed_seconds << " ms\n";

	return 0;

}

int printUsage(std::string errorMessage, char* argv) {
	if(argv != nullptr) {
		std::cout << "\nERROR: " << errorMessage << ": " << argv << std::endl;
	}
	else {
		std::cout << "\nERROR: " << errorMessage << std::endl;	
	}

	printf("\nUsages:\n");
	printf("normalize or substract background:\n\t./normalize -fin Run1.root -fnorm Run2.root -fout Out"
			" -ns 0 -nb 0 -act normalize\n");

	return -1;
}

