// Cluster-trigger-matching
// -------------------------------------------------------------------
// lucian.scharenberg@cern.ch
// 06/Oct/2022 and 04/Oct/2023 and 05/Oct/2023

#include "../../src/DataStructures.h"

using namespace ROOT;
using namespace std;


double myGaus(double *x, double *par) {

    double scale, mean, sigma;

    mean = par[1];
    sigma = par[2];
    scale = par[0];

    return scale * TMath::Exp(- TMath::Power(x[0] - mean, 2) / (2 * TMath::Power(sigma, 2)));

}


pair<vector<pair<double, double>>, vector<pair<double, double>>> matching(vector<pair<double, double>> reference, vector<pair<double, double>> candidate) {

    vector<pair<double, double>> referenceMatched;
    vector<pair<double, double>> candidateMatched;

    int i = 0;
    int k = 0;
    double timeDifference;
    double dtMax = 5000.0;
    
    while (i < reference.size()) {

	if (i % 1000 == 0) {

	    cout << i << endl;

	}

	while (k < candidate.size()) {

	    timeDifference = candidate[k].first - reference[i].first;

	    if (timeDifference <= -dtMax) {

		k += 1;
		continue;

	    }

	    else if (timeDifference > -dtMax && timeDifference < dtMax) {

		referenceMatched.push_back(make_pair(reference[i].first, reference[i].second));
		candidateMatched.push_back(make_pair(candidate[k].first, candidate[k].second));

		k += 1;
		i += 1;
		break;

	    }

	    else if (timeDifference > dtMax) {

		i += 1;
		break;

	    }

	    else if (i == reference.size() || k == candidate.size()) {

		break;

	    }

	}

	if (i == reference.size() || k == candidate.size()) {

	    break;

	}

    }

    return {referenceMatched, candidateMatched};

}


void clusterTriggerMatchingPosition() {

    vector<double> posClusters, posTrigger, timeDifference;
  
    // Open root file
    // Run7 = PosA (central)
    // Run2 = PosB (close to cathode)
    // Run3 = PosC (close to readout)
    TFile *f = new TFile("./run3.root");
    TTree *treeClusters = (TTree*)f->Get("clusters_plane");
    TTree *treeHits = (TTree*)f->Get("hits");

    // Fill the data arrays while applying cuts
    vector<pair<double, double>> events;
    vector<pair<double, double>> clusters;

    char hitsDet;
    double hitsTime;
    unsigned short hitsADC;
    unsigned short hitsPos;
    treeHits->SetBranchAddress("hits.det", &hitsDet);
    treeHits->SetBranchAddress("hits.time", &hitsTime);
    treeHits->SetBranchAddress("hits.adc", &hitsADC);
    treeHits->SetBranchAddress("hits.pos", &hitsPos);

    for (int i = 0; i < treeHits->GetEntries(); i++) {

	treeHits->GetEntry(i);

	if ((int)hitsDet == 1 && hitsPos == 62 && hitsADC > 250.0) {

	    events.push_back(make_pair(hitsTime, hitsPos));

	}

    }

    cout << events.size() << endl;
    
    char clustersDet;
    double clustersTime;
    double clustersPos;
    treeClusters->SetBranchAddress("clusters_plane.det", &clustersDet);
    treeClusters->SetBranchAddress("clusters_plane.time", &clustersTime);
    treeClusters->SetBranchAddress("clusters_plane.pos", &clustersPos);

    for (int i = 0; i < treeClusters->GetEntries(); i++) {

	treeClusters->GetEntry(i);

	if ((int)clustersDet == 0) {

	    clusters.push_back(make_pair(clustersTime, clustersPos));

	}

    }

    // Match the events
    sort(clusters.begin(), clusters.end());
    sort(events.begin(), events.end());

    pair<vector<pair<double, double>>, vector<pair<double, double>>> match = matching(events, clusters);

    for (int i = 0; i < match.first.size(); i++) {
	
	posClusters.push_back(match.second[i].second);
	posTrigger.push_back(match.first[i].second);
	timeDifference.push_back(match.second[i].first - match.first[i].first);
	
    }

    // Plot the data
    gROOT->LoadMacro("./fellowStyle.cpp");
    gROOT->ProcessLine("fellowStyle(1)");
    double w = 1000 * 1.2;
    double h = 1000;
    TCanvas *c = new TCanvas("c", "", w, h);
    c->SetWindowSize(w + (w - c->GetWw()), h + (h - c->GetWh()));
    gROOT->ForceStyle();

    double axisLimitX1 = 0.0;
    double axisLimitX2 = 480.0;
    double axisLimitY1 = -5.0;
    double axisLimitY2 = 5.0;
    double axisLimitZ1 = 0.0;
    double axisLimitZ2 = 20.0;

    TH2D *hist = new TH2D("hist", "", 192, axisLimitX1, axisLimitX2, 200, axisLimitY1, axisLimitY2);

    cout << match.first.size() << endl;
    
    for (int i = 0; i < match.first.size(); i++) {

	hist->Fill(posClusters[i], timeDifference[i] / 1000.0);
	
    }

    hist->GetXaxis()->SetTitle("Position / strips [1.2 mm]");
    hist->GetXaxis()->SetLimits(axisLimitX1, axisLimitX2);
    hist->GetXaxis()->CenterTitle(true);
    hist->GetYaxis()->SetTitle("Time difference / us");
    hist->GetYaxis()->SetRangeUser(axisLimitY1, axisLimitY2);
    hist->GetYaxis()->CenterTitle(true);
    hist->GetZaxis()->SetTitle("Counts");
    hist->GetZaxis()->SetRangeUser(axisLimitZ1, axisLimitZ2);
    hist->GetZaxis()->CenterTitle(true);
    hist->GetZaxis()->SetTickLength(0.0);
    hist->GetZaxis()->SetAxisColor(0);
    hist->Draw("COLZ");
    c->Update();
    c->SaveAs("clusterTriggerMatchingPosition.pdf");
    c->SaveAs("clusterTriggerMatchingPosition.png");
    
    TPaletteAxis *colorbar = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
    double colorbarX1 = 1.0 - (0.27 - 0.021) / 1.2;
    double colorbarX2 = colorbarX1 + 0.05 / 1.2;
    double colorbarY1 = 0.30;
    double colorbarY2 = 0.83;
    colorbar->SetX1NDC(colorbarX1);
    colorbar->SetX2NDC(colorbarX2);
    colorbar->SetY1NDC(colorbarY1);
    colorbar->SetY2NDC(colorbarY2);

    double colorbarAxisX1 = axisLimitX2 + (colorbarX1 - (1.0 - gStyle->GetPadRightMargin())) * (axisLimitX2 - axisLimitX1) / (1.0 - gStyle->GetPadRightMargin() - gStyle->GetPadLeftMargin());
    double colorbarAxisX2 = colorbarAxisX1 + (colorbarX2 - colorbarX1) * (axisLimitX2 - axisLimitX1) / (1.0 - gStyle->GetPadRightMargin() - gStyle->GetPadLeftMargin());
    double colorbarAxisY1 = axisLimitY1 + (colorbarY1 - gStyle->GetPadBottomMargin()) * (axisLimitY2 - axisLimitY1) / (1.0 - gStyle->GetPadTopMargin() - gStyle->GetPadBottomMargin());
    double colorbarAxisY2 = colorbarAxisY1 + (colorbarY2 - colorbarY1) * (axisLimitY2 - axisLimitY1) / (1.0 - gStyle->GetPadTopMargin() - gStyle->GetPadBottomMargin());
    //TGaxis *colorbarAxisLeft = new TGaxis(colorbarAxisX1, colorbarAxisY1, colorbarAxisX1, colorbarAxisY2, axisLimitZ1, axisLimitZ2, 505, "S-");
    //colorbarAxisLeft->SetLabelSize(0.0);
    //colorbarAxisLeft->SetTickLength(0.035 * 1.2);
    //colorbarAxisLeft->Draw();   
    TGaxis *colorbarAxisLeft = new TGaxis(colorbarAxisX1, colorbarAxisY1, colorbarAxisX1, colorbarAxisY2, axisLimitZ1, axisLimitZ2, 100, "S-");
    colorbarAxisLeft->SetTickLength(0.035 * 2.0);
    colorbarAxisLeft->Draw();   
    TGaxis *colorbarAxisRight = new TGaxis(colorbarAxisX2, colorbarAxisY1, colorbarAxisX2, colorbarAxisY2, axisLimitZ1, axisLimitZ2, 505, "S+");
    colorbarAxisRight->SetLabelSize(0.0);
    colorbarAxisRight->SetTickLength(0.035 * 1.2);
    colorbarAxisRight->Draw();
    
    c->Update();
    c->SaveAs("clusterTriggerMatchingPosition.png");
    c->SaveAs("clusterTriggerMatchingPosition.pdf");
  
}
