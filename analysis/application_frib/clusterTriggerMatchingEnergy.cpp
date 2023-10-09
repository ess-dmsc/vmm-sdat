// Cluster-trigger-matching
// -------------------------------------------------------------------
// lucian.scharenberg@cern.ch
// 06/Oct/2022 and 04/Oct/2023

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


void clusterTriggerMatchingEnergy() {

    vector<double> adcClusters, adcTrigger;
  
    // Open root file
    TFile *f = new TFile("./run2.root");
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

	    events.push_back(make_pair(hitsTime, hitsADC));

	}

    }

    cout << events.size() << endl;
    
    char clustersDet;
    double clustersTime;
    unsigned short clustersADC;
    treeClusters->SetBranchAddress("clusters_plane.det", &clustersDet);
    treeClusters->SetBranchAddress("clusters_plane.time", &clustersTime);
    treeClusters->SetBranchAddress("clusters_plane.adc", &clustersADC);

    for (int i = 0; i < treeClusters->GetEntries(); i++) {

	treeClusters->GetEntry(i);

	if ((int)clustersDet == 0) {

	    clusters.push_back(make_pair(clustersTime, clustersADC));

	}

    }

    // Match the events
    sort(clusters.begin(), clusters.end());
    sort(events.begin(), events.end());

    pair<vector<pair<double, double>>, vector<pair<double, double>>> match = matching(events, clusters);

    for (int i = 0; i < match.first.size(); i++) {
	
	adcClusters.push_back(match.second[i].second);
	adcTrigger.push_back(match.first[i].second);
	
    }

    // Plot the data
    gROOT->LoadMacro("./fellowStyle.cpp");
    gROOT->ProcessLine("fellowStyle(0)");
    double w = 1000;
    double h = 1000;
    TCanvas *c = new TCanvas("c", "", w, h);
    c->SetWindowSize(w + (w - c->GetWw()), h + (h - c->GetWh()));
    //gPad->SetLogy();
    gROOT->ForceStyle();

    TH1D *histFull = new TH1D("histFull", "", 150, 0.0, 15.0);
    TH1D *histMatched = new TH1D("histMatched", "", 150, 0.0, 15.0);

    double sizeFull = clusters.size();
    double sizeMatched = match.first.size();
    cout << "FULL: " << sizeFull << ", MATCHED: " << sizeMatched << endl;
    
    for (int i = 0; i < sizeFull; i++) {

	histFull->Fill(clusters[i].second / 1000.0);
	
    }

    for (int i = 0; i < sizeMatched; i++) {

	histMatched->Fill(adcClusters[i] / 1000.0, 100.0);
	
    }

    histFull->SetLineColor(kAzure+2);
    histMatched->SetLineColor(kPink-1);
    histFull->GetXaxis()->SetTitle("Energy loss in detector / a.u.");
    histFull->GetXaxis()->SetLimits(0.0, 15.0);
    histFull->GetXaxis()->CenterTitle(true);
    //histFull->GetYaxis()->SetTitle("Counts");
    //histFull->GetYaxis()->SetRangeUser(1.0, 1.0e6);
    histFull->GetYaxis()->SetTitle("Counts");
    histFull->GetYaxis()->SetRangeUser(0.0, 1.2e4);
    histFull->GetYaxis()->CenterTitle(true);
    histFull->Draw("HIST");
    histMatched->Draw("HIST SAME");

    int legendNentries = 2;
    string longestEntryString = "Triggered (×100)";
    auto *longestEntry = new TLatex(0.0, 0.0, longestEntryString.c_str());
    longestEntry->SetTextSize(0.055);
    longestEntry->SetTextColor(0);
    unsigned int textExtentWidth, textExtentHeight;
    longestEntry->GetTextExtent(textExtentWidth, textExtentHeight, longestEntryString.c_str());
    double legendWidth = (double)textExtentWidth / w + 3 * 0.015 + 0.07; // textwidth + spaces to frame + margin
    cout << legendWidth << endl;
    double legendHeight = legendNentries * 0.07; // 0.055 textsize + 0.015 labeloffset
    double legendX = 0.0 + gStyle->GetPadLeftMargin() + 0.05;
    double legendY = 1.0 - gStyle->GetPadTopMargin() - 0.05 - legendHeight;
    auto legend = new TLegend(legendX, legendY, legendX + legendWidth, legendY + legendHeight, "", "NDC");
    legend->SetTextSize(0.055);
    legend->SetTextFont(42);
    legend->SetBorderSize(0.0);
    legend->SetMargin(0.07 / legendWidth * 1.0 / (1.0 - gStyle->GetPadLeftMargin() - gStyle->GetPadRightMargin()));
    legend->AddEntry(histFull, "All clusters", "l");
    legend->AddEntry(histMatched, "Triggered (#times 100)", "l");
    legend->Draw();
    
    c->Update();
    c->SaveAs("clusterTriggerMatchingEnergy.png");
    c->SaveAs("clusterTriggerMatchingEnergy.pdf");
  
}
