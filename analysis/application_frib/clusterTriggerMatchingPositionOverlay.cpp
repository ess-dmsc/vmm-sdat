// Cluster-trigger-matching
// -------------------------------------------------------------------
// lucian.scharenberg@cern.ch
// 06/Oct/2022 and 04, 05 and 07/Oct/2023

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
    double dtMax = 7000.0;
    
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


void clusterTriggerMatchingPositionOverlay() {

    vector<string> runs = {"./run3.root",
			   "./run7.root",
			   "./run2.root"};
    vector<vector<double>> results;

    for (int r = 0; r < runs.size(); r++) {
    
	vector<double> posClusters, posTrigger, timeDifference;

	TFile *f = new TFile(runs[r].c_str());
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

	results.push_back(timeDifference);

    }


    // Plot the data
    gROOT->LoadMacro("./fellowStyle.cpp");
    gROOT->ProcessLine("fellowStyle(0)");
    double w = 1000;
    double h = 1000;
    TCanvas *c = new TCanvas("c", "", w, h);
    c->SetWindowSize(w + (w - c->GetWw()), h + (h - c->GetWh()));
    gROOT->ForceStyle();

    TH1D *hist1 = new TH1D("hist1", "", 100, -5.0, 5.0);
    TH1D *hist2 = new TH1D("hist2", "", 100, -5.0, 5.0);
    TH1D *hist3 = new TH1D("hist3", "", 100, -5.0, 5.0);

    for (int i = 0; i < results[0].size(); i++) {

	hist1->Fill(results[0][i] / 1000.0);

    }

    for (int i = 0; i < results[1].size(); i++) {

	hist2->Fill(results[1][i] / 1000.0);

    }

    for (int i = 0; i < results[2].size(); i++) {

	hist3->Fill(results[2][i] / 1000.0);

    }   
    
    hist1->SetLineColor(kAzure+2);
    hist2->SetLineColor(kPink-1);
    hist3->SetLineColor(kOrange-3);
    hist1->GetXaxis()->SetTitle("Time difference / us");
    hist1->GetXaxis()->SetLimits(-5.0, 5.0);
    hist1->GetXaxis()->CenterTitle(true);
    hist1->GetYaxis()->SetTitle("Counts");
    hist1->GetYaxis()->SetRangeUser(0.0, 200.0);
    hist1->GetYaxis()->CenterTitle(true);
    hist1->Draw("HIST");
    hist2->Draw("HIST SAME");
    hist3->Draw("HIST SAME");

    TLatex t1;
    t1.SetTextAlign(23);
    t1.SetTextColor(kBlack);
    t1.SetTextFont(42);
    t1.SetTextSize(0.055);
    int binMax = hist1->GetMaximumBin();
    t1.DrawLatexNDC(gStyle->GetPadLeftMargin() + ((hist1->GetXaxis()->GetBinCenter(binMax) + 5.0) / 10.0) * (1.0 - gStyle->GetPadRightMargin() - gStyle->GetPadLeftMargin()),
		    gStyle->GetPadBottomMargin() + (hist1->GetBinContent(binMax) / 200.0) * (1.0 - gStyle->GetPadTopMargin() - gStyle->GetPadBottomMargin()) + 0.05,
		    "15 cm");

    cout << gStyle->GetPadLeftMargin()
	 << " "
	 << gStyle->GetPadBottomMargin()
	 << " "
	 << gStyle->GetPadRightMargin()
	 << " "
	 << gStyle->GetPadTopMargin()
	 << endl;

    cout << binMax
	 << endl;

    cout << hist1->GetXaxis()->GetBinCenter(binMax)
	 << endl;

    cout << hist1->GetBinContent(binMax)
	 << endl;
    
    cout << gStyle->GetPadLeftMargin() + (((hist1->GetXaxis()->GetBinCenter(binMax) + 5.0) / 10.0) * (1.0 - gStyle->GetPadRightMargin() - gStyle->GetPadLeftMargin()))
	 << " "
	 << gStyle->GetPadBottomMargin() + ((hist1->GetBinContent(binMax) / 200.0) * (1.0 - gStyle->GetPadTopMargin() - gStyle->GetPadBottomMargin())) + 0.05
	 << endl;

    TLatex t2;
    t2.SetTextAlign(23);
    t2.SetTextColor(kBlack);
    t2.SetTextFont(42);
    t2.SetTextSize(0.055);
    binMax = hist2->GetMaximumBin();
    t2.DrawLatexNDC(gStyle->GetPadLeftMargin() + ((hist2->GetXaxis()->GetBinCenter(binMax) + 5.0) / 10.0) * (1.0 - gStyle->GetPadRightMargin() - gStyle->GetPadLeftMargin()),
		    gStyle->GetPadBottomMargin() + (hist2->GetBinContent(binMax) / 200.0) * (1.0 - gStyle->GetPadTopMargin() - gStyle->GetPadBottomMargin()) + 0.05,
		     "10 cm");

    TLatex t3;
    t3.SetTextAlign(23);
    t3.SetTextColor(kBlack);
    t3.SetTextFont(42);
    t3.SetTextSize(0.055);
    binMax = hist3->GetMaximumBin();
    t3.DrawLatexNDC(gStyle->GetPadLeftMargin() + ((hist3->GetXaxis()->GetBinCenter(binMax) + 5.0) / 10.0) * (1.0 - gStyle->GetPadRightMargin() - gStyle->GetPadLeftMargin()),
		    gStyle->GetPadBottomMargin() + (hist3->GetBinContent(binMax) / 200.0) * (1.0 - gStyle->GetPadTopMargin() - gStyle->GetPadBottomMargin()) + 0.05,
		    "5 cm");
    
    c->Update();
    c->SaveAs("clusterTriggerMatchingPositionOverlay.png");
    c->SaveAs("clusterTriggerMatchingPositionOverlay.pdf");
  
}
