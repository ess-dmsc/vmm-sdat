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
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include <TObject.h>
#include <TString.h>
#include "DataStructures.h"



int main(int argc, char**argv)
{
 	TString rootFileName = "";
	if(argc != 2)
	{
                rootFileName="../example.root";
	}
	else
	{
		rootFileName = argv[1];
	}

	TFile * f = new TFile(rootFileName);
	
	
	TTree * t = (TTree*)f->Get("events");
	 
	std::vector<ClusterDetector> *  m_clusters_detector = nullptr;
	std::vector<ClusterPlane> *  m_clusters_plane = nullptr;
		
	ClusterDetector theClusterDetector;
	ClusterPlane theClusterPlane;
	t->SetBranchAddress("clusters_plane", &m_clusters_plane);
	t->SetBranchAddress("clusters_detector", &m_clusters_detector);
	
	int numEvents = t->GetEntries();   
	int cnt = 0;
	std::vector<double> vTime;
	double oldTime0 = 0;
	double time0 = 0;
	for (int i = 0; i < numEvents; i++)
	{
		t->GetEntry(i);	
		for( int s = 0; s < m_clusters_detector->size(); s++)
		{   
			theClusterDetector = m_clusters_detector->at(s);
			oldTime0 = time0;
			time0 = (double)m_clusters_detector->at(s).time0;
			if(time0 < oldTime0)
				std::cout << "Time error in event: " << i << ": Time cluster n " << time0 << " - time cluster n-1 " << oldTime0 << std::endl;
			vTime.push_back(time0);
		}
		if(i%100000 == 0)
		{	
			std::cout << "Event " << i << std::endl;
		}

	}
	TCanvas *c1 = new TCanvas("c1","c1");
	c1->cd();
	TH1D *h1 = new TH1D("dt0", "dt0", 10000,0,100000);
	std::cout << "Sorting cluster time vector.." << std::endl;
	 std::sort(begin(vTime), end(vTime));
	 oldTime0 = 0;
	 time0 = 0;
	 double dt0 = 0;
	 for(auto &t: vTime)
	 {
		oldTime0 = time0;
		time0 = t;
		dt0 = time0-oldTime0;
		h1->Fill(dt0);

	 }
	h1->Draw("HIST");		
	c1->SaveAs("dt0.pdf");
}
