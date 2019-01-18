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
		rootFileName="~/data/20190115_Fe55_719uA_gdgem_readouts_20190115-165753_00000_bc_20_tac_100_cxys2_cs1_dt200_mst2_spc500_dp200_HITS.root";
	}
	else
	{
		rootFileName = argv[1];
	}

	TFile * f = new TFile(rootFileName);
	
	
	TTree * t = (TTree*)f->Get("events");
	 
	std::vector<CommonClusterNMX> *  m_clustersXY = nullptr;
	std::vector<ClusterNMX> *  m_clusters = nullptr;
		
	CommonClusterNMX theCommonCluster;
	ClusterNMX theCluster;
	t->SetBranchAddress("clusters", &m_clusters);
	t->SetBranchAddress("clustersXY", &m_clustersXY);
	
	    
	int numEvents = t->GetEntries();
	int cnt = 0;
	for (int i = 0; i < numEvents; i++)
	{
		t->GetEntry(i);	
		for( int s = 0; s < m_clustersXY->size(); s++)
		{   
			theCommonCluster = m_clustersXY->at(s);
			int posX = (int)m_clustersXY->at(s).positionX + 1;
			int posY = (int)m_clustersXY->at(s).positionY + 1;
			std::cout << "Common-Cluster ID " << m_clustersXY->at(s).id << ": x=" << posX << ", y=" << posY << std::endl;
			
		}
	}
		
}