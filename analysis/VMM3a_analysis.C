#define VMM3a_Analysis_cxx
// The class definition in VMM3a_Analysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("VMM3a_Analysis.C")
// root> T->Process("VMM3a_Analysis.C","some options")
// root> T->Process("VMM3a_Analysis.C+")
//


#include "VMM3a_Analysis.h"
#include <TH2.h>
#include <TStyle.h>

void VMM3a_Analysis::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
    
   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(1002210);
   gStyle->SetOptDate();
   gStyle->SetOptFit(1);
  
   cClusterTimeDifference = new TCanvas("cClusterTimeDifference", "cClusterTimeDifference",120,190,700,500);
   hClusterTimeDifference = new TH1D("hClusterTimeDifference","Cluster Time Differences", 1000,-100,100);
}

void VMM3a_Analysis::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

Bool_t VMM3a_Analysis::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);
   
 
  	for (int n=0;  n< clusters_detector_time0.GetSize(); n++) {
    	hClusterTimeDifference->Fill(clusters_detector_time0.At(n) - clusters_detector_time1.At(n));
  	}


   return kTRUE;
}

void VMM3a_Analysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void VMM3a_Analysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
	cClusterTimeDifference->cd();
   	hClusterTimeDifference->Draw();

}