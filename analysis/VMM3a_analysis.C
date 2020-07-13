#define VMM3a_cxx
// The class definition in VMM3a.h has been generated automatically
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
// root> T->Process("VMM3a.C")
// root> T->Process("VMM3a.C","some options")
// root> T->Process("VMM3a.C+")
//


#include "VMM3a.h"
#include <TH2.h>
#include <TStyle.h>

void VMM3a::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(1002210);
   gStyle->SetOptDate();
   gStyle->SetOptFit(1);
   gStyle->SetTimeOffset(dtime.Convert());

  view = TView::CreateView(1);
  view->SetRange(0,0,0,100000,10000,113.5575);

  cXHits = new TCanvas("cXHits", "cXHits",120,190,700,500);
  cXHits->Range(-0.75,-0.75,0.75,0.75);
  cXHits->SetFillColor(0);
  cXHits->SetBorderMode(0);
  cXHits->SetBorderSize(2);
  cXHits->SetTheta(90);
  cXHits->SetPhi(-360);
  cXHits->SetFrameBorderMode(0);

  cXQ = new TCanvas("cXQ", "cXQ",120,190,700,500);
  cXQ->Range(-0.75,-0.75,0.75,0.75);
  cXQ->SetFillColor(0);
  cXQ->SetBorderMode(0);
  cXQ->SetBorderSize(2);
  cXQ->SetTheta(90);
  cXQ->SetPhi(-360);
  cXQ->SetFrameBorderMode(0);

  cXNormalizedCharge = new TCanvas("cXNormalizedCharge", "cXNormalizedCharge",120,190,700,500);
  cXNormalizedCharge->Range(-0.75,-0.75,0.75,0.75);
  cXNormalizedCharge->SetFillColor(0);
  cXNormalizedCharge->SetBorderMode(0);
  cXNormalizedCharge->SetBorderSize(2);
  cXNormalizedCharge->SetTheta(90);
  cXNormalizedCharge->SetPhi(-360);
  cXNormalizedCharge->SetFrameBorderMode(0);

  cXCharge = new TCanvas("cXCharge", "cXCharge",120,190,700,500);
  cXCharge->Range(-0.75,-0.75,0.75,0.75);
  cXCharge->SetFillColor(0);
  cXCharge->SetBorderMode(0);
  cXCharge->SetBorderSize(2);
  cXCharge->SetTheta(90);
  cXCharge->SetPhi(-360);
  cXCharge->SetFrameBorderMode(0);
 
  cXChargeVsTime = new TCanvas("cXChargeVsTime", "cXChargeVsTime",120,190,700,500);
  cXChargeVsTime->Range(-0.75,-0.75,0.75,0.75);
  cXChargeVsTime->SetFillColor(0);
  cXChargeVsTime->SetBorderMode(0);
  cXChargeVsTime->SetBorderSize(2);
  cXChargeVsTime->SetTheta(90);
  cXChargeVsTime->SetPhi(-360);
  cXChargeVsTime->SetFrameBorderMode(0);
 
 
  cXChargeVsPosition = new TCanvas("cXChargeVsPosition", "cXChargeVsPosition",120,190,700,500);
  cXChargeVsPosition->Range(-0.75,-0.75,0.75,0.75);
  cXChargeVsPosition->SetFillColor(0);
  cXChargeVsPosition->SetBorderMode(0);
  cXChargeVsPosition->SetBorderSize(2);
  cXChargeVsPosition->SetTheta(90);
  cXChargeVsPosition->SetPhi(-360);
  cXChargeVsPosition->SetFrameBorderMode(0);
 
  cYHits = new TCanvas("cYHits", "cYHits",120,190,700,500);
  cYHits->Range(-0.75,-0.75,0.75,0.75);
  cYHits->SetFillColor(0);
  cYHits->SetBorderMode(0);
  cYHits->SetBorderSize(2);
  cYHits->SetTheta(90);
  cYHits->SetPhi(-360);
  cYHits->SetFrameBorderMode(0);

  cYQ = new TCanvas("cYQ", "cYQ",120,190,700,500);
  cYQ->Range(-0.75,-0.75,0.75,0.75);
  cYQ->SetFillColor(0);
  cYQ->SetBorderMode(0);
  cYQ->SetBorderSize(2);
  cYQ->SetTheta(90);
  cYQ->SetPhi(-360);
  cYQ->SetFrameBorderMode(0);

  cYNormalizedCharge = new TCanvas("cYNormalizedCharge", "cYNormalizedCharge",120,190,700,500);
  cYNormalizedCharge->Range(-0.75,-0.75,0.75,0.75);
  cYNormalizedCharge->SetFillColor(0);
  cYNormalizedCharge->SetBorderMode(0);
  cYNormalizedCharge->SetBorderSize(2);
  cYNormalizedCharge->SetTheta(90);
  cYNormalizedCharge->SetPhi(-360);
  cYNormalizedCharge->SetFrameBorderMode(0);

  cYCharge = new TCanvas("cYCharge", "cYCharge",120,190,700,500);
  cYCharge->Range(-0.75,-0.75,0.75,0.75);
  cYCharge->SetFillColor(0);
  cYCharge->SetBorderMode(0);
  cYCharge->SetBorderSize(2);
  cYCharge->SetTheta(90);
  cYCharge->SetPhi(-360);
  cYCharge->SetFrameBorderMode(0);
 
  cYChargeVsTime = new TCanvas("cXChargeVsTime", "cXChargeVsTime",120,190,700,500);
  cYChargeVsTime->Range(-0.75,-0.75,0.75,0.75);
  cYChargeVsTime->SetFillColor(0);
  cYChargeVsTime->SetBorderMode(0);
  cYChargeVsTime->SetBorderSize(2);
  cYChargeVsTime->SetTheta(90);
  cYChargeVsTime->SetPhi(-360);
  cYChargeVsTime->SetFrameBorderMode(0);

  cYChargeVsPosition = new TCanvas("cYChargeVsPosition", "cYChargeVsPosition",120,190,700,500);
  cYChargeVsPosition->Range(-0.75,-0.75,0.75,0.75);
  cYChargeVsPosition->SetFillColor(0);
  cYChargeVsPosition->SetBorderMode(0);
  cYChargeVsPosition->SetBorderSize(2);
  cYChargeVsPosition->SetTheta(90);
  cYChargeVsPosition->SetPhi(-360);
  cYChargeVsPosition->SetFrameBorderMode(0);
 

  cXYHits = new TCanvas("cXYHits", "cXYHits",120,190,700,500);
  cXYHits->Range(-0.75,-0.75,0.75,0.75);
  cXYHits->SetFillColor(0);
  cXYHits->SetBorderMode(0);
  cXYHits->SetBorderSize(2);
  cXYHits->SetTheta(90);
  cXYHits->SetPhi(-360);
  cXYHits->SetFrameBorderMode(0);

  cXYCharge = new TCanvas("cXYCharge", "cXYCharge",120,190,700,500);
  cXYCharge->Range(-0.75,-0.75,0.75,0.75);
  cXYCharge->SetFillColor(0);
  cXYCharge->SetBorderMode(0);
  cXYCharge->SetBorderSize(2);
  cXYCharge->SetTheta(90);
  cXYCharge->SetPhi(-360);
  cXYCharge->SetFrameBorderMode(0);

  cXYNormalizedCharge = new TCanvas("cXYNormalizedCharge", "cXYNormalizedCharge",120,190,700,500);
  cXYNormalizedCharge->Range(-0.75,-0.75,0.75,0.75);
  cXYNormalizedCharge->SetFillColor(0);
  cXYNormalizedCharge->SetBorderMode(0);
  cXYNormalizedCharge->SetBorderSize(2);
  cXYNormalizedCharge->SetTheta(90);
  cXYNormalizedCharge->SetPhi(-360);
  cXYNormalizedCharge->SetFrameBorderMode(0);


  //--- histograms ---

  hXHits = new TH1F("hXHits","X Hits", 1280, 0., 1279);
  hXQ = new TH1F("hXQ","Q[X]",1280, 0., 1279);
  hXNormalizedCharge = new TH1F("hXNormalizedCharge","X Normalized Charge", 1280, 0., 1279);
  hXCharge = new TH1F("hXCharge","X Spectrum", 100, 0., 3000);
  hXChargeVsTime = new TH2F("hXChargeVsTime","X Peak Position vs Time",100,0.,1000000,100, 0., 3000);
  hXChargeVsPosition = new TH2F("hXChargeVsPosition","X Charge vs Postion",1280, 0., 1279, 100, 0., 3000);

  hYHits = new TH1F("hYHits","Y Hits",1280, 0., 1279);
  hYQ = new TH1F("hYQ","Q[Y]",1280, 0., 1279);
  hYNormalizedCharge = new TH1F("hYNormalizedCharge","Y Normalized Charge", 1280, 0., 1279);
  hYCharge = new TH1F("hYCharge","Y Spectrum", 100, 0., 3000);
  hYChargeVsTime = new TH2F("hYChargeVsTime","Y Peak Position", 100,0.,1000000,100, 0., 3000);
  hYChargeVsPosition = new TH2F("hYChargeVsPosition","Y Charge vs Postion",1280, 0., 1279, 100, 0., 3000);

  hXYHits = new TH2F("hXYHits","XY Hits", 1280, 0., 1279, 1280, 0., 1279);
  hXYCharge = new TH2F("hXYCharge","XY Charge", 1280, 0., 1279, 1280, 0., 1279);
  hXYNormalizedCharge = new TH2F("hXYNormalizedCharge","XY Normalized Charge", 1280, 0., 1279,1280, 0., 1279);

  hXQCLS = new TH2F("hXQCLS","charge_x:size_x",32,0,32,100,0,3000);
  hYQCLS = new TH2F("hYQCLS","charge_y:size_y",32,0,32,100,0,3000);

  bin_content=0.;

  cClusterTimeDifference = new TCanvas("cClusterTimeDifference", "cClusterTimeDifference",120,190,700,500);
  hClusterTimeDifference = new TH1D("hClusterTimeDifference","Cluster Time Differences", 1000,-100,100);

}

void VMM3a::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t VMM3a::Process(Long64_t entry)
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
   //if(entry<100000){
   if(1){
   if (entry%1000==0) cout << "entry: " << entry << "\r" << flush;
   
	   for (int n=0;  n< clusters_detector_pos0.GetSize(); n++) 
	   	{
   		x_bin = TMath::FloorNint(clusters_detector_pos0.At(n));
   		y_bin = TMath::FloorNint(clusters_detector_pos1.At(n));

		hXCharge->Fill(clusters_detector_adc0.At(n));
   		hYCharge->Fill(clusters_detector_adc1.At(n));

		hXChargeVsTime->Fill(entry,clusters_detector_adc0.At(n));
   		hYChargeVsTime->Fill(entry,clusters_detector_adc1.At(n));

		hXChargeVsPosition->Fill(clusters_detector_pos0.At(n),clusters_detector_adc0.At(n));
   		hYChargeVsPosition->Fill(clusters_detector_pos1.At(n),clusters_detector_adc1.At(n));

   		hXHits->Fill(x_bin);
   		hYHits->Fill(y_bin);
   		hXYHits->Fill(x_bin,y_bin);

   		bin_content = hXQ->GetBinContent(x_bin);
   		hXQ->SetBinContent(x_bin,bin_content+clusters_detector_adc0.At(n));

   		bin_content = hYQ->GetBinContent(y_bin);
   		hYQ->SetBinContent(y_bin,bin_content+clusters_detector_adc1.At(n));

   		bin_content = hXYCharge->GetBinContent(x_bin,y_bin);
   		hXYCharge->SetBinContent(x_bin,y_bin,bin_content+clusters_detector_adc0.At(n)+clusters_detector_adc1.At(n));

   		hXQCLS->Fill(clusters_detector_size0.At(n),clusters_detector_adc0.At(n));
   		hYQCLS->Fill(clusters_detector_size1.At(n),clusters_detector_adc1.At(n));


		hClusterTimeDifference->Fill(clusters_detector_time0.At(n) - clusters_detector_time1.At(n));

        	}
   


   return kTRUE;
   }
   else return kFALSE;
}

void VMM3a::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void VMM3a::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   cout << "graphing...."<< endl;

   cXCharge->cd();
   hXCharge->Draw();

   cXChargeVsTime->cd();
   hXChargeVsTime->Draw("lego2");

   cXChargeVsPosition->cd();
   hXChargeVsPosition->Draw("lego2");

   cYCharge->cd();
   hYCharge->Draw();

   cYChargeVsTime->cd();
   hYChargeVsTime->Draw("lego2");

   cYChargeVsPosition->cd();
   hYChargeVsPosition->Draw("lego2");

   cXYHits->cd();
   hXYHits->Draw("lego2");

   cXYCharge->cd();
   hXYCharge->Draw("lego2");

   cXYNormalizedCharge->cd();
   hXYNormalizedCharge->Divide(hXYCharge,hXYHits,1,1);
   //hXYNormalizedCharge->GetXaxis()->SetRangeUser(-48,48);
   //hXYNormalizedCharge->GetYaxis()->SetRangeUser(-48,48);
   hXYNormalizedCharge->Draw("SURF1");

   cXQ->cd();
   hXQ->Draw();

   cYQ->cd();
   hYQ->Draw();

   cXHits->cd();
   hXHits->Draw();

   cYHits->cd();
   hYHits->Draw();

   cXNormalizedCharge->cd();
   hXNormalizedCharge->Divide(hXQ,hXHits,1,1);
   //hXNormalizedCharge->GetXaxis()->SetRangeUser(-48,48);
   hXNormalizedCharge->Draw();

   cYNormalizedCharge->cd();
   hYNormalizedCharge->Divide(hYQ,hYHits,1,1);
   //hYNormalizedCharge->GetXaxis()->SetRangeUser(-48,48);
   hYNormalizedCharge->Draw();


   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");


   TCanvas *cXQCLS = new TCanvas("cXQCLS", "cXQCLS",375,350,680,520);
   cXQCLS->Range(-1.0,-1.25,1.0,1.25);
   TView *view_1 = TView::CreateView(1);
   view_1->SetRange(0,0,0,20,6800,10652.36);
   cXQCLS->SetFillColor(0);
   cXQCLS->SetBorderMode(0);
   cXQCLS->SetBorderSize(2);
   cXQCLS->SetTheta(57.80679);
   cXQCLS->SetPhi(-405.4827);
   cXQCLS->SetFrameBorderMode(0);
   cXQCLS->cd();

   hXQCLS->SetLineColor(ci);
   hXQCLS->GetXaxis()->SetTitle("size_x");
   hXQCLS->GetXaxis()->SetRange(1,20);
   hXQCLS->GetXaxis()->SetLabelFont(42);
   hXQCLS->GetXaxis()->SetLabelSize(0.035);
   hXQCLS->GetXaxis()->SetTitleSize(0.035);
   hXQCLS->GetXaxis()->SetTitleFont(42);
   hXQCLS->GetYaxis()->SetTitle("charge_x");
   hXQCLS->GetYaxis()->SetRange(1,200);
   hXQCLS->GetYaxis()->SetLabelFont(42);
   hXQCLS->GetYaxis()->SetLabelSize(0.035);
   hXQCLS->GetYaxis()->SetTitleSize(0.035);
   hXQCLS->GetYaxis()->SetTitleFont(42);
   hXQCLS->GetZaxis()->SetLabelFont(42);
   hXQCLS->GetZaxis()->SetLabelSize(0.035);
   hXQCLS->GetZaxis()->SetTitleSize(0.035);
   hXQCLS->GetZaxis()->SetTitleFont(42);
   hXQCLS->Draw("lego2");

   TPaveText *ptx = new TPaveText(0.31,0.94,0.69,0.995,"blNDC");
   ptx->SetName("title");
   ptx->SetBorderSize(0);
   ptx->SetFillColor(0);
   ptx->SetFillStyle(0);
   ptx->SetTextFont(42);
   TText *textx = ptx->AddText("charge_x:size_x");
   ptx->Draw();
   cXQCLS->Modified();
   cXQCLS->SetSelected(cXQCLS);

   TCanvas *cYQCLS = new TCanvas("cYQCLS", "cYQCLS",375,350,680,520);
   cYQCLS->Range(-1.0,-1.25,1.0,1.25);
   TView *view_2 = TView::CreateView(1);
   view_2->SetRange(0,0,0,20,6800,10652.36);
   cYQCLS->SetFillColor(0);
   cYQCLS->SetBorderMode(0);
   cYQCLS->SetBorderSize(2);
   cYQCLS->SetTheta(57.80679);
   cYQCLS->SetPhi(-405.4827);
   cYQCLS->SetFrameBorderMode(0);
   cYQCLS->cd();

   hYQCLS->SetLineColor(ci);
   hYQCLS->GetXaxis()->SetTitle("size_y");
   hYQCLS->GetXaxis()->SetRange(1,20);
   hYQCLS->GetXaxis()->SetLabelFont(42);
   hYQCLS->GetXaxis()->SetLabelSize(0.035);
   hYQCLS->GetXaxis()->SetTitleSize(0.035);
   hYQCLS->GetXaxis()->SetTitleFont(42);
   hYQCLS->GetYaxis()->SetTitle("charge_y");
   hYQCLS->GetYaxis()->SetRange(1,200);
   hYQCLS->GetYaxis()->SetLabelFont(42);
   hYQCLS->GetYaxis()->SetLabelSize(0.035);
   hYQCLS->GetYaxis()->SetTitleSize(0.035);
   hYQCLS->GetYaxis()->SetTitleFont(42);
   hYQCLS->GetZaxis()->SetLabelFont(42);
   hYQCLS->GetZaxis()->SetLabelSize(0.035);
   hYQCLS->GetZaxis()->SetTitleSize(0.035);
   hYQCLS->GetZaxis()->SetTitleFont(42);
   hYQCLS->Draw("lego2");

   TPaveText *pty = new TPaveText(0.31,0.94,0.69,0.995,"blNDC");
   pty->SetName("title");
   pty->SetBorderSize(0);
   pty->SetFillColor(0);
   pty->SetFillStyle(0);
   pty->SetTextFont(42);
   TText *texty = pty->AddText("charge_y:size_y");
   pty->Draw();
   cYQCLS->Modified();
   cYQCLS->SetSelected(cXQCLS);

   cClusterTimeDifference->cd();
   hClusterTimeDifference->Draw();


}
