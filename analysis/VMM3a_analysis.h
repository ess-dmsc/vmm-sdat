//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar  4 14:01:56 2022 by ROOT version 6.24/06
// from TTree clusters_detector/clusters detector
// found on file: example.root
//////////////////////////////////////////////////////////

#ifndef VMM3a_analysis_h
#define VMM3a_analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TView.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class VMM3a_analysis : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   TH2F *hXQCLS;
   TH2F *hYQCLS;

   TH1F *hXHits;
   TH1F *hXQ;
   TH1F *hXNormalizedCharge;
   TH1F *hXCharge;
   TH2F *hXChargeVsTime;
   TH2F *hXChargeVsPosition;

   TH1F *hYHits;
   TH1F *hYQ;
   TH1F *hYNormalizedCharge;
   TH1F *hYCharge;
   TH2F *hYChargeVsTime;
   TH2F *hYChargeVsPosition;

   TH2F *hXYCharge;
   TH2F *hXYNormalizedCharge;
   TH2F *hXYHits;

   TDatime dtime;
   Double_t bin_content;
   Int_t x_bin;
   Int_t y_bin;

   TCanvas *cXQ;
   TCanvas *cYQ;
   TCanvas *cXNormalizedCharge;
   TCanvas *cYNormalizedCharge;
   TCanvas *cXHits;
   TCanvas *cYHits;

   TCanvas *cXCharge;
   TCanvas *cYCharge;
   TCanvas *cXChargeVsTime;
   TCanvas *cYChargeVsTime;
   TCanvas *cXChargeVsPosition;
   TCanvas *cYChargeVsPosition;
   TCanvas *cXYCharge;
   TCanvas *cXYNormalizedCharge;
   TCanvas *cXYHits;
   TView *view;

   TH1D *hClusterTimeDifference;
   TCanvas *cClusterTimeDifference;
   
   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<unsigned int> id = {fReader, "id"};
   TTreeReaderValue<unsigned int> id0 = {fReader, "id0"};
   TTreeReaderValue<unsigned int> id1 = {fReader, "id1"};
   TTreeReaderValue<UChar_t> det = {fReader, "det"};
   TTreeReaderValue<unsigned short> size0 = {fReader, "size0"};
   TTreeReaderValue<unsigned short> size1 = {fReader, "size1"};
   TTreeReaderValue<unsigned short> adc0 = {fReader, "adc0"};
   TTreeReaderValue<unsigned short> adc1 = {fReader, "adc1"};
   TTreeReaderValue<Double_t> pos0 = {fReader, "pos0"};
   TTreeReaderValue<Double_t> pos1 = {fReader, "pos1"};
   TTreeReaderValue<Double_t> pos2 = {fReader, "pos2"};
   TTreeReaderValue<Double_t> time0 = {fReader, "time0"};
   TTreeReaderValue<Double_t> time1 = {fReader, "time1"};
   TTreeReaderValue<Double_t> pos0_utpc = {fReader, "pos0_utpc"};
   TTreeReaderValue<Double_t> pos1_utpc = {fReader, "pos1_utpc"};
   TTreeReaderValue<Double_t> pos2_utpc = {fReader, "pos2_utpc"};
   TTreeReaderValue<Double_t> time0_utpc = {fReader, "time0_utpc"};
   TTreeReaderValue<Double_t> time1_utpc = {fReader, "time1_utpc"};
   TTreeReaderValue<Double_t> pos0_charge2 = {fReader, "pos0_charge2"};
   TTreeReaderValue<Double_t> pos1_charge2 = {fReader, "pos1_charge2"};
   TTreeReaderValue<Double_t> pos2_charge2 = {fReader, "pos2_charge2"};
   TTreeReaderValue<Double_t> time0_charge2 = {fReader, "time0_charge2"};
   TTreeReaderValue<Double_t> time1_charge2 = {fReader, "time1_charge2"};
   TTreeReaderValue<Double_t> pos0_algo = {fReader, "pos0_algo"};
   TTreeReaderValue<Double_t> pos1_algo = {fReader, "pos1_algo"};
   TTreeReaderValue<Double_t> pos2_algo = {fReader, "pos2_algo"};
   TTreeReaderValue<Double_t> time0_algo = {fReader, "time0_algo"};
   TTreeReaderValue<Double_t> time1_algo = {fReader, "time1_algo"};
   TTreeReaderValue<Double_t> dt0 = {fReader, "dt0"};
   TTreeReaderValue<Double_t> dt1 = {fReader, "dt1"};
   TTreeReaderValue<Double_t> delta_plane = {fReader, "delta_plane"};
   TTreeReaderValue<unsigned short> span_cluster0 = {fReader, "span_cluster0"};
   TTreeReaderValue<unsigned short> span_cluster1 = {fReader, "span_cluster1"};
   TTreeReaderValue<unsigned short> max_delta_time0 = {fReader, "max_delta_time0"};
   TTreeReaderValue<unsigned short> max_delta_time1 = {fReader, "max_delta_time1"};
   TTreeReaderValue<unsigned short> max_missing_strip0 = {fReader, "max_missing_strip0"};
   TTreeReaderValue<unsigned short> max_missing_strip1 = {fReader, "max_missing_strip1"};
   TTreeReaderArray<double> strips0 = {fReader, "strips0"};
   TTreeReaderArray<double> times0 = {fReader, "times0"};
   TTreeReaderArray<double> strips1 = {fReader, "strips1"};
   TTreeReaderArray<double> times1 = {fReader, "times1"};
   TTreeReaderArray<double> adcs0 = {fReader, "adcs0"};
   TTreeReaderArray<double> adcs1 = {fReader, "adcs1"};


   VMM3a_analysis(TTree * /*tree*/ =0) { }
   virtual ~VMM3a_analysis() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(VMM3a_analysis,0);

};

#endif

#ifdef VMM3a_analysis_cxx
void VMM3a_analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t VMM3a_analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef VMM3a_analysis_cxx
