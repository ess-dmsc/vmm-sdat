//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun 10 15:09:15 2020 by ROOT version 6.18/04
// from TTree events/vmm3 events
// found on file: example.root
//////////////////////////////////////////////////////////

#ifndef VMM3a_Analysis_h
#define VMM3a_Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <TCanvas.h>


// Headers needed by this particular selector
#include <vector>



class VMM3a_Analysis : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain


   TH1D *hClusterTimeDifference;
   TCanvas *cClusterTimeDifference;

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<unsigned int> hits_id = {fReader, "hits.id"};
   TTreeReaderArray<unsigned int> hits_event = {fReader, "hits.event"};
   TTreeReaderArray<UChar_t> hits_det = {fReader, "hits.det"};
   TTreeReaderArray<UChar_t> hits_plane = {fReader, "hits.plane"};
   TTreeReaderArray<UChar_t> hits_fec = {fReader, "hits.fec"};
   TTreeReaderArray<UChar_t> hits_vmm = {fReader, "hits.vmm"};
   TTreeReaderArray<Double_t> hits_readout_time = {fReader, "hits.readout_time"};
   TTreeReaderArray<Double_t> hits_time = {fReader, "hits.time"};
   TTreeReaderArray<UChar_t> hits_ch = {fReader, "hits.ch"};
   TTreeReaderArray<unsigned short> hits_pos = {fReader, "hits.pos"};
   TTreeReaderArray<unsigned short> hits_bcid = {fReader, "hits.bcid"};
   TTreeReaderArray<unsigned short> hits_tdc = {fReader, "hits.tdc"};
   TTreeReaderArray<unsigned short> hits_adc = {fReader, "hits.adc"};
   TTreeReaderArray<Bool_t> hits_over_threshold = {fReader, "hits.over_threshold"};
   TTreeReaderArray<Float_t> hits_chip_time = {fReader, "hits.chip_time"};
   TTreeReaderArray<unsigned int> clusters_plane_id = {fReader, "clusters_plane.id"};
   TTreeReaderArray<UChar_t> clusters_plane_det = {fReader, "clusters_plane.det"};
   TTreeReaderArray<UChar_t> clusters_plane_plane = {fReader, "clusters_plane.plane"};
   TTreeReaderArray<unsigned short> clusters_plane_size = {fReader, "clusters_plane.size"};
   TTreeReaderArray<unsigned int> clusters_plane_adc = {fReader, "clusters_plane.adc"};
   TTreeReaderArray<Double_t> clusters_plane_time = {fReader, "clusters_plane.time"};
   TTreeReaderArray<Double_t> clusters_plane_time_utpc = {fReader, "clusters_plane.time_utpc"};
   TTreeReaderArray<Double_t> clusters_plane_time_charge2 = {fReader, "clusters_plane.time_charge2"};
   TTreeReaderArray<Double_t> clusters_plane_pos = {fReader, "clusters_plane.pos"};
   TTreeReaderArray<Double_t> clusters_plane_pos_utpc = {fReader, "clusters_plane.pos_utpc"};
   TTreeReaderArray<Double_t> clusters_plane_pos_charge2 = {fReader, "clusters_plane.pos_charge2"};
   TTreeReaderArray<Bool_t> clusters_plane_plane_coincidence = {fReader, "clusters_plane.plane_coincidence"};
   TTreeReaderArray<unsigned short> clusters_plane_max_delta_time = {fReader, "clusters_plane.max_delta_time"};
   TTreeReaderArray<unsigned short> clusters_plane_max_missing_strip = {fReader, "clusters_plane.max_missing_strip"};
   TTreeReaderArray<unsigned short> clusters_plane_span_cluster = {fReader, "clusters_plane.span_cluster"};
   TTreeReaderArray<vector<double>> clusters_plane_strips = {fReader, "clusters_plane.strips"};
   TTreeReaderArray<vector<double>> clusters_plane_times = {fReader, "clusters_plane.times"};
   TTreeReaderArray<unsigned int> clusters_detector_id = {fReader, "clusters_detector.id"};
   TTreeReaderArray<unsigned int> clusters_detector_id0 = {fReader, "clusters_detector.id0"};
   TTreeReaderArray<unsigned int> clusters_detector_id1 = {fReader, "clusters_detector.id1"};
   TTreeReaderArray<UChar_t> clusters_detector_det = {fReader, "clusters_detector.det"};
   TTreeReaderArray<unsigned short> clusters_detector_size0 = {fReader, "clusters_detector.size0"};
   TTreeReaderArray<unsigned short> clusters_detector_size1 = {fReader, "clusters_detector.size1"};
   TTreeReaderArray<unsigned int> clusters_detector_adc0 = {fReader, "clusters_detector.adc0"};
   TTreeReaderArray<unsigned int> clusters_detector_adc1 = {fReader, "clusters_detector.adc1"};
   TTreeReaderArray<Double_t> clusters_detector_pos0 = {fReader, "clusters_detector.pos0"};
   TTreeReaderArray<Double_t> clusters_detector_pos1 = {fReader, "clusters_detector.pos1"};
   TTreeReaderArray<Double_t> clusters_detector_pos2 = {fReader, "clusters_detector.pos2"};
   TTreeReaderArray<Double_t> clusters_detector_time0 = {fReader, "clusters_detector.time0"};
   TTreeReaderArray<Double_t> clusters_detector_time1 = {fReader, "clusters_detector.time1"};
   TTreeReaderArray<Double_t> clusters_detector_pos0_utpc = {fReader, "clusters_detector.pos0_utpc"};
   TTreeReaderArray<Double_t> clusters_detector_pos1_utpc = {fReader, "clusters_detector.pos1_utpc"};
   TTreeReaderArray<Double_t> clusters_detector_pos2_utpc = {fReader, "clusters_detector.pos2_utpc"};
   TTreeReaderArray<Double_t> clusters_detector_time0_utpc = {fReader, "clusters_detector.time0_utpc"};
   TTreeReaderArray<Double_t> clusters_detector_time1_utpc = {fReader, "clusters_detector.time1_utpc"};
   TTreeReaderArray<Double_t> clusters_detector_pos0_charge2 = {fReader, "clusters_detector.pos0_charge2"};
   TTreeReaderArray<Double_t> clusters_detector_pos1_charge2 = {fReader, "clusters_detector.pos1_charge2"};
   TTreeReaderArray<Double_t> clusters_detector_pos2_charge2 = {fReader, "clusters_detector.pos2_charge2"};
   TTreeReaderArray<Double_t> clusters_detector_time0_charge2 = {fReader, "clusters_detector.time0_charge2"};
   TTreeReaderArray<Double_t> clusters_detector_time1_charge2 = {fReader, "clusters_detector.time1_charge2"};
   TTreeReaderArray<Double_t> clusters_detector_dt0 = {fReader, "clusters_detector.dt0"};
   TTreeReaderArray<Double_t> clusters_detector_dt1 = {fReader, "clusters_detector.dt1"};
   TTreeReaderArray<Double_t> clusters_detector_delta_plane = {fReader, "clusters_detector.delta_plane"};
   TTreeReaderArray<unsigned short> clusters_detector_span_cluster0 = {fReader, "clusters_detector.span_cluster0"};
   TTreeReaderArray<unsigned short> clusters_detector_span_cluster1 = {fReader, "clusters_detector.span_cluster1"};
   TTreeReaderArray<unsigned short> clusters_detector_max_delta_time0 = {fReader, "clusters_detector.max_delta_time0"};
   TTreeReaderArray<unsigned short> clusters_detector_max_delta_time1 = {fReader, "clusters_detector.max_delta_time1"};
   TTreeReaderArray<unsigned short> clusters_detector_max_missing_strip0 = {fReader, "clusters_detector.max_missing_strip0"};
   TTreeReaderArray<unsigned short> clusters_detector_max_missing_strip1 = {fReader, "clusters_detector.max_missing_strip1"};
   TTreeReaderArray<vector<double>> clusters_detector_strips0 = {fReader, "clusters_detector.strips0"};
   TTreeReaderArray<vector<double>> clusters_detector_times0 = {fReader, "clusters_detector.times0"};
   TTreeReaderArray<vector<double>> clusters_detector_strips1 = {fReader, "clusters_detector.strips1"};
   TTreeReaderArray<vector<double>> clusters_detector_times1 = {fReader, "clusters_detector.times1"};
 


   VMM3a_Analysis(TTree * /*tree*/ =0) { }
   virtual ~VMM3a_Analysis() { }
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

   ClassDef(VMM3a_Analysis,0);

};

#endif

#ifdef VMM3a_Analysis_cxx
void VMM3a_Analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t VMM3a_Analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef VMM3a_Analysis_cxx
