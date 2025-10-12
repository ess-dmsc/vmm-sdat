// sort_tree.C
// Usage:
//   root -l
// 

#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TError.h>

#include <vector>
#include <utility>
#include <algorithm>
#include <cstdio>
#include <cstring>

// Try to get "clusters_detector"; if missing, pick the first TTree in the file
static TTree* findTree(TFile* f) {
    if (auto* t = dynamic_cast<TTree*>(f->Get("clusters_detector"))) return t;
    TIter nextkey(f->GetListOfKeys());
    while (auto* key = dynamic_cast<TKey*>(nextkey())) {
        if (!strcmp(key->GetClassName(), "TTree"))
            return dynamic_cast<TTree*>(key->ReadObj());
    }
    return nullptr;
}

int sort_tree(const char* infile, const char* outfile, const char* algo, const char* choice, int event_limit)
{
    // --- Open input and find tree
    TFile* fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) { printf("Error: cannot open %s\n", infile); return 1; }
    TTree* tree = findTree(fin);
    if (!tree) { printf("Error: no TTree found in %s\n", infile); return 1; }

    const Long64_t nentries = tree->GetEntries();
    printf("Input tree has %lld entries\n", nentries);

    // --- Pass 1: build (time reference, entry) index using a light reader
    std::vector<std::pair<double, Long64_t>> index;
    index.reserve(nentries);

    {
        TTreeReader rTime(tree);
        TTreeReaderValue<Double_t> time0(rTime, "time0");
        /*
        TTreeReaderValue<Double_t> time1(rTime, "time1");
        TTreeReaderValue<Double_t> time0_charge2(rTime, "time0_charge2");
        TTreeReaderValue<Double_t> time1_charge2(rTime, "time1_charge2");
        TTreeReaderValue<Double_t> time0_utpc(rTime, "time0_utpc");
        TTreeReaderValue<Double_t> time1_utpc(rTime, "time1_utpc");
        TTreeReaderValue<Double_t> time0_algo(rTime, "time0_algo");
        TTreeReaderValue<Double_t> time1_algo(rTime, "time1_algo");         
*/
        Long64_t i = 0;
        while (rTime.Next()) {
        /*
        	if(strcmp(algo, "charge2") == 0) {
        		if(strcmp(choice, "time1") == 0) {
					index.emplace_back(static_cast<double>(*time1_charge2), i);
				}
				else if(strcmp(choice, "mean") == 0) {
					index.emplace_back(static_cast<double>(0.5*(*time0_charge2+ *time1_charge2)), i);
				}
				else {
					index.emplace_back(static_cast<double>(*time0_charge2), i);
				}
			}
        	else if(strcmp(algo, "algo") == 0) {
        	    if(strcmp(choice, "time1") == 0) {
					index.emplace_back(static_cast<double>(*time1_algo), i);
				}
				else if(strcmp(choice, "mean") == 0) {
					index.emplace_back(static_cast<double>(0.5*(*time0_algo+ *time1_algo)), i);
				}
				else {
					index.emplace_back(static_cast<double>(*time0_algo), i);
				}
			}
        	else if(strcmp(algo, "utpc") == 0) {
        	    if(strcmp(choice, "time1") == 0) {
					index.emplace_back(static_cast<double>(*time1_utpc), i);
				}
				else if(strcmp(choice, "mean") == 0) {
					index.emplace_back(static_cast<double>(0.5*(*time0_utpc+ *time1_utpc)), i);
				}
				else {
					index.emplace_back(static_cast<double>(*time0_utpc), i);
				}
        	}
        	else {
				if(strcmp(choice, "time1") == 0) {
					index.emplace_back(static_cast<double>(*time1), i);
				}
				else if(strcmp(choice, "mean") == 0) {
					index.emplace_back(static_cast<double>(0.5*(*time0+ *time1)), i);
				}
				else {
					index.emplace_back(static_cast<double>(*time0), i);
				}
        	}
        	*/
        	index.emplace_back(static_cast<double>(*time0), i);
            ++i;
            if ((i % 100000) == 0) printf("Indexing entry %lld / %lld\r", i, nentries);
        }
        printf("Index built with %lld elements\n", (Long64_t)index.size());
    }

    // --- Sort by time0 (stable: preserves original order for equal time0)
    std::stable_sort(index.begin(), index.end(),
                     [](const auto& a, const auto& b){ return a.first < b.first; });

    // --- Create output file & a fresh flat tree
    TFile* fout = TFile::Open(outfile, "RECREATE");
    if (!fout || fout->IsZombie()) { printf("Error: cannot create %s\n", outfile); return 1; }

    TTree* out = new TTree("clusters_detector", "clusters_detector (sorted by time0)");
    out->SetAutoFlush(30000000);   // ~30 MB
    out->SetAutoSave(300000000);   // ~300 MB

    // --- Output buffers (match your duplication script)
    unsigned int _id, _id0, _id1, _id2;
    UChar_t _det;
    unsigned short _size0, _size1, _size2;
    unsigned short _adc0, _adc1, _adc2;
    Double_t _pos0, _pos1, _pos2;
    Double_t _time0, _time1, _time2;
    Double_t _pos0_utpc, _pos1_utpc, _pos2_utpc;
    Double_t _time0_utpc, _time1_utpc, _time2_utpc;
    Double_t _pos0_charge2, _pos1_charge2, _pos2_charge2;
    Double_t _time0_charge2, _time1_charge2, _time2_charge2;
    Double_t _pos0_algo, _pos1_algo, _pos2_algo;
    Double_t _time0_algo, _time1_algo, _time2_algo;
    Double_t _dt0, _dt1, _dt2;
    Double_t _delta_plane_0_1, _delta_plane_1_2, _delta_plane_0_2;
    unsigned short _span_cluster0, _span_cluster1, _span_cluster2;
    unsigned short _max_delta_time0, _max_delta_time1, _max_delta_time2;
    unsigned short _max_missing_strip0, _max_missing_strip1, _max_missing_strip2;
    std::vector<double> _strips0, _times0, _adcs0;
    std::vector<double> _strips1, _times1, _adcs1;
    std::vector<double> _strips2, _times2, _adcs2;

    // --- Define branches
    out->Branch("id",  &_id);
    out->Branch("id0", &_id0);
    out->Branch("id1", &_id1);
    out->Branch("id2", &_id2);
    out->Branch("det", &_det);
    out->Branch("size0", &_size0);
    out->Branch("size1", &_size1);
    out->Branch("size2", &_size2);
    out->Branch("adc0", &_adc0);
    out->Branch("adc1", &_adc1);
    out->Branch("adc2", &_adc2);
    out->Branch("pos0", &_pos0);
    out->Branch("pos1", &_pos1);
    out->Branch("pos2", &_pos2);
    out->Branch("time0", &_time0);
    out->Branch("time1", &_time1);
    out->Branch("time2", &_time2);
    out->Branch("pos0_utpc", &_pos0_utpc);
    out->Branch("pos1_utpc", &_pos1_utpc);
    out->Branch("pos2_utpc", &_pos2_utpc);
    out->Branch("time0_utpc", &_time0_utpc);
    out->Branch("time1_utpc", &_time1_utpc);
    out->Branch("time2_utpc", &_time2_utpc);
    out->Branch("pos0_charge2", &_pos0_charge2);
    out->Branch("pos1_charge2", &_pos1_charge2);
    out->Branch("pos2_charge2", &_pos2_charge2);
    out->Branch("time0_charge2", &_time0_charge2);
    out->Branch("time1_charge2", &_time1_charge2);
    out->Branch("time2_charge2", &_time2_charge2);
    out->Branch("pos0_algo", &_pos0_algo);
    out->Branch("pos1_algo", &_pos1_algo);
    out->Branch("pos2_algo", &_pos2_algo);
    out->Branch("time0_algo", &_time0_algo);
    out->Branch("time1_algo", &_time1_algo);
    out->Branch("time2_algo", &_time2_algo);
    out->Branch("dt0", &_dt0);
    out->Branch("dt1", &_dt1);
    out->Branch("dt2", &_dt2);
    out->Branch("delta_plane_0_1", &_delta_plane_0_1);
    out->Branch("delta_plane_1_2", &_delta_plane_1_2);
    out->Branch("delta_plane_0_2", &_delta_plane_0_2);
    out->Branch("span_cluster0", &_span_cluster0);
    out->Branch("span_cluster1", &_span_cluster1);
    out->Branch("span_cluster2", &_span_cluster2);
    out->Branch("max_delta_time0", &_max_delta_time0);
    out->Branch("max_delta_time1", &_max_delta_time1);
    out->Branch("max_delta_time2", &_max_delta_time2);
    out->Branch("max_missing_strip0", &_max_missing_strip0);
    out->Branch("max_missing_strip1", &_max_missing_strip1);
    out->Branch("max_missing_strip2", &_max_missing_strip2);
    out->Branch("strips0", &_strips0);
    out->Branch("times0",  &_times0);
    out->Branch("adcs0",   &_adcs0);
    out->Branch("strips1", &_strips1);
    out->Branch("times1",  &_times1);
    out->Branch("adcs1",   &_adcs1);
    out->Branch("strips2", &_strips2);
    out->Branch("times2",  &_times2);
    out->Branch("adcs2",   &_adcs2);

    // --- Full reader for copying (all branches)
    TTreeReader fReader(tree);

    TTreeReaderValue<unsigned int> id(fReader, "id");
    TTreeReaderValue<unsigned int> id0(fReader, "id0");
    TTreeReaderValue<unsigned int> id1(fReader, "id1");
    TTreeReaderValue<unsigned int> id2(fReader, "id2");
    TTreeReaderValue<UChar_t> det(fReader, "det");
    TTreeReaderValue<unsigned short> size0(fReader, "size0");
    TTreeReaderValue<unsigned short> size1(fReader, "size1");
    TTreeReaderValue<unsigned short> size2(fReader, "size2");
    TTreeReaderValue<unsigned short> adc0(fReader, "adc0");
    TTreeReaderValue<unsigned short> adc1(fReader, "adc1");
    TTreeReaderValue<unsigned short> adc2(fReader, "adc2");
    TTreeReaderValue<Double_t> pos0(fReader, "pos0");
    TTreeReaderValue<Double_t> pos1(fReader, "pos1");
    TTreeReaderValue<Double_t> pos2(fReader, "pos2");
    TTreeReaderValue<Double_t> time0(fReader, "time0");
    TTreeReaderValue<Double_t> time1(fReader, "time1");
    TTreeReaderValue<Double_t> time2(fReader, "time2");
    TTreeReaderValue<Double_t> pos0_utpc(fReader, "pos0_utpc");
    TTreeReaderValue<Double_t> pos1_utpc(fReader, "pos1_utpc");
    TTreeReaderValue<Double_t> pos2_utpc(fReader, "pos2_utpc");
    TTreeReaderValue<Double_t> time0_utpc(fReader, "time0_utpc");
    TTreeReaderValue<Double_t> time1_utpc(fReader, "time1_utpc");
    TTreeReaderValue<Double_t> time2_utpc(fReader, "time2_utpc");
    TTreeReaderValue<Double_t> pos0_charge2(fReader, "pos0_charge2");
    TTreeReaderValue<Double_t> pos1_charge2(fReader, "pos1_charge2");
    TTreeReaderValue<Double_t> pos2_charge2(fReader, "pos2_charge2");
    TTreeReaderValue<Double_t> time0_charge2(fReader, "time0_charge2");
    TTreeReaderValue<Double_t> time1_charge2(fReader, "time1_charge2");
    TTreeReaderValue<Double_t> time2_charge2(fReader, "time2_charge2");
    TTreeReaderValue<Double_t> pos0_algo(fReader, "pos0_algo");
    TTreeReaderValue<Double_t> pos1_algo(fReader, "pos1_algo");
    TTreeReaderValue<Double_t> pos2_algo(fReader, "pos2_algo");
    TTreeReaderValue<Double_t> time0_algo(fReader, "time0_algo");
    TTreeReaderValue<Double_t> time1_algo(fReader, "time1_algo");
    TTreeReaderValue<Double_t> time2_algo(fReader, "time2_algo");
    TTreeReaderValue<Double_t> dt0(fReader, "dt0");
    TTreeReaderValue<Double_t> dt1(fReader, "dt1");
    TTreeReaderValue<Double_t> dt2(fReader, "dt2");
    TTreeReaderValue<Double_t> delta_plane_0_1(fReader, "delta_plane_0_1");
    TTreeReaderValue<Double_t> delta_plane_1_2(fReader, "delta_plane_1_2");
    TTreeReaderValue<Double_t> delta_plane_0_2(fReader, "delta_plane_0_2");
    TTreeReaderValue<unsigned short> span_cluster0(fReader, "span_cluster0");
    TTreeReaderValue<unsigned short> span_cluster1(fReader, "span_cluster1");
    TTreeReaderValue<unsigned short> span_cluster2(fReader, "span_cluster2");
    TTreeReaderValue<unsigned short> max_delta_time0(fReader, "max_delta_time0");
    TTreeReaderValue<unsigned short> max_delta_time1(fReader, "max_delta_time1");
    TTreeReaderValue<unsigned short> max_delta_time2(fReader, "max_delta_time2");
    TTreeReaderValue<unsigned short> max_missing_strip0(fReader, "max_missing_strip0");
    TTreeReaderValue<unsigned short> max_missing_strip1(fReader, "max_missing_strip1");
    TTreeReaderValue<unsigned short> max_missing_strip2(fReader, "max_missing_strip2");
    TTreeReaderArray<double> strips0(fReader, "strips0");
    TTreeReaderArray<double> times0(fReader, "times0");
    TTreeReaderArray<double> adcs0(fReader, "adcs0");
    TTreeReaderArray<double> strips1(fReader, "strips1");
    TTreeReaderArray<double> times1(fReader, "times1");
    TTreeReaderArray<double> adcs1(fReader, "adcs1");
    TTreeReaderArray<double> strips2(fReader, "strips2");
    TTreeReaderArray<double> times2(fReader, "times2");
    TTreeReaderArray<double> adcs2(fReader, "adcs2");

    auto copy_and_fill = [&]() {
        _id  = *id;   _id0 = *id0; _id1 = *id1; _id2 = *id2;
        _det = *det;

        _size0 = *size0; _size1 = *size1; _size2 = *size2;
        _adc0  = *adc0;  _adc1  = *adc1;  _adc2  = *adc2;

        _pos0 = *pos0; _pos1 = *pos1; _pos2 = *pos2;
        _time0 = *time0; _time1 = *time1; _time2 = *time2;

        _pos0_utpc = *pos0_utpc; _pos1_utpc = *pos1_utpc; _pos2_utpc = *pos2_utpc;
        _time0_utpc = *time0_utpc; _time1_utpc = *time1_utpc; _time2_utpc = *time2_utpc;

        _pos0_charge2 = *pos0_charge2; _pos1_charge2 = *pos1_charge2; _pos2_charge2 = *pos2_charge2;
        _time0_charge2 = *time0_charge2; _time1_charge2 = *time1_charge2; _time2_charge2 = *time2_charge2;

        _pos0_algo = *pos0_algo; _pos1_algo = *pos1_algo; _pos2_algo = *pos2_algo;
        _time0_algo = *time0_algo; _time1_algo = *time1_algo; _time2_algo = *time2_algo;

        _dt0 = *dt0; _dt1 = *dt1; _dt2 = *dt2;
        _delta_plane_0_1 = *delta_plane_0_1;
        _delta_plane_1_2 = *delta_plane_1_2;
        _delta_plane_0_2 = *delta_plane_0_2;

        _span_cluster0 = *span_cluster0;
        _span_cluster1 = *span_cluster1;
        _span_cluster2 = *span_cluster2;

        _max_delta_time0 = *max_delta_time0;
        _max_delta_time1 = *max_delta_time1;
        _max_delta_time2 = *max_delta_time2;

        _max_missing_strip0 = *max_missing_strip0;
        _max_missing_strip1 = *max_missing_strip1;
        _max_missing_strip2 = *max_missing_strip2;

        _strips0.assign(strips0.begin(), strips0.end());
        _times0.assign(times0.begin(), times0.end());
        _adcs0.assign(adcs0.begin(), adcs0.end());

        _strips1.assign(strips1.begin(), strips1.end());
        _times1.assign(times1.begin(), times1.end());
        _adcs1.assign(adcs1.begin(), adcs1.end());

        _strips2.assign(strips2.begin(), strips2.end());
        _times2.assign(times2.begin(), times2.end());
        _adcs2.assign(adcs2.begin(), adcs2.end());

        out->Fill();
    };

    // --- Pass 2: iterate in sorted order using SetEntry()
    Long64_t filled = 0;
    for (Long64_t i = 0; i < (Long64_t)index.size(); ++i) {
        const Long64_t entry = index[i].second;
        auto st = fReader.SetEntry(entry);
        if (st != TTreeReader::kEntryValid) {
            ::Error("sort_tree", "Invalid entry %lld (status=%d), skipping", entry, int(st));
            continue;
        }
        // optional: reserve to reduce reallocs
        _strips0.reserve(strips0.GetSize()); _times0.reserve(times0.GetSize()); _adcs0.reserve(adcs0.GetSize());
        _strips1.reserve(strips1.GetSize()); _times1.reserve(times1.GetSize()); _adcs1.reserve(adcs1.GetSize());
        _strips2.reserve(strips2.GetSize()); _times2.reserve(times2.GetSize()); _adcs2.reserve(adcs2.GetSize());

        copy_and_fill();

        ++filled;
        if(filled >= event_limit) {
        	break;
        }
    }
	printf("Filled %lld / %lld\r", filled, (Long64_t)index.size());
    // Save
    fout->cd();
    out->Write("", TObject::kOverwrite);
    fout->Close();
    fin->Close();

    printf("\nDone. Wrote %lld sorted entries to %s\n", filled, outfile);
    return 0;
}
