#include <Access.h>
#include <iostream>

Access::Access(std::string fileName, TTree *tree) : fChain(0) {
    if (tree == 0) {
        TFile *f = (TFile *) gROOT->GetListOfFiles()->FindObject(fileName.c_str());
        if (!f || !f->IsOpen()) {
            f = new TFile(fileName.c_str());
        }
        f->GetObject("TOF", tree);
    }
    Init(tree);
}

Long64_t Access::GetEntries() {
    return (Long64_t) fChain->GetEntries();
}

unsigned int * Access::GetWaveform(Long64_t entry) {
    LoadTree(entry);
    GetEntry(entry);

    return &waveform[0];
}

unsigned int *Access::GetNPEOfHits(Long64_t entry) {
    LoadTree(entry);
    GetEntry(entry);
    return &npeOfHits[0];
}

unsigned int *Access::GetTimeOfHitsInNS(Long64_t entry) {
    LoadTree(entry);
    GetEntry(entry);
    return &timeOfHitsInNS[0];
}

unsigned int Access::GetTriggerNHits(Long64_t entry) {
    LoadTree(entry);
    GetEntry(entry);
    return triggerNHits;
}

unsigned int Access::GetTotalNHits(Long64_t entry) {
    LoadTree(entry);
    GetEntry(entry);
    return totalNHits;
}

Access::~Access() {
    if (!fChain)
        return;
    delete fChain->GetCurrentFile();
}

Int_t Access::GetEntry(Long64_t entry) {
    if (!fChain)
        return 0;
    return fChain->GetEntry(entry);
}

Long64_t Access::LoadTree(Long64_t entry) {
    // Set the environment to read one entry
    if (!fChain)
        return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0)
        return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
    }
    return centry;
}

void Access::Init(TTree *tree) {
    if (!tree)
        return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("waveform", waveform, &b_waveform);
    fChain->SetBranchAddress("totalNHits", &totalNHits, &b_totalNHits);
    fChain->SetBranchAddress("triggerNHits", &triggerNHits, &b_triggerNHits);
    fChain->SetBranchAddress("npeOfHits", npeOfHits, &b_npeOfHits);
    fChain->SetBranchAddress("timeOfHitsInNS", timeOfHitsInNS, &b_timeOfHitsInNS);
}

void Access::Test() {
    fChain->GetEntry(0);
    for (unsigned int i = 0; i < 1024; i++) {
        std::cout << "waveform[" << i << "]: " << waveform[i] << "." << std::endl;
    }
}