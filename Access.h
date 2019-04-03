#pragma once

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

class Access {
public:
    Access(std::string fileName, TTree *tree = 0);

    virtual ~Access();

    virtual unsigned int  GetTotalNHits(Long64_t entry);
    virtual unsigned int  GetTriggerNHits(Long64_t entry);
    virtual unsigned int* GetWaveform(Long64_t entry);
    virtual unsigned int* GetNPEOfHits(Long64_t entry);
    virtual unsigned int* GetTimeOfHitsInNS(Long64_t entry);

    void Test();

    virtual Long64_t GetEntries();

    TTree *fChain;
    Int_t fCurrent;

    unsigned int waveform[1024];
    unsigned int totalNHits;
    unsigned int triggerNHits;
    unsigned int npeOfHits[100];
    unsigned int timeOfHitsInNS[100];

    TBranch *b_waveform;
    TBranch *b_totalNHits;
    TBranch *b_triggerNHits;
    TBranch *b_npeOfHits;
    TBranch *b_timeOfHitsInNS;

    virtual Int_t GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void Init(TTree *tree);
};
