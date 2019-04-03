//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  2 13:50:47 2019 by ROOT version 6.16/00
// from TTree TOF/TOF Waveforms
// found on file: ./ToF_Simulation_2_2d5GSs.root
//////////////////////////////////////////////////////////

#ifndef DrawWaveform_h
#define DrawWaveform_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DrawWaveform {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           totalNHits;
   Int_t           triggerNHits;
   Int_t           waveform[1024];
   Int_t           npeOfHits[3];   //[totalNHits]
   Int_t           timeOfHitsInNS[3];   //[totalNHits]

   // List of branches
   TBranch        *b_totalNHits;   //!
   TBranch        *b_triggerNHits;   //!
   TBranch        *b_waveform;   //!
   TBranch        *b_npeOfHits;   //!
   TBranch        *b_timeOfHitsInNS;   //!

   DrawWaveform(TTree *tree=0);
   virtual ~DrawWaveform();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DrawWaveform_cxx
DrawWaveform::DrawWaveform(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./ToF_Simulation_50_2d5GSs.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("./ToF_Simulation_50_2d5GSs.root");
      }
      f->GetObject("TOF",tree);

   }
   Init(tree);
}

DrawWaveform::~DrawWaveform()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DrawWaveform::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DrawWaveform::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DrawWaveform::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("totalNHits", &totalNHits, &b_totalNHits);
   fChain->SetBranchAddress("triggerNHits", &triggerNHits, &b_triggerNHits);
   fChain->SetBranchAddress("waveform", waveform, &b_waveform);
   fChain->SetBranchAddress("npeOfHits", npeOfHits, &b_npeOfHits);
   fChain->SetBranchAddress("timeOfHitsInNS", timeOfHitsInNS, &b_timeOfHitsInNS);
   Notify();
}

Bool_t DrawWaveform::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DrawWaveform::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DrawWaveform::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DrawWaveform_cxx
