#define DrawWaveform_cxx
#include "DrawWaveform.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void DrawWaveform::Loop()
{
//   In a ROOT session, you can do:
//      root> .L DrawWaveform.C
//      root> DrawWaveform t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   double* ticks = new double[1024];
   double* adcs  = new double[1024];
   for(size_t i = 0; i < 1024; i++)
   {
      ticks[i] = 0.4*(double) i;
   }
   

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<10;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      TCanvas* canvas = new TCanvas("canvas", "Canvas", 1200, 900);

      for(size_t i = 0; i < 1024; i++)
      {
         adcs[i] = (double) waveform[i];
      }

      TGraph* graph = new TGraph(1024, ticks, adcs);
      
      // Graph styling
      graph->SetTitle("");
      graph->GetXaxis()->SetTitle("Time (ns)");
      graph->GetYaxis()->SetTitle("ADC");
      graph->GetXaxis()->CenterTitle();
      graph->GetYaxis()->CenterTitle();

      graph->GetXaxis()->SetRangeUser(0, 409.6);
      graph->GetYaxis()->SetRangeUser(2000, 4200);
      
      canvas->cd();
      graph->Draw("AL");
      canvas->SaveAs(Form("Waveform_%i_example.pdf", jentry));
   }
}
