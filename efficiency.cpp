#include <TROOT.h>
#include <TGraph.h>
#include <TCanvas.h>

#include <iostream>

#include <Access.h>
#include <CFDHitFinderAlg.h>

int main(int argc, char* argv[]) {
    std::string npeFile = argv[1];
    std::string fileName = Form("/Users/dphan/Research/novatb/tof-analysis/ToF_Simulation_%s_2d5GSs.root", npeFile.c_str());
    Access* dataAccess = new Access(fileName.c_str());
//    dataAccess->Test();

    std::map<beamlinereco::CFDParams, double> paramSet;
    paramSet[beamlinereco::kADCNBits]                                   = 12;
    paramSet[beamlinereco::kADCDynamicRange]                            = 1;
    paramSet[beamlinereco::kADCOffset]                                  = 0;
    paramSet[beamlinereco::kTimeSamplingInterval]                       = 0.2;
    paramSet[beamlinereco::kNSamplingPoints]                            = 1024;
    paramSet[beamlinereco::kIsWaveformNegativePolarity]                 = 1;
    paramSet[beamlinereco::kCFDThreshold]                               = 0.4;
    paramSet[beamlinereco::kRawHitFinderThresholdInNoiseSigma]          = 6;
    paramSet[beamlinereco::kShortRawHitIgnoringDurationInTicks]         = 10;
    paramSet[beamlinereco::kConsecutiveHitSeperationDurationInTicks]    = 5;
    paramSet[beamlinereco::kGSFilter]                                   = 0;
    paramSet[beamlinereco::kGSFilterWindow]                             = 17;
    paramSet[beamlinereco::kGSFilterDegree]                             = 3;
    paramSet[beamlinereco::kIntergratedWindowFixed]                     = 0;
    paramSet[beamlinereco::kIntergratedWindowLowerLimitIndex]           = 0;
    paramSet[beamlinereco::kIntergratedWindowUpperLimitIndex]           = 1024;

    beamlinereco::CFDHitFinder<double>* cfdHF = new beamlinereco::CFDHitFinder<double>();
    cfdHF->SetParams(paramSet);

    TCanvas* canvas = new TCanvas("canvas", "Canvas", 1200, 1200);
    TH1D* hTimeDiff = new TH1D("hTimeDiff", ";Start Time Reconstruction Error (ns);", 50, -4, 4);
    
    unsigned int sumTrue = 0;
    unsigned int sumTrueWithDarkCurrent = 0;
    unsigned int sumFound = 0;
    for (Long64_t iEvt = 0; iEvt < dataAccess->GetEntries(); iEvt++) {
//    for (Long64_t iEvt = 0; iEvt < 100; iEvt++) {
        unsigned int* waveform;
        if (dataAccess->GetWaveform(iEvt) != NULL) waveform = dataAccess->GetWaveform(iEvt);
        else continue;
        auto nhits = dataAccess->GetTriggerNHits(iEvt);
        auto allhits = dataAccess->GetTotalNHits(iEvt);
        sumTrue += nhits;
        sumTrueWithDarkCurrent += allhits;

        std::vector<uint16_t> wfvector;
        for (unsigned int iADC = 0; iADC < 1024; iADC++) {
            wfvector.push_back((uint16_t)waveform[iADC]);
        }
        cfdHF->SetWaveform(wfvector, 0, 0);
        cfdHF->Go();
        std::map<double, beamlinereco::hit_t<double> > ToFHitCollection = cfdHF->GetHitCollection();
        for (auto itr = ToFHitCollection.begin(); itr != ToFHitCollection.end(); itr++) {
            if (TMath::Abs((*itr).first * 0.2 - 50) <= 10) {
                hTimeDiff->Fill((*itr).first * 0.2 - 50);
            }
        }
        sumFound += ToFHitCollection.size();
    }

    std::cout << "Total number of total hits: " << sumTrueWithDarkCurrent << std::endl;
    std::cout << "Total number of trigger hits: " << sumTrue << std::endl;
    std::cout << "Total number of hits found by CFD: " << sumFound << std::endl;
    std::cout << "Efficiency by CFD: " << (sumFound * 100) / sumTrueWithDarkCurrent << "%" << std::endl;

    canvas->cd();
    hTimeDiff->SetTitle(Form("Trigger NPE = %s", npeFile.c_str()));
    hTimeDiff->Draw("HIST");
    canvas->SaveAs(Form("TimeDiff_%s_2d5GSs.pdf", npeFile.c_str()));

    return 0;
}
