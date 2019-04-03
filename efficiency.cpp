#include <TROOT.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLine.h>

#include <iostream>
#include <stdlib.h>

#include <Access.h>
#include <CFDHitFinderAlg.h>

int main(int argc, char* argv[]) {
    std::string npeFile = argv[1];
    double samplingRateInGSPS = atof(argv[2]);
    std::string samplingRateStr = samplingRateInGSPS == 5 ? "5" : "2d5";
    std::string fileName = Form("/Users/dphan/Research/novatb/tof-analysis/ToF_Simulation_%s_%sGSs.root", npeFile.c_str(), samplingRateStr.c_str());
    Access* dataAccess = new Access(fileName.c_str());
//    dataAccess->Test();

    std::map<beamlinereco::CFDParams, double> paramSet;
    paramSet[beamlinereco::kADCNBits]                                   = 12;
    paramSet[beamlinereco::kADCDynamicRange]                            = 1;
    paramSet[beamlinereco::kADCOffset]                                  = 0;
    paramSet[beamlinereco::kTimeSamplingInterval]                       = 1 / (samplingRateInGSPS);
    paramSet[beamlinereco::kNSamplingPoints]                            = 1024;
    paramSet[beamlinereco::kIsWaveformNegativePolarity]                 = 1;
    paramSet[beamlinereco::kCFDThreshold]                               = 0.4;
    paramSet[beamlinereco::kRawHitFinderThresholdInNoiseSigma]          = 8;
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

    TH1D* hTimeDiff = new TH1D("hTimeDiff", ";Start Time Reconstruction Error (ns);", 100, -4, 4);
    
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

        double baseline = cfdHF->GetPedestal();
        std::map<double, beamlinereco::hit_t<double> > ToFHitCollection = cfdHF->GetHitCollection();

        for (auto & itr : ToFHitCollection) {
            if (TMath::Abs(itr.first * paramSet[beamlinereco::kTimeSamplingInterval] - 50) <= 10) {
                hTimeDiff->Fill(itr.first * paramSet[beamlinereco::kTimeSamplingInterval] - 50);
            }
        }
        sumFound += ToFHitCollection.size();
/*
        if (iEvt < 10) {
            double* ticks = new double[1024];
            double* adcs  = new double[1024];
            for(size_t ik = 0; ik < 1024; ik++) {
                ticks[ik] = paramSet[beamlinereco::kTimeSamplingInterval] * (double) ik;
                adcs[ik] = (double)wfvector.at(ik);
            }

            TCanvas* canvas = new TCanvas(Form("canvas_%i", iEvt), Form("Canvas_%i", iEvt), 1200, 900);
            TGraph* graph = new TGraph(1024, ticks, adcs);

            // Graph styling
            graph->SetTitle("");
            graph->GetXaxis()->SetTitle("Time (ns)");
            graph->GetYaxis()->SetTitle("ADC");
            graph->GetXaxis()->CenterTitle();
            graph->GetYaxis()->CenterTitle();

            graph->GetXaxis()->SetRangeUser(0, 1024 * paramSet[beamlinereco::kTimeSamplingInterval]);
            graph->GetYaxis()->SetRangeUser(2000, 4200);
            graph->SetLineColor(kBlack);

            graph->Draw("AL");

            Double_t x_min = 0;
            Double_t y_min = 2000;
            Double_t x_max = paramSet[beamlinereco::kTimeSamplingInterval] * paramSet[beamlinereco::kNSamplingPoints];
            Double_t y_max = 4200;

            std::cout << "baseline: " << baseline << "." << std::endl;
            
            TLine* line_pedestal = new TLine(x_min, baseline, x_max, baseline);
            line_pedestal->SetLineWidth(1);
            line_pedestal->SetLineColor(kBlue);
            line_pedestal->SetLineStyle(kDashed);
            line_pedestal->Draw();
            canvas->Modified(true);
            canvas->Update();

            std::vector<TLine*> startline;
            for (auto theIter = ToFHitCollection.begin(); theIter != ToFHitCollection.end(); theIter++) {
                unsigned int index = std::distance(ToFHitCollection.begin(), theIter);
                double startLineValue = paramSet[beamlinereco::kTimeSamplingInterval] * (*theIter).first;

                startline.push_back(new TLine(startLineValue, y_min, startLineValue, y_max));
                (startline.at(index))->SetLineColor(kRed - 3);

                (startline.at(index))->Draw();
                canvas->Modified(true);
                canvas->Update();
            }

            canvas->SaveAs(Form("/Users/dphan/Research/novatb/tof-analysis/Waveform_%sGSs/Waveform_%i_%s_%sGSs.pdf", samplingRateStr.c_str(), iEvt, npeFile.c_str(), samplingRateStr.c_str()), "Quiet");
        }
*/
    }

    std::cout << "Total number of total hits: " << sumTrueWithDarkCurrent << std::endl;
    std::cout << "Total number of trigger hits: " << sumTrue << std::endl;
    std::cout << "Total number of hits found by CFD: " << sumFound << std::endl;
    std::cout << "Efficiency by CFD: " << ((double) sumFound * 100) / (double)sumTrueWithDarkCurrent << "%" << std::endl;
/*
    TCanvas* canvas_main = new TCanvas("canvas_main", "Canvas", 1200, 1200);
    canvas_main->cd();
    hTimeDiff->SetTitle(Form("Trigger NPE = %s, Sampling Rate = %1.1fGSPS", npeFile.c_str(), samplingRateInGSPS));
    hTimeDiff->Draw("HIST");
    canvas_main->SaveAs(Form("/Users/dphan/Research/novatb/tof-analysis/TimeDiff/TimeDiff_%s_%1.1fGSs.pdf", npeFile.c_str(), samplingRateInGSPS), "Quiet");
*/
    return 0;
}
