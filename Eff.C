void Eff() {
   double npe[] = {2, 5, 10, 20, 30, 40, 50};
   double eff_5GSPS[]   = {91, 99.3, 100, 100, 100, 100, 100};
   double eff_2d5GSPS[] = {70, 97, 100, 100, 100, 100, 100};
   TCanvas* canvas = new TCanvas("canvas", "Canvas", 1200, 1200);
   TGraph* graph_5GSPS = new TGraph(7, npe, eff_5GSPS);
   TGraph* graph_2d5GSPS = new TGraph(7, npe, eff_2d5GSPS);
   
   // Graph styling
   graph_5GSPS->SetTitle("");
   graph_5GSPS->GetXaxis()->SetTitle("NPE");
   graph_5GSPS->GetYaxis()->SetTitle("Efficiency (%)");
   graph_5GSPS->GetXaxis()->CenterTitle();
   graph_5GSPS->GetYaxis()->CenterTitle();

   graph_5GSPS->SetLineWidth(8);
   graph_5GSPS->SetLineStyle(kSolid);
   graph_5GSPS->SetLineColor(kRed-3);
   graph_5GSPS->SetMarkerStyle(21);
   graph_5GSPS->SetMarkerSize(1);
   graph_5GSPS->SetMarkerColor(kBlack);
   graph_2d5GSPS->SetLineWidth(4);
   graph_2d5GSPS->SetLineStyle(kSolid);
   graph_2d5GSPS->SetLineColor(kBlue-3);
   graph_2d5GSPS->SetMarkerStyle(22);
   graph_2d5GSPS->SetMarkerSize(1);
   graph_2d5GSPS->SetMarkerColor(kBlack);

   gStyle->SetTitleAlign(22);
   gStyle->SetTitleX(.5);
   gStyle->SetTitleY(.95);
   gStyle->SetTitleBorderSize(0);

   gStyle->SetTitleSize(.07, "xyz");
   gStyle->SetTitleOffset(.8, "xyz");
   gStyle->SetTitleOffset(.9, "y");
   gStyle->SetTitleSize(.07, "");
   gStyle->SetTitleOffset(.8, "");
   gStyle->SetLabelSize(.04, "xyz");
   gStyle->SetLabelOffset(.005, "xyz");

   graph_5GSPS->GetYaxis()->SetRangeUser(65, 105);
   
   canvas->cd();
   graph_5GSPS->Draw("APL");
   graph_2d5GSPS->Draw("PL SAME");
   canvas->Modified();
   canvas->Update();

   Double_t x1 = canvas->GetUxmin();
   Double_t y1 = 100;
   Double_t x2 = canvas->GetUxmax();
   Double_t y2 = 100;
   
   TLine* line = new TLine(x1, y1, x2, y2);
   line->SetLineWidth(2);
   line->SetLineColor(kOrange+2);
   line->SetLineStyle(kDashed);
   line->Draw();

   TLegend* legend = new TLegend(0.5, 0.3, 0.7, 0.5);
   legend->SetBorderSize(0);
   legend->AddEntry(graph_5GSPS,"5 GSPS","l");
   legend->AddEntry(graph_2d5GSPS,"2.5 GSPS","l");
   legend->Draw();

   canvas->SaveAs("Efficiency.pdf");
   
   
}