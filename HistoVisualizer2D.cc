#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include <TH2D.h>

void HistoVisualizer2D(const std::string& Detector, const std::string& runName, const std::string& histoName, const std::string& histoTitle, const std::string& xLabel, const std::string& yLabel) {

    TString fileName = Form("%s/AnaFOOT_%s_Decoded_HIT2022_%s.root", Detector.c_str(), Detector.c_str(), runName.c_str());
    TFile *inFile = TFile::Open(fileName.Data());

    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return;
    }

    TH2D *h2 = dynamic_cast<TH2D *>(inFile->Get(histoName.c_str()));
    if (!h2) {
        std::cerr << "Error: Histogram " << histoName << " not found in file " << fileName << std::endl;
        return;
    }

    // Create a canvas to draw the histogram
    TCanvas *c = new TCanvas("c", "2D Histogram", 800, 600);
    c->SetRightMargin(0.15);  // Increase right margin to accommodate the color bar
    gStyle->SetOptStat(1111);  // Display histogram stats (name and entries)
    gStyle->SetStatX(0.8);    // X position of the top-right corner
    gStyle->SetStatY(0.9);     // Y position of the top-right corner
    gPad->SetLogz(1);
    gStyle->SetPalette(1);
    h2->SetTitle(histoTitle.c_str());
    //h2->GetXaxis()->SetRangeUser(6., 14.);
    h2->GetXaxis()->SetTitle(xLabel.c_str());
    h2->GetYaxis()->SetTitle(yLabel.c_str());
    h2->Draw("COLZ");
    c->SaveAs(Form("../../presentazione/%s.png", histoName.c_str()));
    c->WaitPrimitive();
    inFile->Close();
    delete inFile;
    delete c;
    return;


}