#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include <TH1D.h>

void HistoVisualizer(const std::string& Detector, const std::string& runName, const std::string& histoName, const std::string& histoTitle, const std::string& xLabel, const std::string& yLabel) {

    TString fileName = Form("%s/AnaFOOT_%s_Decoded_HIT2022_%s.root", Detector.c_str(), Detector.c_str(), runName.c_str());
    TFile *inFile = TFile::Open(fileName.Data());

    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return;
    }

    TH1D *h1 = dynamic_cast<TH1D *>(inFile->Get(histoName.c_str()));
    if (!h1) {
        std::cerr << "Error: Histogram " << histoName << " not found in file " << fileName << std::endl;
        return;
    }

    // Create a canvas to draw the histogram
    TCanvas *c = new TCanvas("c", "1D Histogram", 800, 600);
    c->SetRightMargin(0.15);  // Increase right margin to accommodate the color bar
    gPad->SetLogy();
    gStyle->SetOptStat(1111);  // Display histogram stats
    gStyle->SetStatX(0.8);    // X position of the top-right corner
    gStyle->SetStatY(0.9);     // Y position of the top-right corner
    gStyle->SetPalette(1);
    h1->SetTitle(histoTitle.c_str());
    //h1->GetXaxis()->SetRangeUser(-0.2, 1.);
    h1->GetXaxis()->SetTitle(xLabel.c_str());
    h1->GetYaxis()->SetTitle(yLabel.c_str());
    h1->Draw();
    c->SaveAs(Form("../../presentazione/%s_%s.png", histoName.c_str(), runName.c_str()));
    c->WaitPrimitive();
    delete c;
    inFile->Close();
    delete inFile;
    return;


}