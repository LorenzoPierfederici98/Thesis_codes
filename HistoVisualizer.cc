#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include <TH1D.h>

TH1D* getHistDir(TFile* file, const TString &dirName, const TString &desiredHistName);

void HistoVisualizer(const std::string& Detector, const std::string& runName, const std::string& histoName, const std::string& histoTitle, const std::string& xLabel) {

    TString fileName = Form("TW/cuts/AnaFOOT_%s_Decoded_HIT2022_%s.root", Detector.c_str(), runName.c_str());
    //TString fileName = Form("MC/TW/AnaFOOT_%s_DecodedMC_HIT2022_MC_%s.root", Detector.c_str(), runName.c_str());
    TFile *inFile = TFile::Open(fileName.Data());

    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return;
    }

    //TH1D *h1 = dynamic_cast<TH1D *>(inFile->Get(histoName.c_str()));
    //TGraph *h1 = dynamic_cast<TGraph *>(inFile->Get(histoName.c_str()));
    //TCanvas *canvas = dynamic_cast<TCanvas *>(inFile->Get(histoName.c_str()));
    //TH1D *h1 = (TH1D*)canvas->FindObject("Charge_Calo_crystalId_1");

    TH1D *h1 = getHistDir(inFile, "ChargeTimeLayerY", histoName);
    if (!h1) {
        std::cerr << "Error: Histogram " << histoName << " not found in file " << fileName << std::endl;
        return;
    }
    // Create a canvas to draw the histogram
    TCanvas *c = new TCanvas("c", "1D Histogram", 800, 600);
    c->SetMargin(0.15, 0.12, 0.15, 0.15); // Left, Right, Bottom, Top margins

    gPad->SetLogy();
    gStyle->SetOptStat(1111);  // Display histogram stats
    gStyle->SetStatX(0.8);    // X position of the top-right corner
    gStyle->SetStatY(0.9);     // Y position of the top-right corner
    gStyle->SetPalette(1);
    h1->SetTitle(histoTitle.c_str());
    //h1->GetXaxis()->SetRangeUser(-0.2, 1.);
    //h1->GetXaxis()->SetRangeUser(0., 0.5);
    h1->GetXaxis()->SetTitle(xLabel.c_str());
    h1->GetYaxis()->SetTitle("Entries");
    h1->SetLineWidth(2);

    gStyle->SetTitleSize(0.07, "T");
    h1->GetXaxis()->SetTitleSize(0.05);  // X-axis title size
    h1->GetYaxis()->SetTitleSize(0.05);  // Y-axis title size
    //h1->SetMinimum(50.);
    h1->Draw();
    c->SaveAs(Form("../../thesis_images/%s_%s.png", histoName.c_str(), runName.c_str()));
    //c->SaveAs(Form("Plots/%s_%s.png", histoName.c_str(), runName.c_str()));
    c->WaitPrimitive();
    delete c;
    inFile->Close();
    delete inFile;
    return;
}

TH1D* getHistDir(TFile* file, const TString &dirName, const TString &desiredHistName) {
    // Change to the desired directory.
    TDirectory* dir = dynamic_cast<TDirectory*>(file->Get(dirName.Data()));
    if (!dir) {
        std::cerr << "Directory " << dirName << " not found in file " << file->GetName() << std::endl;
        file->Close();
        return nullptr;
    }
    
    // Iterate over all keys in the directory.
    TIter next(dir->GetListOfKeys());
    TKey *key;
    TH1D* foundHist = nullptr;
    while ((key = (TKey*) next())) {
        TObject *obj = key->ReadObj();
        TH1D* hist = dynamic_cast<TH1D*>(obj);
        if (!hist) {
            delete obj;
            continue;
        }
        TString hName = hist->GetName();
        // Check if this is the histogram we are looking for.
        if (hName == desiredHistName) {
            foundHist = dynamic_cast<TH1D*>(hist->Clone());
            delete hist;
            break;
        }
        delete hist;
    }
    return foundHist;
}