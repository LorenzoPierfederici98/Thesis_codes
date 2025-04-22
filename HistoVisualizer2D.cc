#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include <TH2D.h>

TH2D* getHist2D(TFile* file, const std::string &desiredHistName);

void DrawBetheBloch();

void HistoVisualizer2D(const std::string& Detector, const std::string& runName, const std::string& histoName, const std::string& histoTitle, const std::string& xLabel, const std::string& yLabel) {

    //TString fileName = Form("%s/AnaFOOT_%s_Decoded_HIT2022_%s.root", Detector.c_str(), Detector.c_str(), runName.c_str());
    TString fileName = Form("Calo/intercalib/AnaFOOT_%s_Decoded_HIT2022_%s.root", Detector.c_str(), runName.c_str());
    TFile *inFile = TFile::Open(fileName.Data());

    TH2D *h2 = getHist2D(inFile, histoName);

    if (!h2) {
        std::cerr << "Error: Histogram " << histoName << " not found in file " << fileName << std::endl;
        return;
    }

    // Create a canvas to draw the histogram
    TCanvas *c = new TCanvas("c", "2D Histogram", 800, 600);
    c->SetMargin(0.15, 0.15, 0.15, 0.15); // Left, Right, Bottom, Top margins
    gStyle->SetOptStat(1111);  // Display histogram stats (name and entries)
    gStyle->SetTitleSize(0.07, "T");
    //gStyle->SetOptStat(0);  // Do not display histogram stats
    gStyle->SetStatX(0.8);    // X position of the top-right corner
    gStyle->SetStatY(0.9);     // Y position of the top-right corner
    gPad->SetLogz(1);
    gStyle->SetPalette(1);
    h2->SetMinimum(1.);
    //h2->GetXaxis()->SetRangeUser(6., 12.);
    //h2->GetYaxis()->SetRangeUser(0., 12.);
    h2->GetXaxis()->SetTitle(xLabel.c_str());
    h2->GetYaxis()->SetTitle(yLabel.c_str());
    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetTitleSize(0.05);
    h2->Draw("COLZ");
    //ex1->Draw();
    //TLegend *leg = new TLegend(0.5, 0.8, 0.85, 0.9);
    //leg->AddEntry(histMyEloss, "My eloss: Q*(1/p0)", "l");
    //leg->AddEntry(histEloss, "SHOE eloss: hit->GetEnergyLoss()", "l");
    //leg->AddEntry(histMC, "MC eloss", "l");
    //leg->AddEntry(h2_Z, "Z", "f");
    //leg->AddEntry(h2_Z1, "Z1", "f");
    //leg->AddEntry(h2_Z2, "Z2", "f");
    //leg->SetTextSize(0.025);
    //leg->Draw();
    h2->SetTitle(histoTitle.c_str());
    //DrawBetheBloch();
    c->SaveAs(Form("Plots/%s_200MeV.png", histoName.c_str()));
    //c->WaitPrimitive();
    inFile->Close();
    delete inFile;
    delete c;
    return;


}

TH2D* getHist2D(TFile* file, const std::string &desiredHistName) {
    if (!file || file->IsZombie()) {
        std::cerr << "Invalid or corrupted file!" << std::endl;
        return nullptr;
    }

    // Attempt to retrieve the histogram by name.
    TH2D* hist = dynamic_cast<TH2D*>(file->Get(desiredHistName.c_str()));
    if (!hist) {
        std::cerr << "Histogram " << desiredHistName << " not found in file " << file->GetName() << std::endl;
        return nullptr;
    }

    // Clone to ensure ownership outside the file.
    TH2D* clonedHist = dynamic_cast<TH2D*>(hist->Clone());
    if (!clonedHist) {
        std::cerr << "Failed to clone histogram " << desiredHistName << std::endl;
        return nullptr;
    }

    return clonedHist;
}

void DrawBetheBloch()
{
    TF1* Bethe_Bloch_Z1 = new TF1("Bethe_Bloch_Z1",
        "[2] * [6] * [3] * [4] * ([1]^2 / pow(([0]/(0.3*x)), 2)) * "
        "(TMath::Log(2 * 0.511 * pow(([0]/(0.3*x)), 2) * (1. / (1 - pow(([0]/(0.3*x)), 2))) / [5]) - pow(([0]/(0.3*x)), 2))",
        6., 12.);

    TF1* Bethe_Bloch_Z2 = new TF1("Bethe_Bloch_Z2",
        "[2] * [6] * [3] * [4] * ([1]^2 / pow(([0]/(0.3*x)), 2)) * "
        "(TMath::Log(2 * 0.511 * pow(([0]/(0.3*x)), 2) * (1. / (1 - pow(([0]/(0.3*x)), 2))) / [5]) - pow(([0]/(0.3*x)), 2))",
        6., 12.);
    
    // d_SC_TW, z, dx, rho, Z/A, I, K
    Bethe_Bloch_Z1->SetParameters(1.439, 1.0, 0.3, 1.023, 0.5417, 64.7E-6, 0.307);
    Bethe_Bloch_Z1->SetLineColor(kRed);
    Bethe_Bloch_Z1->SetLineWidth(2);
    Bethe_Bloch_Z1->Draw("same");
        
    Bethe_Bloch_Z2->SetParameters(1.439, 2.0, 0.3, 1.023, 0.5417, 64.7E-6, 0.307);
    Bethe_Bloch_Z2->SetLineWidth(2);
    Bethe_Bloch_Z1->SetLineColor(kBlue);
    Bethe_Bloch_Z2->Draw("same");

    TLegend *leg = new TLegend(0.2, 0.75, 0.45, 0.85);
    leg->AddEntry(Bethe_Bloch_Z1, "Bethe-Bloch Z1", "l");
    leg->AddEntry(Bethe_Bloch_Z2, "Bethe-Bloch Z2", "l");
    leg->SetTextSize(0.04);
    leg->Draw("same");
    return;
}