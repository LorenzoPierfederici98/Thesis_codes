
#include "DisplayClusterScatter.h"

void DisplayClusterScatter() {
    gStyle->SetOptStat(0);  // Display histogram stats (name and entries)
    std::vector<std::pair<std::string, double>> filesAndEnergies = {
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_100MeV.root", 100},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_140MeV.root", 140},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_180MeV.root", 180},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_200MeV.root", 200},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_220MeV.root", 220}
    };

    for (const auto &[fileName, energy] : filesAndEnergies) {
        TFile *inFile = TFile::Open(fileName.c_str());
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error: Cannot open file " << fileName << std::endl;
            continue;
        }

        for (int crystal_ID = 0; crystal_ID < 63; crystal_ID++) {
            TString histName = Form("ClusterSize_Charge_crystalId_%d", crystal_ID);
            TH2D *h2 = dynamic_cast<TH2D *>(inFile->Get(histName));
            if (!h2) {
                std::cerr << "Error: Histogram " << histName << " not found in file " << fileName << std::endl;
                continue;
            }

            // Create a canvas to draw the histogram
            TCanvas *c = new TCanvas("c", "2D Histogram", 800, 600);
            c->SetRightMargin(0.15);  // Increase right margin to accommodate the color bar
            gPad->SetLogz(1);
            gStyle->SetPalette(1);
            h2->SetTitle(Form("Cluster Size vs. Charge Beam Energy: %.0f MeV | Crystal ID: %d", energy, crystal_ID));
            h2->GetXaxis()->SetTitle("Charge [a.u.]");
            h2->GetYaxis()->SetTitle("Cluster Size");
            h2->Draw("COLZ");  // Draw as a 2D colored histogram

            // Save the histogram as a PNG file
            c->SaveAs(Form("Plots/ClusterCharge/Merged_Histogram_Cluster_Crystal_%d_%.0fMeV.png", crystal_ID, energy));

            delete c;  // Clean up canvas
        }

        inFile->Close();
        delete inFile;  // Clean up input file
    }
}
