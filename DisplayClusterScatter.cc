
#include "DisplayClusterScatter.h"

#include <TFile.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <string>

void DisplayClusterScatter() {
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

            TCanvas *c = new TCanvas("c", "Scatter Plot from 2D Histogram", 800, 600);
            TGraph *scatterPlot = new TGraph();

            int pointIndex = 0;
            // Loop over all bins in the 2D histogram
            for (int binX = 1; binX <= h2->GetNbinsX(); ++binX) {
                for (int binY = 1; binY <= h2->GetNbinsY(); ++binY) {
                    double binContent = h2->GetBinContent(binX, binY);
                    if (binContent > 0) { // Include only non-empty bins
                        double x = h2->GetXaxis()->GetBinCenter(binX);
                        double y = h2->GetYaxis()->GetBinCenter(binY);

                        // Add a point for each occurrence of the bin content
                        for (int i = 0; i < static_cast<int>(binContent); ++i) {
                            scatterPlot->SetPoint(pointIndex++, x, y);
                        }
                    }
                }
            }

            // Customize and draw the scatter plot
            scatterPlot->SetMarkerStyle(20);  // Set marker style
            scatterPlot->SetMarkerSize(0.8);  // Set marker size
            scatterPlot->SetTitle(Form("Scatter Plot Cluster Size-Charge Beam Energy HE: %.0f MeV | Crystal ID: %d", energy, crystal_ID));
            scatterPlot->GetXaxis()->SetTitle("Charge [a.u.]");
            scatterPlot->GetYaxis()->SetTitle("Cluster Size");

            scatterPlot->Draw("AP");
            c->SaveAs(Form("Plots/ClusterCharge/Merged_ScatterPlot_Cluster_Crystal_%d_%.0fMeV.png", crystal_ID, energy));

            delete scatterPlot;
            delete c;
        }

        inFile->Close();
        delete inFile;
    }
}