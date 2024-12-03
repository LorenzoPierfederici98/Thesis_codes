// Macro that plots the 2D histograms of x-y hit coordinates in the calorimeter.
// Histograms are retireved form the AnaLizeFOOT.cc root merged files output.
// To be run with root -l -b -q 'DisplayCaloModules.cc()'

#include "DisplayCaloModules.h"

void DisplayCaloModules() {
    int modules[7] = {7, 6, 5, 4, 3, 2, 1};  // List of modules
    gStyle->SetOptStat(0);  // Display histogram stats (name and entries)

    std::vector<std::pair<std::string, double>> filesAndEnergies = {
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_100MeV.root", 100},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_140MeV.root", 140},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_180MeV.root", 180},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_200MeV.root", 200},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_220MeV.root", 220}
    };
    
    TCanvas* c1 = new TCanvas("c1", "Combined Display", 1200, 800);
    c1->SetRightMargin(0.15);  // Increase right margin to accommodate the color bar
    gPad->SetLogz(1);
    gStyle->SetPalette(1);

    for (const auto &[fileName, beamEnergy] : filesAndEnergies) {
        TFile* inFile = new TFile(fileName.c_str());
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error: Could not open file " << fileName << std::endl;
            return;
        }

        bool firstPlotPos = true;

        for (int moduleID : modules) {
            // Retrieve and draw `hCalMapPos` histogram (x-y hits)
            TH2D* hCalMapPos = (TH2D*)inFile->Get(Form("hCalMapPos_module_%d", moduleID));
            if (hCalMapPos) {
                hCalMapPos->GetXaxis()->SetTitle("X [cm]");
                hCalMapPos->GetYaxis()->SetTitle("Y [cm]");
                hCalMapPos->SetTitle("");

                if (firstPlotPos) {
                    //hCalMapPos->SetTitle(Form("Combined 2D Histograms: x, y Hits and Crystal IDs for Run %d | Beam Energy: %.0f MeV", runNumber, beamEnergy));
                    hCalMapPos->Draw("COLZ");  // Draw first position histogram with axes
                    firstPlotPos = false;
                } else {
                    hCalMapPos->Draw("COLZ SAME");  // Overlay subsequent histograms
                }
            } else {
                std::cerr << "Warning: hCalMapPos for module " << moduleID << " not found in file " << fileName << std::endl;
            }

            // Retrieve and draw `hCalMapCrystalID` histogram (Crystal IDs)
            TH2D* hCalMapCrystalID = (TH2D*)inFile->Get(Form("hCalMapCrystalID_module_%d", moduleID));
            if (hCalMapCrystalID) {
                hCalMapCrystalID->SetMinimum(0);  //the 0 value has to be displayed, it would be neglected instead
                gStyle->SetPaintTextFormat("1.0f");  // Format to display even small values like 0
                hCalMapCrystalID->SetTitle("");
                hCalMapCrystalID->SetMarkerStyle(20);
                hCalMapCrystalID->SetMarkerSize(1.2);  // Adjust marker size for clarity

                hCalMapCrystalID->Draw("TEXT SAME");  // Directly overlay the text histogram
            } else {
                std::cerr << "Warning: hCalMapCrystalID for module " << moduleID << " not found in file " << fileName << std::endl;
            }

        }


        // Use TLatex to add a visible title at the top of the canvas
        TLatex* title = new TLatex();
        title->SetNDC();  // Use normalized device coordinates (from 0 to 1)
        title->SetTextAlign(22);  // Center-align the title
        title->SetTextSize(0.03);  // Adjust text size
        title->DrawLatex(0.5, 0.93, Form("Merged 2D Histograms: x, y Hits and Crystal IDs Beam Energy: %.0f MeV", beamEnergy));
        c1->Modified();
        c1->Update();
        c1->SaveAs(Form("Plots/Merged_CaloModules_Energy_%.0f.png", beamEnergy));

        inFile->Close();
        delete inFile;
    }
}

