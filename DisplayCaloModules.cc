//Macro that plots the 2D histograms of x-y hit coordinates in the calorimeter.
//Histograms are retireved form the AnaLizeFOOT.cc root files output.
//To be run with root -l 'DisplayCaloModules.cc({4742, 4743, 4744, 4745, 4828})'

#include <TFile.h>      // For opening ROOT files
#include <TCanvas.h>    // For creating and handling canvases
#include <TH2D.h>       // For handling 2D histograms
#include <TLegend.h>    // For legends (optional, only if needed)
#include <iostream>     // For input/output operations (std::cerr)

Double_t RetrieveEnergy(int runNumber, TFile* inFile) {
    if (!inFile) {
        std::cerr << "Error: Input file is null for run " << runNumber << std::endl;
        return -1;  // Return a default invalid value
    }

    TObject* obj = inFile->Get("BeamEnergyInfo");
    
    if (!obj) {
        std::cerr << "Error: Could not retrieve 'BeamEnergyInfo' from file for run " << runNumber << std::endl;
        return -1;  // Return a default invalid value
    }

    // Check if the retrieved object is a TObjString
    if (obj->InheritsFrom(TObjString::Class())) {
        TObjString* energyObj = (TObjString*)obj;
        Double_t beamEnergy = energyObj->GetString().Atof();
        return beamEnergy;
    } else {
        std::cerr << "Error: 'BeamEnergyInfo' is not a TObjString for run " << runNumber << std::endl;
        return -1;  // Return a default invalid value
    }
}

std::string ConvertFileNumbersToString(const int& runNumber) {
    std::stringstream ss;
    ss << runNumber;
    return ss.str();
}


void DisplayCaloModules(const vector<int> &fileNumbers) {
    int modules[7] = {7, 6, 5, 4, 3, 2, 1};  // List of modules
    gStyle->SetOptStat(0);  // Display histogram stats (name and entries)
    
    TCanvas* c1 = new TCanvas("c1", "Combined Display", 1200, 800);
    c1->SetRightMargin(0.15);  // Increase right margin to accommodate the color bar
    gPad->SetLogz(1);
    gStyle->SetPalette(1);

    for (int runNumber : fileNumbers) {
        TString filename = Form("Calo/AnaFOOT_Calo_Decoded_HIT2022_%d.root", runNumber);
        TFile* inFile = new TFile(filename);
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        Double_t beamEnergy = RetrieveEnergy(runNumber, inFile);

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
                std::cerr << "Warning: hCalMapPos for module " << moduleID << " not found in file " << filename << std::endl;
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
                std::cerr << "Warning: hCalMapCrystalID for module " << moduleID << " not found in file " << filename << std::endl;
            }

        }

        // Use TLatex to add a visible title at the top of the canvas
        TLatex* title = new TLatex();
        title->SetNDC();  // Use normalized device coordinates (from 0 to 1)
        title->SetTextAlign(22);  // Center-align the title
        title->SetTextSize(0.03);  // Adjust text size
        title->DrawLatex(0.5, 0.93, Form("Combined 2D Histograms: x, y Hits and Crystal IDs for Run %d | Beam Energy: %.0f MeV", runNumber, beamEnergy));
        c1->Modified();
        c1->Update();
        c1->SaveAs(Form("Plots/CombinedCaloModules_Run_%d.png", runNumber));

        inFile->Close();
        delete inFile;
    }
}

