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
    int modules[7] = {7, 6, 5, 4, 3, 2, 1};
    gStyle->SetOptStat(11);  // Display just histo's name and # of entries
    gStyle->SetStatX(0.3);   // Stats box on the left

    for (int runNumber : fileNumbers) {
        // Construct the filename for the specific run
        TString filename = Form("Calo/AnaFOOT_Calo_Decoded_HIT2022_%d.root", runNumber);
        std::string runNumberStr = ConvertFileNumbersToString(runNumber);
        // Open the ROOT file
        TFile* inFile = new TFile(filename);
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        Double_t beamEnergy = RetrieveEnergy(runNumber, inFile);

        // Create two canvases: one for positions and one for crystal IDs
        TCanvas* c1 = new TCanvas("c1", Form("Position Display Run %d", runNumber), 2400, 1200);
        c1->Divide(6, 2, 0.01, 0.09);  // 2 rows and 6 columns for position histograms
        // Add the global title using TLatex
        c1->cd();
        TLatex* latex = new TLatex();
        latex->SetNDC();  // Use normalized device coordinates
        latex->SetTextAlign(22);  // Center the text horizontally
        latex->SetTextSize(0.02);  // Adjust text size to fit properly
        latex->DrawLatex(0.5, 0.95, Form("2D Histograms: x, y hits Run: %s | Beam Energy: %.0f MeV", runNumberStr.c_str(), beamEnergy));

        int customPadOrder[7] = {1, 2, 3, 4, 5, 6, 12};  // Place module 1 in the last pad (12)
        // Loop over all modules and display their respective histograms in subplots
        int index = 0;
        for (int moduleID : modules) {
            // --- For the position plot (hCalMapPos) ---
            int padNumber = customPadOrder[index];  // Get the custom pad number
            c1->cd(padNumber);
            gPad->SetRightMargin(0.15);  // Adjust margin for color palette space
            gPad->SetLeftMargin(0.1);
            gPad->SetTopMargin(0.1);
            gPad->SetBottomMargin(0.15);
            gStyle->SetPalette(1);

            // Retrieve the 2D position histogram for the current module
            TH2D* hCalMapPos = (TH2D*)inFile->Get(Form("hCalMapPos_module_%d", moduleID));
            if (!hCalMapPos) {
                std::cerr << "Error: Could not find hCalMapPos for module " << moduleID << " in file " << filename << std::endl;
                continue;
            }

            int entriesPos = hCalMapPos->GetEntries();
            std::cout << "Drawing position histogram for module " << moduleID << ", Entries: " << entriesPos << std::endl;

            if (entriesPos > 0) {
                hCalMapPos->Draw("COLZ");  // Draw histogram with color palette
                hCalMapPos->GetXaxis()->SetTitle("X");
                hCalMapPos->GetYaxis()->SetTitle("Y");
                //hCalMapPos->SetMarkerStyle(20);
                //hCalMapPos->SetMarkerSize(0.2);
                gPad->SetLogz(1);
            } else {
                std::cout << "Warning: Position histogram for module " << moduleID << " is empty." << std::endl;
            }
            //hCalMapPos->SetDirectory(0);

            // Retrieve the 2D crystal ID histogram for the current module
            TH2D* hCalMapCrystalID = (TH2D*)inFile->Get(Form("hCalMapCrystalID_module_%d", moduleID));
            if (!hCalMapCrystalID) {
                std::cerr << "Error: Could not find hCalMapCrystalID for module " << moduleID << " in file " << filename << std::endl;
                continue;
            }

            int entriesID = hCalMapCrystalID->GetEntries();

            std::cout << "Drawing crystalID histogram for module " << moduleID << ", Entries: " << entriesID << std::endl;
            hCalMapCrystalID->SetMinimum(0);  //the 0 value has to be displayed, it would be neglected instead
            gStyle->SetPaintTextFormat("1.0f");  // Format to display even small values like 0
            if (entriesID > 0) {
                hCalMapCrystalID->Draw("TEXT SAME");
                hCalMapCrystalID->GetXaxis()->SetTitle("X");
                hCalMapCrystalID->GetYaxis()->SetTitle("Y");
                hCalMapCrystalID->SetMarkerStyle(20);
                hCalMapCrystalID->SetMarkerSize(1.2);
            } else {
                std::cout << "Warning: CrystalID histogram for module " << moduleID << " is empty." << std::endl;
            }
            //hCalMapCrystalID->SetDirectory(0);


            index += 1;
        }

        // Update and save the canvases
        c1->cd();
        c1->Modified();
        c1->Update();
        c1->SaveAs(Form("Plots/CaloModules_Run_%d.png", runNumber));

        // Close the input file after processing
        inFile->Close();
        delete inFile;
    }
}

