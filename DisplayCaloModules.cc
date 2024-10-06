//Macro that plots the 2D histograms of x-y hit coordinates in the calorimeter.
//To be run with root -l 'DisplayCaloModules.cc({4742, 4743, 4744, 4745, 4828})'

#include <TFile.h>      // For opening ROOT files
#include <TCanvas.h>    // For creating and handling canvases
#include <TH2D.h>       // For handling 2D histograms
#include <TLegend.h>    // For legends (optional, only if needed)
#include <iostream>     // For input/output operations (std::cerr)

void DisplayCaloModules(const vector<int> &fileNumbers) {
    int nModules = 7;  // Number of modules (adjust if needed)
    gStyle->SetOptStat(11);  // Display just histo's name and # of entries
    gStyle->SetStatX(0.3);  // Stats box on the left

    for(int runNumber : fileNumbers)
    {
        // Construct the filename for the specific run
        TString filename = Form("AnaFOOT_Merge_HIT2022_%d.root", runNumber);

        // Open the ROOT file
        TFile* inFile = new TFile(filename);
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        // Create a canvas and divide it into subplots
        TCanvas* c1 = new TCanvas("c1", Form("Display Run %d", runNumber), 2000, 800);
        c1->Divide(4, 2);  // 2 rows and 4 columns

        // Loop over all modules and display their respective histograms in subplots
        for (int moduleID = 0; moduleID < nModules; ++moduleID) {
            // Move to the appropriate pad (one pad per module)
            c1->cd(moduleID + 1);  // Go to the (moduleID + 1)-th pad

            // Retrieve the 2D histogram for the current module
            TH2D* hCalMapPos = (TH2D*)inFile->Get(Form("hCalMapPos_module_%d", moduleID));
            if (!hCalMapPos) {
                std::cerr << "Error: Could not find hCalMapPos for module " << moduleID << " in file " << filename << std::endl;
                continue;
            }
            // Check if the histogram has entries
            int entries = hCalMapPos->GetEntries();
            std::cout << "Drawing histogram for module " << moduleID << ", Entries: " << entries << std::endl;

            // Only draw if there are entries
            if (entries > 0) {
                hCalMapPos->Draw("COLZ TEXT");  // Draw histogram with color palette
                c1->SetRightMargin(0.5);
            } else {
                std::cout << "Warning: Histogram for module " << moduleID << " is empty." << std::endl;
            }
            hCalMapPos->SetDirectory(0);
        }

        c1->cd();
        c1->Modified();    // Ensure the canvas is marked as modified
        c1->Update();      // Ensure the entire canvas is updated

        c1->SaveAs(Form("Plots/CaloModules_Run_%d.png", runNumber));

        // Close the input file after processing
        inFile->Close();
        delete inFile;
    }
}
