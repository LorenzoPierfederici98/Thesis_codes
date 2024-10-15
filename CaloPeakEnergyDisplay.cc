#include <TFile.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <iostream>
#include <sstream>

std::string ConvertFileNumbersToString(const std::vector<int>& energies) {
    std::stringstream ss;
    for (size_t i = 0; i < energies.size(); ++i) {
        ss << energies[i];
        if (i != energies.size() - 1) {
            ss << ", ";  // Add a comma and space between numbers
        }
    }
    return ss.str();
}

void CaloPeakEnergyDisplay(const std::vector<int> &energies, const int crystalID) {
    TCanvas *canvas = new TCanvas("canvas", "Calo Fit Results vs Energies", 800, 600);
    std::string energiesStr = ConvertFileNumbersToString(energies);
    canvas->SetTitle(Form("Beam Energies HE: %s MeV | Crystal ID: %d", energiesStr.c_str(), crystalID));

    TGraphErrors *graph = new TGraphErrors();  // Create a single graph to accumulate points

    int pointIndex = 0;
    for(int energy : energies) {
        TString filename = Form("Fit_Calo_Crystal_%d_Energy_%dMeV.root", crystalID, energy);
        TFile *inFile = TFile::Open(filename);

        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            continue;
        }

        TFitResult *fitresult = (TFitResult *)inFile->Get("Fit_Results");

        if (!fitresult) {
            std::cerr << Form("Fit result for energy %d not present", energy) << std::endl;
            inFile->Close();
            continue;
        }

        // Retrieve mean charge and error from the fit result
        Double_t mean_charge = fitresult->Parameter(1);
        Double_t mean_charge_error = fitresult->ParError(1);

        // Set the points in the graph
        graph->SetPoint(pointIndex, energy, mean_charge);
        graph->SetPointError(pointIndex, 0, mean_charge_error);

        ++pointIndex;
        inFile->Close();
    }

    // Customize the graph
    graph->SetTitle(Form("Beam Energies HE: %s MeV | Crystal ID: %d", energiesStr.c_str(), crystalID));
    graph->SetMarkerStyle(24);
    graph->SetMarkerSize(1.2);
    graph->GetXaxis()->SetTitle("Primary beam energy [MeV]");
    graph->GetYaxis()->SetTitle("Mean charge");

    // Draw the graph after the loop
    graph->Draw("AP");
    
    // Save the canvas as a PNG image
    canvas->SaveAs(Form("Plots/Calo_FitCharge_Energy_Crystal_%d.png", crystalID));

    delete graph;
    delete canvas;
}
