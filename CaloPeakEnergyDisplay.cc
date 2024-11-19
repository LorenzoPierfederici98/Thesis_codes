//Macro that plots the charge given by the fit performed in AnalyzePeakCrystal.cc in a
//single crystal vs the beam energy. The fit results are stored in root files, named like
//Fit_Calo_Crystal_1_Energy_200MeV.root. To be run with root -l 'CaloPeakEnergyDisplay.cc({180, 200, 220}, 1)'
//{180, 200, 220} being the vector of energy values and 1 the crystalID.

#include <TFile.h>
#include <TF1.h>
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

pair<double, double> RoundMeasurement(double value, double uncertainty){
    //Get the order of magnitude of the uncertainty
    int significantFigures = (int)std::ceil(-std::log10(uncertainty)) + 1;

    //Calculate the rounding factor based on significant figures
    double roundingFactor = std::pow(10, significantFigures);

    //Round the uncertainty and value accordingly
    double roundedUncertainty = std::round(uncertainty * roundingFactor) / roundingFactor;
    double roundedValue = std::round(value * roundingFactor) / roundingFactor;

    return {roundedValue, roundedUncertainty};
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
    graph->GetXaxis()->SetTitle("Beam Energy [MeV]");
    graph->GetYaxis()->SetTitle("Mean Charge [a.u]");
    double x_max = graph->GetXaxis()->GetXmax();

    TF1 *f1 = new TF1("f1", "pol1", 0., x_max);
    TFitResultPtr fitresult = graph->Fit(f1, "S0");

    auto [intercept, sigma_intercept] = RoundMeasurement(fitresult->Parameter(0), fitresult->ParError(0));
    auto [slope, sigma_slope] = RoundMeasurement(fitresult->Parameter(1), fitresult->ParError(1));
    double chi2 = fitresult->Chi2();
    int ndf = fitresult->Ndf();

    TPaveText *fitInfo = new TPaveText(0.15, 0.7, 0.45, 0.85, "NDC"); // coordinates in NDC
    fitInfo->SetFillColor(0);
    fitInfo->AddText(Form("intercept [a.u] = %f#pm %f", intercept, sigma_intercept));
    fitInfo->AddText(Form("slope [a.u / MeV] = %f#pm %f", slope, sigma_slope));
    fitInfo->AddText(Form("#chi^{2} / ndf = %.2f / %d", chi2, ndf));
    fitInfo->SetTextSize(0.02);

    f1->SetParameter(0, intercept);
    f1->SetParameter(1, slope);
    // Draw the graph after the loop
    graph->SetMinimum(0.0);
    graph->GetXaxis()->SetLimits(0., x_max);
    graph->Draw("AP");
    f1->Draw("same");
    fitInfo->Draw();
    gPad->Modified();
    gPad->Update();
    
    // Save the canvas as a PNG image
    canvas->SaveAs(Form("Plots/Calo_FitCharge_Energy_Crystal_%d.png", crystalID));

    delete graph;
    delete canvas;
}
