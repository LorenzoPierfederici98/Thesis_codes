//Macro that plots the charge given by the fit performed in AnalyzePeakCrystal.cc in all
//the available crystals vs the beam energy. The fit results are stored in root files, named like
//Fit_Calo_Crystal_1_Energy_200MeV.root. To be run with root -l 'CaloPeakEnergyDisplay.cc({180, 200, 220})'
//{180, 200, 220} being the vector of energy values, it loops on every crystalID and extracts the fit
//values only on the available IDs.

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

void CaloPeakEnergyDisplay(const std::vector<int> &energies) {
    std::string energiesStr = ConvertFileNumbersToString(energies);

    for (int crystalID = 0; crystalID < 63; crystalID++) {
        TCanvas *canvas = new TCanvas("canvas", "Calo Fit Results vs Energies", 800, 600);

        TGraphErrors *graph = new TGraphErrors();
        int pointIndex = 0;

        for (int energy : energies) {
            TString filename = Form("FitCalo/Fit_Calo_Crystal_%d_Energy_%dMeV.root", crystalID, energy);
            TFile *inFile = TFile::Open(filename);

            if (!inFile || inFile->IsZombie()) {
                std::cerr << "Error: Could not open file " << filename << std::endl;
                if (inFile) inFile->Close(); // Ensure file is closed if it exists
                continue;
            }

            TFitResult *fitresult = (TFitResult *)inFile->Get("Fit_Results");
            if (!fitresult) {
                std::cerr << "Warning: Fit results not found in file " << filename << std::endl;
                inFile->Close();
                continue;
            }

            // Retrieve mean charge and error
            Double_t mean_charge = fitresult->Parameter(1);
            Double_t mean_charge_error = fitresult->ParError(1);

            std::cout << "Crystal " << crystalID << ", Energy " << energy 
                      << " MeV: Mean Charge = " << mean_charge 
                      << " ± " << mean_charge_error << std::endl;

            // Add point to the graph
            graph->SetPoint(pointIndex, energy, mean_charge);
            graph->SetPointError(pointIndex, 0, mean_charge_error);

            ++pointIndex;
            inFile->Close();
        }

        if (graph->GetN() == 0) {
            std::cerr << "No valid points for crystal " << crystalID << std::endl;
            delete graph;
            delete canvas;
            continue;
        }

        // Customize the graph
        graph->SetTitle(Form("Beam Energies HE: %s MeV | Crystal ID: %d", energiesStr.c_str(), crystalID));
        graph->SetMarkerStyle(24);
        graph->SetMarkerSize(1.2);
        graph->GetXaxis()->SetTitle("Beam Energy [MeV]");
        graph->GetYaxis()->SetTitle("Mean Charge [a.u.]");
        graph->SetMinimum(0.0);
        double x_max = graph->GetXaxis()->GetXmax();

        // Perform linear fit
        TF1 *f1 = new TF1("f1", "pol1", 0., graph->GetXaxis()->GetXmax());
        TFitResultPtr fitresult_linear = graph->Fit(f1, "S");

        if (!fitresult_linear.Get() || fitresult_linear->Status() != 0) {
            std::cerr << "Linear fit failed for crystal " << crystalID << std::endl;
        } else {
            auto [intercept, sigma_intercept] = RoundMeasurement(fitresult_linear->Parameter(0), fitresult_linear->ParError(0));
            auto [slope, sigma_slope] = RoundMeasurement(fitresult_linear->Parameter(1), fitresult_linear->ParError(1));
            double chi2 = fitresult_linear->Chi2();
            int ndf = fitresult_linear->Ndf();

            TPaveText *fitInfo = new TPaveText(0.3, 0.7, 0.45, 0.85, "NDC");
            fitInfo->SetFillColor(0);
            fitInfo->AddText(Form("Intercept [a.u.] = %f#pm %f", intercept, sigma_intercept));
            fitInfo->AddText(Form("Slope [a.u. / MeV] = %f#pm %f", slope, sigma_slope));
            fitInfo->AddText(Form("#chi^{2} / ndf = %.2f / %d", chi2, ndf));
            fitInfo->SetTextSize(0.03);

            graph->GetXaxis()->SetLimits(0., x_max);

            f1->SetParameter(0, intercept);
            f1->SetParameter(1, slope);

            // Draw fit info
            graph->Draw("AP");
            f1->Draw("same");
            canvas->cd();
            fitInfo->Draw("same");
            canvas->Modified();
            canvas->Update();
            canvas->SaveAs(Form("Plots/Calo_FitCharge_Energy_Crystal_%d.png", crystalID));

            delete fitInfo;
        }

        delete f1;
        delete graph;
        delete canvas;
    }
}

