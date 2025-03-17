// Macro that collects the charge fit results of AnalyzePeakCrystal.cc for every calorimeter crystal
// from files named like AnaFOOT_Calo_Decoded_HIT2022_100MeV_Fit.root. A charge-beam energy plot is
// built for every crystal and fitted with a 1 parameter (slope) linear function. An intercalibration
// is then performed for those crystals with at least 2 fitted charge peaks out of 5 energy values, by
// computing the ratio of the slopes of a certain crystal with respect to crystal ID 0 (the central crystal).
// Those ratios are used as scaling factors to sum the charge for events with a single cluster of size 2,
// that is with 2 crystals involved. The calibration coefficients for every available crystal are stored
// in a file named SlopeRatios.cal, stored in calib/HIT2022/.

#include "CaloPeakEnergyDisplay.h"

void CaloPeakEnergyDisplay() {

    std::vector<std::pair<std::string, int>> filesAndEnergies = {
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_100MeV_Fit.root", 100},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_140MeV_Fit.root", 140},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_180MeV_Fit.root", 180},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_200MeV_Fit.root", 200},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_220MeV_Fit.root", 220}
    };

    std::vector<int> energies = {100, 140, 180, 200, 220};
    string energiesStr = ConvertFileNumbersToString(energies);
    std::map<int, std::pair<double, double>> slopeRatios;

    double reference_slope = 0.0;
    double reference_error = 0.0;

    for (int crystalID = 0; crystalID < 63; crystalID++) {
        TCanvas *canvas = new TCanvas("canvas", "Calo Fit Results vs Energies", 800, 600);
        TGraphErrors* graph = new TGraphErrors();
        int pointIndex = 0;
        for (const auto& [fileName, energy] : filesAndEnergies) {
            TFile* inFile = TFile::Open(fileName.c_str());
            if (!inFile || inFile->IsZombie()) {
                std::cerr << "Error: Could not open file " << fileName << std::endl;
                continue;
            }

            TFitResult* fitresult = (TFitResult*)inFile->Get(Form("Fit_result_Charge_Calo_crystalId_%d", crystalID));
            if (!fitresult) {
                std::cerr << "Warning: Fit results not found in file " << fileName << std::endl;
                inFile->Close();
                continue;
            }

            double mean_charge = fitresult->Parameter(1);
            double mean_charge_error = fitresult->ParError(1);

            graph->SetPoint(pointIndex, energy, mean_charge);
            graph->SetPointError(pointIndex, 0, mean_charge_error);

            ++pointIndex;
            inFile->Close();
        }

        if (graph->GetN() < 2) {
            std::cout << "Not enough points to fit: " << graph->GetN() << " crystal ID " << crystalID << std::endl;
            delete graph;
            continue;
        }

        TF1* f1 = new TF1("f1", "[0]*x", 0., graph->GetXaxis()->GetXmax());
        TFitResultPtr fitresult_linear = graph->Fit(f1, "S");
        if (fitresult_linear.Get() && fitresult_linear->Status() == 0) {

            //auto [intercept, sigma_intercept] = RoundMeasurement(fitresult_linear->Parameter(0), fitresult_linear->ParError(0));
            auto [slope, sigma_slope] = RoundMeasurement(fitresult_linear->Parameter(0), fitresult_linear->ParError(0));
            double chi2 = fitresult_linear->Chi2();
            int ndf = fitresult_linear->Ndf();

            graph->SetTitle(Form("Beam Energies HE: %s MeV/u | Crystal ID: %d", energiesStr.c_str(), crystalID));
            graph->SetMarkerStyle(24);
            graph->SetMarkerSize(1.2);
            graph->GetXaxis()->SetTitle("Beam Energy [MeV/u]");
            graph->GetYaxis()->SetTitle("Mean Charge [a.u.]");
            graph->SetMinimum(0.0);
            graph->GetXaxis()->SetLimits(0., graph->GetXaxis()->GetXmax());

            //f1->SetParameter(0, intercept);
            f1->SetParameter(0, slope);

            TPaveText *fitInfo = new TPaveText(0.3, 0.7, 0.45, 0.85, "NDC");
            fitInfo->SetFillColor(0);
            //fitInfo->AddText(Form("Intercept [a.u.] = %f#pm %f", intercept, sigma_intercept));
            fitInfo->AddText(Form("Slope [a.u. / MeV] = %f#pm %f", slope, sigma_slope));
            //fitInfo->AddText(Form("#chi^{2} / ndf = %.2f / %d", chi2, ndf));
            fitInfo->SetTextSize(0.03);

            if (crystalID == 0) {
                reference_slope = slope;
                reference_error = sigma_slope;
            } else {
                double ratio = slope / reference_slope;
                double ratio_error = ratio * sqrt(pow(sigma_slope / slope, 2) + pow(reference_error / reference_slope, 2));
                slopeRatios[crystalID] = {ratio, ratio_error};
            }

            // Draw fit info
            graph->Draw("AP");
            f1->Draw("same");
            canvas->cd();
            fitInfo->Draw("same");
            canvas->Modified();
            canvas->Update();
            canvas->SaveAs(Form("Plots/Merged_Calo_FitCharge_Energy_Crystal_%d.png", crystalID));

            delete fitInfo;
        }

        delete f1;
        delete graph;
        delete canvas;
    }

    std::ofstream outFile("SlopeRatios.cal");
    outFile << "#crystalID\t slope ratio(from ID 0)\t err" << std::endl;
    for (const auto& [crystalID, ratioData] : slopeRatios) {
        outFile << crystalID << "\t" << ratioData.first << "\t" << ratioData.second << std::endl;
    }
    outFile.close();
}


std::pair<double, double> RoundMeasurement(double value, double uncertainty) {
    int significantFigures = (int)std::ceil(-std::log10(uncertainty));
    double roundingFactor = std::pow(10, significantFigures);
    double roundedUncertainty = std::round(uncertainty * roundingFactor) / roundingFactor;
    double roundedValue = std::round(value * roundingFactor) / roundingFactor;
    return {roundedValue, roundedUncertainty};
}

std::string ConvertFileNumbersToString(const std::vector<int>& energies) {
    std::stringstream ss;
    for (size_t i = 0; i < energies.size(); ++i) {
        ss << energies[i];
        if (i != energies.size() - 1) {
            ss << ", ";
        }
    }
    return ss.str();
}