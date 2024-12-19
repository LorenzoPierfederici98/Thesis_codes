// Macro thath fits the ToF-raw distribution from the AnalyzeTWChargeTime.cc merged output files.
// The ToF distribution from bars 9 of the two layers are fitted, then the ToF mean and standard 
// deviation are plotted vs the beta value.

#include "AnalyzeTofFit.h"

double HELIUM_MASS = 3755.674;  // MeV/c^2 2m_p + 2m_n

void AnalyzeTofFit() {

    std::vector<std::pair<std::string, double>> filesAndEnergies = {
        {"TW/AnaFOOT_TW_Decoded_HIT2022_100MeV.root", 100},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_140MeV.root", 140},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_180MeV.root", 180},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_200MeV.root", 200},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_220MeV.root", 220}
    };

    std::map<TString, std::vector<std::pair<double, FitResult>>> results;

    for (const auto& [fileName, energy] : filesAndEnergies) {
        ProcessAndPlot(fileName, energy, results);
    }

    PlotResults(results);
}

FitResult SaveFitAndExtractParams(TFile* file, TH1* hist, TF1* fitFunc) {
    file->cd();  // Ensure we are writing to the correct file directory
    TString histName = TString(hist->GetName());
    TString fitName = histName + "_fit";
    hist->Write(fitName, TObject::kOverwrite);

    FitResult result;
    result.mean = fitFunc->GetParameter(1);  // Mean
    result.meanError = fitFunc->GetParError(1);
    result.sigma = fitFunc->GetParameter(2);  // Sigma
    result.sigmaError = fitFunc->GetParError(2);
    return result;
}

void ProcessAndPlot(const std::string& fileName, double energy,
                    std::map<TString, std::vector<std::pair<double, FitResult>>>& results) {
    TFile* inputFile = TFile::Open(fileName.c_str(), "READ");
    TString outputFileName = fileName;
    outputFileName.ReplaceAll(".root", "_Fit.root");
    TFile* outputFile = TFile::Open(outputFileName, "UPDATE");

    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return;
    }

    TString histoNames[] = {"ToF_Bar9_XY"};
    for (const auto& histo : histoNames) {
        TH1* hist = dynamic_cast<TH1*>(inputFile->Get(histo));
        if (hist) {
            TF1* fitFunc = new TF1("fitFunc", "gaus");  // Adjust fit range as needed
            hist->Fit(fitFunc, "Q");

            if (fitFunc) {
                FitResult res = SaveFitAndExtractParams(outputFile, hist, fitFunc);
                double beta = std::sqrt(1 - 1 / std::pow((4 * energy / HELIUM_MASS + 1), 2));
                results[histo].emplace_back(beta, res);
                
                // Debug: Print fit results
                std::cout << "Fit Results for " << histo << " at Energy: " << energy << " Beta: " << beta << std::endl;
                std::cout << "Mean: " << res.mean << " ± " << res.meanError << std::endl;
                std::cout << "Sigma: " << res.sigma << " ± " << res.sigmaError << std::endl;
            }
            delete fitFunc;
        } else {
            std::cerr << "Histogram " << histo << " not found in file: " << fileName << std::endl;
        }
    }
    inputFile->Close();
    outputFile->Close();
}


void PlotResults(const std::map<TString, std::vector<std::pair<double, FitResult>>>& results) {
    TCanvas* canvas = new TCanvas("c", "Fit Results", 800, 600);
    TGraphErrors* meanGraph[1];
    TGraphErrors* sigmaGraph[1];

    int colors[] = {kRed, kBlue, kGreen};  // Different colors for different fit types
    TString fitTypes[] = {"ToF_Bar9_XY"};

    for (int i = 0; i < 1; ++i) {
        meanGraph[i] = new TGraphErrors();
        sigmaGraph[i] = new TGraphErrors();
        meanGraph[i]->SetMarkerColor(colors[i]);
        sigmaGraph[i]->SetMarkerColor(colors[i]);
        meanGraph[i]->SetMarkerStyle(20 + i);
        sigmaGraph[i]->SetMarkerStyle(20 + i);

        int pointIndex = 0;
        auto it = results.find(fitTypes[i]);
        if (it != results.end()) {
            for (const auto& [beta, fitResult] : it->second) {
                meanGraph[i]->SetPoint(pointIndex, beta, fitResult.mean);
                meanGraph[i]->SetPointError(pointIndex, 0, fitResult.meanError);

                sigmaGraph[i]->SetPoint(pointIndex, beta, fitResult.sigma);
                sigmaGraph[i]->SetPointError(pointIndex, 0, fitResult.sigmaError);
                pointIndex++;
            }
        } else {
            std::cerr << "Warning: No results found for " << fitTypes[i] << std::endl;
        }
    }

    canvas->cd();
    meanGraph[0]->SetTitle("Mean ToF vs#beta");
    meanGraph[0]->GetXaxis()->SetTitle("#beta");
    meanGraph[0]->GetYaxis()->SetTitle("Mean ToF [ns]");
    //meanGraph[0]->SetMinimum(0.);
    //double xmean_max = meanGraph[0]->GetXaxis()->GetXmax();
    //meanGraph[0]->GetXaxis()->SetLimits(0., xmean_max);
    meanGraph[0]->Draw("AP");

    TLegend* legend_mean = new TLegend(0.5, 0.6, 0.7, 0.8);  // (x1, y1, x2, y2) in normalized coordinates
    //TLegend* legend_mean = new TLegend(0.2, 0.85, 0.4, 0.7);
    legend_mean->AddEntry(meanGraph[0], "Mean ToF Bar 9(X), Bar 9(Y)", "P");  // Entry for the new graph
    legend_mean->SetBorderSize(0);  // No border
    legend_mean->SetTextSize(0.03);
    legend_mean->Draw();
    canvas->SaveAs("Plots/Merged_Mean_ToF_Fit.png");
    delete legend_mean;

    canvas->Clear();
    sigmaGraph[0]->SetTitle("#sigma ToF vs#beta");
    sigmaGraph[0]->GetXaxis()->SetTitle("#beta");
    sigmaGraph[0]->GetYaxis()->SetTitle("#sigma ToF [ns]");
    sigmaGraph[0]->SetMinimum(0.1);
    //double xsigma_max = sigmaGraph[0]->GetXaxis()->GetXmax();
    //sigmaGraph[0]->GetXaxis()->SetLimits(0., xsigma_max);
    sigmaGraph[0]->Draw("AP");

    //TLegend* legend_sigma = new TLegend(0.5, 0.6, 0.7, 0.8);
    TLegend* legend_sigma = new TLegend(0.2, 0.85, 0.4, 0.7);
    legend_sigma->AddEntry(sigmaGraph[0], "#sigma ToF Bar 9(X), Bar 9(Y)", "P");  // Entry for the new graph
    legend_sigma->SetBorderSize(0);  // No border
    legend_sigma->SetTextSize(0.03);
    legend_sigma->Draw();
    canvas->SaveAs("Plots/Merged_Sigma_ToF_Fit.png");

    delete legend_sigma;
    delete canvas;
}

