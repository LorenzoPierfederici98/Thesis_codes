// Macro that plots the fit results on protons and helium peaks given by the AnalyzeTWFragm.cc macro vs the energies
// retrieved from MC (2 energy loss values per bar and beam energy). The fit charge values vs energy loss values
// are fitted with a 1 parameter linear function. To be run with root -l -b -q 'CalibrateFragm.cc()'

#include "CalibrateFragm.h"

void CalibrateFragm() {
    std::vector<std::pair<std::string, double>> filesAndEnergies = {
        {"TW/AnaFOOT_TW_Decoded_HIT2022_fragm_100MeV_Fit.root", 100},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_fragm_140MeV_Fit.root", 140},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_fragm_200MeV_Fit.root", 200},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_fragm_220MeV_Fit.root", 220}
    };
    std::map<int, double>elossP = {{100, 2.721}, {140, 2.020}, {200, 1.4761}, {220, 1.4267}};  // energy loss for protons, from MC
    std::map<int, double>elossHe = {{100, 9.7638}, {140, 7.2929}, {200, 5.4135}, {220, 5.2352}};  // energy loss for heliums, from MC

    std::map<TString, std::map<int, double>> fitMeansP;  // protons
    std::map<TString, std::map<int, double>> fitErrorsP;
    std::map<TString, std::map<int, double>> fitMeansHe;  // heliums
    std::map<TString, std::map<int, double>> fitErrorsHe;

    for (const auto& [fileName, energy] : filesAndEnergies) {
        ProcessFile(fileName, energy, fitMeansP, fitMeansHe, fitErrorsP, fitErrorsHe);
    }

    PlotFitResultsCombined(fitMeansP, fitErrorsP, fitMeansHe, fitErrorsHe, elossP, elossHe);
}

void ProcessFile(const TString& fileName, 
                 int energy, 
                 std::map<TString, std::map<int, double>>& fitMeansP, 
                 std::map<TString, std::map<int, double>>& fitMeansHe, 
                 std::map<TString, std::map<int, double>>& fitErrorsP, 
                 std::map<TString, std::map<int, double>>& fitErrorsHe) {
    // Open the input ROOT file
    TFile* inputFile = TFile::Open(fileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return;
    }

    // Temporary storage for layer-bar combinations with both fits
    std::map<TString, std::pair<TFitResult *, TFitResult *>> validResults;

    // Loop through all objects in the file
    TIter nextKey(inputFile->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextKey())) {
        TString objectName = key->GetName();

        // Check if the object corresponds to a fit result for Proton or Helium
        bool isProton = objectName.BeginsWith("FitResultP_");
        bool isHelium = objectName.BeginsWith("FitResultHe_");
        if (!isProton && !isHelium) continue;

        // Extract the layer-bar name (e.g., "LayerX_bar9")
        TString layerBarName = objectName;
        layerBarName.Remove(0, objectName.First('_') + 1);

        // Retrieve the fit result object
        TFitResult *fitResult = (TFitResult *)inputFile->Get(objectName);
        if (!fitResult || !fitResult->IsValid()) continue;

        // Store the fit result in the temporary map
        if (isProton) {
            validResults[layerBarName].first = fitResult;
        } else if (isHelium) {
            validResults[layerBarName].second = fitResult;
        }
    }

    // Store valid results for layer-bar combinations with both fits
    for (const auto& [layerBarName, fitPair] : validResults) {
        const TFitResult *fitP = fitPair.first;
        const TFitResult *fitHe = fitPair.second;

        // Only proceed if both Proton and Helium fits are available
        if (!fitP || !fitHe) continue;

        // Extract mean and error for Proton
        double meanP = fitP->Parameter(1);
        double errorP = fitP->ParError(1);

        // Extract mean and error for Helium
        double meanHe = fitHe->Parameter(1);
        double errorHe = fitHe->ParError(1);

        // Store the results in the maps
        fitMeansP[layerBarName][energy] = meanP;
        fitErrorsP[layerBarName][energy] = errorP;
        fitMeansHe[layerBarName][energy] = meanHe;
        fitErrorsHe[layerBarName][energy] = errorHe;
    }

    // Close the input file
    inputFile->Close();
    delete inputFile;
}

void PlotFitResultsCombined(
    const std::map<TString, std::map<int, double>>& fitMeansP,
    const std::map<TString, std::map<int, double>>& fitErrorsP,
    const std::map<TString, std::map<int, double>>& fitMeansHe,
    const std::map<TString, std::map<int, double>>& fitErrorsHe,
    const std::map<int, double>& elossP,
    const std::map<int, double>& elossHe) {

    // Define the energies and layers
    std::vector<int> energies = {100, 140, 200, 220};
    std::vector<TString> layers = {"X", "Y"};

    for (TString layer : layers) {
        for (int bar = 0; bar < 20; ++bar) {
            TString layerBarName = Form("Layer%s_bar%d", layer.Data(), bar);

            TGraphErrors* graphProtons = new TGraphErrors();
            TGraphErrors* graphHeliums = new TGraphErrors();
            TMultiGraph* multiGraph = new TMultiGraph();

            int indexP = 0;
            int indexHe = 0;

            // Add proton points
            for (int energy : energies) {
                if (fitMeansP.count(layerBarName) && fitMeansP.at(layerBarName).count(energy)) {
                    double meanP = fitMeansP.at(layerBarName).at(energy);
                    double errorP = fitErrorsP.at(layerBarName).at(energy);
                    graphProtons->SetPoint(indexP, elossP.at(energy), meanP);
                    graphProtons->SetPointError(indexP, 0, errorP);
                    ++indexP;
                }
            }

            // Add helium points
            for (int energy : energies) {
                if (fitMeansHe.count(layerBarName) && fitMeansHe.at(layerBarName).count(energy)) {
                    double meanHe = fitMeansHe.at(layerBarName).at(energy);
                    double errorHe = fitErrorsHe.at(layerBarName).at(energy);
                    graphHeliums->SetPoint(indexHe, elossHe.at(energy), meanHe);
                    graphHeliums->SetPointError(indexHe, 0, errorHe);
                    if (layerBarName == "LayerX_bar1" && energy == 220) {
                        // Setting the fit point for LayerX_bar9 @220 MeV/u as an outlier.
                        // Points with 0 error are ignored in the fit.
                        graphHeliums->SetPointError(indexHe, 0, 0); 
                    }
                    ++indexHe;
                }
            }

            if (graphProtons->GetN() == 0 && graphHeliums->GetN() == 0) {
                delete graphProtons;
                delete graphHeliums;
                delete multiGraph;
                continue;
            }

            graphProtons->SetMarkerStyle(20);
            graphProtons->SetMarkerColor(kBlue);
            graphHeliums->SetMarkerStyle(20);
            graphHeliums->SetMarkerColor(kRed);

            multiGraph->Add(graphProtons, "P");
            multiGraph->Add(graphHeliums, "P");

            TCanvas* c = new TCanvas(Form("c_%s", layerBarName.Data()), 
                                      Form("Fit Results - %s", layerBarName.Data()), 
                                      800, 600);

            multiGraph->SetTitle(Form("%s - Fit Results", layerBarName.Data()));
            multiGraph->GetXaxis()->SetTitle("Energy loss MC [MeV]");
            multiGraph->GetYaxis()->SetTitle("Mean Charge [a.u.]");
            multiGraph->Draw("A");

            // Combine data for fitting
            TGraphErrors* combinedGraph = new TGraphErrors();
            for (int i = 0; i < graphProtons->GetN(); ++i) {
                double x, y;
                graphProtons->GetPoint(i, x, y);
                combinedGraph->SetPoint(combinedGraph->GetN(), x, y);
                combinedGraph->SetPointError(combinedGraph->GetN() - 1, 0, graphProtons->GetErrorY(i));
            }
            for (int i = 0; i < graphHeliums->GetN(); ++i) {
                double x, y;
                graphHeliums->GetPoint(i, x, y);
                combinedGraph->SetPoint(combinedGraph->GetN(), x, y);
                combinedGraph->SetPointError(combinedGraph->GetN() - 1, 0, graphHeliums->GetErrorY(i));
            }

            // Perform the combined fit
            double fit_min = 1.4267;
            double fit_max = 9.7638;

            //TF1* combinedFit = new TF1("combinedFit", "[0]*x", 0., fit_max);
            TF1* combinedFit = new TF1("combinedFit", "[0]*x/(1 + [1]*x + [2]*x**2)", fit_min, fit_max);
            combinedGraph->Fit(combinedFit);
            combinedFit->SetLineColor(kRed);

            // Draw the fit function
            combinedFit->Draw("same");

            // Add legend
            TLegend* legend = new TLegend(0.15, 0.68, 0.48, 0.85);
            legend->AddEntry(graphProtons, "Protons", "p");
            legend->AddEntry(graphHeliums, "Heliums", "p");
            auto [p0, p0Error] = RoundMeasurement(combinedFit->GetParameter(0), combinedFit->GetParError(0));
            auto [p1, p1Error] = RoundMeasurement(combinedFit->GetParameter(1), combinedFit->GetParError(1));
            auto [p2, p2Error] = RoundMeasurement(combinedFit->GetParameter(2), combinedFit->GetParError(2));
            
            legend->AddEntry((TObject*)0, Form("p_{0} = %s#pm%s", p0.c_str(), p0Error.c_str()), "");
            legend->AddEntry((TObject*)0, Form("p_{1} = %s#pm%s", p1.c_str(), p1Error.c_str()), "");
            legend->AddEntry((TObject*)0, Form("p_{2} = %s#pm%s", p2.c_str(), p2Error.c_str()), "");
            legend->SetTextSize(0.02);
            legend->Draw();

            c->SaveAs(Form("Plots/Fragmentation_%s.png", layerBarName.Data()));

            delete combinedFit;
            delete combinedGraph;
            delete graphProtons;
            delete graphHeliums;
            delete multiGraph;
            delete legend;
            delete c;
        }
    }
}



pair<std::string, std::string> RoundMeasurement(double value, double uncertainty) {
    // Determine the order of magnitude of the uncertainty
    int uncertaintyOrder = (int)std::floor(std::log10(uncertainty));
    double roundingFactor = std::pow(10, uncertaintyOrder);

    // Round uncertainty to 1 significant figure
    double roundedUncertainty = std::round(uncertainty / roundingFactor) * roundingFactor;

    // Adjust the value to match the precision of the uncertainty
    double roundedValue = std::round(value / roundingFactor) * roundingFactor;

    // Determine the number of decimal places to display based on the rounded uncertainty
    int decimalPlaces = std::max(0, -uncertaintyOrder);

    // Format value and uncertainty into strings with the required precision
    std::ostringstream valueStream, uncertaintyStream;
    valueStream << std::fixed << std::setprecision(decimalPlaces) << roundedValue;
    uncertaintyStream << std::fixed << std::setprecision(decimalPlaces) << roundedUncertainty;

    return {valueStream.str(), uncertaintyStream.str()};
}
