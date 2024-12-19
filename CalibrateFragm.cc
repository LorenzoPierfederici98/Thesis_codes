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

void PlotFitResultsCombined(const std::map<TString, std::map<int, double>>& fitMeansP,
                            const std::map<TString, std::map<int, double>>& fitErrorsP,
                            const std::map<TString, std::map<int, double>>& fitMeansHe,
                            const std::map<TString, std::map<int, double>>& fitErrorsHe,
                            std::map<int, double> elossP,
                            std::map<int, double> elossHe) {
    // Define the energies and layers
    std::vector<int> energies = {100, 140, 200, 220};  // Different energies
    std::vector<TString> layers = {"X", "Y"};  // Layers X and Y
    int bar = 9;  // Only for bar 9
    
    // Loop over each energy and layer combination
    for (int energy : energies) {
        for (TString layer : layers) {
            TString layerBarName = Form("Layer%s_bar%d", layer.Data(), bar);

            // Check if the data for both Proton and Helium exist for this energy and layer-bar combination
            if (fitMeansP.find(layerBarName) == fitMeansP.end() || fitMeansHe.find(layerBarName) == fitMeansHe.end())
                continue;

            // Create a graph for Proton and Helium in the same plot
            TGraphErrors *graph = new TGraphErrors();
            // Create a canvas for this energy and layer-bar combination
            TCanvas *c = new TCanvas(Form("c_%s_%d", layerBarName.Data(), energy), 
                                      Form("Fit Results - %s, Energy %d MeV", layerBarName.Data(), energy), 
                                      800, 600);

            int index = 0;

            // Ensure both P and He data exist for the current energy
            if (fitMeansP.at(layerBarName).find(energy) != fitMeansP.at(layerBarName).end() &&
                fitMeansHe.at(layerBarName).find(energy) != fitMeansHe.at(layerBarName).end()) {

                double meanP = fitMeansP.at(layerBarName).at(energy);
                double errorP = fitErrorsP.at(layerBarName).at(energy);

                double meanHe = fitMeansHe.at(layerBarName).at(energy);
                double errorHe = fitErrorsHe.at(layerBarName).at(energy);

                // Add Proton data (Red marker)
                graph->SetPoint(index, elossP[energy], meanP);
                graph->SetPointError(index, 0, errorP);
                graph->SetMarkerColor(kRed);  // Proton color
                index++;

                // Add Helium data (Blue marker)
                graph->SetPoint(index, elossHe[energy], meanHe);
                graph->SetPointError(index, 0, errorHe);
                graph->SetMarkerColor(kBlue);  // Helium color
                index++;
            }

            // Skip plotting if no common points exist
            if (graph->GetN() == 0) {
                delete graph;
                continue;
            }

            // Set graph style and draw
            graph->SetMarkerStyle(20);
            //graph->SetLineColor(kBlack);
            graph->SetTitle(Form("%s beam energy: %d MeV", layerBarName.Data(), energy));
            graph->GetXaxis()->SetTitle("Energy loss MC [MeV]");
            graph->GetYaxis()->SetTitle("Mean Charge [a.u.]");
            graph->Draw("AP");
            graph->SetMinimum(0.);
            //graph->GetXaxis()->SetRangeUser(0., 10.);

            // Fit the graph with the desired function: [0]*x/(1 + [1]*x)
            //TF1* fitFunc = new TF1("fitFunc", "[0]*x/(1 + [1]*x)", 0., 20.);
            TF1* fitFunc = new TF1("fitFunc", "[0]*x");
            graph->Fit(fitFunc, "QR+");
            // Extract the fit parameters: [0] and [1]
            auto [param0, param0Error] = RoundMeasurement(fitFunc->GetParameter(0), fitFunc->GetParError(0));
            //auto [param1, param1Error] = RoundMeasurement(fitFunc->GetParameter(1), fitFunc->GetParError(1));

            //cout << "energy: " << energy << " " << layerBarName << " " << Form("(p0, p1) = (%s, %s)", param0.c_str(), param1.c_str()) << Form(" chisq/ndof: %.1e/%d ", fitFunc->GetChisquare(), fitFunc->GetNDF()) << endl;

            // Add fit results to the plot using TPaveText
            TPaveText *pave = new TPaveText(0.15, 0.7, 0.45, 0.85, "NDC");  // Normalized coordinates
            pave->SetFillColor(0);  // Transparent background
            pave->SetLineColor(0);  // Transparent border
            pave->SetShadowColor(0); // Remove the shadow
            pave->SetTextAlign(12); // Align left
            //pave->AddText(Form("Fit results:"));
            pave->AddText(Form("slope [a.u./MeV] = %s#pm %s", param0.c_str(), param0Error.c_str()));
            pave->SetTextSize(0.035);
            //pave->AddText(Form("p_{0} [a.u./MeV] = %s#pm %s", param0.c_str(), param0Error.c_str()));
            //pave->AddText(Form("p_{1} [1/MeV] = %s#pm %s", param1.c_str(), param1Error.c_str()));
            pave->Draw();

            c->SaveAs(Form("Plots/Fragmentation_%s_%dMeV.png", layerBarName.Data(), energy));

            delete c; // Clean up the canvas
            delete graph; // Clean up the graph
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
