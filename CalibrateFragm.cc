// Macro that plots the fit results on protons and helium peaks given by the AnalyzeTWFragm.cc macro vs the energies
// retrieved from MC (2 energy loss values per bar and beam energy). The fit charge values vs energy loss values
// are fitted with a 1 parameter linear function; these parameters (on per bar) are written in the configuration file.
// To be run with root -l -b -q 'CalibrateFragm.cc()'

#include "CalibrateFragm.h"

void CalibrateFragm() {
    std::vector<std::pair<std::string, int>> filesAndEnergies = {
        {"TW/AnaFOOT_TW_Decoded_HIT2022_fragm_100MeV_Fit.root", 100},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_fragm_140MeV_Fit.root", 140},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_fragm_200MeV_Fit.root", 200},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_fragm_220MeV_Fit.root", 220}
    };
    // Smeared MC values
    // std::map<int, double>elossP = {{100, 2.721}, {140, 2.020}, {200, 1.4761}, {220, 1.4267}};  // energy loss for protons, from MC
    // std::map<int, double>elossHe = {{100, 9.7638}, {140, 7.2929}, {200, 5.4135}, {220, 5.2352}};  // energy loss for heliums, from MC

    std::map<int, double>elossP = {{100, 2.621}, {140, 1.901}, {200, 1.410}, {220, 1.325}};  // energy loss for protons, from MC
    std::map<int, double>elossHe = {{100, 9.7371}, {140, 7.2586}, {200, 5.5314}, {220, 5.1694}};  // energy loss for heliums, from MC

    std::map<TString, std::map<int, double>> fitMeansP;  // protons
    std::map<TString, std::map<int, double>> fitErrorsP;
    std::map<TString, std::map<int, double>> fitMeansHe;  // heliums
    std::map<TString, std::map<int, double>> fitErrorsHe;

    for (const auto& [fileName, energy] : filesAndEnergies) {
        ProcessFile(fileName, energy, fitMeansP, fitMeansHe, fitErrorsP, fitErrorsHe);
    }

    std::map<TString, double> fitValues = PlotFitResultsCombined(fitMeansP, fitErrorsP, fitMeansHe, fitErrorsHe, elossP, elossHe);
    WriteFitValuesOrdered(fitValues);
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
        if (layerBarName == "LayerX_bar9") {
            cout << "meanP: " << meanP << " errorP: " << errorP << endl;
            cout << "meanHe: " << meanHe << " errorHe: " << errorHe << endl;
        }
    }

    // Close the input file
    inputFile->Close();
    delete inputFile;
}

std::map<TString, double> PlotFitResultsCombined(
    const std::map<TString, std::map<int, double>>& fitMeansP,
    const std::map<TString, std::map<int, double>>& fitErrorsP,
    const std::map<TString, std::map<int, double>>& fitMeansHe,
    const std::map<TString, std::map<int, double>>& fitErrorsHe,
    const std::map<int, double>& elossP,
    const std::map<int, double>& elossHe) {

    // Define the energies and layers
    std::vector<int> energies = {100, 140, 200, 220};
    std::vector<TString> layers = {"X", "Y"};

    // Stores the p0 fit values, to write them on an output file
    std::map<TString, double> fitValues;

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

            graphProtons->SetPoint(indexP+1, 0., 0.);

            // Add helium points
            for (int energy : energies) {
                if (fitMeansHe.count(layerBarName) && fitMeansHe.at(layerBarName).count(energy)) {
                    double meanHe = fitMeansHe.at(layerBarName).at(energy);
                    double errorHe = fitErrorsHe.at(layerBarName).at(energy);
                    graphHeliums->SetPoint(indexHe, elossHe.at(energy), meanHe);
                    graphHeliums->SetPointError(indexHe, 0, errorHe);
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

            multiGraph->SetTitle(Form("%s - Fit Results true MC", layerBarName.Data()));
            multiGraph->GetXaxis()->SetTitle("Energy loss MC [MeV]");
            multiGraph->GetYaxis()->SetTitle("Mean Charge [a.u.]");
            multiGraph->SetMinimum(0.);
            //multiGraph->GetXaxis()->SetRangeUser(0., 15.);
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
            //double fit_min = elossP.at(220);  // elossP[220] cannot be used on const, since it allows insertion if the key doesn't exist
            double fit_min = 0.;
            double fit_max = elossHe.at(100);

            double maxX = std::numeric_limits<double>::lowest();
            double maxY = std::numeric_limits<double>::lowest();

            // Loop through all points in the graph to find the true max
            for (int i = 0; i < combinedGraph->GetN(); ++i) {
                double x, y;
                combinedGraph->GetPoint(i, x, y);

                if (x > maxX) maxX = x;
                if (y > maxY) maxY = y;
            }

            // Calculate the initial slope for the fit
            //double initialSlope = (maxY - minY) / (maxX - minX);
            double initialSlope = maxY / maxX;

            TF1* combinedFit = new TF1("combinedFit", "[0]*x", fit_min, fit_max);
            //TF1* combinedFit = new TF1("combinedFit", "[0]*x - [0]*[1]*x**2 + [0]*[1]**2*x**3", fit_min, 2.721);
            combinedFit->SetParameters(initialSlope);  // Linear
            //combinedFit->SetParameters(0.6, -0.09, 0.01);
            //combinedFit->SetParameters(5., 1., 0.6, -2.);  // Birks + linear
            //combinedFit->SetParameters(0.6, -0.5, 0.03);  // MacLaurin Modified Birks
            combinedGraph->Fit(combinedFit, "R");
            combinedFit->SetLineColor(kRed);

            // Draw the fit function
            combinedFit->Draw("same");

            // Add legend
            TLegend* legend = new TLegend(0.15, 0.68, 0.48, 0.85);
            legend->AddEntry(graphProtons, "Protons", "p");
            legend->AddEntry(graphHeliums, "Heliums", "p");
            auto [p0, p0Error] = RoundMeasurement(combinedFit->GetParameter(0), combinedFit->GetParError(0));
            //auto [p1, p1Error] = RoundMeasurement(combinedFit->GetParameter(1), combinedFit->GetParError(1));
            //auto [p2, p2Error] = RoundMeasurement(combinedFit->GetParameter(2), combinedFit->GetParError(2));
            //auto [p3, p3Error] = RoundMeasurement(combinedFit->GetParameter(3), combinedFit->GetParError(3));
            fitValues[layerBarName] = combinedFit->GetParameter(0);
            
            legend->AddEntry((TObject*)0, Form("p_{0} [a.u./MeV] = %s#pm%s", p0.c_str(), p0Error.c_str()), "");
            //legend->AddEntry((TObject*)0, Form("p_{1} = %s#pm%s", p1.c_str(), p1Error.c_str()), "");
            //legend->AddEntry((TObject*)0, Form("p_{2} = %s#pm%s", p2.c_str(), p2Error.c_str()), "");
            //legend->AddEntry((TObject*)0, Form("p_{3} = %s#pm%s", p3.c_str(), p3Error.c_str()), "");
            legend->SetTextSize(0.025);
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
    return fitValues;
}

void WriteFitValuesOrdered(const std::map<TString, double>& fitValues) {
    // Define the output file name inside the function
    const std::string outputFileName = "TATW_Energy_calibration_perBar_HIT2022.cal";

    // Create a vector of pairs to allow sorting
    std::vector<std::tuple<int, TString, double>> sortedValues;

    // Parse the bar ID and store it along with the fit values
    for (const auto& entry : fitValues) {
        const TString& layerBarName = entry.first;
        double p0 = entry.second;

        // Extract the bar ID from the name (e.g., "LayerX_bar0" or "LayerY_bar1")
        int barID = -1;
        bool isLayerX = false;
        if (sscanf(layerBarName.Data(), "LayerX_bar%d", &barID) == 1) {
            isLayerX = true;
        } else if (sscanf(layerBarName.Data(), "LayerY_bar%d", &barID) != 1) {
            std::cerr << "Error parsing bar ID from layerBarName: " << layerBarName.Data() << std::endl;
            continue;
        }

        // Calculate the combined bar ID: 0–19 for Layer Y, 20–39 for Layer X
        int combinedBarID = barID + (isLayerX ? 20 : 0);

        // Add to the sortedValues vector
        sortedValues.emplace_back(combinedBarID, layerBarName, p0);
    }

    // Sort by the combined bar ID
    // The lambda function is used to compare the first element of the tuple (get<0> a and get<0> b) 
    std::sort(sortedValues.begin(), sortedValues.end(),
              [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });

    // Open the output file
    std::ofstream outFile(outputFileName);
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << outputFileName << std::endl;
        return;
    }

    // Write header
    outFile << "#BarId(Pisa)          p0          p1 layer(SHOE)" << std::endl;
    outFile << "#-+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n";

    // Write data
    for (const auto& entry : sortedValues) {
        int combinedBarID = std::get<0>(entry);  // Combined bar ID
        const TString& layerBarName = std::get<1>(entry);
        double p0 = std::get<2>(entry);

        // Determine the layer (0 for Y, 1 for X)
        int layerShoe = layerBarName.Contains("LayerX") ? 1 : 0;

        // Write output (p1 is always 0)
        outFile << std::setw(10) << combinedBarID
                << std::setw(20) << p0
                << std::setw(12) << 0.0
                << std::setw(10) << layerShoe
                << std::endl;
    }

    outFile << "#-+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n";
    outFile.close();
    std::cout << "Fit values written to file: " << outputFileName << std::endl;
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
