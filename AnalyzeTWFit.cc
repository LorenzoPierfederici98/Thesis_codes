#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TFitResult.h>
#include <TObjString.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <cmath>
#include <TSystem.h>
#include <set>

using namespace std;

struct FitData {
    int flag=0;
    set<Double_t> energies;
    map<TString, TH1D*> SumHisto;
    map<TString, Double_t> fitMean;
    map<TString, Double_t> fitMeanErr;
    map<Double_t, vector<int>> runCombinations;
    map<Double_t, vector<pair<TString, TString>>> validCombinations; // To store valid layer-bar combinations
};

void ProcessFile(const TString &fileName, const int runNumber, FitData &heData, FitData &hData, FitData &cData);

void SumHistograms(const TString &fileName, const Double_t energy, const int runNumber, const TString &comboKey, FitData &Data);

pair<Double_t, Double_t> PerformFit(FitData &Data, TString comboKey);

void PlotChargeFitResults(const map<TString, pair<Double_t, Double_t>> &fitCharge, const TString &dataLabel);

void AnalyzeTWFit(const vector<int> &fileNumbers){
    TString baseName = "TW/AnaFOOT_TW_Decoded_HIT2022_";
    TString suffix = ".root";
    TString suffix_fit = "_Fit.root";
    TString fitresult = "FitResult_";
    TString fitcharge = "Charge_";
    vector<TString> layer = {"layer0_", "layer1_"};
    vector<TString> bar = {"bar0", "bar1", "bar2", "bar3", "bar4", "bar5", "bar6", "bar7", "bar8", "bar9", "bar10", "bar11", "bar12", "bar13", "bar14", "bar15", "bar16", "bar17", "bar18", "bar19"};

    // Separate storage for "HE", "H", and "C"
    FitData heData, hData, cData;

    // Process files to populate FitData structures
    for (int runNumber : fileNumbers) {
        TString fileName = baseName + TString::Format("%d", runNumber) + suffix_fit;
        ProcessFile(fileName, runNumber, heData, hData, cData);
    }

    if (heData.flag) {

        map<TString, pair<Double_t, Double_t>> fitChargeHe;
        // Loop over each energy
        for (Double_t energy : heData.energies) {
            // Loop over each layer-bar combination for the current energy
            for (const auto& combination : heData.validCombinations[energy]) {
                TString layer = combination.first; // Extract layer
                TString bar = combination.second;  // Extract bar
                // Create a unique key for this energy-layer-bar combination
                TString comboKey = Form("SumHisto_%g_%s%s", energy, layer.Data(), bar.Data());

                // Sum histograms for this combination across all runs
                for (int runNumber : heData.runCombinations[energy]) {
                    TString fileName = baseName + TString::Format("%d", runNumber) + suffix;
                    SumHistograms(fileName, energy, runNumber, comboKey, heData);
                }
                // Perform the fit for this specific combination
                fitChargeHe[comboKey] = PerformFit(heData, comboKey);
            }
        }

        // Plot the charge fit results
        PlotChargeFitResults(fitChargeHe, "HE");
    }

    if (hData.flag) {

        map<TString, pair<Double_t, Double_t>> fitChargeH;
        // Loop over each energy
        for (Double_t energy : heData.energies) {
            // Loop over each layer-bar combination for the current energy
            for (const auto& combination : heData.validCombinations[energy]) {
                TString layer = combination.first; // Extract layer
                TString bar = combination.second;  // Extract bar
                // Create a unique key for this energy-layer-bar combination
                TString comboKey = Form("SumHisto_%g_%s%s", energy, layer.Data(), bar.Data());

                // Sum histograms for this combination across all runs
                for (int runNumber : heData.runCombinations[energy]) {
                    TString fileName = baseName + TString::Format("%d", runNumber) + suffix;
                    SumHistograms(fileName, energy, runNumber, comboKey, heData);
                }
                // Perform the fit for this specific combination
                fitChargeH[comboKey] = PerformFit(heData, comboKey);
            }
        }

        // Plot the charge fit results
        PlotChargeFitResults(fitChargeH, "H");
    }

    if (cData.flag) {

        map<TString, pair<Double_t, Double_t>> fitChargeC;
        // Loop over each energy
        for (Double_t energy : heData.energies) {
            // Loop over each layer-bar combination for the current energy
            for (const auto& combination : heData.validCombinations[energy]) {
                TString layer = combination.first; // Extract layer
                TString bar = combination.second;  // Extract bar
                // Create a unique key for this energy-layer-bar combination
                TString comboKey = Form("SumHisto_%g_%s%s", energy, layer.Data(), bar.Data());

                // Sum histograms for this combination across all runs
                for (int runNumber : heData.runCombinations[energy]) {
                    TString fileName = baseName + TString::Format("%d", runNumber) + suffix;
                    SumHistograms(fileName, energy, runNumber, comboKey, heData);
                }
                // Perform the fit for this specific combination
                fitChargeC[comboKey] = PerformFit(heData, comboKey);
            }
        }

        // Plot the charge fit results
        PlotChargeFitResults(fitChargeC, "C");
    }

    return;
}

void ProcessFile(const TString &fileName, const int runNumber, FitData &heData, FitData &hData, FitData &cData) {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open the file " << fileName << endl;
        return;
    }

    TObjString *ionInfoObj = (TObjString *)file->Get("IonInfo");
    TObjString *energyObj = (TObjString *)file->Get("BeamEnergyInfo");

    if (!ionInfoObj || !energyObj) {
        cerr << "Error: IonInfo or BeamEnergyInfo not found in " << fileName << endl;
        file->Close();
        delete file;
        return;
    }

    TString ionType = ionInfoObj->GetString();
    Double_t beamEnergy = energyObj->GetString().Atof();

    FitData *data;
    if (ionType == "HE") {
        data = &heData;
    } else if (ionType == "H") {
        data = &hData;
    } else if (ionType == "C") {
        data = &cData;
    } else {
        cerr << "Warning: Unknown ion type " << ionType << " in " << fileName << endl;
        file->Close();
        delete file;
        return;
    }

    (*data).flag=1;
    (*data).energies.insert(beamEnergy);
    (*data).runCombinations[beamEnergy].push_back(runNumber);

    TIter nextKey(file->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextKey())) {
        TString keyName = key->GetName();

        // Check if the key matches the pattern for FitResult
        if (keyName.BeginsWith("FitResult_Charge_layer")) {
            // Extract the layer-bar combination using string parsing
            TString layer = keyName(keyName.Index("layer"), keyName.Index("_bar") - keyName.Index("layer") + 1);
            TString bar = keyName(keyName.Index("bar"), keyName.Length() - keyName.Index("bar"));

            // Store the combination (key ensures uniqueness)
            (*data).validCombinations[beamEnergy].emplace_back(layer, bar);
        }
    }
    file->Close();
    delete file;
}

void SumHistograms(const TString &fileName, const Double_t energy, const int runNumber, const TString &comboKey, FitData &Data) {
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file " << fileName << endl;
        return;
    }

    // Check if energy exists in validCombinations
    if (Data.validCombinations.find(energy) == Data.validCombinations.end()) {
        cerr << "Warning: No valid combinations for energy " << energy << endl;
        file->Close();
        delete file;
        return;
    }

    for (const auto& combination : Data.validCombinations[energy]) {
        TString layer = combination.first;
        TString bar = combination.second;
        TString DirName;

        TString layer_root;
        if (layer == "layer0_") {
            layer_root = "LayerY_";
            DirName = "ChargeTimeLayerY/";
        } else if (layer == "layer1_") {
            layer_root = "LayerX_";
            DirName = "ChargeTimeLayerX/";
        } else {
            cerr << "Invalid layer name" << endl;
            continue;
        }

        TString histoName = Form("%sCharge_%s%s", DirName.Data(), layer_root.Data(), bar.Data());
        TH1D* histo = (TH1D*)file->Get(histoName);
        if (!histo) {
            cerr << "Warning: Histogram " << histoName << " not found in " << fileName << endl;
            continue;
        }

        // Check if the histogram for this layer-bar combination exists in SumHisto
        if (!Data.SumHisto[comboKey]) {
            // Create the histogram if it doesn't exist
            Data.SumHisto[comboKey] = (TH1D*)histo->Clone(Form("SumHisto_%s", comboKey.Data()));
            Data.SumHisto[comboKey]->SetDirectory(0); // Detach from file
        } else {
            // Add the histogram to the existing one for the same layer-bar combination
            Data.SumHisto[comboKey]->Add(histo);
        }
    }

    file->Close();
    delete file;
}

pair<Double_t, Double_t> PerformFit(FitData &Data, TString comboKey) {
    // Check if the key exists in SumHisto
    if (Data.SumHisto.find(comboKey) == Data.SumHisto.end()) {
        cerr << "Error: No summed histogram found for key " << comboKey << endl;
        return {-1.0, -1.0}; // Return default values
    }

    // Retrieve the summed histogram
    TH1D* histo = Data.SumHisto[comboKey];
    if (!histo) {
        cerr << "Error: Summed histogram for key " << comboKey << " is null" << endl;
        return {-1.0, -1.0}; // Return default values
    }

    // Perform the Gaussian fit
    TFitResultPtr fitResult = histo->Fit("gaus", "QS", " ", 2., 10.);
    if (!fitResult->IsValid()) {
        cerr << "Fit failed for key " << comboKey << endl;
        return {-1.0, -1.0}; // Return default values
    }

    // Extract fit parameters
    Double_t mean = fitResult->Parameter(1);
    Double_t meanErr = fitResult->ParError(1);

    // Store the fit results
    Data.fitMean[comboKey] = mean;
    Data.fitMeanErr[comboKey] = meanErr;

    // Print fit results
    cout << "Fit Results for Key " << comboKey << ":\n"
         << "  Mean: " << mean << " ± " << meanErr << endl;

    return {mean, meanErr};
}

void PlotChargeFitResults(const map<TString, pair<Double_t, Double_t>> &fitCharge, const TString &dataLabel) {
    // Group data by layer-bar combinations
    map<TString, vector<pair<Double_t, pair<Double_t, Double_t>>>> groupedData;

    for (const auto& [key, fitResult] : fitCharge) {
        // Extract energy, layer, and bar from the key
        TString energyStr = key(9, key.Index("_", 9) - 9); // Extract energy from "SumHisto_energy_..."
        TString layerBar = key(key.Index("_", 9) + 1, key.Length()); // Extract layer-bar (e.g., "layer0_bar10")

        // Convert energy to Double_t
        Double_t energy = atof(energyStr.Data());

        // Group the data by layer-bar
        groupedData[layerBar].emplace_back(energy, fitResult);
    }

    // Create a TGraphErrors for each layer-bar combination
    for (const auto& [layerBar, data] : groupedData) {
        TGraphErrors *graph = new TGraphErrors();
        int pointIndex = 0;

        // Add all energy-charge pairs to the graph
        for (const auto& [energy, fitResult] : data) {
            Double_t meanCharge = fitResult.first;
            Double_t meanChargeErr = fitResult.second;

            if (meanCharge != -1.0 && meanChargeErr != -1.0) {
                graph->SetPoint(pointIndex, energy, meanCharge);
                graph->SetPointError(pointIndex, 0, meanChargeErr); // No x-error
                pointIndex++;
            }
        }

        // Style the graph
        graph->SetTitle(Form("Sum Fit Charge %s for %s;Beam Energy (MeV);Mean Charge [a.u.]", dataLabel.Data(), layerBar.Data()));
        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(kBlue);
        graph->SetLineColor(kBlue);

        // Adjust axis labels
        graph->GetYaxis()->SetTitle("Mean Charge [a.u.]");
        graph->GetYaxis()->SetTitleOffset(1.8); // Adjust for better visibility
        graph->GetXaxis()->SetTitle("Beam Energy (MeV)");

        // Create a canvas for this combination and draw the graph
        TCanvas *canvas = new TCanvas(Form("SumFitCharge_%s_%s", dataLabel.Data(), layerBar.Data()), Form("Sum Fit Charge %s for %s", dataLabel.Data(), layerBar.Data()), 800, 600);
        canvas->SetLeftMargin(0.18); // Increase left margin to avoid cropping
        canvas->SetBottomMargin(0.12); // Increase bottom margin for clarity

        graph->Draw("AP");

        // Save the plot with the correct name format
        canvas->SaveAs(Form("Plots/SumFitCharge_%s_%s.png", dataLabel.Data(), layerBar.Data()));

        // Clean up memory
        delete graph;
        delete canvas;
    }
}



