#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TObjString.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <cmath>
#include <TSystem.h>

using namespace std;

struct FitData {
    map<Double_t, vector<vector<Double_t>>> sumFitMean;
    map<Double_t, vector<vector<Double_t>>> sumFitMeanErr;
    map<Double_t, int> count;
    vector<Double_t> energies;
    vector<vector<vector<Double_t>>> fitMean;
    vector<vector<vector<Double_t>>> fitMeanErr;
    map<Double_t, vector<pair<int, int>>> validCombinations; // To store valid layer-bar combinations
};

void InitializeFitData(FitData &data, Double_t beamEnergy, int layerlength, int barlength) {
    if (data.sumFitMean.find(beamEnergy) == data.sumFitMean.end()) {
        data.sumFitMean[beamEnergy] = vector<vector<Double_t>>(layerlength, vector<Double_t>(barlength, 0.0));
        data.sumFitMeanErr[beamEnergy] = vector<vector<Double_t>>(layerlength, vector<Double_t>(barlength, 0.0));
        data.count[beamEnergy] = 0;
    }
}

void ProcessFile(const TString &fileName, const TString &fitresult, const TString &fitcharge,
                 const vector<TString> &layer, const vector<TString> &bar, FitData &heData, FitData &hData, FitData &cData) {
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

    InitializeFitData(*data, beamEnergy, layer.size(), bar.size());

    for (size_t i = 0; i < layer.size(); i++) {
        for (size_t j = 0; j < bar.size(); j++) {
            TFitResult *fitresultcharge = (TFitResult *)file->Get(fitresult + fitcharge + layer[i] + bar[j]);
            if (fitresultcharge) {
                Double_t mean = fitresultcharge->Parameter(1);
                Double_t meanErr = fitresultcharge->ParError(1);

                (*data).sumFitMean[beamEnergy][i][j] += mean;
                (*data).sumFitMeanErr[beamEnergy][i][j] += meanErr * meanErr; // Square the error
                (*data).validCombinations[beamEnergy].emplace_back(i, j); // Store valid combinations
            } else {
                cerr << "Warning: Fit result not found for " << layer[i] << bar[j] << " in " << fileName << endl;
            }
        }
    }

    (*data).count[beamEnergy]++;
    file->Close();
    delete file;
}

void ComputeMeans(FitData &data, int layerlength, int barlength) {
    data.fitMean = vector<vector<vector<Double_t>>>(layerlength, vector<vector<Double_t>>(barlength));
    data.fitMeanErr = vector<vector<vector<Double_t>>>(layerlength, vector<vector<Double_t>>(barlength));

    for (const auto &entry : data.sumFitMean) {
        Double_t beamEnergy = entry.first;
        data.energies.push_back(beamEnergy);
        for (int i = 0; i < layerlength; i++) {
            for (int j = 0; j < barlength; j++) {
                if (data.count[beamEnergy] > 0) {
                    data.fitMean[i][j].push_back(data.sumFitMean[beamEnergy][i][j] / data.count[beamEnergy]);
                    data.fitMeanErr[i][j].push_back(sqrt(data.sumFitMeanErr[beamEnergy][i][j]) / data.count[beamEnergy]);
                }
            }
        }
    }
}

void PrintFinalMeans(const FitData &data, const vector<TString> &layer, const vector<TString> &bar, const TString &ionType) {
    if (data.energies.empty()) {
        cout << "No data available for ion type: " << ionType << "\n";
        return; // No data to print
    }

    cout << "Final mean values for ion type: " << ionType << "\n";
    for (size_t i = 0; i < layer.size(); ++i) {
        for (size_t j = 0; j < bar.size(); ++j) {
            cout << "Layer: " << layer[i] << ", Bar: " << bar[j] << "\n";
            for (size_t k = 0; k < data.energies.size(); ++k) {
                if (data.fitMean[i][j].size() > k) {
                    cout << "  Energy: " << data.energies[k]
                         << " MeV, Mean: " << data.fitMean[i][j][k]
                         << ", Error: " << data.fitMeanErr[i][j][k] << "\n";
                }
            }
        }
    }
}

void SavePlots(FitData &data, const vector<TString> &layer, const vector<TString> &bar,
               const TString &type, int markerColor) {
    if (data.energies.empty()) return;

    // Create the directory if it doesn't exist
    TString directory = "Plots";
    gSystem->mkdir(directory, true);  // 'true' allows creating nested directories if needed

    TCanvas *canvas = new TCanvas("canvas", "Fit Results", 800, 600);

    for (const auto &entry : data.validCombinations) {
        Double_t beamEnergy = entry.first;
        for (const auto &[i, j] : entry.second) {
            if (data.fitMean[i][j].empty()) continue;

            TGraphErrors *graph = new TGraphErrors();
            bool validData = false;
            int pointIndex = 0;
            for (size_t k = 0; k < data.energies.size(); ++k) {
                if (!isfinite(data.fitMean[i][j][k]) || !isfinite(data.fitMeanErr[i][j][k]) || data.fitMean[i][j][k] == 0)
                    continue; // Skip invalid data
                graph->SetPoint(pointIndex, data.energies[k], data.fitMean[i][j][k]);
                graph->SetPointError(pointIndex, 0, data.fitMeanErr[i][j][k]);
                validData = true;
                pointIndex++;
            }

            if (validData) {
                TString graphName = TString::Format("Graph_%s_Layer%d_%s", type.Data(), i, bar[j].Data());
                graph->SetName(graphName);
                graph->SetTitle(graphName);
                graph->SetMarkerStyle(24);
                graph->SetMarkerSize(1.2);
                graph->SetMarkerColor(markerColor);
                graph->GetXaxis()->SetTitle("Primary beam energy [MeV]");
                graph->GetYaxis()->SetTitle("Mean charge");
                graph->Draw("AP");

                // Save the plot inside the directory
                TString outputFileName = directory + "/" + TString::Format("plot_%s_Layer%d_%s.png", type.Data(), i, bar[j].Data());
                canvas->SaveAs(outputFileName);
            }
            delete graph;
        }
    }
    delete canvas;
}

void AnalyzeFit(const vector<int> &fileNumbers) {
    // Define file names and paths
    TString baseName = "AnaFOOT_Merge_HIT2022_";
    TString suffix = "_Fit.root";
    TString fitresult = "FitResult_";
    TString fitcharge = "Charge_";
    vector<TString> layer = {"layer0_", "layer1_"};
    vector<TString> bar = {"bar0", "bar1", "bar2", "bar3", "bar4", "bar5", "bar6", "bar7", "bar8", "bar9", "bar10", "bar11", "bar12", "bar13", "bar14", "bar15", "bar16", "bar17", "bar18", "bar19"};

    // Separate storage for "HE", "H", and "C"
    FitData heData, hData, cData;

    // Process each file
    for (int number : fileNumbers) {
        TString fileName = baseName + TString::Format("%d", number) + suffix;
        ProcessFile(fileName, fitresult, fitcharge, layer, bar, heData, hData, cData);
    }

    // Process "HE" IonType if data is present in sumFitMean
    if (!heData.sumFitMean.empty()) {
        ComputeMeans(heData, layer.size(), bar.size());
        PrintFinalMeans(heData, layer, bar, "HE");
        SavePlots(heData, layer, bar, "HE", kBlue);
    }

    // Process "H" IonType if data is present in sumFitMean
    if (!hData.sumFitMean.empty()) {
        ComputeMeans(hData, layer.size(), bar.size());
        PrintFinalMeans(hData, layer, bar, "H");
        SavePlots(hData, layer, bar, "H", kRed);
    }

    // Process "C" IonType if data is present in sumFitMean
    if (!cData.sumFitMean.empty()) {
        ComputeMeans(cData, layer.size(), bar.size());
        PrintFinalMeans(cData, layer, bar, "C");
        SavePlots(cData, layer, bar, "C", kGreen);
    }
}
