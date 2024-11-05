//Macro to analyze the fit results on ToF on the entire TW and on the central bars (8, 9, 10),
//performed by the AnalyzeTWChargeTime.cc macro. A mean ToF value and a mean sigma value from the
//gaussian fit are performed if there are multiple values for a given energy. The mean ToF values and
//sigmas are plotted against the beam energy, both for the entire TW and for the central bars only.
//To be run with e.g. root -l 'AnalyzeTofFit.cc({4742, 4743, 4828})'


#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TObjString.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <iostream>
#include <map>
#include <cmath>

using namespace std;

// Function to compute the mean and error
pair<double, double> ComputeMeanAndError(const vector<double> &values, const vector<double> &errors) {
    double sum = 0.0, squaredErrorSum = 0.0;
    for (size_t i = 0; i < values.size(); ++i) {
        sum += values[i];
        squaredErrorSum += errors[i] * errors[i];
    }
    return {sum / values.size(), sqrt(squaredErrorSum)};
}

// Modified PlotGraph function to include a third graph (Bar9XY)
void PlotGraph(TGraphErrors* graph1, TGraphErrors* graph2, TGraphErrors* graph3, TString canvasTitle, TString graphTitle, TString yAxisTitle, TString legend2, TString legend3, double ylim) {
    TCanvas* c1 = new TCanvas("c1", canvasTitle, 800, 600);
    graph1->SetMarkerStyle(24); graph1->SetMarkerColor(kBlue);
    graph1->SetTitle(graphTitle); 
    graph1->GetXaxis()->SetTitle("Beam Energy [MeV]");
    graph1->GetYaxis()->SetTitle(yAxisTitle + " [ns]"); 
    graph1->GetYaxis()->SetRangeUser(0, ylim);
    graph1->Draw("AP");

    graph2->SetMarkerStyle(24); graph2->SetMarkerColor(kRed);
    graph2->Draw("P same");

    graph3->SetMarkerStyle(24); graph3->SetMarkerColor(kGreen);  // New color for the third graph
    graph3->Draw("P same");

    // Updated legend with entries for all three graphs
    TLegend* legend = new TLegend(0.7, 0.3, 0.9, 0.5);
    legend->AddEntry(graph1, yAxisTitle, "P");
    legend->AddEntry(graph2, legend2, "P");
    legend->AddEntry(graph3, legend3, "P");  // Entry for the new graph
    legend->Draw();

    c1->SaveAs("Plots/" + canvasTitle + ".png");
    delete c1;
}


// Function to analyze TOF fit results
void AnalyzeTofFit(const vector<int> &fileNumbers) {
    TString baseName = "TW/AnaFOOT_TW_Decoded_HIT2022_";
    TString suffix = "_Fit.root";
    TString fitresult_tof = "TW_Fit";
    TString fitresult_tof_centralbars = "TW_Fit_CentralBars";
    TString fitresult_tof_xy = "TW_Fit_Bar9XY"; // New variable for the fit result of bar 9

    map<int, vector<double>> meanTofMap, meanTofCentralBarsMap, meanTofXYMap, sigmaTofMap, sigmaTofCentralBarsMap, sigmaTofXYMap;
    map<int, vector<double>> meanTofErrorMap, meanTofCentralBarsErrorMap, meanTofXYErrorMap, sigmaTofErrorMap, sigmaTofCentralBarsErrorMap, sigmaTofXYErrorMap;
    vector<double> beamEnergies;

    for (int number : fileNumbers) {
        TString fileName = baseName + TString::Format("%d", number) + suffix;
        TFile *inFile = TFile::Open(fileName, "READ");
        if (!inFile || inFile->IsZombie()) continue;

        TObjString *energyObj = (TObjString *)inFile->Get("BeamEnergyInfo");
        TFitResult *fitresult = (TFitResult *)inFile->Get(fitresult_tof);
        TFitResult *fitresult_centralbars = (TFitResult *)inFile->Get(fitresult_tof_centralbars);
        TFitResult *fitresult_xy = (TFitResult *)inFile->Get(fitresult_tof_xy); // Access the new fit result for bar 9

        if (!energyObj || !fitresult || !fitresult_centralbars || !fitresult_xy) {
            inFile->Close(); delete inFile; continue;
        }

        int beamEnergy = energyObj->GetString().Atof();
        if (find(beamEnergies.begin(), beamEnergies.end(), beamEnergy) == beamEnergies.end()) beamEnergies.push_back(beamEnergy);

        meanTofMap[beamEnergy].push_back(fitresult->Parameter(1));
        meanTofErrorMap[beamEnergy].push_back(fitresult->ParError(1));
        meanTofCentralBarsMap[beamEnergy].push_back(fitresult_centralbars->Parameter(1));
        meanTofCentralBarsErrorMap[beamEnergy].push_back(fitresult_centralbars->ParError(1));
        meanTofXYMap[beamEnergy].push_back(fitresult_xy->Parameter(1));  // Add results for bar 9 XY
        meanTofXYErrorMap[beamEnergy].push_back(fitresult_xy->ParError(1));

        sigmaTofMap[beamEnergy].push_back(fitresult->Parameter(2));
        sigmaTofErrorMap[beamEnergy].push_back(fitresult->ParError(2));
        sigmaTofCentralBarsMap[beamEnergy].push_back(fitresult_centralbars->Parameter(2));
        sigmaTofCentralBarsErrorMap[beamEnergy].push_back(fitresult_centralbars->ParError(2));
        sigmaTofXYMap[beamEnergy].push_back(fitresult_xy->Parameter(2));  // Add results for bar 9 XY
        sigmaTofXYErrorMap[beamEnergy].push_back(fitresult_xy->ParError(2));

        inFile->Close();
        delete inFile;
    }

    TGraphErrors* graphMeanTof = new TGraphErrors();
    TGraphErrors* graphMeanTofCentralbars = new TGraphErrors();
    TGraphErrors* graphMeanTofXY = new TGraphErrors(); // New graph for bar 9 XY

    TGraphErrors* graphSigmaTof = new TGraphErrors();
    TGraphErrors* graphSigmaTofCentralbars = new TGraphErrors();
    TGraphErrors* graphSigmaTofXY = new TGraphErrors(); // New graph for bar 9 XY

    int pointIndex = 0;
    for (double energy : beamEnergies) {
        auto meanTof = ComputeMeanAndError(meanTofMap[energy], meanTofErrorMap[energy]);
        auto meanTofCentralBars = ComputeMeanAndError(meanTofCentralBarsMap[energy], meanTofCentralBarsErrorMap[energy]);
        auto meanTofXY = ComputeMeanAndError(meanTofXYMap[energy], meanTofXYErrorMap[energy]); // Calculate for bar 9 XY

        auto sigmaTof = ComputeMeanAndError(sigmaTofMap[energy], sigmaTofErrorMap[energy]);
        auto sigmaTofCentralBars = ComputeMeanAndError(sigmaTofCentralBarsMap[energy], sigmaTofCentralBarsErrorMap[energy]);
        auto sigmaTofXY = ComputeMeanAndError(sigmaTofXYMap[energy], sigmaTofXYErrorMap[energy]); // Calculate for bar 9 XY

        graphMeanTof->SetPoint(pointIndex, energy, meanTof.first);
        graphMeanTof->SetPointError(pointIndex, 0, meanTof.second);
        graphMeanTofCentralbars->SetPoint(pointIndex, energy, meanTofCentralBars.first);
        graphMeanTofCentralbars->SetPointError(pointIndex, 0, meanTofCentralBars.second);
        graphMeanTofXY->SetPoint(pointIndex, energy, meanTofXY.first); // Add point for bar 9 XY
        graphMeanTofXY->SetPointError(pointIndex, 0, meanTofXY.second);

        graphSigmaTof->SetPoint(pointIndex, energy, sigmaTof.first);
        graphSigmaTof->SetPointError(pointIndex, 0, sigmaTof.second);
        graphSigmaTofCentralbars->SetPoint(pointIndex, energy, sigmaTofCentralBars.first);
        graphSigmaTofCentralbars->SetPointError(pointIndex, 0, sigmaTofCentralBars.second);
        graphSigmaTofXY->SetPoint(pointIndex, energy, sigmaTofXY.first); // Add point for bar 9 XY
        graphSigmaTofXY->SetPointError(pointIndex, 0, sigmaTofXY.second);

        pointIndex++;
    }

    PlotGraph(graphMeanTof, graphMeanTofCentralbars, graphMeanTofXY, "Mean_ToF", "Mean ToF vs Beam Energy", "Mean ToF", "Mean ToF Central Bars (8, 9, 10)", "Mean ToF Bar 9(X), Bar 9(Y)", 12.);
    PlotGraph(graphSigmaTof, graphSigmaTofCentralbars, graphSigmaTofXY, "Sigma_ToF", "#sigma ToF vs Beam Energy", "#sigma ToF", "#sigma ToF Central Bars (8, 9, 10)", "#sigma ToF Bar 9(X), Bar 9(Y)", 0.2);
}
