// Macro that plots the energy loss and ToF calibration coefficients for each bar of the 2 layers
// To be run with root -l -b -q 'CalibCoefficients.cc()'
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <iostream>
#include <TF1.h>
#include <TLegend.h>

void PlotElossCoefficients(const std::vector<int>& barID, const std::vector<double>& elossCoeff_error_layerY, const std::vector<double>& elossCoeff_error_layerX);

void PlotTofCoefficients(const std::vector<int>& barID, const std::vector<int>& energies, const std::map<int, std::string>& calibTofFiles);

void CalibCoefficients()
{
    std::vector<int> energies = {100, 140, 200, 220};
    std::vector<int> barID;
    for (int i = 0; i < 20; i++)
    {
        barID.push_back(i);
    }
    std::map<int, std::string> calibTofFiles;
    calibTofFiles[100] = "4766";
    calibTofFiles[140] = "4801";
    calibTofFiles[200] = "4742";
    calibTofFiles[220] = "4828";
    std::vector<double> elossCoeff_error_layerY = {0.002, 0.001, 0.004, 0.003, 0.003, 0.002, 0.0009, 0.0005, 0.00009, 0.00003, 0.00005, 0.0004, 0.001, 0.001, 0.002, 0.003, 0.001, 0.004, 0.005, 0.002};
    std::vector<double> elossCoeff_error_layerX = {0.006, 0.001, 0.004, 0.002, 0.002, 0.001, 0.0008, 0.0004, 0.00003, 0.00003, 0.0001, 0.0005, 0.001, 0.001, 0.002, 0.003, 0.001, 0.003, 0.003, 0.003};

    PlotElossCoefficients(barID, elossCoeff_error_layerY, elossCoeff_error_layerX);
    PlotTofCoefficients(barID, energies, calibTofFiles);
}

void PlotElossCoefficients(const std::vector<int>& barID, 
    const std::vector<double>& elossCoeff_error_layerY, 
    const std::vector<double>& elossCoeff_error_layerX)
{
std::ifstream file("TATW_Energy_Calibration_perBar_4742.cal");
std::vector<double> elossCoeff_layerY, elossCoeff_layerX;

std::string line;
while (std::getline(file, line)) {
if (line.empty() || line[0] == '#') continue;

std::istringstream iss(line);
int barId;
double p0_inverse, p1, layer;

if (!(iss >> barId >> p0_inverse >> p1 >> layer)) continue;

if (layer == 0) {
elossCoeff_layerY.push_back(p0_inverse);
} else {
elossCoeff_layerX.push_back(p0_inverse);
}
}

// Convert bar IDs to doubles for plotting
std::vector<double> x(barID.begin(), barID.end());

TGraphErrors* grY = new TGraphErrors(x.size(), x.data(), elossCoeff_layerY.data(), nullptr, elossCoeff_error_layerY.data());
TGraphErrors* grX = new TGraphErrors(x.size(), x.data(), elossCoeff_layerX.data(), nullptr, elossCoeff_error_layerX.data());

grY->SetMarkerStyle(32);
grY->SetMarkerSize(1.5);
grY->SetMarkerColor(kRed);
grY->SetLineColor(kRed);
grY->SetLineWidth(2);
grY->SetTitle("Charge-energy loss calibration coefficients");

grX->SetMarkerStyle(4);
grX->SetMarkerSize(1.5);
grX->SetLineWidth(2);
grX->SetMarkerColor(kBlack);
grX->SetLineColor(kBlack);

TCanvas* c = new TCanvas("c", "Eloss Coefficients", 800, 600);
c->SetMargin(0.15, 0.12, 0.15, 0.15);
c->SetGrid();
grY->GetXaxis()->SetNdivisions(520);
grY->GetXaxis()->SetTitle("Bar ID");
grY->GetYaxis()->SetTitle("1/p0 [a.u. / MeV]");
grY->GetXaxis()->SetTitleSize(0.05);
grY->GetYaxis()->SetTitleSize(0.05);
gStyle->SetTitleSize(0.07, "T");
grY->Draw("APL");
grX->Draw("PL SAME");

grY->GetXaxis()->SetLimits(-0.5, 20.);
grY->SetMinimum(0.55);

TLegend* legend = new TLegend(0.5, 0.3, 0.65, 0.45);
legend->AddEntry(grY, "Layer Y", "p");
legend->AdcalculateddEntry(grX, "Layer X", "p");
legend->Draw();

c->SaveAs("Plots/ElossCoefficients.png");  // Optional: save plot
}

void PlotTofCoefficients(const std::vector<int>& barID, const std::vector<int>& energies, const std::map<int, std::string>& calibTofFiles)
{
    for (int energy : energies)
    {
        std::string fileName = "TATW_Tof_Calibration_perBar_" + calibTofFiles.at(energy) + ".cal";
        std::ifstream file(fileName);

        std::vector<double> tofCoeff_layerY, tofCoeff_layerX;
        std::vector<double> tofCoeff_error_layerY, tofCoeff_error_layerX;
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            int barId;
            double Delta_t, sigma_t, layer;

            if (!(iss >> barId >> Delta_t >> sigma_t >> layer)) continue;

            if (layer == 0) {
                tofCoeff_layerY.push_back(Delta_t);
                tofCoeff_error_layerY.push_back(sigma_t);
            } else {
                tofCoeff_layerX.push_back(Delta_t);
                tofCoeff_error_layerX.push_back(sigma_t);
            }
        }

        // Convert bar IDs to doubles for plotting
        std::vector<double> x(barID.begin(), barID.end());

        TGraphErrors* grY = new TGraphErrors(x.size(), x.data(), tofCoeff_layerY.data(), nullptr, tofCoeff_error_layerY.data());
        TGraphErrors* grX = new TGraphErrors(x.size(), x.data(), tofCoeff_layerX.data(), nullptr, tofCoeff_error_layerX.data());

        grY->SetMarkerStyle(32);
        grY->SetMarkerSize(1.5);
        grY->SetMarkerColor(kRed);
        grY->SetLineColor(kRed);
        grY->SetLineWidth(2);
        grY->SetTitle(Form("TOF data - TOF MC calibration coefficients @ %d MeV/u beam energy", energy));

        grX->SetMarkerStyle(4);
        grX->SetMarkerSize(1.5);
        grX->SetLineWidth(2);
        grX->SetMarkerColor(kBlack);
        grX->SetLineColor(kBlack);

        TCanvas* c = new TCanvas("c", "Tof Coefficients", 800, 600);
        c->SetMargin(0.15, 0.12, 0.15, 0.15);
        c->SetGrid();
        grY->GetXaxis()->SetNdivisions(520);
        grY->GetXaxis()->SetTitle("Bar ID");
        grY->GetYaxis()->SetTitle("#Delta TOF [ns]");
        grY->GetXaxis()->SetTitleSize(0.05);
        grY->GetYaxis()->SetTitleSize(0.05);
        gStyle->SetTitleSize(0.07, "T");
        grY->Draw("APL");
        grX->Draw("PL SAME");
        grY->GetXaxis()->SetLimits(-0.5, 20.);
        grY->SetMinimum(1.8);

        TLegend* legend = new TLegend(0.5, 0.7, 0.65, 0.8);
        legend->AddEntry(grY, "Layer Y", "p");
        legend->AddEntry(grX, "Layer X", "p");
        legend->SetTextSize(0.04);
        legend->Draw();
        c->SaveAs(Form("Plots/TofCoefficients_%d.png", energy));  // Optional: save plot
    }
}