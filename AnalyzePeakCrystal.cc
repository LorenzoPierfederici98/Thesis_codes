#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TLegend.h>
#include <TObjString.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <TSpectrum.h>


void PrintMeasurement(double value, double uncertainty) {
    // Step 1: Get the order of magnitude of the uncertainty
    int significantFigures = (int)std::ceil(-std::log10(uncertainty)) + 1;

    // Step 2: Calculate the rounding factor based on significant figures
    double roundingFactor = std::pow(10, significantFigures);

    // Step 3: Round the uncertainty and value accordingly
    double roundedUncertainty = std::round(uncertainty * roundingFactor) / roundingFactor;
    double roundedValue = std::round(value * roundingFactor) / roundingFactor;

    // Step 4: Print the results with appropriate precision
    std::cout << std::fixed << std::setprecision(significantFigures) 
              << "Peak position: " << roundedValue << " ± " << roundedUncertainty << std::endl;
}

void AnalyzePeakCrystal(const vector<int> &fileNumbers, const int crystal_ID) {
    TCanvas* c1 = new TCanvas("c1", "Fit Results", 800, 600);

    for (int runNumber : fileNumbers) {
        TString filename = Form("AnaFOOT_Merge_HIT2022_%d.root", runNumber);
        TFile* inFile = new TFile(filename);
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            continue;
        }

        TObjString* energyObj = (TObjString*)inFile->Get("BeamEnergyInfo");
        Double_t beamEnergy = energyObj->GetString().Atof();
        TH1D* Charge_Calo_crystal = (TH1D*)inFile->Get(Form("Charge_Calo_crystalId_%d", crystal_ID));

        if (!Charge_Calo_crystal || Charge_Calo_crystal->GetEntries() == 0) {
            std::cerr << "Warning: Histogram is empty for run " << runNumber << std::endl;
            continue;
        }

        // Draw the histogram first, before fitting
        c1->cd();
        gPad->SetLogy();
        double x_min = 0.4; // Set to your desired min x range
        double x_max = 0.6; // Set to your desired max x range
        Charge_Calo_crystal->GetXaxis()->SetRangeUser(x_min, x_max);
        Charge_Calo_crystal->Draw();  // This will ensure the histogram is drawn

        // Set initial parameters manually for Gaussian fit
        TF1 *gausFit = new TF1("gausFit", "gaus", 0.4, 0.6);
        // Definire altro istogramma zoomato per fare il fit.
        gausFit->SetParameters(100, 0.49, 0.01);  // Example parameters: (amplitude, mean, sigma)
        // Perform the fit and check the status
        TFitResultPtr c = Charge_Calo_crystal->Fit(gausFit, "QS", "", 0.47, 0.52);
        if (c->Status() != 0) {
            std::cerr << "Fit failed for run " << runNumber << std::endl;
        } else {
            // If fit is successful, print the results
            std::cout << "Run: " << runNumber << " Beam energy: " << beamEnergy << " MeV ";
            PrintMeasurement(c->Parameter(1), c->ParError(1));
        }

        // Update and save the canvas
        c1->Update();
        //TString outputImage = Form("fit_run%d_crystal%d.png", runNumber, crystal_ID);
        //c1->SaveAs(outputImage);
        // Pause for the user to view each plot (optional)
        //gPad->WaitPrimitive(); // Uncomment to wait for user input before moving to the next plot
        inFile->Close();
        delete inFile;
    }

    delete c1;
}
