//Macro that computes the sum and fits the 1D charge histograms for every crystal of the calorimeter.
//The histograms are retireved from the AnaLyzeCalo.cc root output files and are summed.
//The user provides the run numbers (all of them must have the same energy), the macro then sums
//all the histograms of every crystalID and of the total calorimeter of every run; the fit on a summed
//histogram for a certain crystalID is performed only if its entries are grater than thresh=0.002 times
//the entries of the total summed histogram of the calorimeter.

//To be run with e.g.  root -l -q 'AnalyzePeakCrystal.cc({4742, 4743, 4743, 4744, 4745, 4828}, x_min, x_max)',
//all the runs having the same energy, x_min and x_max define the range surrounding the peak.
//The fit results are inserted in files namekd like e.g. Fit_Calo_Crystal_1_Energy_200MeV.root, the 1 being the
//crystalID.

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

void PrintMeasurement(double value, double uncertainty);

TFitResultPtr FitPeakWithTSpectrum(TH1D *hist, double threshold);

Double_t RetrieveEnergy(int runNumber, TFile* inFile);

std::tuple<TH1D*, TH1D*, Double_t> SumHistograms(const std::vector<int>& fileNumbers, const TString& histName_total, const TString& histName);

std::string ConvertFileNumbersToString(const std::vector<int>& fileNumbers);

void SaveFitResultsToFile(TCanvas* canvas, TH1D* hist, TFitResultPtr fitResult, const TString& outputFileName, const TString& fileNumbersStr);

void AnalyzePeakCrystal(const vector<int> &fileNumbers, const double x_min, const double x_max) {

    std::string fileNumbersStr = ConvertFileNumbersToString(fileNumbers);

    for (int crystal_ID = 1; crystal_ID < 63; crystal_ID++){
        TCanvas* c1 = new TCanvas("c1", Form("Fit Results"), 800, 600);
        const TString histName = Form("Charge_Calo_crystalId_%d", crystal_ID);
        const TString histName_total = Form("Charge_Calo");
        auto [Charge_Calo_crystal, Charge_Calo, beamEnergy] = SumHistograms(fileNumbers, histName_total, histName);

        if (!Charge_Calo_crystal || Charge_Calo_crystal->GetEntries() == 0 || !Charge_Calo || Charge_Calo->GetEntries() == 0) {
            std::cerr << "Warning: Sum Histogram is empty" << std::endl;
            continue;
        }
        if (beamEnergy == -1){
            std::cerr << "Invalid energy value" << std::endl;
            continue;
        }
        cout << "crystal ID: " << crystal_ID << " sum hist entries / total entries: " << Charge_Calo_crystal->GetEntries()/Charge_Calo->GetEntries() << endl;
        double thresh = 0.002;
        if (Charge_Calo_crystal->GetEntries() > thresh * Charge_Calo->GetEntries()){
            c1->SetTitle(Form("Run Numbers: %s | Beam Energy: %.0f MeV | Crystal ID: %d", fileNumbersStr.c_str(), beamEnergy, crystal_ID));

            TH1D* Charge_Calo_fullrange = (TH1D*)Charge_Calo_crystal->Clone("Charge_Calo_fullrange"); // Clone the full-range histogram
            if (crystal_ID == 4){
                Charge_Calo_crystal->GetXaxis()->SetRangeUser(0.2, 0.4);
            }
            else {
                Charge_Calo_crystal->GetXaxis()->SetRangeUser(x_min, x_max);
            }

            TFitResultPtr c = FitPeakWithTSpectrum(Charge_Calo_crystal, 0.3);
            if (!c.Get()) {
                std::cerr << "Fit failed for run " << std::endl;
                continue;
            } else {
                // If fit is successful, print the results
                std::cout << "Beam energy: " << beamEnergy << " MeV ";
                PrintMeasurement(c->Parameter(1), c->ParError(1));
            }

            // Update and save the canvas
            //c1->Update();
            //TString outputImage = Form("fit_run%d_crystal%d.png", runNumber, crystal_ID);
            //c1->SaveAs(outputImage);
            // Pause for the user to view each plot (optional)
            //gPad->WaitPrimitive(); // Uncomment to wait for user input before moving to the next plot

            const TString outFile = Form("FitCalo/Fit_Calo_Crystal_%d_Energy_%.0fMeV.root", crystal_ID, beamEnergy);
            SaveFitResultsToFile(c1, Charge_Calo_fullrange, c, outFile, fileNumbersStr);
        }
        delete c1;
    }
    return;
}

void PrintMeasurement(double value, double uncertainty) {
    //Get the order of magnitude of the uncertainty
    int significantFigures = (int)std::ceil(-std::log10(uncertainty)) + 1;

    //Calculate the rounding factor based on significant figures
    double roundingFactor = std::pow(10, significantFigures);

    //Round the uncertainty and value accordingly
    double roundedUncertainty = std::round(uncertainty * roundingFactor) / roundingFactor;
    double roundedValue = std::round(value * roundingFactor) / roundingFactor;

    //Print the results with appropriate precision
    std::cout << std::fixed << std::setprecision(significantFigures) 
              << "Peak position: " << roundedValue << " ± " << roundedUncertainty << std::endl;
}

TFitResultPtr FitPeakWithTSpectrum(TH1D *hist, double threshold) {
    //Use TSpectrum to find peaks beyond the threshold
    TSpectrum *spectrum = new TSpectrum(1);  // 1 = maximum number of peaks to search for
    int nPeaks = spectrum->Search(hist, 2, "", threshold); // 2 = sigma for smoothing, threshold to ignore noise
    if (nPeaks == 0) {
        std::cerr << "No peak found beyond threshold " << threshold << std::endl;
        return TFitResultPtr(0);
    }

    // Step 2: Get the peak position from TSpectrum
    double *xPeaks = spectrum->GetPositionX(); // Array of peak positions in x
    double peakPosition = xPeaks[0]; // We assume we're interested in the first peak

    std::cout << "Peak found at x = " << peakPosition << std::endl;

    //Determine the bin where the peak is located
    int binMax = hist->FindBin(peakPosition);

    //Define the fit range (10-bin window around the peak)
    int bin_range = 10;
    int binLow = binMax - bin_range; // 10 bins before the peak
    int binHigh = binMax + bin_range; // 10 bins after the peak

    if (binLow < 1) binLow = 1; // Ensure the range stays within valid bin numbers
    if (binHigh > hist->GetNbinsX()) binHigh = hist->GetNbinsX();

    double xLow = hist->GetBinLowEdge(binLow);
    double xHigh = hist->GetBinLowEdge(binHigh + 1);

    std::cout << "Fitting in the range: [" << xLow << ", " << xHigh << "]" << " (" << bin_range << " bins) " << "around the peak" << std::endl;

    //Perform the Gaussian fit in the determined range
    TFitResultPtr fitResult = hist->Fit("gaus", "S", "", xLow, xHigh);

    if (fitResult->IsValid()) {
        std::cout << "Fit successful" << endl;
    } else {
        std::cout << "Fit did not converge properly." << std::endl;
    }

    delete spectrum;  // Clean up
    return fitResult;
}

Double_t RetrieveEnergy(int runNumber, TFile* inFile) {
    if (!inFile) {
        std::cerr << "Error: Input file is null for run " << runNumber << std::endl;
        return -1;  // Return a default invalid value
    }

    TObject* obj = inFile->Get("BeamEnergyInfo");
    
    if (!obj) {
        std::cerr << "Error: Could not retrieve 'BeamEnergyInfo' from file for run " << runNumber << std::endl;
        return -1;  // Return a default invalid value
    }

    // Check if the retrieved object is a TObjString
    if (obj->InheritsFrom(TObjString::Class())) {
        TObjString* energyObj = (TObjString*)obj;
        Double_t beamEnergy = energyObj->GetString().Atof();
        return beamEnergy;
    } else {
        std::cerr << "Error: 'BeamEnergyInfo' is not a TObjString for run " << runNumber << std::endl;
        return -1;  // Return a default invalid value
    }
}

std::tuple<TH1D*, TH1D*, Double_t> SumHistograms(const std::vector<int>& fileNumbers, const TString& histName_total, const TString& histName) {
    TH1D* summedHist = nullptr;
    TH1D* summedHist_total = nullptr;
    Double_t beamEnergy = 0;  // Store the energy value
    
    // Loop over all file numbers
    for (int runNumber : fileNumbers) {
        TString filename = Form("Calo/AnaFOOT_Calo_Decoded_HIT2022_%d.root", runNumber);
        TFile* inFile = TFile::Open(filename);
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            continue;
        }

        if (beamEnergy == 0) {  // Extract beam energy only from the first valid file
            beamEnergy = RetrieveEnergy(runNumber, inFile);
        }

        TH1D* hist = (TH1D*)inFile->Get(histName);
        TH1D* hist_total = (TH1D*)inFile->Get(histName_total);

        if (!hist || hist->GetEntries() == 0) {
            std::cerr << "Error: Could not find histogram or it is empty " << histName << std::endl;
            inFile->Close();
            continue;
        }

        if (!hist_total || hist_total->GetEntries() == 0) {
            std::cerr << "Error: Could not find total crystals histogram or it is empty " << histName_total << std::endl;
            inFile->Close();
            continue;
        }

        if (!summedHist) {
            summedHist = (TH1D*)hist->Clone();
            summedHist->SetDirectory(0);
        } else {
            summedHist->Add(hist);
        }

        if (!summedHist_total) {
            summedHist_total = (TH1D*)hist_total->Clone();
            summedHist_total->SetDirectory(0);
        } else {
            summedHist_total->Add(hist_total);
        }

        inFile->Close();
    }

    return {summedHist, summedHist_total, beamEnergy};
}


std::string ConvertFileNumbersToString(const std::vector<int>& fileNumbers) {
    std::stringstream ss;
    for (size_t i = 0; i < fileNumbers.size(); ++i) {
        ss << fileNumbers[i];
        if (i != fileNumbers.size() - 1) {
            ss << ", ";  // Add a comma and space between numbers
        }
    }
    return ss.str();
}

void SaveFitResultsToFile(TCanvas* canvas, TH1D* hist, TFitResultPtr fitResult, const TString& outputFileName, const TString& fileNumbersStr) {
    // Create a new ROOT file using TString's Data() method
    TFile* outFile = new TFile(outputFileName.Data(), "RECREATE");
    
    // Check if file creation was successful
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Could not create output file " << outputFileName << std::endl;
        return;
    }

    canvas->cd();
    gPad->SetLogy();
    canvas->Update();
    hist->Write(Form("Full range histogram runs: %s", fileNumbersStr.Data()));
    canvas->Write(outputFileName);         // Save the histogram

    // Write the run numbers as a TObjString
    TObjString runNumbersObj(fileNumbersStr);
    runNumbersObj.Write("RunNumbers");  // Save the run numbers with a label "RunNumbers"
    fitResult->Write(Form("Fit_Results"));

    // Close the ROOT file to ensure everything is saved properly
    outFile->Close();
    
    // Clean up
    delete outFile;

    std::cout << "Fit results, histogram, and run numbers saved to " << outputFileName << std::endl;
}






