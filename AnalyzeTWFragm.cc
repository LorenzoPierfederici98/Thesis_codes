// Macro that fits the charge distribution from the AnalyzeTWChargeTime.cc merged output files (for fragmentation runs).
// A fit is performed with 2 separate gaussians, one for proton and one for helium peaks. The fit limits for the proton
// peaks depend on the beam energy. Only the histograms whose entries are greater than a fraction of the sum of the merged files
// total event number are fitted. The fit results are stored in files name like e.g.
// TW/AnaFOOT_TW_Decoded_HIT2022_140MeV_Fit.root. To be run with root -l -b -q 'AnalyzeTWFragm.cc()'.

#include "AnalyzeTWFragm.h"

// Main analysis function
void AnalyzeTWFragm() {
    std::vector<std::pair<std::string, double>> filesAndEnergies = {
        {"TW/AnaFOOT_TW_Decoded_HIT2022_fragm_100MeV.root", 100},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_fragm_140MeV.root", 140},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_fragm_200MeV.root", 200},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_fragm_220MeV.root", 220}
    };

    std::map<int, double>elossP = {{100, 2.721}, {140, 2.2020}, {200, 1.4761}, {220, 1.4267}};  // energy loss for protons, from MC
    std::map<int, double>elossHe = {{100, 9.7638}, {140, 7.2929}, {200, 5.4135}, {220, 5.2352}};  // energy loss for heliums, from MC

    std::map<TString, std::map<int, double>> fitMeansP;  // protons
    std::map<TString, std::map<int, double>> fitErrorsP;
    std::map<TString, std::map<int, double>> fitMeansHe;  // heliums
    std::map<TString, std::map<int, double>> fitErrorsHe;

    for (const auto& [fileName, energy] : filesAndEnergies) {
        ProcessFile(fileName, energy, fitMeansP, fitMeansHe, fitErrorsP, fitErrorsHe);
    }
}

// Function to sum all "nentries" values in a file
int SumNentries(TFile* file) {
    int totalNentries = 0;
    TIter next(file->GetListOfKeys());
    TKey* key;
    std::regex nentriesRegex("nentries_run_\\d+");

    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        TNamed* nentriesObj = dynamic_cast<TNamed*>(obj);
        if (nentriesObj && std::regex_match(nentriesObj->GetName(), nentriesRegex)) {
            totalNentries += std::stoi(nentriesObj->GetTitle());
        }
        delete obj;
    }
    return totalNentries;
}

// Function to fit histograms in a directory and store results
void FitHistogramsInDirectory(
    TDirectory* dir, 
    double threshold, 
    double energy, 
    std::map<TString, std::map<int, double>>& fitMeansP, 
    std::map<TString, std::map<int, double>>& fitErrorsP, 
    std::map<TString, std::map<int, double>>& fitMeansHe, 
    std::map<TString, std::map<int, double>>& fitErrorsHe, 
    TFile* outputFile
) {
    if (!dir) return;

    TIter next(dir->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        TH1* hist = dynamic_cast<TH1*>(obj);
        if (!hist) {
            delete obj;
            continue;
        }

        TString histName = hist->GetName();
        if (!histName.BeginsWith("Charge_Layer")) {
            delete hist;
            continue;
        }

        // Extract layer and bar identifiers
        TString layerBarCombination = histName(7, histName.Length() - 7);  // e.g. LayerY_bar9
        outputFile->cd();
        // Fit histogram if above the threshold
        if (hist->GetEntries() > threshold) {
            double Proton_xlimit1;
            double Proton_xlimit2;
            if (energy == 100) {
                Proton_xlimit1 = 0.45;
                Proton_xlimit2 = 4.3;
            }
            else if (energy == 140) {
                Proton_xlimit1 = 0.48;
                Proton_xlimit2 = 3.5;
            }
            else if (energy == 200 || energy == 220) {
                Proton_xlimit1 = 0.25;
                Proton_xlimit2 = 2.5;               
            }

            TF1* fitFunc1 = new TF1("fitFunc1", "gaus", Proton_xlimit1, Proton_xlimit2);
            TFitResultPtr fitResult1 = hist->Fit(fitFunc1, "QSR");
            //hist->Draw("E");
            if (fitResult1->Status() == 0) {  // Check first fit success
                double meanChargeP = fitFunc1->GetParameter(1);
                double stdChargeP = fitFunc1->GetParameter(2);
                double meanChargeErrP = fitFunc1->GetParError(1);
                if (meanChargeP > 0 || stdChargeP / meanChargeP < 0.2 || meanChargeP / meanChargeErrP - 1 > 0.5) {
                    fitMeansP[layerBarCombination][energy] = meanChargeP;
                    fitErrorsP[layerBarCombination][energy] = meanChargeErrP;
                    fitResult1->Write(Form("FitResultP_%s", layerBarCombination.Data()));
                }
            }

            TF1* fitFunc2 = new TF1("fitFunc2", "gaus", 3., 13.);
            TFitResultPtr fitResult2 = hist->Fit(fitFunc2, "QSR+");

            if (fitResult2->Status() == 0) {  // Check second fit success
                double meanChargeHe = fitFunc2->GetParameter(1);
                double stdChargeHe = fitFunc2->GetParameter(2);
                double meanChargeErrHe = fitFunc2->GetParError(1);
                if (meanChargeHe > 0 || stdChargeHe / meanChargeHe < 0.2 || meanChargeHe / meanChargeErrHe - 1 > 0.5) {
                    fitMeansHe[layerBarCombination][energy] = meanChargeHe;
                    fitErrorsHe[layerBarCombination][energy] = meanChargeErrHe;
                    fitResult2->Write(Form("FitResultHe_%s", layerBarCombination.Data()));
                }
            }

            // Save the histogram with both fits to the output file
            hist->Write();
            delete fitFunc1;
            delete fitFunc2;
        }

        delete hist;
    }
}

// Function to process a single file and extract fit data
void ProcessFile(
    const std::string& fileName, 
    double energy, 
    std::map<TString, std::map<int, double>>& fitMeansP,
    std::map<TString, std::map<int, double>>& fitMeansHe,
    std::map<TString, std::map<int, double>>& fitErrorsP,
    std::map<TString, std::map<int, double>>& fitErrorsHe
) {
    TFile* file = TFile::Open(fileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        return;
    }

    int totalNentries = SumNentries(file);
    double threshold = 0.009 * totalNentries;

    std::cout << "Processing file: " << fileName << ", Total nentries: " << totalNentries << ", Threshold: " << threshold << std::endl;

    // Create an output ROOT file to store fitted histograms
    TString outputFileName = TString(fileName).ReplaceAll(".root", "_Fit.root");
    TFile* outputFile = TFile::Open(outputFileName, "RECREATE");

    const std::vector<std::string> directories = {"ChargeTimeLayerX", "ChargeTimeLayerY"};
    for (const auto& dirName : directories) {
        TDirectory* dir = dynamic_cast<TDirectory*>(file->Get(dirName.c_str()));
        if (!dir) {
            std::cerr << "Error: Directory " << dirName << " not found in file " << fileName << std::endl;
            continue;
        }
        FitHistogramsInDirectory(dir, threshold, energy, fitMeansP, fitMeansHe, fitErrorsP, fitErrorsHe, outputFile);
    }

    outputFile->Close();
    file->Close();
    delete outputFile;
    delete file;
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
