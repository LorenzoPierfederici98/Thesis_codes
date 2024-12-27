// Macro that fits the charge distribution from the AnalyzeTWChargeTime.cc merged output files (for fragmentation runs).
// A fit is performed with 2 separate gaussians, one for proton and one for helium peaks. The fit limits for the proton
// peaks depend on the beam energy. Only the histograms whose entries are greater than a fraction of the sum of the merged files
// total event number are fitted. The fit results are stored in files name like e.g.
// TW/AnaFOOT_TW_Decoded_HIT2022_140MeV_Fit.root. To be run with root -l -b -q 'AnalyzeTWFragm.cc()'.

#include "AnalyzeTWFragm.h"

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
        TH1D* hist = dynamic_cast<TH1D*>(obj);
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
        hist->Draw();
        // Fit histogram if above the threshold
        if (hist->GetEntries() > threshold) {
            // the bins of the original histograms before thresh_peak_low beyond thresh_peak_high
            // are set to 0 and 2 peaks are searched for inbetween
            double thresh_peak_high;
            double thresh_peak_low = 1.;
            if (energy == 100) {
                thresh_peak_high = 15.;
            }
            else if (energy == 140) {
                thresh_peak_high = 12.;
            }
            else if (energy == 200) {
                thresh_peak_high = 9.;              
            }
            else if (energy == 220) {
                thresh_peak_high = 8.5;
                thresh_peak_low = 0.85;
            }
            cout << endl;
            cout << "beam energy: " << energy << " MeV/u " << layerBarCombination << endl;
            std::pair<TFitResultPtr, TFitResultPtr> fitResults = FitPeaksWithTSpectrum(hist, energy, thresh_peak_low, thresh_peak_high);
            TFitResultPtr fitResult1 = fitResults.first, fitResult2 = fitResults.second;

            if (fitResult1.Get() != nullptr){
                double meanChargeP = fitResult1->Parameter(1);
                double stdChargeP = fitResult1->Parameter(2);
                double meanChargeErrP = fitResult1->Error(1);
                if (meanChargeP > 0 || stdChargeP / meanChargeP < 0.2 || meanChargeP / meanChargeErrP - 1 > 0.5) {
                    fitMeansHe[layerBarCombination][energy] = meanChargeP;
                    fitErrorsHe[layerBarCombination][energy] = meanChargeErrP;
                    fitResult1->Write(Form("FitResultP_%s", layerBarCombination.Data()));
                }
            }

            if (fitResult2.Get() != nullptr) {
                double meanChargeHe = fitResult2->Parameter(1);
                double stdChargeHe = fitResult2->Parameter(2);
                double meanChargeErrHe = fitResult2->Error(1);
                if (meanChargeHe > 0 || stdChargeHe / meanChargeHe < 0.2 || meanChargeHe / meanChargeErrHe - 1 > 0.5) {
                    fitMeansP[layerBarCombination][energy] = meanChargeHe;
                    fitErrorsP[layerBarCombination][energy] = meanChargeErrHe;
                    fitResult2->Write(Form("FitResultHe_%s", layerBarCombination.Data()));
                }
            }

            // Save the histogram with both fits to the output file
            hist->Write();
        }

        delete hist;
    }
}

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


std::pair<TFitResultPtr, TFitResultPtr> FitPeaksWithTSpectrum(TH1D *hist, double energy, double thresh_peak_low, double thresh_peak_high) {

    int nPeaks = 0;

    int binLow = hist->FindBin(thresh_peak_low); // Find the bin corresponding to x = low_thresh
    int binHigh = hist->FindBin(thresh_peak_high);
    int binMax = hist->GetNbinsX(); // Upper limit of the histogram

    // Clone the histogram and set the desired range
    TH1D* histRestricted = (TH1D*)hist->Clone("histRestricted");
    for (int bin = 1; bin < binLow; ++bin) {
        histRestricted->SetBinContent(bin, 0); // Zero out bins below x = thresh
    }
    for (int bin = binHigh; bin < binMax; ++bin) {
        histRestricted->SetBinContent(bin, 0);
    }
    // Use TSpectrum to search for peaks
    TSpectrum spectrum(2); // Max number of peaks to find
    nPeaks = spectrum.Search(histRestricted, 2, "", 0.0015); // Use the restricted histogram
    delete histRestricted;

    cout << "# of peaks: " << nPeaks << endl;
    TFitResultPtr fitResult1, fitResult2;

    double *peakPositions = spectrum.GetPositionX();

    // Sort peaks: the first belongs to protons, the second to heliums
    std::vector<double> sortedPeaks(peakPositions, peakPositions + nPeaks);
    std::sort(sortedPeaks.begin(), sortedPeaks.end());

    if (nPeaks == 0) {
        std::cerr << "No peaks found beyond threshold " << 1. << std::endl;
        return{nullptr, nullptr};
    }

    if (nPeaks >= 1) {
        double x1 = sortedPeaks[0];
        if (x1 < 1.2084 && energy == 100) { // fix LayerX_bar1 at 100 MeV/u
            x1 = 3.268;
        }
        cout << "first peak found at x=" << x1 << endl;
        int binMax1 = hist->FindBin(x1);
        int binLow1 = std::max(1, binMax1 - 4);
        int binHigh1 = std::min(hist->GetNbinsX(), binMax1 + 4);
        TF1* gaus1 = new TF1("gaus1", "gaus", hist->GetBinLowEdge(binLow1), hist->GetBinLowEdge(binHigh1 + 1));
        gaus1->SetParameter(1, x1);
        gaus1->SetParLimits(1, x1 - 0.1, x1 + 0.1);
        fitResult1 = hist->Fit("gaus1", "RS");
        if (fitResult1->Status() != 0){
            cout << "fit 1 failed" << endl;
            fitResult1 = nullptr;
        }
        else {
            hist->GetListOfFunctions()->Add(gaus1);
            gaus1->SetLineColor(kBlue);
            gaus1->Draw("same");
        }
    }

    if (nPeaks == 2) {
        double x2 = sortedPeaks[1];
        if (x2 - sortedPeaks[0] < 2.5)
            x2 = sortedPeaks[0]*4;
        if (x2 < 4.84 && energy == 100) {  // fix LayerX_bar1 at 100 MeV/u
            x2 = 12.5884;
        }
        cout << "second peak found at x=" << x2 << endl;
        int binMax2 = hist->FindBin(x2);
        int binLow2 = std::max(1, binMax2 - 6);
        int binHigh2 = std::min(hist->GetNbinsX(), binMax2 + 6);
        TF1* gaus2 = new TF1("gaus2", "gaus", hist->GetBinLowEdge(binLow2), hist->GetBinLowEdge(binHigh2 + 1));
        gaus2->SetParameter(1, x2);
        gaus2->SetParLimits(1, x2 - 0.1, x2 + 0.1);
        fitResult2 = hist->Fit("gaus2", "RS+");
        if (fitResult2->Status() != 0){
            cout << "fit 2 failed" << endl;
            fitResult2 = nullptr;
        }
        else {
            hist->GetListOfFunctions()->Add(gaus2);
            gaus2->SetLineColor(kRed);
            gaus2->Draw("same");
        }
    }

    return {fitResult1, fitResult2};
}