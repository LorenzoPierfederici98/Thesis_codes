// Macro that fits the charge distribution from the AnalyzeTWChargeTime.cc merged output files (for MC runs).
// A fit is performed with 2 separate gaussians, one for proton and one for helium peaks. The peaks are automatically
// found with TSPectrum in a certain range (between 1. and 12., the bins outside the range are set to 0); the peaks are
// then fitted within a certain bin-range centered around the peak. The fit results are stored in files name like e.g.
// TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_140_Fit.root. To be run with root -l -b -q 'AnalyzeTWMC.cc()'.

#include "AnalyzeTWMC.h"

void AnalyzeTWMC() {

    std::vector<std::pair<std::string, int>> filesAndEnergies = {
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_100.root", 100},
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_140.root", 140},
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_200.root", 200},
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_220.root", 220}
    };

    std::map<int, double> fitMeanP;
    std::map<int, double> fitErrorP;
    std::map<int, double> fitMeanHe;
    std::map<int, double> fitErrorHe;

    for (const auto& [fileName, energy] : filesAndEnergies) {
        ProcessFile(fileName, energy, fitMeanP, fitMeanHe, fitErrorP, fitErrorHe);
    }

}

void ProcessFile(
    const std::string& fileName, 
    int energy, 
    std::map<int, double>& fitMeanP,
    std::map<int, double>& fitMeanHe,
    std::map<int, double>& fitErrorP,
    std::map<int, double>& fitErrorHe
) {
    TFile* file = TFile::Open(fileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        return;
    }
    // Create an output ROOT file to store fitted histograms
    TString outputFileName = TString(fileName).ReplaceAll(".root", "_Fit.root");
    TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
    std::cout << "Processing file: " << fileName << std::endl;
    
    TH1D* hist = dynamic_cast<TH1D*>(file->Get("Eloss_true_ch"));
    if (!hist) {
        std::cerr << "Error: Histogram 'Eloss_true_ch' not found in file " << fileName << std::endl;
        file->Close();
        delete file;
        outputFile->Close();
        delete outputFile;
        return;
    }
    outputFile->cd();
    hist->Draw();

    std::pair<TFitResultPtr, TFitResultPtr> fitResults = FitPeaksWithTSpectrum(hist, energy);
    TFitResultPtr fitResult1 = fitResults.first, fitResult2 = fitResults.second;

    if (fitResult1.Get() != nullptr){
        double meanChargeP = fitResult1->Parameter(1);
        double stdChargeP = fitResult1->Parameter(2);
        double meanChargeErrP = fitResult1->Error(1);
        if (meanChargeP > 0 || stdChargeP / meanChargeP < 0.2 || meanChargeP / meanChargeErrP - 1 > 0.5) {
            fitMeanHe[energy] = meanChargeP;
            fitErrorHe[energy] = meanChargeErrP;
            fitResult1->Write(Form("FitResultP_MC_%dMeV", energy));
        }
    }

    if (fitResult2.Get() != nullptr) {
        double meanChargeHe = fitResult2->Parameter(1);
        double stdChargeHe = fitResult2->Parameter(2);
        double meanChargeErrHe = fitResult2->Error(1);
        if (meanChargeHe > 0 || stdChargeHe / meanChargeHe < 0.2 || meanChargeHe / meanChargeErrHe - 1 > 0.5) {
            fitMeanP[energy] = meanChargeHe;
            fitErrorP[energy] = meanChargeErrHe;
            fitResult2->Write(Form("FitResultHe_MC_%dMeV", energy));
        }
    }
    // Save the histogram with both fits to the output file
    hist->Write();
    outputFile->Close();
    file->Close();
    delete outputFile;
    delete file;
}

std::pair<TFitResultPtr, TFitResultPtr> FitPeaksWithTSpectrum(TH1D *hist, int energy) {

    int nPeaks = 0;
    double thresh_peak_low = 1.;
    double thresh_peak_high = 12.;

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
        cout << "first peak found at x=" << x1 << endl;
        int binMax1 = hist->FindBin(x1);
        int bins_fit_p;
        if (energy == 100){
            bins_fit_p = 3;
        }
        else {
            bins_fit_p = 2;
        }
        int binLow1 = std::max(1, binMax1 - bins_fit_p);
        int binHigh1 = std::min(hist->GetNbinsX(), binMax1 + bins_fit_p);
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
        cout << "second peak found at x=" << x2 << endl;
        int binMax2 = hist->FindBin(x2);
        int binLow2 = std::max(1, binMax2 - 3);
        int binHigh2 = std::min(hist->GetNbinsX(), binMax2 + 3);
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