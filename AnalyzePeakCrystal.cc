// Macro that fits the 1D charge histograms for a given crystal of the calorimeter.
// The histograms are retireved form the AnaLyzeCalo.cc merged output files (the histograms are automatically summed).
// To be run with e.g.  root -l -b -q 'AnalyzePeakCrystal.cc(x_min, x_max, fit_thresh)', -b doesn't display plots. 
// In an interval between x_min and x_max a peak beyond fit_thresh is automatically found with TSpectrum and fitted.
// The fit results are inserted in files namekd like e.g. Calo/AnaFOOT_Calo_Decoded_HIT2022_140MeV_Fit.root which also
// store the crystalID charge histogram and the fit plot restricted in [x_min, x_max].

#include "AnalyzePeakCrystal.h"

void AnalyzePeakCrystal() {
    std::vector<std::pair<std::string, double>> filesAndEnergies = {
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_100MeV.root", 100},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_140MeV.root", 140},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_180MeV.root", 180},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_200MeV.root", 200},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_220MeV.root", 220}
    };

    const double x_min = 0.2;
    const double x_max = 0.6;

    for (const auto &[fileName, energy] : filesAndEnergies) {
        TFile *inFile = TFile::Open(fileName.c_str());
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error: Cannot open file " << fileName << std::endl;
            continue;
        }

        for (int crystal_ID = 1; crystal_ID < 63; crystal_ID++) {
            TCanvas *canvas = new TCanvas(Form("c1_%d", crystal_ID), "Fit Results", 800, 600);
            TString histName = Form("Charge_Calo_crystalId_%d", crystal_ID);
            TString histName_total = "Charge_Calo";

            auto [Charge_Calo_crystal, Charge_Calo] = FindHistograms(inFile, histName_total, histName);

            if (!Charge_Calo_crystal || Charge_Calo_crystal->GetEntries() == 0 || !Charge_Calo || Charge_Calo->GetEntries() == 0) {
                std::cerr << "Warning: Histogram for crystal " << crystal_ID << " is empty or missing." << std::endl;
                delete canvas;
                continue;
            }
            double threshold = 0.002;
            if (Charge_Calo_crystal->GetEntries() > threshold * Charge_Calo->GetEntries()) {
                canvas->SetTitle(Form("Beam Energy: %.0f MeV | Crystal ID: %d", energy, crystal_ID));

                TH1D *Charge_Calo_fullrange = static_cast<TH1D *>(Charge_Calo_crystal->Clone());
                Charge_Calo_crystal->GetXaxis()->SetRangeUser(
                    (energy == 100) ? 0.1 : x_min,
                    (energy == 100) ? 0.35 : x_max
                );

                const double fit_thresh = (energy == 100) ? 0.1 : 0.25;

                TFitResultPtr fitResult = FitPeakWithTSpectrum(Charge_Calo_crystal, fit_thresh);
                if (fitResult.Get() != nullptr) {
                    std::cout << "Beam energy: " << energy << " MeV " << " crystal: " << crystal_ID << " ";
                    double meanCharge = fitResult->Parameter(1);
                    double meanChargeErr = fitResult->ParError(1);
                    double stdCharge = fitResult->Parameter(2);
                    PrintMeasurement(meanCharge, meanChargeErr);
                    if (meanCharge / stdCharge - 1 < 0.1 || stdCharge == 0. || meanCharge / meanChargeErr - 1 < 0.5) {
                        std::cout << "Skipping fit because the mean/std or mean/meanErr charge ratios are close to 1" << std::endl;
                        continue;
                    }
                } else {
                    std::cerr << "Fit failed for crystal " << crystal_ID << std::endl;
                    continue;
                }

                TString outputFileName = fileName;
                outputFileName.ReplaceAll(".root", "_Fit.root");
                SaveFitResultsToFile(canvas, Charge_Calo_fullrange, fitResult, outputFileName);

                delete Charge_Calo_fullrange;
            }

            delete canvas;
        }
        inFile->Close();
        delete inFile;
    }
}


void PrintMeasurement(double value, double uncertainty) {
    int sigFigs = static_cast<int>(std::ceil(-std::log10(uncertainty))) + 1;
    double roundFactor = std::pow(10, sigFigs);
    std::cout << std::fixed << std::setprecision(sigFigs) 
              << "Peak position: " << std::round(value * roundFactor) / roundFactor 
              << " ± " << std::round(uncertainty * roundFactor) / roundFactor << std::endl;
}


TFitResultPtr FitPeakWithTSpectrum(TH1D *hist, double threshold) {
    TSpectrum spectrum(1);
    int nPeaks = spectrum.Search(hist, 2, "", threshold);

    if (nPeaks == 0) {
        std::cerr << "No peak found beyond threshold " << threshold << std::endl;
        return TFitResultPtr(nullptr);
    }

    double peakPos = spectrum.GetPositionX()[0];
    std::cout << "Peak found at x = " << peakPos << std::endl;

    int binMax = hist->FindBin(peakPos);
    int binLow = std::max(1, binMax - 10);
    int binHigh = std::min(hist->GetNbinsX(), binMax + 10);

    std::cout << "Fitting range: [" << hist->GetBinLowEdge(binLow) << ", " << hist->GetBinLowEdge(binHigh + 1) << "]" << std::endl;
    return hist->Fit("gaus", "S", "", hist->GetBinLowEdge(binLow), hist->GetBinLowEdge(binHigh + 1));
}


std::tuple<TH1D *, TH1D *> FindHistograms(TFile *inFile, const TString &histName_total, const TString &histName) {
    if (!inFile || inFile->IsZombie()) return {nullptr, nullptr};

    TH1D *hist = dynamic_cast<TH1D *>(inFile->Get(histName));
    TH1D *hist_total = dynamic_cast<TH1D *>(inFile->Get(histName_total));

    return {hist, hist_total};
}

void SaveFitResultsToFile(TCanvas *canvas, TH1D *hist, TFitResultPtr fitResult, const TString &outputFileName) {
    TFile *outFile = TFile::Open(outputFileName, "UPDATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << outputFileName << std::endl;
        return;
    }

    canvas->cd();
    gPad->SetLogy();
    canvas->Update();
    hist->Write(hist->GetName(), TObject::kOverwrite);
    canvas->Write(Form("Canvas_Crystal_%s", hist->GetName()),  TObject::kOverwrite);
    fitResult->Write(Form("Fit_result_%s", hist->GetName()),  TObject::kOverwrite);
    outFile->Close();
    std::cout << "Saved fit results to " << outputFileName << std::endl;
}