/**
 * @file AnalyzeTWMC.h
 * @brief Header file for processing and fitting charge distribution histograms from ROOT files.
 * 
 * Macro that fits the charge distribution from the AnalyzeTWChargeTime.cc merged output files (for MC runs).
 * A fit is performed with 2 separate gaussians, one for proton and one for helium peaks. The peaks are automatically
 * found with TSpectrum in a certain range (between 1. and 12., the bins outside the range are set to 0); the peaks are
 * then fitted within a certain bin-range centered around the peak. The fit results are stored in files name like e.g.
 * TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_140_Fit.root. To be run with root -l -b -q 'AnalyzeTWMC.cc()'.
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TNamed.h>
#include <TSpectrum.h>
#include <iostream>
#include <vector>
#include <map>
#include <regex>
#endif

/**
 * @brief Processes a ROOT file to fit the charge distribution.
 * 
 * This function processes a ROOT file containing charge distribution data.
 * It performs a fit with two separate Gaussians, one for proton and one for helium peaks.
 * The peaks are automatically found with TSpectrum in a certain range (between 1. and 12.),
 * and then fitted within a certain bin-range centered around the peak.
 * The fit results are stored in files named like e.g. TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_140_Fit.root.
 * 
 * @param fileName The name of the ROOT file to process.
 * @param energy The energy value associated with the file.
 * @param fitMeanP A map to store the mean fit values for protons.
 * @param fitMeanHe A map to store the mean fit values for helium.
 * @param fitErrorP A map to store the fit errors for protons.
 * @param fitErrorHe A map to store the fit errors for helium.
 */
void ProcessFile(
    const std::string& fileName, 
    int energy, 
    std::map<int, double>& fitMeanP,
    std::map<int, double>& fitMeanHe,
    std::map<int, double>& fitErrorP,
    std::map<int, double>& fitErrorHe
);

/**
 * @brief Fits the peaks in a histogram using TSpectrum.
 * 
 * This function uses TSpectrum to search for peaks in a histogram within a specified range.
 * It then fits the found peaks with Gaussian functions. The first peak is assumed to be for protons,
 * and the second peak (if found) is assumed to be for heliums.
 * 
 * @param hist The histogram to fit.
 * @param energy The energy value associated with the histogram.
 * @return A pair of TFitResultPtr, where the first element is the fit result for protons and the second element is the fit result for helium.
 */
std::pair<TFitResultPtr, TFitResultPtr> FitPeaksWithTSpectrum(TH1D *hist, int energy);