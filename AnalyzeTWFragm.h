/**
 * @file AnalyzeTWFragm.h
 * @brief Header file for AnalyzeTWFragm.cc
 * 
 * This file contains function declarations for fitting the charge distribution from the AnalyzeTWChargeTime.cc merged output files
 * (for  fragmentation runs).
 * A fit is performed with 2 separate gaussians, one for proton and one for helium peaks. The fit limits for the proton
 * peaks depend on the beam energy. Only the histograms whose entries are greater than a fraction of the sum of the merged files
 * total event number are fitted (the entries of each file part of the merged output are saved in it and summed).
 * The histogram peaks are automatically found with TSPectrum in a certain range (energy-dependent, the bins outside the range
 * are set to 0 and the peaks are searched for inbetween); the fits are then performed within a certain bin-range centered around
 * the peak, depending on the specific bar and beam energy. The fit results are stored in files named like e.g.
 * TW/AnaFOOT_TW_Decoded_HIT2022_140MeV_Fit.root. To be run with root -l -b -q 'AnalyzeTWFragm.cc()'.
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
 * @brief Sums the number of entries in a ROOT file.
 * 
 * This function iterates over the keys in the provided ROOT file and sums the
 * entries of objects whose names match the pattern "nentries_run_\\d+".
 * 
 * @param file Pointer to the ROOT file.
 * @return Total number of entries.
 */
int SumNentries(TFile* file);

/**
 * @brief Fits histograms in a directory and stores the results.
 * 
 * This function fits the histograms in the specified directory using two separate
 * Gaussian functions, one for proton peaks and one for helium peaks. The fit limits
 * depend on the beam energy. Only histograms with entries greater than a specified
 * threshold are fitted. The fit results are stored in the provided maps.
 * 
 * @param dir Pointer to the directory containing the histograms.
 * @param threshold Minimum number of entries required to fit a histogram.
 * @param energy Beam energy.
 * @param fitMeansP Map to store the mean values of the proton fits.
 * @param fitErrorsP Map to store the errors of the proton fits.
 * @param fitMeansHe Map to store the mean values of the helium fits.
 * @param fitErrorsHe Map to store the errors of the helium fits.
 * @param outputFile Pointer to the output ROOT file to store the fitted histograms.
 */
void FitHistogramsInDirectory(
    TDirectory* dir, 
    double threshold, 
    double energy, 
    std::map<TString, std::map<int, double>>& fitMeansP, 
    std::map<TString, std::map<int, double>>& fitErrorsP, 
    std::map<TString, std::map<int, double>>& fitMeansHe, 
    std::map<TString, std::map<int, double>>& fitErrorsHe, 
    TFile* outputFile
);

/**
 * @brief Processes a ROOT file and fits the histograms.
 * 
 * This function processes the specified ROOT file, sums the number of entries,
 * calculates the threshold, and fits the histograms in the specified directories.
 * The fit results are stored in the provided maps.
 * 
 * @param fileName Name of the ROOT file.
 * @param energy Beam energy.
 * @param fitMeansP Map to store the mean values of the proton fits.
 * @param fitMeansHe Map to store the mean values of the helium fits.
 * @param fitErrorsP Map to store the errors of the proton fits.
 * @param fitErrorsHe Map to store the errors of the helium fits.
 */
void ProcessFile(
    const std::string& fileName, 
    double energy, 
    std::map<TString, std::map<int, double>>& fitMeansP,
    std::map<TString, std::map<int, double>>& fitMeansHe,
    std::map<TString, std::map<int, double>>& fitErrorsP,
    std::map<TString, std::map<int, double>>& fitErrorsHe
);

/**
 * @brief Fits peaks in a histogram using TSpectrum.
 * 
 * This function fits two peaks in the specified histogram using TSpectrum. The peaks
 * are fitted with Gaussian functions, and the fit results are returned.
 * 
 * @param hist Pointer to the histogram.
 * @param energy Beam energy.
 * @param thresh_peak_low Lower threshold for peak search.
 * @param thresh_peak_high Upper threshold for peak search.
 * @param layerBarCombination String representing the layer and bar combination.
 * @return Pair of fit results for the proton and helium peaks.
 */
std::pair<TFitResultPtr, TFitResultPtr> FitPeaksWithTSpectrum(
    TH1D *hist, 
    double energy, 
    double thresh_peak_low, 
    double thresh_peak_high, 
    const TString& layerBarCombination
);