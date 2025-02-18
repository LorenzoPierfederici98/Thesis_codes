/**
 * @file CalibrateFragm.h
 * @brief Header file for the charge-energy calibration.
 * 
 * Macro that plots the fit results on protons and helium peaks given by the AnalyzeTWFragm.cc macro vs the energies
 * loss values retrieved from MC (2 energy loss values per bar and beam energy, given by the AnalyzeTWMC.cc macro).
 * The fit charge values vs energy loss values are fitted with a 1 parameter linear function; these parameters
 * (one per bar) are written in the configuration file.
 * To be run with root -l -b -q 'CalibrateFragm.cc()'.
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
#include <iostream>
#include <vector>
#include <map>
#include <regex>
#endif

/**
 * @brief Processes a ROOT file to extract fit results for protons and helium.
 * 
 * This function processes a ROOT file containing fit results for protons and helium.
 * It extracts the mean and error values for each fit and stores them in the provided maps.
 * 
 * @param fileName The name of the ROOT file to process.
 * @param energy The energy value associated with the file.
 * @param fitMeansP A map to store the mean fit values for protons.
 * @param fitMeansHe A map to store the mean fit values for helium.
 * @param fitErrorsP A map to store the fit errors for protons.
 * @param fitErrorsHe A map to store the fit errors for helium.
 */
void ProcessFile(const TString& fileName, 
                 int energy, 
                 std::map<TString, std::map<int, double>>& fitMeansP, 
                 std::map<TString, std::map<int, double>>& fitMeansHe, 
                 std::map<TString, std::map<int, double>>& fitErrorsP, 
                 std::map<TString, std::map<int, double>>& fitErrorsHe);

/**
 * @brief Plots the combined fit results for protons and helium.
 * 
 * This function plots the combined fit results for protons and helium.
 * It creates graphs for the fit results and performs a combined fit with a linear function.
 * The fit values are stored in a map and returned.
 * 
 * @param fitMeansP A map containing the mean fit values for protons.
 * @param fitErrorsP A map containing the fit errors for protons.
 * @param fitMeansHe A map containing the mean fit values for helium.
 * @param fitErrorsHe A map containing the fit errors for helium.
 * @param elossP A map containing the energy loss values for protons.
 * @param elossHe A map containing the energy loss values for helium.
 * @return A map containing the fit values for each layer-bar combination.
 */
std::map<TString, double> PlotFitResultsCombined(
    const std::map<TString, std::map<int, double>>& fitMeansP,
    const std::map<TString, std::map<int, double>>& fitErrorsP,
    const std::map<TString, std::map<int, double>>& fitMeansHe,
    const std::map<TString, std::map<int, double>>& fitErrorsHe,
    const std::map<int, double>& elossP,
    const std::map<int, double>& elossHe);

/**
 * @brief Writes the fit values to a configuration file.
 * 
 * This function writes the fit values to a configuration file.
 * The fit values are sorted by bar ID and written to the file in a specific format.
 * 
 * @param fitValues A map containing the fit values for each layer-bar combination.
 */
void WriteFitValuesOrdered(const std::map<TString, double>& fitValues);

/**
 * @brief Rounds a measurement value and its uncertainty.
 * 
 * This function rounds a measurement value and its uncertainty to the appropriate number of significant figures.
 * The uncertainty is rounded to 1 significant figure.
 * 
 * @param value The measurement value.
 * @param uncertainty The uncertainty of the measurement.
 * @return A pair of strings containing the rounded value and uncertainty.
 */
pair<std::string, std::string> RoundMeasurement(double value, double uncertainty);