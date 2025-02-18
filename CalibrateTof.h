/**
 * @file CalibrateTof.h
 * @brief Header file for the ToF calibration.
 * 
 * This header file contains the declarations for the functions used in the ToF (Time of Flight) calibration.
 * The calibration is performed by computing the differences between the ToF from data and ToF from MC mean fit values
 * for every bar, as provided by the AnalyzeTofFragm.cc and AnalyzeTofMC.cc macros. These differences are written in
 * four calibration files, one for each energy, along with the corresponding uncertainties (squared sum of the values).
 * To be run with root -l -b -q 'CalibrateTof.cc()'.
 * 
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <tuple>
#include <map>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdio>  // for sscanf
#include "TString.h"
#endif

/**
 * @brief Loads a fit result from a specified file and directory.
 * 
 * @param fileName The name of the file containing the fit result.
 * @param fitResultName The name of the fit result to load.
 * @return A pointer to the loaded TFitResult, or nullptr if an error occurs.
 */
TFitResult* LoadFitResult(const std::string& fileName, const std::string& fitResultName);

/**
 * @brief Extracts mean values and uncertainties from fit results.
 * 
 * @param filesAndEnergies A vector of pairs containing file names and corresponding energies.
 * @param layer The layer identifier (e.g., "X" or "Y").
 * @param bar The bar number.
 * @param means A vector to store the extracted mean values.
 * @param meanErrors A vector to store the extracted mean errors.
 */
void ExtractMeanValues(const std::vector<std::pair<std::string, int>>& filesAndEnergies, const std::string& layer, int bar, std::vector<double>& means, std::vector<double>& meanErrors);

/**
 * @brief Writes the differences between data and MC mean values to calibration files.
 * 
 * This function writes the fit values to configuration files.
 * The fit values are sorted by bar ID and written to the file in a specific format.
 * 
 * @param filesAndEnergies A vector of pairs containing data file names and corresponding energies.
 * @param filesAndEnergiesMC A vector of pairs containing MC file names and corresponding energies.
 * @param layer The layer identifier (e.g., "X" or "Y").
 * @param bar The bar number.
 */
void WriteMeanDifferences(const std::vector<std::pair<std::string, int>>& filesAndEnergies, const std::vector<std::pair<std::string, int>>& filesAndEnergiesMC, const std::string& layer, int bar);



