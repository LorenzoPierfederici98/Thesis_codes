/**
 * @file CaloPeakEnergyDisplay.h
 * @brief Header file for the CaloPeakEnergyDisplay macro.
 * 
 * Macro that collects the charge fit results of AnalyzePeakCrystal.cc for every calorimeter crystal
 * from files named like AnaFOOT_Calo_Decoded_HIT2022_100MeV_Fit.root. A charge-beam energy plot is
 * built for every crystal and fitted with a 1 parameter (slope) linear function. An intercalibration
 * is then performed for those crystals with at least 2 fitted charge peaks out of 5 energy values, by
 * computing the ratio of the slopes of a certain crystal with respect to crystal ID 0 (the central crystal).
 * Those ratios are used as scaling factors to sum the charge for events with a single cluster of size 2,
 * that is with 2 crystals involved. The calibration coefficients for every available crystal are stored
 * in a file named SlopeRatios.cal, stored in calib/HIT2022/.
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <sstream>
#endif

/**
 * @brief Converts a vector of integers representing energies to a comma-separated string.
 * 
 * @param energies A vector of integers representing energies.
 * @return A string containing the comma-separated list of energies.
 */
std::string ConvertFileNumbersToString(const std::vector<int>& energies);

/**
 * @brief Rounds the given value and uncertainty to the appropriate number of significant figures based on the uncertainty.
 * 
 * @param value The value to be rounded.
 * @param uncertainty The uncertainty associated with the value.
 * @return A pair containing the rounded value and rounded uncertainty.
 */
std::pair<double, double> RoundMeasurement(double value, double uncertainty);