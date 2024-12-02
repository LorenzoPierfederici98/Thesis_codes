#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <iostream>
#include <sstream>
#endif

std::string ConvertFileNumbersToString(const std::vector<int>& energies);

pair<double, double> RoundMeasurement(double value, double uncertainty);