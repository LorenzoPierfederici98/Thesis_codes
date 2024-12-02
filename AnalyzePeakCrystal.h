#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TLegend.h>
#include <TObjString.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <TSpectrum.h>
#endif

void PrintMeasurement(double value, double uncertainty);

TFitResultPtr FitPeakWithTSpectrum(TH1D *hist, double threshold);

std::tuple<TH1D*, TH1D*> FindHistograms(TFile *inFile, const TString &histName_total, const TString &histName);

void SaveFitResultsToFile(TCanvas* canvas, TH1D* hist, TFitResultPtr fitResult, const TString& outputFileName);