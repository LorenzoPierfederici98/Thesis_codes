#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TKey.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TH2.h>
#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <string>
#endif

// A small struct to hold fit results
struct FitResult {
    bool success;
    double p0;     // Intercept
    double p1;     // Slope
    double p0err;  // Error on intercept
    double p1err;  // Error on slope
};

std::pair<FitResult, TGraph*> FitScatterPlot(TGraph *graph, int i, int j, double maxX, double maxY);

TH2* Rebin2DHistogram(const TH2* h2, int newXbins, int newYbins);

void ProcessFile(const std::string &fileName, int energy);



