#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>         // For working with ROOT files
#include <TKey.h>          // For iterating over objects in ROOT files
#include <TDirectory.h>    // For handling directories in ROOT files
#include <TH1.h>           // For handling histograms
#include <TF1.h>           // For fitting functions
#include <TGraphErrors.h>  // For plotting graphs with error bars
#include <TCanvas.h>       // For canvas drawing
#include <TString.h>       // For working with ROOT-specific strings
#include <TIterator.h>     // For iterating over lists in ROOT
#include <TLegend.h>       // For adding legends to graphs
#include <iostream>        // For standard input/output
#include <vector>          // For using std::vector containers
#include <map>             // For using std::map containers
#include <utility>         // For std::pair (energy, fit result)
#include <regex>           // For regular expressions (if you're filtering object names)
#include <cmath>           // Required for sqrt
#endif

struct FitResult {
    double mean;
    double meanError;
    double sigma;
    double sigmaError;
};

FitResult SaveFitAndExtractParams(TFile* file, TH1* hist, TF1* fitFunc);

void ProcessAndPlot(const std::string& fileName, double energy,
                    std::map<TString, std::vector<std::pair<double, FitResult>>>& results);

void PlotResults(const std::map<TString, std::vector<std::pair<double, FitResult>>>& results);


