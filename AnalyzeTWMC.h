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

void ProcessFile(
    const std::string& fileName, 
    int energy, 
    std::map<int, double>& fitMeanP,
    std::map<int, double>& fitMeanHe,
    std::map<int, double>& fitErrorP,
    std::map<int, double>& fitErrorHe
);

pair<std::string, std::string> RoundMeasurement(double value, double uncertainty);

std::pair<TFitResultPtr, TFitResultPtr> FitPeaksWithTSpectrum(TH1D *hist, int energy);