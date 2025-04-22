
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
#include <numeric>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TNamed.h>
#include <TSpectrum.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <regex>
#endif

void ProcessFile(const std::string& fileName, const int energy, std::map<int, std::vector<std::pair<int, double>>>& Res, std::map<int, std::vector<std::pair<int, double>>>& Res_err);

std::pair<double, double> RoundMeasurement(double value, double uncertainty);

void PlotResolutions(const std::map<int, std::vector<std::pair<int, double>>>& Res,
    const std::map<int, std::vector<std::pair<int, double>>>& Res_err);