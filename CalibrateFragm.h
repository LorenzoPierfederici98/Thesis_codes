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

void ProcessFile(const TString& fileName, 
                 int energy, 
                 std::map<TString, std::map<int, double>>& fitMeansP, 
                 std::map<TString, std::map<int, double>>& fitMeansHe, 
                 std::map<TString, std::map<int, double>>& fitErrorsP, 
                 std::map<TString, std::map<int, double>>& fitErrorsHe);

void PlotFitResultsCombined(const std::map<TString, std::map<int, double>>& fitMeansP,
                            const std::map<TString, std::map<int, double>>& fitErrorsP,
                            const std::map<TString, std::map<int, double>>& fitMeansHe,
                            const std::map<TString, std::map<int, double>>& fitErrorsHe,
                            std::map<int, double>elossP,
                            std::map<int, double>elossHe);

pair<std::string, std::string> RoundMeasurement(double value, double uncertainty);