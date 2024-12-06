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

int SumNentries(TFile* file);

void FitHistogramsInDirectory(
    TDirectory* dir, 
    double threshold,
    double energy, 
    std::map<TString, std::map<int, double>>& fitMeans, 
    std::map<TString, std::map<int, double>>& fitErrors, 
    TFile* outputFile
);

void ProcessFile(
    const std::string& fileName, 
    double energy, 
    std::map<TString, std::map<int, double>>& fitMeans, 
    std::map<TString, std::map<int, double>>& fitErrors
);

void CreateAndSaveGraph(
    const TString& layerBarCombination, 
    const std::map<int, double>& energiesAndFits, 
    const std::map<TString, std::map<int, double>>& fitErrors,
    std::map<int, double>& stopping_power
);