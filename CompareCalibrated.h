#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include <iostream>
#include <vector>
#include <utility>
#include <string>
#endif

TH1D* getElossByBar(TFile* file, const std::string &dirName, const TString &desiredHistName);

TH1D* getTofByBar(TFile* file, const std::string &dirName, const TString &desiredHistName);

void ElossNormalizedHistograms(TH1* hData, TH1* hMC, const TString &ElossTitle, const TString &ElossCanvasName);

void TofNormalizedHistograms(TH1* hTofData, TH1* hTofMC, const TString &TofTitle, const TString &TofCanvasName);

void ProcessFile(const std::string &dataFile, const std::string &mcFile, int energy, const std::string &layer, int bar);