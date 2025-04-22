#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
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

void ProcessFile(
    const std::string& fileName,
    int energy,
    std::vector<double>& R_Z_vec,
    std::vector<double>& R_Z_err_vec,
    std::vector<double>& meanZ_vec,
    std::vector<double>& meanZ_err_vec,
    std::ofstream& LatexFile_Z
);

void FitHistograms(
    TFile* inFile,
    int energy,
    std::vector<double>& R_Z_vec,
    std::vector<double>& R_Z_err_vec,
    std::vector<double>& meanZ_vec,
    std::vector<double>& meanZ_err_vec,
    std::ofstream& LatexFile_Z
);

TFitResultPtr FitWithTSpectrum(TH1D *hist, int energy);

void WriteZTable(std::ofstream& outFile, int energy, 
    double meanZ, double stdZ, double R_Z,
    double meanErrZ, double stdErrZ, double R_Z_err);

void OpenZFile(std::ofstream& outFile);

void CloseZFile(std::ofstream& outFile);

void PlotResolutionGraphs(
    const std::vector<double>& R_Z_vec,
    const std::vector<double>& R_Z_err_vec
);

pair<std::string, std::string> RoundMeasurement(double value, double uncertainty);

