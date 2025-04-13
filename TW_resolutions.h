
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
 
 void FitHistograms(
     TFile* inFile,
     int energy,
     std::vector<double>& R_He_vec,
    std::vector<double>& R_He_err_vec,
    std::vector<double>& meanEloss_vec,
    std::vector<double>& meanEloss_err_vec,
    std::vector<double>& R_Tof_vec,
    std::vector<double>& R_Tof_err_vec,
    std::vector<double>& meanTof_vec,
    std::vector<double>& meanTof_err_vec,
     std::ofstream& LatexFile_Eloss,
     std::ofstream& LatexFile_Tof
 );
 
 void ProcessFile(
     const std::string& fileName, 
     int energy,
     std::vector<double>& R_He_vec,
    std::vector<double>& R_He_err_vec,
    std::vector<double>& meanEloss_vec,
    std::vector<double>& meanEloss_err_vec,
    std::vector<double>& R_Tof_vec,
    std::vector<double>& R_Tof_err_vec,
    std::vector<double>& meanTof_vec,
    std::vector<double>& meanTof_err_vec,
     std::ofstream& LatexFile_Eloss,
     std::ofstream& LatexFile_Tof
 );

 TFitResultPtr FitWithTSpectrum(TH1D *hist, int energy);

 void WriteElossTable(std::ofstream& outFile, int energy, 
    double meanElossHe, double stdElossHe, double R_He,
    double meanElossErrHe, double stdElossErrHe, double R_He_err);

void WriteTofTable(std::ofstream& outFile, int energy, 
    double meanTof, double meanTofErr,
    double R_Tof, double R_Tof_err);

void OpenElossFile(std::ofstream& outFile);

void OpenTofFile(std::ofstream& outFile);

void CloseElossFile(std::ofstream& outFile);

void CloseTofFile(std::ofstream& outFile);

void PlotResolutionGraphs(
    const std::vector<double>& R_He_vec,
    const std::vector<double>& R_He_err_vec,
    const std::vector<double>& meanEloss_vec,
    const std::vector<double>& meanEloss_err_vec,
    const std::vector<double>& R_Tof_vec,
    const std::vector<double>& R_Tof_err_vec,
    const std::vector<double>& meanTof_vec,
    const std::vector<double>& meanTof_err_vec
);

pair<std::string, std::string> RoundMeasurement(double value, double uncertainty);

