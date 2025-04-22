
#include "Z_resolutions.h"

void Z_resolutions()
{
    std::vector<std::pair<std::string, int>> filesAndEnergies = {
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_100MeV.root", 100},
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_140MeV.root", 140},
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_200MeV.root", 200},
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_220MeV.root", 220}
    };

    TString LatexFileName_Z = TString("ZResolutions.txt");
    std::ofstream LatexFile_Z(LatexFileName_Z.Data());
    OpenZFile(LatexFile_Z);

    std::vector<int> energies;
    std::vector<double> R_Z_vec, R_Z_err_vec, meanZ_vec, meanZ_err_vec;

    for (const auto& [fileName, energy] : filesAndEnergies) {
        energies.push_back(energy);
        ProcessFile(fileName, energy, R_Z_vec, R_Z_err_vec, meanZ_vec, meanZ_err_vec, LatexFile_Z);
    }

    CloseZFile(LatexFile_Z);
    PlotResolutionGraphs(R_Z_vec, R_Z_err_vec);
}

void ProcessFile(
        const std::string& fileName,
        int energy,
        std::vector<double>& R_Z_vec,
        std::vector<double>& R_Z_err_vec,
        std::vector<double>& meanZ_vec,
        std::vector<double>& meanZ_err_vec,
        std::ofstream& LatexFile_Z
)
{
    TFile* file = TFile::Open(fileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        return;
    }

    std::cout << "Processing file: " << std::endl;

    // Create an output ROOT file to store fitted histograms
    TString outputFileName = TString(fileName).ReplaceAll(".root", "_Fit.root");
    TFile* outputFile = TFile::Open(outputFileName, "UPDATE");
    outputFile->cd();

    FitHistograms(file, energy, R_Z_vec, R_Z_err_vec, meanZ_vec, meanZ_err_vec, LatexFile_Z);

    outputFile->Close();
    file->Close();
    delete outputFile;
    delete file;
}

void FitHistograms(
    TFile* inFile,
    int energy,
    std::vector<double>& R_Z_vec,
    std::vector<double>& R_Z_err_vec,
    std::vector<double>& meanZ_vec,
    std::vector<double>& meanZ_err_vec,
    std::ofstream& LatexFile_Z
)
{
    if (!inFile) return;

    TIter next(inFile->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        TH1D* hist = dynamic_cast<TH1D*>(obj);
        if (!hist) {
            delete obj;
            continue;
        }

        TString histName = hist->GetName();
        if (histName == "Z_reco") {       
            cout << endl;
            cout << "beam energy: " << energy << " MeV/u " << endl;
            hist->Draw();

            TFitResultPtr fitResult = FitWithTSpectrum(hist, energy);
            //FitStats statsZ = FitResampling(hist, energy, fitResult->Parameter(1), fitResult->Parameter(2));

            double meanZ = fitResult->Parameter(1);
            double meanErrZ = fitResult->Error(1);
            double stdZ = fitResult->Parameter(2);
            double stdErrZ = fitResult->Error(2);

            double R_Z;
            double R_Z_err;

            if (fitResult.Get() != nullptr) {
                double cov_mu_sigma_Z = fitResult->CovMatrix(1, 2);

                if (meanZ > 0 || stdZ / meanZ < 0.2 || meanZ / meanErrZ - 1 > 0.5) {
                    if (meanErrZ > 0.000045) meanErrZ = 0.0001;
                    if (stdErrZ > 0.000045) stdErrZ = 0.0001;
                    R_Z = stdZ / meanZ;
                    R_Z_err = (1. / meanZ) * sqrt(pow(stdErrZ, 2) + pow(R_Z * meanErrZ, 2));
                    
                    //if (R_Z_err >= 0.0005) R_Z_err = 0.001;
                    auto [R_Z_str, R_Z_err_str] = RoundMeasurement(R_Z, R_Z_err);
                    R_Z = std::stod(R_Z_str);
                    R_Z_err = std::stod(R_Z_err_str);

                    R_Z *= 100.;
                    R_Z_err *= 100.;

                    R_Z_vec.push_back(R_Z);
                    R_Z_err_vec.push_back(R_Z_err);

                    auto [meanZ_str, meanErrZ_str] = RoundMeasurement(meanZ, meanErrZ);
                    auto [stdZ_str, stdErrZ_str] = RoundMeasurement(stdZ, stdErrZ);
                    meanZ = std::stod(meanZ_str);
                    meanErrZ = std::stod(meanErrZ_str);
                    stdZ = std::stod(stdZ_str);
                    stdErrZ = std::stod(stdErrZ_str);

                    meanZ_vec.push_back(meanZ);
                    meanZ_err_vec.push_back(meanErrZ);

                    fitResult->Write(Form("FitResultZ"), TObject::kOverwrite);

                    cout << "Z resolution (%): " << R_Z << " +/- " << R_Z_err << endl;
                    cout << "cov_mu_sigma_Z: " << cov_mu_sigma_Z << endl;

                    WriteZTable(LatexFile_Z, energy,
                        meanZ, stdZ, R_Z,
                        meanErrZ, stdErrZ, R_Z_err);
                }

            }
            // Save tZ histogram with both fits to tZ output file
            hist->Write("", TObject::kOverwrite);
        }
        delete hist;
    }
}

TFitResultPtr FitWithTSpectrum(TH1D *hist, int energy)
{
    int nPeaks = 0;
    TSpectrum spectrum(1);
    nPeaks = spectrum.Search(hist, 1, "", 0.0015);
    cout << "# of peaks: " << nPeaks << endl;
    TFitResultPtr fitResult;

    if (nPeaks == 0)
    {
        std::cerr << "No TOF peaks found" << std::endl;
        return nullptr;
    }
    else
    {
        double peakPosition = spectrum.GetPositionX()[0];
        int binFit;
        if (energy == 100)
        {
            binFit = 10;
        }
        else if (energy == 140)
        {
            binFit = 12;
        }
        else
        {
            binFit = 15;
        }
        int binPeak = hist->FindBin(peakPosition);
        int binLow = std::max(1, binPeak - binFit);
        int binHigh = std::min(hist->GetNbinsX(), binPeak + binFit);
        TF1 *gaus = new TF1("gaus", "gaus", hist->GetBinLowEdge(binLow), hist->GetBinLowEdge(binHigh + 1));
        gaus->SetParameter(1, peakPosition);
        gaus->SetParLimits(1, peakPosition - 0.5, peakPosition + 0.5);
        fitResult = hist->Fit("gaus", "QRS");
        if (fitResult->Status() != 0)
        {
            cout << "fit failed" << endl;
            return nullptr;
        }
        else
        {
            hist->GetListOfFunctions()->Add(gaus);
            gaus->SetLineColor(kRed);
            gaus->Draw("same");
        }
    }
    return fitResult;
}

void WriteZTable(std::ofstream& outFile, int energy, 
    double meanZ, double stdZ, double R_Z,
    double meanErrZ, double stdErrZ, double R_Z_err)
{
    if (!outFile.is_open()) {
    std::cerr << "Error: Output file stream is not open." << std::endl;
    return;
    }

    // Write LaTeX table row with errors in \pm notation
    outFile << energy << " & "
    << meanZ << " $\\pm $ " << meanErrZ << " & "
    << stdZ << " $\\pm $ " << stdErrZ << " & "
    << R_Z << " $\\pm $ " << R_Z_err << " \\\\ \\midrule\n";
}

void OpenZFile(std::ofstream& outFile)
{
    // Write LaTeX table header
    outFile << "\\begin{table}[h]\n"
    << "    \\centering\n"
    << "    \\begin{tabular}{cccc}\n"
    << "        \\toprule\n"
    << "        Beam energy & $\\mu(Z)$ & $\\sigma(Z)$ & $R_{Z}$\\\\ \n"
    << "        {[}MeV/u{]} &  &  & {[}\\%{]}\\\\ \\midrule\n";
}

void CloseZFile(std::ofstream& outFile)
{
    //if (!outFile.is_open()) return;
    // Write LaTeX table footer
    outFile << "        \\bottomrule\n"
    << "    \\end{tabular}\n"
    << "    \\caption{Reconstructed Z resolution results at each beam energy value.}\n"
    << "    \\label{Tab:Z_resolutions}\n"
    << "\\end{table}\n";
    outFile.close();
    std::cout << "LaTeX table written successfully." << std::endl;
}

void PlotResolutionGraphs(
    const std::vector<double>& R_Z_vec,
    const std::vector<double>& R_Z_err_vec
) {
    gStyle->SetTitleSize(0.07, "T");

    // ----------- Canvas 1: R_Z -----------
    auto cZ = new TCanvas("cZ", "Z Resolution vs energy loss", 800, 600);
    cZ->SetGrid();
    cZ->SetMargin(0.15, 0.1, 0.15, 0.1);

    std::vector<double> meanEloss_vec = {9.0379, 7.1191, 5.6532, 5.2572};
    std::vector<double> meanElossErr_vec = {0.0004, 0.0003, 0.0004, 0.0003};

    auto gMeanZ = new TGraphErrors(meanEloss_vec.size(),
                                        meanEloss_vec.data(),
                                        R_Z_vec.data(),
                                        meanElossErr_vec.data(),
                                        R_Z_err_vec.data());

    gMeanZ->SetTitle("Reconstructed Z resolutions;#mu(#Delta E);Resolution [%]");
    gMeanZ->SetMarkerStyle(20);
    gMeanZ->SetMarkerSize(1.5);
    gMeanZ->SetLineWidth(2);
    gMeanZ->GetXaxis()->SetTitleSize(0.05);
    gMeanZ->GetYaxis()->SetTitleSize(0.05);
    gMeanZ->Draw("AP");

    cZ->SaveAs("Plots/Z_Resolution.png");
}

pair<std::string, std::string> RoundMeasurement(double value, double uncertainty) {
    // Determine the order of magnitude of the uncertainty
    int uncertaintyOrder = (int)std::floor(std::log10(uncertainty));
    double roundingFactor = std::pow(10, uncertaintyOrder);

    // Round uncertainty to 1 significant figure
    double roundedUncertainty = std::round(uncertainty / roundingFactor) * roundingFactor;

    // Adjust the value to match the precision of the uncertainty
    double roundedValue = std::round(value / roundingFactor) * roundingFactor;

    // Determine the number of decimal places to display based on the rounded uncertainty
    int decimalPlaces = std::max(0, -uncertaintyOrder);

    // Format value and uncertainty into strings with the required precision
    std::ostringstream valueStream, uncertaintyStream;
    valueStream << std::fixed << std::setprecision(decimalPlaces) << roundedValue;
    uncertaintyStream << std::fixed << std::setprecision(decimalPlaces) << roundedUncertainty;

    return {valueStream.str(), uncertaintyStream.str()};
}