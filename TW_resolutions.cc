// Macro that computes the energy loss (sigma/mu) and TOF resolutions (sigma)
// from the calibrated energy loss histogram and calibrated TOF histogram for
// all the TW bars.

#include "TW_resolutions.h"

void TW_resolutions() {
    std::vector<std::pair<std::string, int>> filesAndEnergies = {
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_100MeV.root", 100},
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_140MeV.root", 140},
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_200MeV.root", 200},
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_220MeV.root", 220}
    };

    TString LatexFileName_Eloss = TString("ElossResolutions.txt");
    TString LatexFileName_Tof = TString("TofResolutions.txt");

    std::ofstream LatexFile_Eloss(LatexFileName_Eloss.Data());
    std::ofstream LatexFile_Tof(LatexFileName_Tof.Data());

    OpenElossFile(LatexFile_Eloss);
    OpenTofFile(LatexFile_Tof);

    std::vector<int> energies;
    std::vector<double> R_He_vec, R_He_err_vec, meanEloss_vec, meanEloss_err_vec;
    std::vector<double> R_Tof_vec, R_Tof_err_vec, meanTof_vec, meanTof_err_vec;

    for (const auto& [fileName, energy] : filesAndEnergies) {
        energies.push_back(energy);
        ProcessFile(fileName, energy, R_He_vec, R_He_err_vec, meanEloss_vec, meanEloss_err_vec, R_Tof_vec, R_Tof_err_vec, meanTof_vec, meanTof_err_vec, LatexFile_Eloss, LatexFile_Tof);
    }

    CloseElossFile(LatexFile_Eloss);
    CloseTofFile(LatexFile_Tof);

    PlotResolutionGraphs(R_He_vec, R_He_err_vec, meanEloss_vec, meanEloss_err_vec, R_Tof_vec, R_Tof_err_vec, meanTof_vec, meanTof_err_vec);

}

// Function to fit histograms in a directory and store results
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
) {
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
        if (histName == "Eloss") {       
            cout << endl;
            cout << "beam energy: " << energy << " MeV/u " << endl;
            hist->Draw();

            TFitResultPtr fitResult = FitWithTSpectrum(hist, energy);

            double meanElossHe;
            double meanElossErrHe;
            double stdElossHe;
            double stdElossErrHe;
            double R_He;
            double R_He_err;

            if (fitResult.Get() != nullptr) {
                meanElossHe = fitResult->Parameter(1);
                meanElossErrHe = fitResult->Error(1);
                stdElossHe = fitResult->Parameter(2);
                stdElossErrHe = fitResult->Error(2);

                double cov_mu_sigma_He = fitResult->CovMatrix(1, 2);
                if (meanElossHe > 0 || stdElossHe / meanElossHe < 0.2 || meanElossHe / meanElossErrHe - 1 > 0.5) {
                    R_He = stdElossHe / meanElossHe;
                    R_He_err = (1. / meanElossHe) * sqrt(pow(stdElossErrHe, 2) + pow(R_He * meanElossErrHe, 2) - 2 * R_He * cov_mu_sigma_He);

                    auto [R_He_str, R_He_err_str] = RoundMeasurement(R_He, R_He_err);
                    R_He = std::stod(R_He_str);
                    R_He_err = std::stod(R_He_err_str);

                    R_He *= 100.;
                    R_He_err *= 100.;

                    R_He_vec.push_back(R_He);
                    R_He_err_vec.push_back(R_He_err);

                    auto [meanElossHe_str, meanElossErrHe_str] = RoundMeasurement(meanElossHe, meanElossErrHe);
                    auto [stdElossHe_str, stdElossErrHe_str] = RoundMeasurement(stdElossHe, stdElossErrHe);
                    meanElossHe = std::stod(meanElossHe_str);
                    meanElossErrHe = std::stod(meanElossErrHe_str);
                    stdElossHe = std::stod(stdElossHe_str);
                    stdElossErrHe = std::stod(stdElossErrHe_str);

                    meanEloss_vec.push_back(meanElossHe);
                    meanEloss_err_vec.push_back(meanElossErrHe);

                    fitResult->Write(Form("FitResultHe"), TObject::kOverwrite);

                    cout << "heliums resolution (%): " << R_He << " +/- " << R_He_err << endl;
                    cout << "cov_mu_sigma_He: " << cov_mu_sigma_He << endl;

                    WriteElossTable(LatexFile_Eloss, energy,
                        meanElossHe, stdElossHe, R_He,
                        meanElossErrHe, stdElossErrHe, R_He_err);
                }

            }
            // Save the histogram with both fits to the output file
            hist->Write("", TObject::kOverwrite);
        }
        else if (histName == "ToF") {
            cout << endl;
            cout << "beam energy: " << energy << " MeV/u " << endl;
            hist->Draw();

            TFitResultPtr fitResult = FitWithTSpectrum(hist, energy);
            if (fitResult.Get() != nullptr) {
                double meanTof = fitResult->Parameter(1);
                double meanTofErr = fitResult->Error(1);
                double stdTof = fitResult->Parameter(2);
                double stdTofErr = fitResult->Error(2);
                double cov_mu_sigma_ToF = fitResult->CovMatrix(1, 2);

                auto [meanTof_str, meanTofErr_str] = RoundMeasurement(meanTof, meanTofErr);
                auto [stdTof_str, stdTofErr_str] = RoundMeasurement(stdTof, stdTofErr);
                meanTof = std::stod(meanTof_str);
                meanTofErr = std::stod(meanTofErr_str);
                stdTof = std::stod(stdTof_str);
                stdTofErr = std::stod(stdTofErr_str);

                // conversion from ns to ps
                stdTof *= 1000.;
                stdTofErr *= 1000.;

                R_Tof_vec.push_back(stdTof);
                R_Tof_err_vec.push_back(stdTofErr);

                meanTof_vec.push_back(meanTof);
                meanTof_err_vec.push_back(meanTofErr);

                cout << "TOF resolution (ps): " << stdTof << " +/- " << stdTofErr << endl;
                cout << "cov_mu_sigma_ToF: " << cov_mu_sigma_ToF << endl;

                WriteTofTable(LatexFile_Tof, energy, meanTof, meanTofErr,
                    stdTof, stdTofErr);

                fitResult->Write(Form("fitResultTOF"), TObject::kOverwrite);
            }

            hist->Write("", TObject::kOverwrite);
        }
        delete hist;
    }
}

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
) {
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

    FitHistograms(file, energy, R_He_vec, R_He_err_vec, meanEloss_vec, meanEloss_err_vec, R_Tof_vec, R_Tof_err_vec, meanTof_vec, meanTof_err_vec, LatexFile_Eloss, LatexFile_Tof);

    outputFile->Close();
    file->Close();
    delete outputFile;
    delete file;
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
        if ((TString(hist->GetName()) == "Eloss") && energy == 220) peakPosition = 5.22;
        int binFit = (TString(hist->GetName()) == "Tof") ? 5 : 7;
        int binPeak = hist->FindBin(peakPosition);
        int binLow = std::max(1, binPeak - binFit);
        int binHigh = std::min(hist->GetNbinsX(), binPeak + binFit);
        TF1 *gaus = new TF1("gaus", "gaus", hist->GetBinLowEdge(binLow), hist->GetBinLowEdge(binHigh + 1));
        gaus->SetParameter(1, peakPosition);
        gaus->SetParLimits(1, peakPosition - 0.5, peakPosition + 0.5);
        fitResult = hist->Fit("gaus", "RS");
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


void FitHistogramWithUniformFluctuation(TH1D* hOriginal, int nTrials = 100) {
    TRandom3 rand(0); // Random number generator
    std::vector<double> meanValues;

    for (int i = 0; i < nTrials; ++i) {
        // Clone and reset the histogram
        TH1* hToy = (TH1D*)hOriginal->Clone(Form("hToy_%d", i));
        hToy->Reset();

        // Fill each bin with uniformly fluctuated bin content
        for (int bin = 1; bin <= hOriginal->GetNbinsX(); ++bin) {
            double n = hOriginal->GetBinContent(bin);
            double delta = rand.Uniform(-std::sqrt(n), std::sqrt(n));
            double newContent = n + delta;
            if (newContent < 0) newContent = 0; // Clamp to avoid negative content
            hToy->SetBinContent(bin, newContent);
        }

        // Fit the toy histogram
        TF1* fitFunc = new TF1("fitFunc", "gaus", hToy->GetXaxis()->GetXmin(), hToy->GetXaxis()->GetXmax());
        if (hToy->Fit(fitFunc, "Q0") == 0) { // Fit successful
            double mean = fitFunc->GetParameter(1);
            meanValues.push_back(mean);
        }

        delete hToy;
        delete fitFunc;
    }

    int nSuccess = meanValues.size();
    if (nSuccess < 2) {
        std::cout << "Not enough successful fits to compute statistics." << std::endl;
        return;
    }

    // Compute mean of the means
    double sum = std::accumulate(meanValues.begin(), meanValues.end(), 0.0);
    double meanOfMeans = sum / nSuccess;

    // Compute the unbiased error on the mean
    double variance = 0.0;
    for (double m_i : meanValues) {
        variance += (m_i - meanOfMeans) * (m_i - meanOfMeans);
    }
    double stddev = std::sqrt(variance / (nSuccess * (nSuccess - 1)));

    // Print results
    std::cout << "Mean of fit means     : " << meanOfMeans << std::endl;
    std::cout << "Uncertainty on the mean: " << stddev << std::endl;
}


void WriteElossTable(std::ofstream& outFile, int energy, 
                     double meanElossHe, double stdElossHe, double R_He,
                     double meanElossErrHe, double stdElossErrHe, double R_He_err)
{
    if (!outFile.is_open()) {
        std::cerr << "Error: Output file stream is not open." << std::endl;
        return;
    }

    // Write LaTeX table row with errors in \pm notation
    outFile << energy << " & "
            << meanElossHe << " $\\pm $ " << meanElossErrHe << " & "
            << stdElossHe << " $\\pm $ " << stdElossErrHe << " & "
            << R_He << " $\\pm $ " << R_He_err << " \\\\ \\midrule\n";
}

void WriteTofTable(std::ofstream& outFile, int energy, 
    double meanTof, double meanTofErr,
    double stdTof, double stdTofErr)
{
if (!outFile.is_open()) {
std::cerr << "Error: Output file stream is not open." << std::endl;
return;
}

// Write LaTeX table row with errors in \pm notation
outFile << energy << " & "
<< meanTof << " $\\pm $ " << meanTofErr << " & "
<< stdTof << " $\\pm $ " << stdTofErr << " \\\\ \\midrule\n";
}

void OpenElossFile(std::ofstream& outFile)
{
    // Write LaTeX table header
    outFile << "\\begin{table}[h]\n"
    << "    \\centering\n"
    << "    \\begin{tabular}{cccc}\n"
    << "        \\toprule\n"
    << "        Beam energy & $\\mu(\\Delta E)$ & $\\sigma(\\Delta E)$ & $R_{\\Delta E}$\\\\ \n"
    << "        {[}MeV/u{]} & {[}MeV{]} & {[}MeV{]} & {[}\\%{]}\\\\ \\midrule\n";
}

void OpenTofFile(std::ofstream& outFile)
{
    // Write LaTeX table header
    outFile << "\\begin{table}[h]\n"
    << "    \\centering\n"
    << "    \\begin{tabular}{ccc}\n"
    << "        \\toprule\n"
    << "        Beam energy & $\\mu(TOF)$ & $R_{TOF}$ $(\\sigma)$\\\\ \n"
    << "        {[}MeV/u{]} & {[}ns{]} & {[}ps{]}\\\\ \\midrule\n";
}

void CloseElossFile(std::ofstream& outFile)
{
    //if (!outFile.is_open()) return;
    // Write LaTeX table footer
    outFile << "        \\bottomrule\n"
    << "    \\end{tabular}\n"
    << "    \\caption{Calibrated $\\Delta E$ resolution results at each beam energy value.}\n"
    << "    \\label{Tab:Eloss_resolutions}\n"
    << "\\end{table}\n";
    outFile.close();
    std::cout << "LaTeX table written successfully." << std::endl;
}

void CloseTofFile(std::ofstream& outFile)
{
    //if (!outFile.is_open()) return;
    // Write LaTeX table footer
    outFile << "        \\bottomrule\n"
    << "    \\end{tabular}\n"
    << "    \\caption{Calibrated TOF resolution results at each beam energy value.}\n"
    << "    \\label{Tab:TOF_resolutions}\n"
    << "\\end{table}\n";
    outFile.close();
    std::cout << "LaTeX table written successfully." << std::endl;
}

void PlotResolutionGraphs(
    const std::vector<double>& R_He_vec,
    const std::vector<double>& R_He_err_vec,
    const std::vector<double>& meanEloss_vec,
    const std::vector<double>& meanEloss_err_vec,
    const std::vector<double>& R_Tof_vec,
    const std::vector<double>& R_Tof_err_vec,
    const std::vector<double>& meanTof_vec,
    const std::vector<double>& meanTof_err_vec
) {
    gStyle->SetTitleSize(0.07, "T");

    // ----------- Canvas 1: R_He vs meanEloss -----------
    auto cEloss = new TCanvas("cEloss", "Helium Resolution vs Energy Loss", 800, 600);
    cEloss->SetGrid();
    cEloss->SetMargin(0.15, 0.1, 0.15, 0.1);

    auto gMeanEloss = new TGraphErrors(meanEloss_vec.size(),
                                        meanEloss_vec.data(),
                                        R_He_vec.data(),
                                        meanEloss_err_vec.data(),
                                        R_He_err_vec.data());

    gMeanEloss->SetTitle("Calibrated energy loss resolutions;#mu(#Delta E) [MeV];Resolution [%]");
    gMeanEloss->SetMarkerStyle(20);
    gMeanEloss->SetMarkerSize(1.5);
    gMeanEloss->SetLineWidth(2);
    gMeanEloss->GetXaxis()->SetTitleSize(0.05);
    gMeanEloss->GetYaxis()->SetTitleSize(0.05);
    gMeanEloss->Draw("AP");

    cEloss->SaveAs("Plots/Eloss_Resolution.png");

    // ----------- Canvas 2: R_Tof vs meanTof -----------
    auto cTof = new TCanvas("cTof", "TOF Resolution vs TOF", 800, 600);
    cTof->SetGrid();
    cTof->SetMargin(0.15, 0.1, 0.15, 0.1);

    auto gRTof = new TGraphErrors(meanTof_vec.size(),
                                   meanTof_vec.data(),
                                   R_Tof_vec.data(),
                                   meanTof_err_vec.data(),
                                   R_Tof_err_vec.data());

    gRTof->SetTitle("Calibrated TOF resolutions;#mu(TOF) [ns];Resolution [ps]");
    gRTof->SetMarkerStyle(20);
    gRTof->SetMarkerSize(1.5);
    gRTof->SetLineWidth(2);
    gRTof->GetXaxis()->SetTitleSize(0.05);
    gRTof->GetYaxis()->SetTitleSize(0.05);
    gRTof->Draw("AP");

    cTof->SaveAs("Plots/TOF_Resolution.png");
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