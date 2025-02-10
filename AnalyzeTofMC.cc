// Macro that fits the ToF (Time of Flight) distributions from the AnalyzeTWFragMC.cc output files (for MC runs).
// The fits are performed with a gaussian function. The fit limits depend on the bar and on the beam energy.
// The histogram peaks are automatically found with TSPectrum; the fits are then performed within a certain bin-range centered around
// the peak, depending on the specific bar and beam energy. The fit results are stored in files name like e.g.
// AnaFOOT_TW_Decoded_HIT2022_140_Fit.root (created if they don't already exist, or overwritten if they do),
// inside of the TofFit directory. To be run with root -l -b -q 'AnalyzeTofMC.cc()'.

#include "AnalyzeTofMC.h"

void AnalyzeTofMC()
{
    std::vector<std::pair<std::string, int>> filesAndEnergies = {
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_100.root", 100},
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_140.root", 140},
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_200.root", 200},
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_220.root", 220}};

    std::map<TString, std::map<int, double>> TofMeansMC;
    std::map<TString, std::map<int, double>> TofErrorsMC;

    for (const auto &[fileName, energy] : filesAndEnergies)
    {
        ProcessFile(fileName, energy, TofMeansMC, TofErrorsMC);
    }
}

// Function to fit histograms and store results
void FitHistogramsInDirectory(
    TDirectory *dir,
    int energy,
    std::map<TString, std::map<int, double>> &TofMeansMC,
    std::map<TString, std::map<int, double>> &TofErrorsMC,
    TFile *outputFile)
{
    if (!dir)
        return;

    TIter next(dir->GetListOfKeys());
    TKey *key;

    while ((key = (TKey *)next()))
    {
        TObject *obj = key->ReadObj();
        TH1D *hist = dynamic_cast<TH1D *>(obj);
        if (!hist)
        {
            delete obj;
            continue;
        }

        TString histName = hist->GetName();
        if (!histName.BeginsWith("ToF_Layer"))
        {
            delete hist;
            continue;
        }

        // Extract layer and bar identifiers
        TString layerBarCombination = histName(4, histName.Length() - 4); // e.g. LayerY_bar9
        cout << endl;
        cout << "beam energy: " << energy << " MeV/u " << layerBarCombination << endl;
        TFitResultPtr fitResult = FitPeaksWithTSpectrum(hist, energy, layerBarCombination);
        if (fitResult.Get() != nullptr)
        {
            TofMeansMC[layerBarCombination][energy] = fitResult->Parameter(1);
            TofErrorsMC[layerBarCombination][energy] = fitResult->Error(1);
            fitResult->Write(Form("fitResultTof_%s", layerBarCombination.Data()), TObject::kOverwrite);
            hist->Write("", TObject::kOverwrite);
        }
        delete hist;
    }
}

void ProcessFile(
    const std::string &fileName,
    int energy,
    std::map<TString, std::map<int, double>> &TofMeansMC,
    std::map<TString, std::map<int, double>> &TofErrorsMC)
{
    TFile *file = TFile::Open(fileName.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        return;
    }

    std::cout << "Processing file: " << fileName << std::endl;

    // Create an output ROOT file to store fitted histograms
    TString outputFileName = TString(fileName).ReplaceAll(".root", "_Fit.root");
    // Update the output file if it already exists
    TFile *outputFile = TFile::Open(outputFileName, "UPDATE");

    TDirectory *fitDir = outputFile->GetDirectory("TofFit");
    if (!fitDir)
    {
        fitDir = outputFile->mkdir("TofFit"); // Create only if it doesn't exist
    }
    fitDir->cd(); // Move into the directory before writing

    const std::vector<std::string> directories = {"TofLayerX", "TofLayerY"};
    for (const auto &dirName : directories)
    {
        TDirectory *dir = dynamic_cast<TDirectory *>(file->Get(dirName.c_str()));
        if (!dir)
        {
            std::cerr << "Error: Directory " << dirName << " not found in file " << fileName << std::endl;
            continue;
        }
        FitHistogramsInDirectory(dir, energy, TofMeansMC, TofErrorsMC, outputFile);
    }

    outputFile->Close();
    file->Close();
    delete outputFile;
    delete file;
}

TFitResultPtr FitPeaksWithTSpectrum(TH1D *hist, int energy, const TString &layerBarCombination)
{
    int layerStart = layerBarCombination.Index("Layer") + 5; // "Layer" is 5 characters long
    int barStart = layerBarCombination.Index("_bar") + 4;    // "_bar" is 4 characters long

    // Extract the layer (e.g., 'X' or 'Y')
    TString layer = layerBarCombination(layerStart, 1); // Extract one character at layerStart

    // Extract the bar number as a string and convert to integer
    TString barString = layerBarCombination(barStart, layerBarCombination.Length() - barStart);
    int barNumber = barString.Atoi();

    // outputFile->cd();
    hist->Draw();
    // Fit histogram
    TSpectrum spectrum(1);
    int nPeaks = spectrum.Search(hist, 1, "", 0.0015);
    cout << "# of peaks: " << nPeaks << endl;
    TFitResultPtr fitResult;
    double peakPosition;
    int binFit;

    if (nPeaks == 0)
    {
        std::cerr << "No peaks found beyond threshold " << 1. << std::endl;
        delete hist;
        return nullptr;
    }
    else
    {
        peakPosition = spectrum.GetPositionX()[0];
        binFit = 3;
        if (energy == 100)
        {
            if (layer == "Y")
            {
                if (barNumber == 0)
                {
                    binFit = 5;
                }
                else if (barNumber == 1 || barNumber == 4) {
                    binFit = 4;
                }
                else if (barNumber == 5 || barNumber == 14) {
                    binFit = 4;
                }
                else if (barNumber == 3 || barNumber == 16 || barNumber == 17 || barNumber == 18)
                {
                    binFit = 5;
                }
                else if (barNumber == 19)
                {
                    binFit = 7;
                }
            }
            else if (layer == "X")
            {
                if (barNumber == 1) {
                    binFit = 6;
                }
                else if (barNumber == 0) {
                    binFit = 5;
                }
                else if (barNumber ==2 || barNumber == 5) {
                    binFit = 4;
                }
                else if (barNumber == 8 || barNumber == 10) {
                    binFit = 2;
                }
                else if (barNumber == 14 || barNumber == 16) {
                    binFit = 4;
                }
                else if (barNumber == 17)
                {
                    peakPosition = 11.7;
                    binFit = 7;
                }
                else if (barNumber == 18) {
                    binFit = 8;
                }
                else if (barNumber == 19) {
                    binFit = 9;
                }
            }
        }
        else if (energy == 140)
        {
            if (layer == "Y")
            {
                if (barNumber == 0)
                {
                    peakPosition = 10.3;
                    binFit = 15;
                }
                else if (barNumber == 2) {
                    binFit = 4;
                }
                else if (barNumber == 3) {
                    binFit = 7;
                }
                else if (barNumber == 4)
                {
                    binFit = 6;
                }
                else if (barNumber == 5) {
                    binFit = 7;
                }
                else if (barNumber == 13 || barNumber == 14 || barNumber == 15) {
                    binFit = 6;
                }
                else if (barNumber == 16)
                {
                    peakPosition = 10.;
                    binFit = 10;
                }
                else if (barNumber == 17)
                {
                    peakPosition = 10.1;
                    binFit = 8;
                }
                else if (barNumber == 18)
                {
                    peakPosition = 10.2;
                    binFit = 10;
                }
                else if (barNumber == 19)
                {
                    peakPosition = 10.1;
                    binFit = 8;
                }
            }
            else if (layer == "X")
            {
                if (barNumber < 3)
                {
                    binFit = 10;
                    if (barNumber == 1)
                    {
                        binFit = 12;
                    }
                    else if (barNumber == 0 || barNumber == 2) {
                        peakPosition = 10.2;
                    }

                }
                else if (barNumber == 3)
                {
                    binFit = 6;
                }
                else if (barNumber == 4) {
                    binFit = 7;
                }
                else if (barNumber == 7) {
                    binFit = 2;
                }
                else if (barNumber == 13 || barNumber == 14 || barNumber == 15)
                {
                    binFit = 5;
                }
                else if (barNumber == 16)
                {
                    peakPosition = 10.1;
                    binFit = 10;
                }
                else if (barNumber == 17) {
                    binFit = 6;
                }
                else if (barNumber == 18) {
                    peakPosition = 10.2;
                    binFit = 10;
                } 
                else if (barNumber == 19)
                {
                    peakPosition = 10.3;
                    binFit = 11;
                }
            }
        }
        else if (energy == 200)
        {
            if (layer == "Y")
            {
                if (barNumber == 0)
                {
                    peakPosition = 8.8;
                    binFit = 12;
                }
                else if (barNumber == 1) {
                    binFit = 5;
                }
                else if (barNumber == 3) {
                    binFit = 7;
                }
                else if (barNumber == 2)
                {
                    peakPosition = 8.6;
                    binFit = 10;
                }
                else if (barNumber == 4) {
                    binFit = 6;
                }
                else if (barNumber == 7) {
                    binFit = 3;
                }
                else if (barNumber == 12) {
                    binFit = 4;
                }
                else if (barNumber == 13) {
                    binFit = 6;
                }
                else if (barNumber == 14)
                {
                    peakPosition = 8.6;
                    binFit = 9;
                }
                else if (barNumber == 15)
                {
                    binFit = 6;
                }
                else if (barNumber == 16)
                {
                    peakPosition = 8.7;
                    binFit = 10;
                }
                else if (barNumber == 17) {
                    peakPosition = 8.9;
                    binFit = 12;
                }
                else if (barNumber == 18)
                {
                    peakPosition = 9.0;
                    binFit = 13;
                }
                else if (barNumber == 19)
                {
                    peakPosition = 8.8;
                    binFit = 12;
                }
            }
            else if (layer == "X")
            {
                if (barNumber == 0)
                {
                    peakPosition = 8.9;
                    binFit = 15;
                }
                else if (barNumber == 1 || barNumber == 2)
                {
                    peakPosition = 8.8;
                    binFit = (barNumber == 1) ? 10 : 9;
                }
                else if (barNumber == 3 || barNumber == 4 || barNumber == 5) {
                    binFit = 6;
                }
                else if (barNumber == 6 || barNumber == 7 || barNumber == 8) {
                    binFit = 3;
                }
                else if (barNumber == 13) {
                    binFit = 7;
                }
                else if (barNumber == 14) {
                    binFit = 5;
                }
                else if (barNumber == 15)
                {
                    binFit = 8;
                }
                else if (barNumber == 16)
                {
                    peakPosition = 8.7;
                    binFit = 10;
                }
                else if (barNumber == 17)
                {
                    peakPosition = 8.9;
                    binFit = 13;
                }
                else if (barNumber == 18)
                {
                    peakPosition = 8.9;
                    binFit = 10;
                }
                else if (barNumber == 19)
                {
                    peakPosition = 9.2;
                    binFit = 15;
                }
            }
        }
        else if (energy == 220) {
            binFit = 4;
            if (layer == "Y") {
                if (barNumber == 0) {
                    peakPosition = 8.6;
                    binFit = 10;
                }
                else if (barNumber == 1) {
                    peakPosition = 8.5;
                    binFit = 10;
                }
                else if (barNumber == 2) {
                    peakPosition = 8.3;
                    binFit = 8;
                }
                else if (barNumber == 3) {
                    binFit = 6;
                }
                else if (barNumber == 4 || barNumber == 14) {
                    binFit = 5;
                }
                else if (barNumber == 8 || barNumber == 10) {
                    binFit = 3;
                }
                else if (barNumber == 13) {
                    binFit = 6;
                }
                else if (barNumber == 15 || barNumber == 16) {
                    binFit = 9;
                }
                else if (barNumber == 17) {
                    peakPosition = 8.4;
                    binFit = 9;
                }
                else if (barNumber == 18) {
                    binFit = 11;
                }
                else if (barNumber == 19) {
                    binFit = 9;
                }
            }
            else if (layer == "X") {
                if (barNumber == 0) {
                    peakPosition = 8.5;
                    binFit = 12;
                }
                else if (barNumber == 1) {
                    peakPosition = 8.45;
                    binFit = 10;
                }
                else if (barNumber == 2) {
                    peakPosition = 8.5;
                    binFit = 10;
                }
                else if (barNumber == 3 || barNumber == 4) {
                    binFit = 5;
                }
                else if (barNumber == 8 || barNumber == 10) {
                    binFit = 3;
                }
                else if (barNumber == 13) {
                    binFit = 6;
                }
                else if (barNumber == 14) {
                    binFit = 5;
                }
                else if (barNumber == 15) {
                    peakPosition = 8.35;
                    binFit = 10;
                }
                else if (barNumber == 16) {
                    binFit = 10;
                }
                else if (barNumber == 17) {
                    binFit = 7;
                }
                else if (barNumber == 18) {
                    binFit = 10;
                }
                else if (barNumber == 19) {
                    peakPosition = 8.7;
                    binFit = 12;
                }
            }
        }

    }
    int binPeak = hist->FindBin(peakPosition);
    int binLow = std::max(1, binPeak - binFit);
    int binHigh = std::min(hist->GetNbinsX(), binPeak + binFit);
    TF1 *gaus = new TF1("gaus", "gaus", hist->GetBinLowEdge(binLow), hist->GetBinLowEdge(binHigh + 1));
    gaus->SetParameter(1, peakPosition);
    gaus->SetParLimits(1, peakPosition - 0.1, peakPosition + 0.1);

    fitResult = hist->Fit("gaus", "RS");
    if (fitResult->Status() != 0)
    {
        cout << "fit failed" << endl;
        delete hist;
        return nullptr;
    }
    else
    {
        hist->GetListOfFunctions()->Add(gaus);
        gaus->Draw("same");
    }

    return fitResult;
}
