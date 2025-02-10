// Macro that fits the ToF (Time of Flight) distributions from the AnalyzeTWFragm.cc merged output files (for fragmentation runs).
// The fits are performed with a gaussian function. The fit limits depend on the bar and on the beam energy.
// The histogram peaks are automatically found with TSPectrum; the fits are then performed within a certain bin-range centered around
// the peak, depending on the specific bar and beam energy. The fit results are stored in files name like e.g.
// AnaFOOT_TW_Decoded_HIT2022_140MeV_Fit.root (created if they don't already exist, or overwritten if they do),
// inside of the TofFit directory. To be run with root -l -b -q 'AnalyzeTofFragm.cc()'.

#include "AnalyzeTofFragm.h"

void AnalyzeTofFragm()
{
    std::vector<std::pair<std::string, int>> filesAndEnergies = {
        {"TW/cuts/AnaFOOT_TW_Decoded_HIT2022_fragm_100MeV.root", 100},
        {"TW/cuts/AnaFOOT_TW_Decoded_HIT2022_fragm_140MeV.root", 140},
        {"TW/cuts/AnaFOOT_TW_Decoded_HIT2022_fragm_200MeV.root", 200},
        {"TW/cuts/AnaFOOT_TW_Decoded_HIT2022_fragm_220MeV.root", 220}};

    std::map<TString, std::map<int, double>> TofMeans;
    std::map<TString, std::map<int, double>> TofErrors;

    for (const auto &[fileName, energy] : filesAndEnergies)
    {
        ProcessFile(fileName, energy, TofMeans, TofErrors);
    }
}

// Function to fit histograms in a directory and store results
void FitHistogramsInDirectory(
    TDirectory *dir,
    int energy,
    std::map<TString, std::map<int, double>> &TofMeans,
    std::map<TString, std::map<int, double>> &TofErrors,
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
        if (fitResult.Get() != nullptr) {
            TofMeans[layerBarCombination][energy] = fitResult->Parameter(1);
            TofErrors[layerBarCombination][energy] = fitResult->Error(1);
            fitResult->Write(Form("fitResultTof_%s", layerBarCombination.Data()), TObject::kOverwrite);
            hist->Write("", TObject::kOverwrite);

        }
        delete hist;
    }
}

void ProcessFile(
    const std::string &fileName,
    int energy,
    std::map<TString, std::map<int, double>> &TofMeans,
    std::map<TString, std::map<int, double>> &TofErrors)
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

    TDirectory* fitDir = outputFile->GetDirectory("TofFit");
    if (!fitDir) {
        fitDir = outputFile->mkdir("TofFit");  // Create only if it doesn't exist
    }
    fitDir->cd();  // Move into the directory before writing

    const std::vector<std::string> directories = {"ToFLayerX", "ToFLayerY"};
    for (const auto &dirName : directories)
    {
        TDirectory *dir = dynamic_cast<TDirectory *>(file->Get(dirName.c_str()));
        if (!dir)
        {
            std::cerr << "Error: Directory " << dirName << " not found in file " << fileName << std::endl;
            continue;
        }
        FitHistogramsInDirectory(dir, energy, TofMeans, TofErrors, outputFile);
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

    //outputFile->cd();
    hist->Draw();
    // Fit histogram
    TSpectrum spectrum(1);
    int nPeaks = spectrum.Search(hist, 1, "", 0.0015);
    cout << "# of peaks: " << nPeaks << endl;
    TFitResultPtr fitResult;

    if (nPeaks == 0)
    {
        std::cerr << "No peaks found beyond threshold " << 1. << std::endl;
        delete hist;
        return nullptr;
    }
    else
    {
        double peakPosition = spectrum.GetPositionX()[0];
        int binFit = 4;
        if (energy == 100)
        {
            if (layer == "X")
            {
                if (barNumber == 0)
                {
                    binFit = 5;
                }
                else if (barNumber == 17)
                {
                    peakPosition = 11.48;
                    binFit = 10;
                }
                else if (barNumber == 18)
                {
                    peakPosition = 12.34;
                    binFit = 10;
                }
                else if (barNumber == 19)
                {
                    peakPosition = 11.9;
                    binFit = 10;
                }
            }
            else if (layer == "Y")
            {
                if (barNumber == 0)
                {
                    peakPosition = 11.7;
                    binFit = 8;
                }
                else if (barNumber == 1)
                {
                    binFit = 8;
                }
                else if (barNumber == 13) {
                    binFit = 3;
                }
                else if (barNumber == 19)
                {
                    peakPosition = 12.5;
                    binFit = 20;
                }
            }
        }
        else if (energy == 140)
        {
            if (layer == "X")
            {
                if (barNumber == 0)
                {
                    peakPosition = 10.1;
                    binFit = 10;
                }
                else if (barNumber == 1 || barNumber == 2)
                {
                    binFit = 10;
                }
                else if (barNumber == 3)
                {
                    binFit = 13;
                }
                else if (barNumber == 4)
                {
                    binFit = 6;
                }
                else if (barNumber == 14)
                {
                    binFit = 5;
                }
                else if (barNumber == 15)
                {
                    binFit = 10;
                }
                else if (barNumber == 17)
                {
                    peakPosition = 9.6;
                    binFit = 11;
                }
                else if (barNumber == 18)
                {
                    peakPosition = 10.5;
                    binFit = 16;
                }
                else if (barNumber == 19)
                {
                    peakPosition = 10.3;
                    binFit = 14;
                }
            }
            else if (layer == "Y")
            {
                if (barNumber == 0)
                {
                    peakPosition = 9.8;
                    binFit = 15;
                }
                else if (barNumber == 1)
                {
                    binFit = 6;
                }
                else if (barNumber == 3)
                {
                    binFit = 7;
                }
                else if (barNumber == 4)
                {
                    binFit = 6;
                }
                else if (barNumber == 13 || barNumber == 14)
                {
                    binFit = 6;
                }
                else if (barNumber == 15)
                {
                    binFit = 7;
                }
                else if (barNumber == 16 || barNumber == 17 || barNumber == 19)
                {
                    binFit = 10;
                    if (barNumber == 17) {
                        binFit = 11;
                    }
                }
                else if (barNumber == 18)
                {
                    peakPosition = 9.5;
                    binFit = 12;
                }
            }
        }

        else if (energy == 200)
        {
            binFit = 5;
            if (layer == "X")
            {
                if (barNumber < 3)
                {
                    binFit = 10;
                }
                if (barNumber == 3)
                {
                    peakPosition = 8.85;
                    binFit = 7;
                }
                else if (barNumber == 15)
                {
                    peakPosition = 8.19;
                    binFit = 7;
                }
                else if (barNumber == 16)
                {
                    binFit = 7;
                }
                else if (barNumber == 17)
                {
                    binFit = 10;
                }
                else if (barNumber == 18)
                {
                    peakPosition = 9.45;
                    binFit = 10;
                }
                else if (barNumber == 19)
                {
                    peakPosition = 9.2;
                    binFit = 10;
                }
            }
            else if (layer == "Y")
            {
                if (barNumber == 0)
                {
                    peakPosition = 8.81;
                    binFit = 10;
                }
                else if (barNumber == 1)
                {
                    peakPosition = 8.7;
                    binFit = 8;
                }
                else if (barNumber == 2)
                {
                    binFit = 7;
                }
                else if (barNumber == 3 || barNumber == 4)
                {
                    binFit = 7;
                }
                else if (barNumber == 17 || barNumber == 18)
                {
                    binFit = 10;
                }
                else if (barNumber == 19)
                {
                    // peakPosition = 9.1;
                    binFit = 10;
                }
            }
        }
        else if (energy == 220)
        {
            binFit = 5;
            if (layer == "X")
            {
                if (barNumber == 0)
                {
                    binFit = 10;
                }
                else if (barNumber == 2)
                {
                    peakPosition = 8.5;
                    binFit = 10;
                }
                else if (barNumber == 3 || barNumber == 4)
                {
                    binFit = 7;
                }
                else if (barNumber == 14)
                {
                    peakPosition = 8.4;
                    binFit = 7;
                }
                else if (barNumber == 13 || barNumber == 15)
                {
                    binFit = 7;
                }
                else if (barNumber == 16)
                {
                    peakPosition = 8.3;
                    binFit = 10;
                }
                else if (barNumber == 18)
                {
                    peakPosition = 9.2;
                    binFit = 10;
                }
                else if (barNumber == 19)
                {
                    binFit = 10;
                }
            }
            else if (layer == "Y")
            {
                if (barNumber == 0)
                {
                    binFit = 10;
                }
                else if (barNumber == 1)
                {
                    peakPosition = 8.5;
                    binFit = 6;
                }
                else if (barNumber == 2) {
                    binFit = 7;
                }
                else if (barNumber == 14)
                {
                    peakPosition = 8.3;
                    binFit = 8;
                }
                else if (barNumber == 16)
                {
                    peakPosition = 8.6;
                    binFit = 10;
                }
                else if (barNumber == 17 || barNumber == 19)
                {
                    binFit = 10;
                }
                else if (barNumber == 18)
                {
                    peakPosition = 7.7;
                    binFit = 10;
                }
            }
        }
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
            delete hist;
            return nullptr;
        }
        else
        {
            hist->GetListOfFunctions()->Add(gaus);
            gaus->Draw("same");
        }
    }
    return fitResult;
}