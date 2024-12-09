// Macro that fits the charge distribution from the AnalyzeTWChargeTime.cc merged output files.
// Only the histograms whose entries are greater than a fraction of the sum of the merged files
// total entries are fitted. The fit results are stored in files name like e.g.
// TW/AnaFOOT_TW_Decoded_HIT2022_140MeV_Fit.root. To be run with root -l -b -q 'AnalyzeTWFit.cc()'.
// The mean charge values given by the fits are then plotted vs the beam energies, one plot for every
// layer-bar combination.

#include "AnalyzeTWFit.h"

double DENSITY = 1.023;  //EJ200 density g/cm^3
double THICKNESS = 0.3;  //slab scintillator thickness cm

// Main analysis function
void AnalyzeTWFit() {
    std::vector<std::pair<std::string, double>> filesAndEnergies = {
        {"TW/AnaFOOT_TW_Decoded_HIT2022_100MeV.root", 100},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_140MeV.root", 140},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_180MeV.root", 180},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_200MeV.root", 200},
        {"TW/AnaFOOT_TW_Decoded_HIT2022_220MeV.root", 220}
    };

    std::map<TString, std::map<int, double>> fitMeans;
    std::map<TString, std::map<int, double>> fitErrors;
    std::map<int, double> stopping_power = {
    {100, 7.245},
    {140, 5.674},
    {180, 4.777},
    {200, 4.458},
    {220, 4.196}
    };  // stopping power for protons, MeV*cm^2/g, to be multiplied by 4, density and thickness to obtain E_loss

    for (const auto& [fileName, energy] : filesAndEnergies) {
        ProcessFile(fileName, energy, fitMeans, fitErrors);
    }

    for (const auto& layerBarEntry : fitMeans) {
        CreateAndSaveGraph(layerBarEntry.first, layerBarEntry.second, fitErrors, stopping_power);
    }
}

// Function to sum all "nentries" values in a file
int SumNentries(TFile* file) {
    int totalNentries = 0;
    TIter next(file->GetListOfKeys());
    TKey* key;
    std::regex nentriesRegex("nentries_run_\\d+");

    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        TNamed* nentriesObj = dynamic_cast<TNamed*>(obj);
        if (nentriesObj && std::regex_match(nentriesObj->GetName(), nentriesRegex)) {
            totalNentries += std::stoi(nentriesObj->GetTitle());
        }
        delete obj;
    }
    return totalNentries;
}

// Function to fit histograms in a directory and store results
void FitHistogramsInDirectory(
    TDirectory* dir, 
    double threshold, 
    double energy, 
    std::map<TString, std::map<int, double>>& fitMeans, 
    std::map<TString, std::map<int, double>>& fitErrors, 
    TFile* outputFile
) {
    if (!dir) return;

    TIter next(dir->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        TH1* hist = dynamic_cast<TH1*>(obj);
        if (!hist) {
            delete obj;
            continue;
        }

        TString histName = hist->GetName();
        if (!histName.BeginsWith("Charge_Layer")) {
            delete hist;
            continue;
        }

        // Extract layer and bar identifiers
        TString layerBarCombination = histName(7, histName.Length() - 7);  // e.g. LayerY_bar9

        // Fit histogram if above the threshold
        if (hist->GetEntries() > threshold) {
            TF1* fitFunc = new TF1("fitFunc", "gaus");
            TFitResultPtr fitResult = hist->Fit(fitFunc, "QS", " ", 2., 13.);

            double meanCharge = fitFunc->GetParameter(1);
            double stdCharge = fitFunc->GetParameter(2);
            double meanChargeErr = fitFunc->GetParError(1);

            if (!meanCharge || meanCharge < 0 || stdCharge / meanCharge > 0.2 || meanCharge / meanChargeErr - 1 < 0.5) {
                std::cout << layerBarCombination << " Energy: " << energy << " MeV: Skipping fit (negative charge or std/charge > 0.2 or charge error too big)" << std::endl;
                continue;
            }
            fitMeans[layerBarCombination][energy] = meanCharge;
            fitErrors[layerBarCombination][energy] = meanChargeErr;

            // Save the histogram with fit to the output file
            outputFile->cd();
            hist->Write();
            fitResult->Write(Form("FitResult_%s", layerBarCombination.Data()));  // Save the TFitResult
            delete fitFunc;
        }

        delete hist;
    }
}

// Function to process a single file and extract fit data
void ProcessFile(
    const std::string& fileName, 
    double energy, 
    std::map<TString, std::map<int, double>>& fitMeans, 
    std::map<TString, std::map<int, double>>& fitErrors
) {
    TFile* file = TFile::Open(fileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        return;
    }

    int totalNentries = SumNentries(file);
    double threshold = 0.009 * totalNentries;

    std::cout << "Processing file: " << fileName << ", Total nentries: " << totalNentries << ", Threshold: " << threshold << std::endl;

    // Create an output ROOT file to store fitted histograms
    TString outputFileName = TString(fileName).ReplaceAll(".root", "_Fit.root");
    TFile* outputFile = TFile::Open(outputFileName, "RECREATE");

    const std::vector<std::string> directories = {"ChargeTimeLayerX", "ChargeTimeLayerY"};
    for (const auto& dirName : directories) {
        TDirectory* dir = dynamic_cast<TDirectory*>(file->Get(dirName.c_str()));
        if (!dir) {
            std::cerr << "Error: Directory " << dirName << " not found in file " << fileName << std::endl;
            continue;
        }
        FitHistogramsInDirectory(dir, threshold, energy, fitMeans, fitErrors, outputFile);
    }

    outputFile->Close();
    file->Close();
    delete outputFile;
    delete file;
}

// Function to create and save a graph for a given layer-bar combination
void CreateAndSaveGraph(
    const TString& layerBarCombination, 
    const std::map<int, double>& energiesAndFits, 
    const std::map<TString, std::map<int, double>>& fitErrors,
    std::map<int, double>& stopping_power
) {
    TGraphErrors* graph = new TGraphErrors();
    graph->SetMinimum(0.0);
    int pointIndex = 0;

    for (const auto& [energy, fitMean] : energiesAndFits) {
        double fitError = fitErrors.at(layerBarCombination).at(energy);
        double e_loss = stopping_power[energy] * 4 * DENSITY * THICKNESS;
        cout << "Beam Energy: " << energy << " MeV " << "E_loss: " << e_loss << " MeV" << endl;
        graph->SetPoint(pointIndex, e_loss, fitMean);
        graph->SetPointError(pointIndex, 0, fitError);
        ++pointIndex;
    }

    double x_max = graph->GetXaxis()->GetXmax();
    graph->GetXaxis()->SetLimits(0., x_max);

    // Perform linear fit
    TF1 *f1 = new TF1("f1", "pol1", 0., x_max);
    TFitResultPtr fitresult_linear = graph->Fit(f1, "S");

    if (!fitresult_linear.Get() || fitresult_linear->Status() != 0) {
        std::cerr << "Linear fit failed for " << layerBarCombination.Data() << std::endl;
    } else {
        auto [intercept, sigma_intercept] = RoundMeasurement(fitresult_linear->Parameter(0), fitresult_linear->ParError(0));
        auto [slope, sigma_slope] = RoundMeasurement(fitresult_linear->Parameter(1), fitresult_linear->ParError(1));
        double chi2 = fitresult_linear->Chi2();
        int ndf = fitresult_linear->Ndf();

        TPaveText *fitInfo = new TPaveText(0.3, 0.7, 0.45, 0.85, "NDC");
        fitInfo->SetFillColor(0);
        fitInfo->AddText(Form("Intercept [a.u.] = %f#pm %f", intercept, sigma_intercept));
        fitInfo->AddText(Form("Slope [a.u. / MeV] = %f#pm %f", slope, sigma_slope));
        fitInfo->AddText(Form("#chi^{2} / ndf = %.2f / %d", chi2, ndf));
        fitInfo->SetTextSize(0.03);

        f1->SetParameter(0, intercept);
        f1->SetParameter(1, slope);

        TCanvas* c = new TCanvas("c", "Fit Results", 800, 600);
        graph->SetTitle(Form("Merged Fit Mean Charge vs Energy Loss for %s", layerBarCombination.Data()));
        graph->GetXaxis()->SetTitle("Energy Loss [MeV]");
        graph->GetYaxis()->SetTitle("Mean Charge [a.u.]");
        c->SetLeftMargin(0.15);
        graph->SetMarkerStyle(24);
        graph->Draw("AP");

        f1->Draw("same");
        c->cd();
        fitInfo->Draw("same");
        c->Modified();
        c->Update();

        c->SaveAs(Form("Plots/Merged_Fit_%s.png", layerBarCombination.Data()));
        delete c;
        delete graph;
        delete fitInfo;
    }
}

pair<double, double> RoundMeasurement(double value, double uncertainty) {
    //Get the order of magnitude of the uncertainty
    int significantFigures = (int)std::ceil(-std::log10(uncertainty)) + 1;

    //Calculate the rounding factor based on significant figures
    double roundingFactor = std::pow(10, significantFigures);

    //Round the uncertainty and value accordingly
    double roundedUncertainty = std::round(uncertainty * roundingFactor) / roundingFactor;
    double roundedValue = std::round(value * roundingFactor) / roundingFactor;

    return {roundedValue, roundedUncertainty};
}
