#include "Calo_resolutions.h"

void Calo_resolutions()
{
    std::vector<std::pair<std::string, int>> filesAndEnergies = {
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_100MeV_Fit.root", 100},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_140MeV_Fit.root", 140},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_180MeV_Fit.root", 180},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_200MeV_Fit.root", 200},
        {"Calo/AnaFOOT_Calo_Decoded_HIT2022_220MeV_Fit.root", 220}
    };

    std::vector<int> energies = {100, 140, 180, 200, 220};
    // resolution values: for each crystal there's a vector of at most 5 values
    // 1 value for each beam energy
    std::map<int, std::vector<std::pair<int, double>>> Res;
    std::map<int, std::vector<std::pair<int, double>>> Res_err;
    // crystals with at least two points in the mean charge vs beam energy plot
    std::vector<int> validCrystals = {0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 15, 16, 18, 20, 21, 23, 24, 26};

    for (const auto &[fileName, energy] : filesAndEnergies)
    {
        ProcessFile(fileName, energy, Res, Res_err);
    }

    PlotResolutions(Res, Res_err);



}

void ProcessFile(const std::string& fileName, const int energy, std::map<int, std::vector<std::pair<int, double>>>& Res, std::map<int, std::vector<std::pair<int, double>>>& Res_err)
{
    TFile* inFile = TFile::Open(fileName.c_str(), "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        return;
    }

    cout << "Processing energy: " << energy << " MeV/u" << endl;

    TIter next(inFile->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next()))
    {
        TObject* obj = key->ReadObj();
        TString objectName = key->GetName();
        if (objectName.BeginsWith("Fit_result"))
        {
            int lastUnderscore = objectName.Last('_');
            TString numberPart = objectName(lastUnderscore + 1, objectName.Length() - lastUnderscore - 1);
            int crystalId = numberPart.Atoi();

            cout << "crystal ID: " << crystalId << endl;

            TFitResult* fitResult = (TFitResult *) inFile->Get(objectName);
            double mean = fitResult->Parameter(1);
            double mean_err = fitResult->ParError(1);
            double std = fitResult->Parameter(2);
            double std_err = fitResult->ParError(2);
            double cov_mean_std = fitResult->GetCovarianceMatrix()(1, 2);

            // resolution in %
            double res = 100. * std / mean;
            double res_err = 100 * (1. / mean) * sqrt(pow(std_err, 2) + pow(res * mean_err, 2));

            auto [res_round, res_err_round] = RoundMeasurement(res, res_err);
            res = res_round;
            res_err = res_err_round;

            cout << "mean: " << mean << " +/- " << mean_err << endl;
            cout << "std: " << std << " +/- " << std_err << endl;
            cout << "covariance: " << cov_mean_std << endl;
            cout << "Resolution (%): " << res << " +/- " << res_err << endl;

            Res[crystalId].emplace_back(energy, res);
            Res_err[crystalId].emplace_back(energy, res_err);

        }
    }
}

std::pair<double, double> RoundMeasurement(double value, double uncertainty) {
    int significantFigures = (int)std::ceil(-std::log10(uncertainty));
    double roundingFactor = std::pow(10, significantFigures);
    double roundedUncertainty = std::round(uncertainty * roundingFactor) / roundingFactor;
    double roundedValue = std::round(value * roundingFactor) / roundingFactor;
    return {roundedValue, roundedUncertainty};
}

void PlotResolutions(const std::map<int, std::vector<std::pair<int, double>>>& Res,
    const std::map<int, std::vector<std::pair<int, double>>>& Res_err)
{
    for (const auto& [crystalId, resVec] : Res)
    {
        const auto& errVec = Res_err.at(crystalId);

        int nPoints = resVec.size();
        std::vector<double> energies(nPoints), resolutions(nPoints), errors(nPoints);

        for (int i = 0; i < nPoints; ++i)
        {
        energies[i] = resVec[i].first;
        resolutions[i] = resVec[i].second;
        errors[i] = errVec[i].second;
        }

        TGraphErrors* graph = new TGraphErrors(nPoints, energies.data(), resolutions.data(), nullptr, errors.data());
        graph->SetTitle(Form("Resolution vs beam energy crystal ID %d;Beam Energy [MeV/u];Resolution [%%]", crystalId));
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.5);
        gStyle->SetTitleSize(0.07, "T");
        graph->GetXaxis()->SetTitleSize(0.05);
        graph->GetYaxis()->SetTitleSize(0.05);

        //TF1* func = new TF1("res_energy", "[0]/sqrt(x)", 100, 220);
        //func->SetParameters(10. * resolutions[0]);
        //graph->Fit(func, "R");
        //func->SetLineColor(kRed);
        //func->SetLineWidth(2);


        TCanvas* c = new TCanvas(Form("c_crystal_%d", crystalId), Form("Crystal %d", crystalId), 800, 600);
        c->SetMargin(0.12, 0.12, 0.15, 0.15); // Left, Right, Bottom, Top margins
        c->SetGrid();
        graph->Draw("AP");
        //func->Draw("SAME");
        c->SaveAs(Form("Plots/CaloRes_crystal_%d.png", crystalId));
    }
}
