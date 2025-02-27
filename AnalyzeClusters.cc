
#include "AnalyzeClusters.h"

void AnalyzeClusters()
{

    std::vector<std::pair<std::string, int>> filesAndEnergies = {
        {"prova/AnaFOOT_Calo_Decoded_HIT2022_100MeV.root", 100},
        {"prova/AnaFOOT_Calo_Decoded_HIT2022_140MeV.root", 140},
        {"prova/AnaFOOT_Calo_Decoded_HIT2022_180MeV.root", 180},
        {"prova/AnaFOOT_Calo_Decoded_HIT2022_200MeV.root", 200},
        {"prova/AnaFOOT_Calo_Decoded_HIT2022_220MeV.root", 220}};

    for (const auto &[fileName, energy] : filesAndEnergies)
    {
        ProcessFile(fileName, energy);
    }
}

// Function to process a ROOT file
void ProcessFile(const std::string &fileName, int energy)
{
    // Open the input file
    TFile *f = TFile::Open(fileName.c_str(), "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return;
    }

    // Create an output file to store the results
    TString outputFileName = TString(fileName).ReplaceAll(".root", "_Fit.root");
    TFile *fout = TFile::Open(outputFileName.Data(), "RECREATE");
    if (!fout || fout->IsZombie())
    {
        std::cerr << "Error creating output file: " << outputFileName << std::endl;
        f->Close();
        return;
    }

    // Loop over crystal pairs
    for (int i = 0; i <= 8; ++i)
    {
        for (int j = i + 1; j <= 8; ++j)
        {
            std::string graphName = Form("scatter_Crystal%d_vs_Crystal%d", i, j);

            // Retrieve scatter plot (TGraph)
            TGraph *graph = dynamic_cast<TGraph *>(f->Get(graphName.c_str()));
            if (!graph)
            {
                std::cout << "Graph " << graphName << " not found, skipping..." << std::endl;
                continue;
            }

            // Get X and Y axis ranges
            double xMin = 0., xMax = 0.55;
            double yMin = 0., yMax = 0.55;
            //graph->ComputeRange(xMin, yMin, xMax, yMax);

            // Create a 2D histogram for visualization
            int nBinsX = 200, nBinsY = 200;
            TH2D *h2 = new TH2D(Form("hist2D_Crystal%d_vs_Crystal%d", i, j),
                                Form("2D Histogram: Crystal %d vs. Crystal %d;Charge Crystal %d;Charge Crystal %d",
                                     i, j, i, j),
                                nBinsX, xMin, xMax, nBinsY, yMin, yMax);

            double maxX = std::numeric_limits<double>::lowest();
            double maxY = std::numeric_limits<double>::lowest();

            double thresh_x, thresh_y;

            if (energy == 100) {
                thresh_x = 0.1;
                thresh_y = 0.2;
            }
            else if (energy == 140) {
                thresh_x = 0.2;
                thresh_y = 0.35;
            }
            else if (energy == 180) {
                thresh_x = 0.27;
                thresh_y = 0.4;
            }
            else if (energy == 200) {
                thresh_x = 0.55;
                thresh_y = 0.55;
            }
            else if (energy == 220) {
                thresh_x = 0.35;
                thresh_y = 0.4;
            }

            // Fill the 2D histogram from the TGraph
            for (int k = 0; k < graph->GetN(); ++k)
            {
                double x, y;
                graph->GetPoint(k, x, y);
                h2->Fill(x, y);
                if (x > maxX && x < thresh_x)
                    maxX = x;
                if (y > maxY && y < thresh_y)
                    maxY = y;
            }

            cout << "maxX: " << maxX << " maxY: " << maxY << endl;

            // Fit the scatter plot
            auto [fitRes, filteredGraph] = FitScatterPlot(graph, i, j, maxX, maxY);
            if (!fitRes.success)
            {
                continue;
            }
            else
            {
                std::cout << "Energy: " << energy
                          << ", ids: " << i << ", " << j
                          << ", p0: " << fitRes.p0 << " ± " << fitRes.p0err
                          << ", p1: " << fitRes.p1 << " ± " << fitRes.p1err
                          << std::endl
                          << std::endl;
            }

            // Create canvas and draw histogram
            TCanvas *c1 = new TCanvas(Form("c1_%s", graphName.c_str()), "2D Histogram and Fit", 800, 600);
            TCanvas *c2 = new TCanvas(Form("c2_%s", graphName.c_str()), "Scatter Plot and Fit", 800, 600);

            TF1 *fitFunc = new TF1(Form("fit_Crystal%d_vs_Crystal%d", i, j), "[0] + [1]*x", 0.02, maxX);
            fitFunc->SetParameters(fitRes.p0, fitRes.p1);
            fitFunc->SetLineColor(kRed);
            fitFunc->SetLineWidth(2);

            TLegend *leg = new TLegend(0.15, 0.75, 0.45, 0.85);
            leg->AddEntry(fitFunc, Form("y = %.2f + (%.2f x)", fitRes.p0, fitRes.p1), "l");

            filteredGraph->SetName(Form("filteredGraph_Crystal%d_vs_Crystal%d", i, j));
            filteredGraph->SetTitle(Form("filtered Graph: Crystal %d vs. Crystal %d",
                i, j, i, j));
            filteredGraph->SetMarkerStyle(20);
            filteredGraph->SetMarkerSize(0.35);
            filteredGraph->SetMarkerColor(kRed); // Red for filtered points
            filteredGraph->SetLineStyle(0);
            filteredGraph->SetLineWidth(0);
            filteredGraph->SetMinimum(0.);
            filteredGraph->SetMaximum(0.55);
            filteredGraph->GetXaxis()->SetLimits(0., 0.55);

            graph->SetMinimum(0.);
            graph->SetMaximum(0.55);
            graph->GetXaxis()->SetLimits(0., 0.55);
            graph->GetXaxis()->SetTitle(Form("Charge crystal ID %d [a.u.]", i));
            graph->GetYaxis()->SetTitle(Form("Charge crystal ID %d [a.u.]", j));
            graph->SetTitle(Form("Charge Scatter Plot (cluster size = 2): Crystal %d vs Crystal %d @ %d MeV/u", i, j, energy));

            fout->cd();

            c1->cd();
            gPad->SetLogz(1);
            gStyle->SetPalette(1);
            h2->GetXaxis()->SetTitle(Form("Charge crystal ID %d [a.u.]", i));
            h2->GetYaxis()->SetTitle(Form("Charge crystal ID %d [a.u.]", j));
            h2->SetMaximum(2000);
            h2->Draw("COLZ"); // "COLZ" for colored 2D plot
            fitFunc->Draw("same");
            leg->Draw();
            c1->Write();

            c2->cd();
            graph->Draw("AP");
            filteredGraph->Draw("P same");
            fitFunc->Draw("same");
            leg->AddEntry(filteredGraph, "Points considered in fit (red)");
            leg->Draw();
            c2->Write();
            

            delete c1;  // Cleanup
            delete c2;
            delete h2; // Cleanup
            delete fitFunc;
        }
    }

    // Close files
    fout->Close();
    f->Close();
}

// Function to fit scatter plots and return fit results
std::pair<FitResult, TGraph*> FitScatterPlot(TGraph *graph, int i, int j, double maxX, double maxY)
{
    FitResult result = {false, 0.0, 0.0, 0.0, 0.0};

    if (!graph)
        return {result, nullptr};

    std::vector<double> xFiltered, yFiltered;

    for (int i = 0; i < graph->GetN(); ++i) {
        double x, y;
        graph->GetPoint(i, x, y);

        // Apply selection condition
        if ((y > (maxY - 0.01 - x * (maxY / maxX))) && (y > 0.02)) {
            xFiltered.push_back(x);
            yFiltered.push_back(y);
        }
    }
    // Create a new graph with filtered points
    TGraph *filteredGraph = new TGraph(xFiltered.size(), xFiltered.data(), yFiltered.data());

    cout << "Intercept init. value: " << maxY << " Slope init. value: " << ((maxX != 0) ? -maxY / maxX : 0) << endl;

    // Define a linear fit function
    TF1 *fitFunc = new TF1(Form("fit_Crystal%d_vs_Crystal%d", i, j), "[0] + [1]*x", 0., maxX);
    fitFunc->SetParameter(0, maxY);
    fitFunc->SetParameter(1, (maxX != 0) ? -maxY / maxX : 0);
    filteredGraph->Fit(fitFunc, "QWR"); // "Q" for quiet mode

    result.success = true;
    result.p0 = fitFunc->GetParameter(0);
    result.p1 = fitFunc->GetParameter(1);
    result.p0err = fitFunc->GetParError(0);
    result.p1err = fitFunc->GetParError(1);

    return {result, filteredGraph};
}


TH2 *Rebin2DHistogram(const TH2 *h2, int newXbins, int newYbins)
{
    if (!h2)
    {
        std::cerr << "Error: Input histogram is null!" << std::endl;
        return nullptr;
    }

    // Get the histogram range
    double xMin = h2->GetXaxis()->GetXmin();
    double xMax = h2->GetXaxis()->GetXmax();
    double yMin = h2->GetYaxis()->GetXmin();
    double yMax = h2->GetYaxis()->GetXmax();

    // Create a new rebinned histogram
    TString newName = h2->GetName();
    newName += "_rebinned";
    TH2D *h2Rebinned = new TH2D(newName, h2->GetTitle(),
                                newXbins, xMin, xMax, newYbins, yMin, yMax);

    // Fill the new histogram with content from the original one
    for (int i = 1; i <= h2->GetNbinsX(); ++i)
    {
        for (int j = 1; j <= h2->GetNbinsY(); ++j)
        {
            double content = h2->GetBinContent(i, j);
            double x = h2->GetXaxis()->GetBinCenter(i);
            double y = h2->GetYaxis()->GetBinCenter(j);
            h2Rebinned->Fill(x, y, content);
        }
    }

    return h2Rebinned;
}
