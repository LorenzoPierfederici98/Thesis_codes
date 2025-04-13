#include <TFile.h>        
#include <TGraph.h>        
#include <TF1.h>       
#include <TCanvas.h>       
#include <iostream>

void DrawLinesInterCalib(int energy, int crystalID) 
{
    std::string fileName = Form("Calo/intercalib/AnaFOOT_Calo_Decoded_HIT2022_%dMeV.root", energy);
    TFile *file = TFile::Open(fileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        return;
    }
    std::string scatterPlotName = Form("scatter_Crystal0_vs_Crystal%d", crystalID);
    TGraph* scatterPlot = dynamic_cast<TGraph*>(file->Get(scatterPlotName.c_str()));
    if (!scatterPlot) {
        std::cerr << "Error: Scatter plot not found with name " << scatterPlotName << " in file " << fileName << std::endl;
        return;
    }

    int nLines = 4;  // number of lines to be plotted
    double slope_ratio;
    if (crystalID == 1)
    {
        slope_ratio = -1.69269;
    }
    else if (crystalID == 6)
    {
        slope_ratio = -1.63262;
    }
    double intercept = 0.13;
    
    TCanvas* canvas = new TCanvas("canvas", "Scatter Plot with Lines", 800, 600);
    canvas->SetMargin(0.15, 0.12, 0.15, 0.15);
    scatterPlot->Draw("AP");  // Draw the scatter plot with axis and points

    scatterPlot->GetXaxis()->SetRangeUser(0., 0.35);

    scatterPlot->GetXaxis()->SetTitle("Charge Crystal ID 0 [a.u.]");
    scatterPlot->GetYaxis()->SetTitle(Form("Charge Crystal ID %d [a.u.]", crystalID));
    scatterPlot->GetXaxis()->SetTitleSize(0.05);
    scatterPlot->GetYaxis()->SetTitleSize(0.05);
    scatterPlot->SetTitle(Form("Cluster size 2: crystal IDs 0 vs %d @ %d MeV/u", crystalID, energy));
    gStyle->SetTitleSize(0.07, "T");
    double thresh_x = (energy == 180) ? 0.26 : 0.29;
    //double thresh_y = 0.45;
    double thresh_y = (energy == 180) ? 0.36 : 0.45;

    double maxX = std::numeric_limits<double>::lowest();
    double maxY = std::numeric_limits<double>::lowest();

    for (int k = 0; k < scatterPlot->GetN(); k++)
    {
        double x, y;
        scatterPlot->GetPoint(k, x, y);
        if (x > maxX && x < thresh_x) maxX = x;
        if (y > maxY && y < thresh_y) maxY = y;
    }

    cout << "maxX: " << maxX << endl;
    cout << "maxY: " << maxY << endl;

    TLegend* legend = new TLegend(0.5, 0.6, 0.9, 0.8);

    for (int i = 0; i < nLines + 1; i++)
    {
        TF1 *line = new TF1(Form("line%d", i), "[0]*x + [1]", 0., 0.3);
        line->SetParameter(0, slope_ratio);
        line->SetParameter(1, intercept);
        line->SetLineColor(kBlack);
        line->SetLineWidth(2);
        line->Draw("same");
        intercept += 0.1;
        if (i == nLines - 1)
        {
            legend->AddEntry(line, Form("Inter-calib. slope = %.3f", slope_ratio), "l");
        }
    }

    double intercept_low = 0.35;
    double intercept_high = 0.4;
    if (energy == 200)
    {
        if (crystalID == 1)
        {
            intercept_low = 0.37;
            intercept_high = 0.43;
        }
        else if (crystalID == 6)
        {
            intercept_low = 0.43;
            intercept_high = 0.6;
        }
    }
    double slope = - maxY / maxX;

    TGraph* filteredGraph = new TGraph();

    // Building the filtered graph to fit around the upper diagonal line
    for (int k = 0; k < scatterPlot->GetN(); k++)
    {
        double x, y;
        scatterPlot->GetPoint(k, x, y);
        double y_low = slope * x + intercept_low;
        double y_high = slope * x + intercept_high;
        if (y > y_low && y < y_high)
        {
            filteredGraph->SetPoint(filteredGraph->GetN(), x, y);
        }

    }

    filteredGraph->SetMarkerColor(kRed);
    filteredGraph->SetMarkerStyle(20);
    filteredGraph->SetMarkerSize(0.3);
    filteredGraph->Draw("P SAME");
    legend->AddEntry(filteredGraph, "Points fitted (red)", "p");

    TF1 *line_fit = new TF1("line_fit", "[0]*x + [1]", 0., maxX);
    line_fit->SetParameter(0, - maxY / maxX);
    line_fit->SetParameter(1, (intercept_low + intercept_high) / 2.);
    //line_max->SetParameter(0, - 0.435 / maxX);
    //line_max->SetParameter(1, 0.46);
    line_fit->SetLineColor(kRed);
    line_fit->SetLineWidth(2);
    filteredGraph->Fit("line_fit", "R");
    legend->AddEntry(line_fit, Form("Fit line slope %.3f", line_fit->GetParameter(0)), "l");

    legend->SetTextSize(0.035);
    legend->Draw();

    TString outputFileName = TString(fileName).ReplaceAll(".root", "_Scatter.root");
    TFile* outputFile = TFile::Open(outputFileName, "UPDATE");
    outputFile->cd();
    canvas->Write(Form("c_Scatter_Line_Crystal0_Crystal%d", crystalID), TObject::kOverwrite);
    canvas->SaveAs(Form("Plots/ScatterFit_Crystal0_vs_Crystal%d_%dMeV.png", crystalID, energy));

    file->Close();
    outputFile->Close();
    delete file;
    delete outputFile;

}