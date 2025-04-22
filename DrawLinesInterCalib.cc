#include <TFile.h>        
#include <TGraph.h>        
#include <TF1.h>       
#include <TCanvas.h>       
#include <iostream>

std::pair<double, double> RoundMeasurement(double value, double uncertainty);

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

    // int nLines = 4;  // number of lines to be plotted
    // double slope_ratio;
    // if (crystalID == 1)
    // {
    //     slope_ratio = -1.69269;
    // }
    // else if (crystalID == 6)
    // {
    //     slope_ratio = -1.63262;
    // }
    // double intercept = 0.13;
    
    TCanvas* canvas = new TCanvas("canvas", "Scatter Plot with Lines", 800, 600);
    canvas->SetMargin(0.15, 0.12, 0.15, 0.15);
    //scatterPlot->SetMarkerSize(0.2);
    scatterPlot->Draw("AP");  // Draw the scatter plot with axis and points

    scatterPlot->GetXaxis()->SetRangeUser(0., 0.35);

    scatterPlot->GetXaxis()->SetTitle("Charge Crystal ID 0 [a.u.]");
    scatterPlot->GetYaxis()->SetTitle(Form("Charge Crystal ID %d [a.u.]", crystalID));
    scatterPlot->GetXaxis()->SetTitleSize(0.05);
    scatterPlot->GetYaxis()->SetTitleSize(0.05);
    scatterPlot->SetTitle(Form("Single clusters of size 2: fit crystal IDs %d vs 0 @ %d MeV/u", crystalID, energy));
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

    TLegend* legend = new TLegend(0.45, 0.54, 0.9, 0.8);

    // for (int i = 0; i < nLines + 1; i++)
    // {
    //     TF1 *line = new TF1(Form("line%d", i), "[0]*x + [1]", 0., 0.3);
    //     line->SetParameter(0, slope_ratio);
    //     line->SetParameter(1, intercept);
    //     line->SetLineColor(kBlack);
    //     line->SetLineWidth(2);
    //     line->Draw("same");
    //     intercept += 0.1;
    //     if (i == nLines - 1)
    //     {
    //         legend->AddEntry(line, Form("Inter-calib. slope = %.3f", slope_ratio), "l");
    //     }
    // }

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
    double m_err_noRound = line_fit->GetParError(0);
    double m_noRound = line_fit->GetParameter(0);
    auto [m, m_err] = RoundMeasurement(m_noRound, m_err_noRound);
    auto [q, q_err] = RoundMeasurement(line_fit->GetParameter(1), line_fit->GetParError(1));
    cout << "Energy: " << energy << " MeV/u " << " crystal ID " << crystalID << " vs 0" << endl;
    cout << "Upper diagonal fit line: m = " << m << " +/- " << m_err << " q = " << q << " +/- " << q_err << endl; 
    legend->AddEntry(line_fit, Form("m_{%d0}^{u} = %.3f#pm %.3f", crystalID, m, m_err), "l");

    // Define slope and point;
    double x0 = 0.1;
    double y0 = 0.2;
    double q_parall = y0 - m * x0;

    // Define the line
    TF1* line_parallel = new TF1("line", "[0]*x + [1]", 0., maxX);
    line_parallel->SetParameter(0, m);
    line_parallel->SetParameter(1, q_parall);
    line_parallel->SetLineColor(kBlack);
    line_parallel->SetLineWidth(2);
    line_parallel->Draw("SAME");  // you can use "SAME" if on top of something else

    // Compute intersection with x-axis
    double x_intersect_noRound = -q_parall / m;
    double x_intersect_err_noRound = y0 * m_err_noRound / (m_noRound * m_noRound);

    auto [x_intersect, x_intersect_err] = RoundMeasurement(x_intersect_noRound, x_intersect_err_noRound);

    cout << "Equalized Q = " << x_intersect << " +/- " << x_intersect_err << endl;

    // Create markers
    TMarker* pointMarker = new TMarker(x0, y0, 20);
    pointMarker->SetMarkerColor(kRed);
    pointMarker->SetMarkerSize(1.5);
    pointMarker->Draw("SAME");

    TMarker* intersectMarker = new TMarker(x_intersect, 0.0, 20);
    intersectMarker->SetMarkerColor(kGreen + 2);
    intersectMarker->SetMarkerSize(1.5);
    intersectMarker->Draw("SAME");

    legend->AddEntry(line_parallel, "Line parallel to the upp. diagonal", "l");
    legend->AddEntry(pointMarker, Form("Q_{%d} = %.3f", crystalID, y0), "p");
    legend->AddEntry(intersectMarker, Form("Q_{%d,eq} = %.4f#pm %.4f", crystalID, x_intersect, x_intersect_err), "p");

    legend->SetTextSize(0.035);
    legend->Draw("SAME");

    TString outputFileName = TString(fileName).ReplaceAll(".root", "_Scatter.root");
    TFile* outputFile = TFile::Open(outputFileName, "UPDATE");
    outputFile->cd();
    canvas->Write(Form("c_Scatter_Line_Crystal%d_Crystal0", crystalID), TObject::kOverwrite);
    canvas->SaveAs(Form("Plots/ScatterFit_Crystal%d_vs_Crystal0_%dMeV.png", crystalID, energy));

    file->Close();
    outputFile->Close();
    delete file;
    delete outputFile;

}

std::pair<double, double> RoundMeasurement(double value, double uncertainty) {
    int significantFigures = (int)std::ceil(-std::log10(uncertainty));
    double roundingFactor = std::pow(10, significantFigures);
    double roundedUncertainty = std::round(uncertainty * roundingFactor) / roundingFactor;
    double roundedValue = std::round(value * roundingFactor) / roundingFactor;
    return {roundedValue, roundedUncertainty};
}