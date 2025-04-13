#include "CompareCalibrated.h"

void CompareCalibrated() {

    // Vectors with data and MC files paired with their corresponding energies.
    std::vector<std::pair<std::string, int>> filesAndEnergiesRaw = {
        {"TW/cuts/AnaFOOT_TW_Decoded_HIT2022_fragm_100MeV.root", 100},
        {"TW/cuts/AnaFOOT_TW_Decoded_HIT2022_fragm_140MeV.root", 140},
        {"TW/cuts/AnaFOOT_TW_Decoded_HIT2022_fragm_200MeV.root", 200},
        {"TW/cuts/AnaFOOT_TW_Decoded_HIT2022_fragm_220MeV.root", 220}
    };

    std::vector<std::pair<std::string, int>> filesAndEnergies = {
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_fragm_100MeV.root", 100},
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_fragm_140MeV.root", 140},
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_fragm_200MeV.root", 200},
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_fragm_220MeV.root", 220}
    };

    std::vector<std::pair<std::string, int>> filesAndEnergiesMC = {
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_100.root", 100},
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_140.root", 140},
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_200.root", 200},
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_220.root", 220}
    };

    std::map<std::string, std::map<int, double>> calibCoeff = extractBarData();

    // Define the two layers.
    std::vector<std::string> layers = {"X", "Y"};

    // Loop over layers and bars (0 to 19).
    for (const auto &layer : layers) {
        for (int bar = 0; bar < 20; bar++) {
            // Loop over each energy (assumed to match between data and MC).
            for (size_t i = 0; i < filesAndEnergies.size(); i++) {
                int energy = filesAndEnergies[i].second;
                std::string dataFileRaw = filesAndEnergiesRaw[i].first;
                std::string dataFile = filesAndEnergies[i].first;
                std::string mcFile   = filesAndEnergiesMC[i].first;

                std::map<int, std::map<int, double>> tofCoeff = extractTofData(energy);
                cout << "Processing energy: " << energy << " MeV/u Layer " << layer.c_str() << " bar " << bar << endl;

                cout << "Calibrated eloss = Q*1/p0, 1/p0 = " << calibCoeff.at(layer).at(bar) << " MeV/a.u." << endl;
                cout << "Calibrated tof = raw tof - Delta_t, Delta_t = " << tofCoeff.at(0).at(bar) << " ns" << endl;
                
                // Process the file combination: retrieve histograms and draw them.
                ProcessFile(dataFileRaw, dataFile, mcFile, energy, layer, bar);
            }
        }
    }
}


// New helper function: getHistogramByBar
// Opens the file, goes to the given directory, and iterates over keys
// to find a histogram whose name exactly equals 'desiredHistName'.
// Only histograms whose names start with "noCuts_Eloss" are considered.
TH1D* getElossByBar(TFile* file, const std::string &dirName, const TString &desiredHistName) {
    // Change to the desired directory.
    TDirectory* dir = dynamic_cast<TDirectory*>(file->Get(dirName.c_str()));
    if (!dir) {
        std::cerr << "Directory " << dirName << " not found in file " << file->GetName() << std::endl;
        file->Close();
        return nullptr;
    }
    
    // Iterate over all keys in the directory.
    TIter next(dir->GetListOfKeys());
    TKey *key;
    TH1D* foundHist = nullptr;
    while ((key = (TKey*) next())) {
        TObject *obj = key->ReadObj();
        TH1D* hist = dynamic_cast<TH1D*>(obj);
        if (!hist) {
            delete obj;
            continue;
        }
        TString hName = hist->GetName();
        // Check if this is the histogram we are looking for.
        if (hName == desiredHistName) {
            foundHist = dynamic_cast<TH1D*>(hist->Clone());
            delete hist;
            break;
        }
        delete hist;
    }
    return foundHist;
}

TH1D* getTofByBar(TFile* file, const std::string &dirName, const TString &desiredHistName) {
    // Change to the desired directory.
    TDirectory* dir = dynamic_cast<TDirectory*>(file->Get(dirName.c_str()));
    if (!dir) {
        std::cerr << "Directory " << dirName << " not found in file " << file->GetName() << std::endl;
        file->Close();
        return nullptr;
    }
    
    // Iterate over all keys in the directory.
    TIter next(dir->GetListOfKeys());
    TKey *key;
    TH1D* foundHist = nullptr;
    while ((key = (TKey*) next())) {
        TObject *obj = key->ReadObj();
        TH1D* hist = dynamic_cast<TH1D*>(obj);
        if (!hist) {
            delete obj;
            continue;
        }
        TString hName = hist->GetName();
        // Check if this is the histogram we are looking for.
        if (hName == desiredHistName) {
            foundHist = dynamic_cast<TH1D*>(hist->Clone());
            delete hist;
            break;
        }
        delete hist;
    }
    return foundHist;
}

// Function to draw two normalized histograms on a canvas.
// hData is drawn in red and hMC in blue.
void ElossNormalizedHistograms(TH1* hDataRaw, TH1* hData, TH1* hMC, const TString &ElossTitle, const TString &ElossCanvasName) {
    TCanvas *c = new TCanvas(ElossCanvasName.Data(), ElossTitle.Data(), 800, 600);
    if (!hData || !hMC) {
        std::cerr << "One or both histograms are null for canvas " << ElossCanvasName.Data() << std::endl;
        delete c;
        return;
    }
    c->SetMargin(0.1, 0.1, 0.15, 0.15);

    hDataRaw->SetLineColor(kBlack);
    hData->SetLineColor(kBlue);
    hMC->SetLineColor(kRed);

    gPad->SetLogy();
    gStyle->SetOptStat(0);

    hMC->Scale(hData->GetMaximum() / hMC->GetMaximum());
    hDataRaw->Scale(hData->GetMaximum() / hDataRaw->GetMaximum());

    double maxY = std::max({hDataRaw->GetMaximum(), hData->GetMaximum(), hMC->GetMaximum()});
    maxY *= 1.5;

    double minY = std::min({hDataRaw->GetMinimum(0), hData->GetMinimum(0), hMC->GetMinimum(0)});
    //minY *= 0.5;

    hData->GetYaxis()->SetRangeUser(minY, maxY);
    hData->GetXaxis()->SetTitle("Energy loss [MeV] - Q [a.u.]");
    hData->GetYaxis()->SetTitle("Entries");
    hData->SetTitleSize(0.1, "T");
    hData->GetXaxis()->SetTitleSize(0.05);
    hData->GetYaxis()->SetTitleSize(0.05);

    hData->SetLineWidth(2);
    hDataRaw->SetLineWidth(2);
    hMC->SetLineWidth(2);
    
    // Draw histograms normalized (without modifying the original bin contents)
    hData->Draw();
    hDataRaw->Draw("HIST SAME");
    hMC->Draw("HIST SAME");
    
    // Build a legend to distinguish the histograms
    TLegend *leg = new TLegend(0.65, 0.7, 0.95, 0.85);
    leg->AddEntry(hData, "SHOE Calib. Eloss", "l");
    leg->AddEntry(hDataRaw, "Raw Eloss #sqrt{Q_{A}Q_{B}}", "l");
    leg->AddEntry(hMC, "MC", "l");
    leg->SetTextSize(0.04);
    leg->Draw();

    TString plotName = ElossCanvasName;
    plotName.ReplaceAll("c_Eloss", "Normalized_Eloss");
    
    c->Update();
    TString outputPath = "Plots/calibrated/" + plotName + ".png";
    c->SaveAs(outputPath.Data());
    delete c;
}


void TofNormalizedHistograms(TH1* hTofDataRaw, TH1* hTofData, TH1* hTofMC, const TString &TofTitle, const TString &TofCanvasName) {
    TCanvas *c = new TCanvas(TofCanvasName.Data(), TofTitle.Data(), 800, 600);
    c->SetMargin(0.1, 0.1, 0.15, 0.15); // Left, Right, Bottom, Top margins
    //c->SetGrid();
    if (!hTofData || !hTofMC) {
        std::cerr << "One or both histograms are null for canvas " << TofCanvasName.Data() << std::endl;
        delete c;
        return;
    }

    hTofDataRaw->SetLineColor(kBlack);
    hTofData->SetLineColor(kBlue);
    hTofMC->SetLineColor(kRed);

    gPad->SetLogy();
    gStyle->SetOptStat(0);

    hTofMC->Scale(hTofData->GetMaximum() / hTofMC->GetMaximum());
    hTofDataRaw->Scale(hTofData->GetMaximum() / hTofDataRaw->GetMaximum());

    double maxY = std::max({hTofDataRaw->GetMaximum(), hTofData->GetMaximum(), hTofMC->GetMaximum()});
    maxY *= 1.5;

    double minY = std::min({hTofDataRaw->GetMinimum(0), hTofData->GetMinimum(0), hTofMC->GetMinimum(0)});
    //minY *= 0.5;

    hTofData->GetYaxis()->SetRangeUser(minY, maxY);
    hTofData->GetXaxis()->SetTitle("TOF [ns]");
    hTofData->GetYaxis()->SetTitle("Entries");
    gStyle->SetTitleSize(0.08, "T");

    hTofData->GetXaxis()->SetTitleSize(0.05);
    hTofData->GetYaxis()->SetTitleSize(0.05);
    // Draw histograms normalized (without modifying the original bin contents)
    hTofData->SetLineWidth(2);
    hTofDataRaw->SetLineWidth(2);
    hTofMC->SetLineWidth(2);

    hTofData->Draw();
    hTofDataRaw->Draw("HIST SAME");
    hTofMC->Draw("HIST SAME");
    
    // Build a legend to distinguish the histograms
    TLegend *leg = new TLegend(0.65, 0.7, 0.95, 0.85);
    leg->AddEntry(hTofData, "SHOE Calib. TOF", "l");
    leg->AddEntry(hTofDataRaw, "Raw TOF T_{bar} - T_{SC}", "l");
    leg->AddEntry(hTofMC, "MC", "l");
    leg->SetTextSize(0.04);
    leg->Draw();

    TString plotName = TofCanvasName;
    plotName.ReplaceAll("c_Tof", "Normalized_ToF");
    
    c->Update();
    TString outputPath = "Plots/calibrated/" + plotName + ".png";
    c->SaveAs(outputPath.Data());
    delete c;
}


void ProcessFile(const std::string &dataFileRaw, const std::string &dataFile, const std::string &mcFile, int energy,
                 const std::string &layer, int bar) {
    // Build directory names based on the layer.
    TString dataDir = TString::Format("ChargeElossLayer%s", layer.c_str());
    TString dataRawDir = TString::Format("ChargeTimeLayer%s", layer.c_str());
    TString mcDir   = TString::Format("ElossLayer%s", layer.c_str());
    TString tofDirData = TString::Format("ToFLayer%s", layer.c_str());
    TString tofDirMC = TString::Format("TofLayer%s", layer.c_str());
    // Construct histogram name.
    TString histNameMyEloss = TString::Format("MyEloss_Layer%s_bar%d", layer.c_str(), bar);
    TString histNameEloss = TString::Format("Eloss_Layer%s_bar%d", layer.c_str(), bar);
    TString histNameElossRaw = TString::Format("Charge_Layer%s_bar%d", layer.c_str(), bar);
    TString histNameTof = TString::Format("ToF_Layer%s_bar%d", layer.c_str(), bar);

    
    TFile *dataRawFilePtr = TFile::Open(dataFileRaw.c_str(), "READ");
    if (!dataRawFilePtr || dataRawFilePtr->IsZombie()) {
        std::cerr << "Error opening data file: " << dataFileRaw << std::endl;
        dataRawFilePtr->Close();
        delete dataRawFilePtr;
        return;
    }

    // Open the ROOT files outside the getHistogramByBar function
    TFile *dataFilePtr = TFile::Open(dataFile.c_str(), "READ");
    if (!dataFilePtr || dataFilePtr->IsZombie()) {
        std::cerr << "Error opening data file: " << dataFile << std::endl;
        dataFilePtr->Close();
        delete dataFilePtr;
        return;
    }

    TFile *mcFilePtr = TFile::Open(mcFile.c_str(), "READ");
    if (!mcFilePtr || mcFilePtr->IsZombie()) {
        std::cerr << "Error opening MC file: " << mcFile << std::endl;
        mcFilePtr->Close();
        delete mcFilePtr;
        return;
    }

    // Retrieve the histograms from the opened files
    TH1D* hData = getElossByBar(dataFilePtr, dataDir.Data(), histNameEloss.Data());
    TH1D* hDataRaw = getElossByBar(dataRawFilePtr, dataRawDir.Data(), histNameElossRaw.Data());
    TH1D* hMC = getElossByBar(mcFilePtr, mcDir.Data(), histNameEloss.Data());

    TH1D *hTofData = getTofByBar(dataFilePtr, tofDirData.Data(), histNameTof.Data());
    TH1D *hTofDataRaw = getTofByBar(dataRawFilePtr, tofDirData.Data(), histNameTof.Data());
    TH1D *hTofMC = getTofByBar(mcFilePtr, tofDirMC.Data(), histNameTof.Data());
    
    // Create a title and canvas name that include layer, bar, and energy.
    TString ElossTitle = TString::Format("Raw, calibrated and MC Eloss layer %s bar %d @ %d MeV/u", layer.c_str(), bar, energy);
    TString ElossCanvasName = TString::Format("c_Eloss_Layer%s_bar%d_%dMeV", layer.c_str(), bar, energy);
    TString TofTitle = TString::Format("Raw, calibrated and MC TOF layer %s bar %d @ %d MeV/u", layer.c_str(), bar, energy);
    TString TofCanvasName = TString::Format("c_Tof_Layer%s_bar%d_%dMeV", layer.c_str(), bar, energy);

    hData->SetTitle(ElossTitle.Data());
    hTofData->SetTitle(TofTitle.Data());
    
    // Only proceed if both histograms were successfully retrieved.
    if (!hData || !hDataRaw || !hMC || !hTofData || !hTofDataRaw || !hTofMC) {
        std::cerr << "Skipping Layer " << layer << " Bar " << bar << " Energy " << energy << " due to missing histogram(s)." << std::endl;
        dataFilePtr->Close();  // Close the files if histograms were not retrieved
        dataRawFilePtr->Close();
        mcFilePtr->Close();
        delete dataFilePtr;
        delete dataRawFilePtr;
        delete mcFilePtr;
        return;
    }

    // Draw the normalized histograms.
    cout << "Drawing histograms" << endl;
    ElossNormalizedHistograms(hDataRaw, hData, hMC, ElossTitle, ElossCanvasName);
    TofNormalizedHistograms(hTofDataRaw, hTofData, hTofMC, TofTitle, TofCanvasName);
    // Clean up: delete the histograms
    delete hData;
    delete hDataRaw;
    delete hMC;
    delete hTofData;
    delete hTofDataRaw;
    delete hTofMC;

    // Close the ROOT files after processing
    dataRawFilePtr->Close();
    dataFilePtr->Close();
    mcFilePtr->Close();
    delete dataRawFilePtr;
    delete dataFilePtr;
    delete mcFilePtr;
}

std::map<std::string, std::map<int, double>> extractBarData() {
    std::map<std::string, std::map<int, double>> barData;
    std::string path = "../shoe/build/Reconstruction/calib/HIT2022/";
    std::string filename = "TATW_Energy_Calibration_perBar_4742.cal"; // File is hardcoded
    std::ifstream file(path + filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << path + filename << std::endl;
        return barData;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue; // Skip headers

        std::istringstream iss(line);
        int barId, p1, shoeLayer;
        double p0;

        if (!(iss >> barId >> p0 >> p1 >> shoeLayer)) continue; // Skip invalid lines

        std::string layer = (shoeLayer == 0) ? "Y" : "X";
        int correctedBar = (shoeLayer == 0) ? barId : barId - 20; // Adjust for X layer

        barData[layer][correctedBar] = p0;
    }

    file.close();
    return barData;
}

std::map<int, std::map<int, double>> extractTofData(int energy) {
    std::map<int, std::map<int, double>> tofData;
    std::string filename;
    std::string path = "../shoe/build/Reconstruction/calib/HIT2022/";
    if (energy == 100) filename = "TATW_Tof_Calibration_perBar_4766.cal";
    else if (energy == 140) filename = "TATW_Tof_Calibration_perBar_4801.cal";
    else if (energy == 200) filename = "TATW_Tof_Calibration_perBar_4742.cal";
    else if (energy == 220) filename = "TATW_Tof_Calibration_perBar_4828.cal";
    std::ifstream file(path + filename);
  
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << path + filename << std::endl;
        return tofData;
    }
  
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue; // Skip headers
  
        std::istringstream iss(line);
        int barId, shoeLayer;
        double Delta_t, sigma_t;
  
        if (!(iss >> barId >> Delta_t >> sigma_t >> shoeLayer)) continue; // Skip invalid lines
  
        int layer = shoeLayer;
        int correctedBar = (shoeLayer == 0) ? barId : barId - 20; // Adjust for X layer
  
        tofData[layer][correctedBar] = Delta_t;
    }
  
    file.close();
    return tofData;
  }

