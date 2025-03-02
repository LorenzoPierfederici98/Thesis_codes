#include "CompareCalibrated.h"

void CompareCalibrated() {

    // Vectors with data and MC files paired with their corresponding energies.
    std::vector<std::pair<std::string, int>> filesAndEnergies = {
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_fragm_100MeV.root", 100},
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_fragm_140MeV.root", 140},
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_fragm_200MeV.root", 200},
        {"calibrated/AnaFOOT_TW_Decoded_HIT2022_fragm_220MeV.root", 220}
    };

    std::vector<std::pair<std::string, int>> filesAndEnergiesMC = {
        {"calibrated/AnaFOOT_TW_DecodedMC_HIT2022_MC_100.root", 100},
        {"calibrated/AnaFOOT_TW_DecodedMC_HIT2022_MC_140.root", 140},
        {"calibrated/AnaFOOT_TW_DecodedMC_HIT2022_MC_200.root", 200},
        {"calibrated/AnaFOOT_TW_DecodedMC_HIT2022_MC_220.root", 220}
    };

    // Define the two layers.
    std::vector<std::string> layers = {"X", "Y"};

    // Loop over layers and bars (0 to 19).
    for (const auto &layer : layers) {
        for (int bar = 0; bar < 20; bar++) {
            // Loop over each energy (assumed to match between data and MC).
            for (size_t i = 0; i < filesAndEnergies.size(); i++) {
                int energy = filesAndEnergies[i].second;
                std::string dataFile = filesAndEnergies[i].first;
                std::string mcFile   = filesAndEnergiesMC[i].first;
                
                // Process the file combination: retrieve histograms and draw them.
                ProcessFile(dataFile, mcFile, energy, layer, bar);
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
        // Consider only histograms beginning with "noCuts_Eloss"
        if (!hName.BeginsWith("noCuts_Eloss")) {
            delete hist;
            continue;
        }
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
        // Consider only histograms beginning with "noCuts_Eloss"
        if (!hName.BeginsWith("noCuts_ToF")) {
            delete hist;
            continue;
        }
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
void ElossNormalizedHistograms(TH1* hData, TH1* hMC, const TString &ElossTitle, const TString &ElossCanvasName) {
    TCanvas *c = new TCanvas(ElossCanvasName.Data(), ElossTitle.Data(), 800, 600);
    if (!hData || !hMC) {
        std::cerr << "One or both histograms are null for canvas " << ElossCanvasName.Data() << std::endl;
        delete c;
        return;
    }

    hData->SetLineColor(kBlue);
    hMC->SetLineColor(kRed);

    gPad->SetLogy();
    
    // Draw histograms normalized (without modifying the original bin contents)
    hData->DrawNormalized();
    hMC->DrawNormalized("SAME");
    
    // Build a legend to distinguish the histograms
    TLegend *leg = new TLegend(0.68, 0.8, 0.78, 0.9);
    leg->AddEntry(hData, "Data", "l");
    leg->AddEntry(hMC, "MC", "l");
    leg->SetTextSize(0.03);
    leg->Draw();

    TString plotName = ElossCanvasName;
    plotName.ReplaceAll("c_Eloss", "Normalized_Eloss");
    
    c->Update();
    TString outputPath = "Plots/calibrated/" + plotName + ".png";
    c->SaveAs(outputPath.Data());
    delete c;
}


void TofNormalizedHistograms(TH1* hTofData, TH1* hTofMC, const TString &TofTitle, const TString &TofCanvasName) {
    TCanvas *c = new TCanvas(TofCanvasName.Data(), TofTitle.Data(), 800, 600);
    if (!hTofData || !hTofMC) {
        std::cerr << "One or both histograms are null for canvas " << TofCanvasName.Data() << std::endl;
        delete c;
        return;
    }

    hTofData->SetLineColor(kBlue);
    hTofMC->SetLineColor(kRed);

    gPad->SetLogy();
    
    // Draw histograms normalized (without modifying the original bin contents)
    hTofData->DrawNormalized();
    hTofMC->DrawNormalized("SAME");
    
    // Build a legend to distinguish the histograms
    TLegend *leg = new TLegend(0.68, 0.8, 0.78, 0.9);
    leg->AddEntry(hTofData, "Data", "l");
    leg->AddEntry(hTofMC, "MC", "l");
    leg->SetTextSize(0.03);
    leg->Draw();

    TString plotName = TofCanvasName;
    plotName.ReplaceAll("c_Tof", "Normalized_ToF");
    
    c->Update();
    TString outputPath = "Plots/calibrated/" + plotName + ".png";
    c->SaveAs(outputPath.Data());
    delete c;
}


void ProcessFile(const std::string &dataFile, const std::string &mcFile, int energy, const std::string &layer, int bar) {
    // Build directory names based on the layer.
    TString dataDir = TString::Format("ChargeElossLayer%s", layer.c_str());
    TString mcDir   = TString::Format("ElossLayer%s", layer.c_str());
    TString tofDirData = TString::Format("ToFLayer%s", layer.c_str());
    TString tofDirMC = TString::Format("TofLayer%s", layer.c_str());
    // Construct histogram name.
    TString histNameEloss = TString::Format("noCuts_Eloss_Layer%s_bar%d", layer.c_str(), bar);
    TString histNameTof = TString::Format("noCuts_ToF_Layer%s_bar%d", layer.c_str(), bar);

    cout << "Processing energy: " << energy << " MeV/u Layer " << layer.c_str() << " bar " << bar << endl;
    
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
    TH1D* hMC = getElossByBar(mcFilePtr, mcDir.Data(), histNameEloss.Data());

    TH1D *hTofData = getTofByBar(dataFilePtr, tofDirData.Data(), histNameTof.Data());
    TH1D *hTofMC = getTofByBar(mcFilePtr, tofDirMC.Data(), histNameTof.Data());
    
    // Create a title and canvas name that include layer, bar, and energy.
    TString ElossTitle = TString::Format("Calibrated eloss layer %s bar %d @ %d MeV/u", layer.c_str(), bar, energy);
    TString ElossCanvasName = TString::Format("c_Eloss_Layer%s_bar%d_%dMeV", layer.c_str(), bar, energy);
    TString TofTitle = TString::Format("Calibrated TOF layer %s bar %d @ %d MeV/u", layer.c_str(), bar, energy);
    TString TofCanvasName = TString::Format("c_Tof_Layer%s_bar%d_%dMeV", layer.c_str(), bar, energy);

    hData->SetTitle(ElossTitle.Data());
    hTofData->SetTitle(TofTitle.Data());
    
    // Only proceed if both histograms were successfully retrieved.
    if (!hData || !hMC || !hTofData || !hTofMC) {
        std::cerr << "Skipping Layer " << layer << " Bar " << bar << " Energy " << energy << " due to missing histogram(s)." << std::endl;
        dataFilePtr->Close();  // Close the files if histograms were not retrieved
        mcFilePtr->Close();
        delete dataFilePtr;
        delete mcFilePtr;
        return;
    }

    // Draw the normalized histograms.
    cout << "Drawing histograms" << endl;
    ElossNormalizedHistograms(hData, hMC, ElossTitle, ElossCanvasName);
    TofNormalizedHistograms(hTofData, hTofMC, TofTitle, TofCanvasName);
    // Clean up: delete the histograms
    delete hData;
    delete hMC;
    delete hTofData;
    delete hTofMC;

    // Close the ROOT files after processing
    dataFilePtr->Close();
    mcFilePtr->Close();
    delete dataFilePtr;
    delete mcFilePtr;
}

