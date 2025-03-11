// Macro that performs the ToF (Time of Flight) calibration by computing the differences
// between the ToF from data and ToF from MC mean fit values for every bar, given by the AnalyzeTofFragm.cc
// and AnalyzeTofMC.cc macros. These differences are written in 4 calibration files, one for each energy,
// along with the corresponding uncertainties (squared sum of the values).
// To be run with root -l -b -q 'CalibrateTof.cc()'.

#include "CalibrateTof.h"

void CalibrateTof() {
    std::vector<std::pair<std::string, int>> filesAndEnergies = {
        {"TW/cuts/AnaFOOT_TW_Decoded_HIT2022_fragm_100MeV_Fit.root", 100},
        {"TW/cuts/AnaFOOT_TW_Decoded_HIT2022_fragm_140MeV_Fit.root", 140},
        {"TW/cuts/AnaFOOT_TW_Decoded_HIT2022_fragm_200MeV_Fit.root", 200},
        {"TW/cuts/AnaFOOT_TW_Decoded_HIT2022_fragm_220MeV_Fit.root", 220}};

    std::vector<std::pair<std::string, int>> filesAndEnergiesMC = {
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_100_Fit.root", 100},
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_140_Fit.root", 140},
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_200_Fit.root", 200},
        {"MC/TW/AnaFOOT_TW_DecodedMC_HIT2022_MC_220_Fit.root", 220}};

    std::map<int, std::vector<std::string>> calibFiles;
    calibFiles[100] = {"4766", "4884"};
    calibFiles[140] = {"4801", "4877"};
    calibFiles[200] = {"4742", "4868"};
    calibFiles[220] = {"4828"};

    for (const std::string& layer : {"Y", "X"}) {
        for (int bar = 0; bar < 20; ++bar) {
            WriteMeanDifferences(calibFiles, filesAndEnergies, filesAndEnergiesMC, layer, bar);
        }
    }
}

// Function to load a fit result from the TofFit directory in a file
TFitResult* LoadFitResult(const std::string& fileName, const std::string& fitResultName) {
    TFile file(fileName.c_str(), "READ");
    if (file.IsZombie()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return nullptr;
    }
    TDirectory* tofFitDir = file.GetDirectory("TofFit");
    if (!tofFitDir) {
        std::cerr << "Error: TofFit directory not found in file: " << fileName << std::endl;
        file.Close();
        return nullptr;
    }
    tofFitDir->cd();  // Change directory context to TofFit
    TFitResult* fitResult = (TFitResult*)tofFitDir->Get(fitResultName.c_str());
    if (!fitResult) {
        std::cerr << "Error retrieving fit result: " << fitResultName << std::endl;
    }
    file.Close();
    return fitResult;
}

// Function to extract mean values and uncertainties from fit results
void ExtractMeanValues(const std::vector<std::pair<std::string, int>>& filesAndEnergies, const std::string& layer, int bar, std::vector<double>& means, std::vector<double>& meanErrors) {
    for (const auto& fileEnergyPair : filesAndEnergies) {
        std::string fitResultName = "fitResultTof_Layer" + layer + "_bar" + std::to_string(bar);
        TFitResult* fitResult = LoadFitResult(fileEnergyPair.first, fitResultName);
        if (fitResult) {
            double meanValue = fitResult->Parameter(1); // Assuming the mean value is the second parameter
            double meanError = fitResult->ParError(1); // Assuming the mean error is the error of the second parameter
            means.push_back(meanValue);
            meanErrors.push_back(meanError);
        }
    }
}

void WriteMeanDifferences(const std::map<int, std::vector<std::string>> calibFiles,
                          const std::vector<std::pair<std::string, int>>& filesAndEnergies,
                          const std::vector<std::pair<std::string, int>>& filesAndEnergiesMC,
                          const std::string& layer, int bar) {
    std::vector<double> meansData, meanErrorsData;
    std::vector<double> meansMC, meanErrorsMC;
    
    // Extract mean values and uncertainties from the data files
    ExtractMeanValues(filesAndEnergies, layer, bar, meansData, meanErrorsData);
    // Extract mean values and uncertainties from the MC files
    ExtractMeanValues(filesAndEnergiesMC, layer, bar, meansMC, meanErrorsMC);
    
    if (meansData.size() != meansMC.size()) {
        std::cerr << "Mismatch in number of points between data and MC for Layer " 
                  << layer << " Bar " << bar << std::endl;
        return;
    }
    
    // For this layer and bar, determine the combined bar ID and layer indicator.
    int barID = bar;  
    bool isLayerX = (layer == "X");
    int combinedBarID = barID + (isLayerX ? 20 : 0);
    int layerShoe = isLayerX ? 1 : 0;
    
    // Loop over each energy point (assuming same ordering in both vectors)
    for (size_t i = 0; i < meansData.size(); ++i) {
        // Compute difference: data mean minus MC mean
        double diff = meansData[i] - meansMC[i];
        // Compute uncertainty on the difference assuming independent errors (quadrature sum)
        double diffError = std::sqrt(meanErrorsData[i]*meanErrorsData[i] + meanErrorsMC[i]*meanErrorsMC[i]);
        
        // Use the energy value from the filesAndEnergies vector
        int energy = filesAndEnergies[i].second;

        std::cout << "Layer " << layer << " Bar " << bar 
                  << ", Energy " << energy << " MeV/u:" << std::endl;
        std::cout << "   Data tof mean [ns] = " << meansData[i] 
                  << "   MC tof mean [ns] = " << meansMC[i] << std::endl;
        std::cout << "   Difference (Data - MC) [ns] = " << diff 
                  << " +- " << diffError << std::endl << std::endl;
        
        for (size_t i = 0; i < calibFiles.at(energy).size(); i++)
        {
            // Build the output file name (one file per energy)
            std::string outputFileName = "TATW_Tof_Calibration_perBar_" + calibFiles.at(energy)[i] + ".cal";
            
            // Determine whether to write the header.
            // We do this by checking if the file already exists.
            bool writeHeader = false;
            {
                std::ifstream infile(outputFileName);
                if (!infile.good()) {
                    writeHeader = true;
                }
                infile.close();
            }
            
            // Open the file in append mode.
            std::ofstream outFile;
            outFile.open(outputFileName, std::ios::app);
            if (!outFile.is_open()) {
                std::cerr << "Error: Unable to open file " << outputFileName << std::endl;
                continue;
            }
            
            // Write header if needed.
            if (writeHeader) {
                outFile << "#BarId(Pisa)          DeltaT                σ(DeltaT)        layer(SHOE)" << std::endl;
                outFile << "#-+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-" << std::endl;
            }
            
            // Write the data line
            outFile << std::setw(10) << combinedBarID
                    << std::setw(20) << diff
                    << std::setw(19) << diffError
                    << std::setw(10) << layerShoe
                    << std::endl;
            
            // If this is the last combination (Layer X and bar 19), append the footer line once.
            if ((layer == "X") && bar == 19) {
                outFile << "#-+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-" << std::endl;
            }
            outFile.close();
        }
    }
}
