// written by mtoppi 05/2023
// to be run in shoe/build/Reconstruction
// [ therein set-up the enviroment:  source ../setupFOOT.sh ]
// and run with (for example): root -l -b -q AnalyzeFOOT.cc++g\(\"../../../../rootfiles/outMC_16O_C_400_1_GSI.root\",1,10,\"\"\)
// sul tier1:  root -l -b -q AnalyzeFOOT.cc++g\(\"/storage/gpfs_data/foot/mtoppi/DataDecoded/CNAO2023/test.root\",0,1000,\"testAnaFOOT\",\"/storage/gpfs_data/foot/mtoppi/OutputMacro/\"\)

#include "AnalyzeCaloClustersInterCalib.h"

// main
void AnalyzeCaloClustersInterCalib(TString infile = "testMC.root", Bool_t isMax = kFALSE, Int_t nev = 10, TString outfile = "AnaFOOT.root", TString outDir = "/Users/marco/FOOT/Analisi/shoe/build/Reconstruction/OutputMacro")

{

    // InitializeContainers();

    TAGroot gTAGroot;

    TFile *inputFile = new TFile(infile.Data());

    TString runInfoName = "runinfo";

    for (TKey *key : ROOT::RangeStaticCast<TKey *>(*inputFile->GetListOfKeys()))
    {
        if (debug)
            std::cout << "key: " << key->GetName() << " points to an object of class: " << key->GetClassName() << '\n';

        TString keyClassName = key->GetClassName();
        if (!keyClassName.CompareTo("TAGrunInfo"))
            runInfoName = key->GetName();
    }

    cout << "runInfoName::" << runInfoName << endl;

    runinfo = (TAGrunInfo *)(inputFile->Get(runInfoName));
    // runinfo=(TAGrunInfo*)(inputFile->Get("runinfo"));
    const TAGrunInfo construninfo(*runinfo);
    gTAGroot.SetRunInfo(construninfo);

    TString expName = runinfo->CampaignName();
    if (expName.EndsWith("/")) // fix a bug present in shoe
        expName.Remove(expName.Length() - 1);

    Int_t runNumber = runinfo->RunNumber();
    TAGrecoManager::Instance(expName);

    TAGrecoManager::GetPar()->FromFile();
    TAGrecoManager::GetPar()->Print();

    runinfo->Print();

    TAGcampaignManager *campManager = new TAGcampaignManager(expName);
    campManager->FromFile();

    // retrieve info from run and geometry
    GetRunAndGeoInfo(campManager, runNumber);

    const Char_t *name = FootActionDscName("TAGactTreeReader");
    TAGactTreeReader *fActReader = new TAGactTreeReader(name);

    SetTreeBranchAddress(fActReader);

    ////////////////////////////

    Int_t pos = outfile.Last('.');
    if (pos > 0)
        outfile = outfile(0, pos);

    Int_t inpos = infile.Last('/');
    TString in = "_" + infile(inpos + 1, infile.Length());
    cout << endl
         << "modified input name::  " << in.Data() << endl;

    outfile.Append(Form("%s", in.Data()));
    cout << outfile.Data() << endl;
    outfile = outDir + "/" + outfile;

    std::cout << "Input file (-in) : " << infile.Data() << std::endl;
    std::cout << "Output file (-out) : " << outfile.Data() << std::endl
              << endl;

    Int_t TG_region = parGeo->GetRegTarget();
    Int_t AirPreTW_region = parGeo->GetRegAirPreTW();
    Int_t AirTW_region = parGeo->GetRegAirTW();

    Int_t firstBarID(0), lastBarID(19);

    Int_t TW_regions[kTWreg] = {twparGeo->GetRegStrip((int)LayerY, firstBarID),
                                twparGeo->GetRegStrip((int)LayerY, lastBarID),
                                twparGeo->GetRegStrip((int)LayerX, firstBarID),
                                twparGeo->GetRegStrip((int)LayerX, lastBarID)}; // min and max rear and front TW regions

    printf("ExpName::%s  TG_region::%d  TW_regions: front=[%d,%d] rear=[%d,%d], AirPreTW::%d, AirTW::%d\n", expName.Data(), TG_region, TW_regions[2], TW_regions[3], TW_regions[0], TW_regions[1], AirPreTW_region, AirTW_region);
    // getchar();

    TFile *fout = new TFile(outfile.Data(), "RECREATE");
    fout->cd();

    BookHistograms();

    TVecPair vPairWrongZ;
    vPairWrongZ.clear();

    map<Int_t, vector<TVector3>> pMap; // primaries
    map<Int_t, vector<TVector3>> nMap; // neutrons
    pMap.clear();
    nMap.clear();

    cout << "Beam Z:: " << parGeo->GetBeamPar().AtomicNumber << " A:: " << parGeo->GetBeamPar().AtomicMass << " ion::" << parGeo->GetBeamPar().Material << endl;
    cout << "Beam Energy:: " << parGeo->GetBeamPar().Energy * 1000 << " MeV/u" << endl;
    cout << "TG center:: " << geoTrafo->GetTGCenter().z() << endl
         << " TG thickness:: " << parGeo->GetTargetPar().Size.z() << " TG material:: " << parGeo->GetTargetPar().Material << endl;

    printf("TW center =  (%2.1f,%2.1f,%2.1f)\n", GetTwCenter().x(), GetTwCenter().y(), GetTwCenter().z());

    cout << "theta angle of TW  acceptance::" << GetMaxAngle() * 180 / TMath::Pi() << endl;

    fActReader->Open(infile);
    gTAGroot.AddRequiredItem(fActReader);

    Int_t nentries = fActReader->NEvents();

    printf("Max available number of Entries in this Tree is::%d\n\n", (int)nentries);

    if (!isMax)
    {
        if (nev < nentries)
        {
            nentries = nev;
            if (nev <= 0)
                printf("number of events to be processed is %d...set a number > 0\n\n", nev);
        }
        else
            printf("Warning! nev (%d) has been set to a value > of maximum number of entries (%d), so take the maximum \n\n", nev, (int)nentries);
    }

    printf("Total Entries to be processed::%d\n\n", (int)nentries);

    string beamEnergyStr;

    if (runNumber == 4723 || runNumber == 4725 || runNumber == 4726 || runNumber == 4628)
    {
        beamEnergyStr = to_string(180.0);
    }
    else if (runNumber == 4727 || runNumber == 4728)
    {
        beamEnergyStr = to_string(140.0);
    }
    else if (runNumber == 4624)
    {
        beamEnergyStr = to_string(110.0);
    }
    else if (runNumber == 4625)
    {
        beamEnergyStr = to_string(130.0);
    }
    else
    {
        beamEnergyStr = to_string(parGeo->GetBeamPar().Energy * 1000);
    }

    TObjString objString(beamEnergyStr.c_str());
    TObjString materialObj(parGeo->GetBeamPar().Material);
    TObjString nentriesObj(to_string(nentries).c_str());

    fout->cd();
    // Check if the object already exists
    objString.Write(Form("BeamEnergyInfo run %d", runNumber));
    materialObj.Write(Form("IonInfo run %d", runNumber));
    nentriesObj.Write(Form("nentries run %d", runNumber));

    Int_t ev = -1;
    Int_t energy = std::stoi(beamEnergyStr);

    std::map<Int_t, Double_t> calibCoeff = extractCrystalData();
    std::map<std::pair<Int_t, Int_t>, TGraph *> scatterPlots;

    // Intercalibration coefficients to align the charge
    // values of crystal 1 and 6 with the one of crystal 0.
    // Valid only for 180 and 200 MeV/u beam energy values.
    // Obtained by fitting the upper diagonal of the crystal ID
    // 1 vs 0 and 6 vs 0 scatter plots.
    Double_t m_0_1, q_0_1;
    Double_t m_0_6, q_0_6;

    // intercalib refers to the ratio of the slopes of crystal IDs
    // 1 and 0 with respect to 0
    Double_t m_0_1_intercalib = -1.69269;
    Double_t m_0_6_intercalib = -1.63262;

    if (energy == 180)
    {
        m_0_1 = -1.382;
        q_0_1 = 0.3811;
        m_0_6 = -1.404;
        q_0_6 = 0.3786;
    }
    else
    {
        m_0_1 = -1.297;
        q_0_1 = 0.4074;
        m_0_6 = -1.481;
        q_0_6 = 0.4464;
    }

    // Loop over the TTree

    gTAGroot.BeginEventLoop();
    while (gTAGroot.NextEvent() && ev != nentries)
    {

        ev++;

        if (debug)
            printf("\n Event: %d\n", ev);
        else if (ev % 10000 == 0)
            printf("Event: %d\n", ev);

        Int_t trigID = -1;
        Int_t nMBplusVETOcounts(-1), nVETOcounts(-1), nMBcounts(-1), nSTcounts(-1), statusMB(-1);

        if (!IncludeMC)
        {
            trigID = wdNtuTrig->GetTriggerID();
            nMBplusVETOcounts = wdNtuTrig->GetTriggersCounter()[kMBplusVeto];
            nVETOcounts = wdNtuTrig->GetTriggersCounter()[kVeto];
            nMBcounts = wdNtuTrig->GetTriggersCounter()[kMB];
            nSTcounts = wdNtuTrig->GetTriggersCounter()[kSTtrig];
            statusMB = wdNtuTrig->GetTriggersStatus()[kMB];
        }

        Int_t nHitsCalo = caNtuHit->GetHitsN(); // number of hits in the calorimeter

        for (int ihit = 0; ihit < nHitsCalo; ihit++)
        {
            TACAhit *hit_ca = caNtuHit->GetHit(ihit);
            if (!hit_ca || !hit_ca->IsValid())
                continue;

            Int_t crystal_id = hit_ca->GetCrystalId();
            Double_t charge_calo = hit_ca->GetCharge();
            Int_t ModuleID = crystal_id / kCrysPerModule;
            TVector3 CaloPosition = hit_ca->GetPosition();

            Charge_Calo_crystal_noCuts[crystal_id]->Fill(charge_calo);

            if (charge_calo > 0.02)
            {
                //Charge_Calo_total->Fill(charge_calo);
                Charge_Calo_crystal[crystal_id]->Fill(charge_calo);
                // Charge_Calo_Module[ModuleID]->Fill(charge_calo);
            }
        }

        Int_t nClusters = caNtuClus->GetClustersN(); // number of clusters
        Clusters_number->Fill(nClusters);

        for (int icluster = 0; icluster < nClusters; icluster++)
        {
            TACAcluster *cluster = caNtuClus->GetCluster(icluster);
            if (!cluster || !cluster->IsValid())
                continue;

            Int_t nClusterHits = cluster->GetHitsN(); // i.e. cluster size
            Clusters_size_noCuts->Fill(nClusterHits);

            // Double_t charge_cluster_total = 0;
            // loop to find the total and minimum charge in the cluster
            Int_t valid_cluster_size = 0;
            Double_t min_charge_cluster = std::numeric_limits<double>::max();
            for (int iclusterhit = 0; iclusterhit < nClusterHits; iclusterhit++)
            {
                TACAhit *hit = cluster->GetHit(iclusterhit);
                if (!hit || !hit->IsValid())
                    continue;
                Double_t charge_cluster = hit->GetCharge();
                if (charge_cluster < min_charge_cluster)
                {
                    min_charge_cluster = charge_cluster;
                }
                if (charge_cluster > 0.02)
                {
                    valid_cluster_size++;
                }
                // charge_cluster_total += charge_cluster;
            }
            if(nClusterHits < 7) minCharge[nClusterHits]->Fill(min_charge_cluster);
            MinCharge_ClusterSize->Fill(nClusterHits, min_charge_cluster);
            Clusters_size->Fill(valid_cluster_size);

            // cluster size 1
            if (nClusterHits == 1 && nClusters == 1)
            {
                //cout << "event: " << ev << "cluster size = 1" << endl;
                for (int iclusterhit = 0; iclusterhit < nClusterHits; iclusterhit++)
                {
                    TACAhit *hit = cluster->GetHit(iclusterhit);
                    if (!hit || !hit->IsValid())
                        continue;

                    Int_t crystal_id = hit->GetCrystalId();
                    //cout << "crystal_id: " << crystal_id << endl;
                    Int_t ModuleID = crystal_id / kCrysPerModule;
                    Double_t charge_clusterHit = hit->GetCharge();
                    if (charge_clusterHit > 0.02)
                    {
                        
                        if (crystal_id == 1 || crystal_id == 6 || crystal_id == 0)
                        {
                            ClusterCharge_Calo_crystal[crystal_id]->Fill(charge_clusterHit);
                            Double_t charge_clusterHit_calibrated;
                            if (crystal_id == 1)
                            {
                                charge_clusterHit_calibrated = (1. / m_0_1) * (charge_clusterHit - q_0_1);
                                ClusterCharge_Calo_Calibrated[crystal_id]->Fill(charge_clusterHit_calibrated);
                            }
                            else if (crystal_id == 6)
                            {
                                charge_clusterHit_calibrated = (1. / m_0_6) * (charge_clusterHit - q_0_6);
                                ClusterCharge_Calo_Calibrated[crystal_id]->Fill(charge_clusterHit_calibrated);
                            }
                        }
                    }
                }
            }

            if (nClusterHits == 2 && nClusters == 1)
            {
                //cout << "event: " << ev << " cluster size = 2" << endl;
                std::vector<int> crystal_ids;
                std::vector<double> charge_values;

                Double_t charge_sum = 0;
                Double_t charge_sum_calibrated = 0;
                // the total charge of a cluster with size 2 is stored
                // only if one of the two crystals is crystal 0
                bool isZero = false;
                Double_t stored_charge = 0;
                Double_t nonCalib_stored_charge = 0;
                Double_t charge_0 = -100.;
                Double_t charge_1 = -100.;
                Double_t charge_6 = -100.;
                for (int iclusterhit = 0; iclusterhit < nClusterHits; iclusterhit++)
                {
                    TACAhit *hit = cluster->GetHit(iclusterhit);
                    if (!hit || !hit->IsValid())
                        continue;
                    // TVector3 ClusterHitPos = hit->GetPosition();
                    Int_t crystal_id = hit->GetCrystalId();
                    Int_t ModuleID = crystal_id / kCrysPerModule;
                    Double_t charge_clusterHit = hit->GetCharge();

                    if (crystal_id == 0 && charge_clusterHit > 0.02)
                    {
                        isZero = true;
                        charge_0 = charge_clusterHit;
                        charge_sum_calibrated += charge_clusterHit;
                        charge_sum += charge_clusterHit;
                        crystal_ids.push_back(0);
                        charge_values.push_back(charge_0);
                    }
                    // The find() method returns an iterator to the matching element if found, or end() if not found.
                    else if ((crystal_id == 1 || crystal_id == 6) && charge_clusterHit > 0.02)
                    {
                        //cout << "crystal_id part of the calibration file" << endl;
                        //cout << "normalizing its charge to crystal 0.." << endl;
                        nonCalib_stored_charge += charge_clusterHit;
                        if (crystal_id == 1)
                        {
                            charge_1 = charge_clusterHit;
                            crystal_ids.push_back(1);
                            charge_values.push_back(charge_1);
                            charge_clusterHit = (1. / m_0_1) * (charge_clusterHit - q_0_1);
                        }
                        else if (crystal_id == 6)
                        {
                            charge_6 = charge_clusterHit;
                            crystal_ids.push_back(6);
                            charge_values.push_back(charge_6);
                            charge_clusterHit = (1. / m_0_6) * (charge_clusterHit - q_0_6);
                        }
                        stored_charge += charge_clusterHit;
                    }
                }

                if (crystal_ids.size() == 2 && isZero && (charge_1 != -100. || charge_6 != -100.))
                {
                    Int_t id1 = crystal_ids[0], id2 = crystal_ids[1];
                    Double_t charge1 = charge_values[0], charge2 = charge_values[1];

                    //cout << "crystal_id1: " << id1 << " charge1: " << charge1 << endl;
                    //cout << "crystal_id2: " << id2 << " charge2: " << charge2 << endl;

                    // Ensure the pair is always stored in a consistent order
                    std::pair<Int_t, Int_t> crystalPair;
                    Double_t chargeLower, chargeHigher;

                    if (id1 < id2)
                    {
                        crystalPair = std::make_pair(id1, id2);
                        chargeLower = charge1;  // charge corresponding to the lower id
                        chargeHigher = charge2; // charge corresponding to the higher id
                    }
                    else
                    {
                        crystalPair = std::make_pair(id2, id1);
                        chargeLower = charge2; // swap the charges to follow the sorted order
                        chargeHigher = charge1;
                    }

                    if (crystalPair.first == 0 && (crystalPair.second == 1 || crystalPair.second == 6))
                    {
                        Correlated_ClusterCharge[crystalPair.first][crystalPair.second]->Fill(chargeLower, chargeHigher);

                        if (crystalPair.first != 0) cout << "first crystal: " << crystalPair.first << endl;

                        if (crystalPair.second == 1 && charge_6 != -100.)
                        {
                            cout << "2nd crystalID = 1, charge_1: " << charge_1 << " charge_6: " << charge_6 << endl;
                        }
                        else if (crystalPair.second == 6 && charge_1 != -100.)
                        {
                            cout << "2nd crystalID = 6, charge_1: " << charge_1 << " charge_6: " << charge_6 << endl;
                        }


                        // Create scatter plot if it doesn't exist
                        if (scatterPlots.find(crystalPair) == scatterPlots.end())
                        {
                            TString graphName = Form("scatter_Crystal%d_vs_Crystal%d", crystalPair.first, crystalPair.second);
                            scatterPlots[crystalPair] = new TGraph();
                            scatterPlots[crystalPair]->SetName(graphName);
                            scatterPlots[crystalPair]->SetTitle(Form("Charge Scatter Plot (cluster size = 2): Crystal %d vs Crystal %d", crystalPair.first, crystalPair.second));
                            scatterPlots[crystalPair]->SetMarkerStyle(20);
                            scatterPlots[crystalPair]->SetMarkerSize(0.35); // Default is 1
                            scatterPlots[crystalPair]->SetMarkerColor(kBlue);
                            scatterPlots[crystalPair]->SetMinimum(0.);
                            scatterPlots[crystalPair]->SetMaximum(0.5);
                            scatterPlots[crystalPair]->SetLineStyle(0);
                            scatterPlots[crystalPair]->SetLineWidth(0);
                        }

                        // Add point to scatter plot
                        scatterPlots[crystalPair]->SetPoint(scatterPlots[crystalPair]->GetN(), chargeLower, chargeHigher);

                        Double_t expected_charge_1 = m_0_1 * charge_0;
                        Double_t expected_charge_1_intercalib = m_0_1_intercalib * charge_0;
                        Double_t expected_charge_6 = m_0_6 * charge_0;
                        Double_t expected_charge_6_intercalib = m_0_6_intercalib * charge_0;

                        charge_sum_calibrated += stored_charge;
                        charge_sum += nonCalib_stored_charge;

                        Charge_Calo_nonCalibrated->Fill(charge_sum);
                        Charge_Calo_Calibrated->Fill(charge_sum_calibrated);

                        if (charge_1 != -100.)
                        {
                            Charge_Calo_Calibrated_q->Fill(charge_0);
                            Charge_Calo_Calibrated_q_intercalib->Fill(charge_0);
                            //Charge_Calo_Calibrated_q_0_1->Fill(charge_0);
                            Charge_Calo_Calibrated_q_intercalib_0_1->Fill(charge_0);
                            Charge_Calo_Calibrated_q->Fill(charge_1 - expected_charge_1);
                            Charge_Calo_Calibrated_q_0_1->Fill(charge_0 - (charge_1 / m_0_1));
                            Charge_Calo_Calibrated_q_intercalib->Fill(charge_1 - expected_charge_1_intercalib);
                            Charge_Calo_Calibrated_q_intercalib_0_1->Fill(charge_1 - expected_charge_1_intercalib);
                        }
                        if (charge_6 != -100.)
                        {
                            Charge_Calo_Calibrated_q->Fill(charge_0);
                            Charge_Calo_Calibrated_q_intercalib->Fill(charge_0);
                            Charge_Calo_Calibrated_q->Fill(charge_6 - expected_charge_6);
                            //Charge_Calo_Calibrated_q_0_6->Fill(charge_0);
                            Charge_Calo_Calibrated_q_intercalib_0_6->Fill(charge_0);
                            Charge_Calo_Calibrated_q_0_6->Fill(charge_0 - (charge_6 / m_0_6));
                            Charge_Calo_Calibrated_q_intercalib->Fill(charge_6 - expected_charge_6_intercalib);
                            Charge_Calo_Calibrated_q_intercalib_0_6->Fill(charge_6 - expected_charge_6_intercalib);
                        }
                    }
                }
            }
        }
    }
    // close for loop over events
    gTAGroot.EndEventLoop();

    SetTitleAndLabels(Charge_Calo_Calibrated, Form("Charge Calo Calibrated (single cluster of size 2) @ %d MeV/u", energy), "Charge [a.u.]", "");
    SetTitleAndLabels(Charge_Calo_nonCalibrated, Form("Non Calibrated Charge Calo (single cluster of size 2) @ %d MeV/u", energy), "Charge [a.u.]", "");

    cout << endl
         << "Job Done!" << endl;

    fout->cd();
    for (auto &entry : scatterPlots)
    {
        entry.second->Write();
    }
    fout->Write();
    fout->Close();

    return;
}

//-----------------------------------------------------------------------------
// Function to fit scatter plots and return fit results
std::pair<FitResult, TGraph*> FitScatterPlot(TGraph *graph, Int_t i, Int_t j, Double_t maxX, Double_t maxY, Double_t energy)
{
    FitResult result = {false, 0.0, 0.0, 0.0, 0.0};

    if (!graph)
        return {result, nullptr};

    std::vector<Double_t> xFiltered, yFiltered;

    // Threhsold to select the points between the two parallel lines
    // with slope -(maxY / maxX) with maximum y value: maxY +- line_threhsold 
    Double_t line_threshold = (energy == 200 || energy == 220) ? 0.02 : 0.01;

    for (int i = 0; i < graph->GetN(); ++i) {
        Double_t x, y;
        graph->GetPoint(i, x, y);

        // Apply selection condition
        if ((y > (maxY - line_threshold - x * (maxY / maxX)) && y < (maxY + 2*line_threshold - x * (maxY/ maxX))) && (y > 0.02)) {
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

void InitializeContainers()
{

    mapTrigHisto.clear();
    mapTrig.clear();

    mapTrig[kMBplusVeto] = "MBplusVeto";
    mapTrig[kVeto] = "Veto";
    mapTrig[kMB] = "MB";
    mapTrig[kSTtrig] = "STtrig";

    return;
}

std::map<Int_t, Double_t> extractCrystalData() {
    std::map<Int_t, Double_t> crystalData;
    std::string filename = "calib/HIT2022/SlopeRatios.cal"; // File is hardcoded
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return crystalData;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue; // Skip headers

        std::istringstream iss(line);
        Int_t crystalId;
        Double_t p0, p0_err;

        if (!(iss >> crystalId >> p0 >> p0_err)) continue; // Skip invalid lines

        crystalData[crystalId]= p0;
    }

    file.close();
    return crystalData;
}

void SetTitleAndLabels(TObject* obj, const char* title, const char* xLabel, const char* yLabel) {
    if (!obj) {
        std::cerr << "Error: Null object passed to SetTitleAndLabels." << std::endl;
        return;
    }

    if (TH1D* hist1D = dynamic_cast<TH1D*>(obj)) {
        hist1D->SetTitle(title);
        hist1D->GetXaxis()->SetTitle(xLabel);
        hist1D->GetYaxis()->SetTitle(yLabel);
    } 
    else if (TH2D* hist2D = dynamic_cast<TH2D*>(obj)) {
        hist2D->SetTitle(title);
        hist2D->GetXaxis()->SetTitle(xLabel);
        hist2D->GetYaxis()->SetTitle(yLabel);
    } 
    else if (TGraph* graph = dynamic_cast<TGraph*>(obj)) {
        graph->SetTitle(title);
        graph->GetXaxis()->SetTitle(xLabel);
        graph->GetYaxis()->SetTitle(yLabel);
        graph->SetLineStyle(0);
        graph->SetLineWidth(0);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.2);
    } 
    else {
        std::cerr << "Error: Object type not supported." << std::endl;
    }
}

//-----------------------------------------------------------------------------

void BookHistograms()
{

    // fpHisSeedMap = new TH1F(Form("msSeedMap%d", 4+1), Form("MSD - seed map for sensor %d", i+1), pGeoMap->GetStripsN(), 0, msdparGeo->GetStripsN());
    // AddHistogram(fpHisSeedMap);

    // fpHisStripMap = new TH1F(Form("msStripMap%d", 4+1), Form("MSD - strip map for sensor %d", i+1), pGeoMap->GetStripsN(), 0, msdparGeo->GetStripsN());
    // AddHistogram(fpHisStripMap);
    int modules[7] = {4, 5, 3, 6, 2, 7, 1}; // module enumeration following the excel file
    std::map<Int_t, Double_t> calibCoeff = extractCrystalData();
    //Charge_Calo_total = new TH1D(Form("Charge_Calo"), Form("Charge_Calo"), 500, -0.2, 1.5);
    Charge_Calo_Calibrated = new TH1D(Form("Charge_Calo_Calibrated"), Form("Charge Calo Calibrated (single cluster of size 2)"), 500, -0.2, 1.5);
    Charge_Calo_nonCalibrated = new TH1D(Form("Charge_Calo_nonCalibrated"), Form("Non-Calibrated Charge Calo (single cluster of size 2)"), 500, -0.2, 1.5);

    Charge_Calo_Calibrated_q = new TH1D(Form("Charge_Calo_Calibrated_q"), Form("Charge Calo Calibrated Intercept (single cluster of size 2)"), 500, -0.2, 1.5);
    Charge_Calo_Calibrated_q_intercalib = new TH1D(Form("Charge_Calo_Calibrated_q_intercalib"), Form("Charge Calo Calibrated Intercept Intercalib. (single cluster of size 2)"), 500, -0.2, 1.5);

    Charge_Calo_Calibrated_q_0_1 = new TH1D(Form("Charge_Calo_Calibrated_q_0_1"), Form("Charge Calo Calibrated Crystal 0 vs 1 Intercept (single cluster of size 2)"), 500, -0.2, 0.5);
    Charge_Calo_Calibrated_q_0_6 = new TH1D(Form("Charge_Calo_Calibrated_q_0_6"), Form("Charge Calo Calibrated Crystal 0 vs 6 Intercept (single cluster of size 2)"), 500, -0.2, 0.5);

    Charge_Calo_Calibrated_q_intercalib_0_1 = new TH1D(Form("Charge_Calo_Calibrated_q_intercalib_0_1"), Form("Charge Calo Calibrated Crystal 0 vs 1 Intercept Intercalib. (single cluster of size 2)"), 500, -0.2, 1.5);
    Charge_Calo_Calibrated_q_intercalib_0_6 = new TH1D(Form("Charge_Calo_Calibrated_q_intercalib_0_6"), Form("Charge Calo Calibrated Crystal 0 vs 6 Intercept Intercalib. (single cluster of size 2)"), 500, -0.2, 1.5);

    Clusters_number = new TH1D(Form("Clusters_number"), Form("Clusters_number"), 100, -1., 5.);
    Clusters_size_noCuts = new TH1D(Form("noCuts_Clusters_size"), Form("No Cuts Clusters_size"), 100, -1., 20.);
    Clusters_size = new TH1D(Form("Clusters_size"), Form("Clusters size"), 100, -1., 20.);

    for (int iclusterHit = 1; iclusterHit < 7; iclusterHit++)
    {
        minCharge[iclusterHit] = new TH1D(Form("MinCharge_ClusterSize_%d", iclusterHit), Form("Min. Charge Cluster Size %d", iclusterHit), 220, -0.5, 1.5);
    }
    MinCharge_ClusterSize = new TH2D(Form("h2D_MinCharge_ClusterSize"), Form("Min.Charge vs Cluster Size"), 100, 0, 6, 220, -0.5, 1.5);

    // Correlated charge for cluster-size = 2, for the central module (ids from 0 to 8)
    for (int icrystal_x = 0; icrystal_x < kCrysPerModule; icrystal_x++)
    {
        for (int icrystal_y = icrystal_x + 1; icrystal_y < kCrysPerModule; icrystal_y++)
        {
            Correlated_ClusterCharge[icrystal_x][icrystal_y] = new TH2D(Form("Correlated_ClusterCharge_ids_%d_%d", icrystal_x, icrystal_y), Form("Cluster Size 2: Charge crystalID %d vs crystalID %d", icrystal_x, icrystal_y), 200, 0., 0.55, 200, 0., 0.55);
            //Correlated_ClusterCharge[icrystal_x][icrystal_y] = new TH2D(Form("Correlated_ClusterCharge_ids_%d_%d", icrystal_x, icrystal_y), Form("Cluster Size 2: Charge crystalID %d vs crystalID %d", icrystal_x, icrystal_y), 100, 0., 0.55, 100, 0., 0.55);
        }
    }

    for (int icrystal = 0; icrystal < kModules * kCrysPerModule; icrystal++)
    {
        Charge_Calo_crystal[icrystal] = new TH1D(Form("Charge_Calo_crystalId_%d", icrystal), Form("Charge_Calo_crystalID_%d", icrystal), 500, -0.2, 1.);
        Charge_Calo_crystal_noCuts[icrystal] = new TH1D(Form("noCuts_Charge_Calo_crystalId_%d", icrystal), Form("No Cuts Charge Calo crystalID %d", icrystal), 500, -0.2, 1.);

        // hClusterSize_Charge[icrystal] = new TH2D(Form("ClusterSize_Charge_crystalId_%d", icrystal), Form("Cluster-Size vs Norm.Charge crystalId %d", icrystal), 100, -0.2, 1.1, 30, 0.5, 6.5);
        if (icrystal == 0 || icrystal == 1 || icrystal == 6)
        {
            ClusterCharge_Calo_crystal[icrystal] = new TH1D(Form("ClusterCharge_Calo_crystalId_%d", icrystal), Form("Single Cluster of Size 1: Non-Calibrated Charge crystalID %d", icrystal), 500, -0.2, 1.1);
            if (icrystal == 1 || icrystal == 6)
                ClusterCharge_Calo_Calibrated[icrystal] = new TH1D(Form("Calibrated_ClusterCharge_Calo_crystalId_%d", icrystal), Form("Single Cluster of Size 1: Calibrated Charge crystalID %d", icrystal), 500, -0.2, 1.1);
        }
    }

    return;
}

//-----------------------------------------------------------------------------

void ProjectTracksOnTw(int Z, TVector3 initPos, TVector3 initP)
{
    // Ogni traccia è una retta nello spazio di eq parametriche (espr. vettoriale) X = initPos + initP * t (parametro t reale) e si interseca con il piano z=z_TW (centro tra i due layer in z del TW) nei punti:

    Int_t x_intTW = initPos.x() + initP.x() / initP.z() * (GetTwCenter().z() - initPos.z());
    Int_t y_intTW = initPos.y() + initP.y() / initP.z() * (GetTwCenter().z() - initPos.z());

    // Select only TW portion of the plan z=z_TW:
    Int_t TwHalfLength = (nBarsPerLayer * twparGeo->GetBarWidth()) / 2;

    return;
}

//-----------------------------------------------------------------------------

Bool_t IsVTregion(int reg)
{

    TString firstVTregName = "VTXE0";
    TString lastVTregName = "VTXP3";

    // first and last VT region
    Int_t VT_regions[kVTreg] = {parGeo->GetCrossReg(firstVTregName),
                                parGeo->GetCrossReg(lastVTregName)};

    bool isvtreg = false;

    if (VT_regions[0] < reg < VT_regions[1])
        isvtreg = true;
    else
        isvtreg = false;

    return isvtreg;
}

//-----------------------------------------------------------------------------

void GetFOOTgeo(TAGcampaignManager *campManager,
                Int_t runNumber)
{

    // TAGgeoTrafo*
    geoTrafo = new TAGgeoTrafo();
    TString parFileName = campManager->GetCurGeoFile(TAGgeoTrafo::GetBaseName(), runNumber);
    geoTrafo->FromFile(parFileName);
    printf("geoTrafo::  %s\n", parFileName.Data());

    // beam
    TAGparaDsc *fpParGeoG = new TAGparaDsc(new TAGparGeo());

    // TAGparGeo*
    parGeo = (TAGparGeo *)fpParGeoG->Object();
    parFileName = campManager->GetCurGeoFile(TAGparGeo::GetBaseName(), runNumber);
    parGeo->FromFile(parFileName.Data());
    printf("parGeo::  %s\n", parFileName.Data());

    // Get the detectors geo files

    // ST
    TAGparaDsc *parGeoST = new TAGparaDsc(new TASTparGeo());

    // TASTparGeo*
    stparGeo = (TASTparGeo *)parGeoST->Object();
    parFileName = campManager->GetCurGeoFile(TASTparGeo::GetBaseName(), runNumber);
    stparGeo->FromFile(parFileName);
    printf("stparGeo::  %s\n", parFileName.Data());

    // BM
    TAGparaDsc *parGeoBm = new TAGparaDsc(new TABMparGeo());

    bmparGeo = (TABMparGeo *)parGeoBm->Object();
    parFileName = campManager->GetCurGeoFile(TABMparGeo::GetBaseName(), runNumber);
    bmparGeo->FromFile(parFileName);
    printf("bmparGeo::  %s\n", parFileName.Data());

    // VTX
    TAGparaDsc *parGeoVtx = new TAGparaDsc(new TAVTparGeo());

    // TAVTparGeo*
    vtparGeo = (TAVTparGeo *)parGeoVtx->Object();
    parFileName = campManager->GetCurGeoFile(TAVTparGeo::GetBaseName(), runNumber);
    vtparGeo->FromFile(parFileName);
    vtxSensorsN = vtparGeo->GetSensorsN();
    printf("vtparGeo::  %s\n", parFileName.Data());

    // IT
    TAGparaDsc *parGeoIt = new TAGparaDsc(new TAITparGeo());

    // TAITparGeo*
    itparGeo = (TAITparGeo *)parGeoIt->Object();
    parFileName = campManager->GetCurGeoFile(TAITparGeo::GetBaseName(), runNumber);
    itparGeo->FromFile(parFileName);
    itSensorsN = itparGeo->GetSensorsN();
    printf("itparGeo::  %s\n", parFileName.Data());

    // MSD
    TAGparaDsc *parGeoMsd = new TAGparaDsc(new TAMSDparGeo());

    // TAMSDparGeo*
    msdparGeo = (TAMSDparGeo *)parGeoMsd->Object();
    parFileName = campManager->GetCurGeoFile(TAMSDparGeo::GetBaseName(), runNumber);
    msdparGeo->FromFile(parFileName);
    msdSensorsN = msdparGeo->GetSensorsN();
    msdStationsN = msdparGeo->GetStationsN();
    printf("msdparGeo::  %s\n", parFileName.Data());

    // TW
    TAGparaDsc *parGeoTW = new TAGparaDsc(new TATWparGeo());

    // TATWparGeo*
    twparGeo = (TATWparGeo *)parGeoTW->Object();
    parFileName = campManager->GetCurGeoFile(TATWparGeo::GetBaseName(), runNumber);
    twparGeo->FromFile(parFileName);
    printf("twparGeo::  %s\n", parFileName.Data());

    // CA
    TAGparaDsc *parGeoCA = new TAGparaDsc(new TACAparGeo());

    // TACAparGeo*
    caparGeo = (TACAparGeo *)parGeoCA->Object();
    parFileName = campManager->GetCurGeoFile(TACAparGeo::GetBaseName(), runNumber);
    caparGeo->FromFile(parFileName);
    printf("caparGeo::  %s\n", parFileName.Data());

    return;
}

//-----------------------------------------------------------------------------

void GetRunAndGeoInfo(TAGcampaignManager *campManager, Int_t runNumber)
{

    IncludeTrk = runinfo->GetGlobalPar().EnableTracking;
    IncludeReg = runinfo->GetGlobalPar().EnableRegionMc;
    IncludeMC = campManager->GetCampaignPar(campManager->GetCurrentCamNumber()).McFlag;
    IncludeDI = campManager->IsDetectorOn("DI");
    IncludeSC = campManager->IsDetectorOn("ST");
    IncludeBM = campManager->IsDetectorOn("BM");
    IncludeTG = campManager->IsDetectorOn("TG");
    IncludeVT = campManager->IsDetectorOn("VT");
    IncludeIT = campManager->IsDetectorOn("IT");
    IncludeMSD = campManager->IsDetectorOn("MSD");
    IncludeTW = campManager->IsDetectorOn("TW");
    IncludeCA = campManager->IsDetectorOn("CA");

    if (!IncludeMC)
    {
        IncludeDAQ = true;
        IncludeWD = true;
    }

    if (IncludeMC && IncludeDAQ)
    {
        cout << "IncludeMC and IncludeDAQ are both true... check your input file and the configuration files, this program will be ended" << endl;
        return;
    }

    cout << endl;
    cout << "  Include DAQ:: " << IncludeDAQ << endl;
    cout << "  Include WD:: " << IncludeWD << endl;
    cout << "  Include MC:: " << IncludeMC << endl
         << endl;

    if (debug)
        campManager->Print();

    // global FOOT geometry
    GetFOOTgeo(campManager, runNumber);

    return;
}

//-----------------------------------------------------------------------------

void SetTreeBranchAddress(TAGactTreeReader *actTreeReader)
{

    // blocco traccia MC
    mcNtuPart = new TAMCntuPart();
    TAGdataDsc *mcPart = new TAGdataDsc(mcNtuPart);
    actTreeReader->SetupBranch(mcPart);

    TString name(mcNtuPart->ClassName());
    const char *branch = TAGnameManager::GetBranchName(name);
    // cout<<branch<<endl;
    printf("%s\n", TAGnameManager::GetBranchName(name).Data());

    // blocco crossings MC
    mcNtuRegion = new TAMCntuRegion();
    TAGdataDsc *mcRegion = new TAGdataDsc(mcNtuRegion);
    actTreeReader->SetupBranch(mcRegion);
    cout << TAGnameManager::GetBranchName(mcNtuRegion->ClassName()) << endl;

    // ST
    stNtuHit = new TASTntuHit(); // blocco hit ST reco
    TAGdataDsc *stHit = new TAGdataDsc(stNtuHit);
    actTreeReader->SetupBranch(stHit);
    cout << TAGnameManager::GetBranchName(stNtuHit->ClassName()) << endl;

    if (IncludeMC)
    {
        stMcNtuHit = new TAMCntuHit();
        TAGdataDsc *stMcHit = new TAGdataDsc(FootDataDscMcName(kST), stMcNtuHit);
        actTreeReader->SetupBranch(stMcHit, FootBranchMcName(kST));
        // actTreeReader->SetupBranch(stMcHit);
        cout << TAGnameManager::GetBranchMcName(kST) << endl;
    }

    // TW
    twNtuHit = new TATWntuHit(); // blocco hit TW reco
    TAGdataDsc *twHit = new TAGdataDsc(twNtuHit);
    actTreeReader->SetupBranch(twHit);
    cout << TAGnameManager::GetBranchName(twNtuHit->ClassName()) << endl;

    twNtuPoint = new TATWntuPoint();
    TAGdataDsc *twPoint = new TAGdataDsc(twNtuPoint);
    actTreeReader->SetupBranch(twPoint);
    cout << TAGnameManager::GetBranchName(twNtuPoint->ClassName()) << endl;

    if (IncludeMC)
    {
        twMcNtuHit = new TAMCntuHit(); // blocco hit TW MC
        // TAGdataDsc* twMcHit    = new TAGdataDsc(FootBranchMcName(kTW),twMcNtuHit);
        // actTreeReader->SetupBranch(twMcHit);
        TAGdataDsc *twMcHit = new TAGdataDsc(FootDataDscMcName(kTW), twMcNtuHit);
        actTreeReader->SetupBranch(twMcHit, FootBranchMcName(kTW));
        cout << TAGnameManager::GetBranchMcName(kTW) << endl;
    }

    // BM
    bmNtuHit = new TABMntuHit();
    TAGdataDsc *bmHit = new TAGdataDsc(bmNtuHit);
    actTreeReader->SetupBranch(bmHit);
    cout << TAGnameManager::GetBranchName(bmNtuHit->ClassName()) << endl;

    if (IncludeTrk)
    {
        bmNtuTrack = new TABMntuTrack();
        TAGdataDsc *bmTrack = new TAGdataDsc(bmNtuTrack);
        actTreeReader->SetupBranch(bmTrack);
        cout << TAGnameManager::GetBranchName(bmNtuTrack->ClassName()) << endl;
    }
    if (IncludeMC)
    {
        bmMcNtuHit = new TAMCntuHit();
        TAGdataDsc *bmMcHit = new TAGdataDsc(FootDataDscMcName(kBM), bmMcNtuHit);
        // TAGdataDsc* bmMcHit    = new TAGdataDsc(FootBranchMcName(kBM),bmMcNtuHit);
        actTreeReader->SetupBranch(bmMcHit, FootBranchMcName(kBM));
        // TAGdataDsc* bmMcHit    = new TAGdataDsc(FootBranchMcName(kBM),bmMcNtuHit);
        // actTreeReader->SetupBranch(bmMcHit);
        cout << TAGnameManager::GetBranchMcName(kBM) << endl;
    }

    // VT
    if (IncludeTrk)
    {

        vtxNtuVertex = new TAVTntuVertex();
        TAGdataDsc *vtVertex = new TAGdataDsc(vtxNtuVertex);
        actTreeReader->SetupBranch(vtVertex);

        vtxNtuTrack = new TAVTntuTrack();
        TAGdataDsc *vtTrack = new TAGdataDsc(vtxNtuTrack);
        actTreeReader->SetupBranch(vtTrack);
    }

    // Int_t sensorsN = vtparGeo->GetSensorsN();
    vtxNtuCluster = new TAVTntuCluster(vtxSensorsN);
    // vtxNtuCluster = new TAVTntuCluster(sensorsN);
    TAGdataDsc *vtCluster = new TAGdataDsc(vtxNtuCluster);
    actTreeReader->SetupBranch(vtCluster);

    if (IncludeMC)
    {
        vtMcNtuHit = new TAMCntuHit();
        // TAGdataDsc* vtMcHit    = new TAGdataDsc(FootBranchMcName(kVTX),vtMcNtuHit);
        // actTreeReader->SetupBranch(vtMcHit);
        TAGdataDsc *vtMcHit = new TAGdataDsc(FootDataDscMcName(kVTX), vtMcNtuHit);
        actTreeReader->SetupBranch(vtMcHit, FootBranchMcName(kVTX));
    }

    // IT
    itNtuClus = new TAITntuCluster(itSensorsN);
    TAGdataDsc *itClus = new TAGdataDsc(itNtuClus);
    actTreeReader->SetupBranch(itClus);

    if (IncludeMC)
    {
        itMcNtuHit = new TAMCntuHit();
        TAGdataDsc *itMcHit = new TAGdataDsc(FootDataDscMcName(kITR), itMcNtuHit);
        actTreeReader->SetupBranch(itMcHit, FootBranchMcName(kITR));
    }

    // MSD
    msdNtuClus = new TAMSDntuCluster(msdSensorsN);
    TAGdataDsc *msdClus = new TAGdataDsc(msdNtuClus);
    actTreeReader->SetupBranch(msdClus);

    msdNtuPoint = new TAMSDntuPoint(msdStationsN);
    TAGdataDsc *msdPoint = new TAGdataDsc(msdNtuPoint);
    actTreeReader->SetupBranch(msdPoint);

    msdNtuHit = new TAMSDntuHit(msdSensorsN);
    TAGdataDsc *msdHit = new TAGdataDsc(msdNtuHit);
    actTreeReader->SetupBranch(msdHit);

    msdNtuRaw = new TAMSDntuRaw(msdSensorsN);
    TAGdataDsc *msdRaw = new TAGdataDsc(msdNtuRaw);
    actTreeReader->SetupBranch(msdRaw);

    if (IncludeMC)
    {
        msdMcNtuHit = new TAMCntuHit();
        TAGdataDsc *msdMcHit = new TAGdataDsc(FootDataDscMcName(kMSD), msdMcNtuHit);
        actTreeReader->SetupBranch(msdMcHit, FootBranchMcName(kMSD));
    }

    // CA
    caNtuHit = new TACAntuHit();
    TAGdataDsc *caHit = new TAGdataDsc(caNtuHit);
    actTreeReader->SetupBranch(caHit);

    caNtuClus = new TACAntuCluster();
    TAGdataDsc *caClus = new TAGdataDsc(caNtuClus);
    actTreeReader->SetupBranch(caClus);

    if (IncludeMC)
    {
        caMcNtuHit = new TAMCntuHit();
        TAGdataDsc *caMcHit = new TAGdataDsc(FootDataDscMcName(kCAL), caMcNtuHit);
        actTreeReader->SetupBranch(caMcHit, FootBranchMcName(kCAL));
        // actTreeReader->SetupBranch(caMcHit);
    }

    // DAQ
    if (IncludeDAQ && IncludeWD)
    {
        tgNtuEvent = new TAGntuEvent();
        TAGdataDsc *tgEvent = new TAGdataDsc(tgNtuEvent);
        actTreeReader->SetupBranch(tgEvent);

        wdNtuTrig = new TAWDntuTrigger();
        TAGdataDsc *wdTrg = new TAGdataDsc(wdNtuTrig);
        actTreeReader->SetupBranch(wdTrg);
        // tree->SetBranchAddress(TAWDntuTrigger::GetBranchName(), &wdNtuTrig);
    }

    return;
}
