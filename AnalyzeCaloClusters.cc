// written by mtoppi 05/2023
// to be run in shoe/build/Reconstruction
// [ therein set-up the enviroment:  source ../setupFOOT.sh ]
// and run with (for example): root -l -b -q AnalyzeFOOT.cc++g\(\"../../../../rootfiles/outMC_16O_C_400_1_GSI.root\",1,10,\"\"\)
// sul tier1:  root -l -b -q AnalyzeFOOT.cc++g\(\"/storage/gpfs_data/foot/mtoppi/DataDecoded/CNAO2023/test.root\",0,1000,\"testAnaFOOT\",\"/storage/gpfs_data/foot/mtoppi/OutputMacro/\"\)

#include "AnalyzeCaloClusters.h"

// main
void AnalyzeCaloClusters(TString infile = "testMC.root", Bool_t isMax = kFALSE, Int_t nev = 10, TString outfile = "AnaFOOT.root", TString outDir = "/Users/marco/FOOT/Analisi/shoe/build/Reconstruction/OutputMacro")

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

    // Loop over the TTree
    gTAGroot.BeginEventLoop();
    Int_t ev = -1;
    std::map<std::pair<int, int>, TGraph *> scatterPlots;
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

            if (charge_calo > 0.02)
            {
                Charge_Calo_total->Fill(charge_calo);
                Charge_Calo_crystal[crystal_id]->Fill(charge_calo);
                // Charge_Calo_Module[ModuleID]->Fill(charge_calo);
                hCalMapPos[ModuleID]->Fill(CaloPosition.X(), CaloPosition.Y());
                double valueToSet = (crystal_id == 0) ? 0.0001 : static_cast<double>(crystal_id);
                hCalMapCrystalID[ModuleID]->SetBinContent(
                    hCalMapCrystalID[ModuleID]->GetXaxis()->FindBin(CaloPosition.X()),
                    hCalMapCrystalID[ModuleID]->GetYaxis()->FindBin(CaloPosition.Y()),
                    valueToSet);
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
            Clusters_size->Fill(nClusterHits);

            // Double_t charge_cluster_total = 0;
            //  loop to find the total charge in the cluster
            Int_t valid_cluster_size = 0;
            for (int iclusterhit = 0; iclusterhit < nClusterHits; iclusterhit++)
            {
                TACAhit *hit = cluster->GetHit(iclusterhit);
                if (!hit || !hit->IsValid())
                    continue;
                Double_t charge_cluster = hit->GetCharge();
                if (charge_cluster > 0.02) {
                    valid_cluster_size++;
                }
                // charge_cluster_total += charge_cluster;
            }
            Clusters_size_noCuts->Fill(valid_cluster_size);

            // cluster size 1
            if (nClusterHits == 1)
            {
                cout << "event: " << ev << "cluster size = 1" << endl;
                for (int iclusterhit = 0; iclusterhit < nClusterHits; iclusterhit++)
                {
                    TACAhit *hit = cluster->GetHit(iclusterhit);
                    if (!hit || !hit->IsValid())
                        continue;

                    Int_t crystal_id = hit->GetCrystalId();
                    cout << "crystal_id: " << crystal_id << endl;
                    Int_t ModuleID = crystal_id / kCrysPerModule;
                    TVector3 CaloPosition = hit->GetPosition();
                    Double_t charge_clusterHit = hit->GetCharge();

                    // if (crystal_id == 0) scatterPlot->SetPoint(pointIndex++, charge_cluster, nClusterHits);
                    if (charge_clusterHit > 0.02)
                    {
                        ClusterCharge_Calo_crystal[crystal_id]->Fill(charge_clusterHit);
                        // hClusterSize_Charge[crystal_id]->Fill(charge_clusterHit, nClusterHits);
                    }
                }
            }
            if (nClusterHits == 2 && nClusters == 1)
            {
                cout << "event: " << ev << " cluster size = 2" << endl;
                std::vector<int> crystal_ids;
                std::vector<double> charge_values;
                for (int iclusterhit = 0; iclusterhit < nClusterHits; iclusterhit++)
                {
                    TACAhit *hit = cluster->GetHit(iclusterhit);
                    if (!hit || !hit->IsValid())
                        continue;
                    // TVector3 ClusterHitPos = hit->GetPosition();
                    Int_t crystal_id = hit->GetCrystalId();
                    Int_t ModuleID = crystal_id / kCrysPerModule;
                    Double_t charge_clusterHit = hit->GetCharge();
                    
                    if (charge_clusterHit > 0.02)
                    {
                        crystal_ids.push_back(crystal_id);
                        charge_values.push_back(charge_clusterHit);
                    }
                }

                if (crystal_ids.size() == 2)
                {
                    int id1 = crystal_ids[0], id2 = crystal_ids[1];
                    double charge1 = charge_values[0], charge2 = charge_values[1];

                    cout << "crystal_id1: " << id1 << " charge1: " << charge1 << endl;
                    cout << "crystal_id2: " << id2 << " charge2: " << charge2 << endl;

                    // Ensure the pair is always stored in a consistent order
                    std::pair<int, int> crystalPair;
                    double chargeLower, chargeHigher;

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

                    if (crystalPair.first < 9 && crystalPair.second < 9)
                    {
                        Correlated_ClusterCharge[crystalPair.first][crystalPair.second]->Fill(chargeLower, chargeHigher);

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
                            scatterPlots[crystalPair]->GetXaxis()->SetLimits(0., 0.55);
                            scatterPlots[crystalPair]->GetXaxis()->SetRangeUser(0., 0.55);
                            scatterPlots[crystalPair]->SetMinimum(0.);
                            scatterPlots[crystalPair]->SetMaximum(0.55);
                            scatterPlots[crystalPair]->SetLineStyle(0);
                            scatterPlots[crystalPair]->SetLineWidth(0);
                        }

                        // Add point to scatter plot
                        scatterPlots[crystalPair]->SetPoint(scatterPlots[crystalPair]->GetN(), chargeLower, chargeHigher);
                    }
                }
            }
        }
    }
    // close for loop over events
    gTAGroot.EndEventLoop();

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

//-----------------------------------------------------------------------------

void BookHistograms()
{

    // fpHisSeedMap = new TH1F(Form("msSeedMap%d", 4+1), Form("MSD - seed map for sensor %d", i+1), pGeoMap->GetStripsN(), 0, msdparGeo->GetStripsN());
    // AddHistogram(fpHisSeedMap);

    // fpHisStripMap = new TH1F(Form("msStripMap%d", 4+1), Form("MSD - strip map for sensor %d", i+1), pGeoMap->GetStripsN(), 0, msdparGeo->GetStripsN());
    // AddHistogram(fpHisStripMap);
    int modules[7] = {4, 5, 3, 6, 2, 7, 1}; // module enumeration following the excel file

    Charge_Calo_total = new TH1D(Form("Charge_Calo"), Form("Charge_Calo"), 500, -0.2, 1.5);

    Clusters_number = new TH1D(Form("Clusters_number"), Form("Clusters_number"), 100, -1., 5.);
    Clusters_size = new TH1D(Form("Clusters_size"), Form("Clusters_size"), 100, -1., 20.);
    Clusters_size_noCuts = new TH1D(Form("Clusters_size_noCuts"), Form("Clusters size no Cuts"), 100, -1., 20.);

    // Correlated charge for cluster-size = 2, for the central module (ids from 0 to 8)
    for (int icrystal_x = 0; icrystal_x < kCrysPerModule; icrystal_x++)
    {
        for (int icrystal_y = icrystal_x + 1; icrystal_y < kCrysPerModule; icrystal_y++)
        {
            Correlated_ClusterCharge[icrystal_x][icrystal_y] = new TH2D(Form("Correlated_ClusterCharge_ids_%d_%d", icrystal_x, icrystal_y), Form("Cluster Size 2: Charge crystalID %d vs crystalID %d", icrystal_x, icrystal_y), 200, 0., 0.55, 200, 0., 0.55);
        }
    }

    for (int icrystal = 0; icrystal < kModules * kCrysPerModule; icrystal++)
    {
        Charge_Calo_crystal[icrystal] = new TH1D(Form("Charge_Calo_crystalId_%d", icrystal), Form("Charge_Calo_crystalID_%d", icrystal), 500, -0.2, 1.);

        // hClusterSize_Charge[icrystal] = new TH2D(Form("ClusterSize_Charge_crystalId_%d", icrystal), Form("Cluster-Size vs Norm.Charge crystalId %d", icrystal), 100, -0.2, 1.1, 30, 0.5, 6.5);
        ClusterCharge_Calo_crystal[icrystal] = new TH1D(Form("ClusterCharge_Calo_crystalId_%d", icrystal), Form("Cluster Size 1: Charge crystalID %d", icrystal), 500, -0.2, 1.1);
    }
    for (int imodule = 0; imodule < kModules; imodule++)
    {
        hCalMapPos[imodule] = new TH2D(Form("hCalMapPos_module_%d", modules[imodule]), Form("hCalMapPos_module_%d", modules[imodule]), 27, -27., 27., 11, -11., 11.);
        hCalMapCrystalID[imodule] = new TH2D(Form("hCalMapCrystalID_module_%d", modules[imodule]), Form("hCalMapCrystalID_module_%d", modules[imodule]), 27, -27., 27., 11, -11., 11.);
        // Charge_Calo_Module[imodule] = new TH1D(Form("Charge_Calo_Module_%d", modules[imodule]), Form("Charge_Calo_Module_%d", modules[imodule]), 200, 0., 1.);
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
