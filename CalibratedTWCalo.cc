// written by mtoppi 05/2023
// to be run in shoe/build/Reconstruction
// [ therein set-up the enviroment:  source ../setupFOOT.sh ]
// and run with (for example): root -l -b -q AnalyzeFOOT.cc++g\(\"../../../../rootfiles/outMC_16O_C_400_1_GSI.root\",1,10,\"\"\)
// sul tier1:  root -l -b -q AnalyzeFOOT.cc++g\(\"/storage/gpfs_data/foot/mtoppi/DataDecoded/CNAO2023/test.root\",0,1000,\"testAnaFOOT\",\"/storage/gpfs_data/foot/mtoppi/OutputMacro/\"\)

#include "CalibratedTWCalo.h"

// main
void CalibratedTWCalo(TString infile = "testMC.root", Bool_t isMax = kFALSE, Int_t nev = 10, TString outfile = "AnaFOOT.root", TString outDir = "/Users/marco/FOOT/Analisi/shoe/build/Reconstruction/OutputMacro")

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
  TDirectory *DirChargeElossLayerX = fout->mkdir("ChargeElossLayerX");
  TDirectory *DirChargeElossLayerY = fout->mkdir("ChargeElossLayerY");
  TDirectory *DirToFLayerX = fout->mkdir("ToFLayerX");
  TDirectory *DirToFLayerY = fout->mkdir("ToFLayerY");
  BookHistograms(DirChargeElossLayerX, DirChargeElossLayerY, DirToFLayerX, DirToFLayerY);

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
  objString.Write(Form("BeamEnergyInfo run %d", runNumber));
  materialObj.Write(Form("IonInfo run %d", runNumber));
  nentriesObj.Write(Form("nentries run %d", runNumber));

  // Int_t pointIndex = 0;
  // TCanvas *c = new TCanvas("c", "Scatter Plot", 800, 600);
  // scatter plot of beta (x) vs dE/r*dx (y) [MeV * cm^2 / g]
  // TGraph *betaEloss = new TGraph();

  TCanvas *c = new TCanvas("c", "Scatter Plot", 800, 600);
  TGraph *scatterPlot_bar9_layerX = new TGraph();
  TGraph *scatterPlot_bar9_layerY = new TGraph();
  TGraph *scatterPlot_bar9_layerX_filtered = new TGraph();
  TGraph *scatterPlot_bar9_layerY_filtered = new TGraph();
  Int_t pointIndex_bar9_layerX = 0;
  Int_t pointIndex_bar9_layerY = 0;
  Int_t pointIndex_bar9_layerX_filtered = 0;
  Int_t pointIndex_bar9_layerY_filtered = 0;

  Double_t d_SC_TW = 1.;  // distance between SC and TW, in m
  Double_t bar_density = 1.023;  // density of the bars in g/cm^3
  Double_t bar_thickness = 0.3;  // thickness of the bars in cm

  Int_t energy = std::stoi(beamEnergyStr);
  Int_t ev = -1;

  // Loop over the TTree to build the ampl-charge scatterplot
  // to be fitted with a linear function, to discard the pileup hits
  // which are far from the charge-ampl fit line
  gTAGroot.BeginEventLoop();

  while (gTAGroot.NextEvent() && ev != nentries)
  {
    ev++;
    Int_t nHitsX = twNtuHit->GetHitN((Int_t)LayerX);
    Int_t nHitsY = twNtuHit->GetHitN((Int_t)LayerY);

    if (debug)
      cout << " TWhits X::" << nHitsX << " Y::" << nHitsY << endl;

    for (int ihitX = 0; ihitX < nHitsX; ihitX++)
    {
      TATWhit *hitX = twNtuHit->GetHit(ihitX, (Int_t)LayerX);
      Double_t barX = hitX->GetBar();
      Double_t QAX = hitX->GetChargeChA();
      Double_t QBX = hitX->GetChargeChB();
      Double_t QBarX = sqrt(QAX * QBX);
      Double_t amplAX = hitX->GetAmplitudeChA();
      Double_t amplBX = hitX->GetAmplitudeChB();
      Double_t amplX = sqrt(amplAX * amplBX);
      if (barX == 9 && hitX->IsValid())
      {
        if (std::isnan(amplX) || std::isnan(QBarX) || std::isinf(amplX) || std::isinf(QBarX)) continue;
        scatterPlot_bar9_layerX->SetPoint(pointIndex_bar9_layerX++, amplX, QBarX);
      }
    }
    for (int ihitY = 0; ihitY < nHitsY; ihitY++)
    {
      TATWhit *hitY = twNtuHit->GetHit(ihitY, (Int_t)LayerY);
      Double_t barY = hitY->GetBar();
      Double_t QAY = hitY->GetChargeChA();
      Double_t QBY = hitY->GetChargeChB();
      Double_t QBarY = sqrt(QAY * QBY);
      Double_t amplAY = hitY->GetAmplitudeChA();
      Double_t amplBY = hitY->GetAmplitudeChB();
      Double_t amplY = sqrt(amplAY * amplBY);
      if (barY == 9 && hitY->IsValid())
      {
        if (std::isnan(amplY) || std::isnan(QBarY) || std::isinf(amplY) || std::isinf(QBarY)) continue;
        scatterPlot_bar9_layerY->SetPoint(pointIndex_bar9_layerY++, amplY, QBarY);
      }
    }
  }

  gTAGroot.EndEventLoop();

  cout << "Ended first event loop" << endl;

  int nPointsX = scatterPlot_bar9_layerX->GetN();
  int nPointsY = scatterPlot_bar9_layerY->GetN();

  if (nPointsX == 0 || nPointsY == 0) {
      std::cerr << "Error: No points in scatter plots!" << std::endl;
  }

  TF1 *fitFunc_layerX = new TF1("fitFunc_layerX", "[0]*x", 0., 0.99);
  TF1 *fitFunc_layerY = new TF1("fitFunc_layerY", "[0]*x", 0., 0.99);

  fitFunc_layerX->SetParameter(0, 35.);
  fitFunc_layerY->SetParameter(0, 35.);

  TFitResultPtr fitX = scatterPlot_bar9_layerX->Fit(fitFunc_layerX, "SQROBR");
  TFitResultPtr fitY = scatterPlot_bar9_layerY->Fit(fitFunc_layerY, "SQROBR");

  Double_t slopeX = fitFunc_layerX->GetParameter(0);
  Double_t slopeY = fitFunc_layerY->GetParameter(0);

  cout << "slopeX: " << slopeX << " slopeY: " << slopeY << endl;

  if (fitX->IsValid() && !fitY->IsValid())
  {
    cout << "fit layerY failed, setting slopeY = slopeX" << endl;
    slopeY = slopeX;
  }
  else if (fitY->IsValid() && !fitX->IsValid())
  {
    cout << "fit layerX failed, setting slopeX = slopeY" << endl;
    slopeX = slopeY;
  }
  else if (!fitX->IsValid() && !fitY->IsValid())
  {
    cout << "both fit failed, setting slopes to 33." << endl;
    slopeX = 33.;
    slopeY = 33.;
  }

  // needed for the second loop
  ev = -1;
  fActReader->Open(infile);
  gTAGroot.AddRequiredItem(fActReader);

  std::map<Int_t, Double_t> calibCoeff = extractCrystalData();

  gTAGroot.BeginEventLoop();

  // Second loop to build the charge-eloss histograms selecting events
  // with 1 valid hit (defined with a threshold on charge values) on both
  // layers and discarding pileup events (charge > charge_threshold from the fit line)
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

    if (debug)
      cout << "trigID::" << trigID << " nMBplusVeto::" << nMBplusVETOcounts << " nVETO::" << nVETOcounts << " nMB::" << nMBcounts << " nSTcounts::" << nSTcounts << "  statusMB::" << statusMB << endl;

    Int_t nHitsX = twNtuHit->GetHitN((Int_t)LayerX);
    Int_t nHitsY = twNtuHit->GetHitN((Int_t)LayerY);

    hHits_X->Fill(nHitsX);
    hHits_Y->Fill(nHitsY);

    Int_t nValidHitsX = 0;
    Int_t nValidHitsY = 0;

    // Track which hitY has already been counted
    std::vector<bool> countedY(nHitsY, false);

    if (debug)
      cout << " TWhits X::" << nHitsX << " Y::" << nHitsY << endl;

    // Indexes to select the valid hits and fill the histograms
    // when there's 1 valid hit on both X and Y layers
    Int_t hitNumber_X;
    Int_t hitNumber_Y;

    Double_t charge_threshold = (energy == 220) ? 0.5 : 0.7;

    for (int ihitX = 0; ihitX < nHitsX; ihitX++)
    {

      TATWhit *hitX = twNtuHit->GetHit(ihitX, (Int_t)LayerX);
      Double_t posBarY = twparGeo->GetBarPosition((Int_t)LayerX, hitX->GetBar())[1];
      Double_t barX = hitX->GetBar();
      Double_t amplAX = hitX->GetAmplitudeChA();
      Double_t amplBX = hitX->GetAmplitudeChB();
      Double_t amplX = sqrt(amplAX * amplBX);
      Double_t QAX = hitX->GetChargeChA();
      Double_t QBX = hitX->GetChargeChB();
      Double_t QBarX = sqrt(QAX * QBX);
      Double_t expected_chargeX = slopeX * amplX;

      if (hitX->IsValid() && (barX != 9 || (barX == 9 && (fabs(QBarX - expected_chargeX) < charge_threshold))) && ((energy == 100 && QBarX > 1.) ||
          (energy == 140 && QBarX > 0.95) ||
          (energy == 200 && QBarX > 0.8) ||
          (energy == 220 && QBarX > 0.7)))
      {
        nValidHitsX++;
        hitNumber_X = ihitX;
        if (barX == 9)
        {
          scatterPlot_bar9_layerX_filtered->SetPoint(pointIndex_bar9_layerX_filtered++, amplX, QBarX);
        }
      }

      for (int ihitY = 0; ihitY < nHitsY; ihitY++)
      {

        TATWhit *hitY = twNtuHit->GetHit(ihitY, (Int_t)LayerY);
        Double_t posBarX = twparGeo->GetBarPosition((Int_t)LayerY, hitY->GetBar())[0];
        Double_t barY = hitY->GetBar();
        Double_t amplAY = hitY->GetAmplitudeChA();
        Double_t amplBY = hitY->GetAmplitudeChB();
        Double_t amplY = sqrt(amplAY * amplBY);
        Double_t QAY = hitY->GetChargeChA();
        Double_t QBY = hitY->GetChargeChB();
        Double_t QBarY = sqrt(QAY * QBY);
        Double_t expected_chargeY = slopeY * amplY;

        if (hitY->IsValid() && !countedY[ihitY] && (barY != 9 || (barY == 9 && (fabs(QBarY - expected_chargeY) < charge_threshold))) && ((energy == 100 && QBarY > 1.) ||
          (energy == 140 && QBarY > 0.95) ||
          (energy == 200 && QBarY > 0.8) ||
          (energy == 220 && QBarY > 0.7)))
        {
          nValidHitsY++;
          countedY[ihitY] = true; // Mark this hitY as counted
          hitNumber_Y = ihitY;
          if (barY == 9)
          {
            scatterPlot_bar9_layerY_filtered->SetPoint(pointIndex_bar9_layerY_filtered++, amplY, QBarY);
          }
        }
      }
    }

    if (ev % 10000 == 0)
    {
      cout << "Energy [MeV/u]: " << energy << endl;
      cout << "nHitsX: " << nHitsX << " nHitsY: " << nHitsY << endl;
      cout << "nValidHitsX: " << nValidHitsX << " nValidHitsY: " << nValidHitsY << endl;
    }

    h_nValidHits_X->Fill(nValidHitsX);
    h_nValidHits_Y->Fill(nValidHitsY);

    if (nValidHitsX == 1 && nValidHitsY == 1)
    {
      TATWhit *hitX = twNtuHit->GetHit(hitNumber_X, (Int_t)LayerX);
      Int_t barX = hitX->GetBar();
      Double_t posAlongX = hitX->GetPosition();
      Double_t QAX = hitX->GetChargeChA();
      Double_t QBX = hitX->GetChargeChB();
      Double_t QBarX = sqrt(QAX * QBX);
      Double_t tofX = hitX->GetToF();
      Double_t elossX = hitX->GetEnergyLoss();
      Double_t ZX = hitX->GetChargeZ();

      TATWhit *hitY = twNtuHit->GetHit(hitNumber_Y, (Int_t)LayerY);
      Int_t barY = hitY->GetBar();
      Double_t posAlongY = hitY->GetPosition();
      Double_t QAY = hitY->GetChargeChA();
      Double_t QBY = hitY->GetChargeChB();
      Double_t QBarY = sqrt(QAY * QBY);
      Double_t tofY = hitY->GetToF();
      Double_t elossY = hitY->GetEnergyLoss();
      Double_t ZY = hitY->GetChargeZ();

      Bar_ID_X->Fill(barX);
      Bar_ID_Y->Fill(barY);

      Int_t nClusters = caNtuClus->GetClustersN(); // number of clusters

      if (ZX == ZY)
      {
        // selecting a single cluster
        if (nClusters == 1)
        {
          TACAcluster *cluster = caNtuClus->GetCluster(0);
          if (cluster && cluster->IsValid())
          {
              Int_t nClusterHits = cluster->GetHitsN(); // i.e. cluster size
              if (nClusterHits == 1)
              {
                Z_clusterSize1->Fill(ZX);
                Z_clusterSize1->Fill(ZY);

                dE_vs_tof_clusterSize1->Fill(tofX, elossX);
                dE_vs_tof_clusterSize1->Fill(tofY, elossY);

                TACAhit *hit = cluster->GetHit(0);
                Int_t crystal_id = hit->GetCrystalId();
                Double_t charge_clusterHit = hit->GetCharge();
                if (hit->IsValid() && charge_clusterHit > 0.02 && (crystal_id == 0 || calibCoeff.find(crystal_id) != calibCoeff.end()))
                {
                  Double_t charge_filling = (crystal_id == 0) ? charge_clusterHit : charge_clusterHit / calibCoeff.at(crystal_id);
                  if (ZX == 1)
                  {
                    CS1_Calo_Calibrated_Z1[crystal_id]->Fill(charge_filling);
                  }
                  else if (ZX == 2)
                  {
                    CS1_Calo_Calibrated_Z2[crystal_id]->Fill(charge_filling);
                  }
                  
                }

              }
              else if (nClusterHits == 2)
              {
                Z_clusterSize2->Fill(ZX);
                Z_clusterSize2->Fill(ZY);
                dE_vs_tof_clusterSize2->Fill(tofX, elossX);
                dE_vs_tof_clusterSize2->Fill(tofY, elossY);

                Double_t charge_sum_calibrated = 0;
                bool isZero = false;
                Double_t stored_charge = 0;
                for (int iclusterhit = 0; iclusterhit < nClusterHits; iclusterhit++)
                {
                  TACAhit *hit = cluster->GetHit(iclusterhit);
                  if (hit->IsValid())
                  {
                    Int_t crystal_id = hit->GetCrystalId();
                    Double_t charge_clusterHit = hit->GetCharge();
                    
                    if (crystal_id == 0 && charge_clusterHit > 0.02)
                    {
                        isZero = true;
                        charge_sum_calibrated += charge_clusterHit;
                    }
                    // The find() method returns an iterator to the matching element if found, or end() if not found.
                    else if ((calibCoeff.find(crystal_id) != calibCoeff.end()) && charge_clusterHit > 0.02)
                    {
                    
                        charge_clusterHit = charge_clusterHit / calibCoeff.at(crystal_id);
                        stored_charge += charge_clusterHit;
                    }
                  }

                }
                if (isZero)
                {
                    charge_sum_calibrated += stored_charge;
                    if (ZX == 1)
                    {
                      CS2_Calo_Calibrated_Z1->Fill(charge_sum_calibrated);
                    }
                    else if (ZX == 2)
                    {
                      CS2_Calo_Calibrated_Z2->Fill(charge_sum_calibrated);
                    }
                }
              }
          }
        }
      }

    }
  }

  gTAGroot.EndEventLoop();

  //SetTitleAndLabels(betaEloss, Form("#beta vs mass stopping power @ %d MeV/u beam energy", energy), "#beta", "#frac{dE}{#rho dx} [MeV cm^{2} g^{-1}]");
  //betaEloss->Write("ScatterPlot_beta_vs_dE");
  //SetTitleAndLabels(beta_vs_dE, Form("#beta vs mass stopping power @ %d MeV/u beam energy", energy), "#beta", "#frac{dE}{#rho dx} [MeV cm^{2} g^{-1}]");

  cout << endl
       << "Job Done!" << endl;

  fout->cd();
  fout->Write();
  fout->Close();

  return;
}

//-----------------------------------------------------------------------------
std::map<Int_t, std::map<Int_t, Double_t>> extractBarData() {
    std::map<Int_t, std::map<Int_t, Double_t>> barData;
    std::string filename = "calib/HIT2022/TATW_Energy_Calibration_perBar_4742.cal"; // File is hardcoded
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return barData;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue; // Skip headers

        std::istringstream iss(line);
        Int_t barId, p1, shoeLayer;
        Double_t p0;

        if (!(iss >> barId >> p0 >> p1 >> shoeLayer)) continue; // Skip invalid lines

        Int_t layer = shoeLayer;
        Int_t correctedBar = (shoeLayer == 0) ? barId : barId - 20; // Adjust for X layer

        barData[layer][correctedBar] = p0;
    }

    file.close();
    return barData;
}

std::map<Int_t, std::map<Int_t, Double_t>> extractTofData(Int_t energy) {
  std::map<Int_t, std::map<Int_t, Double_t>> tofData;
  std::string filename;
  if (energy == 100) filename = "calib/HIT2022/TATW_Tof_Calibration_perBar_4766.cal";
  else if (energy == 140) filename = "calib/HIT2022/TATW_Tof_Calibration_perBar_4801.cal";
  else if (energy == 200) filename = "calib/HIT2022/TATW_Tof_Calibration_perBar_4742.cal";
  else if (energy == 220) filename = "calib/HIT2022/TATW_Tof_Calibration_perBar_4828.cal";
  std::ifstream file(filename);

  if (!file.is_open()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return tofData;
  }

  std::string line;
  while (std::getline(file, line)) {
      if (line.empty() || line[0] == '#') continue; // Skip headers

      std::istringstream iss(line);
      Int_t barId, shoeLayer;
      Double_t Delta_t, sigma_t;

      if (!(iss >> barId >> Delta_t >> sigma_t >> shoeLayer)) continue; // Skip invalid lines

      Int_t layer = shoeLayer;
      Int_t correctedBar = (shoeLayer == 0) ? barId : barId - 20; // Adjust for X layer

      tofData[layer][correctedBar] = Delta_t;
  }

  file.close();
  return tofData;
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

void AdjustHistoRange(TH1D *Histo)
{
  Histo->GetXaxis()->SetRangeUser(Histo->GetBinLowEdge(Histo->FindFirstBinAbove()),
                                  Histo->GetBinLowEdge(Histo->FindLastBinAbove() + 1));
  Histo->GetXaxis()->SetRangeUser(Histo->GetBinLowEdge(Histo->FindFirstBinAbove()),
                                  Histo->GetBinLowEdge(Histo->FindLastBinAbove() + 1));
  Histo->GetXaxis()->SetRangeUser(Histo->GetBinLowEdge(Histo->FindFirstBinAbove()),
                                  Histo->GetBinLowEdge(Histo->FindLastBinAbove() + 1));
  return;
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

void BookHistograms(TDirectory *DirChargeElossLayerX, TDirectory *DirChargeElossLayerY,
                    TDirectory *DirToFLayerX, TDirectory *DirToFLayerY)
{

  // fpHisSeedMap = new TH1F(Form("msSeedMap%d", 4+1), Form("MSD - seed map for sensor %d", i+1), pGeoMap->GetStripsN(), 0, msdparGeo->GetStripsN());
  // AddHistogram(fpHisSeedMap);

  // fpHisStripMap = new TH1F(Form("msStripMap%d", 4+1), Form("MSD - strip map for sensor %d", i+1), pGeoMap->GetStripsN(), 0, msdparGeo->GetStripsN());
  // AddHistogram(fpHisStripMap);

  std::map<Int_t, Double_t> calibCoeff = extractCrystalData();

  Bar_ID_X = new TH1D("BarID_LayerX", "BarID LayerX (1 valid hit on both layers)", 200, 0, 19);
  Bar_ID_Y = new TH1D("BarID_LayerY", "BarID LayerY (1 valid hit on both layers)", 200, 0, 19);

  hHits_X = new TH1D("Hits_LayerX", "Number of Hits LayerX", 100, 0, 10);
  hHits_Y = new TH1D("Hits_LayerY", "Number of Hits LayerY", 100, 0, 10);
  h_nValidHits_X = new TH1D("nValidHits_LayerX", "Number of Valid Hits LayerX", 100, 0, 10);
  h_nValidHits_Y = new TH1D("nValidHits_LayerY", "Number of Valid Hits LayerY", 100, 0, 10);

  Z_clusterSize1 = new TH1D("Z_clusterSize1", "Z cluster size 1", 50, 0., 5.);
  Z_clusterSize2 = new TH1D("Z_clusterSize2", "Z cluster size 2", 50, 0., 5.);

  dE_vs_tof_clusterSize1 = new TH2D("dE_vs_tof_clusterSize1", "dE vs TOF Cluster Size 1", 1500, 5., 20., 2500, -5., 20.); // 0.01~ns/bin - 0.01 MeV/bin
  dE_vs_tof_clusterSize2 = new TH2D("dE_vs_tof_clusterSize2", "dE vs TOF Cluster Size 2", 1500, 5., 20., 2500, -5., 20.); // 0.01~ns/bin - 0.01 MeV/bin

  CS2_Calo_Calibrated_Z1 = new TH1D(Form("CS2_Calo_Calibrated_Z1"), Form("Calibrated Calo Charge - Cluster Size 2 - Z1"), 200, -0.5, 1.5);
  CS2_Calo_Calibrated_Z2 = new TH1D(Form("CS2_Calo_Calibrated_Z2"), Form("Calibrated Calo Charge - Cluster Size 2 - Z2"), 200, -0.5, 1.5);

  for (int icrystal = 0; icrystal < kModules * kCrysPerModule; icrystal++)
  {
    if (icrystal == 0 || calibCoeff.find(icrystal) != calibCoeff.end())
    {
      CS1_Calo_Calibrated_Z1[icrystal] = new TH1D(Form("CS1_Calo_Calibrated_Z1_%d", icrystal), Form("Calibrated Calo Charge - Cluster Size 1 - Z1 - Crystal ID %d", icrystal), 200, -0.5, 1.5);
      CS1_Calo_Calibrated_Z2[icrystal] = new TH1D(Form("CS1_Calo_Calibrated_Z2_%d", icrystal), Form("Calibrated Calo Charge - Cluster Size 1 - Z2 - Crystal ID %d", icrystal), 200, -0.5, 1.5);
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
