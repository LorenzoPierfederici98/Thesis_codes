// written by mtoppi 05/2023
// to be run in shoe/build/Reconstruction
// [ therein set-up the enviroment:  source ../setupFOOT.sh ]
// and run with (for example): root -l -b -q AnalyzeFOOT.cc++g\(\"../../../../rootfiles/outMC_16O_C_400_1_GSI.root\",1,10,\"\"\)
// sul tier1:  root -l -b -q AnalyzeFOOT.cc++g\(\"/storage/gpfs_data/foot/mtoppi/DataDecoded/CNAO2023/test.root\",0,1000,\"testAnaFOOT\",\"/storage/gpfs_data/foot/mtoppi/OutputMacro/\"\)

#include "AnalyzeTWFragMC.h"

// main
void AnalyzeTWFragMC(TString infile = "testMC.root", Bool_t isMax = kFALSE, Int_t nev = 10, TString outfile = "AnaFOOT.root", TString outDir = "/Users/marco/FOOT/Analisi/shoe/build/Reconstruction/OutputMacro")

{

  InitializeContainers();

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
  TDirectory *DirElossLayerX = fout->mkdir("ElossLayerX");
  TDirectory *DirElossLayerY = fout->mkdir("ElossLayerY");
  TDirectory *DirTofLayerX = fout->mkdir("TofLayerX");
  TDirectory *DirTofLayerY = fout->mkdir("TofLayerY");
  //fileBinTwCalib.open(Form("%s/calibTwInput_%s_run%d.dat", outDir.Data(), expName.Data(), runNumber), ios::out | ios::binary); // apre il file binario in SCRITTURA

  BookHistograms(DirElossLayerX, DirElossLayerY, DirTofLayerX, DirTofLayerY);

  TVecPair vPairWrongZ;
  vPairWrongZ.clear();

  map<Int_t, vector<TVector3>> pMap; // primaries
  map<Int_t, vector<TVector3>> nMap; // neutrons
  pMap.clear();
  nMap.clear();

  cout << "Beam Z:: " << parGeo->GetBeamPar().AtomicNumber << " A:: " << parGeo->GetBeamPar().AtomicMass << " ion::" << parGeo->GetBeamPar().Material << endl;

  cout << "TG center:: " << geoTrafo->GetTGCenter().z() << endl
       << " TG thickness:: " << parGeo->GetTargetPar().Size.z() << " TG material:: " << parGeo->GetTargetPar().Material << endl;

  printf("TW center =  (%2.1f,%2.1f,%2.1f)\n", GetTwCenter().x(), GetTwCenter().y(), GetTwCenter().z());

  cout << "theta angle of TW  acceptance::" << GetMaxAngle() * 180 / TMath::Pi() << endl;

  fActReader->Open(infile);
  gTAGroot.AddRequiredItem(fActReader);

  Int_t nentries = fActReader->NEvents();

  Int_t energy = static_cast<Int_t>(parGeo->GetBeamPar().Energy * 1000);  // MeV/u

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

  // Loop over the TTree
  gTAGroot.BeginEventLoop();

  Int_t ev = -1;

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

    /**
    if (calibTw)
    {
      if (!IncludeMC)
        FillBinaryForTwCalib(ev, nentries, runNumber);
      else
        FillBinaryForTwCalib_MC(ev, nentries, runNumber);
    }
    */

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

    for (int ihitX = 0; ihitX < nHitsX; ihitX++)
    {

      TATWhit *hitX = twNtuHit->GetHit(ihitX, (Int_t)LayerX);
      Double_t posAlongX = hitX->GetPosition(); // it provides the X coordinate
      Double_t posBarY = twparGeo->GetBarPosition((Int_t)LayerX, hitX->GetBar())[1];
      Double_t elossX = hitX->GetAmplitudeChA();
      Int_t barX = hitX->GetBar();

      // if ((energy == 100 && elossX > 1.5 && elossX < 13.) ||
      //     (energy == 140 && elossX > 1.1 && elossX < 11.) ||
      //     (energy == 200 && elossX > 0.9 && elossX < 9.) ||
      //     (energy == 220 && elossX > 0.8 && elossX < 8.5))
      if ((energy == 100 && elossX > 1.5) ||
          (energy == 140 && elossX > 1.1) ||
          (energy == 200 && elossX > 0.9) ||
          (energy == 220 && elossX > 0.8))    
          {
            nValidHitsX++;
            hitNumber_X = ihitX;
          }

      for (int ihitY = 0; ihitY < nHitsY; ihitY++)
      {

        TATWhit *hitY = twNtuHit->GetHit(ihitY, (Int_t)LayerY);
        Double_t posAlongY = hitY->GetPosition(); // it provides the Y coordinate
        Double_t posBarX = twparGeo->GetBarPosition((Int_t)LayerY, hitY->GetBar())[0];
        Double_t elossY = hitY->GetAmplitudeChA();
        Int_t barY = hitY->GetBar();

        if (!countedY[ihitY] && ((energy == 100 && elossY > 1.5) ||
          (energy == 140 && elossY > 1.1) ||
          (energy == 200 && elossY > 0.9) ||
          (energy == 220 && elossY > 0.8))) {
            nValidHitsY++;
            countedY[ihitY] = true; // Mark this hitY as counted
            hitNumber_Y = ihitY;
          }

        hTwMapPos->Fill(posAlongX, posAlongY);

        if (debug)
        {
          cout << " alongX::" << posAlongX << " alongY::" << posAlongY << " BarX::" << posBarX << " BarY::" << posBarY << endl;
          cout << " BarX_ID::::" << hitX->GetBar() << " BarY_ID::" << hitY->GetBar() << endl;
        }

        if (debug)
        {
          cout << " hitX_dE::" << hitX->GetEnergyLoss() << " hitY_dE::" << hitY->GetEnergyLoss() << endl;
          cout << " hitX_tof::" << hitX->GetToF() << " hitY_tof::" << hitY->GetToF() << endl;
          // cout<<" hitX_Z::"<<hitX->GetChargeZ()<<" hitY_Z::"<<hitY->GetChargeZ()<<endl;
        }
      }
    }

    hValidHits_X->Fill(nValidHitsX);
    hValidHits_Y->Fill(nValidHitsY);

    if (nValidHitsX == 1 && nValidHitsY == 1)
    {
      TATWhit *hitX = twNtuHit->GetHit(hitNumber_X, (Int_t)LayerX);
      Int_t barX = hitX->GetBar();
      Double_t elossX = hitX->GetAmplitudeChA();
      Double_t tofX = hitX->GetToF();

      TATWhit *hitY = twNtuHit->GetHit(hitNumber_Y, (Int_t)LayerY);
      Int_t barY = hitY->GetBar();
      Double_t elossY = hitY->GetAmplitudeChA();
      Double_t tofY = hitY->GetToF();

      if (ev % 10000 == 0) {
        cout << "1 hit in each layer" << endl;
        cout << "Energy [MeV/u]: " << energy << " elossX: " << elossX << " elossY: " << elossY << endl;
      }

      eloss_true_perBar[(Int_t)LayerX][barX]->Fill(elossX);
      eloss_true_perBar[(Int_t)LayerY][barY]->Fill(elossY);
      hToF[(Int_t)LayerX][barX]->Fill(tofX);
      hToF[(Int_t)LayerY][barY]->Fill(tofY);

      eloss_true_ch->Fill(elossX);
      eloss_true_ch->Fill(elossY);

    }

    if (nHitsX == 1 && nHitsY == 1)
    {
      TATWhit *hitX = twNtuHit->GetHit(0, (Int_t)LayerX);
      Double_t posAlongX = hitX->GetPosition(); // it provides the X coordinate
      Double_t posBarY = twparGeo->GetBarPosition((Int_t)LayerX, hitX->GetBar())[1];

      TATWhit *hitY = twNtuHit->GetHit(0, (Int_t)LayerY);
      Double_t posAlongY = hitY->GetPosition(); // it provides the Y coordinate
      Double_t posBarX = twparGeo->GetBarPosition((Int_t)LayerY, hitY->GetBar())[0];

      hResX_1Cross->Fill(posAlongX - posBarX);
      hResY_1Cross->Fill(posAlongY - posBarY);
      hTwMapPos_1Cross->Fill(posAlongX, posAlongY);
    }



    Int_t nHits = twNtuHit->GetHitN();
    if (debug)
      cout << " TWhits::" << nHits << endl;

    for (int ihit = 0; ihit < nHits; ihit++)
    {

      TATWhit *hit = twNtuHit->GetHit(ihit);
      if (!hit || !hit->IsValid())
        continue;

      Int_t bar = hit->GetBar();
      Int_t layer = hit->GetLayer();
      Int_t NmcTrk = hit->GetMcTracksN();
      Double_t eloss_ch = hit->GetAmplitudeChA();  // for MC it provides the true energy loss
      Double_t tof = hit->GetToF();
      Int_t Z = hit->GetChargeZ();

      if (debug)
        printf("twhit::%d  %s  bar::%d  Z::%d  eloss::%f  NmcTrk::%d\n", ihit, LayerName[(TLayer)layer].data(), bar, Z, eloss_ch, NmcTrk);

      dE_vs_tof_perBar[layer][bar]->Fill(tof, eloss_ch);
      
      noCuts_eloss_true_perBar[layer][bar]->Fill(eloss_ch);
      noCuts_hToF[layer][bar]->Fill(tof);
      eloss_true_nocuts->Fill(eloss_ch);

      //dE_vs_tof_perBar[layer][bar]->Fill(tof, eloss_ch);
      //dE_vs_tof[layer]->Fill(tof, eloss_ch);
      //heloss_all->Fill(eloss_ch);

      //eloss_true_ch->Fill(eloss_ch);
      //hToF[layer][bar]->Fill(tof);

      //if (layer == (Int_t)LayerX && bar >= CentralBarsID[0] && bar <= CentralBarsID[2])
      //{
      //  if (!IncludeMC)
      //    mapTrigHisto[(TrigID)trigID]->Fill(eloss);
      //}

      if (layer == (Int_t)LayerX)
      {
        Double_t posAlongX = hit->GetPosition();
        Double_t posBarY = twparGeo->GetBarPosition(layer, bar)[1];
        hTwPos[layer]->Fill(posAlongX, posBarY);
      }
      else if (layer == (Int_t)LayerY)
      {
        Double_t posAlongY = hit->GetPosition();
        Double_t posBarX = twparGeo->GetBarPosition(layer, bar)[0];
        hTwPos[layer]->Fill(posBarX, posAlongY);
      }

      for (int imctr = 0; imctr < NmcTrk; imctr++)
      {

        Int_t mcId = hit->GetMcIndex(imctr);
        Int_t mcTrkId = hit->GetMcTrackIdx(imctr);

        if (debug)
          cout << "mcId::" << mcId << "  mcTrkId::" << mcTrkId << endl;
      }
    }

    // Int_t NvtClus = vtxNtuCluster->GetClustersN(0);

    // for(int iclus=0; iclus<NvtClus; iclus++) {

    //   TAVTcluster*  vtClus =vtxNtuCluster->GetCluster(0,iclus);
    //   TVector3 clusPos = vtClus->GetPositionG();
    //   cout<<clusPos.X()<<"  "<<clusPos.Y()<<"  "<<clusPos.Z()<<endl;

    // }

    // Int_t NmsdHits = msdNtuHit->GetStripsN(4);

    // for(int istrip=0; istrip<NmsdHits; istrip++) {

    //   TAMSDntuHit*  msdHit =msdNtuHit->GetStrip(4,istrip);
    //   // TVector3 clusPos = msdClus->GetPositionG();
    //   // cout<<clusPos.X()<<"  "<<clusPos.Y()<<"  "<<clusPos.Z()<<endl;
    //   cout<<istrip<<"  charge::"<<msdHit->GetCharge()<<"  Strip::"<<msdHit->GetStrip()<<" stripId::"<<msdHit->GetSensorId()<<"  view::"<<msdHit->GetView()<<endl;

    //   charge[msdHit->GetStrip()]+=msdHit->GetCharge();
    // }

    // hChargeMSD->SetBinContent(msdHit->GetStrip(),charge);

    // Int_t NmsdClus = msdNtuCluster->GetClustersN(4);

    // for(int iclus=0; iclus<NmsdClus; iclus++) {

    //   TAMSDcluster*  msdClus =msdNtuCluster->GetCluster(4,iclus);
    //   TVector3 clusPos = msdClus->GetPositionG();
    //   cout<<clusPos.X()<<"  "<<clusPos.Y()<<"  "<<clusPos.Z()<<endl;

    // }
    /**
    int NtwPoints = twNtuPoint->GetPointsN();
    if (debug)
      cout << NtwPoints << " TW points" << endl;

    for (int ipt = 0; ipt < NtwPoints; ipt++)
    {

      if (debug)
        cout << "pointID::" << ipt << endl;

      Int_t cnt_all(0), cnt_good_1trk(0), cnt_wrong_1trk(0);

      TATWpoint *twpnt = twNtuPoint->GetPoint(ipt);
      TATWhit *hitX = twpnt->GetRowHit();
      TATWhit *hitY = twpnt->GetColumnHit();
      Int_t Z_hitX = hitX->GetChargeZ();
      Int_t Z_hitY = hitY->GetChargeZ();
      Int_t Z = twpnt->GetChargeZ();

      Int_t barX = hitX->GetBar();
      Int_t barY = hitY->GetBar();

      double eloss(-1);
      if (twpnt->GetMainLayer() == (Int_t)LayerX)
        eloss = hitX->GetEnergyLoss();
      else if (twpnt->GetMainLayer() == (Int_t)LayerY)
        eloss = hitY->GetEnergyLoss();

      Double_t tof = twpnt->GetMeanTof();
      // TVector3 glPntPos = twpnt->GetPositionGlb();

      TATWhit *hit, *other_hit;
      if (twpnt->GetMainLayer() == (Int_t)LayerX)
      {
        hit = (TATWhit *)hitX;
        other_hit = (TATWhit *)hitY;
      }
      else
      {
        hit = (TATWhit *)hitY;
        other_hit = (TATWhit *)hitX;
      }

      TVector3 pointLocPos = twpnt->GetPositionG();
      hTwMapPos_TWpntBin->Fill(pointLocPos.X(), pointLocPos.Y());
      hTwMapPos_TWpnt->Fill(hitX->GetPosition(), hitY->GetPosition());

      if (Z > 0 && Z < GetZbeam() + 1)
      {
        hTwMapPos_TWpnt_Z[Z - 1]->Fill(hitX->GetPosition(), hitY->GetPosition());
        if (twpnt->GetMainLayer() == (Int_t)LayerX)
          hResX[Z - 1]->Fill(pointLocPos.X() - hit->GetPosition());
        if (twpnt->GetMainLayer() == (Int_t)LayerY)
          hResY[Z - 1]->Fill(pointLocPos.Y() - hit->GetPosition());
      }
    }
    */

  } // close for loop over events

  gTAGroot.EndEventLoop();
  //eloss_true_ch->Add(eloss_true_X, eloss_true_Y);
  cout << endl
       << "Job Done!" << endl;

  fout->cd();
  fout->Write();
  fout->Close();

  //fileBinTwCalib.close();

  if (readBinTwCalibFile)
  {

    int pos1 = 1;
    int pos2 = 2;

    Int_t runN(-1);
    Int_t nEntries(-1);

    ifstream fileBinario("calibTwInput.dat", ios::in | ios::binary); // apre il file binario in LETTURA

    if (fileBinario.is_open())
    {
      fileBinario.seekg((pos1 - 1) * sizeof(Int_t), ios::beg); // si posiziona per leggere il centesimo intero
      fileBinario.read((char *)&runN, sizeof(runN));           // legge il centesimo intero
      fileBinario.seekg((pos2 - 1) * sizeof(Int_t), ios::beg); // si posiziona per leggere il centesimo intero
      fileBinario.read((char *)&nEntries, sizeof(nEntries));   // legge il centesimo intero
      fileBinario.close();
    }

    cout << "runN::" << runN << " entries::" << nEntries << endl;

    // Int_t barID_r;
    // Int_t layerID_r;
    // Double_t rawTof_r;
    // Double_t rawEnergy_r;
    // fileBin.read((char*) &layerID_r, sizeof(Int_t)); //scrive i byte della variabile n nel file binario
    // fileBin.read((char*) &barID_r, sizeof(Int_t)); //scrive i byte della variabile n nel file binario
    // fileBin.read((char*) &rawTof_r, sizeof(Double_t)); //scrive i byte della variabile n nel file binario
    // fileBin.read((char*) &rawEnergy_r, sizeof(Double_t)); //scrive i byte della variabile n nel file binario
    // cout<<"lay::"<<layerID_r<<" barID_r::"<<barID_r<<" rawTof_r::"<<rawTof_r<<" rawEnergy_r::"<<rawEnergy_r<<endl;
  }

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

void BookHistograms(TDirectory *DirElossLayerX, TDirectory *DirElossLayerY, TDirectory *DirTofLayerX, TDirectory *DirTofLayerY)
{

  // fpHisSeedMap = new TH1F(Form("msSeedMap%d", 4+1), Form("MSD - seed map for sensor %d", i+1), pGeoMap->GetStripsN(), 0, msdparGeo->GetStripsN());
  // AddHistogram(fpHisSeedMap);

  // fpHisStripMap = new TH1F(Form("msStripMap%d", 4+1), Form("MSD - strip map for sensor %d", i+1), pGeoMap->GetStripsN(), 0, msdparGeo->GetStripsN());
  // AddHistogram(fpHisStripMap);

  eloss_true_ch = new TH1D("Eloss_true_ch", "Eloss true MC ch", 120, -2., 20.);
  eloss_true_nocuts = new TH1D("Eloss_true_ch_nocuts", "Eloss true MC ch no cuts", 120, -2., 20.);
  //eloss_1hit_sum = new TH1D("Eloss_1hit_sum", "Eloss true MC 1 hit", 120, -2., 20.);

  hHits_X = new TH1D("Hits_LayerX", "Number of Hits LayerX", 100, 0, 10);
  hHits_Y = new TH1D("Hits_LayerY", "Number of Hits LayerY", 100, 0, 10);
  hValidHits_X = new TH1D("nValidHits_LayerX", "Number of Valid Hits LayerX", 100, 0, 10);
  hValidHits_Y = new TH1D("nValidHits_LayerY", "Number of Valid Hits LayerY", 100, 0, 10);

  for (int ilay = 0; ilay < kLayers; ilay++)
  {
    for (int ibar = 0; ibar < (int)nBarsPerLayer; ibar++)
    {
      //eloss_true_MC[ilay][ibar] = new TH1D(Form("Eloss_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), Form("Eloss true MC %s bar%d", LayerName[(TLayer)ilay].data(), ibar), 120, -2., 20.);
      hToF[ilay][ibar] = new TH1D(Form("ToF_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), Form("ToF %s bar %d", LayerName[(TLayer)ilay].data(), ibar), 140, 6., 20.);
      noCuts_hToF[ilay][ibar] = new TH1D(Form("noCuts_ToF_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), Form("No cuts ToF %s bar %d", LayerName[(TLayer)ilay].data(), ibar), 140, 6., 20.);
      eloss_true_perBar[ilay][ibar] = new TH1D(Form("Eloss_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), Form("Eloss true MC %s bar%d", LayerName[(TLayer)ilay].data(), ibar), 120, -2., 20.);
      noCuts_eloss_true_perBar[ilay][ibar] = new TH1D(Form("noCuts_Eloss_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), Form("No cuts Eloss true MC %s bar%d", LayerName[(TLayer)ilay].data(), ibar), 120, -2., 20.);
      if (ilay == (Int_t)LayerX) {
        hToF[ilay][ibar]->SetDirectory(DirTofLayerX);
        eloss_true_perBar[ilay][ibar]->SetDirectory(DirElossLayerX);
        noCuts_eloss_true_perBar[ilay][ibar]->SetDirectory(DirElossLayerX);
        noCuts_hToF[ilay][ibar]->SetDirectory(DirTofLayerX);
      }
      else {
        hToF[ilay][ibar]->SetDirectory(DirTofLayerY);
        eloss_true_perBar[ilay][ibar]->SetDirectory(DirElossLayerY);
        noCuts_eloss_true_perBar[ilay][ibar]->SetDirectory(DirElossLayerY);
        noCuts_hToF[ilay][ibar]->SetDirectory(DirTofLayerY);
      }

      dE_vs_tof_perBar[ilay][ibar] = new TH2D(Form("dE_vs_tof_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), Form("dE_vs_tof_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), 2500, 5., 30., 20, 0., 2.); // 10~ps/bin - 0.1 MeV/bin
    }
    //dE_vs_tof[ilay] = new TH2D(Form("dE_vs_tof_%s", LayerName[(TLayer)ilay].data()), Form("dE_vs_tof_%s", LayerName[(TLayer)ilay].data()), 25000, 5., 30., 1200, 0., 120.); // 1~ps/bin - 0.1 MeV/bin

    if (ilay == (Int_t)LayerX)
      hTwPos[ilay] = new TH2D(Form("hTwPos_%s", LayerName[(TLayer)ilay].data()), Form("hTwPos_%s", LayerName[(TLayer)ilay].data()), 220, -22., 22., 20, -20., 20.); // 2 mm/bin - 2 cm/bin
    else if (ilay == (Int_t)LayerY)
      hTwPos[ilay] = new TH2D(Form("hTwPos_%s", LayerName[(TLayer)ilay].data()), Form("hTwPos_%s", LayerName[(TLayer)ilay].data()), 20, -20., 20., 220, -22., 22.); // 2 cm/bin - 2 mm/bin
  }

  //heloss_all = new TH1D("Eloss_all", "Eloss_all", 1200, 0., 120.);

  for (auto itr = mapTrig.begin(); itr != mapTrig.end(); ++itr)
  {
    TrigID trig_id = itr->first;
    TString TrigName = itr->second;
    cout << "TrigID::" << itr->first << " TrigName::" << TrigName << endl;
    // heloss[itrig] = new TH1D(Form("Eloss_%s",TrigName.Data()),Form("Eloss_%s",TrigName.Data()),1200,0.,120.);
    mapTrigHisto[trig_id] = new TH1D(Form("Eloss_%s", TrigName.Data()), Form("Eloss_%s", TrigName.Data()), 1200, 0., 120.);
  }

  hTwMapPos_TWpntBin = new TH2D("hTwMapPos_TWpntBin", "hTwMapPos_TWpntBin", 22, -22., 22., 22, -22., 22.); // 2 cm/bin - 2 cm/bin
  hTwMapPos_TWpnt = new TH2D("hTwMapPos_TWpnt", "hTwMapPos_TWpnt", 220, -22., 22., 220, -22., 22.);        // 2 mm/bin - 2 mm/bin

  hTwMapPos = new TH2D("hTwMapPos", "hTwMapPos", 220, -22., 22., 220, -22., 22.); // 2 mm/bin - 2 mm/bin

  hTwMapPos_1Cross = new TH2D("hTwMapPos_1Cross", "hTwMapPos_1Cross", 220, -22., 22., 220, -22., 22.); // 2 mm/bin - 2 mm/bin
  hResX_1Cross = new TH1D("hResX_1Cross", "hResX_1Cross", 400, -20., 20.);
  hResY_1Cross = new TH1D("hResY_1Cross", "hResY_1Cross", 400, -20., 20.);

  for (int ichg = 0; ichg < GetZbeam(); ichg++)
  {

    hTwMapPos_TWpnt_Z[ichg] = new TH2D(Form("hTwMapPos_TWpnt_Z%d", ichg + 1), Form("hTwMapPos_TWpnt_Z%d", ichg + 1), 220, -22., 22., 220, -22., 22.); // 2 mm/bin - 2 mm/bin

    hResX[ichg] = new TH1D(Form("hResX_Z%d", ichg + 1), Form("hResX_Z%d", ichg + 1), 400, -20., 20.);
    hResY[ichg] = new TH1D(Form("hResY_Z%d", ichg + 1), Form("hResY_Z%d", ichg + 1), 400, -20., 20.);
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
void FillBinaryForTwCalib(Int_t eventN, Int_t nEntries, Int_t runN)
{

  if (eventN == 0)
    cout << "Fill Binaries File for Data" << endl;

  if (twNtuHit->GetHitN() < 1)
  {
    Info("FillBinaryForTwCalib()", "N of TW hits is::%d...skip event\n", twNtuHit->GetHitN());
    return;
  }

  Int_t nHitsX = twNtuHit->GetHitN((Int_t)LayerX);
  Int_t nHitsY = twNtuHit->GetHitN((Int_t)LayerY);

  if (nHitsX == 0 || nHitsY == 0)
  {
    Info("FillBinaryForTwCalib()", "N of TW hitsX is::%d and hitsY is::%d...skip event\n", nHitsX, nHitsY);
    return;
  }

  map<Int_t, TATWhit *> fmapHitX;
  map<Int_t, TATWhit *> fmapHitY;
  fmapHitX.clear();
  fmapHitY.clear();

  if (isTwScan)
  {

    for (int idx = 0; idx < nHitsX; idx++)
    {

      TATWhit *hitX = twNtuHit->GetHit(idx, LayerX);

      fmapHitX[idx] = hitX;
    }

    for (int idy = 0; idy < nHitsY; idy++)
    {

      TATWhit *hitY = twNtuHit->GetHit(idy, LayerY);

      fmapHitY[idy] = hitY;
    }

    for (auto it1 = fmapHitX.begin(); it1 != fmapHitX.end(); ++it1)
    {

      Int_t idx = it1->first;
      TATWhit *hitX = it1->second;
      Int_t layerX = hitX->GetLayer();
      Int_t barX = hitX->GetBar();

      if (barX != 9)
        continue;

      for (auto it2 = fmapHitY.begin(); it2 != fmapHitY.end(); ++it2)
      {

        Int_t idx = it2->first;
        TATWhit *hitY = it2->second;
        Int_t layerY = hitY->GetLayer();
        Int_t barY = hitY->GetBar();

        Double_t rawTof = hitY->GetTime();

        Double_t chargeA = hitY->GetChargeChA();
        Double_t chargeB = hitY->GetChargeChB();

        Double_t rawEnergy = sqrt(chargeA * chargeB);

        dE_vs_tof_perBar[layerY][barY]->Fill(rawTof, rawEnergy);
      }
    }

    for (auto it1 = fmapHitY.begin(); it1 != fmapHitY.end(); ++it1)
    {

      Int_t idx = it1->first;
      TATWhit *hitY = it1->second;
      Int_t layerY = hitY->GetLayer();
      Int_t barY = hitY->GetBar();

      if (barY != 9)
        continue;

      for (auto it2 = fmapHitX.begin(); it2 != fmapHitX.end(); ++it2)
      {

        Int_t idx = it2->first;
        TATWhit *hitX = it2->second;
        Int_t layerX = hitX->GetLayer();
        Int_t barX = hitX->GetBar();

        Double_t rawTof = hitX->GetTime();

        Double_t chargeA = hitX->GetChargeChA();
        Double_t chargeB = hitX->GetChargeChB();

        Double_t rawEnergy = sqrt(chargeA * chargeB);

        dE_vs_tof_perBar[layerX][barX]->Fill(rawTof, rawEnergy);
      }
    }
  }

  Int_t nHits = twNtuHit->GetHitN();

  if (fileBinTwCalib.is_open())
  {

    // if(eventN==0){
    //   fileBinTwCalib.write((char*) &runN, sizeof(Int_t)); //scrive i byte della variabile n nel file binario
    //   fileBinTwCalib.write((char*) &nEntries, sizeof(Int_t)); //scrive i byte della variabile n nel file binario
    // }
    fileBinTwCalib.write((char *)&eventN, sizeof(Int_t)); // scrive i byte della variabile n nel file binario
    fileBinTwCalib.write((char *)&nHits, sizeof(Int_t));  // scrive i byte della variabile n nel file binario
  }

  for (int ihit = 0; ihit < nHits; ihit++)
  {

    TATWhit *hit = twNtuHit->GetHit(ihit);

    Int_t barID = hit->GetBar();
    Int_t layerID = hit->GetLayer();

    Double_t rawTof = hit->GetTime();

    Double_t chargeA = hit->GetChargeChA();
    Double_t chargeB = hit->GetChargeChB();

    Double_t rawEnergy = sqrt(chargeA * chargeB);

    if (!isTwScan)
    {
      dE_vs_tof_perBar[layerID][barID]->Fill(rawTof, rawEnergy);
      dE_vs_tof[layerID]->Fill(rawTof, rawEnergy);
    }

    if (fileBinTwCalib.is_open())
    {

      // fileBinTwCalib<<setw(12)<<"#Layer"<<setw(12)<<"barID "<<setw(12)<<"rawTof"<<setw(12)<<"rawEnergy"<<endl;
      // fileBinTwCalib<<setw(12)<<layerID<<setw(12)<<barID<<setw(12)<<rawTof<<setw(12)<<rawEnergy<<endl;
      // fileBinTwCalib<<"#-+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-"<<endl;

      fileBinTwCalib.write((char *)&layerID, sizeof(Int_t)); // scrive i byte della variabile n nel file binario
      fileBinTwCalib.write((char *)&barID, sizeof(Int_t));   // scrive i byte della variabile n nel file binario

      fileBinTwCalib.write((char *)&rawTof, sizeof(Double_t));    // scrive i byte della variabile n nel file binario
      fileBinTwCalib.write((char *)&rawEnergy, sizeof(Double_t)); // scrive i byte della variabile n nel file binario
    }
  }

  return;
}
//-----------------------------------------------------------------------------
void FillBinaryForTwCalib_MC(Int_t eventN, Int_t nEntries, Int_t runN)
{

  if (eventN == 0)
    cout << "Fill Binaries File for MC" << endl;

  if (twNtuHit->GetHitN() < 1)
  {
    Info("FillBinaryForTwCalib()", "N of TW hits is::%d...skip event\n", twNtuHit->GetHitN());
    return;
  }

  Int_t nHitsX = twNtuHit->GetHitN((Int_t)LayerX);
  Int_t nHitsY = twNtuHit->GetHitN((Int_t)LayerY);

  if (nHitsX == 0 || nHitsY == 0)
  {
    Info("FillBinaryForTwCalib()", "N of TW hitsX is::%d and hitsY is::%d...skip event\n", nHitsX, nHitsY);
    return;
  }

  Int_t nHits = twNtuHit->GetHitN();
  Int_t nHits_cut(0);
  for (int ihit = 0; ihit < nHits; ihit++)
  {

    TATWhit *hit = twNtuHit->GetHit(ihit);

    Int_t barID = hit->GetBar();
    Int_t layerID = hit->GetLayer();

    if (stNtuHit->GetHitsN() < 1)
    {
      Info("FillBinaryForTwCalib()", "N of ST hits is::%d...skip event\n", stNtuHit->GetHitsN());
      return;
    }

    if (hit->GetMcTracksN() == 1)
    {

      Int_t trkMcId = hit->GetMcTrackIdx(0);

      TAMCpart *mctrk = mcNtuPart->GetTrack(trkMcId);
      int mothId = mctrk->GetMotherID();
      int reg = mctrk->GetRegion(); // region where track is generated

      if (reg == parGeo->GetRegTarget() && mothId == kPrimaryID)
      {
        nHits_cut++;
      }
    }
  }

  if (fileBinTwCalib.is_open())
  {

    if (nHits_cut > 0)
    {
      fileBinTwCalib.write((char *)&eventN, sizeof(Int_t));    // scrive i byte della variabile n nel file binario
      fileBinTwCalib.write((char *)&nHits_cut, sizeof(Int_t)); // scrive i byte della variabile n nel file binario
      if (debug)
        cout << "eventN::" << eventN << "  nHits::" << nHits << "  nHits_cut::" << nHits_cut << endl;
    }
  }

  for (int ihit = 0; ihit < nHits; ihit++)
  {

    TATWhit *hit = twNtuHit->GetHit(ihit);

    if (stNtuHit->GetHitsN() < 1)
    {
      Info("FillBinaryForTwCalib()", "N of ST hits is::%d...skip event\n", stNtuHit->GetHitsN());
      return;
    }

    Int_t barID = hit->GetBar();
    Int_t layerID = hit->GetLayer();

    Double_t timeSTtrue = stNtuHit->GetHit(0)->GetCharge(); // in GetCharge() stored the true ST time in ns
    Double_t elossTWtrue = hit->GetAmplitudeChA();          // in GetAmplitudeChA() stored the true TW eloss: in MeV it has a meaning only for noPileUp production
    Double_t timeTWtrue = hit->GetAmplitudeChB() * TAGgeoTrafo::PsToNs();
    ; // in GetAmplitudeChB() stored the true TW time: it has a meaning only for noPileUp production

    Double_t trueTof = (timeTWtrue - timeSTtrue);

    if (hit->GetMcTracksN() == 1)
    {

      Int_t trkMcId = hit->GetMcTrackIdx(0);

      TAMCpart *mctrk = mcNtuPart->GetTrack(trkMcId);
      int mothId = mctrk->GetMotherID();
      int reg = mctrk->GetRegion(); // region where track is generated

      if (debug)
        cout << "trkId::" << trkMcId << " Eloss::" << elossTWtrue << " Tof::" << trueTof << " timeST::" << timeSTtrue << " timeTW::" << timeTWtrue << " mothId::" << mothId << " reg::" << reg << endl;

      if (reg == parGeo->GetRegTarget() && mothId == kPrimaryID)
      {

        if (debug)
          cout << "in TG-->trkId::" << trkMcId << " Eloss::" << elossTWtrue << " Tof::" << trueTof << " mothId::" << mothId << " reg::" << reg << endl;

        dE_vs_tof_perBar[layerID][barID]->Fill(trueTof, elossTWtrue);
        dE_vs_tof[layerID]->Fill(trueTof, elossTWtrue);
        heloss_all->Fill(elossTWtrue);

        if (fileBinTwCalib.is_open())
        {

          fileBinTwCalib.write((char *)&layerID, sizeof(Int_t)); // scrive i byte della variabile n nel file binario
          fileBinTwCalib.write((char *)&barID, sizeof(Int_t));   // scrive i byte della variabile n nel file binario

          fileBinTwCalib.write((char *)&trueTof, sizeof(Double_t));     // scrive i byte della variabile n nel file binario
          fileBinTwCalib.write((char *)&elossTWtrue, sizeof(Double_t)); // scrive i byte della variabile n nel file binario

          if (debug)
            cout << "fill-->layerID::" << layerID << " barID::" << barID << " Eloss::" << elossTWtrue << " Tof::" << trueTof << endl;
        }
      }
    }
  }

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
