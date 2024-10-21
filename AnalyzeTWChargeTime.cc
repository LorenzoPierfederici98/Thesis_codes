// written by mtoppi 05/2023
// to be run in shoe/build/Reconstruction
// [ therein set-up the enviroment:  source ../setupFOOT.sh ]
// and run with (for example): root -l -b -q AnalyzeFOOT.cc++g\(\"../../../../rootfiles/outMC_16O_C_400_1_GSI.root\",1,10,\"\"\)
// sul tier1:  root -l -b -q AnalyzeFOOT.cc++g\(\"/storage/gpfs_data/foot/mtoppi/DataDecoded/CNAO2023/test.root\",0,1000,\"testAnaFOOT\",\"/storage/gpfs_data/foot/mtoppi/OutputMacro/\"\)

#include "AnalyzeTWChargeTime.h"

// main
void AnalyzeTWChargeTime(TString infile = "testMC.root", Bool_t isMax = kFALSE, Int_t nev = 10, TString outfile = "AnaFOOT.root", TString outDir = "/Users/marco/FOOT/Analisi/shoe/build/Reconstruction/OutputMacro")

{

  //InitializeContainers();

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

  const char *insertStr = "_Fit";
  const char *dotPos = strrchr(outfile, '.');

  size_t prefixLength = dotPos - outfile; // Length of string before ".root"
  size_t originalLength = strlen(outfile);
  size_t insertLength = strlen(insertStr);
  size_t newLength = originalLength + insertLength;
  char *fitFilename = new char[newLength + 1];

  strncpy(fitFilename, outfile, prefixLength);
  // Append the insert string
  strcpy(fitFilename + prefixLength, insertStr);
  // Append the rest of the original string (from ".root" onward)
  strcpy(fitFilename + prefixLength + insertLength, dotPos);

  TFile *fitFile = new TFile(fitFilename, "RECREATE");
  delete[] fitFilename;
  string beamEnergyStr = to_string(parGeo->GetBeamPar().Energy * 1000);
  TObjString objString(beamEnergyStr.c_str());
  TObjString materialObj(parGeo->GetBeamPar().Material);
  fout->cd();
  objString.Write("BeamEnergyInfo");
  fitFile->cd();                     // Ensure the correct file is active
  objString.Write("BeamEnergyInfo"); // Write the primary beam energy and type of particle in the fit file
  materialObj.Write("IonInfo");

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

    Int_t nHitsX = twNtuHit->GetHitN((Int_t)LayerX);
    Int_t nHitsY = twNtuHit->GetHitN((Int_t)LayerY);

    
    if (debug)
      cout << " TWhits X::" << nHitsX << " Y::" << nHitsY << endl;

    for (int ihitX = 0; ihitX < nHitsX; ihitX++)
    {

      TATWhit *hitX = twNtuHit->GetHit(ihitX, (Int_t)LayerX);
      Double_t posAlongX = hitX->GetPosition(); // it provides the X coordinate
      Double_t posBarY = twparGeo->GetBarPosition((Int_t)LayerX, hitX->GetBar())[1];
      Double_t barX = hitX->GetBar();
      if (hitX->GetEnergyLoss() > 0)
        Bar_ID_X->Fill(barX);
        PosX->Fill(posAlongX);

      for (int ihitY = 0; ihitY < nHitsY; ihitY++)
      {

        TATWhit *hitY = twNtuHit->GetHit(ihitY, (Int_t)LayerY);
        Double_t posAlongY = hitY->GetPosition(); // it provides the Y coordinate
        Double_t posBarX = twparGeo->GetBarPosition((Int_t)LayerY, hitY->GetBar())[0];
        Double_t barY = hitY->GetBar();
        if (hitY->GetEnergyLoss() > 0)
          Bar_ID_Y->Fill(barY);
          PosY->Fill(posAlongX);
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

      if (!hit)
        continue;

      Int_t bar = hit->GetBar();
      Int_t layer = hit->GetLayer();
      Int_t NmcTrk = hit->GetMcTracksN();
      Double_t eloss = hit->GetEnergyLoss();
      Double_t tof = hit->GetToF();
      Int_t Z = hit->GetChargeZ();
      Double_t chargeA = hit->GetChargeChA();
      Double_t chargeB = hit->GetChargeChB();
      if (chargeA < 0)
        chargeA = 0;
      if (chargeB < 0)
        chargeB = 0;
      Double_t chargeBar = sqrt(chargeA * chargeB);
      Double_t timeA = hit->GetTimeChA();
      Double_t timeB = hit->GetTimeChB();
      Double_t timeBar = 0.5 * (timeA + timeB);

      if (debug)
        printf("twhit::%d  %s  bar::%d  Z::%d  eloss::%f  NmcTrk::%d\n", ihit, LayerName[(TLayer)layer].data(), bar, Z, eloss, NmcTrk);

      if (!calibTw)
      {
        // dE_vs_tof_perBar[layer][bar]->Fill(tof,eloss);
        dE_vs_tof[layer]->Fill(tof, eloss);
        hToF->Fill(tof);
        if (bar == 8 || bar == 9 || bar == 10) {
          hToF_CentralBars->Fill(tof);
        }
        heloss_all->Fill(eloss);
        ChargeA_perBar[layer][bar]->Fill(chargeA);
        ChargeB_perBar[layer][bar]->Fill(chargeB);
        Charge_perBar[layer][bar]->Fill(chargeBar);
        TimeA_perBar[layer][bar]->Fill(timeA);
        TimeB_perBar[layer][bar]->Fill(timeB);
        Time_perBar[layer][bar]->Fill(timeBar);

        ChargeA_perBar[layer][bar]->GetXaxis()->SetRangeUser(ChargeA_perBar[layer][bar]->GetBinLowEdge(ChargeA_perBar[layer][bar]->FindFirstBinAbove()), 
                                   ChargeA_perBar[layer][bar]->GetBinLowEdge(ChargeA_perBar[layer][bar]->FindLastBinAbove() + 1));
        ChargeB_perBar[layer][bar]->GetXaxis()->SetRangeUser(ChargeB_perBar[layer][bar]->GetBinLowEdge(ChargeB_perBar[layer][bar]->FindFirstBinAbove()), 
                                   ChargeB_perBar[layer][bar]->GetBinLowEdge(ChargeB_perBar[layer][bar]->FindLastBinAbove() + 1));
        Charge_perBar[layer][bar]->GetXaxis()->SetRangeUser(Charge_perBar[layer][bar]->GetBinLowEdge(Charge_perBar[layer][bar]->FindFirstBinAbove()), 
                                   Charge_perBar[layer][bar]->GetBinLowEdge(Charge_perBar[layer][bar]->FindLastBinAbove() + 1));

        TimeA_perBar[layer][bar]->GetXaxis()->SetRangeUser(TimeA_perBar[layer][bar]->GetBinLowEdge(TimeA_perBar[layer][bar]->FindFirstBinAbove()), 
                                   TimeA_perBar[layer][bar]->GetBinLowEdge(TimeA_perBar[layer][bar]->FindLastBinAbove() + 1));
        TimeB_perBar[layer][bar]->GetXaxis()->SetRangeUser(TimeB_perBar[layer][bar]->GetBinLowEdge(TimeB_perBar[layer][bar]->FindFirstBinAbove()), 
                                   TimeB_perBar[layer][bar]->GetBinLowEdge(TimeB_perBar[layer][bar]->FindLastBinAbove() + 1));
        Time_perBar[layer][bar]->GetXaxis()->SetRangeUser(Time_perBar[layer][bar]->GetBinLowEdge(Time_perBar[layer][bar]->FindFirstBinAbove()), 
                                   Time_perBar[layer][bar]->GetBinLowEdge(Time_perBar[layer][bar]->FindLastBinAbove() + 1));
      }


      if (layer == (Int_t)LayerX)
      {
        Double_t posAlongX = hit->GetPosition();
        if (bar == 9){
          PosX_Bar9->Fill(posAlongX);
        }
        Double_t posBarY = twparGeo->GetBarPosition(layer, bar)[1];
        hTwPos[layer]->Fill(posAlongX, posBarY);
      }
      else if (layer == (Int_t)LayerY)
      {
        Double_t posAlongY = hit->GetPosition();
        if (bar == 9){
          PosY_Bar9->Fill(posAlongY);
        }
        Double_t posBarX = twparGeo->GetBarPosition(layer, bar)[0];
        hTwPos[layer]->Fill(posBarX, posAlongY);
      }
    }
  }

  gTAGroot.EndEventLoop();

  TFitResultPtr t = hToF->Fit("gaus", "QS");

  if (t->IsValid())
      {
        t->Write(Form("TW_Fit"));
      }
  else
      {
        std::cerr << "Fit TW failed" << std::endl;
      }

  TFitResultPtr t_c = hToF_CentralBars->Fit("gaus", "QS");

  if (t_c->IsValid())
      {
        t_c->Write(Form("TW_Fit_CentralBars"));
      }
  else
      {
        std::cerr << "Fit TW Central Bars failed" << std::endl;
      }
  

  const Double_t fractionThreshold = 0.009; // Fraction of total entries required (e.g., 0.1 for 10%)
  Int_t fitThreshold = static_cast<Int_t>(fractionThreshold * nentries);

  for (int ilay = 0; ilay < kLayers; ilay++)
  {
    for (int ibar = 0; ibar < 20; ibar++) // Loop through all bars from 0 to 19
    {
      // Get the number of entries in the histograms
      int chargeEntries = Charge_perBar[ilay][ibar]->GetEntries();
      int timeEntries = Time_perBar[ilay][ibar]->GetEntries();

      // Check if Charge_perBar histogram has sufficient entries
      if (chargeEntries >= fitThreshold)
      {
        TFitResultPtr c = Charge_perBar[ilay][ibar]->Fit("gaus", "QS", " ", 2., 10.);

        if (c->IsValid())
        {
          c->Write(Form("FitResult_Charge_layer%d_bar%d", ilay, ibar));
        }
        else
        {
          std::cerr << "Fit failed for Charge in layer " << ilay << ", bar " << ibar << std::endl;
        }
      }

      // Check if Time_perBar histogram has sufficient entries
      // if (timeEntries >= fitThreshold)
      // {
      //   TFitResultPtr t = Time_perBar[ilay][ibar]->Fit("gaus", "QS");

      //   if (t->IsValid())
      //   {
      //     t->Write(Form("FitResult_Time_layer%d_bar%d", ilay, ibar));
      //   }
      //   else
      //   {
      //     std::cerr << "Fit failed for Time in layer " << ilay << ", bar " << ibar << std::endl;
      //   }
      // }
    }
  }

  fitFile->Close();
  delete fitFile;
  cout << endl
       << "Job Done!" << endl;

  fout->cd();
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
  int modules[7] = {4, 5, 3, 6, 2, 7, 1};  // module enumeration following the excel file

  PosX = new TH1D("Hit_Pos_LayerX", "Hit_Pos_LayerX", 220, -22., 22.);  //2 mm/bin
  PosY = new TH1D("Hit_Pos_LayerY", "Hit_Pos_LayerY", 220, -22., 22.);
  PosX_Bar9 = new TH1D("Hit_Pos_LayerX_Bar9", "Hit_Pos_LayerX_Bar9", 220, -22., 22.);
  PosY_Bar9 = new TH1D("Hit_Pos_LayerY_Bar9", "Hit_Pos_LayerY_Bar9", 220, -22., 22.);
  Bar_ID_X = new TH1D("BarID_LayerX", "BarID_LayerX", 200, 0, 19);
  Bar_ID_Y = new TH1D("BarID_LayerY", "BarID_LayerY", 200, 0, 19);

  hToF = new TH1D("ToF_TW", "ToF_TW", 220, 6., 12.);
  hToF_CentralBars = new TH1D("ToF_TW_CentralBars", "ToF_TW_CentralBars", 220, 6., 12.);

  for (int ilay = 0; ilay < kLayers; ilay++)
  {

    dE_vs_tof[ilay] = new TH2D(Form("dE_vs_tof_%s", LayerName[(TLayer)ilay].data()), Form("dE_vs_tof_%s", LayerName[(TLayer)ilay].data()), 25000, 5., 30., 1200, 0., 120.); // 1~ps/bin - 0.1 MeV/bin

    for (int ibar = 0; ibar < (int)nBarsPerLayer; ibar++)
    {
      // dE_vs_tof_perBar[ilay][ibar] = new TH2D(Form("dE_vs_tof_%s_bar%d",LayerName[(TLayer)ilay].data(),ibar),Form("dE_vs_tof_%s_bar%d",LayerName[(TLayer)ilay].data(),ibar),25000,5.,30.,1200,0.,120.);  // 1~ps/bin - 0.1 MeV/bin
      ChargeA_perBar[ilay][ibar] = new TH1D(Form("Charge_ChA_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), Form("Charge_ChA_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), 120, -2., 12.);
      ChargeB_perBar[ilay][ibar] = new TH1D(Form("Charge_ChB_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), Form("Charge_ChB_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), 120, -2., 12.);
      Charge_perBar[ilay][ibar] = new TH1D(Form("Charge_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), Form("Charge_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), 120, -2., 12.);
      TimeA_perBar[ilay][ibar] = new TH1D(Form("Time_ChA_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), Form("Time_ChA_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), 150, -100., 400.);
      TimeB_perBar[ilay][ibar] = new TH1D(Form("Time_ChB_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), Form("Time_ChB_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), 150, -100., 400.);
      Time_perBar[ilay][ibar] = new TH1D(Form("Time_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), Form("Time_%s_bar%d", LayerName[(TLayer)ilay].data(), ibar), 150, -100., 400.);
    }
    if (ilay == (Int_t)LayerX)
      hTwPos[ilay] = new TH2D(Form("hTwPos_%s", LayerName[(TLayer)ilay].data()), Form("hTwPos_%s", LayerName[(TLayer)ilay].data()), 220, -22., 22., 20, -20., 20.); // 2 mm/bin - 2 cm/bin
    else if (ilay == (Int_t)LayerY)
      hTwPos[ilay] = new TH2D(Form("hTwPos_%s", LayerName[(TLayer)ilay].data()), Form("hTwPos_%s", LayerName[(TLayer)ilay].data()), 20, -20., 20., 220, -22., 22.); // 2 cm/bin - 2 mm/bin
  }

  heloss_all = new TH1D("Eloss_all", "Eloss_all", 1200, 0., 120.);

  hTwMapPos_TWpntBin = new TH2D("hTwMapPos_TWpntBin", "hTwMapPos_TWpntBin", 22, -22., 22., 22, -22., 22.); // 2 cm/bin - 2 cm/bin
  hTwMapPos_TWpnt = new TH2D("hTwMapPos_TWpnt", "hTwMapPos_TWpnt", 220, -22., 22., 220, -22., 22.);        // 2 mm/bin - 2 mm/bin

  hTwMapPos = new TH2D("hTwMapPos", "hTwMapPos", 220, -22., 22., 220, -22., 22.); // 2 mm/bin - 2 mm/bin
  hTwMapPos_Bar9 = new TH2D("hTwMapPos_Bar9", "hTwMapPos_Bar9", 220, -22., 22., 220, -22., 22.);

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
