#if !defined(__CINT__) || defined(__MAKECINT__)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <bitset>
#include <string>
#include <vector>
#include <map>

#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TObjString.h>
#include <TVector3.h>
#include <TH1.h>
#include <TROOT.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TDirectory.h>
#include <TLinearFitter.h>

#include "TAGcampaignManager.hxx"
#include "TAGactTreeReader.hxx"
#include "TAGrecoManager.hxx"
#include "TAGnameManager.hxx"
#include "TAGgeoTrafo.hxx"
#include "TAGionisMaterials.hxx"
#include "TAGrunInfo.hxx"
#include "TAGroot.hxx"
#include "TAGparGeo.hxx"

#include "TASTparGeo.hxx"
#include "TASTntuHit.hxx"

#include "TABMparGeo.hxx"
#include "TABMntuHit.hxx"
#include "TABMntuHit.hxx"
#include "TABMntuTrack.hxx"

#include "TAVTparGeo.hxx"
#include "TAVTntuHit.hxx"
#include "TAVThit.hxx"
#include "TAVTntuCluster.hxx"
#include "TAVTcluster.hxx"
#include "TAVTntuVertex.hxx"
#include "TAVTntuTrack.hxx"
#include "TAVTtrack.hxx"
#include "TAVTntuVertex.hxx"

#include "TAITparGeo.hxx"
#include "TAITntuCluster.hxx"

#include "TAMSDparGeo.hxx"
#include "TAMSDntuCluster.hxx"
#include "TAMSDcluster.hxx"
#include "TAMSDntuTrack.hxx"
#include "TAMSDntuPoint.hxx"
#include "TAMSDntuHit.hxx"
#include "TAMSDntuRaw.hxx"
#include "TAMSDhit.hxx"
#include "TAMSDtrack.hxx"
#include "TAMSDparameters.hxx"

#include "TATWparGeo.hxx"
#include "TATWntuHit.hxx"
#include "TATWntuPoint.hxx"
#include "TATWparameters.hxx"

#include "TACAparGeo.hxx"
#include "TACAntuHit.hxx"
#include "TACAntuCluster.hxx"

#include "TAMCntuEvent.hxx"
#include "TAMCntuPart.hxx"
#include "TAMCntuRegion.hxx"
#include "TAMCntuHit.hxx"

#include "TAGntuEvent.hxx"

#include "TAWDntuTrigger.hxx"

#include "TAGntuGlbTrack.hxx"


#endif

// static TAGrunInfo*   runinfo;  //runinfo
TAGrunInfo*   runinfo;  //runinfo

static TAGgeoTrafo*  geoTrafo; //geometry
static TAGparGeo* parGeo; //beam info
static TASTparGeo* stparGeo;
static TABMparGeo* bmparGeo;
static TAVTparGeo* vtparGeo;
static TAITparGeo* itparGeo;
static TAMSDparGeo* msdparGeo;
static TATWparGeo* twparGeo;
static TACAparGeo* caparGeo;

static TASTntuHit *stNtuHit;
static TAMCntuHit *stMcNtuHit;
static TABMntuHit*  bmNtuHit;
static TABMntuTrack*  bmNtuTrack;
static TAMCntuHit *bmMcNtuHit;
static TAVTntuVertex* vtxNtuVertex;
static TAVTntuCluster *vtxNtuCluster;
static TAVTntuTrack *vtxNtuTrack;
static TAMCntuHit *vtMcNtuHit;
static TAMSDntuCluster *msdNtuClus;
static TAMSDntuTrack *msdNtuTrack;
static TAMSDntuPoint *msdNtuPoint;
static TAMSDntuHit *msdNtuHit;
static TAMSDntuRaw *msdNtuRaw;
static TAMCntuHit *msdMcNtuHit;
static TAITntuCluster *itNtuClus;
static TAMCntuHit *itMcNtuHit;
static TATWntuHit *twNtuHit;
static TATWntuPoint *twNtuPoint;
static TAMCntuHit *twMcNtuHit;
static TACAntuHit *caNtuHit;
static TACAntuCluster *caNtuClus;
static TAMCntuHit *caMcNtuHit;

static TAMCntuPart *mcNtuPart;
static TAMCntuRegion *mcNtuRegion;
static TAGntuGlbTrack *glbntutrk;
static TAGntuEvent* tgNtuEvent;
static TAWDntuTrigger *wdNtuTrig;

static  Int_t vtxSensorsN;
static  Int_t itSensorsN;
static  Int_t msdSensorsN;
static  Int_t msdStationsN;

static int  IncludeTrk; // from runinfo
static int  IncludeReg; // from runinfo
static int  IncludeTOE; // from runinfo

//other global parameters taken from campaign manager
static int  IncludeMC;
static int  IncludeDI;
static int  IncludeSC;
static int  IncludeBM;
static int  IncludeTG;
static int  IncludeVT;
static int  IncludeIT;
static int  IncludeMSD;
static int  IncludeTW;
static int  IncludeCA;
static int  IncludeDAQ;
static int  IncludeWD;

Bool_t debug              = false;
Bool_t calibTw            = false;
Bool_t isTwScan           = true;
Bool_t readBinTwCalibFile = false;

enum{kCharges=8,kLayers=2,kBars=20};  //TW
enum{kModules=7, kCrysPerModule=9};  //Calorimeter
// enum{kCharges=8,kLayers=2,kCentralBars=3};
enum{kVTreg=2,kTWreg=4};
enum FlukaVar {kPrimaryID=0,kNeutronFlukaId=8};
enum TrigID {kTrigsN=4,kMBplusVeto=0,kVeto=1,kMB=40,kSTtrig=42};


typedef std::map<TrigID,TString> TMapTrig;
typedef std::map<TrigID,TH1D*> TMapTrigHisto;
TMapTrig mapTrig;
TMapTrigHisto mapTrigHisto;

const Int_t maxTGx = 10;
const Int_t maxTGy = 10;

const TString ElementName[kCharges+1] = {"H","He","Li","Be","B","C","N","O"};
const TString ElementName_with_n[kCharges+1] = {"n","H","He","Li","Be","B","C","N","O"};

//const Int_t  CentralBarsID[kCentralBars] = {8,9,10}; // for both the layers 

typedef std::map<Int_t,std::vector<Int_t> > TMapMatch;
typedef std::vector<std::pair<Int_t,Int_t> > TVecPair;

// TH1F*           fpHisSeedMap;    ///< seed map
// TH1F*           fpHisStripMap;   ///< strip map
// TH1F*           fpHisSeedMap[MaxPlane];    ///< seed map
// TH1F*           fpHisStripMap[MaxPlane];   ///< strip map

TH1D *Eloss;
TH1D *Eloss_noCuts;
TH1D *Tof;
TH1D *Tof_noCuts;
TH2D *dE_vs_tof;
TH2D *dE_vs_tof_Z1;
TH2D *dE_vs_tof_Z2;
TH2D *beta_vs_dE;
TH1D *Eloss_perBar[kLayers][nBarsPerLayer];
TH1D *Eloss_perBar_noCuts[kLayers][nBarsPerLayer];
TH1D *My_eloss_noCuts[kLayers][nBarsPerLayer];
TH1D *Charge_perBar[kLayers][nBarsPerLayer];
TH1D *Charge_perBar_noCuts[kLayers][nBarsPerLayer];
TH1D *hToF[kLayers][nBarsPerLayer];
TH1D *hToF_noCuts[kLayers][nBarsPerLayer];
TH1D *PosX;  //bar hit position in layerX
TH1D *PosY;  //bar hit position in layerY
TH1D *Bar_ID_X;  //bar ID for a given layer
TH1D *Bar_ID_Y;
TH1D *hHits_X;
TH1D *h_nValidHits_X;
TH1D *h_nValidHits_Y;
TH1D *hHits_Y;
TH2D *hTwMapPos;
TH1D *Z_reco;


void  InitializeContainers();
void  BookHistograms(TDirectory* DirChargeTimeLayerX, TDirectory* DirChargeTimeLayerY, 
                    TDirectory* DirToFLayerX, TDirectory* DirToFLayerY);
void  GetFOOTgeo(TAGcampaignManager* camp_manager, Int_t run_number);
void  GetRunAndGeoInfo(TAGcampaignManager* campManager, Int_t runNumber);
void  SetTreeBranchAddress(TAGactTreeReader *treeReader);
void  ProjectTracksOnTw(int Z, TVector3 init_pos, TVector3 init_p);
void  LoopOverMCtracks(Int_t Emin, Int_t Emax, Bool_t isnotrig);
void  AdjustHistoRange(TH1D *Histo);
void SetTitleAndLabels(TObject* obj, const char* title, const char* xLabel, const char* yLabel);
Double_t CalculateZ(Double_t dE, Double_t beta);
std::map<Int_t, std::map<Int_t, Double_t>> extractBarData();
std::map<Int_t, std::map<Int_t, Double_t>> extractTofData(Int_t energy);

Bool_t IsVTregion(int reg);

ofstream fileBinTwCalib;

inline Int_t GetZbeam() {return parGeo->GetBeamPar().AtomicNumber;}   
// TW center in global ref frame
inline TVector3 GetTwCenter() {return geoTrafo->GetTWCenter();}   
// TW theta angle geometrical acceptance
inline Double_t GetMaxAngle() {return TMath::ATan(((nBarsPerLayer*twparGeo->GetBarWidth())/2-TMath::Abs(GetTwCenter().y())-maxTGy)/GetTwCenter().z()); } // rad   
