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
#include <TSpectrum.h>


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

TH2D *dE_vs_tof[kLayers];
TH2D *dE_vs_tof_perBar[kLayers][nBarsPerLayer];
TH1D *ChargeA_perBar[kLayers][nBarsPerLayer];
TH1D *ChargeB_perBar[kLayers][nBarsPerLayer];
TH1D *Charge_perBar[kLayers][nBarsPerLayer];  //sqrt(QA*QB)
TH1D *TimeA_perBar[kLayers][nBarsPerLayer];
TH1D *TimeB_perBar[kLayers][nBarsPerLayer];
TH1D *Time_perBar[kLayers][nBarsPerLayer];  //0.5*(TA + TB)
TH1D *Charge_Calo_total;  //charge in all calo
TH1D *Charge_Calo_crystal[kModules * kCrysPerModule];  //charge per crystal id

// TH2D *dE_vs_tof_perBar[kLayers][kBars];
TH1D *heloss_all;
// TH1D *heloss[kTrigsN];;
TH2D *hTwPos[kLayers];
TH2D *hTwMapPos;
TH2D *hCalMapPos[kModules];  //2D histogram of x, y positions in the calorimeter
TH2D *hTwMapPos_TWpntBin;
TH2D *hTwMapPos_TWpnt;
TH2D *hTwMapPos_TWpnt_Z[kCharges];
TH1D *hResX[kCharges];
TH1D *hResY[kCharges];
TH2D *hTwMapPos_1Cross;
TH1D *hResX_1Cross;
TH1D *hResY_1Cross;


void  InitializeContainers();
void  BookHistograms();
void  GetFOOTgeo(TAGcampaignManager* camp_manager, Int_t run_number);
void  GetRunAndGeoInfo( TAGcampaignManager* campManager, Int_t runNumber);
void  SetTreeBranchAddress(TAGactTreeReader *treeReader);
void  ProjectTracksOnTw(int Z, TVector3 init_pos, TVector3 init_p);
void  LoopOverMCtracks(Int_t Emin, Int_t Emax, Bool_t isnotrig);
void  FillBinaryForTwCalib(Int_t eventN, Int_t nentries, Int_t runNumber);
void  FillBinaryForTwCalib_MC(Int_t eventN, Int_t nentries, Int_t runNumber);

Bool_t IsVTregion(int reg);

ofstream fileBinTwCalib;

inline Int_t GetZbeam() {return parGeo->GetBeamPar().AtomicNumber;}   
// TW center in global ref frame
inline TVector3 GetTwCenter() {return geoTrafo->GetTWCenter();}   
// TW theta angle geometrical acceptance
inline Double_t GetMaxAngle() {return TMath::ATan(((nBarsPerLayer*twparGeo->GetBarWidth())/2-TMath::Abs(GetTwCenter().y())-maxTGy)/GetTwCenter().z()); } // rad   


