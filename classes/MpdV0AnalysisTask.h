#ifndef MPDV0ANALYSISTASK_H
#define MPDV0ANALYSISTASK_H

#include "MpdAnalysisTask.h"
#include "FairMCEventHeader.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "MpdHelix.h"
#include "MpdParticle.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanTrack.h"
#include "MpdV0track.h"
//#include "MpdTpcSectorGeo.h"
#include "MpdPid.h" 
#include "MpdKalmanFilter.h"
#include "MpdZdcDigi.h"

class MpdTpcKalmanTrack;
class MpdV0AnalysisTask : public MpdAnalysisTask {

public:
   MpdV0AnalysisTask() {}
   MpdV0AnalysisTask(const char *name, const char *outputName = "taskName");
   virtual ~MpdV0AnalysisTask() {}


   virtual void UserInit();
   virtual void ProcessEvent(MpdAnalysisEvent &event);
   virtual void Finish();

   void setOutFile(std::string filename = "histos.root") { mOutFile = filename; }

protected:

   bool selectEvent(MpdAnalysisEvent &event); // vertex 
   bool selectTrackpr(MpdTrack *track); //choose nice track
   bool selectTrackpi(MpdTrack *track); //choose nice track
   //functions
   Int_t GetMultiplicity(MpdAnalysisEvent &event, Double_t ptCut, Double_t etaCut, Double_t dcaCut, Int_t nhitsCut, Bool_t isMC);
   Double_t GetProjection(TLorentzVector La, TLorentzVector pr, TVector3 Proy);
   Double_t GetDeltaPhi(TLorentzVector La, TLorentzVector pr, TVector3 Epv);
   Int_t GetCentClass(Int_t multiplicityE);
   bool BuildLambda(MpdTrack *trp, MpdTpcKalmanTrack *ktrp, MpdTrack *trn, MpdTpcKalmanTrack *ktrn, MpdV0track &v0can);
   Double_t EventPlaneZDC(MpdAnalysisEvent &event);
   MpdHelix MakeHelix(const MpdKalmanTrack &tr) const ;
   Double_t Getctau(TLorentzVector La, Double_t distance);   

private:

   // event properties
   float cen;
   int   anaBin;

   bool     isInitialized = false;
   int      mCenBin       = 0;
   int      mZvtxBin      = 0;
   int      mRPBin        = 0;
   int      mixBin;

   TVector3 mPrimaryVertex;
   static constexpr short nMixEventZ    = 10; //(V) number of bins in z direction
   static constexpr short nMixEventCent = 10; //(V) number of bins of centrality
   static constexpr short nMixEventRP   = 1;  //(V) number of bins of Reaction Plane orientation

   MpdKalmanFilter* mKF          = nullptr ;
   MpdPid *mPID                  = nullptr ;
   MpdKalmanHit mKHit;


  const Double_t cut_pt = 0.15; // default: 0.15 GeV/c
  const Double_t cut_eta = 0.5; // default: 0.5
  const Int_t cut_nhits = 16;   // default: 16
  const Double_t dca_cut = 0.5; // default: 0.5 cm


   std::string mOutFile = "histos.root" ;

   FairMCEventHeader *MCHeader;
   MpdEvent *fTDstEvent;
   TClonesArray *fTMCTracks;
   TClonesArray *fTVtx;
   TClonesArray  *fTpcKlTracks;
   TClonesArray  *fTMpdGlobalTracks = nullptr;
   TClonesArray  *fTZDCDigits;

  
  //objects
   std::vector<MpdV0track> mV0;  // v0 array? 
   vector<MpdParticle*> mPartV0;			 //

   //Histograms
   
   // General QA  ==>>> taken from example class
   //
   TH1F *mhEvents       = nullptr;
   TH1F *mhVertex       = nullptr;
   TH1F *mhCentrality   = nullptr;
   TH1F *mhMultiplicity = nullptr;
   //
   TH2F *mhdEdx;
   TH2F *mhdEdxAss;
   TH2F *mhdEdxvsmass2;
   TH2F *mhdEdxvsmass2Ass;
   TH2F *mhdEdxpi;
   TH2F *mhdEdxpiAss;
   TH2F *mhdEdxvsmass2pi;
   TH2F *mhdEdxvsmass2piAss;

   TProfile *mhdEdxvspAssProf;
   TProfile *mhdEdxvsppiAssProf;
   //
   //
   TH1F *fhistPtMC;

   TH1F *fhRefMultMC;
   TH2F *fhBvsRefMultMC;

   TProfile *fhCentvsB;
   TProfile *fhCentvsBbis;

   TH2F *fhArmPodMC;
   TH2F *fhArmPodMCcut;

   TH1F *hMassV0MC[2]={0,0};
   TH1F *hMassV0MCcut[2]={0,0};
   TH1F *fhProjection[2][3][10];
   TH1F *fhProjectioncut[2][3][10];
   TH1F *fhDeltaPhi[2][10];
   TH1F *fhDeltaPhicut[2][10];
   TH1F *fhDeltaPhiRec[2][10];
   TH1F *fhDeltaPhiAss[2][10];
   TH1F *fhProjLocRec[2][10];
   TH1F *fhProjLocAss[2][10];

   TH1F *fhDeltaPsiRec[2][10];
   TH1F *fhDeltaPsiAss[2][10];
   TH1F *fhProjPsiRec[2][10];
   TH1F *fhProjPsiAss[2][10];

   TH1F *fhPolarization[2][3][10];
   TH1F *fhPolarizationcut[2][3][10];

   TH1F *fhtest;
   TH1F *fhtest2;

   TH1F *fhchi2La;
   TH1F *fhptLa;
   TH1F *fhrConv;
   TH1F *fhAsym1;
   TH1F *fhAsym2;

   TH1F *fhchi2LaAss;
   TH1F *fhptLaAss;
   TH1F *fhrConvAss;
   TH1F *fhAsymAss1;
   TH1F *fhAsymAss2;

   //Reconstructed tracks
   TH1F *fhRefMult;
   TH2F *fhBvsRefMult;

   TH2F *fhistPsiRecvsPsiMC;
   TProfile *fhResolutionEP;
   TH2F *fhistPsiRecvsPsiMCbis;
   TProfile *fhResolutionEPbis;



   TH1F *hInvMass;
   TH2F *fhArmPod;
   TH2F *fhDecayctau;
   TH2F *fhCosInvMass;
   TH2F *fhDCAv0InvMass;
   TH2F *fhDCAposInvMass;
   TH2F *fhDCAnegInvMass;
   TH1F *fhDeltaPhiRecLa;
   TH1F *fhProjLocRecLa;

   //Associated tracks
   TH1F *hInvMassAss;
   TH2F *fhArmPodAss;
   TH2F *fhDecayctauAss;
   TH2F *fhCosInvMassAss;
   TH2F *fhDCAv0InvMassAss;
   TH2F *fhDCAposInvMassAss;
   TH2F *fhDCAnegInvMassAss;
   TH1F *fhDeltaPhiAssLa;
   TH1F *fhProjLocAssLa;



   ClassDef(MpdV0AnalysisTask,1);
};
#endif














