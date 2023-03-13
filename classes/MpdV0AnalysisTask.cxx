#include "MpdV0AnalysisTask.h"

#include "MpdMCTrack.h"
#include "MpdVertex.h"
#include "MpdParticle.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdKalmanTrack.h"



#include <iostream>
#include <fstream>
using std::cout;
using std::endl;

ClassImp(MpdV0AnalysisTask);

MpdV0AnalysisTask::MpdV0AnalysisTask(const char *name, const char *outputName) : MpdAnalysisTask(name,outputName)
//	MCHeader(0),
//	fTDstEvent(0),
//	fTMCTracks(0)
{
   fOutputList = nullptr;
}
//_____________________________________________________________________________________

//_____________________________________________________________________________________
void MpdV0AnalysisTask::UserInit()
{
	const Char_t *Hyperon[2]={"Lambda","LambdaBar"};
	const Char_t *Ref[3]={"P","n","L"};
	const Char_t *Eje[3]={"x","y","z"};

	fOutputList = new TList();
	fOutputList->SetOwner(kTRUE);

	TH1::AddDirectory(kFALSE);

Double_t n2pi = 2*TMath::Pi();


	   // General QA
   mhEvents = new TH1F("hEvents", "Number of events", 10, 0., 10.);
   fOutputList->Add(mhEvents);
   mhVertex = new TH1F("hVertex", "Event vertex distribution", 100, -200., 200.);
   fOutputList->Add(mhVertex);
   mhCentrality = new TH1F("hCentrality", "Centrality distribution", 100, 0., 100.);
   fOutputList->Add(mhCentrality);
   mhMultiplicity = new TH1F("hMultiplicity", "Multiplicity distribution", 2000, -0.5, 1999.5);
   fOutputList->Add(mhMultiplicity);


	fhistPtMC = new TH1F("fhistPtMC","MC p_{T} distribution; p_{T}(GeV/c); dN/dp_{T}",400,0,10);
	fOutputList->Add(fhistPtMC);

	fhRefMultMC =new TH1F("hRefMultMC","hRefMultMC;N_{ch}; dN_{Ev}/dN_{ch}",2500,0,2500);
	fOutputList->Add(fhRefMultMC);

	fhBvsRefMultMC = new TH2F("hBvsRefMultMC","hBvsRefMultMC;N_{ch};b(fm)",2500,0,2500,200,0.,20.);
	fOutputList->Add(fhBvsRefMultMC);

	fhCentvsB = new TProfile("hCentvsB","hCentvsB",10,0,100,0,20);
	fOutputList->Add(fhCentvsB);

        fhArmPodMC = new TH2F("hArmPodMC","hArmPodMC; #alpha; p_{T}",100,-1.0,1.0,100,0,0.25);
	fOutputList->Add(fhArmPodMC);
        fhArmPodMCcut = new TH2F("hArmPodMCcut","hArmPodMC cut; #alpha; p_{T}",100,-1.0,1.0,100,0,0.25);
	fOutputList->Add(fhArmPodMCcut);

 	for(Int_t ip=0;ip<2;++ip){
 	 hMassV0MC[ip]=new TH1F(Form("hMassMC_%s",Hyperon[ip]),Form("hMassMC_%s",Hyperon[ip]),100,1.08,1.12);
	 fOutputList->Add(hMassV0MC[ip]);
 	 hMassV0MCcut[ip]=new TH1F(Form("hMassMCcut_%s",Hyperon[ip]),Form("hMassMCcut_%s",Hyperon[ip]),100,1.08,1.12);
	 fOutputList->Add(hMassV0MCcut[ip]);

	for(Int_t k=0;k<10;++k){
 	fhDeltaPhi[ip][k]=new TH1F(Form("hdeltaphi%s_%d",Hyperon[ip],k),Form("hdeltaphi%s_%d;dN/d #Delta #phi; #Delta #phi",Hyperon[ip],k),100,0,2*TMath::Pi());
	fOutputList->Add(fhDeltaPhi[ip][k]);
 	fhDeltaPhicut[ip][k]=new TH1F(Form("hdeltaphicut%s_%d",Hyperon[ip],k),Form("hdeltaphicut%s_%d;dN/d #Delta #phi; #Delta #phi",Hyperon[ip],k),100,0,2*TMath::Pi());
	fOutputList->Add(fhDeltaPhicut[ip][k]);
	}

 	for(Int_t j=0;j<3;++j){
	for(Int_t k=0;k<10;++k){
	fhProjection[ip][j][k]=new TH1F(Form("hproj%s_%s_%d",Hyperon[ip],Ref[j],k),Form("hproj%s_%s_%d;dN/d cos(#sigma); cos(#sigma)",Hyperon[ip],Ref[j],k),100,-1,1);
	 fOutputList->Add(fhProjection[ip][j][k]);
	fhProjectioncut[ip][j][k]=new TH1F(Form("hprojcut%s_%s_%d",Hyperon[ip],Ref[j],k),Form("hprojcut%s_%s_%d;dN/d cos(#sigma); cos(#sigma)",Hyperon[ip],Ref[j],k),100,-1,1);
	 fOutputList->Add(fhProjectioncut[ip][j][k]);
	}
	}

 	for(Int_t j=0;j<3;++j){
	for(Int_t k=0;k<10;++k){
	fhPolarization[ip][j][k]=new TH1F(Form("hpol%s_%s_%d",Hyperon[ip],Eje[j],k),Form("hpol%s_%s_%d;dN/dP; P_%s",Hyperon[ip],Eje[j],k,Eje[j]),100,-1,1);
	 fOutputList->Add(fhPolarization[ip][j][k]);
	fhPolarizationcut[ip][j][k]=new TH1F(Form("hpolcut%s_%s_%d",Hyperon[ip],Eje[j],k),Form("hpolcut%s_%s_%d;dN/dP; P_%s",Hyperon[ip],Eje[j],k,Eje[j]),100,-1,1);
	 fOutputList->Add(fhPolarizationcut[ip][j][k]);
	}
	}

 	}

	fhtest = new TH1F("fhtest","MC or rec variable distribution; variable; dN/dvariable",1000,-5,200);
	fOutputList->Add(fhtest);

	fhtest2 = new TH1F("fhtest2","MC variable distribution; variable; dN/dvariable",1000,-5,200);
	fOutputList->Add(fhtest2);

	//Reconstructed tracks

        fhchi2La = new TH1F("hchi2La","hchi2La",100,0,30);
	fOutputList->Add(fhchi2La);
   	fhptLa = new TH1F("hptLa","hPtLa",100,-1,9);
	fOutputList->Add(fhptLa);
	fhrConv = new TH1F("hrConv","hrConv",100,0,150);
	fOutputList->Add(fhrConv);
	fhAsym1 = new TH1F("h1asym1y2","h1asym1y2",100,0,2);
	fOutputList->Add(fhAsym1);
	fhAsym2 = new TH1F("h2asym1y2","h2asym1y2",100,0,2);
	fOutputList->Add(fhAsym2);

	//

        fhistPsiRecvsPsiMC = new TH2F("histEPrecvsEPmc","EP_{rec} vs EP_{MC},#Psi_{rec};#Psi_{MC}",100,0,n2pi,100,0,n2pi);
	fOutputList->Add(fhistPsiRecvsPsiMC);

	fhResolutionEP = new TProfile("hResolutionEP","hResolutionEP;centrality;Resolution",10,0,100,-2,2);
	fOutputList->Add(fhResolutionEP);
	//
	mhdEdx = new TH2F("hdEdx", "dEdx", 200, 0, 5., 1000, 0, 1000);
        fOutputList->Add(mhdEdx);

	mhdEdxvsmass2 = new TH2F("hdEdxvsmass2", "dEdx vs mass2", 200, 0, 3., 1000, 0, 1000);
        fOutputList->Add(mhdEdxvsmass2);

	mhdEdxpi = new TH2F("hdEdxpi", "dEdx", 200, 0, 5., 1000, 0, 1000);
        fOutputList->Add(mhdEdxpi);

	mhdEdxvsmass2pi = new TH2F("hdEdxvsmass2pi", "dEdx vs mass2", 500, 0, 1., 1000, 0, 1000);
        fOutputList->Add(mhdEdxvsmass2pi);

	//
	fhRefMult =new TH1F("hRefMult","hRefMult;N_{ch}; dN_{Ev}/dN_{ch}",2500,0,2500);
	fOutputList->Add(fhRefMult);

	fhBvsRefMult = new TH2F("hBvsRefMult","hBvsRefMult;N_{ch};b(fm)",2500,0,2500,200,0.,20.);
	fOutputList->Add(fhBvsRefMult);


 	hInvMass=new TH1F("hInvMass","hInvMass",100,1.08,1.2);
	fOutputList->Add(hInvMass);

        fhArmPod = new TH2F("hArmPod","hArmPod; #alpha; p_{T}",100,-1.0,1.0,100,0,0.25);
	fOutputList->Add(fhArmPod);

	fhDecayctau = new TH2F("hDecayctau","hDecayctau",100,-1.0,15,100,-1,15);
	fOutputList->Add(fhDecayctau);


	fhCosInvMass = new TH2F("hCosInvMass","hCosInvMass",100,1.08,1.17,100,-1.0,1.0);
	fOutputList->Add(fhCosInvMass);

	fhDCAv0InvMass = new TH2F("hDCAv0InvMass","hDCAv0InvMass",100,1.08,1.17,100,0,5);
	fOutputList->Add(fhDCAv0InvMass);

	fhDCAposInvMass = new TH2F("hDCAposInvMass","hDCAposInvMass",100,1.08,1.17,100,0,5);
	fOutputList->Add(fhDCAposInvMass);

	fhDCAnegInvMass = new TH2F("hDCAnegInvMass","hDCAnegInvMass",100,1.08,1.17,100,0,5);
	fOutputList->Add(fhDCAnegInvMass);

	fhDeltaPhiRecLa=new TH1F("hdeltaphiRecLa","hdeltaphiRecLa;dN/d #Delta #phi; #Delta #phi",100,0,2*TMath::Pi());
	fOutputList->Add(fhDeltaPhiRecLa);


	fhProjLocRecLa=new TH1F("hprojLocRecLa","hprojLocRecLa;dN/d cos(#sigma); cos(#sigma)",100,-1,1);
	 fOutputList->Add(fhProjLocRecLa);

 	for(Int_t ip=0;ip<2;++ip){
	for(Int_t k=0;k<10;++k){

	fhDeltaPhiRec[ip][k]=new TH1F(Form("hdeltaphiRec%s_%d",Hyperon[ip],k),Form("hdeltaphiRec%s_%d;dN/d #Delta #phi; #Delta #phi",Hyperon[ip],k),100,0,2*TMath::Pi());
	fOutputList->Add(fhDeltaPhiRec[ip][k]);

	fhProjLocRec[ip][k]=new TH1F(Form("hprojLocRec%s_%d",Hyperon[ip],k),Form("hprojLocRec%s_%d;dN/d cos(#sigma); cos(#sigma)",Hyperon[ip],k),100,-1,1);
	 fOutputList->Add(fhProjLocRec[ip][k]);

	fhDeltaPsiRec[ip][k]=new TH1F(Form("hdeltapsiRec%s_%d",Hyperon[ip],k),Form("hdeltapsiRec%s_%d;dN/d #Delta #psi; #Delta #psi",Hyperon[ip],k),100,0,2*TMath::Pi());
	fOutputList->Add(fhDeltaPsiRec[ip][k]);

	fhProjPsiRec[ip][k]=new TH1F(Form("hprojPsiRec%s_%d",Hyperon[ip],k),Form("hprojPsiRec%s_%d;dN/d cos(#sigma); cos(#sigma)",Hyperon[ip],k),100,-1,1);
	 fOutputList->Add(fhProjPsiRec[ip][k]);

	}
	}



	//Associated tracks
	//
	
        fhchi2LaAss = new TH1F("hchi2LaAss","hchi2LaAss",100,0,30);
	fOutputList->Add(fhchi2LaAss);
   	fhptLaAss = new TH1F("hptLaAss","hPtLaAss",100,-1,9);
	fOutputList->Add(fhptLaAss);
	fhrConvAss = new TH1F("hrConvAss","hrConvAss",100,0,150);
	fOutputList->Add(fhrConvAss);
	fhAsymAss1 = new TH1F("hasym1y2Ass1","hasym1y2Ass1",100,0,2);
	fOutputList->Add(fhAsymAss1);
	fhAsymAss2 = new TH1F("hasym1y2Ass2","hasym1y2Ass2",100,0,2);
	fOutputList->Add(fhAsymAss2);
	//
		mhdEdxAss = new TH2F("hdEdxAss", "dEdxAss", 200, 0, 5., 1000, 0, 1000);
        fOutputList->Add(mhdEdxAss);

		mhdEdxvspAssProf = new TProfile("hdEdxvspAssProf", "dEdxvspAss", 200, 0, 5., 0, 1000);
        fOutputList->Add(mhdEdxvspAssProf);


	mhdEdxvsmass2Ass = new TH2F("hdEdxvsmass2Ass", "dEdx vs mass2", 200, 0, 3., 1000, 0, 1000);
        fOutputList->Add(mhdEdxvsmass2Ass);

	mhdEdxpiAss = new TH2F("hdEdxpiAss", "dEdx", 200, 0, 5., 1000, 0, 1000);
        fOutputList->Add(mhdEdxpiAss);

	mhdEdxvsppiAssProf = new TProfile("hdEdxvsppiAssProf", "dEdx", 200, 0, 5., 0, 1000);
        fOutputList->Add(mhdEdxvsppiAssProf);

	mhdEdxvsmass2piAss = new TH2F("hdEdxvsmass2piAss", "dEdx vs mass2", 500, 0, 1., 1000, 0, 1000);
        fOutputList->Add(mhdEdxvsmass2piAss);

	//
	
 	hInvMassAss=new TH1F("hInvMassAss","hInvMass",100,1.08,1.2);
	fOutputList->Add(hInvMassAss);

        fhArmPodAss = new TH2F("hArmPodAss","hArmPodAss; #alpha; p_{T}",100,-1.0,1.0,100,0,0.25);
	fOutputList->Add(fhArmPodAss);

	fhDecayctauAss = new TH2F("hDecayctauAss","hDecayctauAss;d_{measured};d_{mass hypothesis}",100,-1.0,15,100,-1,15);
	fOutputList->Add(fhDecayctauAss);


	fhCosInvMassAss = new TH2F("hCosInvMassAss","hCosInvMassAss",100,1.08,1.17,100,-1.0,1.0);
	fOutputList->Add(fhCosInvMassAss);

	fhDCAv0InvMassAss = new TH2F("hDCAv0InvMassAss","hDCAv0InvMassAss",100,1.08,1.17,100,0,5);
	fOutputList->Add(fhDCAv0InvMassAss);

	fhDCAposInvMassAss = new TH2F("hDCAposInvMassAss","hDCAposInvMassAss",100,1.08,1.17,100,0,5);
	fOutputList->Add(fhDCAposInvMassAss);

	fhDCAnegInvMassAss = new TH2F("hDCAnegInvMassAss","hDCAnegInvMassAss",100,1.08,1.17,100,0,5);
	fOutputList->Add(fhDCAnegInvMassAss);
 
	fhDeltaPhiAssLa=new TH1F("hdeltaphiAssLa","hdeltaphiAssLa;dN/d #Delta #phi; #Delta #phi",100,0,2*TMath::Pi());
	fOutputList->Add(fhDeltaPhiAssLa);

	fhProjLocAssLa=new TH1F("hprojLocAssLa","hprojLocAssLa;dN/d cos(#sigma); cos(#sigma)",100,-1,1);
	 fOutputList->Add(fhProjLocAssLa);


 	for(Int_t ip=0;ip<2;++ip){
	for(Int_t k=0;k<10;++k){
 	fhDeltaPhiAss[ip][k]=new TH1F(Form("hdeltaphiAss%s_%d",Hyperon[ip],k),Form("hdeltaphiAss%s_%d;dN/d #Delta #phi; #Delta #phi",Hyperon[ip],k),100,0,2*TMath::Pi());
	fOutputList->Add(fhDeltaPhiAss[ip][k]);

	fhProjLocAss[ip][k]=new TH1F(Form("hprojLocAss%s_%d",Hyperon[ip],k),Form("hprojLocAss%s_%d;dN/d cos(#sigma); cos(#sigma)",Hyperon[ip],k),100,-1,1);
	 fOutputList->Add(fhProjLocAss[ip][k]);

 	fhDeltaPsiAss[ip][k]=new TH1F(Form("hdeltapsiAss%s_%d",Hyperon[ip],k),Form("hdeltapsiAss%s_%d;dN/d #Delta #psi; #Delta #psi",Hyperon[ip],k),100,0,2*TMath::Pi());
	fOutputList->Add(fhDeltaPsiAss[ip][k]);

	fhProjPsiAss[ip][k]=new TH1F(Form("hprojPsiAss%s_%d",Hyperon[ip],k),Form("hprojPsiAss%s_%d;dN/d cos(#sigma); cos(#sigma)",Hyperon[ip],k),100,-1,1);
	 fOutputList->Add(fhProjPsiAss[ip][k]);

	}
	}


}
//____________________________________________________________________________________
void MpdV0AnalysisTask::ProcessEvent(MpdAnalysisEvent &event)
{
    if (!isInitialized) {
      mKF = MpdKalmanFilter::Instance();
      mKHit.SetType(MpdKalmanHit::kFixedR);
      mPID          = new MpdPid(2.0, 2.0, 9.2, 1.0, "NSIG", "CFHM", "pikapr"); //https://git.jinr.ru/nica/mpdroot/-/tree/dev/core/mpdPid#MpdPid_2_2
			//MpdPid(Double_t sigM, Double_t sigE, Double_t E, Double_t C, TString generator, TString tracking, TString nSigPart);
      isInitialized = true;
   }
 
   if(!selectEvent(event)){
   return;
   }

 Double_t impb=-99.0;
 Double_t PsiEP=-99.0;

 Double_t PsizdcEP=-99.0;


 MCHeader = event.fMCEventHeader;
 impb = MCHeader->GetB();
 PsiEP = MCHeader->GetRotZ();


 fTMCTracks = event.fMCTrack; //branches name defined in MpdAnalysisManager
 Int_t nmctracks=fTMCTracks->GetEntriesFast();


 Int_t multiplicityMC =  GetMultiplicity(event,cut_pt,cut_eta,dca_cut,cut_nhits,kTRUE);
 fhRefMultMC->Fill(multiplicityMC);
 fhBvsRefMultMC->Fill(multiplicityMC,impb);
 Int_t multiplicity =  GetMultiplicity(event,cut_pt,cut_eta,dca_cut,cut_nhits,kFALSE);
 fhRefMult->Fill(multiplicity);
 fhBvsRefMult->Fill(multiplicity,impb);

 Int_t index_cent= -99;
 index_cent = GetCentClass(multiplicity);

 Double_t centrality = 10*index_cent + 5;

 fhCentvsB->Fill(centrality,impb);

 PsizdcEP = EventPlaneZDC(event);
 
 Double_t reszdc = TMath::Cos(PsizdcEP - PsiEP);

 fhResolutionEP->Fill(centrality,reszdc);
 fhistPsiRecvsPsiMC->Fill(PsizdcEP,PsiEP);



 //cout << "N of MC tracks = " << nmctracks << endl;

 for (Int_t i = 0; i < nmctracks; i++)
  {
    MpdMCTrack *MCtrack = (MpdMCTrack*) fTMCTracks->UncheckedAt(i);
    Double_t ptmc=MCtrack->GetPt();
    fhistPtMC->Fill(ptmc);

    if(MCtrack->GetMotherId()==-1)continue; // remove primaries and keep only secondaries
    Int_t ppdg = MCtrack->GetPdgCode();
          if((ppdg==2212) || (ppdg == 211)){
      auto mctrackmom = (MpdMCTrack*) fTMCTracks->UncheckedAt(MCtrack->GetMotherId());

       Int_t mompdg = mctrackmom->GetPdgCode();
       if(TMath::Abs(mompdg)!=3122)continue;
           for(Int_t iTrn=0; iTrn<nmctracks; iTrn++){
      if(iTrn == i)continue;
      auto mctrackn = (MpdMCTrack*) fTMCTracks->UncheckedAt(iTrn);

    if (mctrackn->GetMotherId()!= MCtrack->GetMotherId())continue;
      Int_t npdg=mctrackn->GetPdgCode();

            if((npdg==-2212) || (npdg == -211)){
//    cout << "Mother ID = " << mctrackn->GetMotherId() << ",      Mother pdg = " << mompdg << ",       pdg code pos = "<< ppdg<< ",       pdg code neg = "<< npdg<<  endl;
      // positive daughter 
       TVector3 P(MCtrack->GetPx(),MCtrack->GetPy(),MCtrack->GetPz());
       Double_t masspos=MCtrack->GetMass();
       Double_t Enpos=MCtrack->GetEnergy();
       TLorentzVector PosL(P,Enpos);
     // negative daughter 
       TVector3 Pn(mctrackn->GetPx(),mctrackn->GetPy(),mctrackn->GetPz());
       Double_t massneg=mctrackn->GetMass();
       Double_t Enneg=mctrackn->GetEnergy();


       Double_t etamc=0.5*TMath::Log((P.Mag() + MCtrack->GetPz())/(P.Mag() - MCtrack->GetPz()+1.e-13));
       Double_t etamcn=0.5*TMath::Log((Pn.Mag() + mctrackn->GetPz())/(Pn.Mag() - mctrackn->GetPz()+1.e-13));
       Double_t ptmcn=mctrackn->GetPt();

     // Lambda Reco  
       Double_t laInvmass=TMath::Sqrt(masspos*masspos + massneg*massneg +2*(Enpos*Enneg - P.Dot(Pn)));
      if(mompdg ==  3122)hMassV0MC[0]->Fill(laInvmass);
      if(mompdg == -3122)hMassV0MC[1]->Fill(laInvmass);


      TVector3 PLa=Pn+P;
      Double_t EnLa=Enpos+Enneg;
      TLorentzVector HypL(PLa,EnLa);

      //ArmPod variables
 
      Float_t       lQlNeg = Pn.Dot(PLa)/PLa.Mag();
      Float_t       lQlPos = P.Dot(PLa)/PLa.Mag();
      Float_t       alpha = (lQlPos - lQlNeg)/(lQlPos + lQlNeg);
      Float_t       ptarm = Pn.Perp(PLa);
      
      fhArmPodMC->Fill(alpha, ptarm);
//    Double_t ctau=Getctau(HypL);
//    fhtest2->Fill(ctau);

      // Angular momentum
      TVector3 mL(TMath::Sin(PsiEP),-TMath::Cos(PsiEP),0); 
      TVector3 PolL(mctrackmom->GetPolar(0),mctrackmom->GetPolar(1),mctrackmom->GetPolar(2));

      Double_t ptLa=TMath::Sqrt(TMath::Power(PLa.X(),2) + TMath::Power(PLa.Y(),2));
      TVector3 nL(-PLa.Y()/ptLa,PLa.X()/ptLa,0);

      Double_t adPL=GetProjection(HypL,PosL,PolL);
      Double_t adnL=GetProjection(HypL,PosL,nL);
      Double_t admL=GetProjection(HypL,PosL,mL);

      //if(PolL.Mag()!=0.)
      //fhtest->Fill(GetDeltaPhi(HypL,PosL,mL));

      Double_t weight_pol = mctrackmom->GetWeight(); // if PHSD data
      //Double_t weight_pol = 0.4; // if UrQMD data

	      if(mompdg == 3122){
	      	fhPolarization[0][0][index_cent]->Fill(weight_pol*PolL.X());
	      	fhPolarization[0][1][index_cent]->Fill(weight_pol*PolL.Y());
	      	fhPolarization[0][2][index_cent]->Fill(weight_pol*PolL.Z());
	      	if(PolL.Mag()!=0.)fhProjection[0][0][index_cent]->Fill(adPL);
	      	fhProjection[0][1][index_cent]->Fill(adnL);
	      	fhProjection[0][2][index_cent]->Fill(admL);
		fhDeltaPhi[0][index_cent]->Fill(GetDeltaPhi(HypL,PosL,mL));
		if( 
				(etamc < cut_eta) && (etamcn < cut_eta) && (ptmc > cut_pt) && (ptmcn > cut_pt)
				){
	      		fhPolarizationcut[0][0][index_cent]->Fill(weight_pol*PolL.X());
	      		fhPolarizationcut[0][1][index_cent]->Fill(weight_pol*PolL.Y());
	      		fhPolarizationcut[0][2][index_cent]->Fill(weight_pol*PolL.Z());
	      		if(PolL.Mag()!=0.)fhProjectioncut[0][0][index_cent]->Fill(adPL);
	      		fhProjectioncut[0][1][index_cent]->Fill(adnL);
	      		fhProjectioncut[0][2][index_cent]->Fill(admL);
			fhDeltaPhicut[0][index_cent]->Fill(GetDeltaPhi(HypL,PosL,mL));
                        hMassV0MCcut[0]->Fill(laInvmass);
                        fhArmPodMCcut->Fill(alpha, ptarm);
		}
	      }
	      if(mompdg == -3122){
	      	fhPolarization[1][0][index_cent]->Fill(weight_pol*PolL.X());
	      	fhPolarization[1][1][index_cent]->Fill(weight_pol*PolL.Y());
	      	fhPolarization[1][2][index_cent]->Fill(weight_pol*PolL.Z());
	      	if(PolL.Mag()!=0.)fhProjection[1][0][index_cent]->Fill(adPL);
	      	fhProjection[1][1][index_cent]->Fill(adnL);
	      	fhProjection[1][2][index_cent]->Fill(admL);
		fhDeltaPhi[1][index_cent]->Fill(GetDeltaPhi(HypL,PosL,mL));
		if(
				(etamc < cut_eta) && (etamcn < cut_eta) && (ptmc > cut_pt) && (ptmcn > cut_pt)
				){
	      		fhPolarizationcut[1][0][index_cent]->Fill(weight_pol*PolL.X());
	      		fhPolarizationcut[1][1][index_cent]->Fill(weight_pol*PolL.Y());
	      		fhPolarizationcut[1][2][index_cent]->Fill(weight_pol*PolL.Z());
	      		if(PolL.Mag()!=0.)fhProjectioncut[1][0][index_cent]->Fill(adPL);
	      		fhProjectioncut[1][1][index_cent]->Fill(adnL);
	      		fhProjectioncut[1][2][index_cent]->Fill(admL);
			fhDeltaPhicut[1][index_cent]->Fill(GetDeltaPhi(HypL,PosL,mL));
                        hMassV0MC[1]->Fill(laInvmass);
                        fhArmPodMCcut->Fill(alpha, ptarm);
		}
	      }
            } // id negative daughter
	   } // second loop

          } // id positive daughter 
  }//end loop MC Tracks


 //Vertex info
fTVtx = event.fVertex;

//PrimaryVtx point
Int_t nVert = fTVtx->GetEntriesFast();
MpdVertex *vtx = (MpdVertex*) fTVtx->First();

Double_t posVtxX = vtx->GetX();
Double_t posVtxY = vtx->GetY();
Double_t posVtxZ = vtx->GetZ();

Double_t vertice=TMath::Sqrt(posVtxX*posVtxX + posVtxY*posVtxY);
TVector3 primvtx;
vtx->Position(primvtx);

     // fhtest2->Fill(vertice);


 //Reconstructed tracks 

fTDstEvent = event.fMPDEvent;
fTMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();
//TClonesArray *MpdGlobalTracks = (TClonesArray*)fTDstEvent->GetGlobalTracks();
Int_t ntracks=fTMpdGlobalTracks->GetEntriesFast();

fTpcKlTracks = event.fTPCKalmanTrack;

//variables for reconstructed tracks
vector<MpdParticle*> v0particle;
MpdV0track v0f;

 for (Int_t i = 0; i < ntracks-1; i++){
	 //if(i>20)continue;
	 MpdTrack *trackp = (MpdTrack*)fTMpdGlobalTracks->UncheckedAt(i);
      if (!selectTrackpr(trackp)) {
         continue;
      }

//	 if(trackp->GetChi2()< -8)continue;
         MpdTpcKalmanTrack *trp = (MpdTpcKalmanTrack*)fTpcKlTracks->UncheckedAt(i);

	 Double_t pt  = trackp->GetPt();
         Double_t eta = trackp->GetEta();
         Int_t nhits  = trackp->GetNofHits();

for (int iTrn=i+1; iTrn<ntracks; iTrn++){
	 MpdTrack *trackn = (MpdTrack*)fTMpdGlobalTracks->UncheckedAt(iTrn);
      if (!selectTrackpi(trackn)) {
         continue;
      }

//	 if(trackn->GetChi2()< -8)continue;
         MpdTpcKalmanTrack *trn = (MpdTpcKalmanTrack*)fTpcKlTracks->UncheckedAt(iTrn);
 
	 Double_t ptn  = trackn->GetPt();
         Double_t etan = trackn->GetEta();
         Int_t nhitsn  = trackn->GetNofHits();

      //opposite charge tracks
      int chargep, chargen;
      if (pt < 0)
	      chargep =  1;
      else
	      chargep = -1;
      if(chargep < 0) continue;
      if (ptn < 0)
	      chargen =  1;
      else
	      chargen = -1;
      if(chargen > 0) continue;
   //reject same sign pairs
   //   if(chargep*chargen > 0) continue;
     
//********************************************************************************************************
 
  MpdTpcKalmanTrack trCorK1(*trp);
  MpdHelix helixp = MakeHelix(trCorK1);
  MpdParticle       pr(trCorK1, 0);
  if(chargep > 0)pr.SetPdg(2212);
  if(chargep < 0)pr.SetPdg(-2212);
	    //masspos=0.938272;
  pr.SetMass(0.938272);

//********************************************************************************************************
  MpdTpcKalmanTrack trCorK2(*trn);
  MpdHelix helixn = MakeHelix(trCorK2);
  MpdParticle       pi(trCorK2, 0);
  if(chargen > 0)pi.SetPdg(211);
  if(chargep < 0)pi.SetPdg(-211);
  pi.SetMass(0.139570); // if I don't add the number appears a Break segmentation violation in MpdParticle::SetMass(double) ()

  //v0 daughters
  mPartV0.clear();
  mPartV0.emplace_back(&pr);
  mPartV0.emplace_back(&pi);

  MpdParticle V0La;
   float       chi2La = TMath::Abs(V0La.BuildMother(mPartV0));
   float       ptLa   = V0La.Pt();

  fhchi2La->Fill(chi2La);
  fhptLa->Fill(ptLa);

TVector3 v0(V0La.Getx()(0,0), V0La.Getx()(1,0), V0La.Getx()(2,0));
v0 -= primvtx;
//Double_t decay = v0.Mag();

   float rConv = TMath::Sqrt(pow(V0La.Getx()(0, 0), 2) + pow(V0La.Getx()(1, 0), 2));
      fhtest->Fill(rConv);
      fhrConv->Fill(rConv);
//   if (rConv < 5. || rConv > 120)continue; // Minimal and Maximal conversion radius (to exclude Dalitz and poorly reconstructed tracks respectively) check if values apply for V0s Dalitz plot three body decay but radius? could be

   float angle;
   angle =v0.Angle(V0La.Momentum3());

   std::pair<float, float> paths = helixp.pathLengths(helixn);
   TVector3                p1    = helixp.at(paths.first);
   TVector3                p2    = helixn.at(paths.second);
   p1 -= p2;
   float dist = p1.Mag(); // Closest distance between daughters
   if(dist > 2.8)continue;

   float massV0La = V0La.GetMass();

   if(massV0La < 1.08 && massV0La > 1.2) continue;

     // propagate trCorK1,trCorK2 to conversion point
  MpdKalmanHit hitTmp;
   hitTmp.SetType(MpdKalmanHit::kFixedR);
   hitTmp.SetPos(trCorK1.GetPos());
   trCorK1.SetParamNew(*trCorK1.GetParam());
   trCorK1.SetPos(trCorK1.GetPosNew());
   trCorK1.ReSetWeight();
   //  TMatrixDSym w = *trCorK1.GetWeight(); // save current weight matrix
   mKHit.SetPos(rConv);
   if (!mKF->PropagateToHit(&trCorK1, &mKHit, kFALSE, kFALSE))continue;
   trCorK1.SetDirection(MpdKalmanTrack::kInward);
   TVector3 m1 = trCorK1.Momentum3();

   hitTmp.SetPos(trCorK2.GetPos());
   trCorK2.SetParamNew(*trCorK2.GetParam());
   trCorK2.SetPos(trCorK2.GetPosNew());
   trCorK2.ReSetWeight();
   TMatrixDSym w = *trCorK1.GetWeight(); // save current weight matrix
   mKHit.SetPos(rConv);
   if (!mKF->PropagateToHit(&trCorK2, &mKHit, kFALSE, kFALSE))continue;
   trCorK2.SetDirection(MpdKalmanTrack::kInward);
   TVector3 m2 = trCorK2.Momentum3();

   // Asymmetry cut
   float asym1 = m1.Mag() / V0La.Momentum();
   float asym2 = m2.Mag() / V0La.Momentum();

fhAsym1->Fill(asym1);
fhAsym2->Fill(asym2);


  TVector3 pca_p;
  Double_t s_p = helixp.pathLength(primvtx);
  pca_p = helixp.at(s_p);
  pca_p -= primvtx;
  Double_t dcapos = pca_p.Mag();

  TVector3 pca_n;
  Double_t s_n = helixn.pathLength(primvtx);
  pca_n = helixn.at(s_p);
  pca_n -= primvtx;
  Double_t dcaneg = pca_n.Mag();


v0particle.clear();

v0particle.emplace_back(&pr);
v0particle.emplace_back(&pi);

MpdParticle lambPart;
Double_t chi2 = TMath::Abs (lambPart.BuildMother(v0particle));


Double_t cosangle = TMath::Cos(angle);

      TVector3 Pmpos(trackp->GetPx(),trackp->GetPy(),trackp->GetPz());
      TVector3 Pmneg(trackn->GetPx(),trackn->GetPy(),trackn->GetPz());
      TVector3 PLa=Pmpos+Pmneg;

      Double_t masspos=0;
      Double_t massneg=0;
  
     //ArmPod variables
 
      Float_t       lQlNeg = Pmneg.Dot(PLa)/PLa.Mag();
      Float_t       lQlPos = Pmpos.Dot(PLa)/PLa.Mag();
      Float_t       alpha = (lQlPos - lQlNeg)/(lQlPos + lQlNeg);
      Float_t       ptarm = Pmneg.Perp(PLa);

    //Lambda reco
	    masspos=0.938272;
	    massneg=0.139570;

      Double_t Enpos=TMath::Sqrt(masspos*masspos + Pmpos.Mag2());
      Double_t Enneg=TMath::Sqrt(massneg*massneg + Pmneg.Mag2());

      Double_t laInvmass=TMath::Sqrt(masspos*masspos + massneg*massneg +2*(Enpos*Enneg - Pmpos.Dot(Pmneg)));

Double_t decay = v0.Mag();
Double_t EnLa=TMath::Sqrt(PLa.Mag2() + 1.245456);
TLorentzVector LambdaLV(PLa,EnLa);
TLorentzVector ProtonLV(Pmpos,Enpos);
TVector3 mL(TMath::Sin(PsiEP),-TMath::Cos(PsiEP),0);
TVector3 mLR(TMath::Sin(PsizdcEP),-TMath::Cos(PsizdcEP),0); 
Double_t ptLam=TMath::Sqrt(PLa.X()*PLa.X() + PLa.Y()*PLa.Y());
TVector3 nL(-PLa.Y()/ptLam,PLa.X()/ptLam,0);


Double_t ctau=Getctau(LambdaLV,decay);
Double_t deltaphipr=GetDeltaPhi(LambdaLV,ProtonLV,mL);
Double_t deltapsipr=GetDeltaPhi(LambdaLV,ProtonLV,mLR);


Double_t adnL=GetProjection(LambdaLV,ProtonLV,nL);
Double_t admLrec=GetProjection(LambdaLV,ProtonLV,mLR);


  if(dist < 0.5){ 
  //if(dcapos < 0.5)continue;
  //if(dcaneg < 0.5)continue;
  //if((asym1 < 0.6 || asym1 > 1.0) && (asym2 > 0.4 || asym2 < 0.05))continue;
if(
  (ptLa > 0.005  ) &&
  (chi2La < 10   ) &&
  (angle < 0.102 ) &&
  (dcapos > 0.15 ) &&
  (dcaneg > 0.15 ) &&
  (asym1 > 0.6   ) && 
  (asym1 < 1.0   ) && 
  (asym2 > 0.05  ) && 
  (asym2 < 0.4   )  
		 ){

      hInvMass->Fill(laInvmass);
      fhArmPod->Fill(alpha, ptarm);
      fhCosInvMass->Fill(laInvmass,cosangle);
      fhDCAv0InvMass->Fill(laInvmass,dist);
      fhDCAposInvMass->Fill(laInvmass,dcapos);
      fhDCAnegInvMass->Fill(laInvmass,dcaneg);
      fhDecayctau->Fill(decay,ctau);
		fhDeltaPhiRec[0][index_cent]->Fill(deltaphipr);
		fhDeltaPsiRec[0][index_cent]->Fill(deltapsipr);
		fhDeltaPhiRecLa->Fill(deltaphipr);
		fhProjLocRec[0][index_cent]->Fill(adnL);
		fhProjPsiRec[0][index_cent]->Fill(admLrec);
		fhProjLocRecLa->Fill(adnL);
	}
      }//cut on dcabetweenV0

	//Get MC ID
      auto mctrackp = (MpdMCTrack*)fTMCTracks->UncheckedAt(trackp->GetID());

      auto mctrackn = (MpdMCTrack*)fTMCTracks->UncheckedAt(trackn->GetID());
      Int_t pdgpos = mctrackp->GetPdgCode();
      if(mctrackp->GetMotherId()==-1)continue; // remove primaries and keep only secondaries
      if(mctrackn->GetMotherId()==-1)continue; // remove primaries and keep only secondaries
      if (mctrackn->GetMotherId()!= mctrackp->GetMotherId())continue;
      Int_t pdgneg = mctrackn->GetPdgCode();
	//if(!(pdgpos == 2212 || pdgpos == 211))continue;
	//if(!(pdgneg == -211 || pdgneg == -2212))continue;
     auto mctrackmomp = (MpdMCTrack*) fTMCTracks->UncheckedAt(mctrackp->GetMotherId());
     auto mctrackmomn = (MpdMCTrack*) fTMCTracks->UncheckedAt(mctrackn->GetMotherId());

     Int_t mompdg = mctrackmomp->GetPdgCode();
     if(TMath::Abs(mompdg)!=3122)continue;
     //End MC ID
      hInvMassAss->Fill(laInvmass);
      fhArmPodAss->Fill(alpha, ptarm);
      fhCosInvMassAss->Fill(laInvmass,cosangle);
      fhDCAv0InvMassAss->Fill(laInvmass,dist);
      fhDCAposInvMassAss->Fill(laInvmass,dcapos);
      fhDCAnegInvMassAss->Fill(laInvmass,dcaneg);
      fhDecayctauAss->Fill(decay,ctau);
      fhtest2->Fill(admLrec);
		fhDeltaPhiAss[0][index_cent]->Fill(deltaphipr);
		fhDeltaPsiAss[0][index_cent]->Fill(deltapsipr);
		fhDeltaPhiAssLa->Fill(deltaphipr);
		fhProjLocAss[0][index_cent]->Fill(adnL);
		fhProjPsiAss[0][index_cent]->Fill(admLrec);
		fhProjLocAssLa->Fill(adnL);

  fhchi2LaAss->Fill(chi2La);
  fhptLaAss->Fill(ptLa);
  fhrConvAss->Fill(rConv);
fhAsymAss1->Fill(asym1);
fhAsymAss2->Fill(asym2);

  }//second loop end
 }//end loop reconstructed tracks

}
//___________________________________________________________________________________
void MpdV0AnalysisTask::Finish()
{
}
//___________________________________________________________________________________
bool MpdV0AnalysisTask::selectEvent(MpdAnalysisEvent &event)
{
   mhEvents->Fill(0.5);
   // first test if event filled?
   if (!event.fVertex) { // if even vertex not filled, skip event
      return false;
   }
   // Vertex z coordinate
   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);
   mhVertex->Fill(mPrimaryVertex.Z());
   float mZvtxCut = 50.;//vertex cut to be implemented in add class
   if (fabs(mPrimaryVertex.Z()) > mZvtxCut) { 
      return false;
   }
   mZvtxBin = 0.5 * (mPrimaryVertex.Z() / mZvtxCut + 1) * nMixEventZ;
   if (mZvtxBin < 0) mZvtxBin = 0;
   if (mZvtxBin >= nMixEventZ) mZvtxBin = nMixEventZ - 1;
   mhEvents->Fill(1.5);

   float cen = event.getCentrTPC();
   mCenBin   = (cen / 100.) * nMixEventCent; // very rough
   if (mCenBin < 0) mCenBin = 0;
   if (mCenBin >= nMixEventCent) mCenBin = nMixEventCent - 1;

   // Multiplicity
   fTMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();
   int ntr          = fTMpdGlobalTracks->GetEntriesFast();
   mhMultiplicity->Fill(ntr);

   // Centrality
   mhCentrality->Fill(mCenBin);
   // mCenBin = 0;
   mhEvents->Fill(2.5);

   // ZCD vs TPC (pileup?)
   mhEvents->Fill(3.5);

   // Eventplane  TODO
   mRPBin = 0;
   mhEvents->Fill(4.5);

   mixBin = (mCenBin + 1) * (mZvtxBin + 1) * (mRPBin + 1);
   // cout<<"Mixing bin: "<<mixBin<<" = "<<mCenBin<<" "<<mZvtxBin<<" "<<mRPBin<<endl;

   return true;
}
//___________________________________________________________________________________
bool MpdV0AnalysisTask::selectTrackpr(MpdTrack *track)
{

	 Double_t pt  = track->GetPt();
         Double_t eta = track->GetEta();
         Int_t nhits  = track->GetNofHits();


	 if(pt > 0 ) return false; // reject negative charge

      if (TMath::Abs(pt) < cut_pt) return false;
      if (TMath::Abs(eta) > cut_eta) return false;
      if (nhits < cut_nhits) return false;

      bool isGoodPID;
      if (track->GetTofFlag()== 2 || track->GetTofFlag()==6){
      isGoodPID = mPID->FillProbs(TMath::Abs(pt)*TMath::CosH(eta),track->GetdEdXTPC()*6.036e-3,track->GetTofMass2(),1);
      } else {
      isGoodPID = mPID->FillProbs(TMath::Abs(pt)*TMath::CosH(eta),track->GetdEdXTPC()*6.036e-3,1);
      }
      if (isGoodPID && (mPID->GetProbPr() < 0.75)) {
	      return false;
      }
      float dEdx = track->GetdEdXTPC();
      float tofmass2=track->GetTofMass2();
      mhdEdx->Fill(TMath::Abs(pt)*TMath::CosH(eta), dEdx); // | p | = p_T cosh(η)
      mhdEdxvsmass2->Fill(track->GetTofMass2(), dEdx);
      long int prim1 = track->GetID();
      if (
		      (abs((static_cast<MpdMCTrack *>(fTMCTracks->At(prim1)))->GetPdgCode()) == 2212)		      
		      ){
	      mhdEdxvsmass2Ass->Fill(track->GetTofMass2(),dEdx);
      	      mhdEdxAss->Fill(TMath::Abs(pt)*TMath::CosH(eta),dEdx);
              mhdEdxvspAssProf->Fill(TMath::Abs(pt)*TMath::CosH(eta), dEdx); // | p | = p_T cosh(η)
      }
/*
      float pionC, electronC, kaonC, protonC, maxprobC;
      int pidpdg, charge;
 	 pionC = track->GetPidProbPion(); 
	 electronC = track->GetPidProbElectron();
	 protonC = track->GetPidProbProton();
	 kaonC = track->GetPidProbKaon();
 	 Double_t ProbsC[] = {0,electronC, pionC, kaonC, protonC};
 	 maxprobC =  TMath::MaxElement(5, ProbsC);
      int idmax3 = TMath::LocMax(5,ProbsC);
*/

      if(tofmass2 > 0.046722 && tofmass2 < 1.22858 && // gaussianfit +4sigma at low dedx; previous values // (0.6804,1.0793)
        (track->GetTofFlag() == 2 || track->GetTofFlag() == 6)
		      ){
	      return true;
      }else {
              return false;
      }

/*
      if ((fabs(dEdx) < 3.0) && (fabs(mPID->GetNsigmaToBetheBloch(MpdPidUtils::kProton)) < 3.0 &&
      //if ((fabs(dEdx) < 3.0) && (fabs(Beta_sigma(track->GetTofBeta(), sqrt(pow(pt, 2) + pow(track->GetPz(), 2)))) < 3.0 &&
        (track->GetTofFlag() == 2 || track->GetTofFlag() == 6))) {
      return true;
   } else {
      return false;
   }*/
return false;
}
//______________________________________________________________________________________________
bool MpdV0AnalysisTask::selectTrackpi(MpdTrack *track)
{

	 Double_t pt  = track->GetPt();
         Double_t eta = track->GetEta();
         Int_t nhits  = track->GetNofHits();

	 if(pt < 0 ) return false; // reject positive charge

      if (TMath::Abs(pt) < cut_pt) return false;
      if (TMath::Abs(eta) > cut_eta) return false;
      if (nhits < cut_nhits) return false;

      bool isGoodPID;
      if (track->GetTofFlag()== 2 || track->GetTofFlag()==6){
      isGoodPID = mPID->FillProbs(TMath::Abs(pt)*TMath::CosH(eta),track->GetdEdXTPC()*6.036e-3,track->GetTofMass2(),-1);
      } else {
      isGoodPID = mPID->FillProbs(TMath::Abs(pt)*TMath::CosH(eta),track->GetdEdXTPC()*6.036e-3,-1);
      }
      if (isGoodPID && (mPID->GetProbPi() < 0.75)) {
	      return false;
      }
      float dEdx = track->GetdEdXTPC();
      float tofmass2=track->GetTofMass2();
      mhdEdxpi->Fill(TMath::Abs(pt)*TMath::CosH(eta), dEdx);
      mhdEdxvsmass2pi->Fill(track->GetTofMass2(), dEdx);
      long int prim1 = track->GetID();
      if (
		      (abs((static_cast<MpdMCTrack *>(fTMCTracks->At(prim1)))->GetPdgCode()) == 211)
		      ){
	      mhdEdxvsmass2piAss->Fill(track->GetTofMass2(),dEdx);
      	      mhdEdxpiAss->Fill(TMath::Abs(pt)*TMath::CosH(eta),dEdx);
              mhdEdxvsppiAssProf->Fill(TMath::Abs(pt)*TMath::CosH(eta), dEdx); // | p | = p_T cosh(η)
      }
/*
      float pion, electron, kaon, proton, maxprob, fProbCut = 0.6;
      float pionT, electronT, kaonT, protonT, maxprobT;
      float pionC, electronC, kaonC, protonC, maxprobC;
      int pidpdg, charge;

  pion = track->GetTPCPidProbPion(); 
  electron = track->GetTPCPidProbElectron();
  proton = track->GetTPCPidProbProton();
  kaon = track->GetTPCPidProbKaon();
  Double_t Probs[] = {0,electron, pion, kaon, proton};
  maxprob =  TMath::MaxElement(5, Probs);
  pidpdg = mPID->GetMaxProb();

  pionT = track->GetTOFPidProbPion(); 
  electronT = track->GetTOFPidProbElectron();
  protonT = track->GetTOFPidProbProton();
  kaonT = track->GetTOFPidProbKaon();
  Double_t ProbsT[] = {0,electronT, pionT, kaonT, protonT};
  maxprobT =  TMath::MaxElement(5, ProbsT);

  pionC = track->GetPidProbPion(); 
  electronC = track->GetPidProbElectron();
  protonC = track->GetPidProbProton();
  kaonC = track->GetPidProbKaon();
  Double_t ProbsC[] = {0,electronC, pionC, kaonC, protonC};
  maxprobC =  TMath::MaxElement(5, ProbsC);

int idmax1 = TMath::LocMax(5,Probs);
int idmax2 = TMath::LocMax(5,ProbsT);
int idmax3 = TMath::LocMax(5,ProbsC);
*/

      if(tofmass2 > 0.005755 && tofmass2 < 0.03123 &&
        (track->GetTofFlag() == 2 || track->GetTofFlag() == 6) //&&
		      ){
	      return true;
      }else {
              return false;
      }
  
/*
      if ((fabs(dEdx) < 3.0) && (fabs(mPID->GetNsigmaToBetheBloch(MpdPidUtils::kProton)) < 3.0 &&
      //if ((fabs(dEdx) < 3.0) && (fabs(Beta_sigma(track->GetTofBeta(), sqrt(pow(pt, 2) + pow(track->GetPz(), 2)))) < 3.0 &&
        (track->GetTofFlag() == 2 || track->GetTofFlag() == 6))) {
      return true;
   } else {
      return false;
   }*/
return false;
}
//___________________________________________________________________________________

//___________________________________________________________________________________
Int_t MpdV0AnalysisTask::GetMultiplicity(MpdAnalysisEvent &event,Double_t ptCut, Double_t etaCut, Double_t dcaCut, Int_t nhitsCut, Bool_t isMC) {
//Int_t MpdV0AnalysisTask::GetMultiplicity(MpdAnalysisEvent &event) {

	Int_t refMult = -10;
       
if(isMC==kTRUE){

 fTMCTracks = event.fMCTrack;
 Int_t nmctracks=fTMCTracks->GetEntriesFast();

 for (Int_t i = 0; i < nmctracks; i++)
  {
    MpdMCTrack *MCtrack = (MpdMCTrack*) fTMCTracks->UncheckedAt(i);
    Double_t ptmc=MCtrack->GetPt();
    TVector3 P(MCtrack->GetPx(),MCtrack->GetPy(),MCtrack->GetPz());
    Double_t etamc=0.5*TMath::Log((P.Mag() + MCtrack->GetPz())/(P.Mag() - MCtrack->GetPz()+1.e-13));
    Int_t abspdg=TMath::Abs(MCtrack->GetPdgCode());


      if(MCtrack->GetMotherId()!=-1) continue;
      if (ptmc < ptCut) continue;
      if(TMath::Abs(etamc) > etaCut) continue;
      if(abspdg == 211 || abspdg == 321 || abspdg == 2212)    refMult++;

  }
}
else{

fTDstEvent = event.fMPDEvent;
TClonesArray *MpdGlobalTracks = (TClonesArray*)fTDstEvent->GetGlobalTracks();
Int_t ntracks=MpdGlobalTracks->GetEntriesFast();
 for (Int_t i = 0; i < ntracks; i++){
	 MpdTrack *track = (MpdTrack*)MpdGlobalTracks->UncheckedAt(i);

	 Double_t pt  = track->GetPt();
	 Double_t eta = track->GetEta();
	 Int_t nhits  = track->GetNofHits();

      if (TMath::Abs(pt) < ptCut) continue;
      if (TMath::Abs(eta) > etaCut) continue;
      if (nhits < nhitsCut) continue;
      if (TMath::Sqrt(TMath::Power(track->GetDCAX(),2) + TMath::Power(track->GetDCAY(),2) + TMath::Power(track->GetDCAZ(),2)) > dcaCut) continue;
      refMult++;
  }
}
return refMult;
}
//_______________________________________________________________________________
Double_t MpdV0AnalysisTask::GetProjection(TLorentzVector La, TLorentzVector pr, TVector3 Proy){

	Double_t projection = -2.0; 

TVector3 PLa=La.Vect();
Double_t EnLa=La.E();

TVector3 P=pr.Vect();
Double_t Enpos=pr.E();

      TVector3 lBetta = (1./EnLa)*PLa;
      Double_t MagBetta = lBetta.Mag();
      Double_t lgamma  = 1.0/TMath::Sqrt(1.0 - (MagBetta*MagBetta));

//      Double_t ldotLL = (PLa.Dot(PLa))/EnLa;
//      TVector3 lPLaRF = PLa + ((lgamma-1.0)/::pow(MagBetta,2))*ldotLL*lBetta - lgamma*lBetta*EnLa;

      Double_t ldotPB = (PLa.Dot(P))/EnLa;
      TVector3 lPRF = P + ((lgamma-1.0)/::pow(MagBetta,2))*ldotPB*lBetta - lgamma*lBetta*Enpos;

      projection = TMath::Cos(lPRF.Angle(Proy));
//projection = EnLa;
//projection = lPLaRF.Mag();
return projection;
}
//_______________________________________________________________________________
Double_t MpdV0AnalysisTask::GetDeltaPhi(TLorentzVector La, TLorentzVector pr, TVector3 EPv){

Double_t deltaphi = -99.0; 

TVector3 PLa=La.Vect();
Double_t EnLa=La.E();

TVector3 P=pr.Vect();
Double_t Enpos=pr.E();

      TVector3 lBetta = (1./EnLa)*PLa;
      Double_t MagBetta = lBetta.Mag();
      Double_t lgamma  = 1.0/TMath::Sqrt(1.0 - (MagBetta*MagBetta));

      Double_t ldotPB = (PLa.Dot(P))/EnLa;
      TVector3 lPRF = P + ((lgamma-1.0)/::pow(MagBetta,2))*ldotPB*lBetta - lgamma*lBetta*Enpos;

     //proton azimutal angle
     Double_t phip=TMath::ATan2(-lPRF.Y(),-lPRF.X()) + TMath::Pi();
     Double_t PEL=TMath::ATan2(-EPv.X(),EPv.Y()) + TMath::Pi();

    // deltaphi = phip;
     deltaphi = PEL - phip;
     if(deltaphi < 0)deltaphi = deltaphi + 2*TMath::Pi();


return deltaphi;
}
//_______________________________________________________________________________
Int_t MpdV0AnalysisTask::GetCentClass(Int_t multiplicityE){ // Fit results from MC-Glauber method(Primary tracks) 
							   // https://mpdforum.jinr.ru/uploads/short-url/fuR0zSFXPfAh4DDxlCE5nlemk5m.pdf
	Int_t centclass;
if(multiplicityE >= 115){
	centclass = 0;
}
else if(multiplicityE >= 83 && multiplicityE < 115){
	centclass = 1;
}
else if(multiplicityE >= 59 && multiplicityE < 83){
	centclass = 2;
}
else if(multiplicityE >= 41 && multiplicityE < 59){
	centclass = 3;
}
else if(multiplicityE >= 28 && multiplicityE < 41){
	centclass = 4;
}
else if(multiplicityE >= 18 && multiplicityE < 28){
	centclass = 5;
}
else if(multiplicityE >= 11 && multiplicityE < 18){
	centclass = 6;
}
else if(multiplicityE >= 6 && multiplicityE < 11){
	centclass = 7;
}
else if(multiplicityE >= 3 && multiplicityE < 6){
	centclass = 8;
}
else if(multiplicityE >= 1 && multiplicityE < 3){
	centclass = 9;
}
else {
	centclass = 10;
}

return centclass;
}
//______________________________________________________________________________________
bool MpdV0AnalysisTask::BuildLambda(MpdTrack *trp, MpdTpcKalmanTrack *ktrp, MpdTrack *trn, MpdTpcKalmanTrack *ktrn, MpdV0track &v0can)
{
   // reject same sign pairs
   if (trp->GetPt() * trn->GetPt() > 0) return false;
cout<<"valor de C pos = "<< trp->GetCharge() << ",   valor de C neg = "<< trn->GetCharge() <<endl;

  MpdTpcKalmanTrack trCorK1(*ktrp);
  MpdHelix helixp = MakeHelix(trCorK1);
  MpdParticle       pos1(trCorK1);

  MpdTpcKalmanTrack trCorK2(*ktrn);
  MpdHelix helixn = MakeHelix(trCorK2);
  MpdParticle       neg2(trCorK2);

mPartV0.clear();

return true;
}
//______________________________________________________________________________________
Double_t MpdV0AnalysisTask::EventPlaneZDC(MpdAnalysisEvent &event){

	Double_t zdcEP = 0;

/*************************************************************/
	// position of ZDC modules
	
Int_t nModulesZDC=90;
Float_t ZDC_energy_mpd[nModulesZDC];
Float_t phi_angle_of_module[nModulesZDC];
for (long int j = 0; j < 90; ++j){
                        ZDC_energy_mpd[j] = 0;
                }
        for (int i = 0; i < 2; ++i)
        {
                Int_t x_axis_switch;
                if (i == 0) x_axis_switch = 1;
                else if (i == 1) x_axis_switch = -1;

                for (Int_t j = 0; j < nModulesZDC/2; ++j)
                {
                        Double_t y = 0, x = 0;

                        if ((j>=0) && (j<=4))
                        {
                                y = 45., x = (j-2)*15.;
                                phi_angle_of_module[j + i*nModulesZDC/2] = TMath::ATan2(y,x_axis_switch*x);
                        }
                        else if ((j>=5) && (j<=39))
                        {
                                y = (3-(j+2)/7)*15, x = (3-(j+2)%7)*15;
                                phi_angle_of_module[j + i*nModulesZDC/2] = TMath::ATan2(y,x_axis_switch*x);
                        }
                        else if ((j>=40) && (j<=44))
                        {
                                y = -45. , x = (j-42)*15.;
                                phi_angle_of_module[j + i*nModulesZDC/2] = TMath::ATan2(y,x_axis_switch*x);
                        }
                }
        }
/*****************************************************************************/

 fTZDCDigits = event.fZDCDigit; 
 Int_t nZdc = fTZDCDigits->GetEntriesFast();

for(Int_t j=0;j < nZdc; ++j){
                MpdZdcDigi* ZDCHit = (MpdZdcDigi*) fTZDCDigits->At(j);
                Int_t detector_ID = ZDCHit->GetDetectorID();//1,2
                Int_t module_ID = ZDCHit->GetModuleID();//1-45
                Double_t energyhit = ZDCHit->GetELoss();
                ZDC_energy_mpd[ (detector_ID - 1)*45 + module_ID - 1] += energyhit;
}

Double_t Qx_zdc_p=0, Qx_zdc_m=0;
Double_t Qy_zdc_p=0, Qy_zdc_m=0;
for(Int_t j=0; j<45; ++j){
Qx_zdc_p += ZDC_energy_mpd[j]*TMath::Cos(phi_angle_of_module[j]);
Qy_zdc_p += ZDC_energy_mpd[j]*TMath::Sin(phi_angle_of_module[j]);
}
for(Int_t j=45; j<90; ++j){
Qx_zdc_m += ZDC_energy_mpd[j]*TMath::Cos(phi_angle_of_module[j]);
Qy_zdc_m += ZDC_energy_mpd[j]*TMath::Sin(phi_angle_of_module[j]);
}

Double_t QQxzdc1 = Qx_zdc_p - Qx_zdc_m;
Double_t QQyzdc1 = Qy_zdc_p - Qy_zdc_m;

  zdcEP = TMath::ATan2(-QQyzdc1,-QQxzdc1) + TMath::Pi();

  return zdcEP;

}
//______________________________________________________________________________________
MpdHelix MpdV0AnalysisTask::MakeHelix(const MpdKalmanTrack &tr) const
{

  Double_t r = tr.GetPosNew();
  Double_t phi = tr.GetParam(0) / r;
  Double_t x = r * TMath::Cos(phi);
  Double_t y = r * TMath::Sin(phi);
  Double_t dip = tr.GetParam(3);
  Double_t cur = 0.3 * 0.01 * 5 / 10; // 5 kG
  cur *= TMath::Abs (tr.GetParam(4));
  TVector3 o(x, y, tr.GetParam(1));
  Int_t h = (Int_t) TMath::Sign(1.1,tr.GetParam(4));
  MpdHelix helix(cur, dip, tr.GetParam(2)-TMath::PiOver2()*h, o, h);
  return helix;
}
//______________________________________________________________________________________
Double_t MpdV0AnalysisTask::Getctau(TLorentzVector La, Double_t distance)
{

	// distance travel in lab frame d = betta*gamma*ctau ctau=7.89cm for Lambda
Double_t ct = 0;

TVector3 PLa=La.Vect();
Double_t EnLa=La.E();

      TVector3 lBetta = (1./EnLa)*PLa;
      Double_t MagBetta = lBetta.Mag();
      Double_t lgamma  = 1.0/TMath::Sqrt(1.0 - (MagBetta*MagBetta));

ct = distance/(MagBetta * lgamma);
//ct = MagBetta * lgamma * 7.89;

return ct;
}
