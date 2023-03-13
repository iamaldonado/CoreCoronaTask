# CoreCorona simply task to get &Lambda polarization

This working example, is based on the lego trains introduced by V.Riabov [presentation](https://indico.jinr.ru/event/3391/contributions/18477/attachments/13910/23277/AnalysisFramework_RiabovV.pdf) at MPD Cross-PWG Meeting. It implements the measurement of &Lambda polarization within the core-corona approach.

It employs the classes MpdAnalysisManager, MpdAnalysisEvent, and MpdAnalysisTask commited by D. Peresunko originally. It runs in the tag v22.12.22.

## MpdAnalysisEvent Class

Contains references to all branches in the DST file for this event and may contain extra global variables of interest (centrality, event plane, T<sub>0</sub>, n-sigma matching ... etc.). In the actual version reads

```ruby
EventHeader
Vertex
MPDEvent
MCEventHeader
MCTrack
EMCCluster
TPCKalmanTrack
TOFHit
TOFMatching
ZdcDigi
ZDCEloss1Value
ZDCEloss2Value
ZDCEloss1Histo
ZDCEloss2Histo
```

## MpdAnalysisManager Class

Process Input File and allow to add task to the analysis


## MpdAnalysisTas - Basic Structure

Tasks called by the MpdAnalysis Manager should be derived from MpdAnalysisTask. It consists in two files, its own header (.h) and implementation file (.cxx). In the following lines the basic structure is described.

- ** Constructors and destructors

It contains the standard  C++ constructors and destructors, called each time an Instance of the class is created or deleted, they are defined in the header file .h

```ruby
   MpdAnalysisTask() = default;
   MpdAnalysisTask(const char *name, const char *outputName = "taskName");
   virtual ~MpdAnalysisTask(); // Destructor
```

It also implements the following methods that needs to be re-defined in each class

- ** Initialization of objects
Users should prepare objects to fill (histograms, trees, etc.)
```ruby
	void UserInit();
```
- ** Main - Execution and process of analysis
Method is called for each event, data are provided by container MpdAnalysisEvent
```ruby
	void ProcessEvent(MpdAnalysisEvent &event);
```
- ** End
Method is called when scan is finished but class data are not written yet
```ruby
	void Finish();
```
- ** Output
It define the output file and the objects to be stored in it. In the example is a TList with TH1F and TH2F histograms.

```ruby
     void setOutFile(std::string filename = "histos.root") { mOutFile = filename; }
  private:
     std::string mOutFile = "histos.root" ;
```
## Common methods for the different examples (pairKK, photons)

We implement the common functions to select only events and tracks that satisfy certain parameters: z vertex position, &eta, p<sub>T</sub>, number of hits, PID, ...

### bool selectEvent(MpdAnalysisEvent &event)

Is defined in the event header as a protected member class, it plot histograms like multiplicity and z vertex position that should be defined also in the header file as is shown in the following lines

<details>
<summary>Click me</summary>

```ruby
protected:

bool selectEvent(MpdAnalysisEvent &event);  

private:

   // event properties
   bool     isInitialized = false;
   int      mCenBin       = 0;
   int      mZvtxBin      = 0;
   int      mRPBin        = 0;
   int      mixBin;

   TVector3 mPrimaryVertex;
   static constexpr short nMixEventZ    = 10; //(V) number of bins in z direction
   static constexpr short nMixEventCent = 10; //(V) number of bins of centrality
   static constexpr short nMixEventRP   = 1;  //(V) number of bins of Reaction Plane orientation

   // General QA  ==>>> taken from example class
   //
   TH1F *mhEvents       = nullptr;
   TH1F *mhVertex       = nullptr;
   TH1F *mhCentrality   = nullptr;
   TH1F *mhMultiplicity = nullptr;
```
</details>

The Implementation is written in file .cxx after finish() method:


<details>
<summary>Click me</summary>

```ruby
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
```
</details>


### bool selectTrack(MpdTrack &track)

It select events that pass cuts on |&eta|, p<sub>T</sub>, nhits and dca and applies PID selection with [MpdPID](https://git.jinr.ru/nica/mpdroot/-/tree/dev/core/mpdPid) class. 
Is defined in the header file by:

```ruby
 bool selectTrackpr(MpdTrack *track);
 
   TH2F *mhdEdx;
   TH2F *mhdEdxAss;
   TH2F *mhdEdxvsmass2;
   TH2F *mhdEdxvsmass2Ass;
```
It plots PID histograms, as is shown in the figures for p and &pi<sup>-</sup> 
<image 
src="/figures/protondedx.jpg" src="/figures/piondedx.jpg" alt="Energy loss proton and pion">


<details>
<summary>Click me</summary>

```ruby
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
      	      mhdEdxAss->Fill(TMath::Abs(pt)*TMath::CosH(eta),dEdx);// | p | = p_T cosh(η)
      }

      if(tofmass2 > 0.046722 && tofmass2 < 1.22858 && // gaussianfit +4sigma at low dedx; previous values // (0.6804,1.0793)
        (track->GetTofFlag() == 2 || track->GetTofFlag() == 6)
		      ){
	      return true;
      }else {
              return false;
      }

return false;
}
```
</details>








________________________________________________
MpdPID is initialized in the member class **ProcessEvent** with

To select the events, this function is called in the member class **ProcessEvent** with:

```ruby
if(!selectEvent(event)){
   return;
   }
```
```ruby
if (!isInitialized) {
      mPID          = new MpdPid(2.0, 2.0, 9.2, 1.0, "NSIG", "CFHM", "pikapr");
      isInitialized = true;
   }
```

 









 

