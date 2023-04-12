//         Funtion to check centrality files
//__________________________________________________________________
bool CheckFileExist(TString fileName){
    gSystem->ExpandPathName(fileName);
    if (gSystem->AccessPathName(fileName.Data()) == true)
    {
        cout<<endl<<"no specified file: "<<fileName<<endl;
        return false;
    }                

    return true;
}
//____________________________________________________________________
void RunAnalyses(int nEvents = -1, TString inFileList = "list10.txt"){

  gROOT->LoadMacro("mpdloadlibs.C"); // 
  gROOT->ProcessLine("mpdloadlibs()"); // you can add here the functions directly as appear in mpdloadlibs file

  TStopwatch timer;
  timer.Start();

  ProcInfo_t proc;
  MemInfo_t memory;

   MpdAnalysisManager man("ManagerAnal", nEvents) ;
if (!CheckFileExist(inFileList)) return;
	man.InputFileList(inFileList) ;
   	man.ReadBranches("*") ; 
   	man.SetOutput("histos.root") ;
   
  	MpdCentralityAll pCentr("pCentr","pCentra") ;
   	man.AddTask(&pCentr) ;
   
   MpdEventPlaneAll pEP("pEP","pEPa") ;
   man.AddTask(&pEP) ;

   MpdV0AnalysisTask MpdV0AnalysisTask("V0","V0a") ;
   man.AddTask(&MpdV0AnalysisTask) ;

   man.Process() ;


  gSystem->GetProcInfo(&proc);
  cout << " User CPU time: " << proc.fCpuUser << "seconds, Memory: resident " << proc.fMemResident << "KB, virtual " << proc.fMemVirtual << "KB"<< endl;

  gSystem->GetMemInfo(&memory);
  cout << " Total RAM: "<< memory.fMemTotal << "MB, Used RAM: " << memory.fMemUsed << " MB, Free RAM: " << memory.fMemFree << endl; 

   timer.Stop();
   Double_t rtime = timer.RealTime(), ctime = timer.CpuTime();
   printf("RealTime=%f seconds, CpuTime=%f seconds\n", rtime, ctime);
   cout << "Macro finished successfully." << endl;

ofstream myfile;//file that stores time and memory consumption you can comment
myfile.open("/scratch1/maldonado/test/info.txt");
myfile<< "Time and memory usage" << "\n";

  myfile << " User CPU time: " << proc.fCpuUser << "seconds, Memory: resident " << proc.fMemResident << "KB, virtual " << proc.fMemVirtual << "KB"<< "\n";
  myfile << " Total RAM: "<< memory.fMemTotal << "MB, Used RAM: " << memory.fMemUsed << " MB, Free RAM: " << memory.fMemFree << "\n"; 
  myfile << "RealTime = " << rtime << " seconds, CpuTime = " << ctime << " seconds" << "\n";

  myfile.close();
}
