void RunAnalyses(){

  gROOT->LoadMacro("mpdloadlibs.C"); // looks unnecessary, but it is to load at list MpdParticle
  gROOT->ProcessLine("mpdloadlibs()");

  TStopwatch timer;
  timer.Start();

  ProcInfo_t proc;
  MemInfo_t memory;

   MpdAnalysisManager man("ManagerAnal") ;
 //man.InputFileList("list.txt") ;
// man.InputFileList("onefile.txt") ;
man.InputFileList("listshort.txt") ;
   man.ReadBranches("*") ; 
   man.SetOutput("histos.root") ;
   
   MpdCentralityAll pCentr("pCentr","pCentr") ;
   man.AddTask(&pCentr) ;
   
   MpdConvPi0 pDef("pi0Def","ConvDef") ; //name, parametes file
   man.AddTask(&pDef) ;

   MpdPairKK pKK("pKK","pKK") ;
   man.AddTask(&pKK) ;


   MpdEPAnalysisTask MpdEPAnalysisTask("EP","EP") ;
   man.AddTask(&MpdEPAnalysisTask) ;

   MpdV0AnalysisTask MpdV0AnalysisTask("V0","V0") ;
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

ofstream myfile;
myfile.open("info.txt");
myfile<< "Time and memory usage" << "\n";

  myfile << " User CPU time: " << proc.fCpuUser << "seconds, Memory: resident " << proc.fMemResident << "KB, virtual " << proc.fMemVirtual << "KB"<< "\n";
  myfile << " Total RAM: "<< memory.fMemTotal << "MB, Used RAM: " << memory.fMemUsed << " MB, Free RAM: " << memory.fMemFree << "\n"; 
  myfile << "RealTime = " << rtime << " seconds, CpuTime = " << ctime << " seconds" << "\n";

  myfile.close();
}
