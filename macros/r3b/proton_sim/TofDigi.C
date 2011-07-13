void TofDigi(){

  // Input files
   TString inFile="/lustre/land/jmargan/simulation/Ox_1p/r3bsim.root";
   TString parFile="/lustre/land/jmargan/simulation/Ox_1p/r3bpar.root";

  // Output file         
  TString outFile ="/lustre/land/jmargan/simulation/Ox_1p/histo_1p.root";


  
  // ----  Load libraries   -------------------------------------------------
  gROOT->LoadMacro("$VMCWORKDIR/gconfig/basiclibs.C");
  basiclibs();

  gSystem->Load("libGenVector");

  gSystem->Load("libGeoBase");
  gSystem->Load("libParBase");


  gSystem->Load("libBase");
  gSystem->Load("libMCStack");
  gSystem->Load("libField");
  gSystem->Load("libGen");

  //----  Load R3B specific libraries ---------------------------------------
  gSystem->Load("libR3Bbase");
  gSystem->Load("libR3BGen");
  gSystem->Load("libR3BPassive");
  gSystem->Load("libR3BData");
  gSystem->Load("libR3BDch");
  
  gSystem->Load("libR3BTof");
  gSystem->Load("libR3BmTof");
  gSystem->Load("libR3BTra");
  gSystem->Load("libR3BGfi");


  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  timer.Start();
  // ------------------------------------------------------------------------

  // -----   Digitization run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  fRun->SetInputFile(inFile);
  fRun->SetOutputFile(outFile);



  // Verbosity Mode level
  // (0=quiet, 1=event level, 2=track level, 3=debug)
  Int_t iVerbose = 1;
  //Connect the Digitization  Task

  R3BTofDigitizer* Tof  = new R3BTofDigitizer();

        

  fRun->AddTask(Tof);  

    

  // Runtime DataBase info
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo*  parIo1 = new FairParRootFileIo();
  parIo1->open(parFile.Data());
  rtdb->setFirstInput(parIo1);
  rtdb->setOutput(parIo1);
  rtdb->saveOutput();

  // Load the Root Geometry
  fRun->LoadGeometry();
  
  // Number of events to process
  Int_t nEvents = 100;
  
  // -----   Intialise and run   --------------------------------------------
  fRun->Init();
  fRun->Run(0, nEvents);

  // -----   Finish   -------------------------------------------------------
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Macro finished succesfully." << endl;
  cout << "Output file writen:  "    << outFile << endl;
  cout << "Parameter file writen " << parFile << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
  cout << endl;
  // ------------------------------------------------------------------------


}