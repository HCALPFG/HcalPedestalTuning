// ------------------------------------------------------------------------------------
//  ROOT macro that produces a pedestal table from a PFG ntuple
//
//  Author : Jae Hyeok Yoo (jae.hyeok.yoo@cern.ch)
//  Written on 06/01/2015 
//  Last update on 09/28/2017
// ------------------------------------------------------------------------------------
//  
// Pre-requisite :
//
//   You should have the PFG ntuple for the Run from which you want to do a measurement. 
//   Instruction on how to make PFG ntuples can be found at
//   https://github.com/HCALPFG/HcalTupleMaker/tree/PFG-CMSSW_9_0_X
//
// Caution : 
//
//   If the pedestal table file already exists, and you want to generate a new one 
//   with the same name, remove the existing file or change its name because new 
//   contents will be "appeneded" to the existing file, i.e., existing file is not overwritten. 
//
// Usage : 
//
//   $ root -b  
//   root> .L HCALPedestalTableMaker.C++ 
//   root> HCALPedestalTableMaker("PFGntuple.root")
//    
// -----------------------------------------------------------------------------------
// 
// There are three options to extract the mean and the RMS of a pedestal distribution.
// Note that option 2 is the validated and the default option. 
//
//   - option 0           : Gaussian fit  
//   - option 1           : GetMean() and GetRMS() in TH1
//   - option 2 (default) : manual calculation in the range of 0 - max+6  
//

// 
// Indices of channels in the subdetectors 
//
//  HBHE -----------------------
//      IEta    = -29 - 29
//      IPhi    =  1 - 72 
//      Depth   =  1,2,3 
//  HO -------------------------
//      IEta    = -15 - 15 
//      IPhi    =  1 - 72 
//      Depth   =  4 
//  HF ------------------------
//      IEta    = -41--29, 29-41
//      IPhi    =  3,5,7,9,...,25 
//      Depth   =  1,2 

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip> // for setw()

#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TString.h"
#include "TStyle.h"

// In order to use vector of vectors : vector<vector<data type> >
// ACLiC makes dictionary for this
// [ref] http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=10236&p=44117#p44117
#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<int> >+;
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

using namespace std;

bool DRAWPLOTS    = false;   // draw plots or not (make "Fig" directory first before turning this on)
bool VERBOSE      = false;  // print out mean +/- sigma for each channel or not
bool doEffective  = true;  // Effective pedestal for Phase 1 HE 

//
// Conversion from decimal to hexadecimal 
//
TString DeciToHexa(int dec)
{   
  TString hexa = TString::Itoa(dec,16);
  hexa.ReplaceAll("a","A");
  hexa.ReplaceAll("b","B");
  hexa.ReplaceAll("c","C");
  hexa.ReplaceAll("d","D");
  hexa.ReplaceAll("e","E");
  hexa.ReplaceAll("f","F");
  return hexa; 
}

const char* GetDetName(int Subdet) 
{ 
  const char* DetName;
  if(Subdet==1) DetName = "HB"; 
  if(Subdet==2) DetName = "HE"; 
  if(Subdet==3) DetName = "HO"; 
  if(Subdet==4) DetName = "HF"; 
  if(Subdet==5) DetName = "HF"; // For QIE10, use HF 
  if(Subdet==8) DetName = "HE"; // Foe QIE11, use HE 
  return DetName;
}

bool isBadChannel(int ieta, int iphi, int depth, TString sub)
{ 
  bool isbadchannel = false;
  if(ieta==-16 && iphi==7  && depth==1 && sub=="HB") isbadchannel=true;
  if(ieta==-15 && iphi==7  && depth==1 && sub=="HB") isbadchannel=true;
  if(ieta==-13 && iphi==7  && depth==1 && sub=="HB") isbadchannel=true;
  if(ieta==-6  && iphi==2  && depth==4 && sub=="HO") isbadchannel=true;
  if(ieta==-6  && iphi==3  && depth==4 && sub=="HO") isbadchannel=true;
  if(ieta==-6  && iphi==4  && depth==4 && sub=="HO") isbadchannel=true;
  if(ieta==5   && iphi==5  && depth==4 && sub=="HO") isbadchannel=true;
  if(ieta==5   && iphi==6  && depth==4 && sub=="HO") isbadchannel=true;
  if(ieta==5   && iphi==7  && depth==4 && sub=="HO") isbadchannel=true;
  if(ieta==15  && iphi==11 && depth==4 && sub=="HO") isbadchannel=true;
  if(ieta==15  && iphi==12 && depth==4 && sub=="HO") isbadchannel=true;
  if(ieta==15  && iphi==13 && depth==4 && sub=="HO") isbadchannel=true;
  if(ieta==-29 && iphi==31 && depth==1 && sub=="HE") isbadchannel=true;
  if(ieta==-28 && iphi==31 && depth==1 && sub=="HE") isbadchannel=true;
  if(ieta==-18 && iphi==31 && depth==2 && sub=="HE") isbadchannel=true;
  if(ieta==-17 && iphi==31 && depth==2 && sub=="HE") isbadchannel=true;
  if(ieta==-16 && iphi==31 && depth==4 && sub=="HE") isbadchannel=true;
  if(ieta==-29 && iphi==31 && depth==2 && sub=="HE") isbadchannel=true;
  return isbadchannel;
}

//
void HCALPedestalTableMakerSubdet(TString rootfile="../../results.root", TString SubDetName="HB", int option=2) 
{ 
  // do not print out too many errors ... 
  gErrorIgnoreLevel=kError+1;

  cout << "[HCAL Pedestal table maker] Running option " << option << " for " << SubDetName << endl; 

  // fit pannel display option
  gStyle->SetOptFit(1011);
  gStyle->SetOptStat(111111); 

  //
  // Get the tree from the PFG ntuple 
  //
  TFile *f = TFile::Open(rootfile, "READ");
  TDirectory* dir = f->GetDirectory("hcalTupleTree");
  dir->cd();
  TTree *tree = (TTree*)dir->Get("tree");

  //
  // Set up branch address
  //

  // event,ls and run  
  UInt_t   event_ = 0;
  tree->SetBranchAddress("event", &event_);
  UInt_t   ls_ = 0;
  tree->SetBranchAddress("ls", &ls_);
  UInt_t   run_ = 0;
  tree->SetBranchAddress("run", &run_);

  // HBHE
  vector<int>   *HBHEDigiRawID_ = 0;
  tree->SetBranchAddress("HBHEDigiRawID", &HBHEDigiRawID_);
  vector<int>   *HBHEDigiSubdet_ = 0;
  tree->SetBranchAddress("HBHEDigiSubdet", &HBHEDigiSubdet_);
  vector<int>   *HBHEDigiIEta_ = 0;
  tree->SetBranchAddress("HBHEDigiIEta", &HBHEDigiIEta_);
  vector<int>   *HBHEDigiIPhi_ = 0;
  tree->SetBranchAddress("HBHEDigiIPhi", &HBHEDigiIPhi_);
  vector<int>   *HBHEDigiDepth_ = 0;
  tree->SetBranchAddress("HBHEDigiDepth", &HBHEDigiDepth_);
  vector<vector<int> >   *HBHEDigiCapID_ = 0;
  tree->SetBranchAddress("HBHEDigiCapID", &HBHEDigiCapID_);
  vector<vector<float> >   *HBHEDigiNomFC_ = 0; // linearlized ADC count
  tree->SetBranchAddress("HBHEDigiNomFC", &HBHEDigiNomFC_);
  vector<vector<float> >   *HBHEDigiADC_ = 0; // unlinearlized ADC count
  tree->SetBranchAddress("HBHEDigiADC", &HBHEDigiADC_);

  // HO  
  vector<int>   *HODigiRawID_ = 0;
  tree->SetBranchAddress("HODigiRawID", &HODigiRawID_);
  vector<int>   *HODigiSubdet_ = 0;
  tree->SetBranchAddress("HODigiSubdet", &HODigiSubdet_);
  vector<int>   *HODigiIEta_ = 0;
  tree->SetBranchAddress("HODigiIEta", &HODigiIEta_);
  vector<int>   *HODigiIPhi_ = 0;
  tree->SetBranchAddress("HODigiIPhi", &HODigiIPhi_);
  vector<int>   *HODigiDepth_ = 0;
  tree->SetBranchAddress("HODigiDepth", &HODigiDepth_);
  vector<vector<int> >   *HODigiCapID_ = 0;
  tree->SetBranchAddress("HODigiCapID", &HODigiCapID_);
  vector<vector<float> >   *HODigiNomFC_ = 0; // linearlized ADC count
  tree->SetBranchAddress("HODigiNomFC", &HODigiNomFC_);
  vector<vector<float> >   *HODigiADC_ = 0; // unlinearlized ADC count
  tree->SetBranchAddress("HODigiADC", &HODigiADC_);

  // HF 
  vector<int>   *HFDigiRawID_ = 0;
  tree->SetBranchAddress("HFDigiRawID", &HFDigiRawID_);
  vector<int>   *HFDigiSubdet_ = 0;
  tree->SetBranchAddress("HFDigiSubdet", &HFDigiSubdet_);
  vector<int>   *HFDigiIEta_ = 0;
  tree->SetBranchAddress("HFDigiIEta", &HFDigiIEta_);
  vector<int>   *HFDigiIPhi_ = 0;
  tree->SetBranchAddress("HFDigiIPhi", &HFDigiIPhi_);
  vector<int>   *HFDigiDepth_ = 0;
  tree->SetBranchAddress("HFDigiDepth", &HFDigiDepth_);
  vector<vector<int> >   *HFDigiCapID_ = 0;
  tree->SetBranchAddress("HFDigiCapID", &HFDigiCapID_);
  vector<vector<float> >   *HFDigiNomFC_ = 0; // linearlized ADC count
  tree->SetBranchAddress("HFDigiNomFC", &HFDigiNomFC_);
  vector<vector<float> >   *HFDigiADC_ = 0; // unlinearlized ADC count
  tree->SetBranchAddress("HFDigiADC", &HFDigiADC_);

  // QIE10 
  vector<int>   *QIE10DigiRawID_ = 0;
  tree->SetBranchAddress("QIE10DigiRawID", &QIE10DigiRawID_);
  vector<int>   *QIE10DigiSubdet_ = 0;
  tree->SetBranchAddress("QIE10DigiSubdet", &QIE10DigiSubdet_);
  vector<int>   *QIE10DigiIEta_ = 0;
  tree->SetBranchAddress("QIE10DigiIEta", &QIE10DigiIEta_);
  vector<int>   *QIE10DigiIPhi_ = 0;
  tree->SetBranchAddress("QIE10DigiIPhi", &QIE10DigiIPhi_);
  vector<int>   *QIE10DigiDepth_ = 0;
  tree->SetBranchAddress("QIE10DigiDepth", &QIE10DigiDepth_);
  vector<vector<int> >   *QIE10DigiCapID_ = 0;
  tree->SetBranchAddress("QIE10DigiCapID", &QIE10DigiCapID_);
  vector<vector<float> >   *QIE10DigiADC_ = 0; // unlinearlized ADC count
  tree->SetBranchAddress("QIE10DigiADC", &QIE10DigiADC_);

  // QIE11 (HEP17) 
  vector<int>   *QIE11DigiRawID_ = 0;
  tree->SetBranchAddress("QIE11DigiRawID", &QIE11DigiRawID_);
  vector<int>   *QIE11DigiSubdet_ = 0;
  tree->SetBranchAddress("QIE11DigiSubdet", &QIE11DigiSubdet_);
  vector<int>   *QIE11DigiIEta_ = 0;
  tree->SetBranchAddress("QIE11DigiIEta", &QIE11DigiIEta_);
  vector<int>   *QIE11DigiIPhi_ = 0;
  tree->SetBranchAddress("QIE11DigiIPhi", &QIE11DigiIPhi_);
  vector<int>   *QIE11DigiDepth_ = 0;
  tree->SetBranchAddress("QIE11DigiDepth", &QIE11DigiDepth_);
  vector<vector<int> >   *QIE11DigiCapID_ = 0;
  tree->SetBranchAddress("QIE11DigiCapID", &QIE11DigiCapID_);
  vector<vector<int> >   *QIE11DigiADC_ = 0; // unlinearlized ADC count
  tree->SetBranchAddress("QIE11DigiADC", &QIE11DigiADC_);
  vector<vector<float> >   *QIE11DigiFC_ = 0;
  tree->SetBranchAddress("QIE11DigiFC", &QIE11DigiFC_);

  // 
  // Define histograms for each channel 
  //  - One channel has 4 capacitors, so there are four plots per channel
  //  - Unlearized ADC count goes from 0 to 127, so there are 128 bins 
  //    and the range is from -0.5 to 127.5 
  //    

  // number of indices in eta, phi, depth
  int nieta = 83;
  int niphi = 72;
  int ndepth = 4; if(SubDetName=="QIE11") ndepth = 7; 
  TH1F *h1_ADC[nieta][niphi][ndepth][4]; // the last dimention is capid 
  for(int ieta=0; ieta<nieta; ieta++) 
  { 
    for(int iphi=0; iphi<niphi; iphi++) 
    {
      for(int idepth=0; idepth<ndepth; idepth++) 
      { 
        for(int icap=0; icap<4; icap++)  
        {
          h1_ADC[ieta][iphi][idepth][icap] = new 
            TH1F( Form("h1_ADC_ieta%s_iphi%i_depth%i_cap%i",
                  (ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),
                  iphi+1, //SubDetName=="QIE11"?iphi:iphi+1, 
                  (idepth+1),icap),
                Form("h1_ADC_ieta%s_iphi%i_depth%i_cap%i",
                  (ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),
                  iphi+1, //SubDetName=="QIE11"?iphi:iphi+1, 
                  (idepth+1),icap),
                128, -0.5, 127.5);
                //30, -0.5, 29.5);
          h1_ADC[ieta][iphi][idepth][icap]->Sumw2();
        }
      }
    }
  }

  //
  // Define and initialize arrays to be used to make a text file 
  //
  float ADC_mean[nieta][niphi][ndepth][4];    // mean of ADC count
  float ADC_sigma[nieta][niphi][ndepth][4];   // sigma of ADC count 
  int DetId[nieta][niphi][ndepth];            // Id for channel : it is decimal in the ntuple but to be converted into Heximal  
  int Subdet[nieta][niphi][ndepth];           // Id for subdetectors : HB=1, HE=2, HO=3, HF=4, QIE10=5, CRF=8
  for(int ieta=0; ieta<nieta; ieta++) 
  { 
    for(int iphi=0; iphi<niphi; iphi++) 
    {
      for(int idepth=0; idepth<ndepth; idepth++) 
      { 
        DetId[ieta][iphi][idepth] = -999.; 
        Subdet[ieta][iphi][idepth] = -999.; 
        for(int icap=0; icap<4; icap++)  
        {
          ADC_mean[ieta][iphi][idepth][icap] = -999.; 
          ADC_sigma[ieta][iphi][idepth][icap] = -999.; 
        }
      }
    }
  }

  //
  // Loop over events 
  //
  unsigned int nentries = (Int_t)tree->GetEntries();
  cout << "[HCAL Pedestal table maker] The number of entries is: " << nentries << endl;

  // main event loop
  for(unsigned int ievent = 0; ievent<nentries; ievent++) 
  //for(unsigned int ievent = 0; ievent<100; ievent++) 
  {
    tree->GetEntry(ievent); 

    // Progress indicator 
    if(ievent%100==0) cout << "[HCAL Pedestal table maker] Processed " << ievent << " out of " << nentries << " events" << endl; 

    // Fill HBHE
    if(SubDetName=="HB" || SubDetName=="HE") 
    { 
      for(unsigned int i=0; i<HBHEDigiRawID_->size(); i++) 
      {
        if(SubDetName=="HB" && HBHEDigiSubdet_->at(i)!=1) continue;
        if(SubDetName=="HE" && HBHEDigiSubdet_->at(i)!=2) continue;

        int ieta =  HBHEDigiIEta_->at(i);
        int iphi =  HBHEDigiIPhi_->at(i);
        int idepth =  HBHEDigiDepth_->at(i);

        DetId[ieta+41][iphi-1][idepth-1] = HBHEDigiRawID_->at(i);  
        Subdet[ieta+41][iphi-1][idepth-1] = HBHEDigiSubdet_->at(i);  

        for(unsigned int its=0; its<HBHEDigiCapID_->at(i).size(); its++)  
        { 
          h1_ADC[ieta+41][iphi-1][idepth-1][HBHEDigiCapID_->at(i).at(its)]->Fill(HBHEDigiADC_->at(i).at(its)); 
        }
      }
    } 

    // Fill HO
    if(SubDetName=="HO") 
    {
      for(unsigned int i=0; i<HODigiRawID_->size(); i++) 
      {
        int ieta =  HODigiIEta_->at(i);
        int iphi =  HODigiIPhi_->at(i);
        int idepth =  HODigiDepth_->at(i);

        DetId[ieta+41][iphi-1][idepth-1] = HODigiRawID_->at(i);  
        Subdet[ieta+41][iphi-1][idepth-1] = HODigiSubdet_->at(i);  

        for(unsigned int its=0; its<HODigiCapID_->at(i).size(); its++)  
        { 
          h1_ADC[ieta+41][iphi-1][idepth-1][HODigiCapID_->at(i).at(its)]->Fill(HODigiADC_->at(i).at(its)); 
        }
      }
    } 

    // Fill HF
    if(SubDetName=="HF") 
    {
      for(unsigned int i=0; i<HFDigiRawID_->size(); i++) 
      {
        int ieta =  HFDigiIEta_->at(i);
        int iphi =  HFDigiIPhi_->at(i);
        int idepth =  HFDigiDepth_->at(i);

        DetId[ieta+41][iphi-1][idepth-1] = HFDigiRawID_->at(i);  
        Subdet[ieta+41][iphi-1][idepth-1] = HFDigiSubdet_->at(i);  

        for(unsigned int its=0; its<HFDigiCapID_->at(i).size(); its++)  
        { 
          h1_ADC[ieta+41][iphi-1][idepth-1][HFDigiCapID_->at(i).at(its)]->Fill((HFDigiADC_->at(i).at(its))); 
        }
      }
    } 

    // Fill QIE10 
    if(SubDetName=="QIE10") 
    {
      for(unsigned int i=0; i<QIE10DigiRawID_->size(); i++) 
      {
        int ieta =  QIE10DigiIEta_->at(i);
        int iphi =  QIE10DigiIPhi_->at(i);
        int idepth =  QIE10DigiDepth_->at(i);

        DetId[ieta+41][iphi-1][idepth-1] = QIE10DigiRawID_->at(i);  
        Subdet[ieta+41][iphi-1][idepth-1] = QIE10DigiSubdet_->at(i);  

        for(unsigned int its=0; its<QIE10DigiCapID_->at(i).size(); its++)   // There are 3 TSs in local runs 
        { 
          h1_ADC[ieta+41][iphi-1][idepth-1][QIE10DigiCapID_->at(i).at(its)]->Fill(QIE10DigiADC_->at(i).at(its)); 
        }
      }
    } 

    // Fill QIE11 
    if(SubDetName=="QIE11") 
    {
      for(unsigned int i=0; i<QIE11DigiRawID_->size(); i++) 
      {
        int ieta =  QIE11DigiIEta_->at(i);
        int iphi =  QIE11DigiIPhi_->at(i);
        int idepth =  QIE11DigiDepth_->at(i);

        DetId[ieta+41][iphi-1][idepth-1] = QIE11DigiRawID_->at(i);  
        Subdet[ieta+41][iphi-1][idepth-1] = QIE11DigiSubdet_->at(i);  

        for(unsigned int its=0; its<QIE11DigiCapID_->at(i).size(); its++)
        { 
          h1_ADC[ieta+41][iphi-1][idepth-1][QIE11DigiCapID_->at(i).at(its)]->Fill(QIE11DigiADC_->at(i).at(its));  
        }
      }
    } 


  } //for(unsigned int ievent = 0; ievent<nentries; ievent++) 

  // 
  // Extract mean and sigma 
  // 
  cout << endl; 
  cout << " ........................................................................................  " << endl; 
  cout << " ........................... Extraction of mean and sigma ...............................  " << endl; 
  cout << " ........................................................................................  " << endl; 
  cout << endl; 

  for(int ieta=0; ieta<nieta; ieta++) 
  { 
    for(int iphi=0; iphi<niphi; iphi++) 
    {
      for(int idepth=0; idepth<ndepth; idepth++) 
      { 
        if( h1_ADC[ieta][iphi][idepth][0]->Integral()==0 ) continue;
        if( Subdet[ieta][iphi][idepth]==-999. ) continue; 

        for(int icap=0; icap<4; icap++)  
        { 
          if(VERBOSE) 
          { 
            cout << "[HCAL Pedestal table maker] For ieta, iphi, depth, icap = ";
            cout << (ieta-41) <<  ", " << (iphi+1) << ", " << (idepth+1) << ", " << icap << endl;
            cout << "[HCAL Pedestal table maker]   pedestal = " << h1_ADC[ieta][iphi][idepth][icap]->GetMean() << " +/- " 
              << h1_ADC[ieta][iphi][idepth][icap]->GetRMS() << endl;  
          } 

          // Gaussian fit  
          if(option==0) 
          { 
            float FitRangeBegin = h1_ADC[ieta][iphi][idepth][icap]->GetMean()-3*h1_ADC[ieta][iphi][idepth][icap]->GetRMS(); 
            float FitRangeEnd   = h1_ADC[ieta][iphi][idepth][icap]->GetMean()+3*h1_ADC[ieta][iphi][idepth][icap]->GetRMS(); 
            TF1 *gfit = new TF1("gfit","gaus",FitRangeBegin,FitRangeEnd); 
            h1_ADC[ieta][iphi][idepth][icap]->Fit("gfit","R");  
            ADC_mean[ieta][iphi][idepth][icap]   = gfit->GetParameter(1);
            ADC_sigma[ieta][iphi][idepth][icap]  = gfit->GetParameter(2); 
            //delete gfit; 
          } 
          // No fit : TH1 GetMean() and GetRMS()  
          else if(option==1) 
          {
            ADC_mean[ieta][iphi][idepth][icap] = h1_ADC[ieta][iphi][idepth][icap]->GetMean();   
            ADC_sigma[ieta][iphi][idepth][icap] = h1_ADC[ieta][iphi][idepth][icap]->GetRMS(); 
          }
          // !! THIS IS THE DRFAULT METHOD !!
          // No fit : manual calculation of mean and RMS
          // Note that range is bin 1 to bin max  
          // where max is bin peak + 3  
          else if(option==2) 
          {
            // Get the integration range 
            int from,to,max=h1_ADC[ieta][iphi][idepth][icap]->GetBinContent(1),maxi=1;
            for(int i=1;i<=128;i++)
            { 
              if(h1_ADC[ieta][iphi][idepth][icap]->GetBinContent(i)>max)
              { 
                max=h1_ADC[ieta][iphi][idepth][icap]->GetBinContent(i); 
                maxi=i;
              }
            } 
            from=1;  
            if(Subdet[ieta][iphi][idepth]==8 && doEffective) to=256; // For effective PED for QIE11 (Phase1 HE)  
            else to=maxi+3; 
            if(Subdet[ieta][iphi][idepth]==8 && to>256) to=256;  
            else if(to>128) to=128; 

            // peak at ADC=0
            if(maxi==1)
            { 
              if(Subdet[ieta][iphi][idepth]==8) to=256;
              else to=128;
            }

            // Calculate mean and width 
            double Sum=0,nSum=0;
            for(int i=from;i<=to;i++)
            {
              Sum+=(i-1)*h1_ADC[ieta][iphi][idepth][icap]->GetBinContent(i);
              nSum+=h1_ADC[ieta][iphi][idepth][icap]->GetBinContent(i);
            }
            ADC_mean[ieta][iphi][idepth][icap]=Sum/nSum; 
            Sum=0;
            for(int i=from;i<=to;i++) Sum+=h1_ADC[ieta][iphi][idepth][icap]->GetBinContent(i)*((i-1)-ADC_mean[ieta][iphi][idepth][icap])*((i-1)-ADC_mean[ieta][iphi][idepth][icap]);
            ADC_sigma[ieta][iphi][idepth][icap]=TMath::Sqrt(Sum/nSum); 
          }
        }
      }
    }
  }

  // 
  // Drawing : pedestal distribution per channel 
  // 
  // 
  if(DRAWPLOTS)
  {
    cout << endl; 
    cout << " ........................................................................................  " << endl; 
    cout << " ..................................... Drawing ..........................................  " << endl; 
    cout << " ........................................................................................  " << endl; 
    cout << endl; 

    for(int ieta=0; ieta<nieta; ieta++) 
    { 
      for(int iphi=0; iphi<niphi; iphi++) 
      {
        for(int idepth=0; idepth<ndepth; idepth++) 
        {  
          // Draw plots if the mean is too large or too small 
          // You can define other criteria such as large or small RMS  
          bool DRAWCHANNEL=false;
          float highADC, lowADC;
          float highADCsigma, lowADCsigma;
          if(Subdet[ieta][iphi][idepth]==3) { highADC=11; lowADC=7;}                    
          else { highADC=6; lowADC=1.5;}                    
          if(ADC_mean[ieta][iphi][idepth][0]>highADC || ADC_mean[ieta][iphi][idepth][0]<lowADC) DRAWCHANNEL=true;
          if(ADC_mean[ieta][iphi][idepth][1]>highADC || ADC_mean[ieta][iphi][idepth][1]<lowADC) DRAWCHANNEL=true;
          if(ADC_mean[ieta][iphi][idepth][2]>highADC || ADC_mean[ieta][iphi][idepth][2]<lowADC) DRAWCHANNEL=true;
          if(ADC_mean[ieta][iphi][idepth][3]>highADC || ADC_mean[ieta][iphi][idepth][3]<lowADC) DRAWCHANNEL=true;
          if(Subdet[ieta][iphi][idepth]==3) { highADCsigma=1.2; lowADCsigma=0.4;}                     
          else { highADCsigma=2; lowADCsigma=0.2;}                     
          if(ADC_sigma[ieta][iphi][idepth][0]>highADCsigma || ADC_sigma[ieta][iphi][idepth][0]<lowADCsigma) DRAWCHANNEL=true;
          if(ADC_sigma[ieta][iphi][idepth][1]>highADCsigma || ADC_sigma[ieta][iphi][idepth][1]<lowADCsigma) DRAWCHANNEL=true;
          if(ADC_sigma[ieta][iphi][idepth][2]>highADCsigma || ADC_sigma[ieta][iphi][idepth][2]<lowADCsigma) DRAWCHANNEL=true;
          if(ADC_sigma[ieta][iphi][idepth][3]>highADCsigma || ADC_sigma[ieta][iphi][idepth][3]<lowADCsigma) DRAWCHANNEL=true;
/*          
          for(int i=0; i<4; i++)
          {
            if((int)h1_ADC[ieta][iphi][idepth][i]->Integral()%2000!=0) DRAWCHANNEL=true;
          } 
   
          if(((ieta-41)==32 && (iphi+1)==15 && (idepth+1)==3) ||
             ((ieta-41)==-33 && (iphi+1)==71 && (idepth+1)==2) ||
             ((ieta-41)==33 && (iphi+1)==45 && (idepth+1)==2) ||
             ((ieta-41)==-31 && (iphi+1)==35 && (idepth+1)==4) ||
             ((ieta-41)==31 && (iphi+1)==15 && (idepth+1)==3)
            ) DRAWCHANNEL=true;
*/
          if(!DRAWCHANNEL) continue;

          // Skip if the histogram is empty (= channel does not exis in the ntuple) 
          if( h1_ADC[ieta][iphi][idepth][0]->Integral()==0 ) continue;
          if( Subdet[ieta][iphi][idepth]==-999. ) continue; 

          // Canvas for each channel
          TCanvas *c = new TCanvas("c", "c", 800, 800); 
          c->Divide(2,2);  

          c->cd(1); c->cd(1)->SetLogy(1); h1_ADC[ieta][iphi][idepth][0]->Draw("hist e"); 
          c->cd(2); c->cd(2)->SetLogy(1); h1_ADC[ieta][iphi][idepth][1]->Draw("hist e"); 
          c->cd(3); c->cd(3)->SetLogy(1); h1_ADC[ieta][iphi][idepth][2]->Draw("hist e"); 
          c->cd(4); c->cd(4)->SetLogy(1); h1_ADC[ieta][iphi][idepth][3]->Draw("hist e"); 

          c->Print(Form("Fig/ADC_ieta%s_iphi%i_depth%i_%s_option%i.pdf",(ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),(iphi+1),(idepth+1),GetDetName(Subdet[ieta][iphi][idepth]),option)); 
          delete c;
        }
      }
    }
  } 

  // 
  // Print pedestal table file 
  // 
  cout << endl; 
  cout << " ........................................................................................  " << endl; 
  cout << " .......................... Printing pedestal table .....................................  " << endl; 
  cout << " ........................................................................................  " << endl; 
  cout << endl; 

  // File name depending on option
  TString PedestalTable = rootfile;
  PedestalTable.ReplaceAll(".root",".txt");
  PedestalTable.ReplaceAll("../","");
  if(option==0) PedestalTable = "PedestalTable_option0_"+PedestalTable;
  if(option==1) PedestalTable = "PedestalTable_option1_"+PedestalTable;
  if(option==2) PedestalTable = "PedestalTable_"+PedestalTable;

  cout << "[HCAL Pedestal table maker] Printing pedestal table file : " << PedestalTable.Data() << endl;

  // Open file 
  ofstream fout(PedestalTable.Data(), ios_base::app | ios_base::out);

  // Printing header
  if(SubDetName=="HB")
  {   
    fout << "#U ADC  << this is the unit" << endl;
    fout <<
      setw(1) <<  "#"   <<
      setw(16) << "eta" <<
      setw(16) << "phi" <<
      setw(16) << "dep" <<
      setw(16) << "det" << 
      setw(9) << "cap0" << 
      setw(9) << "cap1" << 
      setw(9) << "cap2" << 
      setw(9) << "cap3" << 
      setw(10) << "widthcap0" << 
      setw(10) << "widthcap1" << 
      setw(10) << "widthcap2" << 
      setw(10) << "widthcap3" << 
      setw(11) << "DetId"  
      << endl; 
  }

  // Printing table 
  for(int iphi=0; iphi<niphi; iphi++) 
  { 
    for(int ieta=0; ieta<nieta; ieta++) 
    {
      for(int idepth=0; idepth<ndepth; idepth++) 
      { 
        if(ADC_mean[ieta][iphi][idepth][0] == -999 ) continue;
        if(ADC_mean[ieta][iphi][idepth][1] == -999 ) continue;
        if(ADC_mean[ieta][iphi][idepth][2] == -999 ) continue;
        if(ADC_mean[ieta][iphi][idepth][3] == -999 ) continue;

        if(isBadChannel(ieta-41, iphi+1, idepth+1, SubDetName)) continue;

        fout <<
          setw(17) << (ieta-41)   <<
          setw(16) << iphi+1 <<
          setw(16) << (idepth+1)  <<
          setw(16) << GetDetName(Subdet[ieta][iphi][idepth])   << 
          setw(9) << Form("%.5f", ADC_mean[ieta][iphi][idepth][0]) << 
          setw(9) << Form("%.5f", ADC_mean[ieta][iphi][idepth][1]) << 
          setw(9) << Form("%.5f", ADC_mean[ieta][iphi][idepth][2]) << 
          setw(9) << Form("%.5f", ADC_mean[ieta][iphi][idepth][3]) << 
          setw(9) << Form("%.5f", ADC_sigma[ieta][iphi][idepth][0]) << 
          setw(9) << Form("%.5f", ADC_sigma[ieta][iphi][idepth][1]) << 
          setw(9) << Form("%.5f", ADC_sigma[ieta][iphi][idepth][2]) << 
          setw(9) << Form("%.5f", ADC_sigma[ieta][iphi][idepth][3]) << 
          setw(11) << DeciToHexa(DetId[ieta][iphi][idepth])
          << endl; 
      }
    }
  }

  fout.close();
  f->Close();

  cout << "[HCAL Pedestal table maker] Done with " << SubDetName << endl;
  cout << " ........................................................................................  " << endl; 

}




void HCALPedestalTableMakerZDC(TString rootfile="../../results.root", int option=2)
{ 
  // File name depending on option
  TString PedestalTable = rootfile;
  PedestalTable.ReplaceAll(".root",".txt");
  PedestalTable.ReplaceAll("../","");
  if(option==0) PedestalTable = "PedestalTable_option0_"+PedestalTable;
  if(option==1) PedestalTable = "PedestalTable_option1_"+PedestalTable;
  if(option==2) PedestalTable = "PedestalTable_"+PedestalTable;

  cout << "[HCAL Pedestal table maker] Printing pedestal table file : " << PedestalTable.Data() << endl;
  cout << "[HCAL Pedestal table maker] Adding ZDC channels ... " << endl;

  // Open file
  ofstream fout(PedestalTable.Data(), ios_base::app | ios_base::out);
  // Write ZDC channels 
  fout << "               -1               1             -99          ZDC_EM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000011" << endl;
  fout << "               -1               2             -99          ZDC_EM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000012" << endl;
  fout << "               -1               3             -99          ZDC_EM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000013" << endl;
  fout << "               -1               4             -99          ZDC_EM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000014" << endl;
  fout << "               -1               5             -99          ZDC_EM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000015" << endl;
  fout << "               -1               1             -99         ZDC_HAD  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000021" << endl;
  fout << "               -1               2             -99         ZDC_HAD  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000022" << endl;
  fout << "               -1               3             -99         ZDC_HAD  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000023" << endl;
  fout << "               -1               4             -99         ZDC_HAD  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000024" << endl;
  fout << "               -1               1             -99         ZDC_LUM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000031" << endl;
  fout << "               -1               2             -99         ZDC_LUM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000032" << endl;
  fout << "                1               1             -99          ZDC_EM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000051" << endl;
  fout << "                1               2             -99          ZDC_EM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000052" << endl;
  fout << "                1               3             -99          ZDC_EM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000053" << endl;
  fout << "                1               4             -99          ZDC_EM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000054" << endl;
  fout << "                1               5             -99          ZDC_EM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000055" << endl;
  fout << "                1               1             -99         ZDC_HAD  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000061" << endl;
  fout << "                1               2             -99         ZDC_HAD  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000062" << endl;
  fout << "                1               3             -99         ZDC_HAD  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000063" << endl;
  fout << "                1               4             -99         ZDC_HAD  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000064" << endl;
  fout << "                1               1             -99         ZDC_LUM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000071" << endl;
  fout << "                1               2             -99         ZDC_LUM  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000   54000072" << endl;

  fout.close();
}

void HCALPedestalTableMakerMissingCh(TString rootfile="../../results.root", int option=2)
{ 
  // File name depending on option
  TString PedestalTable = rootfile;
  PedestalTable.ReplaceAll(".root",".txt");
  PedestalTable.ReplaceAll("../","");
  if(option==0) PedestalTable = "PedestalTable_option0_"+PedestalTable;
  if(option==1) PedestalTable = "PedestalTable_option1_"+PedestalTable;
  if(option==2) PedestalTable = "PedestalTable_"+PedestalTable;

  cout << "[HCAL Pedestal table maker] Printing pedestal table file : " << PedestalTable.Data() << endl;
  cout << "[HCAL Pedestal table maker] Adding missing channels ... " << endl;

  // Open file
  ofstream fout(PedestalTable.Data(), ios_base::app | ios_base::out);
  // Write missing channels 
  fout << "              -16               7               1              HB  3.08350  3.13850  3.22650  3.49300  0.66184  0.69233  0.68935  0.70954   43104007" << endl;
  fout << "              -15               7               1              HB  3.31875  2.54675  2.87550  3.03000  0.73800  0.72754  0.72007  0.74101   43103C07" << endl;
  fout << "              -13               7               1              HB  3.44300  3.05125  3.22975  3.35300  0.68392  0.66115  0.67599  0.67037   43103407" << endl;
  fout << "               -6               2               4              HO  8.83071  9.05575  9.10280  8.99975  0.57068  0.55601  0.53882  0.53583   47401802" << endl;
  fout << "               -6               3               4              HO  9.33083  9.06056  9.19260  9.27107  0.62614  0.58325  0.58617  0.59677   47401803" << endl;
  fout << "               -6               4               4              HO  9.23181  9.17797  8.90468  9.18409  0.67241  0.68169  0.63677  0.60399   47401804" << endl;
  fout << "                5               5               4              HO  8.64291  9.21241  8.95773  8.85671  0.70366  0.65954  0.62941  0.68472   47481405" << endl;
  fout << "                5               6               4              HO  9.32908  9.36909  8.91544  9.01526  0.73356  0.69313  0.68389  0.63468   47481406" << endl;
  fout << "                5               7               4              HO  9.14429  9.43118  9.51014  8.93047  0.69466  0.66893  0.68496  0.66135   47481407" << endl;
  fout << "               15              11               4              HO  9.01124  8.80490  8.70990  8.75854  0.77310  0.77459  0.77766  0.77456   47483C0B" << endl;
  fout << "               15              12               4              HO  9.17137  8.85672  8.65863  8.56180  0.72147  0.80455  0.80457  0.70899   47483C0C" << endl;
  fout << "               15              13               4              HO  8.79438  8.69312  8.69490  8.68659  0.76303  0.77146  0.79822  0.78791   47483C0D" << endl;
  fout << "              -29              31               1              HE  4.68019  4.74444  4.70787  4.38100  0.51177  0.47945  0.49223  0.51134   4510741F" << endl;
  fout << "              -28              31               1              HE  4.73106  4.47875  4.32544  4.37388  0.49192  0.52635  0.49158  0.52342   4510701F" << endl;
  fout << "              -18              31               2              HE  4.53075  4.81050  4.28950  4.66463  0.52386  0.43542  0.48212  0.50623   4520481F" << endl;
  fout << "              -17              31               2              HE  4.59381  4.78719  4.80419  4.75875  0.51522  0.44331  0.43081  0.47135   4520441F" << endl;
  fout << "              -16              31               4              HE  4.29944  4.70856  4.78312  4.79725  0.48324  0.49232  0.44913  0.43848   4540401F" << endl;
  fout << "              -29              31               2              HE  4.61338  4.78094  4.77037  4.26088  0.51614  0.45450  0.45951  0.46497   4520741F" << endl;
  fout.close();
}



//
void HCALPedestalTableMaker(TString rootfile="../../HcalTupleMaker_ped_321028.root") 
{
  HCALPedestalTableMakerSubdet(rootfile, "HB", 2);
  HCALPedestalTableMakerSubdet(rootfile, "HO", 2);
  HCALPedestalTableMakerSubdet(rootfile, "QIE10", 2); // HF
  HCALPedestalTableMakerSubdet(rootfile, "QIE11", 2); // HE
  HCALPedestalTableMakerZDC(rootfile, 2);
  HCALPedestalTableMakerMissingCh(rootfile, 2);
}
