// ------------------------------------------------------------------------------------
//  ROOT macro that produces a pedestal table from a PFG ntuple
//
//  Author : Jae Hyeok Yoo (jae.hyeok.yoo@cern.ch)
//  Written on 01 June 2015 
// ------------------------------------------------------------------------------------
//  
// Pre-requisite :
//
//   You should have the PFG ntuple for the Run from which you want to do a measurement. 
//   Instruction on how to make PFG ntuples can be found here : FIXME link here 
//
// Caution : 
//
//   If the pedestal table file already exists, and you want to generate a new one 
//   with the same name, remove the existing file or change its name because new 
//   contents will be "appeneded" to the existing file, i.e., existing file not overwritten. 
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
//   - option 2 (default) : manual calculation 
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

bool DRAWPLOTS  = false;  // draw plots or not (make "Fig" directory first before turning this on)
bool VERBOSE    = false;  // print out mean +/- sigma for each channel or not

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
    if(Subdet==8) DetName = "CRF"; 
    return DetName;
}

//
void HCALPedestalTableMakerSubdet(TString rootfile="../../results.root", TString SubDetName="HB", int option=2) 
{ 

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

    // CRF (QIE11) 
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
    vector<vector<float> >   *QIE11DigiADC_ = 0; // unlinearlized ADC count
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
    int ndepth = 4;
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
                                SubDetName=="QIE11"?iphi:iphi+1,
                                (idepth+1),icap),
                          Form("h1_ADC_ieta%s_iphi%i_depth%i_cap%i",
                                (ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),
                                SubDetName=="QIE11"?iphi:iphi+1,
                                (idepth+1),icap),
                          //128, -0.5, 127.5);
                          //16, -0.5, 15.5);
                          100, -0.5, 200);
                    h1_ADC[ieta][iphi][idepth][icap]->Sumw2();
                }
            }
        }
    }
    
    //
    // Define and initialize arrays to be used to make text file 
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
    // Loop over entries
    //
    unsigned int nentries = (Int_t)tree->GetEntries();
    cout << "[HCAL Pedestal table maker] The number of entries is: " << nentries << endl;

    // main event loop
    for(unsigned int ievent = 0; ievent<nentries; ievent++) 
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

                for(unsigned int icap=0; icap<8; icap++)  
                { 
                    h1_ADC[ieta+41][iphi-1][idepth-1][HBHEDigiCapID_->at(i).at(icap)]->Fill(HBHEDigiADC_->at(i).at(icap)); 
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

                for(unsigned int icap=0; icap<8; icap++)  
                { 
                    h1_ADC[ieta+41][iphi-1][idepth-1][HODigiCapID_->at(i).at(icap)]->Fill(HODigiADC_->at(i).at(icap)); 
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

                for(unsigned int icap=0; icap<8; icap++)  
                { 
                    h1_ADC[ieta+41][iphi-1][idepth-1][HFDigiCapID_->at(i).at(icap)]->Fill((HFDigiADC_->at(i).at(icap))); 
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

                for(unsigned int icap=0; icap<8; icap++)  
                { 
                    h1_ADC[ieta+41][iphi-1][idepth-1][QIE10DigiCapID_->at(i).at(icap)]->Fill(QIE10DigiADC_->at(i).at(icap)); 
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

                for(unsigned int icap=0; icap<8; icap++)  
                { 
                    //h1_ADC[ieta+41][iphi-1][idepth-1][QIE11DigiCapID_->at(i).at(icap)]->Fill(QIE11DigiADC_->at(i).at(icap)); 
                    h1_ADC[ieta+41][iphi-1][idepth-1][QIE11DigiCapID_->at(i).at(icap)]->Fill(QIE11DigiFC_->at(i).at(icap)); 
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
                    // No fit : manual calculation of mean and RMS  
                    // Snippet came from DQM code : link 
                    // THIS IS THE DRFAULT METHOD
                    else if(option==2) 
                    {
                        double Sum=0,nSum=0;
                        int from,to,max=h1_ADC[ieta][iphi][idepth][icap]->GetBinContent(1),maxi=0;
                        for(int i=1;i<=128;i++)
                        { 
                            if(h1_ADC[ieta][iphi][idepth][icap]->GetBinContent(i)>max)
                            { 
                                max=h1_ADC[ieta][iphi][idepth][icap]->GetBinContent(i); 
                                maxi=i;
                            }
                        }
                        from=1; 
                        to=maxi+6; 
                        if(to>128) to=128;
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
                    bool DRAWCHANNEL=false;
                    int highADC, lowADC;
                    if(idepth<3) { highADC=5; lowADC=1; }                    
                    else { highADC=11; lowADC=7; }                    
                    if(ADC_mean[ieta][iphi][idepth][0]>highADC || ADC_mean[ieta][iphi][idepth][0]<lowADC) DRAWCHANNEL=true;
                    if(ADC_mean[ieta][iphi][idepth][1]>highADC || ADC_mean[ieta][iphi][idepth][1]<lowADC) DRAWCHANNEL=true;
                    if(ADC_mean[ieta][iphi][idepth][2]>highADC || ADC_mean[ieta][iphi][idepth][2]<lowADC) DRAWCHANNEL=true;
                    if(ADC_mean[ieta][iphi][idepth][3]>highADC || ADC_mean[ieta][iphi][idepth][3]<lowADC) DRAWCHANNEL=true;
                    if(ADC_sigma[ieta][iphi][idepth][0]>1.5 || ADC_sigma[ieta][iphi][idepth][0]<0.5) DRAWCHANNEL=true;
                    if(ADC_sigma[ieta][iphi][idepth][1]>1.5 || ADC_sigma[ieta][iphi][idepth][1]<0.5) DRAWCHANNEL=true;
                    if(ADC_sigma[ieta][iphi][idepth][2]>1.5 || ADC_sigma[ieta][iphi][idepth][2]<0.5) DRAWCHANNEL=true;
                    if(ADC_sigma[ieta][iphi][idepth][3]>1.5 || ADC_sigma[ieta][iphi][idepth][3]<0.5) DRAWCHANNEL=true;
                    //if(!DRAWCHANNEL) continue;

                    // ieta,iphi= 35+41,47
                    // ieta,iphi= 34+41,39
                    // ieta,iphi= 30+41,39
                    // ieta,iphi= -18+41,18
                    //if(!((ieta==76 && iphi==46) || (ieta==75 && iphi==38) || (ieta==71 && iphi==38) || (ieta==23 && iphi==17) 
                    //      || (ieta==76 && iphi==48) || (ieta==76 && iphi==44))) continue;

                    if( h1_ADC[ieta][iphi][idepth][0]->Integral()==0 ) continue;
                    if( Subdet[ieta][iphi][idepth]==-999. ) continue; 
                   
                    if(SubDetName=="QIE11") 
                    {
                        for(int icap=0; icap<4; icap++)
                        {
                            h1_ADC[ieta][iphi][idepth][icap]->SetTitle(
                                    Form("h1_fC_ieta%s_iphi%i_depth%i_cap%i",
                                        (ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),
                                        SubDetName=="QIE11"?iphi:iphi+1,
                                        8,icap));
                            h1_ADC[ieta][iphi][idepth][icap]->SetName(
                                    Form("h1_fC_ieta%s_iphi%i_depth%i_cap%i",
                                        (ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),
                                        SubDetName=="QIE11"?iphi:iphi+1,
                                        8,icap));
                        }
                    }

                    // Canvas for each channel
                    TCanvas *c = new TCanvas("c", "c", 800, 800); 
                    c->Divide(2,2);  

                    c->cd(1); c->cd(1)->SetLogy(1); h1_ADC[ieta][iphi][idepth][0]->Draw("hist e"); 
                    c->cd(2); c->cd(2)->SetLogy(1); h1_ADC[ieta][iphi][idepth][1]->Draw("hist e"); 
                    c->cd(3); c->cd(3)->SetLogy(1); h1_ADC[ieta][iphi][idepth][2]->Draw("hist e"); 
                    c->cd(4); c->cd(4)->SetLogy(1); h1_ADC[ieta][iphi][idepth][3]->Draw("hist e"); 

                    //c->Print(Form("Fig/ADC_ieta%s_iphi%i_depth%i_%s_option%i.C",(ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),(iphi+1),(idepth+1),GetDetName(Subdet[ieta][iphi][idepth]),option)); 
                    c->Print(Form("Fig/%s_ieta%s_iphi%i_depth%i_%s_option%i.pdf",SubDetName=="QIE11"?"fC":"ADC",(ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),(SubDetName=="QIE11"?iphi:iphi+1),(idepth+1),GetDetName(Subdet[ieta][iphi][idepth]),option)); 
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
    if(SubDetName=="HB" || SubDetName=="QIE10" || SubDetName=="QIE11") 
    {   
        if(SubDetName=="QIE11") fout << "#U fC  << this is the unit" << endl;
        else fout << "#U ADC  << this is the unit" << endl;
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

    // printing table 
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

                fout <<
                    setw(17) << (ieta-41)   <<
                    setw(16) << (SubDetName=="QIE11"?iphi:iphi+1) <<
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

//
//void HCALPedestalTableMaker(TString rootfile="../../outputFile_run279721.root") 
void HCALPedestalTableMaker(TString rootfile="outputFile_USC_282314.root") 
{
//    HCALPedestalTableMakerSubdet(rootfile, "HB", 2);
//    HCALPedestalTableMakerSubdet(rootfile, "HE", 2);
//    HCALPedestalTableMakerSubdet(rootfile, "HO", 2);
//    HCALPedestalTableMakerSubdet(rootfile, "HF", 2);
//    HCALPedestalTableMakerSubdet(rootfile, "QIE10", 2);
    HCALPedestalTableMakerSubdet(rootfile, "QIE11", 1);
}
