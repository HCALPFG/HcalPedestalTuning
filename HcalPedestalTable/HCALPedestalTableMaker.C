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
    return DetName;
}

//
void HCALPedestalTableMakerSubdet(TString rootfile="../HcalNtuples_239995_Viktor.root", TString SubDet="HB", int option=2) 
{ 

    cout << "[HCAL Pedestal table maker] Running option " << option << " for " << SubDet << endl; 

    // fit pannel display option
    gStyle->SetOptFit(1011);

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
                    TH1F( Form("h1_ADC_ieta%s_iphi%i_depth%i_cap%i",(ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),(iphi+1),(idepth+1),icap),
                          Form("h1_ADC_ieta%s_iphi%i_depth%i_cap%i",(ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),(iphi+1),(idepth+1),icap),
                          128, -0.5, 127.5);
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
    int Subdet[nieta][niphi][ndepth];           // Id for subdetectors : HB=1, HE=2, HO=3, HF=4
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
        if(SubDet=="HB" || SubDet=="HE") 
        { 
            for(unsigned int i=0; i<HBHEDigiRawID_->size(); i++) 
            {
                if(SubDet=="HB" && HBHEDigiSubdet_->at(i)!=1) continue;
                if(SubDet=="HE" && HBHEDigiSubdet_->at(i)!=2) continue;
                
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
        if(SubDet=="HO") 
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
        if(SubDet=="HF") 
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
                    h1_ADC[ieta+41][iphi-1][idepth-1][HFDigiCapID_->at(i).at(icap)]->Fill(HFDigiADC_->at(i).at(icap)); 
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
                    if( h1_ADC[ieta][iphi][idepth][0]->Integral()==0 ) continue;
                    if( Subdet[ieta][iphi][idepth]==-999. ) continue; 

                    // Canvas for each channel
                    TCanvas *c = new TCanvas("c", "c", 800, 800); 
                    c->Divide(2,2);  

                    c->cd(1); h1_ADC[ieta][iphi][idepth][0]->Draw(); 
                    c->cd(2); h1_ADC[ieta][iphi][idepth][1]->Draw(); 
                    c->cd(3); h1_ADC[ieta][iphi][idepth][2]->Draw(); 
                    c->cd(4); h1_ADC[ieta][iphi][idepth][3]->Draw(); 

                    c->Print(Form("Fig/ADC_ieta%s_iphi%i_depth%i_%s_option%i.C",(ieta>=41?Form("%i",ieta-41):Form("m%i",41-ieta)),(iphi+1),(idepth+1),GetDetName(Subdet[ieta][iphi][idepth]),option)); 
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
    if(SubDet=="HB") 
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
                    setw(16) << (iphi+1)    <<
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
    
    cout << "[HCAL Pedestal table maker] Done with " << SubDet << endl;
    cout << " ........................................................................................  " << endl; 
    
}

//
void HCALPedestalTableMaker(TString rootfile="../HcalNtuples_244174_Viktor.root") 
{
    HCALPedestalTableMakerSubdet(rootfile, "HB", 2);
    HCALPedestalTableMakerSubdet(rootfile, "HE", 2);
    HCALPedestalTableMakerSubdet(rootfile, "HO", 2);
    HCALPedestalTableMakerSubdet(rootfile, "HF", 2);
}
