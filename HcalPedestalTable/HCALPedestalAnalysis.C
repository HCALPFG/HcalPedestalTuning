// ------------------------------------------------------------------------------------
//  ROOT macro that analyzes 
//      - a pedestal table
//      - comparison between two pedestal tables 
//
//  Author : Jae Hyeok Yoo (jae.hyeok.yoo@cern.ch)
//  Written on 04/28/2016
//  Last update on 09/28/2017
//
// The previous pedestal tables can be found at 
//    https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalPedestalsTags2011
// ------------------------------------------------------------------------------------



#include <iostream>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TColor.h"
#include "TLatex.h"

//
// h2 cosmetics
//
void h2cosmetic(TH2D* &h2, char* title, TString Xvar="", TString Yvar="", TString Zvar="Events/bin", float Max=4, float Min=2)
{
    h2->SetTitle(title);
    h2->SetXTitle(Xvar);
    h2->SetYTitle(Yvar);
    h2->SetZTitle(Zvar);
    h2->SetMaximum(Max);
    h2->SetMinimum(Min);
    h2->SetStats(0);
}


void HCALPedestalTableAnalysis(TString PedTable="PedestalTable_results.txt") 
{ 

    // [HB HE HO HF QIE11][cap 0 1 2 3] 
    TH1D        *h1_mean[5][4];
    TH1D        *h1_width[5][4];
    for(int isub=0; isub<5; isub++)
    {
        TString subdet="HB"; if(isub==1) subdet="HE"; if(isub==2) subdet="HO"; if(isub==3) subdet="HF"; if(isub==4) subdet="QIE11";
        for(int cap=0; cap<4; cap++)
        {
            //if(isub==3) { Max=11; Min=7; }
            h1_mean[isub][cap]  = new TH1D(Form("h1_mean_%s_cap%i",subdet.Data(),cap),   Form("h1_mean_%s_cap%i",subdet.Data(),cap),      50,0,16);
            h1_width[isub][cap]  = new TH1D(Form("h1_width_%s_cap%i",subdet.Data(),cap),   Form("h1_width_%s_cap%i",subdet.Data(),cap),   50,0,2);
        }
    }
    // [depth 0-6][cap 0 1 2 3] 
    TH2D        *h2_mean[7][4];
    TH2D        *h2_width[7][4];
    for(int depth=0; depth<7; depth++)
    {
        for(int cap=0; cap<4; cap++)
        {
            float Max=10;
            float Min=0; 
            h2_mean[depth][cap]  = new TH2D(Form("h2_mean_depth%i_cap%i",(depth+1),cap),   Form("h2_mean_depth%i_cap%i",(depth+1),cap),      83,-41.5,41.5,72,0.5,72.5);
            h2_width[depth][cap] = new TH2D(Form("h2_width_depth%i_cap%i",(depth+1),cap),   Form("h2_width_depth%i_cap%i",(depth+1),cap),      83,-41.5,41.5,72,0.5,72.5);
            h2cosmetic(h2_mean[depth][cap],Form("Mean of pedestal in depth=%i and capid=%i",(depth+1),cap),"ieta", "iphi", "mean ADC", Max, Min);
            h2cosmetic(h2_width[depth][cap],Form("Width of pedestal in depth=%i and capid=%i",(depth+1),cap),"ieta", "iphi", "width ADC", 2, 0);
        }
    }

    string line;
    int ieta, iphi, depth;
    double cap0, cap1, cap2, cap3, widthcap0, widthcap1, widthcap2, widthcap3;
    string det;
    ifstream fin(PedTable.Data()); 
    if(fin.is_open()) { 
        while(fin.good()){

            // get a line from fin
            getline(fin, line);

            if( line.find("#")!=string::npos ) continue;

            // Store each element in the line to the defined variables
            stringstream stream(line);
            stream >> ieta >> iphi >> depth >> det >> cap0 >> cap1 >> cap2 >> cap3 >> widthcap0 >> widthcap1 >> widthcap2 >> widthcap3;  
            
            if( !fin.good() ) continue;
            
            if(det=="HB") { h1_mean[0][0]->Fill(cap0); h1_mean[0][1]->Fill(cap1); h1_mean[0][2]->Fill(cap2); h1_mean[0][3]->Fill(cap3);}
            if(det=="HE") { h1_mean[1][0]->Fill(cap0); h1_mean[1][1]->Fill(cap1); h1_mean[1][2]->Fill(cap2); h1_mean[1][3]->Fill(cap3);}
            if(det=="HO") { h1_mean[2][0]->Fill(cap0); h1_mean[2][1]->Fill(cap1); h1_mean[2][2]->Fill(cap2); h1_mean[2][3]->Fill(cap3);}
            if(det=="HF") { h1_mean[3][0]->Fill(cap0); h1_mean[3][1]->Fill(cap1); h1_mean[3][2]->Fill(cap2); h1_mean[3][3]->Fill(cap3);}
            if(det=="QIE11") { h1_mean[4][0]->Fill(cap0); h1_mean[4][1]->Fill(cap1); h1_mean[4][2]->Fill(cap2); h1_mean[4][3]->Fill(cap3);}
            
            if(det=="HB") { h1_width[0][0]->Fill(widthcap0); h1_width[0][1]->Fill(widthcap1); h1_width[0][2]->Fill(widthcap2); h1_width[0][3]->Fill(widthcap3);}
            if(det=="HE") { h1_width[1][0]->Fill(widthcap0); h1_width[1][1]->Fill(widthcap1); h1_width[1][2]->Fill(widthcap2); h1_width[1][3]->Fill(widthcap3);}
            if(det=="HO") { h1_width[2][0]->Fill(widthcap0); h1_width[2][1]->Fill(widthcap1); h1_width[2][2]->Fill(widthcap2); h1_width[2][3]->Fill(widthcap3);}
            if(det=="HF") { h1_width[3][0]->Fill(widthcap0); h1_width[3][1]->Fill(widthcap1); h1_width[3][2]->Fill(widthcap2); h1_width[3][3]->Fill(widthcap3);}
            if(det=="QIE11") { h1_width[4][0]->Fill(widthcap0); h1_width[4][1]->Fill(widthcap1); h1_width[4][2]->Fill(widthcap2); h1_width[4][3]->Fill(widthcap3);}

            h2_mean[depth-1][0]->SetBinContent(ieta+42,iphi,cap0);
            h2_mean[depth-1][1]->SetBinContent(ieta+42,iphi,cap1);
            h2_mean[depth-1][2]->SetBinContent(ieta+42,iphi,cap2);
            h2_mean[depth-1][3]->SetBinContent(ieta+42,iphi,cap3);
            h2_width[depth-1][0]->SetBinContent(ieta+42,iphi,widthcap0);
            h2_width[depth-1][1]->SetBinContent(ieta+42,iphi,widthcap1);
            h2_width[depth-1][2]->SetBinContent(ieta+42,iphi,widthcap2);
            h2_width[depth-1][3]->SetBinContent(ieta+42,iphi,widthcap3);

        }
    }
   

    TCanvas *c_mean_2d[7], *c_width_2d[7], *c_mean_1d[5], *c_width_1d[5];
    // 2d 
    for(int depth=0; depth<7; depth++)
    { 
        c_mean_2d[depth] = new TCanvas(Form("c_mean_2d_depth%i", depth+1), Form("c_mean_2d_depth%i", depth+1), 800,600);
        c_mean_2d[depth]->Divide(2,2);
        c_mean_2d[depth]->cd(1); h2_mean[depth][0]->Draw("colz");
        c_mean_2d[depth]->cd(2); h2_mean[depth][1]->Draw("colz");
        c_mean_2d[depth]->cd(3); h2_mean[depth][2]->Draw("colz");
        c_mean_2d[depth]->cd(4); h2_mean[depth][3]->Draw("colz");
        c_mean_2d[depth]->Print(Form("Fig/c_mean_2d_depth%i.pdf", depth+1));

        c_width_2d[depth] = new TCanvas(Form("c_width_2d_depth%i", depth+1), Form("c_width_2d_depth%i", depth+1), 800,600);
        c_width_2d[depth]->Divide(2,2);
        c_width_2d[depth]->cd(1); h2_width[depth][0]->Draw("colz");
        c_width_2d[depth]->cd(2); h2_width[depth][1]->Draw("colz");
        c_width_2d[depth]->cd(3); h2_width[depth][2]->Draw("colz");
        c_width_2d[depth]->cd(4); h2_width[depth][3]->Draw("colz");
        c_width_2d[depth]->Print(Form("Fig/c_width_2d_depth%i.pdf", depth+1));
    } 
    // 1d 
    for(int isub=0; isub<5; isub++)
    {
        TString subdet="HB"; if(isub==1) subdet="HE"; if(isub==2) subdet="HO"; if(isub==3) subdet="HF";if(isub==4) subdet="QIE11"; 
        c_mean_1d[isub] = new TCanvas(Form("c_mean_1d_%s", subdet.Data()), Form("c_mean_1d_%s", subdet.Data()), 800,600);
        c_mean_1d[isub]->Divide(2,2);
        c_mean_1d[isub]->cd(1); c_mean_1d[isub]->cd(1)->SetLogy(1); h1_mean[isub][0]->Draw("hist");
        c_mean_1d[isub]->cd(2); c_mean_1d[isub]->cd(2)->SetLogy(1); h1_mean[isub][1]->Draw("hist");
        c_mean_1d[isub]->cd(3); c_mean_1d[isub]->cd(3)->SetLogy(1); h1_mean[isub][2]->Draw("hist");
        c_mean_1d[isub]->cd(4); c_mean_1d[isub]->cd(4)->SetLogy(1); h1_mean[isub][3]->Draw("hist");
        c_mean_1d[isub]->Print(Form("Fig/c_mean_1d_%s.pdf", subdet.Data()));
        
        c_width_1d[isub] = new TCanvas(Form("c_width_1d_%s", subdet.Data()), Form("c_width_1d_%s", subdet.Data()), 800,600);
        c_width_1d[isub]->Divide(2,2);
        c_width_1d[isub]->cd(1); c_width_1d[isub]->cd(1)->SetLogy(1); h1_width[isub][0]->Draw("hist");
        c_width_1d[isub]->cd(2); c_width_1d[isub]->cd(2)->SetLogy(1); h1_width[isub][1]->Draw("hist");
        c_width_1d[isub]->cd(3); c_width_1d[isub]->cd(3)->SetLogy(1); h1_width[isub][2]->Draw("hist");
        c_width_1d[isub]->cd(4); c_width_1d[isub]->cd(4)->SetLogy(1); h1_width[isub][3]->Draw("hist");
        c_width_1d[isub]->Print(Form("Fig/c_width_1d_%s.pdf", subdet.Data()));
    }

    fin.close();
}

void HCALPedestalCompareAnalysis(TString CompareFile="compare.root") 
{ 
    //
    TFile* comparefile = TFile::Open(CompareFile);

    //
    // 1D
    //

    // 1D plots 
    TH1D *h1[2][5][4]; // [2]: mean, width / [4]: HB,HE,HO,HF / [4]: capid
    TH1D *h1_capsum[2][5]; // [2]: mean, width / [4]: HB,HE,HO,HF

    for(int var=0; var<2; var++)
    { 
        for(int subdet=0; subdet<5; subdet++) 
        {   
            // 
            TString Var="hm"; if(var==1) Var="hr"; 
            TString SubDet="HB"; if(subdet==1) SubDet="HE"; 
                                 if(subdet==2) SubDet="HO"; 
                                 if(subdet==3) SubDet="HF"; 
                                 if(subdet==4) SubDet="QIE11"; 
            
            TCanvas *c = new TCanvas("c","c",800,600);
            c->Divide(2,2); 

            for(int capid=0; capid<4; capid++) 
            { 
                h1[var][subdet][capid] = (TH1D*)comparefile->Get(Form("%s/%s%i_diff", SubDet.Data(), Var.Data(), capid));
                // 
                c->cd(capid+1);
                c->cd(capid+1)->SetLogy(1);
                h1[var][subdet][capid]->SetStats(0); 
                h1[var][subdet][capid]->SetTitle(Form("Difference of %s: %s capid=%i", SubDet.Data(), Var.Data(), capid));
                h1[var][subdet][capid]->Draw("hist"); 

                TLatex *tex_mean = new TLatex(0.7,0.8,Form("Mean = %.3f", h1[var][subdet][capid]->GetMean()));
                tex_mean->SetNDC();
                tex_mean->SetTextSize(0.06);
                tex_mean->SetLineWidth(2);
                tex_mean->SetTextAlign(22);
                
                TLatex *tex_rms = new TLatex(0.7,0.7,Form("RMS = %.3f", h1[var][subdet][capid]->GetRMS()));
                tex_rms->SetNDC();
                tex_rms->SetTextSize(0.06);
                tex_rms->SetLineWidth(2);
                tex_rms->SetTextAlign(22); 

                tex_mean->Draw();
                tex_rms->Draw();
            } 
            TString plotname  = Form("%s_%s_diff.pdf", SubDet.Data(), Var.Data());
            cout <<  "Wrigint this plot: " << plotname << endl;
            c->Print("Fig/Compare/1D/"+plotname);
            delete c; 

            if(var==0) 
            { 
                TCanvas *c_capsum = new TCanvas("c_capsum","c_capsum",400,300);
                h1_capsum[var][subdet] = (TH1D*)comparefile->Get(Form("%s/%s_diff", SubDet.Data(), Var.Data()));
                
                c_capsum->cd(1);
                c_capsum->cd(1)->SetLogy(1);
                h1_capsum[var][subdet]->SetStats(0); 
                h1_capsum[var][subdet]->SetTitle(Form("Difference of %s: %s", SubDet.Data(), Var.Data()));
                h1_capsum[var][subdet]->Draw("hist"); 

                TLatex *tex_mean = new TLatex(0.7,0.8,Form("Mean = %.3f", h1_capsum[var][subdet]->GetMean()));
                tex_mean->SetNDC();
                tex_mean->SetTextSize(0.06);
                tex_mean->SetLineWidth(2);
                tex_mean->SetTextAlign(22);
                
                TLatex *tex_rms = new TLatex(0.7,0.7,Form("RMS = %.3f", h1_capsum[var][subdet]->GetRMS()));
                tex_rms->SetNDC();
                tex_rms->SetTextSize(0.06);
                tex_rms->SetLineWidth(2);
                tex_rms->SetTextAlign(22); 

                tex_mean->Draw();
                tex_rms->Draw();
            
                TString plotname_capsum  = Form("%s_%s_capsum_diff.pdf", SubDet.Data(), Var.Data());
                cout <<  "Wrigint this plot: " << plotname_capsum << endl;
                c_capsum->Print("Fig/Compare/1D/"+plotname_capsum);
                delete c_capsum; 
            } 
        }
    }
    
    // 
    // 2D
    // 
    
    // 2D plots
    TH2D *h2[2][5][4][4]; // [2]: mean, width / [4]: HB,HE,HO,HF / [4]: depth / [4]: capid
  
    // from: https://root.cern.ch/phpBB3/viewtopic.php?t=18987
    const UInt_t Number = 5;
    Double_t Red[Number]    = { 1.00, 0.00, 1.00, 1.00, 0.00};
    Double_t Green[Number]  = { 0.00, 1.00, 1.00, 0.00, 0.00};
    Double_t Blue[Number]   = { 0.00, 0.00, 1.00, 1.00, 1.00};
    Double_t Length[Number] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
    Int_t nb=50;
    TColor::CreateGradientColorTable(Number,Length,Blue,Green,Red,nb);

    for(int var=0; var<2; var++)
    { 
        for(int subdet=0; subdet<5; subdet++) 
        {   
            for(int depth=0; depth<4; depth++) 
            { 
                // 
                if(subdet==0 && depth>1) continue;
                if(subdet==1 && depth>2) continue;
                if(subdet==2 && depth<3) continue;
                if(subdet==3 && depth>1) continue;
                   
                // 
                TString Var="rm"; if(var==1) Var="rr"; 
                TString SubDet="HB"; if(subdet==1) SubDet="HE"; 
                                     if(subdet==2) SubDet="HO"; 
                                     if(subdet==3) SubDet="HF"; 
                                     if(subdet==4) SubDet="QIE11"; 
                
                TCanvas *c = new TCanvas("c","c",800,600);
                c->Divide(2,2);
                
                for(int capid=0; capid<4; capid++) 
                {

                    // Check first if the histogram exists
                    cout << Form("Analyzing %s %s capid=%i depth=%i", SubDet.Data(), Var.Data(), capid, depth+1) << endl; 

                    // 
                    h2[var][subdet][depth][capid] = (TH2D*)comparefile->Get(Form("%s/2D/h2%s%i_diff_%i", SubDet.Data(), Var.Data(), capid, depth+1));
                    
                    // 
                    c->cd(capid+1);
                    h2[var][subdet][depth][capid]->SetContour(nb);
                    h2[var][subdet][depth][capid]->SetMaximum(1);
                    h2[var][subdet][depth][capid]->SetMinimum(-1);
                    h2[var][subdet][depth][capid]->SetStats(0);
                    h2[var][subdet][depth][capid]->SetTitle(Form("Difference of %s: %s depth=%i capid=%i", SubDet.Data(), Var.Data(), depth+1, capid));
                    h2[var][subdet][depth][capid]->Draw("colz");

                } 
                
                TString plotname  = Form("%s_2D_%s_depth%i_diff.pdf", SubDet.Data(), Var.Data(), depth+1);
                cout <<  "Wrigint this plot: " << plotname << endl;
                c->Print("Fig/Compare/2D/"+plotname);
                delete c;
            }
        }
    }
} 

void HCALPedestalAnalysis() 
{ 
    //
    // Analysis: visualization of a ped table
    //
    HCALPedestalTableAnalysis("PedestalTable_outputFile_run293777_run293775HEP17.txt"); 
    
    //
    // Analysis: compare two tables 
    //
    // Need to run ped_compare.py first to get the input file that contains histroms
    // This function visualizes them  
    //HCALPedestalCompareAnalysis("output_NewvsDB.root"); 

}
