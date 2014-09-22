void plot() {
  
  TFile *file = TFile::Open("alberto_histo_220426.root");
  //  TFile *file = TFile::Open("alberto_histo_217874.root");
  //  TFile* file = TFile::Open("217384_histo.root");
  TDirectoryFile* dir = (TDirectoryFile*)file->Get("pedTuner");
  dir->cd();

  TTree* tree = (TTree*)dir->Get("tree");

  HO_1D->Draw();
  //  c1->Print("HO_Ped.gif");

  /*
  HO_1D0->Draw();
  c1->Print("HO_Ped0.gif");
  HO_1D1->Draw();
  c1->Print("HO_Ped1.gif");
  HO_1D2->Draw();
  c1->Print("HO_Ped2.gif");
  HO_1D3->Draw();
  c1->Print("HO_Ped3.gif");
  */
  
  //  TString cut = "subdet==3&&ped<8.8";
  //  TString var = "ped:dac";
  //  mytree->Draw(var,cut,"BOX");
  //  c1->Print("ped_vs_dac.gif");


}

void plot_HBHE_before() {
  
  TFile *file = TFile::Open("alberto_histo_218117.root");
  TDirectoryFile* dir = (TDirectoryFile*)file->Get("pedTuner");
  dir->cd();

  TTree* tree = (TTree*)dir->Get("tree");

  int rebin= 4;

  HB_1D->Rebin(rebin);
  HB_1D->Draw();
  c1->Print("HB_Ped_before.gif");

  HEP_1D->Rebin(rebin);
  HEP_1D->Draw();
  c1->Print("HEP_Ped_before.gif");

  HEM_1D->Rebin(rebin);
  HEM_1D->Draw();
  c1->Print("HEM_Ped_before.gif");

  HFP_1D->Rebin(rebin);
  HFP_1D->Draw();
  c1->Print("HFP_Ped_before.gif");

  //  HFM_1D->Rebin(rebin);
  //  HFM_1D->Draw();
  //  c1->Print("HFM_Ped_before.gif");


  
  //  TString cut = "subdet==3&&ped<8.8";
  //  TString var = "ped:dac";
  //  mytree->Draw(var,cut,"BOX");
  //  c1->Print("ped_vs_dac.gif");
}


void plot_HBHE_after() {
  
  TFile *file = TFile::Open("alberto_histo_218513.root");
  //  TFile *file = TFile::Open("alberto_histo_218469.root");
  TDirectoryFile* dir = (TDirectoryFile*)file->Get("pedTuner");
  dir->cd();

  TTree* tree = (TTree*)dir->Get("tree");

  int rebin= 4;

  HB_1D->Rebin(rebin);
  HB_1D->Draw();
  c1->Print("HB_Ped_after.gif");

  /*
  HBP_1D->Rebin(rebin);
  HBP_1D->Draw();
  c1->Print("HBP_Ped_after.gif");

  HBM_1D->Rebin(rebin);
  HBM_1D->Draw();
  c1->Print("HBM_Ped_after.gif");

  HEP_1D->Rebin(rebin);
  HEP_1D->Draw();
  c1->Print("HEP_Ped_after.gif");

  HEM_1D->Rebin(rebin);
  HEM_1D->Draw();
  c1->Print("HEM_Ped_after.gif");

  HFP_1D->Rebin(rebin);
  HFP_1D->Draw();
  c1->Print("HFP_Ped_after.gif");
  */
  //  HFM_1D->Rebin(rebin);
  //  HFM_1D->Draw();
  //  c1->Print("HFM_Ped_after.gif");

}
