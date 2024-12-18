#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));


double Start=38;
double End=198;
int bin=80;

void Large_ZtoQQ_msd(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //gStyle->SetOptStat(1000000001);
  //gStyle->SetOptFit(1111);

  TH1::AddDirectory(kFALSE); 
  
  //std::string Y_string{"Events / "+ std::to_string(((End/s-Start/s)/bin))};
  std::string Y_string{"Events / 2 GeV"};
  const char * Y_title{Y_string.c_str()};

  //histogams signal region

  //PN-MD softdropmass histograms
  TH1F *PN_Data_s; 
  PN_Data_s= new TH1F("PNhist_s", "x", bin, Start, End);
  PN_Data_s->GetXaxis()->SetTitle("m_{SD} [GeV]");    
  PN_Data_s->GetYaxis()->SetTitle(Y_title);
  PN_Data_s->SetTitle("");
  PN_Data_s->Sumw2();  
  

 
  TChain *T_Data = new TChain("passedEvents");  // name of the tree is the argument
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT800toInf_2023C.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT800toInf_2023D.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT600to800_2023C.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT600to800_2023D.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT400to600_2023C.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT400to600_2023D.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT200to400_2023C.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT200to400_2023D.root");
 

  
  //definition of the variables and branches 
  Float_t PNMD_XbbVsQCD;//score of the leading AK8 jet
  Double_t lumiWeight;
  Double_t AK8_sdmJet0;  
    
  T_Data->SetBranchAddress("AK8_sdmJet0", &AK8_sdmJet0);
  T_Data->SetBranchAddress("PNMD_XbbVsQCD", &PNMD_XbbVsQCD);
  T_Data->SetBranchAddress("lumiWeight", &lumiWeight);
    
  Int_t nentries = (Int_t)T_Data->GetEntries();

  cout<<nentries<<" entries"<<endl;
  for (int i = 0; i < nentries; i++){
    T_Data->GetEntry(i);
    if(i%10000 == 0){cout<<i<<"/"<<nentries<<" entries done"<<endl;}
    //Double_t AK8_sdmJet0 = AK8PuppiJets_rawsoftdropmass->at(ptIdx->at(0));
    //for(int iPN=0; iPN<PN_DIM;iPN++){
    if(AK8_sdmJet0>Start && AK8_sdmJet0<End ){PN_Data_s->Fill(AK8_sdmJet0, lumiWeight);}
  } 



  
  
  //saving histograms in a Tfile
  cout<<"writing root files"<<endl;
  TFile *myfile;
  myfile = TFile::Open("Large_ZtoQQ_histos.root","RECREATE");
  myfile->cd();


    PN_Data_s->Write();
  

  cout<<"root files completed"<<endl;
  myfile->Close();
  delete myfile;
  

  //PN_Data_s[0]->Draw("hist");
  
}
