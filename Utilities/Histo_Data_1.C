#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));


double Start=50;
double End=150;
int bin=50;
int bin1 = 16;
double l=74;
double u=106;
double s = 250;

static const std::vector<float> PN_scores {0.768733, 0.936303, 0.977057};

int PN_DIM =  PN_scores.size(); //number of PN-MD score bins

static const std::vector<float> DDX_scores {0.022193, 0.0939224, 0.238822};

int DDX_DIM =  DDX_scores.size(); //number of DDX score bins

void Histo_Data_1(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //gStyle->SetOptStat(1000000001);
  //gStyle->SetOptFit(1111);

  TH1::AddDirectory(kFALSE); 
  
  std::string Y_string{"Events / "+ std::to_string(((End/s-Start/s)/bin))};
  const char * Y_title{Y_string.c_str()};

  //histogams sidebands

  //PN-MD softdropmass histograms
  TH1F *PN_Data[PN_DIM];
  
  for(int iPN=0; iPN<PN_DIM;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)};
    const char * nome{nomestringa.c_str()};    
    PN_Data[iPN]= new TH1F(nome, "x", bin, Start/s, End/s);
    PN_Data[iPN]->GetXaxis()->SetTitle("x");    
    PN_Data[iPN]->GetYaxis()->SetTitle(Y_title);
    if(iPN==0){PN_Data[iPN]->SetTitle("PN-MD_BBvsQCD loose WP");}
    if(iPN==1){PN_Data[iPN]->SetTitle("PN-MD_BBvsQCD medium WP");}
    if(iPN==2){PN_Data[iPN]->SetTitle("PN-MD_BBvsQCD tight WP");}
    PN_Data[iPN]->Sumw2();  
  }

  //DDX softdropmass histograms
  TH1F *DDX_Data[PN_DIM];
  
  for(int iPN=0; iPN<DDX_DIM;iPN++){
    std::string nomestringa{"DDXhist"+std::to_string(iPN)};
    const char * nome{nomestringa.c_str()};    
    DDX_Data[iPN]= new TH1F(nome, "x", bin, Start/s, End/s);
    DDX_Data[iPN]->GetXaxis()->SetTitle("x");    
    DDX_Data[iPN]->GetYaxis()->SetTitle(Y_title);
    if(iPN==0){DDX_Data[iPN]->SetTitle("DDX_BBvsQCD loose WP");}
    if(iPN==1){DDX_Data[iPN]->SetTitle("DDX_BBvsQCD medium WP");}
    if(iPN==2){DDX_Data[iPN]->SetTitle("DDX_BBvsQCD tight WP");}
    DDX_Data[iPN]->Sumw2();  
  }

  //histogams signal region

  //PN-MD softdropmass histograms
  TH1F *PN_Data_s[PN_DIM];
  
  for(int iPN=0; iPN<PN_DIM;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};    
    PN_Data_s[iPN]= new TH1F(nome, "x", bin1, l/s, u/s);
    PN_Data_s[iPN]->GetXaxis()->SetTitle("x");    
    PN_Data_s[iPN]->GetYaxis()->SetTitle(Y_title);
    if(iPN==0){PN_Data_s[iPN]->SetTitle("PN-MD_BBvsQCD loose WP");}
    if(iPN==1){PN_Data_s[iPN]->SetTitle("PN-MD_BBvsQCD medium WP");}
    if(iPN==2){PN_Data_s[iPN]->SetTitle("PN-MD_BBvsQCD tight WP");}
    PN_Data_s[iPN]->Sumw2();  
  }

  //DDX softdropmass histograms
  TH1F *DDX_Data_s[PN_DIM];
  
  for(int iPN=0; iPN<DDX_DIM;iPN++){
    std::string nomestringa{"DDXhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};    
    DDX_Data_s[iPN]= new TH1F(nome, "x", bin1, l/s, u/s);
    DDX_Data_s[iPN]->GetXaxis()->SetTitle("x");    
    DDX_Data_s[iPN]->GetYaxis()->SetTitle(Y_title);
    if(iPN==0){DDX_Data_s[iPN]->SetTitle("DDX_BBvsQCD loose WP");}
    if(iPN==1){DDX_Data_s[iPN]->SetTitle("DDX_BBvsQCD medium WP");}
    if(iPN==2){DDX_Data_s[iPN]->SetTitle("DDX_BBvsQCD tight WP");}
    DDX_Data_s[iPN]->Sumw2();  
  }
  
  TChain *T_Data = new TChain("passedEvents");  // name of the tree is the argument
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/Data/RunC.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/Data/RunD.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/Data/RunE.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/Data/RunF.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/Data/RunG.root");
  

  //definition of the variables and branches
  Double_t AK8_sdmJet0, sdm_Higher_DDX_XbbVsQCD, sdm_Higher_PN_MD_XbbVsQCD;
  Float_t AK8_ptJet0, PNMD_XbbVsQCD, DDX_XbbVsQCD;//pt, softdropmass and scores of the leading AK8 jet
  Float_t Higher_PN_MD_XbbVsQCD;//BBvsQCD and softdropmass of the leading_PNMD AK8 jet
  Float_t Higher_DDX_XbbVsQCD;//BBvsQCD and softdropmass of the leading_DDX AK8 jet
  
  T_Data->SetBranchAddress("AK8_ptJet0", &AK8_ptJet0);
  T_Data->SetBranchAddress("AK8_sdmJet0", &AK8_sdmJet0);
  T_Data->SetBranchAddress("PNMD_XbbVsQCD", &PNMD_XbbVsQCD);
  T_Data->SetBranchAddress("DDX_XbbVsQCD", &DDX_XbbVsQCD);
  T_Data->SetBranchAddress("Higher_PN_MD_XbbVsQCD", &Higher_PN_MD_XbbVsQCD);
  T_Data->SetBranchAddress("sdm_Higher_PN_MD_XbbVsQCD", &sdm_Higher_PN_MD_XbbVsQCD);
  T_Data->SetBranchAddress("Higher_DDX_XbbVsQCD", &Higher_DDX_XbbVsQCD);
  T_Data->SetBranchAddress("sdm_Higher_DDX_XbbVsQCD", &sdm_Higher_DDX_XbbVsQCD);  
  
  
    
  Int_t nentries = (Int_t)T_Data->GetEntries();
  
  cout<<nentries<<" entries"<<endl;

  for (int i = 0; i < nentries; i++){    
    T_Data->GetEntry(i);     
    if(i%10000 == 0){cout<<i<<"/"<<nentries<<" entries done"<<endl;}
    
    for(int iPN=0; iPN<PN_DIM;iPN++){
      if(sdm_Higher_PN_MD_XbbVsQCD>l && sdm_Higher_PN_MD_XbbVsQCD<u && Higher_PN_MD_XbbVsQCD>PN_scores.at(iPN)){
	PN_Data_s[iPN]->Fill(sdm_Higher_PN_MD_XbbVsQCD/s);
      }
      if(((sdm_Higher_PN_MD_XbbVsQCD>Start && sdm_Higher_PN_MD_XbbVsQCD<l) || (sdm_Higher_PN_MD_XbbVsQCD >u && sdm_Higher_PN_MD_XbbVsQCD<End)) && Higher_PN_MD_XbbVsQCD>PN_scores.at(iPN)){
	PN_Data[iPN]->Fill(sdm_Higher_PN_MD_XbbVsQCD/s);
      }
    } 

    for(int iPN=0; iPN<DDX_DIM;iPN++){
      if(sdm_Higher_DDX_XbbVsQCD>l && sdm_Higher_DDX_XbbVsQCD<u && Higher_DDX_XbbVsQCD>DDX_scores.at(iPN)){
	DDX_Data_s[iPN]->Fill(sdm_Higher_DDX_XbbVsQCD/s);
      }
      if(((sdm_Higher_DDX_XbbVsQCD>Start && sdm_Higher_DDX_XbbVsQCD<l) || (sdm_Higher_DDX_XbbVsQCD >u && sdm_Higher_DDX_XbbVsQCD<End)) && Higher_DDX_XbbVsQCD>DDX_scores.at(iPN)){
	DDX_Data[iPN]->Fill(sdm_Higher_DDX_XbbVsQCD/s);
      }
    } 
    
  }

  
  //saving histograms in a Tfile

  TFile *myfile;
  myfile = TFile::Open("Data_histos_1.root","RECREATE");
  myfile->cd();

  for(int i=0; i<PN_DIM;i++){
    PN_Data[i]->Write();
    PN_Data_s[i]->Write();
  }

  for(int i=0; i<DDX_DIM;i++){
    DDX_Data[i]->Write();
    DDX_Data_s[i]->Write();
  }

  myfile->Close();
  delete myfile;
    
}
