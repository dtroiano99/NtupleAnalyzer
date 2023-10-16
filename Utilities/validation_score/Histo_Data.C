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

static const std::vector<float> PN_scores {0.83226, 0.913256, 0.949388, 0.968723, 0.981101, 0.989198, 0.993106, 0.996553, 1};

int PN_DIM =  PN_scores.size()-1; //number of PN-MD score bins

static const std::vector<float> DDX_scores {0.0365136, 0.0753757,0.127416,0.195188, 0.286469, 0.414133, 0.581553, 0.770378,1};

int DDX_DIM =  DDX_scores.size()-1; //number of DDX score bins

void Histo_Data(){
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
    float start = floor(PN_scores.at(iPN) * 10000) / 10000;
    std::stringstream ss1;
    ss1 <<start;
    std::string str1 = ss1.str();
    float end = floor(PN_scores.at(iPN+1) * 10000) / 10000;
    std::stringstream ss2;
    ss2 <<end;
    std::string str2 = ss2.str();
    std::string nomestringa_title{str1+" < PN-MD_BBvsQCD #leq "+ str2};
    const char * title{nomestringa_title.c_str()}; 
    PN_Data[iPN]->SetTitle(title);
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
    float start = floor(DDX_scores.at(iPN) * 10000) / 10000;
    std::stringstream ss1;
    ss1 <<start;
    std::string str1 = ss1.str();
    float end = floor(DDX_scores.at(iPN+1) * 10000) / 10000;
    std::stringstream ss2;
    ss2 <<end;
    std::string str2 = ss2.str();
    std::string nomestringa_title{str1+" < DDX_BBvsQCD #leq "+ str2};
    const char * title{nomestringa_title.c_str()}; 
    DDX_Data[iPN]->SetTitle(title);    
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
    float start = floor(PN_scores.at(iPN) * 10000) / 10000;
    std::stringstream ss1;
    ss1 <<start;
    std::string str1 = ss1.str();
    float end = floor(PN_scores.at(iPN+1) * 10000) / 10000;
    std::stringstream ss2;
    ss2 <<end;
    std::string str2 = ss2.str();
    std::string nomestringa_title{str1+" < PN-MD_BBvsQCD #leq "+ str2};
    const char * title{nomestringa_title.c_str()}; 
    PN_Data_s[iPN]->SetTitle(title);    
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
    float start = floor(DDX_scores.at(iPN) * 10000) / 10000;
    std::stringstream ss1;
    ss1 <<start;
    std::string str1 = ss1.str();
    float end = floor(DDX_scores.at(iPN+1) * 10000) / 10000;
    std::stringstream ss2;
    ss2 <<end;
    std::string str2 = ss2.str();
    std::string nomestringa_title{str1+" < DDX_BBvsQCD #leq "+ str2};
    const char * title{nomestringa_title.c_str()}; 
    DDX_Data_s[iPN]->SetTitle(title);   
    DDX_Data_s[iPN]->Sumw2();  
  }

  TChain *T_Data = new TChain("passedEvents");  // name of the tree is the argument
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/Data/RunC.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/Data/RunD.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/Data/RunE.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/Data/RunF.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/Data/RunG.root");
  

  //definition of the variables and branches
  Double_t AK8_sdmJet0;
  Float_t AK8_ptJet0, PNMD_XbbVsQCD, DDX_XbbVsQCD;//pt, softdropmass and scores of the leading AK8 jet
  
  T_Data->SetBranchAddress("AK8_ptJet0", &AK8_ptJet0);
  T_Data->SetBranchAddress("AK8_sdmJet0", &AK8_sdmJet0);
  T_Data->SetBranchAddress("PNMD_XbbVsQCD", &PNMD_XbbVsQCD);
  T_Data->SetBranchAddress("DDX_XbbVsQCD", &DDX_XbbVsQCD);
    
  Int_t nentries = (Int_t)T_Data->GetEntries();

  cout<<nentries<<" entries"<<endl;
  for (int i = 0; i < nentries; i++){      
    T_Data->GetEntry(i); 
    if(i%1000 == 0){cout<<i<<"/"<<nentries<<" entries done"<<endl;}
    for(int iPN=0; iPN<PN_DIM;iPN++){
      if(AK8_sdmJet0>l && AK8_sdmJet0<u && PNMD_XbbVsQCD>PN_scores.at(iPN) && PNMD_XbbVsQCD<=PN_scores.at(iPN+1)){PN_Data_s[iPN]->Fill(AK8_sdmJet0/s);}
       if(((AK8_sdmJet0>Start && AK8_sdmJet0<l) || (AK8_sdmJet0 >u && AK8_sdmJet0<End)) && PNMD_XbbVsQCD>PN_scores.at(iPN) && PNMD_XbbVsQCD<=PN_scores.at(iPN+1)){PN_Data[iPN]->Fill(AK8_sdmJet0/s);}
    } 

    for(int iPN=0; iPN<DDX_DIM;iPN++){
      if(AK8_sdmJet0>l && AK8_sdmJet0<u && DDX_XbbVsQCD>DDX_scores.at(iPN) && DDX_XbbVsQCD<=DDX_scores.at(iPN+1)){DDX_Data_s[iPN]->Fill(AK8_sdmJet0/s);}
      if(((AK8_sdmJet0>Start && AK8_sdmJet0<l) || (AK8_sdmJet0 >u && AK8_sdmJet0<End)) && DDX_XbbVsQCD>DDX_scores.at(iPN) && DDX_XbbVsQCD<=DDX_scores.at(iPN+1)){DDX_Data[iPN]->Fill(AK8_sdmJet0/s);}
    } 
    
  }

  
  //saving histograms in a Tfile
  cout<<"writing root file"<<endl;

  TFile *myfile;
  myfile = TFile::Open("Data_histos.root","RECREATE");
  myfile->cd();

  for(int i=0; i<PN_DIM;i++){
    PN_Data[i]->Write();
    PN_Data_s[i]->Write();
  }

  for(int i=0; i<DDX_DIM;i++){
    DDX_Data[i]->Write();
    DDX_Data_s[i]->Write();
  }
  cout<<"root file completed"<<endl;

  myfile->Close();
  delete myfile;
  
}
