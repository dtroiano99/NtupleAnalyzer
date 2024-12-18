#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>

//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));

double Start=38;
double End=198;
int bin=160;
int bin1 = 56;
double l=70;
double u=126;
double s = 200;


static const std::vector<float> PN_scores {0.641, 0.875,0.957, 0.988,  1};

int PN_DIM =  PN_scores.size()-1; //number of PN-MD score bins

void Histo_MC(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //gStyle->SetOptStat(1000000001);
  //gStyle->SetOptFit(1111);

  TH1::AddDirectory(kFALSE); 
  
  std::string Y_string{"Events / "+ std::to_string(((End/s-Start/s)/bin))};
  const char * Y_title{Y_string.c_str()};

  //histogams signal region

  //PN-MD softdropmass histograms
  TH1F *PN_Zbb_s[PN_DIM];
  
  for(int iPN=0; iPN<PN_DIM;iPN++){
    std::string nomestringa{"Zbbhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};    
    PN_Zbb_s[iPN]= new TH1F(nome, "x", bin1, l/s, u/s);
    PN_Zbb_s[iPN]->GetXaxis()->SetTitle("m_{SD} / 200 GeV");    
    PN_Zbb_s[iPN]->GetYaxis()->SetTitle(Y_title);
    float start = floor(PN_scores.at(iPN) * 10000) / 10000;
    std::stringstream ss1;
    ss1 <<start;
    std::string str1 = ss1.str();
    float end = floor(PN_scores.at(iPN+1) * 10000) / 10000;
    std::stringstream ss2;
    ss2 <<end;
    std::string str2 = ss2.str();
    std::string nomestringa_title{str1+" < PN-MD_{BBvsQCD} #leq "+ str2};
    const char * title{nomestringa_title.c_str()}; 
    PN_Zbb_s[iPN]->SetTitle(title);
    PN_Zbb_s[iPN]->Sumw2();  
  }

  TH1F *PN_Zcc_s[PN_DIM];
  
  for(int iPN=0; iPN<PN_DIM;iPN++){
    std::string nomestringa{"Zcchist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};    
    PN_Zcc_s[iPN]= new TH1F(nome, "x", bin1, l/s, u/s);
    PN_Zcc_s[iPN]->GetXaxis()->SetTitle("m_{SD} / 200 GeV");    
    PN_Zcc_s[iPN]->GetYaxis()->SetTitle(Y_title);
    float start = floor(PN_scores.at(iPN) * 10000) / 10000;
    std::stringstream ss1;
    ss1 <<start;
    std::string str1 = ss1.str();
    float end = floor(PN_scores.at(iPN+1) * 10000) / 10000;
    std::stringstream ss2;
    ss2 <<end;
    std::string str2 = ss2.str();
    std::string nomestringa_title{str1+" < PN-MD_{BBvsQCD} #leq "+ str2};
    const char * title{nomestringa_title.c_str()}; 
    PN_Zcc_s[iPN]->SetTitle(title);
    PN_Zcc_s[iPN]->Sumw2();  
  }

  TH1F *PN_Zll_s[PN_DIM];
  
  for(int iPN=0; iPN<PN_DIM;iPN++){
    std::string nomestringa{"Zllhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};    
    PN_Zll_s[iPN]= new TH1F(nome, "x", bin1, l/s, u/s);
    PN_Zll_s[iPN]->GetXaxis()->SetTitle("m_{SD} / 200 GeV");    
    PN_Zll_s[iPN]->GetYaxis()->SetTitle(Y_title);
    float start = floor(PN_scores.at(iPN) * 10000) / 10000;
    std::stringstream ss1;
    ss1 <<start;
    std::string str1 = ss1.str();
    float end = floor(PN_scores.at(iPN+1) * 10000) / 10000;
    std::stringstream ss2;
    ss2 <<end;
    std::string str2 = ss2.str();
    std::string nomestringa_title{str1+" < PN-MD_{BBvsQCD} #leq "+ str2};
    const char * title{nomestringa_title.c_str()}; 
    PN_Zll_s[iPN]->SetTitle(title);
    PN_Zll_s[iPN]->Sumw2();  
  }

  float BB[PN_DIM];
  float CC[PN_DIM];
  float QQ[PN_DIM];

for(int iPN=0; iPN<PN_DIM;iPN++){
    BB[iPN]=0;
    CC[iPN]=0;
    QQ[iPN]=0;
  }
 
  TChain *T_Data = new TChain("passedEvents");  // name of the tree is the argument
  T_Data->Add("../../new_wind/pt_validation/Full_ZtoQQ.root");
  
  //definition of the variables and branches
  Double_t AK8_sdmJet0;
  Float_t PNMD_XbbVsQCD;//softdropmass and score of the leading AK8 jet
  Double_t lumiWeight;

  std::vector<int> *quark_flavour = 0;
  
  T_Data->SetBranchAddress("AK8_sdmJet0", &AK8_sdmJet0);
  T_Data->SetBranchAddress("PNMD_XbbVsQCD", &PNMD_XbbVsQCD);
  T_Data->SetBranchAddress("lumiWeight", &lumiWeight);

  T_Data->SetBranchAddress("quark_flavour",&quark_flavour);
    
  Int_t nentries = (Int_t)T_Data->GetEntries();
  //Int_t nentries = 10000;
  
  cout<<nentries<<" entries"<<endl;
  for (int i = 0; i < nentries; i++){
    T_Data->GetEntry(i);
    if(i%10000 == 0){cout<<i<<"/"<<nentries<<" entries done"<<endl;}
    /*
      cout<<"Event number "<<i<<endl;
      for(int k=0;k<quark_flavour->size();k++){
      cout<<quark_flavour->at(k)<<endl;
    }
    */
    if(AK8_sdmJet0>l && AK8_sdmJet0<u){
      
      for(int iPN=0; iPN<PN_DIM;iPN++){
	
	if(PNMD_XbbVsQCD>PN_scores.at(iPN) && PNMD_XbbVsQCD<=PN_scores.at(iPN+1)){
	  if(abs(quark_flavour->at(0))==5){
	    BB[iPN]=BB[iPN]+lumiWeight;
	    PN_Zbb_s[iPN]->Fill(AK8_sdmJet0/s, lumiWeight);
	  }
	  if(abs(quark_flavour->at(0))==4){
	    CC[iPN]=CC[iPN]+lumiWeight;
	    PN_Zcc_s[iPN]->Fill(AK8_sdmJet0/s, lumiWeight);
	  }
	  if(abs(quark_flavour->at(0))<4){
	    QQ[iPN]=QQ[iPN]+lumiWeight;
	    PN_Zll_s[iPN]->Fill(AK8_sdmJet0/s, lumiWeight);
	  }
	}
      }
    } 

  }
    
    for(int iPN=0; iPN<PN_DIM;iPN++){
      cout<<"score region # "<<iPN<<endl;
      cout<<"BB events: "<<BB[iPN]<<"; total fraction: "<<100*BB[iPN]/(BB[iPN]+CC[iPN]+QQ[iPN])<<"%"<<endl;
      cout<<"CC events: "<<CC[iPN]<<"; total fraction: "<<100*CC[iPN]/(BB[iPN]+CC[iPN]+QQ[iPN])<<"%"<<endl;
      cout<<"QQ events: "<<QQ[iPN]<<"; total fraction: "<<100*QQ[iPN]/(BB[iPN]+CC[iPN]+QQ[iPN])<<"%"<<endl;
    }


  
  
  //saving histograms in a Tfile
  cout<<"writing root files"<<endl;
  TFile *myfile;
  myfile = TFile::Open("MC_Z2_histos.root","RECREATE");
  myfile->cd();

  for(int i=0; i<PN_DIM;i++){
    PN_Zbb_s[i]->Write();
    PN_Zcc_s[i]->Write();
    PN_Zll_s[i]->Write();
  }

  cout<<"root files completed"<<endl;
  myfile->Close();
  delete myfile;
  
}
