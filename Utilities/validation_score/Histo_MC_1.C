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

static const std::vector<float> PN_scores {0.768733,0.885883,0.936303, 0.961613,0.977057,0.986697,0.992354,0.996177, 1};

int PN_DIM =  PN_scores.size()-1; //number of PN-MD score bins

static const std::vector<float> DDX_scores {0.022193,0.0502867,0.0939224,0.153532, 0.238822,0.356848,0.527736,0.738001,1};

int DDX_DIM =  DDX_scores.size()-1; //number of DDX score bins

void Histo_MC_1(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //gStyle->SetOptStat(1000000001);
  //gStyle->SetOptFit(1111);

  TH1::AddDirectory(kFALSE); 
  
  std::string Y_string{"Events / "+ std::to_string(((End/s-Start/s)/bin))};
  const char * Y_title{Y_string.c_str()};

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
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/MC/ZJetsToQQ_HT200to400.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/MC/ZJetsToQQ_HT400to600.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/MC/ZJetsToQQ_HT600to800.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/MC/ZJetsToQQ_HT800toInf.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/MC/ZJetsToQQ_HT200to400_postEE.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/MC/ZJetsToQQ_HT400to600_postEE.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/MC/ZJetsToQQ_HT600to800_postEE.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/new_variable/MC/ZJetsToQQ_HT800toInf_postEE.root");
  
  //definition of the variables and branches
  Double_t AK8_sdmJet0, sdm_Higher_DDX_XbbVsQCD, sdm_Higher_PN_MD_XbbVsQCD;
  Float_t AK8_ptJet0, PNMD_XbbVsQCD, DDX_XbbVsQCD;//pt, softdropmass and scores of the leading AK8 jet
  Float_t Higher_PN_MD_XbbVsQCD;//BBvsQCD and softdropmass of the leading_PNMD AK8 jet
  Float_t Higher_DDX_XbbVsQCD;//BBvsQCD and softdropmass of the leading_DDX AK8 jet
  Double_t lumiWeight;
  Int_t idx_BestZqq_AK8;
  
  T_Data->SetBranchAddress("AK8_ptJet0", &AK8_ptJet0);
  T_Data->SetBranchAddress("AK8_sdmJet0", &AK8_sdmJet0);
  T_Data->SetBranchAddress("PNMD_XbbVsQCD", &PNMD_XbbVsQCD);
  T_Data->SetBranchAddress("DDX_XbbVsQCD", &DDX_XbbVsQCD);
  T_Data->SetBranchAddress("Higher_PN_MD_XbbVsQCD", &Higher_PN_MD_XbbVsQCD);
  T_Data->SetBranchAddress("sdm_Higher_PN_MD_XbbVsQCD", &sdm_Higher_PN_MD_XbbVsQCD);
  T_Data->SetBranchAddress("Higher_DDX_XbbVsQCD", &Higher_DDX_XbbVsQCD);
  T_Data->SetBranchAddress("sdm_Higher_DDX_XbbVsQCD", &sdm_Higher_DDX_XbbVsQCD);
  T_Data->SetBranchAddress("lumiWeight", &lumiWeight);
  T_Data->SetBranchAddress("idx_BestZqq_AK8", &idx_BestZqq_AK8);

  Int_t nentries = (Int_t)T_Data->GetEntries();

  cout<<nentries<<" entries"<<endl;

  for (int i = 0; i < nentries; i++){      
    T_Data->GetEntry(i); 
    if(i%10000 == 0){cout<<i<<"/"<<nentries<<" entries done"<<endl;}
    //select only Z->bb events
    if(idx_BestZqq_AK8 != -2){
      for(int iPN=0; iPN<PN_DIM;iPN++){
	if(sdm_Higher_PN_MD_XbbVsQCD>l && sdm_Higher_PN_MD_XbbVsQCD<u && Higher_PN_MD_XbbVsQCD>PN_scores.at(iPN) && Higher_PN_MD_XbbVsQCD<=PN_scores.at(iPN+1)){PN_Data_s[iPN]->Fill(sdm_Higher_PN_MD_XbbVsQCD/s, lumiWeight);}
      } 

      for(int iPN=0; iPN<DDX_DIM;iPN++){
	if(sdm_Higher_DDX_XbbVsQCD>l && sdm_Higher_DDX_XbbVsQCD<u && Higher_DDX_XbbVsQCD>DDX_scores.at(iPN) && Higher_DDX_XbbVsQCD<=DDX_scores.at(iPN+1)){DDX_Data_s[iPN]->Fill(sdm_Higher_DDX_XbbVsQCD/s, lumiWeight);}
      }
    }
  }
  //saving histograms in a Tfile
  cout<<"writing root files"<<endl;
  TFile *myfile;
  myfile = TFile::Open("MC_histos_1.root","RECREATE");
  myfile->cd();

  for(int i=0; i<PN_DIM;i++){
    PN_Data_s[i]->Write();
  }

  for(int i=0; i<DDX_DIM;i++){
    DDX_Data_s[i]->Write();
  }
  cout<<"root files completed"<<endl;
  myfile->Close();
  delete myfile;
   
}
