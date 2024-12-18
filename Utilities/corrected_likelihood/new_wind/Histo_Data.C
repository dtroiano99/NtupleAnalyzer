#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));


double Start=48;
double End=152;
int bin=52;
int bin1 = 16;
double l=80;
double u=112;
double s = 250;

//static const std::vector<float> PN_scores {0.641, 0.875};
static const std::vector<float> PN_scores {0.641, 0.791, 0.875, 0.926, 0.957, 0.977, 0.988, 0.995, 1};

int PN_DIM =  PN_scores.size()-1; //number of PN-MD score bins

void Histo_Data(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //gStyle->SetOptStat(1000000001);
  //gStyle->SetOptFit(1111);

  TH1::AddDirectory(kFALSE); 
  
  std::string Y_string{"Events / "+ std::to_string(((End/s-Start/s)/bin))};
  const char * Y_title{Y_string.c_str()};

  //PN-MD softdropmass histograms

  //histogams sidebands
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

  //histogams signal region
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

  TChain *T_Data = new TChain("passedEvents");  // name of the tree is the argument
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022C_corrected.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022D_corrected.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022E_corrected.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022F0_corrected.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022F1_corrected.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022F2_corrected.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022G_corrected.root");


  //definition of the variables and branches
  Double_t AK8_sdmJet0; Float_t PNMD_XbbVsQCD;//softdropmass and score of the leading AK8 jet

  

  bool ispreEE;

  ROOT::VecOps::RVec<float> *AK4_eta = 0;
  ROOT::VecOps::RVec<float> *AK4_phi = 0;
  ROOT::VecOps::RVec<bool> *AK4_isloose = 0;
  
  T_Data->SetBranchAddress("AK8_sdmJet0", &AK8_sdmJet0);
  T_Data->SetBranchAddress("PNMD_XbbVsQCD", &PNMD_XbbVsQCD);
  
  T_Data->SetBranchAddress("ispreEE", &ispreEE);
  
  T_Data->SetBranchAddress("AK4PuppiJets_eta",&AK4_eta);
  T_Data->SetBranchAddress("AK4PuppiJets_phi",&AK4_phi);
  T_Data->SetBranchAddress("AK4PuppiJets_isloose",&AK4_isloose);
  



    
  Int_t nentries = (Int_t)T_Data->GetEntries();

  cout<<nentries<<" entries"<<endl;
  for (int i = 0; i < nentries; i++){      
    T_Data->GetEntry(i);
    if(i%10000 == 0){cout<<i<<"/"<<nentries<<" entries done"<<endl;} 

    for(int iPN=0; iPN<PN_DIM;iPN++){
    
      if(PNMD_XbbVsQCD>PN_scores.at(iPN) && PNMD_XbbVsQCD<=PN_scores.at(iPN+1)){

	if(AK8_sdmJet0>l && AK8_sdmJet0<u ){

	  bool jetvetomap_ = false;
	  
	  TFile *f1;
	  if(ispreEE){
	    f1 = TFile::Open("/lustre/cms/store/user/dtroiano/Commissioning/Utilities/Summer22_23Sep2023_RunCD_v1.root","read");
	  }
	  else{
	    f1 = TFile::Open("/lustre/cms/store/user/dtroiano/Commissioning/Utilities/Summer22EE_23Sep2023_RunEFG_v1.root","read");
	  }
	  
	  TH2D *h_jetvetomap = (TH2D*)f1->Get("jetvetomap");
	  
	  
	  for (UInt_t u4 = 0; u4 < AK4_eta->size(); ++u4) {
	    //cout<<AK4_isloose->at(u4)<<endl;
	    if(AK4_isloose->at(u4) ){
	      double veto = h_jetvetomap->GetBinContent(h_jetvetomap->FindBin(AK4_eta->at(u4), AK4_phi->at(u4)));
	      if(veto>0){
		jetvetomap_ = true;
		break;
	      }
	    }
	  }
	  
	  delete h_jetvetomap;
	  f1->Close();
	  delete f1;
	  
	  if(jetvetomap_ == false){PN_Data_s[iPN]->Fill(AK8_sdmJet0/s);}
	  
	}
	if((AK8_sdmJet0>Start && AK8_sdmJet0<l) || (AK8_sdmJet0 >u && AK8_sdmJet0<End)){
	  
	  bool jetvetomap_ = false;
	  
	  TFile *f1;
	  if(ispreEE){
	    f1 = TFile::Open("/lustre/cms/store/user/dtroiano/Commissioning/Utilities/Summer22_23Sep2023_RunCD_v1.root","read");
	  }
	  else{
	    f1 = TFile::Open("/lustre/cms/store/user/dtroiano/Commissioning/Utilities/Summer22EE_23Sep2023_RunEFG_v1.root","read");
	  }
	  
	  TH2D *h_jetvetomap = (TH2D*)f1->Get("jetvetomap");
	  
	  
	  for (UInt_t u4 = 0; u4 < AK4_eta->size(); ++u4) {
	    //cout<<AK4_isloose->at(u4)<<endl;
	    if(AK4_isloose->at(u4) ){
	      double veto = h_jetvetomap->GetBinContent(h_jetvetomap->FindBin(AK4_eta->at(u4), AK4_phi->at(u4)));
	      if(veto>0){
		jetvetomap_ = true;
		break;
	      }
	    }
	  }
	  
	  delete h_jetvetomap;
	  f1->Close();
	  delete f1;
	  
	  if(jetvetomap_ == false){PN_Data[iPN]->Fill(AK8_sdmJet0/s);}
	} 
      }
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
  cout<<"root file completed"<<endl;
  
  myfile->Close();
  delete myfile;
  
}
