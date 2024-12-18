#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));

static const std::vector<float> PN_scores0 {0.641, 0.791, 0.875, 0.926, 0.957, 0.977, 0.988, 0.995, 1};

static const std::vector<float> PN_scores {0.641, 0.875, 0.957, 0.988, 1};

//static const std::vector<float> PN_scores0 {0.641, 0.791, 0.875};

//static const std::vector<float> PN_scores {0.641, 0.875};

int PN_DIM =  PN_scores.size()-1; //number of PN-MD score bins

void histos(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //gStyle->SetGridStyle(3);
  //gStyle->SetGridWidth(3);

  //gStyle->SetOptStat(1000000001);
  //gStyle->SetOptFit(111);

  TH1::AddDirectory(kFALSE); 
  //histogams signal region

  

  //Data

  TFile *f1;
  f1 = TFile::Open("Data_histos.root","read");

  //PN-MD softdropmass histograms
  TH1F *PN_Data[PN_DIM];
  
  for(int iPN=0; iPN<PN_scores0.size()-1;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};
    if(iPN%2 ==0){
      int i = (int)iPN/2;
      PN_Data[i] = (TH1F*)f1->Get(nome);
      float start = floor(PN_scores.at(i) * 10000) / 10000;
      std::stringstream ss1;
      ss1 <<start;
      std::string str1 = ss1.str();
      float end = floor(PN_scores.at(i+1) * 10000) / 10000;
      std::stringstream ss2;
      ss2 <<end;
      std::string str2 = ss2.str();
      std::string nomestringa_title{str1+" < PN-MD_BBvsQCD #leq "+ str2};
      const char * title{nomestringa_title.c_str()};
      PN_Data[i]->SetTitle(title);
    }
    else{
      int i = (int)(iPN-1)/2;
      TH1F * h1 = (TH1F*)f1->Get(nome);
      PN_Data[i]->Add(h1);
    }
   
  } 
  
  f1->Close();
  delete f1;

  
   //draw raw histograms
 /* 
  TString fname1= "./Plots/Hist.pdf";
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);  
  canv->cd();
  canv->Print(fname1+"[");
  for(int i=0;i<PN_DIM;i++){       
    canv->cd();
    gPad->SetLeftMargin(0.15);
    PN_Data[i]->Draw("hist");
    canv->Print(fname1);
    canv->Clear();
  }
  canv->Print(fname1+"]"); 
   
*/
  
  
  //MC

  TFile *f2;
  f2 = TFile::Open("MC_ZtoQQ_histos.root","read");

  //PN-MD softdropmass histograms
  TH1F *PN_MC[PN_DIM];
  
  for(int iPN=0; iPN<PN_scores0.size()-1;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};
    if(iPN%2 ==0){
      int i = (int)iPN/2;
      PN_MC[i] = (TH1F*)f2->Get(nome);
      float start = floor(PN_scores.at(i) * 10000) / 10000;
      std::stringstream ss1;
      ss1 <<start;
      std::string str1 = ss1.str();
      float end = floor(PN_scores.at(i+1) * 10000) / 10000;
      std::stringstream ss2;
      ss2 <<end;
      std::string str2 = ss2.str();
      std::string nomestringa_title{str1+" < PN-MD_BBvsQCD #leq "+ str2};
      const char * title{nomestringa_title.c_str()};
      PN_MC[i]->SetTitle(title);
    }
    else{
      int i = (int)(iPN-1)/2;
      TH1F * h1 = (TH1F*)f2->Get(nome);
      PN_MC[i]->Add(h1);
    }
    
  }
  
  f2->Close();
  delete f2;

  TFile *myFile;
  myFile = TFile::Open("./rootfiles/Histos.root","RECREATE");
  for(int iPN=0; iPN<PN_scores.size()-1;iPN++){
    myFile->cd();
    PN_Data[iPN]->Write();
    PN_MC[iPN]->Write();
  }
  myFile->Close();
  delete myFile;
}
