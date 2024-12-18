#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));

static const std::vector<float> PN_scores {0.641, 0.875, 0.957, 0.988, 1};

static const std::vector<int> rebinning {4,8,8,8};

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
  
  for(int iPN=0; iPN<PN_scores.size()-1;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};    
    PN_Data[iPN] = (TH1F*)f1->Get(nome);    
    PN_Data[iPN]->Rebin(rebinning.at(iPN));
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
  
  for(int iPN=0; iPN<PN_scores.size()-1;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};
    PN_MC[iPN] = (TH1F*)f2->Get(nome);
    PN_MC[iPN]->Rebin(rebinning.at(iPN));    
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
