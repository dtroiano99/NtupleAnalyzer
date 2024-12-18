#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));

static const std::vector<float> PN_scores {0.588, 0.856, 0.954,  0.988,  1};

static const std::vector<int> rebinning {8,8,8,8};

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
  f2 = TFile::Open("MC_Z2_histos.root","read");

  //PN-MD softdropmass histograms
  TH1F *PN_MC[PN_DIM];
  TH1F *PN_MCcc[PN_DIM];
  TH1F *PN_MCll[PN_DIM];
  
  for(int iPN=0; iPN<PN_scores.size()-1;iPN++){
    std::string nomestringa{"Zbbhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};
    PN_MC[iPN] = (TH1F*)f2->Get(nome);
    PN_MC[iPN]->Rebin(rebinning.at(iPN));
 
    std::string nomestringacc{"Zcchist"+std::to_string(iPN)+"_s"};
    const char * nomecc{nomestringacc.c_str()};
    PN_MCcc[iPN] = (TH1F*)f2->Get(nomecc);
    PN_MCcc[iPN]->Rebin(rebinning.at(iPN)); 

    std::string nomestringall{"Zllhist"+std::to_string(iPN)+"_s"};
    const char * nomell{nomestringall.c_str()};
    PN_MCll[iPN] = (TH1F*)f2->Get(nomell);
    PN_MCll[iPN]->Rebin(rebinning.at(iPN)); 
  }
  
  f2->Close();
  delete f2;

  TFile *myFile;
  myFile = TFile::Open("./rootfiles/Histos2.root","RECREATE");
  for(int iPN=0; iPN<PN_scores.size()-1;iPN++){
    myFile->cd();
    PN_Data[iPN]->Write();
    PN_MC[iPN]->Write();
    PN_MCcc[iPN]->Write();
    PN_MCll[iPN]->Write();
  }
  myFile->Close();
  delete myFile;
}
