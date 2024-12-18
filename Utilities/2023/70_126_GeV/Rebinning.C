#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>

void Rebinning(){

  gStyle->SetOptFit(1);

  TH1::AddDirectory(kFALSE); 
  
  TFile *f;
  f = TFile::Open("Data_histos.root","read");

  static const std::vector<float> PN_scores {0.588, 0.856, 0.954,  0.988,  1};
  int PN_DIM = PN_scores.size()-1;
  static const std::vector<int> rebinning {2,4,8};

  for(int ir=0; ir<3;ir++){

  TH1F *PN_Data[PN_DIM];
  
  for(int i=0; i<PN_DIM;i++){
    std::string nomestringa{"PNhist"+std::to_string(i)};
    const char * nome{nomestringa.c_str()};    
    PN_Data[i] = (TH1F*)f->Get(nome);
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
    PN_Data[i]->GetXaxis()->SetTitle("m_{SD} / 200 GeV");
    PN_Data[i]->Rebin(rebinning.at(ir));
    std::string Y_string{"Events / "+ std::to_string(PN_Data[i]->GetBinWidth(1))};
    const char * Y_title{Y_string.c_str()};
    PN_Data[i]->GetYaxis()->SetTitle(Y_title);
  }

  //std::string fname{"./Plots/Data_sb_"+std::to_string(rebinning.at(ir))+"GeV.pdf"};
  //const char * fname1{nomestringa.c_str()};
  TString fname1= "./Plots/Data_sb"+std::to_string(rebinning.at(ir))+"GeV.pdf";
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
  
  }
  f->Close();
  delete f;

}
