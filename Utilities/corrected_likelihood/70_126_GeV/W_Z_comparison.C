#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include "/lustrehome/dtroiano1/tdrstyle.C"
#include "/lustrehome/dtroiano1/CMS_lumi.C"


int bin1 = 56;
double l=70;
double u=126;
double s = 200;

static const std::vector<int> rebinning {4,8,8,8};

void W_Z_comparison(){
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);

  //gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1);

  TH1::AddDirectory(kFALSE);

  TFile *f1;
  f1 = TFile::Open("MC_ZtoQQ_histos.root","read");

  TFile *f2;
  f2 = TFile::Open("Wjets_histos.root","read");
  
  TString fname1= "./Plots/W_Z_comparison.pdf";
  //auto canv1 = new TCanvas("Canvas1", "Canvas1", 10000, 6000);
  int W = 800;
  int H = 600;
  int H_ref = 600;
  int W_ref = 800;


  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  TCanvas* canv1 = new TCanvas("C","C",50,50,800,600);
  canv1->SetFillColor(0);
  canv1->SetBorderMode(0);
  canv1->SetFrameFillStyle(0);
  canv1->SetFrameBorderMode(0);
  canv1->SetLeftMargin( L/W );
  canv1->SetRightMargin( R/W );
  canv1->SetTopMargin( 0 );
  canv1->SetBottomMargin( B/H );
  canv1->SetTickx(0);
  canv1->SetTicky(0);
  canv1->cd();
  canv1->Print(fname1+"[");  

  for(int i=0;i<rebinning.size();i++){

    std::string nomestringa{"PNhist"+std::to_string(i)+"_s"};
    const char * nome{nomestringa.c_str()};

    TH1F *Zjets = (TH1F*)f1->Get(nome);
    TH1F *Wjets = (TH1F*)f2->Get(nome);
    Zjets->Rebin(rebinning.at(i));
    Wjets->Rebin(rebinning.at(i));

    Wjets->SetFillColor(2);
    Wjets->SetLineColor(2);
    Wjets->SetFillStyle(3335);

    Zjets->SetFillColor(4);
    Zjets->SetLineColor(4);
    Zjets->SetFillStyle(3353);
    
    //Plot part

    Zjets->SetTitle("");
    
    canv1->cd();
    gPad->SetLeftMargin(0.15);
    canv1->Clear();    
    

    
    float min = Wjets->GetMinimum()*0.8;
    Zjets->SetMinimum(0);
    float max = Zjets->GetMaximum()*1.25;
    if(i <1){max = Wjets->GetMaximum()*1.1;}
    Zjets->SetMaximum(max);
        
    Zjets->GetXaxis()->SetTitle("m_{SD} / 200 GeV");
    Zjets->GetYaxis()->SetLabelSize(0.05);
    Zjets->GetYaxis()->SetTitleSize(0.06);
    Zjets->GetXaxis()->SetLabelSize(0.05);
    Zjets->GetXaxis()->SetTitleSize(0.05);
    Zjets->GetYaxis()->SetTitleOffset(1.3);
    canv1->cd();
    canv1->SetTopMargin(0.1);
    std::string Ystring{"Events / "+std::to_string(Wjets->GetBinWidth(1))};
    const char * Yname{Ystring.c_str()};
    Zjets->GetYaxis()->SetTitle(Yname);
    Zjets->Draw("hist");
    Wjets->Draw("histsame");
    auto legend1 = new TLegend(0.65,0.65,0.89,0.89);
    legend1->AddEntry(Wjets,"W #rightarrow q' #bar{q}","f");
    legend1->AddEntry(Zjets,"Z #rightarrow q #bar{q}","f");
    legend1->SetBorderSize(0);
    legend1->Draw("");

    int iPeriod = 5;
    int iPos = 0;
    CMS_lumi( canv1, iPeriod, iPos );

    
    
    canv1->Print(fname1);
    
  }//close score loop

  canv1->Print(fname1+"]");

}

  
