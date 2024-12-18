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

void Data_W_ratio(){
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);

  //gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1);

  TH1::AddDirectory(kFALSE);

  TFile *f1;
  f1 = TFile::Open("Data_histos.root","read");

  TFile *f2;
  f2 = TFile::Open("Wjets_histos.root","read");

  TString fname1= "./Plots/Data_W_ratio.pdf";
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

    std::string nomestringa{"PNhist"+std::to_string(i)};
    const char * nome{nomestringa.c_str()};

    TH1F *Data = (TH1F*)f1->Get(nome);
    TH1F *Wjets = (TH1F*)f2->Get(nome);
    Data->Rebin(rebinning.at(i));
    Wjets->Rebin(rebinning.at(i));

    Data->SetLineColor(1);
    Data->SetMarkerStyle(8);

    Wjets->SetFillColor(41);
    Wjets->SetLineColor(41);

    //Plot part

    Wjets->SetTitle("");

    canv1->cd();
    gPad->SetLeftMargin(0.15);
    canv1->Clear();

    auto hRatio= (TH1F*)Wjets->Clone("");
    hRatio->GetXaxis()->SetLabelSize(0.07);
    hRatio->GetXaxis()->SetTitleSize(0.06);
    hRatio->GetYaxis()->SetLabelSize(0.07);
    hRatio->GetYaxis()->SetTitleSize(0.06);

     auto p1 = new TPad{"p1", "QCD vs qcd", 0.05, 0.25, 0.95, 0.95};//0.05, 0.3, 0.95, 0.95
    p1->SetBottomMargin(0); // Tolgo il margine bianco attorno al pad
    p1->SetTopMargin(0.1);
    //p1->SetTicky(2);
    p1->SetLeftMargin(0.15);
    auto p2 = new TPad{"p2", "Ratio", 0.05, 0, 0.95, 0.25};//0.05, 0, 0.95, 0.3
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0.4);
    p2->SetLeftMargin(0.15);
    //p2->SetTicky(2);
    p1->Draw();
    p2->Draw();


    p1->cd();
    p1->SetLogy();
    float min = Wjets->GetMinimum()*0.8;
    Wjets->SetMinimum(min+0.2);
    float max = Data->GetMaximum()*1.55;
    Wjets->SetMaximum(max);
    Wjets->GetYaxis()->SetLabelSize(0.06);
    Wjets->GetYaxis()->SetTitleSize(0.06);
    Wjets->GetYaxis()->SetTitleOffset(1.1);  
    std::string Ystring{"Events / "+std::to_string(Wjets->GetBinWidth(1))};
    const char * Yname{Ystring.c_str()};
    Wjets->GetYaxis()->SetTitle(Yname);
    Wjets->Draw("hist");
    Data->Draw("pesame");  

    auto legend1 = new TLegend(0.36,0.4,0.53,0.65);
    legend1->AddEntry(Data,"Data","lpe");
    legend1->AddEntry(Wjets,"W #rightarrow q' #bar{q} + jets","f");
    legend1->SetBorderSize(0);
    legend1->Draw("");

    p2->cd();
    p2->SetGridy();
    gPad->SetGridy();

    hRatio->Divide(Data);
    hRatio->SetMarkerStyle(8);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);
    hRatio->SetTitle("");

    hRatio->GetXaxis()->SetLabelSize(0.16);
    hRatio->GetXaxis()->SetTitleSize(0.16);
    hRatio->GetYaxis()->SetTitleSize(0.16);
    hRatio->GetYaxis()->SetLabelSize(0.14);
    hRatio->GetYaxis()->SetTitle("W/data ratio");
    hRatio->GetYaxis()->SetTitleOffset(0.42);
    //hRatio->GetYaxis()->SetRangeUser(-1.1,1.1);
    hRatio->GetYaxis()->SetNdivisions(6);
    hRatio->Draw("lpe");

    int iPeriod = 5;
    int iPos = 0;
    CMS_lumi( p1, iPeriod, iPos );



    canv1->Print(fname1);
    

  }//close score loop

  canv1->Print(fname1+"]");

}
