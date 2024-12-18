#include "/lustrehome/dtroiano1/tdrstyle.C"
#include "/lustrehome/dtroiano1/CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"

void Fit_plot(){
  //setTDRStyle();

  int W = 800;
  int H = 600;
  int H_ref = 600; 
  int W_ref = 800; 


  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  TCanvas* canv = new TCanvas("C","C",50,50,800,600);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);

  TFile* file;
  file = TFile::Open("Data_histos.root","read");
  TH1F* h = (TH1F*)file->Get("PNhist3");  
  h->Rebin(8);
  h->SetTitle("");

  int histLineColor = 1;
  int histFillColor = 1;
  float markerSize  = 5.0;

  
  TFile* file1;
  file1 = TFile::Open("rootfiles/Fit_Cebysev20.root","read");
  TF1* f = (TF1*)file1->Get("3parameters/fun3");
  TFitResult* r = (TFitResult*) file1->Get("3parameters/TFitResult-PNhist3-fun3");
  r->Print();
  
  int n_ = 2;
    
  float x1_l = 0.92;
  float y1_l = 0.60;
  
  float dx_l = 0.30;
  float dy_l = 0.18;
  float x0_l = x1_l-dx_l;
  float y0_l = y1_l-dy_l;

  h->SetDirectory(0);
  h->SetLineColor(histLineColor);
  h->SetFillColor(histFillColor);
  h->SetMarkerStyle(8);
  h->SetMarkerSize(1.3);

  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);


  h->Draw("pe");
  f->Draw("same");
  
  auto legend = new TLegend(0.28,0.7,0.55,0.8);
  legend->SetBorderSize(0);
  legend->AddEntry(h,"2022 data","lep");
  legend->Draw("");

  h->SetMinimum(100);
  h->SetMaximum(550);
  auto l1 = new TLine(126,h->GetMinimum(),126,h->GetMaximum());
  auto l2 = new TLine(70,h->GetMinimum(),70,h->GetMaximum());
  l1->SetLineColor(2);
  l2->SetLineColor(2);
  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  l1->Draw();
  l2->Draw();
  int iPeriod = 5;
  int iPos = 0;

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.SetTextAlign(13);  //align at top
  //latex.DrawLatex(.2,.9,"K_{S}");
  //latex.DrawLatex(.3,.9,"K^{*0}");

  std::string cstring{"#chi^{2} / ndf = "+std::to_string(r->Chi2())+" / "+std::to_string(r->Ndf())};
  const char * cname{cstring.c_str()};
  latex.DrawLatex(.7,.9,cname);
  for(int yyy=0;yyy<4;yyy++){
    std::string string{"c_{"+std::to_string(yyy)+"} = "+std::to_string(r->Parameter(yyy))+" #pm "+std::to_string(r->ParError(yyy))};
    const char * name{string.c_str()};
    latex.DrawLatex(.7,.9-(0.05*(yyy+1)),name);
  }
   
  CMS_lumi( canv, iPeriod, iPos );  

     

}
