#include "/lustrehome/dtroiano1/tdrstyle.C"
#include "/lustrehome/dtroiano1/CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"

void Score_ZtoBB_plot(){
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
  file = TFile::Open("Score_ZtoBB_hist.root","read");
  TH1F* h = (TH1F*)file->Get("Lpt_bAK8_PNMD_BBvsQCD");  

  int histLineColor = 4;
  int histFillColor = 4;
  float markerSize  = 1.0;

  
  
  
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
  h->SetFillStyle(1001);

  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitle("PNet-MD_{bbvsQCD}");

  h->Draw("hist");
  
  auto legend = new TLegend(0.3,0.7,0.58,0.8);
  legend->SetBorderSize(0);
  legend->AddEntry(h,"Z #rightarrow q#bar{q} (MC)","f");
  legend->Draw("");

  h->SetMinimum(0);
  h->SetMaximum(120);
  auto l1 = new TLine(0.875,h->GetMinimum(),0.875,h->GetMaximum());
  auto l2 = new TLine(0.957,h->GetMinimum(),0.957,h->GetMaximum());
  auto l3 = new TLine(0.988,h->GetMinimum(),0.988,h->GetMaximum());
  l1->SetLineColor(2);
  l2->SetLineColor(2);
  l3->SetLineColor(2);
  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  l3->SetLineWidth(2);
  l1->Draw();
  l2->Draw();
  l3->Draw();
  int iPeriod = 5;
  int iPos = 0;
  CMS_lumi( canv, iPeriod, iPos );

  /*
  TLatex *   tex = new TLatex(40,1800,"sideband");
  tex->SetTextColor(2);
  tex->SetTextSize(0.05);
  tex->Draw();

  TLatex *   tex1 = new TLatex(100,1800,"sideband");
  tex1->SetTextColor(2);
  tex1->SetTextSize(0.05);
  tex1->Draw();

  TLatex *   tex2 = new TLatex(100,1000,"signal region");
  tex2->SetTextColor(2);
  tex2->SetTextSize(0.05);
  tex2->Draw();
  */
  

}
