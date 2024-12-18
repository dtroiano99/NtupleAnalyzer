#include "/lustrehome/dtroiano1/tdrstyle.C"
#include "/lustrehome/dtroiano1/CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"

void Large_msd_plot(){
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
  file = TFile::Open("Large_ZtoQQ_histos.root","read");
  TH1F* h = (TH1F*)file->Get("PNhist0_s");  

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
  auto fcolor = new TColor;
  Color_t color = fcolor->GetColor("#5790fc");
  h->SetLineColor(color);
  h->SetFillColor(color);
  //Int_t ci = TColor::GetFreeColorIndex();
  //TColor *color = new TColor(ci, 0.1, 0.2, 0.3);	
  //h->SetFillStyle(3353);

  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);


  h->Draw("hist");
  
  auto legend = new TLegend(0.62,0.5,0.88,0.6);
  legend->SetBorderSize(0);
  legend->AddEntry(h,"Z #rightarrow q#bar{q}","f");
  legend->Draw("");

  h->SetMinimum(200);
  h->SetMaximum(2600);
  auto l1 = new TLine(126,h->GetMinimum(),126,h->GetMaximum());
  auto l2 = new TLine(70,h->GetMinimum(),70,h->GetMaximum());
  l1->SetLineColor(2);
  l2->SetLineColor(2);
  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  l1->Draw();
  l2->Draw();
  int iPeriod = 6;
  int iPos = 0;
  CMS_lumi( canv, iPeriod, iPos );

  TLatex *   tex = new TLatex(40,1800,"sideband");
  //TLatex tex;
  //tex.SetTextAlign(13);
  //tex.DrawLatex(.15,.82,"sideband");
  tex->SetTextColor(2);
  tex->SetTextSize(0.05);
  tex->Draw();

  TLatex *   tex1 = new TLatex(100,1800,"sideband");
  //TLatex tex1;
  //tex1.SetTextAlign(13);
  //tex1.DrawLatex(.75,.82,"sideband");
  tex1->SetTextColor(2);
  tex1->SetTextSize(0.05);
  tex1->Draw();

  TLatex *   tex2 = new TLatex(100,1000,"signal region");
  //TLatex tex2;
  //tex2.SetTextAlign(13);
  //tex1.DrawLatex(.45,.82,"signal region");
  tex2->SetTextColor(2);
  tex2->SetTextSize(0.05);
  tex2->Draw();

  

}
