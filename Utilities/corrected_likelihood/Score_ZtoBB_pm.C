#include "/lustrehome/dtroiano1/tdrstyle.C"
#include "/lustrehome/dtroiano1/CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"

void Score_ZtoBB_pm(){
  //setTDRStyle();

  static const std::vector<float> PN_scores {0.641, 0.875,0.957, 0.988,  1};

  int PN_DIM =  PN_scores.size()-1;

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
  auto fcolor = new TColor;
  Color_t color = fcolor->GetColor("#5790fc");
  h->SetLineColor(color);
  h->SetFillColor(color);
  //h->SetLineColor(histLineColor);
  //h->SetFillColor(histFillColor);
  h->SetFillStyle(1001);

  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitle("PNet-MD_{bbvsQCD}");

   h->SetMinimum(0);
  h->SetMaximum(120);

  auto legend = new TLegend(0.3,0.5,0.58,0.8);
  legend->SetBorderSize(0);
  for(int j=0;j< PN_DIM;j++){    
    TH1F* h1 = (TH1F*)h->Clone("");
    std::string hstring{"h"+std::to_string(j)};
    const char * hname{hstring.c_str()};
    h1->SetName(hname);
    h1->SetFillColorAlpha(color, 1.25-1/j);
    for(int k=1; k< h1->GetNbinsX()+1;k++){
      float xData = h1->GetXaxis()->GetBinCenter(k);
      float ent =  h1->GetBinContent(k);
      if (xData>PN_scores.at(j) && xData<=PN_scores.at(j+1)){
	h1->SetBinContent(k,ent);
      }
      else{
	h1->SetBinContent(k,0);
      }
    }
    h1->GetXaxis()->SetRangeUser(0.641,1);
    if(j==0){
      auto fcolor1 = new TColor;
      Color_t color1 = fcolor1->GetColor("#5790fc");
      h1->SetLineColor(color1);
      h1->SetFillColor(color1);
      legend->AddEntry(h1,"4^{th} highest score  region","f");
    }
    if(j==1){
      auto fcolor1 = new TColor;
      Color_t color1 = fcolor1->GetColor("#f89c20");
      h1->SetLineColor(color1);
      h1->SetFillColor(color1);
      legend->AddEntry(h1,"3^{rd} highest score  region","f");
    }
    if(j==2){
      auto fcolor1 = new TColor;
      Color_t color1 = fcolor1->GetColor("#e42536");
      h1->SetLineColor(color1);
      h1->SetFillColor(color1);
      legend->AddEntry(h1,"2^{nd} highest score  region","f");
    }
    if(j==3){
      auto fcolor1 = new TColor;
      Color_t color1 = fcolor1->GetColor("#964a8b");
      h1->SetLineColor(color1);
      h1->SetFillColor(color1);
      legend->AddEntry(h1,"highest score  region","f");
    }
    if(j==0){
      h1->Draw("hist");
      TLatex tex; tex.SetNDC();
      //tex = new TLatex(40,1800,"Z #rightarrow b#bar{b} ");
      tex.SetTextColor(1);
      tex.SetTextSize(0.09);
      tex.DrawLatex(0.2,0.3,"Z #rightarrow b#bar{b} ");

    }
    else{h1->Draw("histsamea");}
    
  }
  gPad->Update();
   gPad->RedrawAxis();
  
  //legend->AddEntry("h1","Z #rightarrow q#bar{q} (MC)","f");
  legend->Draw("");
  canv->Update();

 
  auto l1 = new TLine(0.875,h->GetMinimum(),0.875,h->GetMaximum());
  auto l2 = new TLine(0.957,h->GetMinimum(),0.957,h->GetMaximum());
  auto l3 = new TLine(0.988,h->GetMinimum(),0.988,h->GetMaximum());
  l1->SetLineColor(15);
  l2->SetLineColor(15);
  l3->SetLineColor(15);
  l1->SetLineWidth(1);
  l2->SetLineWidth(1);
  l3->SetLineWidth(1);
  l1->SetLineStyle(3);
  l2->SetLineStyle(3);
  l3->SetLineStyle(3);
  l1->Draw();
  l2->Draw();
  l3->Draw();

  int iPeriod = 6;
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
