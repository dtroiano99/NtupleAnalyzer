#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>


static const std::vector<int> idx {5,5,4,3};

void syst_wind(){

  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(1001);
  gStyle->SetOptFit(111);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *myfile;
  myfile = TFile::Open("rootfiles/window.root","RECREATE");

  TFile* file;
  file = TFile::Open("./rootfiles/Syst_wind_QCD.root","read");

  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);
  canv->cd();
  TString fname= "./Plots/Window.pdf";
  canv->Print(fname+"[");


  for(int i=0;i<idx.size();i++){
    canv->Clear();
    //take the QCD histogram with the sysystematics
    std::string QCDstring{"QCD"+std::to_string(i)};
    const char * QCDname{QCDstring.c_str()};
    TH1F * hQCD= (TH1F*)file->Get(QCDname);
    std::string Pstring{"QCD"+std::to_string(i)+"_p8"};
    const char * Pname{Pstring.c_str()};
    TH1F * hsystP= (TH1F*)file->Get(Pname);
    std::string Mstring{"QCD"+std::to_string(i)+"_m8"};
    const char * Mname{Mstring.c_str()};
    TH1F * hsystM= (TH1F*)file->Get(Mname);

    //new histograms
    std::string string1{"QCD_systup"+std::to_string(i)};
    const char * name1{string1.c_str()};
    TH1F* QCD_up = (TH1F*)hQCD->Clone(name1);
    std::string string2{"QCD_systdown"+std::to_string(i)};
    const char * name2{string2.c_str()};
    TH1F* QCD_down = (TH1F*)hQCD->Clone(name2);

    for(int k=1; k< hQCD->GetNbinsX()+1;k++){
      float e = hQCD->GetBinContent(k);
      float p = hsystP->GetBinContent(k);
      float m = hsystM->GetBinContent(k);

      float diff1 = abs(e-p);
      float diff2 = abs(e-m);

      float maximum = max(diff1,diff2);

      QCD_up->SetBinContent(k,e+maximum);
      QCD_down->SetBinContent(k,e-maximum);

    }

    auto p1 = new TPad{"p1", "QCD with sistematic",0.05, 0.3, 0.95, 0.95};
    p1->SetBottomMargin(0);
    //p1->SetTicky(2);
    p1->SetLeftMargin(0.105);
    auto p2 = new TPad{"p2", "relative variation", 0.05, 0, 0.95, 0.3};
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0.25);
    p2->SetLeftMargin(0.105);
    //p2->SetTicky(2);
    p1->Draw();
    p2->Draw();

    p1->cd();
    hQCD->GetYaxis()->SetRangeUser(QCD_down->GetMinimum()*0.9,QCD_up->GetMaximum()*1.3);

    hQCD->Draw("hist");

    QCD_up->SetLineStyle(9);   
    QCD_down->SetLineStyle(9);
    QCD_down->SetLineColor(6);
    QCD_up->SetLineColor(3);

    QCD_up->Draw("histsame");
    QCD_down->Draw("histsame");

    myfile->cd();
    hQCD->Write();
    QCD_up->Write();
    QCD_down->Write();

    p2->cd();
    p2->SetGridy();
    gPad->SetGridy();

    TH1F* hQCD_1 = (TH1F*)hQCD->Clone("");
    TH1F* QCD_up_1 = (TH1F*)QCD_up->Clone("");
    TH1F* QCD_down_1 = (TH1F*)QCD_down->Clone("");
  
    hQCD_1->Divide(hQCD);
    QCD_up_1->Divide(hQCD);
    QCD_down_1->Divide(hQCD);

    hQCD_1->SetTitle("");
    hQCD_1->GetYaxis()->SetRangeUser(0.8,1.2);
    hQCD_1->GetYaxis()->SetNdivisions(9);
    hQCD_1->GetYaxis()->SetTitle("QCD #pm syst./ QCD");
    hQCD_1->GetXaxis()->SetTitle("x");
    //hQCD_1->GetYaxis()->SetTitleSize(0.1);
    hQCD_1->GetXaxis()->SetTitleSize(0.08);

    hQCD_1->Draw("hist");
    QCD_up_1->Draw("histsame");
    QCD_down_1->Draw("histsame");
    p2->RedrawAxis();
    canv->ForceUpdate();

    canv->Print(fname);
   


  }

  canv->Print(fname+"]");

}
