#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>


static const std::vector<int> idx {3,3,3,3};

void QCD_and_syst(){

  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(1001);
  gStyle->SetOptFit(111);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *myfile;
  myfile = TFile::Open("rootfiles/QCD.root","RECREATE");
  
  TFile* file;
  file = TFile::Open("./rootfiles/Syst_QCD.root","read");

  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);
  canv->cd();
  TString fname= "./Plots/QCD.pdf";
  canv->Print(fname+"[");

  for(int i=0;i<idx.size();i++){
    //take the QCD histogram with the sysystematics
    std::string QCDstring{"QCD"+std::to_string(i)};
    const char * QCDname{QCDstring.c_str()};
    TH1F * hQCD= (TH1F*)file->Get(QCDname);
    std::string RMSstring{"RMS"+std::to_string(i)};
    const char * RMSname{RMSstring.c_str()};
    TH1F * hsystRMS= (TH1F*)file->Get(RMSname);
    std::string Fitstring{"Fit"+std::to_string(i)};
    const char * Fitname{Fitstring.c_str()};    
    TH1F * hsystFit= (TH1F*)file->Get(Fitname);

    //new histograms 
    std::string string1{"QCD_systup_rms"+std::to_string(i)};
    const char * name1{string1.c_str()};
    TH1F* QCD_uprms = (TH1F*)hQCD->Clone(name1);
    std::string string2{"QCD_systdown_rms"+std::to_string(i)};
    const char * name2{string2.c_str()};
    TH1F* QCD_downrms = (TH1F*)hQCD->Clone(name2);
    std::string string3{"QCD_systup_fit"+std::to_string(i)};
    const char * name3{string3.c_str()};
    TH1F* QCD_upfit = (TH1F*)hQCD->Clone(name3);
    std::string string4{"QCD_systdown_fit"+std::to_string(i)};
    const char * name4{string4.c_str()};
    TH1F* QCD_downfit = (TH1F*)hQCD->Clone(name4);

    QCD_uprms->Add(hsystRMS,1);
    QCD_downrms->Add(hsystRMS,-1);
    QCD_upfit->Add(hsystFit,1);
    QCD_downfit->Add(hsystFit,-1);

    float max = QCD_uprms->GetMaximum();
    if(max<QCD_upfit->GetMaximum()){max=QCD_upfit->GetMaximum();}

    float min = QCD_downrms->GetMinimum();
    if(min<QCD_downfit->GetMinimum()){min=QCD_downfit->GetMinimum();}

    auto p1 = new TPad{"p1", "QCD with sistematic",0.05, 0.3, 0.95, 0.95};
    p1->SetBottomMargin(0);
    p1->SetTicky(2);
    p1->SetLeftMargin(0.105);
    auto p2 = new TPad{"p2", "relative variation", 0.05, 0, 0.95, 0.3};
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0.25);
    p2->SetLeftMargin(0.105);
    p2->SetTicky(2);
    p1->Draw();
    p2->Draw();

    p1->cd();
    hQCD->GetYaxis()->SetRangeUser(min,max);
   
    hQCD->Draw("hist");

    QCD_uprms->SetLineStyle(9);
    QCD_uprms->SetLineColor(2);
    QCD_downrms->SetLineStyle(9);
    QCD_downrms->SetLineColor(6);
    QCD_upfit->SetLineStyle(9);
    QCD_upfit->SetLineColor(3);
    QCD_downfit->SetLineStyle(9);
    QCD_downfit->SetLineColor(8);
    
    
    

    

    QCD_uprms->Draw("histsame");
    QCD_downrms->Draw("histsame");
    QCD_upfit->Draw("histsame");
    QCD_downfit->Draw("histsame");

    myfile->cd();
    hQCD->Write();
    QCD_uprms->Write();
    QCD_downrms->Write();
    QCD_upfit->Write();
    QCD_downfit->Write();

    auto legend1 = new TLegend(0.69,0.58,0.8,0.89);
    legend1->AddEntry(hQCD,"nominal QCD","l");
    legend1->AddEntry(QCD_uprms,"QCD systup_rms","l");
    legend1->AddEntry(QCD_downrms,"QCD systdown_rms","l");
    legend1->AddEntry(QCD_upfit,"QCD systup_fit","l");
    legend1->AddEntry(QCD_downfit,"QCD systdown_fit","l");
    legend1->SetBorderSize(0);

    if(i==0){legend1->Draw();}
    
    p2->cd();
    p2->SetGridy();
    gPad->SetGridy();

    TH1F* hQCD_1 = (TH1F*)hQCD->Clone("");
    TH1F* QCD_uprms_1 = (TH1F*)QCD_uprms->Clone("");
    TH1F* QCD_downrms_1 = (TH1F*)QCD_downrms->Clone("");
    TH1F* QCD_upfit_1 = (TH1F*)QCD_upfit->Clone("");
    TH1F* QCD_downfit_1 = (TH1F*)QCD_downfit->Clone("");

    hQCD_1->Divide(hQCD);
    QCD_uprms_1->Divide(hQCD);
    QCD_downrms_1->Divide(hQCD);
    QCD_upfit_1->Divide(hQCD);
    QCD_downfit_1->Divide(hQCD);

    
    hQCD_1->SetTitle("");
    hQCD_1->GetYaxis()->SetRangeUser(0.985,1.015);
    hQCD_1->GetYaxis()->SetNdivisions(9);
    hQCD_1->GetYaxis()->SetTitle("QCD #pm syst./ QCD");
    hQCD_1->GetXaxis()->SetTitle("x");
    //hQCD_1->GetYaxis()->SetTitleSize(0.1);
    hQCD_1->GetXaxis()->SetTitleSize(0.08);
    
    hQCD_1->Draw("hist");
    QCD_uprms_1->Draw("histsame");
    QCD_downrms_1->Draw("histsame");
    QCD_upfit_1->Draw("histsame");
    QCD_downfit_1->Draw("histsame");
    p2->RedrawAxis();
    canv->ForceUpdate();

    canv->Print(fname);
    canv->Clear();

    
  }//closes loop on the score bin

  canv->Print(fname+"]");

}
