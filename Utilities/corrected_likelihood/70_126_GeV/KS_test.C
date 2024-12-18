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

static const std::vector<int> idx {5,5,4,3};

static const std::vector<int> rebinning {4,8,8,8};

void KS_test(){
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);

  //gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1);

  TH1::AddDirectory(kFALSE);

  TFile *f1;
  f1 = TFile::Open("Data_histos.root","read");

  TFile* file;
  file = TFile::Open("./rootfiles/Fit_Cebysev20.root","read");

  TString fname1= "./Plots/Sidebands_comparison.pdf";
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
    Data->Rebin(rebinning.at(i));

    //take the best function fitting data
    std::string funstring{std::to_string(idx.at(i))+"parameters/fun"+std::to_string(i)};
    const char * funname{funstring.c_str()};
    TF1 *f;
    f = (TF1*) file->Get(funname);
    std::string rstring{std::to_string(idx.at(i))+"parameters/TFitResult-PNhist"+std::to_string(i)+"-fun"+std::to_string(i)};
    const char * rname{rstring.c_str()};
    TFitResult* r = (TFitResult*) file->Get(rname);
    float xwidth = rebinning.at(i)*((u-l)/bin1)/s;


    TH1F *QCD = (TH1F*)Data->Clone("");

    for(int k=1; k< Data->GetNbinsX()+1;k++){
	float xData = QCD->GetXaxis()->GetBinCenter(k);
	float yQCD = (f->Integral(xData-(xwidth/2),xData+(xwidth/2)))/xwidth;
	float erryQCD = (f->IntegralError(xData-(xwidth/2),xData+(xwidth/2),r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()))/xwidth;
	if(Data->GetBinContent(k)==0){
	      QCD->SetBinContent(k,0);
	      QCD->SetBinError(k,0);
	}
	else{
	      QCD->SetBinContent(k,yQCD);
              QCD->SetBinError(k,erryQCD);
	}
    }
    
    cout<<"score region # "<<i<<endl;
    Double_t KS = Data->KolmogorovTest(QCD,"");//The returned function value is the probability of test (much less than one means NOT compatible)
    cout<<"KS probability: "<<KS<<endl;


    Data->SetLineColor(1);
    Data->SetMarkerStyle(8);

    QCD->SetFillColor(2);
    QCD->SetLineColor(2);

    QCD->SetTitle("");

    canv1->cd();
    gPad->SetLeftMargin(0.15);
    canv1->Clear();

    float min = QCD->GetMinimum()*0.8;
    QCD->SetMinimum(min+0.2);
    float max = Data->GetMaximum()*1.1;
    QCD->SetMaximum(max);
    QCD->GetYaxis()->SetLabelSize(0.06);
    QCD->GetYaxis()->SetTitleSize(0.06);
    QCD->GetYaxis()->SetTitleOffset(1.3);
    QCD->GetXaxis()->SetTitle("m_{SD} / 200 GeV");
    QCD->GetXaxis()->SetLabelSize(0.05);
    QCD->GetXaxis()->SetTitleSize(0.05);
    canv1->SetTopMargin(0.1);
    std::string Ystring{"Events / "+std::to_string(QCD->GetBinWidth(1))};
    const char * Yname{Ystring.c_str()};
    QCD->GetYaxis()->SetTitle(Yname);
    QCD->Draw("hist");
    Data->Draw("pesame");

    auto legend1 = new TLegend(0.36,0.4,0.53,0.65);
    legend1->AddEntry(Data,"Data","lpe");
    legend1->AddEntry(QCD,"QCD","f");
    legend1->SetBorderSize(0);
    legend1->Draw("");

    int iPeriod = 5;
    int iPos = 0;
    CMS_lumi( canv1, iPeriod, iPos );



    canv1->Print(fname1);
    


  }//close score loop

  canv1->Print(fname1+"]");

}
