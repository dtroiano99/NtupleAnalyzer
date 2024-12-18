#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

int bin1 = 56;
double l=70;
double u=126;
double s = 200;

static const std::vector<float> PN_scores {0.641, 0.875, 0.957, 0.988, 1};

static const std::vector<int> rebinning {4,8,8,8};

static const std::vector<int> idx {5,5,4,3};
static const std::vector<int> idx_p8 {5,5,4,4};
static const std::vector<int> idx_m8 {5,4,4,3};

void QCD_est(){

  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(1001);
  gStyle->SetOptFit(111);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile* file;
  file = TFile::Open("../rootfiles/Fit_Cebysev20.root","read");
  TFile* fp;
  fp = TFile::Open("./rootfiles/Fit_Cebysev20_p8.root","read");
  TFile* fm;
  fm = TFile::Open("./rootfiles/Fit_Cebysev20_m8.root","read");

  TH1F *QCD[idx.size()];
  TH1F *QCD_p8[idx.size()];
  TH1F *QCD_m8[idx.size()];
  for(int i=0;i<idx.size();i++){
    std::string QCDstring{"QCD"+std::to_string(i)};
    const char * QCDname{QCDstring.c_str()};
    std::string QCDstring1{"QCD"+std::to_string(i)+"_p8"};
    const char * QCDname1{QCDstring1.c_str()};
    std::string QCDstring2{"QCD"+std::to_string(i)+"_m8"};
    const char * QCDname2{QCDstring2.c_str()};
    QCD[i]= new TH1F(QCDname, "", bin1, l/s, u/s);
    QCD_p8[i]= new TH1F(QCDname1, "", bin1, l/s, u/s);
    QCD_m8[i]= new TH1F(QCDname2, "", bin1, l/s, u/s);
    QCD[i]->Rebin(rebinning.at(i));
    QCD_p8[i]->Rebin(rebinning.at(i));
    QCD_m8[i]->Rebin(rebinning.at(i));
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
    QCD[i]->SetTitle(title);
    QCD_p8[i]->SetTitle(title);
    QCD_m8[i]->SetTitle(title);
    QCD[i]->GetYaxis()->SetTitle("QCD events");
  }

  TFile *myfile;
  myfile = TFile::Open("rootfiles/Syst_wind_QCD.root","RECREATE");


  //loop on the score region
  for(int j=0;j<idx.size();j++){

    //take the nominal function fitting data
    std::string funstring{std::to_string(idx.at(j))+"parameters/fun"+std::to_string(j)};
    const char * funname{funstring.c_str()};
    TF1 *f;
    f = (TF1*) file->Get(funname);
    float xwidth = rebinning.at(j)*((u-l)/bin1)/s;

    std::string funstring1{std::to_string(idx_p8.at(j))+"parameters/fun"+std::to_string(j)};
    const char * funname1{funstring1.c_str()};
    TF1 *f1;
    f1 = (TF1*) fp->Get(funname1);

    std::string funstring2{std::to_string(idx_m8.at(j))+"parameters/fun"+std::to_string(j)};
    const char * funname2{funstring2.c_str()};
    TF1 *f2;
    f2 = (TF1*) fm->Get(funname2);

    for(int k=1; k< QCD[j]->GetNbinsX()+1;k++){
      float xData = QCD[j]->GetXaxis()->GetBinCenter(k);
      float yQCD = (f->Integral(xData-(xwidth/2),xData+(xwidth/2)))/xwidth;
      float yQCD1 = (f1->Integral(xData-(xwidth/2),xData+(xwidth/2)))/xwidth;
      float yQCD2 = (f2->Integral(xData-(xwidth/2),xData+(xwidth/2)))/xwidth;

      QCD[j]->SetBinContent(k,yQCD);
      QCD_p8[j]->SetBinContent(k,yQCD1);
      QCD_m8[j]->SetBinContent(k,yQCD2);
      
    }
    
    myfile->cd();
    QCD[j]->Write();
    QCD_p8[j]->Write();
    QCD_m8[j]->Write();
  }//close loop on score region 
}

  
