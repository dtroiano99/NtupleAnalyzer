#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "Math/GSLIntegrator.h" 
#include "Math/SpecFuncMathMore.h"
#include "/lustrehome/dtroiano1/tdrstyle.C"
#include "/lustrehome/dtroiano1/CMS_lumi.C"

static const std::vector<float> PN_scores {0.588, 0.856, 0.954,  0.988,  1};

float Rms(std::vector<float> v){
  float sum1 = 0;
  for(int y =0;y<v.size();y++){
    sum1 = sum1 + v.at(y);
  }
  float avg = sum1 / v.size(); 
  float sum=0;
  for(int ii= 0;ii<v.size();ii++){
    sum= sum  + (v.at(ii)-avg)*(v.at(ii)-avg);
  }
  float rms=sqrt(sum/(v.size()-1));

  return rms;
}



void QCD_2(){

  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(1001);
  gStyle->SetOptFit(111);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  int bin1 = 56;
  float l[]={70, 62, 78, 54, 86, 62, 78, 70, 70};
  float u[]={126,134,118,142,110,126,126,134,118};
  //int l[]={70, 62, 78,  62, 78, 70, 70};
  //int u[]={126,134,118,126,126,134,118};
  double s = 200;
  
  static const std::vector<int> rebinning {8,8,8,8};
  
  static const std::vector< std::vector<int> > idx {
    {4,4,3,3},
      {4,3,3,3},
	{4,4,3,3},
	  {4,3,3,3},
   	    {4,3,3,3},
	      {4,3,3,3},
		{4,4,3,3},
		  {4,3,3,3},
		    {4,3,3,3}
		    };
  /*static const std::vector< std::vector<int> > idx {
    {5,5,4,3},
      {5,4,4,3},
	{5,5,4,4},
	  {5,4,4,3},
	    {5,5,4,3},
	      {5,4,4,3},
		{5,4,4,4}
		};*/
  
  

  TFile *myfile;
  myfile = TFile::Open("QCD.root","RECREATE");

  TString fname= "./Plots/New_QCD.pdf";
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
  canv->SetTopMargin( 0 );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);
  canv->cd();
  canv->Print(fname+"[");

  for(int i=0;i<rebinning.size();i++){
    cout<<"score region #"<<i<<endl;
    //cout<<endl;
    std::string QCDstring1{"QCD"+std::to_string(i)+"_p"};
    const char * QCDname1{QCDstring1.c_str()};
    std::string QCDstring2{"QCD"+std::to_string(i)};
    const char * QCDname2{QCDstring2.c_str()};
    TH1F* frame= new TH1F(QCDname1, "", bin1, l[0]/s, u[0]/s);
    TH1F* QCD= new TH1F(QCDname2, "", bin1, l[0]/s, u[0]/s);
    frame->Rebin(rebinning.at(i));
    QCD->Rebin(rebinning.at(i));
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
    frame->SetTitle("");
    QCD->SetTitle(title);

    //sytematic plots
    std::string QCD_fit_upstring{"QCD_fit_up"+std::to_string(i)};
    const char * QCD_fit_upname{QCD_fit_upstring.c_str()};
    std::string QCD_fit_dwstring{"QCD_fit_dw"+std::to_string(i)};
    const char * QCD_fit_dwname{QCD_fit_dwstring.c_str()};
    std::string QCD_rms_upstring{"QCD_rms_up"+std::to_string(i)};
    const char * QCD_rms_upname{QCD_rms_upstring.c_str()};
    std::string QCD_rms_dwstring{"QCD_rms_dw"+std::to_string(i)};
    const char * QCD_rms_dwname{QCD_rms_dwstring.c_str()};
    TH1F* QCD_fit_up= new TH1F(QCD_fit_upname, "", bin1, l[0]/s, u[0]/s);
    TH1F* QCD_fit_dw= new TH1F(QCD_fit_dwname, "", bin1, l[0]/s, u[0]/s);
    TH1F* QCD_rms_up= new TH1F(QCD_rms_upname, "", bin1, l[0]/s, u[0]/s);
    TH1F* QCD_rms_dw= new TH1F(QCD_rms_dwname, "", bin1, l[0]/s, u[0]/s);
    QCD_fit_up->SetTitle(title);
    QCD_fit_dw->SetTitle(title);
    QCD_rms_up->SetTitle(title);
    QCD_rms_dw->SetTitle(title);
    QCD_fit_up->Rebin(rebinning.at(i));
    QCD_fit_dw->Rebin(rebinning.at(i));
    QCD_rms_up->Rebin(rebinning.at(i));
    QCD_rms_dw->Rebin(rebinning.at(i));

    TH1F* hratio = (TH1F*)frame->Clone("");      

    TF1* f[idx.size()];
    TFitResult* r[idx.size()];
    float xwidth = rebinning.at(i)*((u[0]-l[0])/bin1)/s;
    
    

    for(int z=0; z< (sizeof(l) / sizeof(l[0]));z++){

      TFile* file1;
      int Start = (int) l[z];
      int End = (int) u[z];
      std::string filestring{"./rootfiles/Fit_Cebysev20_"+std::to_string(Start)+"_"+std::to_string(End)+".root"};
      const char * filename{filestring.c_str()};
      file1 = TFile::Open(filename,"read");
      TF1 *f1;
      std::string funstring1{std::to_string(idx.at(z).at(i))+"parameters/fun"+std::to_string(i)};
      const char * funname1{funstring1.c_str()};
      f1 = (TF1*) file1->Get(funname1);
      f[z] = f1;
    

      std::string rstring{std::to_string(idx.at(z).at(i))+"parameters/TFitResult-PNhist"+std::to_string(i)+"-fun"+std::to_string(i)};
      const char * rname{rstring.c_str()};
      TFitResult* r1 = (TFitResult*) file1->Get(rname); 
      float yQCD = (f[z]->Integral(l[0]/s,u[0]/s))/xwidth;
      float ErryQCD = (f1->IntegralError(l[0]/s,u[0]/s,r1->GetParams(), r1->GetCovarianceMatrix().GetMatrixArray()))/xwidth;
      r[z] = r1;
      
      file1->Close();
      delete file1;
      //delete f1;
    }
    //std::string funstring{std::to_string(idx.at(0).at(i))+"parameters/fun"+std::to_string(i)};
    //const char * funname{funstring.c_str()};
    //TF1 *f;
    //f = (TF1*) file->Get(funname);
    

    

    for(int k=1; k< QCD->GetNbinsX()+1;k++){
      
      float xData = QCD->GetXaxis()->GetBinCenter(k);
      

      vector<float> v;
      vector<float> err;
      //cout<<"bin #"<<k<<endl;
      for(int z=0; z< (sizeof(l) / sizeof(l[0]));z++){

	
        float yQCD = (f[z]->Integral(xData-(xwidth/2),xData+(xwidth/2)))/xwidth;
	float ErryQCD = (f[z]->IntegralError(xData-(xwidth/2),xData+(xwidth/2),r[z]->GetParams(), r[z]->GetCovarianceMatrix().GetMatrixArray()))/xwidth;
	
	
	v.push_back(yQCD);
	err.push_back(ErryQCD);
	//cout<<"diff: "<<abs(yQCD-yQCD0)<<endl;
	
      }//end loop on windows
      //cout<<endl;
      float sum = 0;
      float sum1 = 0;
      for(int y =0;y<v.size();y++){
	sum = sum + v.at(y);
	sum1 = sum1 + err.at(y);
      }
      float avg = sum / v.size(); 
      float avg1 = sum1 / v.size();   

      QCD->SetBinContent(k,avg);
      //QCD->SetBinError(k,Rms(v));
      QCD->SetBinError(k,0);

      QCD_fit_up->SetBinContent(k,avg+avg1);
      QCD_fit_dw->SetBinContent(k,avg-avg1);
      QCD_rms_up->SetBinContent(k,avg+Rms(v));
      QCD_rms_dw->SetBinContent(k,avg-Rms(v));
   

      //cout<<"bin #"<<k<<" rms = "<<Rms(v)<<" average "<<avg<<endl;

      hratio->SetBinContent(k,Rms(v)/avg);

    }//end loop on QCD bins

    
    canv->cd();
    gPad->SetLeftMargin(0.15);
    canv->Clear();

    auto p1 = new TPad{"p1", "Dati vs mc", 0.05, 0.25, 0.95, 0.95};//0.05, 0.3, 0.95, 0.95
    p1->SetBottomMargin(0); // Tolgo il margine bianco attorno al pad
    p1->SetTopMargin(0.1);
    //p1->SetTicky(2);
    p1->SetLeftMargin(0.15);
    auto p2 = new TPad{"p2", "Rapporto", 0.05, 0, 0.95, 0.25};//0.05, 0, 0.95, 0.3
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0.4);
    p2->SetLeftMargin(0.15);
    //p2->SetTicky(2);
    p1->Draw();
    p2->Draw();
    
    p1->cd();

    frame->GetYaxis()->SetLabelSize(0.06);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetTitle("Events");
    frame->GetYaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetRangeUser(QCD->GetMinimum()*0.7,QCD->GetMaximum()*1.45);
    frame->Draw("");
    QCD->SetLineColor(41);
    QCD->SetFillColor(41);
    QCD->SetFillStyle(3325);
    QCD->Draw("histsame");

    auto legend1 = new TLegend(0.2,0.58,0.5,0.88);
    legend1->AddEntry(QCD,"new QCD","f");

    for(int z=0; z< (sizeof(l) / sizeof(l[0]));z++){
      f[z]->SetLineColor(1 +z);
      f[z]->Draw("lsame");
      std::string Lstring{"fit function betweeen "+std::to_string(l[z])+ " and "+std::to_string(u[z])+" GeV"};
      const char * Lname{Lstring.c_str()};
      legend1->AddEntry(f[z],Lname,"l");

    }
    legend1->SetBorderSize(0);
    //if(i==0){legend1->Draw();}
    legend1->Draw();
    p2->cd();
    p2->SetGridy();
    gPad->SetGridy();

    hratio->GetXaxis()->SetTitle("m_{SD} / 200 GeV");
    hratio->GetYaxis()->SetTitle("err QCD / QCD");
    hratio->SetMarkerStyle(8);
    hratio->SetMinimum(0);
    hratio->SetMarkerColor(kBlack);
    hratio->GetXaxis()->SetLabelSize(0.16);
    hratio->GetXaxis()->SetTitleSize(0.16);
    hratio->GetYaxis()->SetTitleSize(0.16);
    hratio->GetYaxis()->SetLabelSize(0.14);
    hratio->GetYaxis()->SetTitleOffset(0.38);
    //hratio->GetYaxis()->SetRangeUser(0.9,1.1);
    hratio->GetYaxis()->SetNdivisions(6);

    hratio->Draw("P");

    int iPeriod = 5;
    int iPos = 0;
    CMS_lumi( p1, iPeriod, iPos );

    canv->Update();
    //legend1->Draw();
    canv->Print(fname);
    myfile->cd();
    QCD->Write();
    QCD_fit_up->Write();
    QCD_fit_dw->Write();
    QCD_rms_up->Write();
    QCD_rms_dw->Write();
    
    cout<<QCD->Integral()<<endl;
    cout<<QCD_fit_up->Integral()<<endl;
    cout<<QCD_fit_dw->Integral()<<endl;
    cout<<QCD_rms_up->Integral()<<endl;
    cout<<QCD_rms_dw->Integral()<<endl;

    cout<<endl;

  }//close loop on score region

  canv->Print(fname+"]");
  myfile->Close();
  delete myfile;
  
}
