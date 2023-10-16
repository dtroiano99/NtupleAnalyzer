#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));


double Start=50;
double End=150;
int bin=50;
int bin1 = 16;
double l=74;
double u=106;
double s = 250;

static const std::vector<float> PN_scores {0.768733, 0.936303, 0.977057};

int PN_DIM =  PN_scores.size(); //number of PN-MD score bins

static const std::vector<float> DDX_scores {0.022193, 0.0939224, 0.238822};

int DDX_DIM =  DDX_scores.size(); //number of DDX score bins

static const std::vector<string> functions {"pol","Cebysev","DF"/*,"MDF","PPF", "PEF", "UA2F"*/};

static const std::vector<vector<int>> idx 
{{5,5,4},//indexes of polynomials
    {8,8,4}, //indexes of Cebysev polynomials
      {2,1,1}/*, //indexes of DF
	{3,3,3}, //indexes of MDF
	  {5,4,4},//indexes of PPF 
	    {3,3,3},//indexes of PEF
	    {4,4,4}*/};//indexes of UA2F

void drawratio(TCanvas *C, TH1F *MC, TH1F *DD, TLegend *legend,  int u, int j){
  TLatex tex; tex.SetNDC(); //begin inscription
  gPad->SetLeftMargin(0.15);
  TH1F *hRatio = (TH1F*)DD->Clone("ratio");
  hRatio->GetXaxis()->SetLabelSize(0.07);
  hRatio->GetXaxis()->SetTitleSize(0.06);
  hRatio->GetYaxis()->SetLabelSize(0.07);
  hRatio->GetYaxis()->SetTitleSize(0.06);
  for(int i=1; i< DD->GetNbinsX()+1;i++){
    float yMC = MC->GetBinContent(i);
    if(yMC > 0){      
      float yDD = DD->GetBinContent(i);
      float erryDD = DD->GetBinError(i);
      hRatio->SetBinContent(i,yDD/yMC);
      hRatio->SetBinError(i,erryDD/yMC);
    }
  }
  auto p1 = new TPad{"p1", "Dati vs mc", 0.05, 0.3, 0.95, 0.95};
  p1->SetBottomMargin(0); // Tolgo il margine bianco attorno al pad
  //p1->SetTickx(2);
  //p1->SetTicky(2);
  auto p2 = new TPad{"p2", "Rapporto", 0.05, 0, 0.95, 0.3};
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.25);
  p2->SetLeftMargin(0.1);
  //p2->SetTicky(2);
  p1->Draw();
  p2->Draw();

  p1->cd();
  TH1F* frame = (TH1F*)MC->Clone("");
  frame->SetMinimum(0);
  float max;
  if(MC->GetMaximum()>DD->GetMaximum()){max = MC->GetMaximum() + 0.2*MC->GetMaximum();}
  else{max = DD->GetMaximum()+0.2*DD->GetMaximum();}
  frame->SetMaximum(max);
  frame->GetXaxis()->SetLabelSize(4);
  DD->GetXaxis()->SetLabelSize(4);
  frame->Draw("histe");
  DD->Draw("PEsame");
  legend->Draw();
  if(j==0){tex.DrawLatex(0.25,0.8,"#scale[0.8]{Polynomial}");}
  if(j==1){tex.DrawLatex(0.25,0.8,"#scale[1]{Cebysev polynomial}");}
  if(j==2){tex.DrawLatex(0.25,0.8,"#scale[1]{Dijet family function}");}
  tex.DrawLatex(0.75,0.92,"#scale[0.9]{39.2 fb^{-1} (13.6 TeV)}");


  p2->cd();
  TH1F* ratioframe = (TH1F*)hRatio->Clone("");
  ratioframe->SetMinimum(0);
  ratioframe->SetMaximum(2.5);
  ratioframe->GetYaxis()->SetTitle("Data/MC");
  ratioframe->SetTitle("");
  ratioframe->Draw("X0 E1 p ");
  auto l = new TLine(DD->GetXaxis()->GetXmin(),1,DD->GetXaxis()->GetXmax(),1);
  //l->SetLineWidth(8);
  l->SetLineColor(kBlack);
  l->Draw();
}

void QCD_DDX_Data_1(){
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);

  gStyle->SetOptStat(1000000001);
  gStyle->SetOptFit(1111);

  TH1::AddDirectory(kFALSE); 
  //histogams signal region

  //Data

  TFile *f1;
  f1 = TFile::Open("Data_histos_1.root","read");

  //PN-MD softdropmass histograms
  TH1F *PN_Data[PN_DIM];
  
  for(int iPN=0; iPN<PN_DIM;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};
    PN_Data[iPN] = (TH1F*)f1->Get(nome);
  }
  
  //DDX softdropmass histograms
  TH1F *DDX_Data[DDX_DIM];
  
  for(int iPN=0; iPN<DDX_DIM;iPN++){
    std::string nomestringa{"DDXhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};    
    DDX_Data[iPN]= (TH1F*)f1->Get(nome);
  }

  f1->Close();
  delete f1;

  //MC

  TFile *f2;
  f2 = TFile::Open("MC_histos_1.root","read");

  //PN-MD softdropmass histograms
  TH1F *PN_MC[PN_DIM];
  
  for(int iPN=0; iPN<PN_DIM;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};
    PN_MC[iPN] = (TH1F*)f2->Get(nome);
  }
  
  //DDX softdropmass histograms
  TH1F *DDX_MC[DDX_DIM];
  
  for(int iPN=0; iPN<DDX_DIM;iPN++){
    std::string nomestringa{"DDXhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};    
    DDX_MC[iPN]= (TH1F*)f2->Get(nome);
  }

  f2->Close();
  delete f2;

  //Data-driven procedure

  ofstream myfile;
  myfile.open ("Integral_1.txt", std::ofstream::out | std::ofstream::trunc);

  TString fname= "./Plots/comparison_DDX_1.pdf";
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);  
  canv->cd();
  canv->Print(fname+"[");
 
  TH1F *Zqq_dd[DDX_DIM];
  for(int j=0; j<functions.size();j++){
    myfile << functions.at(j)<< endl;
    cout<<functions.at(j)<<endl;

    //obtain the data-driven Zqq
    for(int i=0;i<DDX_DIM;i++){   
      Zqq_dd[i]= new TH1F;
      std::string QCDstring{"Zqq data-driven "+std::to_string(i)};
      const char * QCDname{QCDstring.c_str()};
      Zqq_dd[i] = (TH1F*)DDX_Data[i]->Clone(QCDname);
      //take the best function fitting data 
      std::string filestring{"./rootfiles_1/DDX_"+ functions.at(j)+".root"};
      const char * filename{filestring.c_str()};
      TFile* file;
      file = TFile::Open(filename,"read");
      
      std::string funstring{std::to_string(idx.at(j).at(i))+"parameters/fun"+std::to_string(i)};
      const char * funname{funstring.c_str()};
      TF1 *f;
      f = (TF1*) file->Get(funname);

      std::string rstring{std::to_string(idx.at(j).at(i))+"parameters/TFitResult-DDXhist"+std::to_string(i)+"-fun"+std::to_string(i)};
      const char * rname{rstring.c_str()};
      TFitResult* r = (TFitResult*) file->Get(rname);
      
      float QCD_ =0;
      float e_QCD_ =0;
      float e_MCZqq_ =0;
      float e_DDZqq_ =0;
      
      //fill the data-driven Zqq histograms
      for(int k=1; k< DDX_Data[i]->GetNbinsX()+1;k++){

	float xData = DDX_Data[i]->GetXaxis()->GetBinCenter(k);
	float xwidth = DDX_Data[i]->GetXaxis()->GetBinWidth(k);
	
	float yData = DDX_Data[i]->GetBinContent(k);
	float erryData = sqrt(yData);
	//float yQCD = f->Eval(xData);
	float yQCD = (f->Integral(xData-(xwidth/2),xData+(xwidth/2)))/xwidth;
	QCD_ = QCD_ +yQCD;
	//float erryQCD = Err( f, xData);
	auto covMatrix = r->GetCovarianceMatrix();
	std::cout << "Covariance matrix from the fit ";
	covMatrix.Print();
	float erryQCD = (f->IntegralError(xData-(xwidth/2),xData+(xwidth/2),r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()))/xwidth;	   	
	e_QCD_ = e_QCD_ + erryQCD*erryQCD;	
	float yZqq = yData - yQCD;	
	float erryZqq = sqrt(erryData*erryData+erryQCD*erryQCD);
	e_DDZqq_ = e_DDZqq_ + erryZqq*erryZqq ;
	
	if(yZqq>0){
	  Zqq_dd[i]->SetBinContent(k,yZqq);
	  Zqq_dd[i]->SetBinError(k,erryZqq);      
	}
	else{
	  Zqq_dd[i]->SetBinContent(k,0);
	  Zqq_dd[i]->SetBinError(k,erryZqq);      
	}

	float yMCZqq = DDX_MC[i]->GetBinContent(k);
	e_MCZqq_ = e_MCZqq_ + PN_MC[i]->GetBinError(k)*PN_MC[i]->GetBinError(k);
      }
      cout<<((u-l)/bin1)/s<<endl;
      float xwidth = ((u-l)/bin1)/s;
      float yQCD = (f->Integral(l/s,u/s))/xwidth;
      float erryQCD = (f->IntegralError(l/s,u/s,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()))/xwidth;	   	
      myfile << "score region #"<<i<<": Data = "<<DDX_Data[i]->Integral()<<" || MC Zqq = "<< DDX_MC[i]->Integral()<<"#pm"<<sqrt(e_MCZqq_)<<" || QCD = "<<QCD_<<"#pm"<<sqrt(e_QCD_)<<" || Data-driven Zqq = "<<DDX_Data[i]->Integral()-QCD_<<"#pm"<<sqrt(e_DDZqq_)<<" ||total QCD = "<<yQCD<<"#pm"<<sqrt(erryQCD)<<  endl;
      DDX_MC[i]->SetFillColor(2);
      DDX_MC[i]->SetLineColor(3);
      Zqq_dd[i]->SetMarkerStyle(20);
      //Zqq_dd[i]->Sumw2();
      
      file->Close();
      delete file;
    }//close loop on DDX score region

    //plot histogrmas    

    auto legend = new TLegend(0.65,0.6,0.83,0.85);
    legend->AddEntry(DDX_MC[0],"MC Z#rightarrow b #bar{b}","fe");
    legend->AddEntry(Zqq_dd[0],"DD Z#rightarrow b #bar{b}","ep");
    legend->SetBorderSize(0);

    canv->cd();    
    for(int ij=0; ij<DDX_DIM; ij++){
      if(ij==0){
	std::string nomestringa_title{functions.at(j)+": DDX_BBvsQCD loose WP"};
	const char * title{nomestringa_title.c_str()};
	DDX_MC[ij]->SetTitle(title);
      }
      if(ij==1){
	std::string nomestringa_title{functions.at(j)+": DDX_BBvsQCD medium WP"};
	const char * title{nomestringa_title.c_str()};
	DDX_MC[ij]->SetTitle(title);
      }
      if(ij==2){
	std::string nomestringa_title{functions.at(j)+": DDX_BBvsQCD tight WP"};
	const char * title{nomestringa_title.c_str()};
	DDX_MC[ij]->SetTitle(title);
      }      
      drawratio(canv,DDX_MC[ij], Zqq_dd[ij],legend,   ij,j);
      //Zqq_dd[ij]->SetMinimum(0);     
    
      canv->Print(fname);
      canv->Clear();
      canv->cd();
    }

    
  }//close loop on function

  myfile.close();
  canv->Print(fname+"]");  
}
