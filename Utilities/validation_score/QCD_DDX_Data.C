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

static const std::vector<float> PN_scores {0.83226, 0.913256, 0.949388, 0.968723, 0.981101, 0.989198, 0.993106, 0.996553, 1};

int PN_DIM =  PN_scores.size()-1; //number of PN-MD score bins

static const std::vector<float> DDX_scores {0.0365136, 0.0753757,0.127416,0.195188, 0.286469, 0.414133, 0.581553, 0.770378,1};

int DDX_DIM =  DDX_scores.size()-1; //number of DDX score bins

static const std::vector<string> functions {"pol","Cebysev","DF","MDF"};

static const std::vector<vector<int>> idx 
{{3,3,2},//indexes of polynomials
    {3,3,2}, //indexes of Cebysev polynomials
      {3,3,2}, //indexes of DF
	{3,3,2}, //indexes of MDF
	  {5,4,3},//indexes of PPF 
	    {5,3,2},//indexes of PEF
	      {5,5,3}};//indexes of UA2F

void drawratio(TCanvas *C, TH1F *MC, TH1F *DD, TLegend *legend,  int u){
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

void QCD_PN_Data(){
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);

  gStyle->SetOptStat(1000000001);
  gStyle->SetOptFit(1111);

  TH1::AddDirectory(kFALSE); 
  //histogams signal region

  //Data

  TFile *f1;
  f1 = TFile::Open("Data_histos.root","read");

  //PN-MD softdropmass histograms
  TH1F *PN_Data[PN_DIM];
  
  for(int iPN=0; iPN<PN_DIM;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};
    PN_Data[iPN] = (TH1F*)f1->Get(nome);
  }
  
  //DDX softdropmass histograms
  TH1F *DDX_Data[PN_DIM];
  
  for(int iPN=0; iPN<DDX_DIM;iPN++){
    std::string nomestringa{"DDXhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};    
    DDX_Data[iPN]= (TH1F*)f1->Get(nome);
  }

  f1->Close();
  delete f1;

  //MC

  TFile *f2;
  f2 = TFile::Open("MC_histos.root","read");

  //PN-MD softdropmass histograms
  TH1F *PN_MC[PN_DIM];
  
  for(int iPN=0; iPN<PN_DIM;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};
    PN_MC[iPN] = (TH1F*)f2->Get(nome);
  }
  
  //DDX softdropmass histograms
  TH1F *DDX_MC[PN_DIM];
  
  for(int iPN=0; iPN<DDX_DIM;iPN++){
    std::string nomestringa{"DDXhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};    
    DDX_MC[iPN]= (TH1F*)f2->Get(nome);
  }

  f2->Close();
  delete f2;

  //PN-MD_BBvsQCD score histogram

  const Int_t NBINS = PN_DIM;
  //Double_t edges[] = {0.83226, 0.913256, 0.949388, 0.968723, 0.981101, 0.989198, 0.993106, 0.996553, 1};
  Double_t edges[NBINS+1];
  for(int iPN=0; iPN<PN_DIM+1;iPN++){
    edges[iPN]=PN_scores.at(iPN);    
  }
  TH1* h = new TH1F("hPN-MD","",NBINS,edges);

  //Data-driven procedure

  ofstream myfile;
  myfile.open ("Integral.txt", std::ofstream::out | std::ofstream::trunc);

  TString fname= "./Plots/comparison_PN.pdf";
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);  
  canv->cd();
  canv->Print(fname+"[");
 
  TH1F *Zqq_dd[PN_DIM];
  for(int j=0; j<functions.size();j++){
    myfile << functions.at(j)<< endl;
    cout<<functions.at(j)<<endl;

    //obtain the data-driven Zqq
    for(int i=0;i<PN_DIM;i++){   
      Zqq_dd[i]= new TH1F;
      std::string QCDstring{"Zqq data-driven "+std::to_string(i)};
      const char * QCDname{QCDstring.c_str()};
      Zqq_dd[i] = (TH1F*)PN_Data[i]->Clone(QCDname);
      //take the best function fitting data 
      std::string filestring{"./rootfiles/PN_"+ functions.at(j)+".root"};
      const char * filename{filestring.c_str()};
      TFile* file;
      file = TFile::Open(filename,"read");
      
      std::string funstring{std::to_string(idx.at(j).at(i))+"parameters/fun"+std::to_string(i)};
      const char * funname{funstring.c_str()};
      TF1 *f;
      f = (TF1*) file->Get(funname);

      std::string rstring{std::to_string(idx.at(j).at(i))+"parameters/TFitResult-PNhist"+std::to_string(i)+"-fun"+std::to_string(i)};
      const char * rname{rstring.c_str()};
      TFitResult* r = (TFitResult*) file->Get(rname);
      
      float QCD_ =0;
      float e_QCD_ =0;
      float e_MCZqq_ =0;
      float e_DDZqq_ =0;
      
      //fill the data-driven Zqq histograms
      for(int k=1; k< PN_Data[i]->GetNbinsX()+1;k++){

	float xData = PN_Data[i]->GetXaxis()->GetBinCenter(k);
	float xwidth = PN_Data[i]->GetXaxis()->GetBinWidth(k);
	
	float yData = PN_Data[i]->GetBinContent(k);
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

	float yMCZqq = PN_MC[i]->GetBinContent(k);
	e_MCZqq_ = e_MCZqq_ + yMCZqq;
      }
      float xwidth = 0.008;
      float yQCD = (f->Integral(l/s,u/s))/xwidth;
      float erryQCD = (f->IntegralError(l/s,u/s,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()))/xwidth;	   	
      myfile << "score region #"<<i<<": Data = "<<PN_Data[i]->Integral()<<" || MC Zqq = "<< PN_MC[i]->Integral()<<"#pm"<<sqrt(e_MCZqq_)<<" || QCD = "<<QCD_<<"#pm"<<sqrt(e_QCD_)<<" || Data-driven Zqq = "<<PN_Data[i]->Integral()-QCD_<<"#pm"<<sqrt(e_DDZqq_)<<" ||total QCD = "<<yQCD<<"#pm"<<sqrt(erryQCD)<<  endl;
      PN_MC[i]->SetFillColor(2);
      PN_MC[i]->SetLineColor(28);
      Zqq_dd[i]->SetMarkerStyle(20);
      //Zqq_dd[i]->Sumw2();
      
      file->Close();
      delete file;
    }//close loop on PN score region

    //plot histogrmas    

    auto legend = new TLegend(0.65,0.6,0.83,0.85);
    legend->AddEntry(PN_MC[0],"MC Z#rightarrow q #bar{q}","f");
    legend->AddEntry(Zqq_dd[0],"DD Z#rightarrow q #bar{q}","ep");
    legend->SetBorderSize(0);

    canv->cd();    
    for(int ij=0; ij<PN_DIM; ij++){
      if(ij==0){
	std::string nomestringa_title{functions.at(j)+": PN-MD_BBvsQCD loose WP"};
	const char * title{nomestringa_title.c_str()};
	PN_MC[ij]->SetTitle(title);
      }
      if(ij==1){
	std::string nomestringa_title{functions.at(j)+": PN-MD_BBvsQCD medium WP"};
	const char * title{nomestringa_title.c_str()};
	PN_MC[ij]->SetTitle(title);
      }
      if(ij==2){
	std::string nomestringa_title{functions.at(j)+": PN-MD_BBvsQCD tight WP"};
	const char * title{nomestringa_title.c_str()};
	PN_MC[ij]->SetTitle(title);
      }      
      drawratio(canv, PN_MC[ij], Zqq_dd[ij],legend,   ij);
      //Zqq_dd[ij]->SetMinimum(0);     
    
      canv->Print(fname);
      canv->Clear();
      canv->cd();
    }

    
  }//close loop on function

  myfile.close();
  canv->Print(fname+"]");  
}
