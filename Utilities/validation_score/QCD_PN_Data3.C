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

static const std::vector<string> functions {"Cebysev"};

static const std::vector<float> PN_scores0 {0.83226, 0.913256, 0.949388, 0.968723, 0.981101, 0.989198, 0.993106, 0.996553, 1};

static const std::vector<float> PN_scores {0.83226, 0.949388, 0.981101, 0.993106,  1};

int PN_DIM =  PN_scores.size()-1; //number of PN-MD score bins

static const std::vector<vector<int>> idx 
{{3,3,2,1}};//indexes of Cebysev polynomials


static const std::vector<vector<float>> syst
{{0,0,0,0}};//systematics of Cebysev polynomials
/*
static const std::vector<vector<float>> syst
{{701.036,331.445,143.977,45.8929}};//systematics of Cebysev polynomials
*/

  
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
  if(j==0){tex.DrawLatex(0.25,0.8,"#scale[1]{Cebysev polynomial}");}
  if(j==1){tex.DrawLatex(0.25,0.8,"#scale[1]{Dijet family function}");}
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

float systematic(TF1 *f, TFitResult* r){
  Int_t N = f->GetNpar();
  float xwidth = ((u-l)/bin1)/s;
  float original =(f->Integral(l/s,u/s))/xwidth;
  float systm = 0;
  for(int i=0;i<N;i++){
    //TF1 * g = f->Copy();
    std::string nomestringa{"g"+std::to_string(i)};
    const char * nome{nomestringa.c_str()};
    TF1 *g = (TF1*)f->Clone(nome);
    float p= r->Parameter(i);
    float pe= r->ParError(i);
    cout<<"p = "<<p<<" pm "<<pe<<endl;
    g->SetParameter(i,p-pe);
    float integral1 = (g->Integral(l/s,u/s))/xwidth;
    float diff1 = abs(integral1-original);
    if(diff1>systm){
      systm = diff1;
    }
    g->SetParameter(i,p+pe);
    float integral2 = (g->Integral(l/s,u/s))/xwidth;
    float diff2 = abs(integral2-original);
    if(diff2>systm){
      systm = diff2;
    }
  }
  return systm;
}

void QCD_PN_Data3(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //gStyle->SetGridStyle(3);
  //gStyle->SetGridWidth(3);

  //gStyle->SetOptStat(1000000001);
  //gStyle->SetOptFit(111);

  TH1::AddDirectory(kFALSE); 
  //histogams signal region

  //Data

  TFile *f1;
  f1 = TFile::Open("Data_histos.root","read");

  //PN-MD softdropmass histograms
  TH1F *PN_Data[PN_DIM];
  
  for(int iPN=0; iPN<PN_scores0.size()-1;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};
    if(iPN%2 ==0){
      int i = (int)iPN/2;
      PN_Data[i] = (TH1F*)f1->Get(nome);
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
      PN_Data[i]->SetTitle(title);
    }
    else{
      int i = (int)(iPN-1)/2;
      TH1F * h1 = (TH1F*)f1->Get(nome);
      PN_Data[i]->Add(h1);
    }
  } 
  
  f1->Close();
  delete f1;

  
  
  
  //MC

  TFile *f2;
  f2 = TFile::Open("MC_histos.root","read");

  //PN-MD softdropmass histograms
  TH1F *PN_MC[PN_DIM];
  
  for(int iPN=0; iPN<PN_scores0.size()-1;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};
    if(iPN%2 ==0){
      int i = (int)iPN/2;
      PN_MC[i] = (TH1F*)f2->Get(nome);
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
      PN_MC[i]->SetTitle(title);
    }
    else{
      int i = (int)(iPN-1)/2;
      TH1F * h1 = (TH1F*)f2->Get(nome);
      PN_MC[i]->Add(h1);
    }
  }
  
  f2->Close();
  delete f2;

  //draw raw histograms
  /*
  TString fname1= "./Plots/Hist.pdf";
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);  
  canv->cd();
  canv->Print(fname1+"[");
  for(int i=0;i<PN_DIM;i++){       
    canv->cd();
    gPad->SetLeftMargin(0.15);
    PN_MC[i]->Draw("hist");
    canv->Print(fname1);
    canv->Clear();
  }
  canv->Print(fname1+"]"); 
  */ 

  //PN-MD_BBvsQCD score histograms

  TH1F *score_MC; TH1F *score_QCD[functions.size()]; TH1F *score_DATA;

  const Int_t NBINS = PN_DIM;
  Double_t edges[NBINS+1];
  for(int iPN=0; iPN<PN_DIM+1;iPN++){
    edges[iPN]=PN_scores.at(iPN);    
  }
  score_MC= new TH1F("PN-MD_MC","",NBINS,edges);
  score_DATA= new TH1F("PN-MD_DATA","",NBINS,edges);
  score_MC->GetXaxis()->SetTitle("PN-MD_BBvsQCD score"); 
  score_DATA->GetXaxis()->SetTitle("PN-MD_BBvsQCD score");
  score_MC->GetYaxis()->SetTitle("Events");  
  score_DATA->GetYaxis()->SetTitle("Events");
  score_MC->SetLineColor(4);
  score_MC->SetFillColor(4);
  score_DATA->SetMarkerStyle(20);
  
  //Data-driven procedure

  ofstream myfile;
  myfile.open ("Integral.txt", std::ofstream::out | std::ofstream::trunc);

  TString fname= "./Plots/comparison_PN2.pdf";
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);  
  canv->cd();
  canv->Print(fname+"[");
 
  TH1F *Zqq_dd[PN_DIM];
  for(int j=0; j<functions.size();j++){
    cout<<functions.at(j)<<endl;
    myfile << functions.at(j)<< endl;   

    //definition of the score histograms
    const char * SCname{functions.at(j).c_str()};
    score_QCD[j]= new TH1F("PN-MD_QCD",SCname,NBINS,edges);
    score_QCD[j]->GetXaxis()->SetTitle("PN-MD_BBvsQCD score");
    score_QCD[j]->GetYaxis()->SetTitle("Events");
    score_QCD[j]->SetLineColor(2);
    score_QCD[j]->SetFillColor(2);

    //obtain the data-driven Zqq
    for(int i=0;i<PN_DIM;i++){   
      Zqq_dd[i]= new TH1F;
      std::string QCDstring{"Zqq data-driven "+std::to_string(i)};
      const char * QCDname{QCDstring.c_str()};
      Zqq_dd[i] = (TH1F*)PN_Data[i]->Clone(QCDname);
      //take the best function fitting data 
      std::string filestring{"./rootfiles/PN_"+ functions.at(j)+"2.root"};
      const char * filename{filestring.c_str()};
      TFile* file;
      file = TFile::Open(filename,"read");
      
      std::string funstring{std::to_string(idx.at(j).at(i))+"parameters/fun"+std::to_string(i)};
      const char * funname{funstring.c_str()};
      TF1 *f;
      f = (TF1*) file->Get(funname);

      std::string rstring{std::to_string(idx.at(j).at(i))+"parameters/TFitResult-PNhist"+std::to_string(i*2)+"-fun"+std::to_string(i)};
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
	if(k==1){
	  std::cout << "Covariance matrix from the fit "<<endl; 
	  covMatrix.Print();
	}
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
	e_MCZqq_ = e_MCZqq_ + PN_MC[i]->GetBinError(k)*PN_MC[i]->GetBinError(k);
      }

      //Fill the score histograms
      Float_t valueMC;
      Double_t valueErrorMC;
      valueMC= (Float_t)PN_MC[i]->IntegralAndError(1, PN_MC[i]->GetNbinsX(), valueErrorMC, "");
      Float_t valueDATA;
      Double_t valueErrorDATA;
      valueDATA= (Float_t)PN_Data[i]->IntegralAndError(1, PN_Data[i]->GetNbinsX(), valueErrorDATA, "");
      if(j==0){
	score_MC->SetBinContent(i+1,valueMC);
	score_MC->SetBinError(i+1,valueErrorMC); 
	score_DATA->SetBinContent(i+1,valueDATA);
	score_DATA->SetBinError(i+1,valueErrorDATA); 
      }
      cout<<((u-l)/bin1)/s<<endl;
      //float erryQCD_syst = syst.at(j).at(i);
      float erryQCD_syst = systematic(f, r);
      float xwidth = ((u-l)/bin1)/s;
      float yQCD = (f->Integral(l/s,u/s))/xwidth;
      float erryQCD_stat = (f->IntegralError(l/s,u/s,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()))/xwidth;      
      float erryQCD = sqrt(erryQCD_stat*erryQCD_stat+erryQCD_syst*erryQCD_syst);
      score_QCD[j]->SetBinContent(i+1,yQCD);
      score_QCD[j]->SetBinError(i+1,erryQCD);
      myfile << "score region #"<<i<<": Data = "<<valueDATA<<"#pm"<<valueErrorDATA<<" || MC Zqq = "<< valueMC<<"#pm"<<valueErrorMC<<" || QCD = "<<QCD_<<"#pm"<<sqrt(e_QCD_)<<" || Data-driven Zqq = "<<PN_Data[i]->Integral()-QCD_<<"#pm"<<sqrt(e_DDZqq_)<<" ||total QCD = "<<yQCD<<"#pm"<<erryQCD<<" ("<<erryQCD_stat<<"+"<<erryQCD_syst<<")"<<  endl;
      PN_MC[i]->SetFillColor(2);
      PN_MC[i]->SetLineColor(3);
      //PN_MC[i]->SetLineWidth(2);
      Zqq_dd[i]->SetMarkerStyle(20);
      //Zqq_dd[i]->Sumw2();
      
      file->Close();
      delete file;
    }//close loop on PN score region

    //plot histogrmas    

    auto legend = new TLegend(0.65,0.6,0.83,0.85);
    legend->AddEntry(PN_MC[0],"MC Z#rightarrow b #bar{b}","fe");
    legend->AddEntry(Zqq_dd[0],"DD Z#rightarrow b #bar{b}","ep");
    legend->SetBorderSize(0);

    canv->cd();    
    for(int ij=0; ij<PN_DIM; ij++){      
      drawratio(canv, PN_MC[ij], Zqq_dd[ij],legend,   ij, j);      
      //Zqq_dd[ij]->SetMinimum(0);     
    
      canv->Print(fname);
      canv->Clear();
      canv->cd();
    }

    
  }//close loop on function

  myfile.close();
  canv->Print(fname+"]");  

  //draw score histograms
  TLatex tex1; tex1.SetNDC();
  TString fname1= "./Plots/PNScore3.pdf";
  auto canv1 = new TCanvas("Canvas1", "Canvas1", 2000, 1500);  
  canv1->cd();
  canv1->Print(fname1+"[");
  for(int i=0;i<functions.size();i++){       
    canv1->cd();
    gPad->SetLeftMargin(0.15);
    
    TH1F *hSM = (TH1F*)score_MC->Clone("");
    TH1F *hShadow = (TH1F*)score_MC->Clone("");
    TH1F *hRatio = (TH1F*)score_MC->Clone("");
    TH1F *hDD = (TH1F*)score_MC->Clone("");
    hRatio->GetXaxis()->SetLabelSize(0.07);
    hRatio->GetXaxis()->SetTitleSize(0.06);
    hRatio->GetYaxis()->SetLabelSize(0.07);
    hRatio->GetYaxis()->SetTitleSize(0.06);
    for(int k=1; k< hSM->GetNbinsX()+1;k++){
      float xSM =  score_MC->GetBinContent(k);
      float errxSM = score_MC->GetBinError(k);
      float xDD = score_DATA->GetBinContent(k)-score_QCD[i]->GetBinContent(k);
      float errxDD =sqrt(score_DATA->GetBinError(k)*score_DATA->GetBinError(k)+score_QCD[i]->GetBinError(k)*score_QCD[i]->GetBinError(k));

      if(xSM>0){
	hSM->SetBinContent(k,xSM);
	hSM->SetBinError(k,errxSM);
	hRatio->SetBinContent(k,xDD/xSM);
	hRatio->SetBinError(k,errxDD/xSM);
	hShadow->SetBinContent(k,1);
	hShadow->SetBinError(k,errxSM/xSM);
	hDD->SetBinContent(k,xDD);
	hDD->SetBinError(k,errxDD);
      }
    }

    
    hDD->SetLineColor(1);
    auto p1 = new TPad{"p1", "Dati vs mc", 0.05, 0.3, 0.95, 0.95};//0.05, 0.3, 0.95, 0.95
    p1->SetBottomMargin(0); // Tolgo il margine bianco attorno al pad
    //p1->SetTickx(2);
    p1->SetTicky(2);
    p1->SetLeftMargin(0.105);
    auto p2 = new TPad{"p2", "Rapporto", 0.05, 0., 0.95, 0.3};//0.05, 0, 0.95, 0.3
    p2->SetTopMargin(0);    
    p2->SetBottomMargin(0.25);
    p2->SetLeftMargin(0.105);
    p2->SetTicky(2);
    p1->Draw();
    p2->Draw();

    p1->cd();
    auto legend1 = new TLegend(0.4,0.58,0.65,0.89);
    legend1->AddEntry(score_MC,"MC Z#rightarrow b #bar{b}","f");  
    legend1->AddEntry(hDD,"DD Z#rightarrow b #bar{b}","ep");
    legend1->AddEntry(hSM,"Standard Model","f");
    legend1->SetBorderSize(0);
    //p1->SetLogy();    
    score_MC->GetYaxis()->SetTitle("Events");  
    score_MC->SetMaximum(2000);
    score_MC->SetMinimum(200);
    score_MC->GetYaxis()->SetLabelSize(0.04);
    score_MC->GetYaxis()->SetTitleSize(0.04);
    score_MC->Draw("hist");
    hDD->Draw("PEsame");
    hSM->SetFillColor(kBlack);
    hSM->SetLineColor(1);
    hSM->SetFillStyle(3004);
    //hSM->SetLineColor(kOrange + 4);
    hSM->SetLineWidth(0);
    hSM->SetMarkerStyle(6);
    hSM->SetMarkerSize(0);
    hSM->Draw("e0 e2 L  same ");
    legend1->Draw();
    tex1.DrawLatex(0.75,0.92,"#scale[0.9]{39.2 fb^{-1} (13.6 TeV)}");

    p2->cd();
    p2->SetGridy();
    gPad->SetGridy();
    //p2->SetFillStyle(4000);
    hShadow->SetFillColor(kBlack);
    hShadow->SetFillStyle(3004);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);
    hShadow->GetYaxis()->SetTitle("Data/SM");
    hRatio->SetTitle("");

    hShadow->GetXaxis()->SetLabelSize(0.08);
    hShadow->GetXaxis()->SetTitleSize(0.08);
    hShadow->GetYaxis()->SetTitleSize(0.14);
    hShadow->GetYaxis()->SetLabelSize(0.08);
    hShadow->GetYaxis()->SetTitleOffset(0.4);
    hShadow->GetYaxis()->SetRangeUser(0.66,1.3);
    hShadow->GetYaxis()->SetNdivisions(9);
    //hShadow->GetYaxis()->SetLabelOffset(999.);
    hShadow->SetMarkerStyle(6);
    hShadow->SetMarkerSize(0);
    
    hShadow->SetLineColor(kOrange + 4);
    hShadow->SetLineWidth(8);  
    // hShadow->GetXaxis()->SetTitleSize(0.1);
    //hShadow->GetYaxis()->SetTitleSize(0.12);
    hShadow->Draw("e0 e2 L");
    hRatio->Draw("PE same");
    p2->RedrawAxis();
    canv1->ForceUpdate();

    auto l = new TLine(score_DATA->GetXaxis()->GetXmin(),1,score_DATA->GetXaxis()->GetXmax(),1);
    l->SetLineColor(kRed);
    l->Draw();
    
    canv1->Print(fname1);
    canv1->Clear();
  }
  canv1->Print(fname1+"]");

  /*
  //draw score histograms
  TLatex tex1; tex1.SetNDC();
  TString fname1= "./Plots/PNScore2.pdf";
  auto canv1 = new TCanvas("Canvas1", "Canvas1", 2000, 1500);  
  canv1->cd();
  canv1->Print(fname1+"[");
  for(int i=0;i<functions.size();i++){       
    canv1->cd();
    gPad->SetLeftMargin(0.15);
    THStack *ts_MC= new THStack;      
    ts_MC->Add(score_MC);
    ts_MC->Add(score_QCD[i]);
    ts_MC->SetHistogram(score_QCD[i]);
    ts_MC->SetTitle(score_QCD[i]->GetTitle());
    TH1F *hSM = (TH1F*)score_MC->Clone("");
    TH1F *hShadow = (TH1F*)score_MC->Clone("");
    TH1F *hRatio = (TH1F*)score_MC->Clone("");
    hRatio->GetXaxis()->SetLabelSize(0.07);
    hRatio->GetXaxis()->SetTitleSize(0.06);
    hRatio->GetYaxis()->SetLabelSize(0.07);
    hRatio->GetYaxis()->SetTitleSize(0.06);
    for(int k=1; k< hSM->GetNbinsX()+1;k++){
      float xSM =  score_MC->GetBinContent(k)+score_QCD[i]->GetBinContent(k);
      float errxSM = sqrt(score_MC->GetBinError(k)*score_MC->GetBinError(k)+score_QCD[i]->GetBinError(k)*score_QCD[i]->GetBinError(k));
      float xData = score_DATA->GetBinContent(k);
      float errxData = score_DATA->GetBinError(k);

      if(xSM>0){
	hSM->SetBinContent(k,xSM);
	hSM->SetBinError(k,errxSM);
	hRatio->SetBinContent(k,xData/xSM);
	hRatio->SetBinError(k,errxData/xSM);
	hShadow->SetBinContent(k,1);
	hShadow->SetBinError(k,errxSM/xSM);
      }
    }

    
    score_DATA->SetLineColor(1);
    auto p1 = new TPad{"p1", "Dati vs mc", 0.05, 0.45, 0.95, 0.95};//0.05, 0.3, 0.95, 0.95
    p1->SetBottomMargin(0); // Tolgo il margine bianco attorno al pad
    //p1->SetTickx(2);
    p1->SetTicky(2);
     p1->SetLeftMargin(0.105);
    auto p2 = new TPad{"p2", "Rapporto", 0.05, 0.33, 0.95, 0.45};//0.05, 0, 0.95, 0.3
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0);
    //p2->SetBottomMargin(0.25);
    p2->SetLeftMargin(0.105);
    p2->SetTicky(2);
    auto p3 = new TPad{"p3", "Z over Data", 0.05, 0, 0.95, 0.33};
    p3->SetTopMargin(0);
    p3->SetBottomMargin(0.25);
    p3->SetLeftMargin(0.105);
    p3->SetTicky(2);
    p1->Draw();
    p2->Draw();
    p3->Draw();

    p1->cd();
    auto legend1 = new TLegend(0.69,0.58,0.8,0.89);
    legend1->AddEntry(score_MC,"Z#rightarrow b #bar{b}","f");
    legend1->AddEntry(score_QCD[0],"QCD","f");
    legend1->AddEntry(score_DATA,"Data","ep");
    legend1->AddEntry(hSM,"Standard Model","f");
    legend1->SetBorderSize(0);
    p1->SetLogy();    
    ts_MC->SetMaximum(300000);
    ts_MC->SetMinimum(800);
    ts_MC->GetYaxis()->SetLabelSize(0.06);
    ts_MC->GetYaxis()->SetTitleSize(0.06);
    ts_MC->Draw("hist");
    score_DATA->Draw("PEsame");
    hSM->SetFillColor(kBlack);
    hSM->SetLineColor(1);
    hSM->SetFillStyle(3004);
    //hSM->SetLineColor(kOrange + 4);
    hSM->SetLineWidth(0);
    hSM->SetMarkerStyle(6);
    hSM->SetMarkerSize(0);
    hSM->Draw("e0 e2 L  same ");
    legend1->Draw();
    tex1.DrawLatex(0.75,0.92,"#scale[0.9]{39.2 fb^{-1} (13.6 TeV)}");

    p2->cd();
    p2->SetGridy();
    gPad->SetGridy();
    //p2->SetFillStyle(4000);
    hShadow->SetFillColor(kBlack);
    hShadow->SetFillStyle(3004);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);
    hShadow->GetYaxis()->SetTitle("Data/SM");
    hRatio->SetTitle("");

    hShadow->GetXaxis()->SetLabelSize(0.08);
    hShadow->GetXaxis()->SetTitleSize(0.08);
    hShadow->GetYaxis()->SetTitleSize(0.14);
    hShadow->GetYaxis()->SetLabelSize(0.08);
    hShadow->GetYaxis()->SetTitleOffset(0.4);
    hShadow->GetYaxis()->SetRangeUser(0.96,1.03);
    hShadow->GetYaxis()->SetNdivisions(9);
    //hShadow->GetYaxis()->SetLabelOffset(999.);
    hShadow->SetMarkerStyle(6);
    hShadow->SetMarkerSize(0);
    
    hShadow->SetLineColor(kOrange + 4);
    hShadow->SetLineWidth(8);  
    // hShadow->GetXaxis()->SetTitleSize(0.1);
    //hShadow->GetYaxis()->SetTitleSize(0.12);
    hShadow->Draw("e0 e2 L");
    hRatio->Draw("X0 E1 p same");
    p2->RedrawAxis();
    canv1->ForceUpdate();

    auto l = new TLine(score_DATA->GetXaxis()->GetXmin(),1,score_DATA->GetXaxis()->GetXmax(),1);
    l->SetLineColor(kRed);
    l->Draw();

    p3->cd();
    TH1F *score_MC1 = (TH1F*)score_MC->Clone("");
    score_MC1->Divide(hSM);
    score_MC1->GetYaxis()->SetTitle("Zbb/SM");
    score_MC1->GetXaxis()->SetTitleSize(0.1);
    score_MC1->GetYaxis()->SetTitleSize(0.08);
    score_MC1->Draw("hist");
    
    canv1->Print(fname1);
    canv1->Clear();
  }
  canv1->Print(fname1+"]");
   */ 
}
