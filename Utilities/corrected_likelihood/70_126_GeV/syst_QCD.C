#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

double Start=38;
double End=198;
int bin=160;
int bin1 = 56;
double l=70;
double u=126;
double s = 200;
int max_int =50;

static const std::vector<float> PN_scores {0.641, 0.875, 0.957, 0.988, 1};

static const std::vector<int> rebinning {4,8,8,8};

static const std::vector<int> idx {5,5,4,3};

float Rms(std::vector<float> v, float m){
  float sum=0;
  for(int i= 0;i<v.size();i++){
    sum= sum  + (v.at(i)-m)*(v.at(i)-m);
  }
  float rms=sqrt(sum/v.size());

  return rms;
}

string Cebysev(int n){

  stringstream ss;

  string T[10];
  T[0] = "1";
  T[1] = "x";
  T[2] = "2*x^(2)-1";
  T[3] = "4*x^(3)-3*x";
  T[4] = "8*x^(4)-8*x^(2)+1";
  T[5] = "16*x^(5)-20*x^(3)+5*x";
  T[6] = "32*x^(6)-48*x^(4)+18*x^(2)-1";
  T[7] = "64*x^(7)-112*x^(5)+56*x^(3)-7*x";
  T[8] = "128*x^(8)-256*x^(6)+160*x^(4)-32*x^(2)+1";
  T[9] = "256*x^(9)-576*x^(7)+432*x^(5)-120*x^(3)+9*x";

  for(int i=0; i<n+1; i++){
    string s2;
    if(i==0){s2="[0]";}
    else{s2="+["+ std::to_string(i)+ "]*("+T[i]+")";}
    ss <<s2;
  }
  std::string str2 = ss.str();

  return str2;
}

void toy(TH1F *Data, TF1 *f,TH1F *original ){
  Double_t ymax = f->GetMaximum(Start/s, End/s);
  
  TRandom3 *eventGenerator = new TRandom3();
  //TRandomRanluxpp *eventGenerator = new TRandomRanluxpp();
  //TRandomMixMax256 *eventGenerator = new TRandomMixMax256();
  //TRandomRanlux48 *eventGenerator = new TRandomRanlux48();
  Double_t x1 = 1;
  eventGenerator->SetSeed(0);  
  while(Data->GetEntries()<=original->GetEntries()){     
    Double_t r= eventGenerator->Uniform(x1);//random number between 0 and x1
    Double_t x = (Start +r*(End-Start))/s;
    Double_t y= eventGenerator->Uniform(ymax);//random number between 0 and ymax   
    if((x<l/s || x>u/s) && y<f->Eval(x)){
      Data->Fill(x);
    }    
  }
}

TF1 * fit(TH1F *Data, int k){
  int PN_idx = idx.at(k)+1;
  vector<float> par(PN_idx);
 
  if(k==0){par[0]=  4500;}
  if(k==1){par[0]=  3200;}
  if(k==2){par[0]=  1100;}
  if(k==3){par[0]=  280;}

  TF1 *flist[PN_idx];

  for(int jOR=0; jOR<PN_idx; jOR++){
    string function = Cebysev(jOR);
    const char * namefit{function.c_str()};

    std::string funstring{std::to_string(k)+"fun"+std::to_string(jOR)};
    const char * funname{funstring.c_str()};
    TF1 *fun1=new TF1(funname, namefit,Data->GetXaxis()->GetXmin(), Data->GetXaxis()->GetXmax());
    if(jOR==0){
      fun1->SetParameter(0,par.at(0));
      fun1->SetParLimits(0, 0, 200*par.at(0));
    }

    for(int jj=0; jj<jOR; jj++){
      fun1->SetParameter(jj,par.at(jj));
    }
    
    TFitResultPtr r = Data->Fit(funname,"S0","", Data->GetXaxis()->GetXmin(),Data->GetXaxis()->GetXmax() );
    for(int jjj=0; jjj<jOR+1; jjj++){
    par[jjj]   = r->Parameter(jjj);// retrieve the value for the parameter jjj
    //parerr[jjj]   = r->ParError(jjj); // retrieve the error for the parameter jjj
    }
       
    flist[jOR]=fun1;
  }//saved all fit functions

  par.clear();
  return flist[idx.at(k)];
}

void syst_QCD(){

  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(1001);
  gStyle->SetOptFit(111);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile* file;
  file = TFile::Open("./rootfiles/Fit_Cebysev20.root","read");

  //histogams sidebands

  //PN-MD softdropmass histograms
  TH1F *PN_Data[idx.size()];

  TFile* hfile;
  hfile = TFile::Open("Data_histos.root","read");
  for(int iPN=0; iPN<PN_scores.size()-1;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)};
    const char * nome{nomestringa.c_str()};    
    PN_Data[iPN] = (TH1F*)hfile->Get(nome);    
    PN_Data[iPN]->Rebin(rebinning.at(iPN));
  }
  hfile->Close();
  delete hfile;

  TFile *myfile;
  myfile = TFile::Open("rootfiles/Syst_QCD.root","RECREATE");

  //create hystograms of QCD and  systematics
  TH1F *systFit[idx.size()];
  TH1F *QCD[idx.size()];
  TH1F *systRMS[idx.size()];
  for(int i=0;i<idx.size();i++){
    std::string QCDstring{"QCD"+std::to_string(i)};
    const char * QCDname{QCDstring.c_str()};
    std::string RMSstring{"RMS"+std::to_string(i)};
    const char * RMSname{RMSstring.c_str()};
    std::string Fitstring{"Fit"+std::to_string(i)};
    const char * Fitname{Fitstring.c_str()};
    systRMS[i]= new TH1F(RMSname, "", bin1, l/s, u/s);
    QCD[i]= new TH1F(QCDname, "", bin1, l/s, u/s);
    systFit[i]= new TH1F(Fitname, "", bin1, l/s, u/s);
    systRMS[i]->Rebin(rebinning.at(i));
    QCD[i]->Rebin(rebinning.at(i));
    systFit[i]->Rebin(rebinning.at(i));
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
    systRMS[i]->SetTitle(title);
    systFit[i]->SetTitle(title);
    systRMS[i]->GetYaxis()->SetTitle("Systematic RMS");
    systFit[i]->GetYaxis()->SetTitle("Systematic Fit");
    QCD[i]->GetYaxis()->SetTitle("QCD events");
  }

  TString fname= "./Plots/Random.pdf";//simulated toy and original function
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);
  canv->cd();
  canv->Print(fname+"[");

  //loop on the score region
  for(int j=0;j<idx.size();j++){

    std::vector<std::vector<float>> yields;

    //take the best function fitting data     
    std::string funstring{std::to_string(idx.at(j))+"parameters/fun"+std::to_string(j)};
    const char * funname{funstring.c_str()};      
    TF1 *f;
    f = (TF1*) file->Get(funname);
    std::string rstring{std::to_string(idx.at(j))+"parameters/TFitResult-PNhist"+std::to_string(j)+"-fun"+std::to_string(j)};
    const char * rname{rstring.c_str()};
    TFitResult* r = (TFitResult*) file->Get(rname);
    float xwidth = rebinning.at(j)*((u-l)/bin1)/s;
    //toy simulation
    for(int yo =0; yo<max_int;yo++){
      
      //toy histo
      TH1F *Data = new TH1F("toy", "Simulated Data", bin, Start/s, End/s);
      std::string Y_string{"Events "};
      const char * Y_title{Y_string.c_str()};
      //PN_Data[iPN]->Rebin(rebinning.at(iPN));
      Data->GetXaxis()->SetTitle("x");
      Data->GetYaxis()->SetTitle(Y_title);
      Data->Rebin(rebinning.at(j));
      float start = floor(PN_scores.at(j) * 10000) / 10000;
      std::stringstream ss1;
      ss1 <<start;
      std::string str1 = ss1.str();
      float end = floor(PN_scores.at(j+1) * 10000) / 10000;
      std::stringstream ss2;
      ss2 <<end;
      std::string str2 = ss2.str();
      std::string nomestringa_title{"Simulated Data: "+str1+" < PN-MD_BBvsQCD #leq "+ str2};
      const char * title{nomestringa_title.c_str()};
      Data->SetTitle(title);
	
      toy(Data, f, PN_Data[j]);           
      //Data->Scale((f->Integral(Start/s, l/s)+f->Integral(u/s, End/s))/(xwidth*Data->Integral(1,Data->GetNbinsX())));       
	
      TF1 *g=fit(Data, j );
      g->SetLineColor(3);
      cout<<"number of parameters: "<<g->GetNpar()<<endl;
    
      canv->cd();
      //canv->SetLogy();
      gPad->SetLeftMargin(0.15);
      Data->Draw("histE");
      f->Draw("same");
      //f->Draw("");
      //g->Draw("same"); 
      canv->Print(fname);
      canv->Clear();

      std::vector<float> values;
      for(int k=1; k< QCD[j]->GetNbinsX()+1;k++){
	float xxData = QCD[j]->GetXaxis()->GetBinCenter(k);
	float yyQCD = (g->Integral(xxData-(xwidth/2),xxData+(xwidth/2)))/xwidth;
	values.push_back(yyQCD);
      }
      yields.push_back(values);
      delete g;
      delete Data;
    }//finishes the loop iteration for the toys generation 
  
    auto covMatrix = r->GetCovarianceMatrix();
    std::cout << "Covariance matrix from the fit "<<endl;
    covMatrix.Print();
    for(int k=1; k< QCD[j]->GetNbinsX()+1;k++){
      //cout<<"aaaaaaaaaaaaaaa!"<<endl;
      float xData = QCD[j]->GetXaxis()->GetBinCenter(k);
      float yQCD = (f->Integral(xData-(xwidth/2),xData+(xwidth/2)))/xwidth;
      float erryQCD = (f->IntegralError(xData-(xwidth/2),xData+(xwidth/2),r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()))/xwidth;
      std::vector<float> values;
      for(int yo =0; yo<max_int;yo++){
	  //cout<<yields.at(yo).at(k)<<endl;
	values.push_back(yields.at(yo).at(k-1));
      }
      float RMS = Rms(values,yQCD);      
      systFit[j]->SetBinContent(k,erryQCD);
      systRMS[j]->SetBinContent(k,RMS);
      QCD[j]->SetBinContent(k,yQCD);
      QCD[j]->SetBinError(k,0);
      //cout<<RMS<<endl;
    }
    //cout<<"end"<<endl;
    myfile->cd();
    systFit[j]->Write();
    systRMS[j]->Write();
    QCD[j]->Write();
    //outfile<<"Root mean square (RMS): "<<Rms(yields)<<" ("<<100*Rms(yields)/integral<<"%)"<<endl;
  }//closes loop on the score bin

  canv->Print(fname+"]");
  for(int y=0;y<idx.size();y++){
    float xwidth = rebinning.at(y)*((u-l)/bin1)/s;
    cout<<xwidth<<endl;
  }
  /*
  TString fname1= "./Plots/SystQCD.pdf";
  auto canv1 = new TCanvas("Canvas1", "Canvas1", 1500, 600);
  canv1->cd();
  canv1->Print(fname1+"[");
  for(int j=0;j<idx.size();j++){
     canv1->cd();
     PN_Syst[j]->Draw("hist");
     canv1->Print(fname1);
     canv1->Clear();
  }
  canv->Print(fname1+"]");
  */
}
