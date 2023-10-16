#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

double Start=50;
double End=150;
int bin=50;
int bin1 = 16;
double l=74;
double u=106;
double s = 250;

static const std::vector<float> PN_scores {0.83226, 0.949388, 0.981101, 0.993106,  1};

static const std::vector<int> idx {3,3,2,1};

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

void toy(TH1F *Data, TF1 *f){
  Double_t ymax = f->GetMaximum(Start/s, End/s);
  float integral = (f->Integral(Start/s, l/s)+f->Integral(u/s, End/s))/(((u-l)/bin1)/s);

  TRandom3 *eventGenerator = new TRandom3();
  //TRandomRanluxpp *eventGenerator = new TRandomRanluxpp();
  //TRandomMixMax256 *eventGenerator = new TRandomMixMax256();
  //TRandomRanlux48 *eventGenerator = new TRandomRanlux48();
  Double_t x1 = 1;
  eventGenerator->SetSeed(0);  
  while(Data->GetEntries()<10000*integral /*&& Data->GetEntries() < 2.26047e+08*/){     
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
 
  if(k==0){par[0]=  4600;}
  if(k==1){par[0]=  1700;}
  if(k==2){par[0]=  650;}
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
  }//saved all 10 fit functions

  par.clear();
  return flist[idx.at(k)];
}

void random_g(){

  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(1001);
  gStyle->SetOptFit(111);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile* file;
  file = TFile::Open("./rootfiles/PN_Cebysev2.root","read");
 
  TString fname= "./Plots/Random_p.pdf";
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);
  canv->cd();
  canv->Print(fname+"[");

  for(int j=0;j<idx.size();j++){
    //toy histo
    TH1F *Data = new TH1F("toy", "Simulated Data", bin, Start/s, End/s);
    std::string Y_string{"Events / "+ std::to_string(((u/s-l/s)/bin1))};
    const char * Y_title{Y_string.c_str()};
    Data->GetXaxis()->SetTitle("x");
    Data->GetYaxis()->SetTitle(Y_title);
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

    //take the best function fitting data     
    std::string funstring{std::to_string(idx.at(j))+"parameters/fun"+std::to_string(j)};
    const char * funname{funstring.c_str()};      
    TF1 *f;
    f = (TF1*) file->Get(funname);
  
    toy(Data, f);   
    float xwidth = ((u-l)/bin1)/s;
    cout<<(f->Integral(Start/s, l/s)+f->Integral(u/s, End/s))/xwidth<<endl;
    //cout<<f->Integral(l/s, u/s)<<endl
    cout<<Data->Integral(1,Data->GetNbinsX())<<endl;
    cout<<"   "<<endl;    
    Data->Scale((f->Integral(Start/s, l/s)+f->Integral(u/s, End/s))/(xwidth*Data->Integral(1,Data->GetNbinsX())));       
    
    TF1 *g=fit(Data, j );
    g->SetLineColor(3);
    cout<<"number of parameters: "<<g->GetNpar()<<endl;
    
    canv->cd();
    //canv->SetLogy();
    gPad->SetLeftMargin(0.15);

    //Data->ResetStats();
    Data->Draw("histE");
    f->Draw("same");
    //f->Draw("");
    //g->Draw("same"); 
    canv->Print(fname);
    canv->Clear();
    cout<<"Difference of the integrals: "<<(f->Integral(l/s, u/s) - g->Integral(l/s, u/s))/xwidth<<endl;
    delete g;
    delete Data;
  }

  canv->Print(fname+"]");

}
