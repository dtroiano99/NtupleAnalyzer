#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));


double Start=38;
double End=198;
int bin=160;
int bin1 = 56-16;
double l=70+8;
double u=126-8;
double s = 200;


static const std::vector<string> functions {"Cebysev"};

int PN_idx = 10; //number possibile fitting coefficinets

static const std::vector<float> PN_scores {0.641, 0.875, 0.957, 0.988, 1};

static const std::vector<int> rebinning {4,8,8,8};


int PN_DIM =  PN_scores.size()-1; //number of PN-MD score bins

float  Fisher(float Chi_a, float Chi_b, int a, int b, int N){

  if(!(b>a)){
    cout<<"Error: second index must be greater than the fisrt one"<<endl;
    float to= -5000000;

    return to;
  }

  float diffchi= Chi_a-Chi_b;
  float num = float(diffchi/(b-a));
  float den = float(Chi_b/(N-b));
  float Fab= num/den;

  float CL = ROOT::Math::fdistribution_cdf_c(Fab,b-a, N-b);

  return CL;

}

string Cebysev(int n){

  stringstream ss;

  string T[PN_idx];
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

void Addiction(TH1F* h1, TH1F* h2, int n){

  int j =0;
  for(int k=1; k< h1->GetNbinsX()+1;k++){

    float e = h1->GetBinContent(k);    

    if(e==0){
      j = j+1;      
      if(j<n+1 || (j>h2->GetNbinsX()-n && j<h2->GetNbinsX()+1)){
	cout<<j<<endl;
	h1->SetBinContent(k,h2->GetBinContent(j));
	h1->SetBinError(k,h2->GetBinError(j));
      }
    }
  }
      

}

void QCD_p8(){


  TH1::AddDirectory(kFALSE);

  TFile *f;
  f = TFile::Open("../Data_histos.root","read");


  //histogams sidebands

  //PN-MD softdropmass histograms
  TH1F *PN_Data[PN_DIM];


  for(int i=0; i<PN_scores.size()-1;i++){
    std::string nomestringa{"PNhist"+std::to_string(i)};
    const char * nome{nomestringa.c_str()};
    PN_Data[i] = (TH1F*)f->Get(nome);
    std::string nomestringa2{"PNhist"+std::to_string(i)+"_s"};
    const char * nome2{nomestringa2.c_str()};
    TH1F* h = (TH1F*)f->Get(nome2);
    Addiction(PN_Data[i], h, 8);
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
    PN_Data[i]->GetXaxis()->SetTitle("m_{SD} / 200 GeV");
    PN_Data[i]->Rebin(rebinning.at(i));
    std::string Y_string{"Events / "+ std::to_string(PN_Data[i]->GetBinWidth(1))};
    const char * Y_title{Y_string.c_str()};
    PN_Data[i]->GetYaxis()->SetTitle(Y_title);
  }

  f->Close();
  delete f;

    //draw raw histograms
  /*
  TString fname1= "./Plots/Data_sb_p8.pdf";
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);
  canv->cd();
  canv->Print(fname1+"[");
  for(int i=0;i<PN_DIM;i++){
    canv->cd();
    canv->Clear();
    gPad->SetLeftMargin(0.15);
    PN_Data[i]->Draw("hist");
    canv->Print(fname1);    
  }
  canv->Print(fname1+"]");
  */

  //Fit procedure

  TLatex tex; tex.SetNDC(); //begin inscription

  vector<vector<float>> Chi2( PN_DIM,vector<float>(PN_idx,-10000));
  vector<vector<float>> fisher( PN_DIM,vector<float>(PN_idx-1,-1));
  vector<vector<float>> Chi2_N( PN_DIM,vector<float>(PN_idx,-10000));
  vector<float> par(PN_idx);
  vector<float> parerr(PN_idx);



  TString fname= "./Plots/Fit_Data20_p8.pdf";
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);
  canv->cd();
  canv->Print(fname+"[");

  ofstream outfile;
  outfile.open ("results_p8.txt", std::ofstream::out | std::ofstream::trunc);

  for(int j=0; j<functions.size();j++){
    cout<<functions.at(j)<<endl;
    std::string filestring1{"rootfiles/Fit_"+functions.at(j)+"20_p8.root"};
    const char * filename1{filestring1.c_str()};
    TFile *myfile;
    myfile = TFile::Open(filename1,"RECREATE");
    //create directories in each root file for each order of the fitting function
    for(int ju=0; ju<PN_idx; ju++){
      std::string dirstring{std::to_string(ju)+"parameters"};
      const char * dirname{dirstring.c_str()};
      myfile->mkdir(dirname);
    }

    for(int i=0;i<PN_DIM;i++){
      int np=(int)((bin-bin1)/rebinning.at(i));
      cout<< np<<endl;
      canv->cd();
      //canv->SetLogy();
      gPad->SetLeftMargin(0.15);

       if(i==0){par[0]=  4500;}
       if(i==1){par[0]=  3200;}
       if(i==2){par[0]=  1100;}
       if(i==3){par[0]=  280;}

       for(int jOR=0; jOR<PN_idx; jOR++){

         std::string dirstring{std::to_string(jOR)+"parameters"};
         const char * dirname{dirstring.c_str()};
         myfile->cd(dirname);
         string function;
         if(j==0){function = Cebysev(jOR);}
         const char * namefit{function.c_str()};

         std::string funstring{"fun"+std::to_string(i)};
	 
	 const char * funname{funstring.c_str()};
         TF1 *fun1=new TF1(funname, namefit,PN_Data[i]->GetXaxis()->GetXmin(), PN_Data[i]->GetXaxis()->GetXmax());
         if(jOR==0){
                fun1->SetParameter(0,par.at(0));
                fun1->SetParLimits(0, 0, 200*par.at(0));
         }

         for(int jj=0; jj<jOR; jj++){
           fun1->SetParameter(jj,par.at(jj));
         }

         TFitResultPtr r = PN_Data[i]->Fit(funname,"S","", PN_Data[i]->GetXaxis()->GetXmin(),PN_Data[i]->GetXaxis()->GetXmax() );
         for(int jjj=0; jjj<jOR+1; jjj++){
           par[jjj]   = r->Parameter(jjj);// retrieve the value for the parameter jjj
           parerr[jjj]   = r->ParError(jjj); // retrieve the error for the parameter jjj
         }

         float chi2   = r->Chi2();
         Chi2_N[i][jOR]=  chi2 /(np-1-jOR);
         Chi2[i][jOR]=  chi2 ;
         if(jOR>0){
           fisher[i][jOR-1]=Fisher(Chi2.at(i).at(jOR-1), Chi2.at(i).at(jOR), jOR, jOR+1,np);
         }

         if(j==0){tex.DrawLatex(0.4,0.4,"#scale[0.8]{Cebysev polynomial}");}
         if(j==1){tex.DrawLatex(0.4,0.4,"#scale[0.8]{Dijet family function}");}
         float start = floor(PN_scores.at(i) * 10000) / 10000;
         std::stringstream ss1;
         ss1 <<jOR+1;
         std::string str1 = ss1.str();
         std::string nomestringa_title{"#scale[0.8]{"+str1+" parameters}"};
         const char * title{nomestringa_title.c_str()};
         tex.DrawLatex(0.4,0.3,title);
         //r->Print();
         //PN_Data[i]->Draw("A");
         canv->Update();
        TPaveText *pt = new TPaveText(1.240187e-08,-2.552966e-10,
                                 1.957123e-08,4.661017e-11,"br");

   pt->AddText(Form("First range %g #pm %g",r->Parameter(0),r->ParError(0)));
   pt->Draw();
         //TPaveText *st = (TPaveText*)r->FindObject("stats");
         //st->SetX1NDC(.15);
         //st->SetX2NDC(.5);
         canv->Update();
         canv->Print(fname);

         r->Write();
         fun1->Write();

         canv->Clear();

       } //close the cicle on the orders


    }//close cicle on the scores

    myfile->Close();
    delete myfile;
    cout<<functions.at(j)<<endl;

    outfile << functions.at(j)<< endl;
    outfile<<""<<endl;
    outfile<<"Chi2 normalized"<<endl;
    outfile<<Chi2.size()<<endl;
    for (int i = 0; i < Chi2_N.size(); i++) {
      for (int j = 0; j < Chi2_N.at(i).size(); j++){
        outfile<< Chi2_N.at(i).at(j) << " ";
      }
      outfile<< endl;
    }
    outfile<<""<<endl;


    float alpha = 0.05;

    outfile<<"Best # of parameters"<<endl;
    outfile<<""<<endl;
    for (int i = 0; i < fisher.size(); i++) {
      int n=0;
        int np=(int)((bin-bin1)/rebinning.at(i));
        cout<< np<<endl;
      for (int j = 1; j < fisher.at(i).size(); j++){
        float fis = Fisher(Chi2.at(i).at(n), Chi2.at(i).at(j), n+1, j+1, np);
        if(fis<alpha){n=j;}
    }
      outfile <<"score region #"<<i<<": "<<n<<" order"<< endl;
      }
    outfile<<""<<endl;
  }//close cicle on the functions

  canv->Print(fname+"]");


}
