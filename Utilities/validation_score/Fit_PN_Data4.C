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

static const std::vector<string> functions {"Cebysev","DF"};

int PN_idx = 10; //number possibile fitting coefficinets

static const std::vector<float> PN_scores0 {0.83226, 0.913256, 0.949388, 0.968723, 0.981101, 0.989198, 0.993106, 0.996553, 1};

static const std::vector<float> PN_scores {0.83226, 0.913256, 0.949388, 0.968723,0.981101, 0.993106,  1};

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


string DF(int n){
  stringstream ss;

  if(n==0){
    string r= "[0]";
    return r;
  }

  if(n==1){
    string r= "[0]*(1-x)^([1])";
    return r;
  }


  
  for(int i=2; i<n+1; i++){
    string s2;
    if(i==2){s2="[2]";}
    else{s2="+["+ std::to_string(i)+ "]*log(x)^("+std::to_string(i-2)+ ")";}
    ss <<s2;
  }
  std::string str2 = ss.str();

  string t= "[0]*(1-x)^([1])*x^(-(" + str2 + "))";

  return t;

}


void Fit_PN_Data4(){
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);

  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);

  TH1::AddDirectory(kFALSE); 
  
  TFile *f;
  f = TFile::Open("Data_histos.root","read");
  

  //histogams sidebands

  //PN-MD softdropmass histograms
  TH1F *PN_Data[PN_DIM];  
  
  for(int iPN=0; iPN<PN_scores0.size()-1;iPN++){
    cout<<iPN<<endl;
    std::string nomestringa{"PNhist"+std::to_string(iPN)};
    const char * nome{nomestringa.c_str()};    
    if(iPN<4){
      PN_Data[iPN] = (TH1F*)f->Get(nome);
      float start = floor(PN_scores.at(iPN) * 10000) / 10000;
      std::stringstream ss1;
      ss1 <<start;
      std::string str1 = ss1.str();
      float end = floor(PN_scores.at(iPN+1) * 10000) / 10000;
      std::stringstream ss2;
      ss2 <<end;
      std::string str2 = ss2.str();
      std::string nomestringa_title{str1+" < PN-MD_BBvsQCD #leq "+ str2};
      const char * title{nomestringa_title.c_str()};
      PN_Data[iPN]->SetTitle(title);
    }
    else{
      if(iPN%2 ==0){
	int i = (int)(iPN -((iPN-4)/2));
	cout<<"  "<<i<<endl;
	PN_Data[i] = (TH1F*)f->Get(nome);
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
	int i = (int)(iPN -1 -(iPN-5)/2);
	cout<<"  "<<i<<endl;
	TH1F * h1 = (TH1F*)f->Get(nome);
	PN_Data[i]->Add(h1);
      }
    }
  }
  
  f->Close();
  delete f;
  
  /*
  //draw raw histograms
  
  TString fname1= "./Plots/Hist4.pdf";
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);  
  canv->cd();
  canv->Print(fname1+"[");
  for(int i=0;i<PN_DIM;i++){       
    canv->cd();
    gPad->SetLeftMargin(0.15);
    PN_Data[i]->Draw("hist");
    canv->Print(fname1);
    canv->Clear();
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

  int np=(int)(bin-bin1);
  cout<< np<<endl;  

  TString fname= "./Plots/PN_Fit4.pdf";
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);  
  canv->cd();
  canv->Print(fname+"[");

  ofstream outfile;
  outfile.open ("results.txt", std::ofstream::out | std::ofstream::trunc);

  for(int j=0; j<functions.size();j++){
    cout<<functions.at(j)<<endl;
    std::string filestring1{"rootfiles/PN_"+functions.at(j)+"4.root"};
    const char * filename1{filestring1.c_str()};
    TFile *myfile;
    myfile = TFile::Open(filename1,"RECREATE");

    myfile->cd();
    //create directories in each root file for each order of the fitting function
    for(int ju=0; ju<PN_idx; ju++){
      std::string dirstring{std::to_string(ju)+"parameters"};
      const char * dirname{dirstring.c_str()};
      myfile->mkdir(dirname);
    }

    for(int i=0;i<PN_DIM;i++){   
       canv->cd();
       //canv->SetLogy(); 
       gPad->SetLeftMargin(0.15);

       if(i==0){par[0]=  3000;}
       if(i==1){par[0]=  1700;}
       if(i==2){par[0]=  1100;}
       if(i==3){par[0]=  700;}
       if(i==4){par[0]=  650;}      
       if(i==5){par[0]=  280;}
     
       for(int jOR=0; jOR<PN_idx; jOR++){
	
	 std::string dirstring{std::to_string(jOR)+"parameters"};
	 const char * dirname{dirstring.c_str()};
	 myfile->cd(dirname);
	 string function;	
	 if(j==0){function = Cebysev(jOR);}
	 if(j==1){function = DF(jOR);}	 
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

	 if(j==0){tex.DrawLatex(0.7,0.4,"#scale[0.8]{Cebysev polynomial}");}
	 if(j==1){tex.DrawLatex(0.7,0.4,"#scale[0.8]{Dijet family function}");}
	 float start = floor(PN_scores.at(i) * 10000) / 10000;
	 std::stringstream ss1;
	 ss1 <<jOR+1;
	 std::string str1 = ss1.str();
	 std::string nomestringa_title{"#scale[0.8]{"+str1+" parameters}"};
	 const char * title{nomestringa_title.c_str()};
	 tex.DrawLatex(0.7,0.3,title);
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
