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

int PN_idx = 10; //number possibile fitting coefficinets

static const std::vector<string> functions {"pol","Cebysev","DF","MDF"};

static const std::vector<float> PN_scores {0.83226, 0.913256, 0.949388, 0.968723, 0.981101, 0.989198, 0.993106, 0.996553, 1};

int PN_DIM =  PN_scores.size()-1; //number of PN-MD score bins

static const std::vector<float> DDX_scores {0.0365136, 0.0753757,0.127416,0.195188, 0.286469, 0.414133, 0.581553, 0.770378,1};

int DDX_DIM =  DDX_scores.size()-1; //number of DDX score bins


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

string MDF(int n){
  stringstream ss;

  if(n==0){
    string r= "[0]";
    return r;
  }

  if(n==1){
    string r= "[0]*(1-x)^([1]/3)";
    return r;
  }


  
  for(int i=2; i<n+1; i++){
    string s2;
    if(i==2){s2="[2]";}
    else{s2="+["+ std::to_string(i)+ "]*log(x)^("+std::to_string(i-2)+ ")";}
    ss <<s2;
  }
  std::string str2 = ss.str();

  string t= "[0]*(1-x)^([1]/3)*x^(-(" + str2 + "))";

  return t;

}


string pol(int n){
  stringstream ss;
  
  for(int i=0; i<n+1; i++){
    string s;
    if(i==0){s="[0]";}
    else{s="+["+ std::to_string(i)+ "]*x^("+ std::to_string(i)+ ")";}
    ss <<s;
  }
std::string str2 = ss.str();
  return str2;

}

void Fit_DDX_Data(){
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);

  gStyle->SetOptStat(1000000001);
  gStyle->SetOptFit(1111);

  TH1::AddDirectory(kFALSE); 
  
  TFile *f;
  f = TFile::Open("Data_histos.root","read");
  

  //histogams sidebands

  //PN-MD softdropmass histograms
  TH1F *PN_Data[PN_DIM];
  
  for(int iPN=0; iPN<PN_DIM;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)};
    const char * nome{nomestringa.c_str()};
    PN_Data[iPN] = (TH1F*)f->Get(nome);
  }
  
  //DDX softdropmass histograms
  TH1F *DDX_Data[PN_DIM];
  
  for(int iPN=0; iPN<DDX_DIM;iPN++){
    std::string nomestringa{"DDXhist"+std::to_string(iPN)};
    const char * nome{nomestringa.c_str()};    
    DDX_Data[iPN]= (TH1F*)f->Get(nome);
  }

  f->Close();
  delete f;

  //Fit procedure

  vector<vector<float>> Chi2( DDX_DIM,vector<float>(PN_idx,-10000));
  vector<vector<float>> fisher( DDX_DIM,vector<float>(PN_idx-1,-1));
  vector<vector<float>> Chi2_N( DDX_DIM,vector<float>(PN_idx,-10000));
  vector<float> par(PN_idx);
  vector<float> parerr(PN_idx);

  int np=(int)(bin-bin1);
  cout<< np<<endl;  

  TString fname= "./Plots/DDX_Fit.pdf";
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);  
  canv->cd();
  canv->Print(fname+"[");

  ofstream outfile;
  outfile.open ("results.txt", std::ofstream::out | std::ofstream::trunc);

  for(int j=0; j<functions.size();j++){
    std::string filestring1{"rootfiles/DDX_"+functions.at(j)+".root"};
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

    for(int i=0;i<DDX_DIM;i++){   
       canv->cd();
       canv->SetLogy(); 
       gPad->SetLeftMargin(0.15);

       if(i==0){par[0]=  7000;}
       if(i==1){par[0]=  2600;}
       if(i==2){par[0]=  900;}

       for(int jOR=0; jOR<PN_idx; jOR++){
	
	 std::string dirstring{std::to_string(jOR)+"parameters"};
	 const char * dirname{dirstring.c_str()};
	 myfile->cd(dirname);
	 string function;
	 if(j==0){function = pol(jOR);}
	 if(j==1){function = Cebysev(jOR);}
	 if(j==2){function = DF(jOR);}
	 if(j==3){function = MDF(jOR);}
	 const char * namefit{function.c_str()};

	 std::string funstring{"fun"+std::to_string(i)};
	 const char * funname{funstring.c_str()};
	 TF1 *fun1=new TF1(funname, namefit,DDX_Data[i]->GetXaxis()->GetXmin(), DDX_Data[i]->GetXaxis()->GetXmax());   
	 if(jOR==0){
 		fun1->SetParameter(0,par.at(0));
		fun1->SetParLimits(0, 0, 200*par.at(0));
	 }
       
	 for(int jj=0; jj<jOR; jj++){
	   fun1->SetParameter(jj,par.at(jj));
	 }

	 TFitResultPtr r = DDX_Data[i]->Fit(funname,"S","", DDX_Data[i]->GetXaxis()->GetXmin(),DDX_Data[i]->GetXaxis()->GetXmax() );
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
