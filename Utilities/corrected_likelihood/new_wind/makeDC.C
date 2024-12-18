#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

void makeDC_PN(int bin = 0){

  TFile *outFile;
  std::string rootstring{"/lustre/home/dtroiano1/CMSSW_10_6_12/src/Zqq/NtupleAnalyzer/Utilities/corrected_likelihood/new_wind/rootfiles/Bin"+std::to_string(bin)+".root"};
  const char * rootname{rootstring.c_str()};
  outFile = TFile::Open(rootname,"RECREATE");

  std::string datacardstring{"Datacard_bin"+std::to_string(bin)+".txt"};
  const char * datacardname{datacardstring.c_str()};
  ofstream datacard;
  datacard.open (datacardname, std::ofstream::out | std::ofstream::trunc);
  datacard << "imax 1"<< endl;
  datacard << "jmax 1"<< endl;
  datacard << "kmax *"<< endl;
  datacard << "---------------"<< endl;
  datacard << "shapes * * "<<rootname<<" $PROCESS $PROCESS_$SYSTEMATIC"<< endl;
  datacard << "---------------"<< endl;
  datacard << "bin bin1"<< endl;
  
  TFile *f1;
  f1 = TFile::Open("rootfiles/Histos.root","read");

  std::string datastring{"PNhist"+std::to_string(2*bin)+"_s;1"};
  const char * dataname{datastring.c_str()};

  TH1F *Data = (TH1F*)f1->Get(dataname);
  Data->SetName("data_obs");
  float dataI =  Data->Integral(1,Data->GetNbinsX());
  int   Datarate = (int)dataI;

  outFile->cd();
  Data->Write();

  datacard << "observation "<<Datarate<< endl;
  datacard << "------------------------------"<< endl;
  datacard << left<< setw(15) << "bin"<< setw(15) << "bin1"<< setw(15) << "bin1"  << endl;
  datacard << left<< setw(15) << "process"<< setw(15) << "Zbb"<< setw(15) << "QCD"  << endl;
  datacard << left<< setw(15) << "process"<< setw(15) << "0"<< setw(15) << "1"  << endl;
  
  std::string MCstring{"PNhist"+std::to_string(2*bin)+"_s;2"};
  const char * MCname{MCstring.c_str()};
  TH1F *MC = (TH1F*)f1->Get(MCname);
  MC->SetName("Zbb");
  float MCI =  MC->Integral(1,MC->GetNbinsX());
  //int   MCrate = (int)MCI;

  outFile->cd();
  MC->Write();

  //cout<<"Zbb Entries: "<<MC->GetEntries()<<" -> "<<MC->GetMean()<<endl;

  f1->Close();
  delete f1;

  TFile *f2;
  f2 = TFile::Open("rootfiles/QCD.root","read");

  std::string QCDstring{"QCD"+std::to_string(bin)};
  const char * QCDname{QCDstring.c_str()};
  TH1F *QCD = (TH1F*)f2->Get(QCDname);
  QCD->SetName("QCD");
  float QCDI =  QCD->Integral(1,QCD->GetNbinsX());
  //int   QCDrate = (int)QCDI;
  datacard << left<< setw(15) << "rate"<< setw(15) << MCI<< setw(15) << QCDI  << endl;  
  datacard << "--------------------------------"<< endl;

  std::string rmsstring{"rms_s"+std::to_string(bin)};
  std::string fitstring{"fit_s"+std::to_string(bin)};

  std::string QCDfitupstring{"QCD_systup_fit"+std::to_string(bin)};
  const char * QCDfitupname{QCDfitupstring.c_str()};
  TH1F *QCD_fitup = (TH1F*)f2->Get(QCDfitupname);
  std::string QCDfitup_string{"QCD_"+fitstring+"Up"};
  const char * QCDfitup_name{QCDfitup_string.c_str()}; 
  QCD_fitup->SetName(QCDfitup_name);

  std::string QCDfitdownstring{"QCD_systdown_fit"+std::to_string(bin)};
  const char * QCDfitdownname{QCDfitdownstring.c_str()};
  TH1F *QCD_fitdown = (TH1F*)f2->Get(QCDfitdownname);
  std::string QCDfitdown_string{"QCD_"+fitstring+"Down"};
  const char * QCDfitdown_name{QCDfitdown_string.c_str()};
  QCD_fitdown->SetName(QCDfitdown_name);

  std::string QCDrmsupstring{"QCD_systup_rms"+std::to_string(bin)};
  const char * QCDrmsupname{QCDrmsupstring.c_str()};
  TH1F *QCD_rmsup = (TH1F*)f2->Get(QCDrmsupname);
  std::string QCDrmsup_string{"QCD_"+rmsstring+"Up"};
  const char * QCDrmsup_name{QCDrmsup_string.c_str()};
  QCD_rmsup->SetName(QCDrmsup_name);

  std::string QCDrmsdownstring{"QCD_systdown_rms"+std::to_string(bin)};
  const char * QCDrmsdownname{QCDrmsdownstring.c_str()};
  TH1F *QCD_rmsdown = (TH1F*)f2->Get(QCDrmsdownname);
  std::string QCDrmsdown_string{"QCD_"+rmsstring+"Down"};
  const char * QCDrmsdown_name{QCDrmsdown_string.c_str()};
  QCD_rmsdown->SetName(QCDrmsdown_name);

  datacard << left<< setw(20) << fitstring<< setw(15) <<"shape"<< setw(15)<< "-"<< setw(15) <<"1"<< endl; 
  datacard << left<< setw(20) << rmsstring<< setw(15) <<"shape"<< setw(15)<< "-"<< setw(15) << "1"  << endl; 
  datacard << left<< setw(20) << "lumi"<< setw(15) <<"lnN"<< setw(15)<< "1.022"<< setw(15) << "-"  << endl; 

  outFile->cd();
  QCD->Write();
  QCD_fitup->Write();
  QCD_fitdown->Write();
  QCD_rmsup->Write();
  QCD_rmsdown->Write();
  /*
  cout<<"QCD Entries: "<<QCD->GetEntries()<<" -> "<<QCD->GetMean()<<endl;
  cout<<"fit up: "<<QCD_fitup->GetEntries()<<" -> "<<QCD_fitup->GetMean()<<endl;
  cout<<"fit down: "<<QCD_fitdown->GetEntries()<<" -> "<<QCD_fitdown->GetMean()<<endl;
  cout<<"rms up: "<<QCD_rmsup->GetEntries()<<" -> "<<QCD_rmsup->GetMean()<<endl;
  cout<<"rms down: "<<QCD_rmsdown->GetEntries()<<" -> "<<QCD_rmsdown->GetMean()<<endl;
  */

  f2->Close();
  delete f2;
  
  static const std::vector<string> sources {"AbsoluteFlavMap","AbsoluteMPFBias","AbsoluteSample","AbsoluteScale","AbsoluteStat","FlavorQCD","Fragmentation","PileUpDataMC","PileUpPtBB","PileUpPtEC1","PileUpPtEC2","PileUpPtHF","PileUpPtRef","RelativeBal","RelativeFSR","RelativeJEREC1","RelativeJEREC2","RelativeJERHF","RelativePtBB","RelativePtEC1","RelativePtEC2","RelativePtHF","RelativeSample","RelativeStatEC","RelativeStatFSR","RelativeStatHF","SinglePionECAL","SinglePionHCAL","TimePtEta"};

  static const std::vector<string> variation {"Down","Up"};
  TFile *f3;
  f3 = TFile::Open("rootfiles/Syst_Zbb.root","read");

  for(int v=0; v<variation.size();v++){
    for(int s=0; s<sources.size();s++){
      std::string Zsyststring{"Zbb_"+variation.at(v)+"_"+sources.at(s)+"_"+std::to_string(bin)};
      const char * Zsystname{Zsyststring.c_str()};
      TH1F *Zbb_syst = (TH1F*)f3->Get(Zsystname);

      std::string Zsyststring1{"Zbb_"+sources.at(s)+variation.at(v)};
      const char * Zsystname1{Zsyststring1.c_str()};      
      Zbb_syst->SetName(Zsystname1);

      //cout<<"Entries: "<<Zbb_syst->GetEntries()<<" -> "<<Zbb_syst->GetMean()<<endl;
      if(v==0){
	datacard << left<< setw(20) << sources.at(s)<< setw(15) <<"shape"<< setw(15)<< "1"<< setw(15) << "-"  << endl;
      }
      outFile->cd();
      Zbb_syst->Write();
      
    }
  }

  f3->Close();
  delete f3;
  
  outFile->Close();
  delete outFile;  
  
  datacard << left<< setw(20) << "bin1"<< setw(15) <<"autoMCStats"<< setw(15)<< "0.01"<< setw(15) << "1"  << endl;
  //datacard << "bin1 autoMCStats 0.01 1" << endl;


  datacard.close();

//  delete datacard; 

}

void makeDC(){
	for(int i=0;i<4;i++){
	  makeDC_PN(i);
	}	
}

