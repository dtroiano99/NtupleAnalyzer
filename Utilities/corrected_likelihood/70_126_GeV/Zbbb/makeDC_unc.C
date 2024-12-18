#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

void makeDC_PN(int bin = 0){
  cout<<bin<<endl;
  TFile *outFile;
  std::string rootstring{"/lustre/home/dtroiano1/CMSSW_10_6_12/src/Zqq/NtupleAnalyzer/Utilities/corrected_likelihood/70_126_GeV/Zbbb/rootfiles/Bin"+std::to_string(bin)+"_unc.root"};
  const char * rootname{rootstring.c_str()};
  outFile = TFile::Open(rootname,"RECREATE");

  std::string datacardstring{"Datacard_bin"+std::to_string(bin)+"_unc.txt"};
  const char * datacardname{datacardstring.c_str()};
  ofstream datacard;
  datacard.open (datacardname, std::ofstream::out | std::ofstream::trunc);
  datacard << "imax 1"<< endl;
  datacard << "jmax 3"<< endl;
  datacard << "kmax *"<< endl;
  datacard << "---------------"<< endl;
  datacard << "shapes * * "<<rootname<<" $PROCESS $PROCESS_$SYSTEMATIC"<< endl;
  datacard << "---------------"<< endl;
  datacard << "bin bin1"<< endl;
  
  TFile *f1;
  f1 = TFile::Open("rootfiles/Histos2.root","read");

  std::string datastring{"PNhist"+std::to_string(bin)+"_s"};
  const char * dataname{datastring.c_str()};

  TH1F *Data = (TH1F*)f1->Get(dataname);
  Data->SetName("data_obs");
  float dataI =  Data->Integral(1,Data->GetNbinsX());
  int   Datarate = (int)dataI;

  outFile->cd();
  Data->Write();

  datacard << "observation "<<Datarate<< endl;
  datacard << "------------------------------"<< endl;
  datacard << left<< setw(15) << "bin"<< setw(15) << "bin1"<< setw(15) << "bin1" << setw(15) << "bin1" << setw(15) << "bin1" << endl;
  datacard << left<< setw(15) << "process"<< setw(15) << "Zbb"<< setw(15) << "Zcc"<< setw(15)<< "Zll"<< setw(15) << "QCD"  << endl;
  datacard << left<< setw(15) << "process"<< setw(15) << "-2"<< setw(15) << "-1"<< setw(15) << "0"<< setw(15) << "1"  << endl;
  
  std::string MCstring{"Zbbhist"+std::to_string(bin)+"_s"};
  const char * MCname{MCstring.c_str()};
  TH1F *MC = (TH1F*)f1->Get(MCname);
  MC->SetName("Zbb");
  float MCI =  MC->Integral(1,MC->GetNbinsX());

  std::string MCccstring{"Zcchist"+std::to_string(bin)+"_s"};
  const char * MCccname{MCccstring.c_str()};
  TH1F *MCcc = (TH1F*)f1->Get(MCccname);
  MCcc->SetName("Zcc");
  float MCccI =  MCcc->Integral(1,MCcc->GetNbinsX());

  std::string MCllstring{"Zllhist"+std::to_string(bin)+"_s"};
  const char * MCllname{MCllstring.c_str()};
  TH1F *MCll = (TH1F*)f1->Get(MCllname);
  MCll->SetName("Zll");
  float MCllI =  MCll->Integral(1,MCll->GetNbinsX());

  outFile->cd();
  MC->Write();
  MCcc->Write();
  MCll->Write();

  //cout<<"Zbb Entries: "<<MC->GetEntries()<<" -> "<<MC->GetMean()<<endl;

  f1->Close();
  delete f1;

  TFile *f2;
  f2 = TFile::Open("../rootfiles/QCD.root","read");

  std::string QCDstring{"QCD"+std::to_string(bin)};
  const char * QCDname{QCDstring.c_str()};
  TH1F *QCD = (TH1F*)f2->Get(QCDname);
  QCD->SetName("QCD");
  float QCDI =  QCD->Integral(1,QCD->GetNbinsX());
  //int   QCDrate = (int)QCDI;
  datacard << left<< setw(15) << "rate"<< setw(15) << MCI<< setw(15) << MCccI<< setw(15) << MCllI<< setw(15) << QCDI  << endl;  
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

  cout<<endl;
  cout<<"fit: +"<<100*(QCD_fitup->Integral(1,QCD_fitup->GetNbinsX())-QCDI)/QCDI<<"%  -("<<100*(QCD_fitdown->Integral(1,QCD_fitdown->GetNbinsX())-QCDI)/QCDI<<")%"<<endl;

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

  cout<<"rms: +"<<100*(QCD_rmsup->Integral(1,QCD_rmsup->GetNbinsX())-QCDI)/QCDI<<"%  -"<<100*(QCD_rmsdown->Integral(1,QCD_rmsdown->GetNbinsX())-QCDI)/QCDI<<"%"<<endl;

  datacard << left<< setw(20) << fitstring<< setw(15) <<"shape"<< setw(15)<< "-"<< setw(15)<< "-"<< setw(15)<< "-"<<setw(15) <<"1"<< endl; 
  datacard << left<< setw(20) << rmsstring<< setw(15) <<"shape"<< setw(15)<< "-"<< setw(15)<< "-"<< setw(15)<< "-"<<setw(15) <<"1"<< endl; 
  datacard << left<< setw(20) << "lumi"<< setw(15) <<"lnN"<< setw(15)<< "1.022"<< setw(15) <<"1.022"<< setw(15)<<"1.022"<< setw(15) << "-"  << endl; 

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
  f3 = TFile::Open("rootfiles/Syst_Zqq2.root","read");
  cout<<endl;
  for(int v=0; v<variation.size();v++){
    for(int s=0; s<sources.size();s++){
      std::string Zsyststring{"Zbb_"+variation.at(v)+"_"+sources.at(s)+"_"+std::to_string(bin)};
      const char * Zsystname{Zsyststring.c_str()};
      TH1F *Zbb_syst = (TH1F*)f3->Get(Zsystname);

      std::string Zsyststring1{"Zbb_"+sources.at(s)+"_s"+std::to_string(bin)+variation.at(v)};
      const char * Zsystname1{Zsyststring1.c_str()};      
      Zbb_syst->SetName(Zsystname1);

      cout<<"Zbb"<<sources.at(s)<<"   "<<variation.at(v)<<"   "<<100*(Zbb_syst->Integral(1,Zbb_syst->GetNbinsX())-MCI)/MCI<<"%"<<endl;
      if(v==0){
	datacard << left<<setw(20)<<sources.at(s)+"_s"+std::to_string(bin)<< setw(15) <<"shape"<< setw(15)<< "1"<< setw(15)<< "1"<< setw(15) << "1"<<setw(15) << "-"  << endl;
      }     
      

      std::string Zccsyststring{"Zcc_"+variation.at(v)+"_"+sources.at(s)+"_"+std::to_string(bin)};
      const char * Zccsystname{Zccsyststring.c_str()};
      TH1F *Zcc_syst = (TH1F*)f3->Get(Zccsystname);

      std::string Zccsyststring1{"Zcc_"+sources.at(s)+"_s"+std::to_string(bin)+variation.at(v)};
      const char * Zccsystname1{Zccsyststring1.c_str()};      
      Zcc_syst->SetName(Zccsystname1);

      cout<<"Zcc"<<sources.at(s)<<"   "<<variation.at(v)<<"   "<<100*(Zcc_syst->Integral(1,Zcc_syst->GetNbinsX())-MCccI)/MCccI<<"%"<<endl;      

      std::string Zllsyststring{"Zll_"+variation.at(v)+"_"+sources.at(s)+"_"+std::to_string(bin)};
      const char * Zllsystname{Zllsyststring.c_str()};
      TH1F *Zll_syst = (TH1F*)f3->Get(Zllsystname);

      std::string Zllsyststring1{"Zll_"+sources.at(s)+"_s"+std::to_string(bin)+variation.at(v)};
      const char * Zllsystname1{Zllsyststring1.c_str()};      
      Zll_syst->SetName(Zllsystname1);

      cout<<"Zll"<<sources.at(s)<<"   "<<variation.at(v)<<"   "<<100*(Zll_syst->Integral(1,Zll_syst->GetNbinsX())-MCllI)/MCllI<<"%"<<endl;      

      outFile->cd();
      Zbb_syst->Write();
      Zcc_syst->Write();
      Zll_syst->Write();
      
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

void makeDC_unc(){
	for(int i=0;i<4;i++){
	  makeDC_PN(i);
	}	
}

