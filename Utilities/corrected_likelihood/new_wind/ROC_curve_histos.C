#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));

double l=70;
double u=126;

double Start=0;
double End=1;
int bin=1000;


void ROC_curve_histos(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);  

  //gStyle->SetOptStat(1000000001);
  //gStyle->SetOptFit(1111);

  TH1::AddDirectory(kFALSE); 
  
  //PN-MD softdropmass histograms

  //histogams sidebands
  TH1F *Data = new TH1F("Data","Data", bin, Start, End);
  TH1F *Zqq = new TH1F("Zqq","Zqq", bin, Start, End);

  TChain *T_Data = new TChain("passedEvents");  // name of the tree is the argument
/*
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022C_corrected.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022D_corrected.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022E_corrected.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022F0_corrected.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022F1_corrected.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022F2_corrected.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/Data/RUN2022G_corrected.root");
*/
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/Data/RUNC_v3.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/Data/RUND_v3.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/Data/RUNE_v3.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/Data/RUNF0_v3.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/Data/RUNF1_v3.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/Data/RUNF2_v3.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/Data/RUNF3_v3.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/Data/RUNG_v3.root");

  //definition of the variables and branches
  Double_t AK8_sdmJet0_d; Float_t PNMD_XbbVsQCD_d;//softdropmass and score of the leading AK8 jet

  

  bool ispreEE;

  ROOT::VecOps::RVec<float> *AK4_eta = 0;
  ROOT::VecOps::RVec<float> *AK4_phi = 0;
  ROOT::VecOps::RVec<bool> *AK4_isloose = 0;
  
  T_Data->SetBranchAddress("AK8_sdmJet0", &AK8_sdmJet0_d);
  T_Data->SetBranchAddress("PNMD_XbbVsQCD", &PNMD_XbbVsQCD_d);
  
  T_Data->SetBranchAddress("ispreEE", &ispreEE);
  
  T_Data->SetBranchAddress("AK4PuppiJets_eta",&AK4_eta);
  T_Data->SetBranchAddress("AK4PuppiJets_phi",&AK4_phi);
  T_Data->SetBranchAddress("AK4PuppiJets_isloose",&AK4_isloose);
  



    
  Int_t nentries = (Int_t)T_Data->GetEntries();

  cout<<nentries<<" entries"<<endl;
  for (int i = 0; i < nentries; i++){      
    T_Data->GetEntry(i);
    if(i%10000 == 0){cout<<i<<"/"<<nentries<<" entries done"<<endl;} 
    
    if(AK8_sdmJet0_d>l && AK8_sdmJet0_d<u ){

	  bool jetvetomap_ = false;
	  
	  TFile *f1;
	  if(ispreEE){
	    f1 = TFile::Open("/lustre/cms/store/user/dtroiano/Commissioning/Utilities/Summer22_23Sep2023_RunCD_v1.root","read");
	  }
	  else{
	    f1 = TFile::Open("/lustre/cms/store/user/dtroiano/Commissioning/Utilities/Summer22EE_23Sep2023_RunEFG_v1.root","read");
	  }
	  
	  TH2D *h_jetvetomap = (TH2D*)f1->Get("jetvetomap");
	  
	  
	  for (UInt_t u4 = 0; u4 < AK4_eta->size(); ++u4) {
	    //cout<<AK4_isloose->at(u4)<<endl;
	    if(AK4_isloose->at(u4) ){
	      double veto = h_jetvetomap->GetBinContent(h_jetvetomap->FindBin(AK4_eta->at(u4), AK4_phi->at(u4)));
	      if(veto>0){
		jetvetomap_ = true;
		break;
	      }
	    }
	  }
	  
	  delete h_jetvetomap;
	  f1->Close();
	  delete f1;
	  
	  if(jetvetomap_ == false){Data->Fill(PNMD_XbbVsQCD_d);}
	  
    }

  }


  TChain *T_Zqq = new TChain("passedEvents");  // name of the tree is the argument
/*
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT800toInf_2022_CD.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT800toInf_2022_E.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT800toInf_2022_FG.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT600to800_2022_CD.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT600to800_2022_E.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT600to800_2022_FG.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT400to600_2022_CD.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT400to600_2022_E.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT400to600_2022_FG.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT200to400_2022_CD.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT200to400_2022_E.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT200to400_2022_FG.root");
*/
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/MC/Zto2Q_200to400_2022_CDv3.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/MC/Zto2Q_200to400_2022_Ev3.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/MC/Zto2Q_200to400_2022_FGv3.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/MC/Zto2Q_400to600_2022_CDv3.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/MC/Zto2Q_400to600_2022_Ev3.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/MC/Zto2Q_400to600_2022_FGv3.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/MC/Zto2Q_600to800_2022_CDv3.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/MC/Zto2Q_600to800_2022_Ev3.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/MC/Zto2Q_600to800_2022_FGv3.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/MC/Zto2Q_800toInf_2022_CDv3.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/MC/Zto2Q_800toInf_2022_Ev3.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/v3/MC/Zto2Q_800toInf_2022_FGv3.root");



  //definition of the variables and branches
  Double_t AK8_sdmJet0_z; Float_t PNMD_XbbVsQCD_z, AK8_ptJet0_z;//softdropmass and score of the leading AK8 jet
  Double_t lumiWeight;

  Int_t idx_BestZqq_AK8;

  std::vector<float> *AK8_pt = 0;
  
  T_Zqq->SetBranchAddress("AK8_sdmJet0", &AK8_sdmJet0_z);
  T_Zqq->SetBranchAddress("PNMD_XbbVsQCD", &PNMD_XbbVsQCD_z);
  T_Zqq->SetBranchAddress("AK8_ptJet0", &AK8_ptJet0_z);
  T_Zqq->SetBranchAddress("lumiWeight", &lumiWeight);
  T_Zqq->SetBranchAddress("idx_BestZqq_AK8", &idx_BestZqq_AK8);
  T_Zqq->SetBranchAddress("AK8PuppiJets_pt",&AK8_pt);
    
  Int_t Nentries = (Int_t)T_Zqq->GetEntries();

  cout<<Nentries<<" entries"<<endl;
  for (int i = 0; i < Nentries; i++){      
    T_Zqq->GetEntry(i);
    if(i%10000 == 0){cout<<i<<"/"<<Nentries<<" entries done"<<endl;} 
    
    if(idx_BestZqq_AK8>-1){
      Float_t BestZbbAK8_pt = AK8_pt->at(idx_BestZqq_AK8);
      if(AK8_sdmJet0_z>l && AK8_sdmJet0_z<u && AK8_ptJet0_z==BestZbbAK8_pt){
	
	Zqq->Fill(PNMD_XbbVsQCD_z, lumiWeight);
      }
	  
    }

  }


  //saving histograms in a Tfile
  cout<<"writing root file"<<endl;

  TFile *myfile;
  myfile = TFile::Open("PN_histosv3Zbb_70_126.root","RECREATE");
  myfile->cd();

  Data->Write();
  Zqq->Write();
  
  cout<<"root file completed"<<endl;
  
  myfile->Close();
  delete myfile;
  
}
