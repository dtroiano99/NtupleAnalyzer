#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));

double Start=38;
double End=198;
int bin=160;
int bin1 = 56;
double l=70;
double u=126;
double s = 200;

static const std::vector<float> PN_scores {0.641, 0.875, 0.957, 0.988, 1};

static const std::vector<int> rebinning {4,8,8,8};

int PN_DIM =  PN_scores.size()-1; //number of PN-MD score bins


void syst_Wjets(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //gStyle->SetOptStat(1000000001);
  //gStyle->SetOptFit(1111);

  TH1::AddDirectory(kFALSE); 
  
  std::string Y_string{"Events / "};
  const char * Y_title{Y_string.c_str()};
 
  TChain *T_Data = new TChain("passedEvents");  // name of the tree is the argument
  //T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/Wjets/Wjets_HT200to400_2022_CD.root"); //is empty
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/Wjets/Wjets_HT200to400_2022_E.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/Wjets/Wjets_HT200to400_2022_FG.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/Wjets/Wjets_HT400to600_2022_CD.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/Wjets/Wjets_HT400to600_2022_E.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/Wjets/Wjets_HT400to600_2022_FG.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/Wjets/Wjets_HT600to800_2022_CD.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/Wjets/Wjets_HT600to800_2022_E.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/Wjets/Wjets_HT600to800_2022_FG.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/Wjets/Wjets_HT800toInf_2022_CD.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/Wjets/Wjets_HT800toInf_2022_E.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/Wjets/Wjets_HT800toInf_2022_FG.root");
  
  //definition of the variables and branches
  Double_t AK8_sdmJet0; Float_t PNMD_XbbVsQCD;//softdropmass and score of the leading AK8 jet
  Double_t lumiWeight;

  map<string, double> *AK8PuppiJets_Lpt_softdropmass_up=0;
  map<string, double> *AK8PuppiJets_Lpt_softdropmass_down=0;

  T_Data->SetBranchAddress("AK8_sdmJet0", &AK8_sdmJet0);
  T_Data->SetBranchAddress("PNMD_XbbVsQCD", &PNMD_XbbVsQCD);
  T_Data->SetBranchAddress("lumiWeight", &lumiWeight);
  T_Data->SetBranchAddress("AK8PuppiJets_Lpt_softdropmass_Up", &AK8PuppiJets_Lpt_softdropmass_up);
  T_Data->SetBranchAddress("AK8PuppiJets_Lpt_softdropmass_Down", &AK8PuppiJets_Lpt_softdropmass_down);


  //histogams signal region

  //PN-MD softdropmass histograms
  TH1F *PN_Wjets[PN_DIM];
  TH1F *PN_Wjets_u[PN_DIM][1000];
  TH1F *PN_Wjets_d[PN_DIM][1000];
  int SYS_DIM = 0; // number of systematics

  T_Data->GetEntry(0);
  
  for(int iPN=0; iPN<PN_DIM;iPN++){
    std::string namestringa{"Whist"+std::to_string(iPN)+"_s"};
    const char * name{namestringa.c_str()};
    PN_Wjets[iPN]= new TH1F(name, "x", bin1, l/s, u/s);
    PN_Wjets[iPN]->Sumw2();
    int k = -1;
    float start = floor(PN_scores.at(iPN) * 10000) / 10000;
    std::stringstream ss1;
    ss1 <<start;
    std::string str1 = ss1.str();
    float end = floor(PN_scores.at(iPN+1) * 10000) / 10000;
    std::stringstream ss2;
    ss2 <<end;
    std::string str2 = ss2.str();
    std::string nomestringa_title{str1+" < PN-MD_{BBvsQCD} #leq "+ str2};
    const char * title{nomestringa_title.c_str()}; 
    for (map<string, double>::iterator i = AK8PuppiJets_Lpt_softdropmass_up->begin(); i != AK8PuppiJets_Lpt_softdropmass_up->end(); i++) {
                        k = k+1;
                        SYS_DIM = SYS_DIM +1;
                        std::string nomestringa_u{"Wjets_Up_"+i->first+"_"+std::to_string(iPN)};
                        const char * nome_u{nomestringa_u.c_str()};    
                        PN_Wjets_u[iPN][k]= new TH1F(nome_u, "x", bin1, l/s, u/s);
                        PN_Wjets_u[iPN][k]->GetXaxis()->SetTitle("m_{SD} / 200 GeV");    
                        PN_Wjets_u[iPN][k]->GetYaxis()->SetTitle(Y_title);
                        PN_Wjets_u[iPN][k]->SetTitle(title);
                        PN_Wjets_u[iPN][k]->Sumw2();
                        std::string nomestringa_d{"Wjets_Down_"+i->first+"_"+std::to_string(iPN)};
                        const char * nome_d{nomestringa_d.c_str()};    
                        PN_Wjets_d[iPN][k]= new TH1F(nome_d, "x", bin1, l/s, u/s);
                        PN_Wjets_d[iPN][k]->GetXaxis()->SetTitle("m_{SD} / 200 GeV");    
                        PN_Wjets_d[iPN][k]->GetYaxis()->SetTitle(Y_title);
                        PN_Wjets_d[iPN][k]->SetTitle(title);
                        PN_Wjets_d[iPN][k]->Sumw2();
    }
  }
  SYS_DIM = SYS_DIM/PN_DIM;
  //cout<<SYS_DIM<<endl;
    
  Int_t nentries = (Int_t)T_Data->GetEntries();
  //Int_t nentries = 100;

  cout<<nentries<<" entries"<<endl;
  for (int i = 0; i < nentries; i++){
    T_Data->GetEntry(i);
    if(i%10000 == 0){cout<<i<<"/"<<nentries<<" entries done"<<endl;}
    
    for(int iPN=0; iPN<PN_DIM;iPN++){

      if(AK8_sdmJet0>l && AK8_sdmJet0<u && PNMD_XbbVsQCD>PN_scores.at(iPN) && PNMD_XbbVsQCD<=PN_scores.at(iPN+1)){
          PN_Wjets[iPN]->Fill(AK8_sdmJet0/s,lumiWeight);
        }

      int k = -1;
	
      for (map<string, double>::iterator i = AK8PuppiJets_Lpt_softdropmass_up->begin(); i != AK8PuppiJets_Lpt_softdropmass_up->end(); i++) {
	k = k +1;
	if(i->second>l && i->second<u && PNMD_XbbVsQCD>PN_scores.at(iPN) && PNMD_XbbVsQCD<=PN_scores.at(iPN+1)){PN_Wjets_u[iPN][k]->Fill((i->second)/s, lumiWeight);}
      }
      
      int yk = -1;
      for (map<string, double>::iterator i = AK8PuppiJets_Lpt_softdropmass_down->begin(); i != AK8PuppiJets_Lpt_softdropmass_down->end(); i++) {
	yk = yk +1;
	if(i->second>l && i->second<u && PNMD_XbbVsQCD>PN_scores.at(iPN) && PNMD_XbbVsQCD<=PN_scores.at(iPN+1)){PN_Wjets_d[iPN][yk]->Fill((i->second)/s, lumiWeight);}
      }
    
    }
    
  }
  
  //saving histograms in a Tfile
  cout<<"writing root files"<<endl;
  TFile *myfile;
  myfile = TFile::Open("rootfiles/Syst_Wjets.root","RECREATE");
  myfile->cd();

  for(int i=0; i<PN_DIM;i++){
    PN_Wjets[i]->Rebin(rebinning.at(i));
    PN_Wjets[i]->Write();
    for(int k=0; k< SYS_DIM; k++){ 
      PN_Wjets_u[i][k]->Rebin(rebinning.at(i)); 
      PN_Wjets_d[i][k]->Rebin(rebinning.at(i));
      PN_Wjets_d[i][k]->Write();
      PN_Wjets_u[i][k]->Write();
    }
  }

  cout<<"root files completed"<<endl;
  myfile->Close();
  delete myfile;
    
}
