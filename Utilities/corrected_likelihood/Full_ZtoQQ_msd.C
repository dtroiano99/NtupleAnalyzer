#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));


double Start=38;
double End=198;
int bin=80;

//static const std::vector<float> PN_scores {0.641, 1};
static const std::vector<float> PN_scores {0.641, 0.875,  0.957, 0.988, 1};

int PN_DIM =  PN_scores.size()-1; //number of PN-MD score bins

void Full_ZtoQQ_msd(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //gStyle->SetOptStat(1000000001);
  //gStyle->SetOptFit(1111);

  TH1::AddDirectory(kFALSE); 
  
  std::string Y_string{"Events / 2 GeV"};
  const char * Y_title{Y_string.c_str()};

  //histogams signal region

  //PN-MD softdropmass histograms
  TH1F *PN_cor_s[PN_DIM];
  TH1F *PN_raw_s[PN_DIM];
/*  
  for(int iPN=0; iPN<PN_DIM;iPN++){
    std::string nomestringa{"PNhist"+std::to_string(iPN)+"_s"};
    const char * nome{nomestringa.c_str()};    
    PN_cor_s[iPN]= new TH1F(nome, "corrected m_{SD}", bin, Start, End);
    PN_cor_s[iPN]->GetXaxis()->SetTitle("m_{SD}");    
    PN_cor_s[iPN]->GetYaxis()->SetTitle(Y_title);
    std::string nomestringa1{"PNhist_raw"+std::to_string(iPN)+"_s"};
    const char * nome1{nomestringa1.c_str()};    
    PN_raw_s[iPN]= new TH1F(nome1, "raw m_{SD}", bin, Start, End);
    PN_raw_s[iPN]->GetXaxis()->SetTitle("m_{SD}");    
    PN_raw_s[iPN]->GetYaxis()->SetTitle(Y_title);
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
    PN_cor_s[iPN]->SetTitle(title);
    PN_cor_s[iPN]->Sumw2(); 
    PN_raw_s[iPN]->SetTitle(title);
    PN_raw_s[iPN]->Sumw2();  
  }

 
  TChain *T_Data = new TChain("passedEvents");  // name of the tree is the argument
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT800toInf_2022_CD.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT800toInf_2022_E.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT800toInf_2022_FG.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT600to800_2022_CD.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT600to800_2022_E.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT600to800_2022_FG.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT400to600_2022_CD.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT400to600_2022_E.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT400to600_2022_FG.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT200to400_2022_CD.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT200to400_2022_E.root");
  T_Data->Add("/lustre/cms/store/user/dtroiano/Commissioning/corrected/MC/ZtoQQ/ZtoQQ_HT200to400_2022_FG.root");

  
  //definition of the variables and branches
  Double_t AK8_sdmJet0;
  Float_t PNMD_XbbVsQCD;//softdropmass and score of the leading AK8 jet
  Double_t lumiWeight;

  ROOT::VecOps::RVec<double> *AK8PuppiJets_rawsoftdropmass=0;
  vector<int> *ptIdx =0;
  
  T_Data->SetBranchAddress("AK8_sdmJet0", &AK8_sdmJet0);
  T_Data->SetBranchAddress("PNMD_XbbVsQCD", &PNMD_XbbVsQCD);
  T_Data->SetBranchAddress("lumiWeight", &lumiWeight);
  T_Data->SetBranchAddress("AK8PuppiJets_rawsoftdropmass", &AK8PuppiJets_rawsoftdropmass);
  T_Data->SetBranchAddress("ptIndex", &ptIdx);
    
  Int_t nentries = (Int_t)T_Data->GetEntries();

  cout<<nentries<<" entries"<<endl;
  for (int i = 0; i < nentries; i++){
    T_Data->GetEntry(i);
    if(i%10000 == 0){cout<<i<<"/"<<nentries<<" entries done"<<endl;}

    Double_t AK8_raw_sdmJet0 = AK8PuppiJets_rawsoftdropmass->at(ptIdx->at(0));
 
    for(int iPN=0; iPN<PN_DIM;iPN++){
      if(AK8_sdmJet0>Start && AK8_sdmJet0<End && PNMD_XbbVsQCD>PN_scores.at(iPN) && PNMD_XbbVsQCD<=PN_scores.at(iPN+1)){PN_cor_s[iPN]->Fill(AK8_sdmJet0, lumiWeight);}
      if(AK8_raw_sdmJet0>Start && AK8_raw_sdmJet0<End && PNMD_XbbVsQCD>PN_scores.at(iPN) && PNMD_XbbVsQCD<=PN_scores.at(iPN+1)){PN_raw_s[iPN]->Fill(AK8_raw_sdmJet0, lumiWeight);}
    } 



  }
  //saving histograms in a Tfile
  cout<<"writing root files"<<endl;
  TFile *myfile;
  myfile = TFile::Open("ZtoQQ_full_msd_histos.root","RECREATE");
  myfile->cd();

  for(int i=0; i<PN_DIM;i++){
    PN_raw_s[i]->Write();
    PN_cor_s[i]->Write();
  }

  cout<<"root files completed"<<endl;
  myfile->Close();
  delete myfile;
*/
  TFile *myfile;
  myfile = TFile::Open("ZtoQQ_full_msd_histos.root","read");

  TString fname= "./Plots/Full_ZtoQQ_sdm_comparison.pdf";
  auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);
  canv->cd();
  canv->Print(fname+"[");

  for(int i=0; i<PN_scores.size()-1;i++){

    std::string nomestringa{"PNhist"+std::to_string(i)+"_s;2"};
    const char * nome{nomestringa.c_str()};
    PN_cor_s[i] = (TH1F*)myfile->Get(nome);
    
    std::string nomestringa1{"PNhist"+std::to_string(i)+"_s;1"};
    const char * nome1{nomestringa1.c_str()};
    PN_raw_s[i]= (TH1F*)myfile->Get(nome1);

    cout<<"Score bin # "<<i<<endl;

    cout<<"raw softdrop mass"<<endl;
    cout<<"Mean:   "<<PN_raw_s[i]->GetMean()<<" GeV"<<endl;
    cout<<"StdDev: "<<PN_raw_s[i]->GetStdDev()<<" GeV"<<endl;

    cout<<endl;

    cout<<"corrected softdrop mass"<<endl;
    cout<<"Mean:   "<<PN_cor_s[i]->GetMean()<<" GeV"<<endl;
    cout<<"StdDev: "<<PN_cor_s[i]->GetStdDev()<<" GeV"<<endl;
    cout<<endl;

    PN_cor_s[i]->SetLineColor(2);
    PN_raw_s[i]->SetLineColor(4);
    PN_cor_s[i]->SetFillColor(2);
    PN_raw_s[i]->SetFillColor(4);
    PN_cor_s[i]->SetFillStyle(3335);
    PN_raw_s[i]->SetFillStyle(3353);


    canv->cd();
    canv->Clear();
    PN_cor_s[i]->SetMaximum(1.15*PN_cor_s[i]->GetMaximum());
    PN_cor_s[i]->SetMinimum(0.8*PN_cor_s[i]->GetMinimum());
    PN_cor_s[i]->Draw("hist");
    PN_raw_s[i]->Draw("histsame");

    auto legend1 = new TLegend(0.15,0.7,0.3,0.85);
    legend1->AddEntry(PN_raw_s[i],"raw m_{SD}","f");
    legend1->AddEntry(PN_cor_s[i],"corrected m_{SD}","f");
    legend1->SetBorderSize(0);
    legend1->Draw();
    TLatex tex1;
    tex1.SetNDC();
    tex1.DrawLatex(0.75,0.92,"#scale[0.9]{34.4 fb^{-1} (13.6 TeV)}");

    canv->Print(fname);


  }
  canv->Print(fname+"]");  
}
