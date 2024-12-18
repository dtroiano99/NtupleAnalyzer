#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));

void score_match(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1::AddDirectory(kFALSE); 

  double Start=0;
  double End=1;
  int bin=1000;

  //histograms of the leading-pt AK8 jets which minimize the distances from the two b-quarks
  TH1F *Lpt_bAK8_PNMD_BBvsQCD = new TH1F("Lpt_bAK8_PNMD_BBvsQCD","Lpt_bAK8_PNMD_BBvsQCD", bin, Start, End);
  //TH1F *Lpt_bAK8_PNMD_BBvsQCD = new TH1F("Lpt_bAK8_PNMD_BBvsQCD","Lpt_bAK8_PNMD_BBvsQCD", 144,0.64, End);
  Lpt_bAK8_PNMD_BBvsQCD->SetTitle("");
  Lpt_bAK8_PNMD_BBvsQCD->GetXaxis()->SetTitle("PN-MD_BBvsQCD");
  std::string Y_string{"Events / "+ std::to_string(((End-Start)/bin))};
  const char * Y_title{Y_string.c_str()};
  Lpt_bAK8_PNMD_BBvsQCD->GetYaxis()->SetTitle(Y_title);
  Lpt_bAK8_PNMD_BBvsQCD->Sumw2();

  TChain *T_Zqq = new TChain("passedEvents");  // name of the tree is the argument
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT800toInf_2023C.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT800toInf_2023D.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT600to800_2023C.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT600to800_2023D.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT400to600_2023C.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT400to600_2023D.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT200to400_2023C.root");
  T_Zqq->Add("/lustre/cms/store/user/dtroiano/Commissioning/2023/MC/ZtoQQ/ZtoQQ_HT200to400_2023D.root");

  //definition of the variables and branches
  Float_t AK8_ptJet0, PNMD_XbbVsQCD;//pt and score of the leading AK8 jet
  Double_t lumiWeight;

  Int_t idx_BestZqq_AK8;

  std::vector<float> *AK8_pt = 0;

  

  T_Zqq->SetBranchAddress("AK8_ptJet0", &AK8_ptJet0);
  T_Zqq->SetBranchAddress("PNMD_XbbVsQCD", &PNMD_XbbVsQCD);
  T_Zqq->SetBranchAddress("lumiWeight", &lumiWeight);

  T_Zqq->SetBranchAddress("idx_BestZqq_AK8", &idx_BestZqq_AK8);
    
  T_Zqq->SetBranchAddress("AK8PuppiJets_pt",&AK8_pt);

  Int_t nentries = (Int_t)T_Zqq->GetEntries();

  for (int i = 0; i < nentries; i++){      
    T_Zqq->GetEntry(i); 
    if(i%10000 == 0){cout<<i<<"/"<<nentries<<" entries done"<<endl;}

    if(idx_BestZqq_AK8>-1){
      Float_t BestZbbAK8_pt, BestZbbPNMD_XbbVsQCD, BestZbbDDX_XbbVsQCD;//pt and scores of the AK8 jet which minimizes the distances from the two b-quarks

      BestZbbAK8_pt = AK8_pt->at(idx_BestZqq_AK8);

      if(AK8_ptJet0==BestZbbAK8_pt){
	Lpt_bAK8_PNMD_BBvsQCD->Fill(PNMD_XbbVsQCD, lumiWeight);
      }
    }
    
  }


  //plot histograms

  //WPs evaluation
  const Int_t nq = 10;
  Double_t yq[nq];  // position where to compute the quantiles in [0,1]
  Double_t xq[nq];  // array to contain the quantiles in [min,max]
  for (Int_t i=0;i<nq;i++) {yq[i] = Float_t(i+1)/nq;}

  Lpt_bAK8_PNMD_BBvsQCD->GetQuantiles(nq,xq,yq);
  cout<<"PN-MD_BBvsQCD"<<endl;
  for(int i=0;i<nq;i++){
    cout<<xq[i]<<endl;
    cout<<yq[i]<<endl;
  }

  TCanvas *Z = new TCanvas("CanvasZ", "Z", 800, 600);
  Z->cd();
  //Z->SetLogy();

  Lpt_bAK8_PNMD_BBvsQCD->SetLineColor(4);
  Lpt_bAK8_PNMD_BBvsQCD->SetFillColor(4);
  Lpt_bAK8_PNMD_BBvsQCD->SetFillStyle(3354);

  Z->cd();
  Lpt_bAK8_PNMD_BBvsQCD->Draw("hist");
  Z->Print("./Plots/Best_Lpt_Zqq_scores.pdf");

  cout<<"writing root files"<<endl;
  TFile *myfile;
  myfile = TFile::Open("Score_ZtoBB_hist.root","RECREATE");
  myfile->cd();

  Lpt_bAK8_PNMD_BBvsQCD->Write();
  

  cout<<"root files completed"<<endl;
  myfile->Close();
  delete myfile;

  
}
