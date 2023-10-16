#include <vector>
//#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <cmath>
#include <algorithm>
//float dr1 = deltaR(GENlep_eta->at(r),GENlep_phi->at(r), lep_eta->at(j), lep_phi->at(j));

/*
double Start=50;
double End=150;
int bin=50;
double l=74;
double u=106;
*/

void m_Zbb(){
  //gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1::AddDirectory(kFALSE); 

  double Start=74;
  double End=106;
  int bin=16;

  //histograms of the AK8 jets which pass the loose working point
  TH1F *l_DDX_BBvsQCD = new TH1F("loose WP DDX_BBvsQCD", "loose WP DDX_BBvsQCD", bin, Start, End);
  //l_DDX_BBvsQCD->SetTitle("");
  l_DDX_BBvsQCD->GetXaxis()->SetTitle("m_{SD}");
  std::string Y_string{"Events / "+ std::to_string(((End-Start)/bin))};
  const char * Y_title{Y_string.c_str()};
  l_DDX_BBvsQCD->GetYaxis()->SetTitle(Y_title);
  l_DDX_BBvsQCD->Sumw2();
  TH1F *l_PNMD_BBvsQCD = new TH1F("loose WP PNMD_BBvsQCD","loose WP PNMD_BBvsQCD", bin, Start, End);
  //l_PNMD_BBvsQCD->SetTitle("");
  l_PNMD_BBvsQCD->GetXaxis()->SetTitle("m_{SD}");
  l_PNMD_BBvsQCD->GetYaxis()->SetTitle(Y_title);
  l_PNMD_BBvsQCD->Sumw2();

  //histograms of the AK8 jets which pass the medium working point
  TH1F *m_DDX_BBvsQCD = new TH1F("medium WP DDX_BBvsQCD", "medium WP DDX_BBvsQCD", bin, Start, End);
  //m_DDX_BBvsQCD->SetTitle("");
  m_DDX_BBvsQCD->GetXaxis()->SetTitle("m_{SD}");  
  m_DDX_BBvsQCD->GetYaxis()->SetTitle(Y_title);
  m_DDX_BBvsQCD->Sumw2();
  TH1F *m_PNMD_BBvsQCD = new TH1F("medium WP PNMD_BBvsQCD","medium WP PNMD_BBvsQCD", bin, Start, End);
  //m_PNMD_BBvsQCD->SetTitle("");
  m_PNMD_BBvsQCD->GetXaxis()->SetTitle("m_{SD}");
  m_PNMD_BBvsQCD->GetYaxis()->SetTitle(Y_title);
  m_PNMD_BBvsQCD->Sumw2();

  //histograms of the AK8 jets which pass the tight working point
  TH1F *t_DDX_BBvsQCD = new TH1F("tight WP DDX_BBvsQCD", "tight WP DDX_BBvsQCD", bin, Start, End);
  //t_DDX_BBvsQCD->SetTitle("");
  t_DDX_BBvsQCD->GetXaxis()->SetTitle("m_{SD}");  
  t_DDX_BBvsQCD->GetYaxis()->SetTitle(Y_title);
  t_DDX_BBvsQCD->Sumw2();
  TH1F *t_PNMD_BBvsQCD = new TH1F("tight WP PNMD_BBvsQCD","tight WP PNMD_BBvsQCD", bin, Start, End);
  //t_PNMD_BBvsQCD->SetTitle("");
  t_PNMD_BBvsQCD->GetXaxis()->SetTitle("m_{SD}");
  t_PNMD_BBvsQCD->GetYaxis()->SetTitle(Y_title);
  t_PNMD_BBvsQCD->Sumw2();

  TChain *T_Zqq = new TChain("passedEvents");  // name of the tree is the argument
  T_Zqq->Add("/lustre/home/dtroiano1/CMSSW_10_6_12/src/Zqq/NtupleAnalyzer/scripts/ZJetsToQQ_HT200to400.root");
  T_Zqq->Add("/lustre/home/dtroiano1/CMSSW_10_6_12/src/Zqq/NtupleAnalyzer/scripts/ZJetsToQQ_HT400to600.root");
  T_Zqq->Add("/lustre/home/dtroiano1/CMSSW_10_6_12/src/Zqq/NtupleAnalyzer/scripts/ZJetsToQQ_HT600to800.root");
  T_Zqq->Add("/lustre/home/dtroiano1/CMSSW_10_6_12/src/Zqq/NtupleAnalyzer/scripts/ZJetsToQQ_HT800toInf.root");

  //definition of the variables and branches
  Double_t AK8_sdmJet0;
  Float_t AK8_ptJet0, PNMD_XbbVsQCD, DDX_XbbVsQCD;//pt, softdropmass and scores of the leading AK8 jet
  Double_t lumiWeight;

  Int_t idx_BestZqq_AK8;

  std::vector<float> *AK8_pt = 0;
  std::vector<float> *AK8_msd = 0;
  std::vector<float> *probXbb = 0;
  std::vector<float> *probQCDbb_md = 0;
  std::vector<float> *probQCDcc_md = 0;
  std::vector<float> *probQCDb_md = 0;
  std::vector<float> *probQCDc_md = 0;
  std::vector<float> *probQCDo_md = 0;  
  std::vector<float> *probHbb_ddx = 0;
  

  T_Zqq->SetBranchAddress("AK8_ptJet0", &AK8_ptJet0);
  T_Zqq->SetBranchAddress("AK8_sdmJet0", &AK8_sdmJet0);
  T_Zqq->SetBranchAddress("PNMD_XbbVsQC", &PNMD_XbbVsQCD);
  T_Zqq->SetBranchAddress("DDX_XbbVsQCD", &DDX_XbbVsQCD);
  T_Zqq->SetBranchAddress("lumiWeight", &lumiWeight);

  T_Zqq->SetBranchAddress("idx_BestZqq_AK8", &idx_BestZqq_AK8);
    
  T_Zqq->SetBranchAddress("AK8PuppiJets_pt",&AK8_pt);
  T_Zqq->SetBranchAddress("AK8PuppiJets_softdropmass",&AK8_msd);
  T_Zqq->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probXbb",&probXbb);
  T_Zqq->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb",&probQCDbb_md);
  T_Zqq->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc",&probQCDcc_md);
  T_Zqq->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probQCDb",&probQCDb_md);
  T_Zqq->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probQCDc",&probQCDc_md);
  T_Zqq->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers",&probQCDo_md);
  T_Zqq->SetBranchAddress("jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb",&probHbb_ddx);

  float loosePN = 0.822862;
  float mediumPN = 0.949899;
  float tightPN = 0.981114;

  float looseDDX = 0.0358491;
  float mediumDDX = 0.121852;
  float tightDDX = 0.280009;

  Int_t nentries = (Int_t)T_Zqq->GetEntries();

  for (int i = 0; i < nentries; i++){      
    T_Zqq->GetEntry(i); 

    //select only Z->bb events
    if(idx_BestZqq_AK8 != -2){
      
      if(PNMD_XbbVsQCD > loosePN){l_PNMD_BBvsQCD->Fill(AK8_sdmJet0, lumiWeight);}
      if(PNMD_XbbVsQCD > mediumPN){m_PNMD_BBvsQCD->Fill(AK8_sdmJet0, lumiWeight);}
      if(PNMD_XbbVsQCD > tightPN){t_PNMD_BBvsQCD->Fill(AK8_sdmJet0, lumiWeight);}

      if(DDX_XbbVsQCD > looseDDX){l_DDX_BBvsQCD->Fill(AK8_sdmJet0, lumiWeight);}
      if(DDX_XbbVsQCD > mediumDDX){m_DDX_BBvsQCD->Fill(AK8_sdmJet0, lumiWeight);}
      if(DDX_XbbVsQCD > tightDDX){t_DDX_BBvsQCD->Fill(AK8_sdmJet0, lumiWeight);}

    }
    
  }

  //plot histograms

  //histograms for the PN-MD_BBvsQCD score
  
  TCanvas *PN = new TCanvas("CanvasPN", "PN", 800, 600);
  PN->cd();
  //PN->SetLogy();

  PN->Print("./Plots/PN-MD_threshold.pdf[");

  l_PNMD_BBvsQCD->Draw("histe");
  PN->Print("./Plots/PN-MD_threshold.pdf");
  PN->Clear();
  PN->cd();

  m_PNMD_BBvsQCD->Draw("histe");
  PN->Print("./Plots/PN-MD_threshold.pdf");
  PN->Clear();
  PN->cd();
  
  t_PNMD_BBvsQCD->Draw("histe");
  PN->Print("./Plots/PN-MD_threshold.pdf");
  PN->Clear();
  PN->cd();
 
  PN->Print("./Plots/PN-MD_threshold.pdf]");


  TCanvas *DDX = new TCanvas("CanvasDDX", "DDX", 800, 600);
  DDX->cd();
  //DDX->SetLogy();

  DDX->Print("./Plots/DDX_threshold.pdf[");

  l_DDX_BBvsQCD->Draw("histe");
  DDX->Print("./Plots/DDX_threshold.pdf");
  DDX->Clear();
  DDX->cd();

  m_DDX_BBvsQCD->Draw("histe");
  DDX->Print("./Plots/DDX_threshold.pdf");
  DDX->Clear();
  DDX->cd();
  
  t_DDX_BBvsQCD->Draw("histe");
  DDX->Print("./Plots/DDX_threshold.pdf");
  DDX->Clear();
  DDX->cd();
 
  DDX->Print("./Plots/DDX_threshold.pdf]");
  
  //histograms of the all AK8 jets for PN-MD_BBvsQCD
  /*
  TCanvas *PNMD = new TCanvas("CanvasPNMD", "PNMD", 800, 600);
  PNMD->cd();
  PNMD->SetLogy();

  PNMD_BBvsQCD->SetLineColor(2);
  PNMD_BBvsQCD->SetFillColor(2);
  PNMD_BBvsQCD->SetFillStyle(3354);
  Lpt_PNMD_BBvsQCD->SetLineColor(1);
  Lpt_PNMD_BBvsQCD->SetFillColor(1);
  Lpt_PNMD_BBvsQCD->SetFillStyle(3335);
  bAK8_PNMD_BBvsQCD->SetLineColor(4);
  bAK8_PNMD_BBvsQCD->SetFillColor(4);
  bAK8_PNMD_BBvsQCD->SetFillStyle(3744);
  Lpt_bAK8_PNMD_BBvsQCD->SetLineColor(8);
  Lpt_bAK8_PNMD_BBvsQCD->SetFillColor(8);
  Lpt_bAK8_PNMD_BBvsQCD->SetFillStyle(3636);
  

  float max_PN=PNMD_BBvsQCD->GetMaximum() ;
  if(max_PN < Lpt_PNMD_BBvsQCD->GetMaximum()){max_PN=Lpt_PNMD_BBvsQCD->GetMaximum();}
  if(max_PN < bAK8_PNMD_BBvsQCD->GetMaximum()){max_PN=bAK8_PNMD_BBvsQCD->GetMaximum();}
  if(max_PN < Lpt_bAK8_PNMD_BBvsQCD->GetMaximum()){max_PN=Lpt_bAK8_PNMD_BBvsQCD->GetMaximum();}

  float min_PN=PNMD_BBvsQCD->GetMinimum() ;
  if(min_PN < Lpt_PNMD_BBvsQCD->GetMinimum()){min_PN=Lpt_PNMD_BBvsQCD->GetMinimum();}
  if(min_PN < bAK8_PNMD_BBvsQCD->GetMinimum()){min_PN=bAK8_PNMD_BBvsQCD->GetMinimum();}
  if(min_PN < Lpt_bAK8_PNMD_BBvsQCD->GetMinimum()){min_PN=Lpt_bAK8_PNMD_BBvsQCD->GetMinimum();}

  PNMD_BBvsQCD->SetMaximum(1.1*max_PN);
  PNMD_BBvsQCD->SetMinimum(0.9*min_PN);

  TH1F *frame_PNMD = new TH1F("frame_PNMD","frame_PNMD", bin, Start, End);
  frame_PNMD->SetTitle("");
  frame_PNMD->GetXaxis()->SetTitle("PN-MD_BBvsQCD");
  frame_PNMD->GetYaxis()->SetTitle("Normalized to 1");
  frame_PNMD->SetMaximum(1);
  frame_PNMD->SetMinimum(0.0004);

  frame_PNMD->Draw("hist");
  PNMD_BBvsQCD->DrawNormalized("histsame");
  Lpt_PNMD_BBvsQCD->DrawNormalized("histsame");
  Lpt_bAK8_PNMD_BBvsQCD->DrawNormalized("histsame");
  bAK8_PNMD_BBvsQCD->DrawNormalized("histsame");

  
  auto legend1 = new TLegend(0.3,0.6,0.58,0.8);  
  legend1->AddEntry( PNMD_BBvsQCD,"match","f");
  legend1->AddEntry( Lpt_PNMD_BBvsQCD, "L_pt && no-match","f");
  legend1->AddEntry( bAK8_PNMD_BBvsQCD,"no-L_pt && match","f");
  legend1->AddEntry( Lpt_bAK8_PNMD_BBvsQCD,"L_pt && match","f");
  legend1->SetBorderSize(0);
  legend1->Draw(); 
  
  
  PNMD->Print("./Plots/PNMD_scores_log.png");
  */

  //histograms of the all AK8 jets for DDX_BBvsQCD
  /*
  TCanvas *DDX = new TCanvas("CanvasDDX", "DDX", 800, 600);
  DDX->cd();
  DDX->SetLogy();

  DDX_BBvsQCD->SetLineColor(2);
  DDX_BBvsQCD->SetFillColor(2);
  DDX_BBvsQCD->SetFillStyle(3354);
  Lpt_DDX_BBvsQCD->SetLineColor(1);
  Lpt_DDX_BBvsQCD->SetFillColor(1);
  Lpt_DDX_BBvsQCD->SetFillStyle(3335);
  bAK8_DDX_BBvsQCD->SetLineColor(4);
  bAK8_DDX_BBvsQCD->SetFillColor(4);
  bAK8_DDX_BBvsQCD->SetFillStyle(3744);
  Lpt_bAK8_DDX_BBvsQCD->SetLineColor(8);
  Lpt_bAK8_DDX_BBvsQCD->SetFillColor(8);
  Lpt_bAK8_DDX_BBvsQCD->SetFillStyle(3636);
  
  TH1F *frame_DDX = new TH1F("frame_DDX","frame_DDX", bin, Start, End);
  frame_DDX->SetTitle("");
  frame_DDX->GetXaxis()->SetTitle("DDX_BBvsQCD");
  frame_DDX->GetYaxis()->SetTitle("Normalized to 1");
  frame_DDX->SetMaximum(1);
  frame_DDX->SetMinimum(0.0004);

  frame_DDX->Draw("hist");
  DDX_BBvsQCD->DrawNormalized("histsame");
  Lpt_DDX_BBvsQCD->DrawNormalized("histsame");
  Lpt_bAK8_DDX_BBvsQCD->DrawNormalized("histsame");
  bAK8_DDX_BBvsQCD->DrawNormalized("histsame");

  
  auto legend2 = new TLegend(0.3,0.6,0.58,0.8);  
  legend2->AddEntry( DDX_BBvsQCD,"match","f");
  legend2->AddEntry( Lpt_DDX_BBvsQCD, "L_pt && no-match","f");
  legend2->AddEntry( bAK8_DDX_BBvsQCD,"no-L_pt && match","f");
  legend2->AddEntry( Lpt_bAK8_DDX_BBvsQCD,"L_pt && match","f");
  legend2->SetBorderSize(0);
  legend2->Draw(); 
    
  DDX->Print("./Plots/DDX_scores_log.png");
  */


 
}
