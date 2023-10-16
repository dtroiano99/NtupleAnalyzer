"""
classes to handle ntuplizer output, skimming and adding new variables for AK8 tagger analysis
"""
import ROOT
import sys,os
from Handler import Handler

class Analysis_Ak8Tagger(Handler):


	def ApplyAK8Selection(self):

		#### All the filter functions can be defined here
		ROOT.gInterpreter.Declare("""
		bool TriggerMatch(ROOT::VecOps::RVec<string> &hlt_name, const ROOT::VecOps::RVec<int> &hlt_decision)
		{
			bool HLT_pass = false;
			for(int hk=0; hk<hlt_name.size(); hk++) {
				TString hltName = hlt_name.at(hk);
				if( hltName.Contains("PFHT1050_v") && hlt_decision.at(hk) == 1) {HLT_pass = true;}
				if( !HLT_pass && hltName.Contains("PFJet500") && hlt_decision.at(hk) == 1) {HLT_pass = true;}
				if( !HLT_pass && hltName.Contains("AK8PFJet500") && hlt_decision.at(hk) == 1) {HLT_pass = true;}
				if( !HLT_pass && hltName.Contains("AK8PFJet400_TrimMass30") && hlt_decision.at(hk) == 1) {HLT_pass = true;}
				if( !HLT_pass && hltName.Contains("AK8PFJet420_TrimMass30") && hlt_decision.at(hk) == 1) {HLT_pass = true;}
				if( !HLT_pass && hltName.Contains("AK8PFHT800_TrimMass50") && hlt_decision.at(hk) == 1) {HLT_pass = true;}
			}

			return HLT_pass;
		}

		bool noElectron(ROOT::VecOps::RVec<double> &Ele_pt, ROOT::VecOps::RVec<double> &Ele_eta, ROOT::VecOps::RVec<bool> &Ele_passID)
		{
			bool noEle = true;

			for (UInt_t iEl = 0; iEl < Ele_pt.size(); ++iEl) {
				if(Ele_pt.at(iEl)>20 && abs(Ele_eta.at(iEl))<2.4 && Ele_passID.at(iEl)){
					noEle = false;
					break;
				}
			}
			return noEle;
		}

		bool noMuon(ROOT::VecOps::RVec<double> &Muon_pt, ROOT::VecOps::RVec<double> &Muon_eta, ROOT::VecOps::RVec<double> &Muon_Iso){
			bool noMuon = true;

			for (UInt_t iMu = 0; iMu < Muon_pt.size(); ++iMu) {
				if(Muon_pt.at(iMu)>20 && abs(Muon_eta.at(iMu))<2.4 && Muon_Iso.at(iMu)<0.4){
					noMuon = false;
					break;
				}
			}
			return noMuon;
            
		}

		bool noBtagAK4(ROOT::VecOps::RVec<float> &AK4_pt, ROOT::VecOps::RVec<float> &AK4_eta, ROOT::VecOps::RVec<float> &AK4_phi, ROOT::VecOps::RVec<float> &AK4_DeepCSVTags_b, ROOT::VecOps::RVec<float> &AK4_DeepCSVTags_bb, const float  &AK8_eta, const float  &AK8_phi){
			bool noBtagAK4 = true;

			for (UInt_t i4 = 0; i4 < AK4_pt.size(); ++i4) {
				if(AK4_pt.at(i4)>30 && abs(AK4_eta.at(i4))<2.4 && AK4_DeepCSVTags_b.at(i4)+ AK4_DeepCSVTags_bb.at(i4) >0.4184){
					float dr = sqrt((AK8_eta-AK4_eta.at(i4))*(AK8_eta-AK4_eta.at(i4))+(AK8_phi-AK4_phi.at(i4))*(AK8_phi-AK4_phi.at(i4))); 
                                        if(dr > 0.8){
                                        	noBtagAK4 = false;
                    				break;
                                        } 
				}
			}
			return noBtagAK4;
            
		} 

		""")

		#### All function to get new variables here
		ROOT.gInterpreter.Declare("""
		#include <bits/stdc++.h>
                

		std::vector<int> idxOrder(ROOT::VecOps::RVec<float> &vec)
		{
			std::vector<int> index;			
			std::vector<double> tmpVec;			
			if(vec.size()>0){
				for (unsigned int i=0; i<vec.size(); i++) tmpVec.push_back(vec.at(i));
				std::sort(tmpVec.begin(),tmpVec.end(),greater<float>());
				
				while (tmpVec.size()>0){
				   for (unsigned int j=0; j<vec.size();j++){  
				       if(vec.at(j)==tmpVec.at(0)){
					if(!std::count(index.begin(), index.end(), j)){
						index.push_back(j);
						break;
					}
				       }
				   }
				   tmpVec.erase(tmpVec.begin());
				}
			}else{index.push_back(-1);}
			return index;
			
		}

                int idxBestPN_MD_BBvsQCD(ROOT::VecOps::RVec<float> &ak8_probXbb, ROOT::VecOps::RVec<float> &ak8_probQCDbb, ROOT::VecOps::RVec<float> &ak8_probQCDcc, ROOT::VecOps::RVec<float> &ak8_probQCDb, ROOT::VecOps::RVec<float> &ak8_probQCDc, ROOT::VecOps::RVec<float> &ak8_probQCDothers)
                {
                        float max_s =ak8_probXbb.at(0)/(ak8_probXbb.at(0)+ak8_probQCDbb.at(0)+ak8_probQCDcc.at(0)+ak8_probQCDb.at(0)+ak8_probQCDc.at(0)+ak8_probQCDothers.at(0));
                        int idx_b= 0;
                        for (unsigned int i=1; i<ak8_probXbb.size(); i++)
                        {
                                float s =  ak8_probXbb.at(i)/(ak8_probXbb.at(i)+ak8_probQCDbb.at(i)+ak8_probQCDcc.at(i)+ak8_probQCDb.at(i)+ak8_probQCDc.at(i)+ak8_probQCDothers.at(i));
                                if(s>max_s)
                                {
                                        max_s = s;
                                        idx_b = i;
                                }
                        }

                        return idx_b;
                }
               
                int idxBestZqq(ROOT::VecOps::RVec<int> &flavour, ROOT::VecOps::RVec<float> &quark_eta, ROOT::VecOps::RVec<float> &quark_phi, ROOT::VecOps::RVec<float> &ak8_eta, ROOT::VecOps::RVec<float> &ak8_phi)
                {
                        if(flavour.size() != 2){
                                //cout<<"no 2 quarks"<<endl;
                                return -4;
                        }

                        if(abs(flavour.at(0))!=5){
                                return -2;
                        }

                        float min_DR0 = 0.8;
                        int idx_q0= -3;
                        for (unsigned int i=0; i<ak8_eta.size(); i++)
                        {
                                float dr0 = sqrt((quark_eta.at(0)-ak8_eta.at(i))*(quark_eta.at(0)-ak8_eta.at(i))+(quark_phi.at(0)-ak8_phi.at(i))*(quark_phi.at(0)-ak8_phi.at(i))); 
                                if(dr0<min_DR0)
                                {
                                        min_DR0 = dr0;
                                        idx_q0 = i;
                                }
                        }
  
                        float min_DR1 = 0.8;
                        int idx_q1= -3;
                        for (unsigned int i=0; i<ak8_eta.size(); i++)
                        {
                                float dr1 = sqrt((quark_eta.at(1)-ak8_eta.at(i))*(quark_eta.at(1)-ak8_eta.at(i))+(quark_phi.at(1)-ak8_phi.at(i))*(quark_phi.at(1)-ak8_phi.at(i))); 
                                if(dr1<min_DR1)
                                {
                                        min_DR1 = dr1;
                                        idx_q1 = i;
                                }
                        }

                        if(idx_q0 != idx_q1)
                        {
                        cout<<"quarks associated to different AK8 jets: "<<idx_q0<<" && "<<idx_q1<<endl;
                        return -1;
                        }

                        return idx_q0;
                }


		""")



		self.rdframe = self.rdframe.Filter("TriggerMatch(Trigger_hltname,Trigger_hltdecision)","TriggerMatch")\
		.Filter("noElectron(Ele_pt,Ele_eta,Ele_isPassID) && noMuon(Muon_pt,Muon_eta,Muon_PF_Iso_R04)","LeptonVeto")\
		.Filter("AK8PuppiJets_pt.size()>1","AK8nJet>2")\
		.Define("ptIndex","std::vector<int> ptIndex = idxOrder(AK8PuppiJets_pt); return ptIndex;")\
		.Define("DDXIndex","std::vector<int> DDXIndex = idxOrder(jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb); return DDXIndex;")\
		.Define("AK8_ptJet0","AK8PuppiJets_pt[ptIndex[0]]")\
		.Define("AK8_etaJet0","AK8PuppiJets_eta[ptIndex[0]]")\
                .Define("AK8_phiJet0","AK8PuppiJets_phi[ptIndex[0]]")\
		.Define("AK8_sdmJet0","AK8PuppiJets_softdropmass[ptIndex[0]]")\
		.Define("AK8_ptJet1","AK8PuppiJets_pt[ptIndex[1]]")\
		.Define("AK8_etaJet1","AK8PuppiJets_eta[ptIndex[1]]")\
                .Define("AK8_phiJet1","AK8PuppiJets_phi[ptIndex[1]]")\
		.Define("AK8_sdmJet1","AK8PuppiJets_softdropmass[ptIndex[1]]")\
                .Filter("noBtagAK4(AK4PuppiJets_pt,AK4PuppiJets_eta,AK4PuppiJets_phi,jet_pfDeepCSVAK4JetTags_probb,jet_pfDeepCSVAK4JetTags_probbb,AK8_etaJet0,AK8_phiJet0)","BtagAK4Veto")\
		.Define("PN_ZbbVsQCD","jet_pfParticleNetJetTags_probZbb.at(ptIndex[0])/(jet_pfParticleNetJetTags_probZbb.at(ptIndex[0])+jet_pfParticleNetJetTags_probQCDbb.at(ptIndex[0])+jet_pfParticleNetJetTags_probQCDcc.at(ptIndex[0])+jet_pfParticleNetJetTags_probQCDb.at(ptIndex[0])+jet_pfParticleNetJetTags_probQCDc.at(ptIndex[0])+jet_pfParticleNetJetTags_probQCDothers.at(ptIndex[0]))")\
		.Define("PNMD_XbbVsQCD","jet_pfMassDecorrelatedParticleNetJetTags_probXbb.at(ptIndex[0])/(jet_pfMassDecorrelatedParticleNetJetTags_probXbb.at(ptIndex[0])+jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb.at(ptIndex[0])+jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc.at(ptIndex[0])+jet_pfMassDecorrelatedParticleNetJetTags_probQCDb.at(ptIndex[0])+jet_pfMassDecorrelatedParticleNetJetTags_probQCDc.at(ptIndex[0])+jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers.at(ptIndex[0]))")\
		.Define("DDX_XbbVsQCD","jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb.at(ptIndex[0])")\
		.Define("Higher_DDX_XbbVsQCD","jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb.at(DDXIndex[0])")\
		.Define("sdm_Higher_DDX_XbbVsQCD","AK8PuppiJets_softdropmass.at(DDXIndex[0])")\
                .Filter("AK8_ptJet0>450. && abs(AK8_etaJet0)<2.4","AK8Jet0 (pt>450 GeV && eta<2.4)")\
                .Filter("AK8_ptJet1>200. && abs(AK8_etaJet1)<2.4","AK8Jet1 (pt>200 GeV && eta<2.4)")\
                .Define("idx_Higher_PN_MD_XbbVsQCD","idxBestPN_MD_BBvsQCD(jet_pfMassDecorrelatedParticleNetJetTags_probXbb, jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb, jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc, jet_pfMassDecorrelatedParticleNetJetTags_probQCDb, jet_pfMassDecorrelatedParticleNetJetTags_probQCDc,jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers)")\
                .Define("Higher_PN_MD_XbbVsQCD","jet_pfMassDecorrelatedParticleNetJetTags_probXbb.at(idx_Higher_PN_MD_XbbVsQCD)/(jet_pfMassDecorrelatedParticleNetJetTags_probXbb.at(idx_Higher_PN_MD_XbbVsQCD)+jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb.at(idx_Higher_PN_MD_XbbVsQCD)+jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc.at(idx_Higher_PN_MD_XbbVsQCD)+jet_pfMassDecorrelatedParticleNetJetTags_probQCDb.at(idx_Higher_PN_MD_XbbVsQCD)+jet_pfMassDecorrelatedParticleNetJetTags_probQCDc.at(idx_Higher_PN_MD_XbbVsQCD)+jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers.at(idx_Higher_PN_MD_XbbVsQCD))")\
                .Define("sdm_Higher_PN_MD_XbbVsQCD","AK8PuppiJets_softdropmass.at(idx_Higher_PN_MD_XbbVsQCD)")\
                .Define("idx_BestZqq_AK8","idxBestZqq(quark_flavour, quark_eta, quark_phi, AK8PuppiJets_eta, AK8PuppiJets_phi)")
