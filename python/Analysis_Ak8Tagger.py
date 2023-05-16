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
		""")



		self.rdframe = self.rdframe.Filter("TriggerMatch(Trigger_hltname,Trigger_hltdecision)","TriggerMatch")\
		.Filter("noElectron(Ele_pt,Ele_eta,Ele_isPassID) && noMuon(Muon_pt,Muon_eta,Muon_PF_Iso_R04)","LeptonVeto")\
		.Filter("AK8PuppiJets_pt.size()>1","AK8nJet>2")\
		.Define("ptIndex","std::vector<int> ptIndex = idxOrder(AK8PuppiJets_pt); return ptIndex;")\
		.Define("AK8_ptJet0","AK8PuppiJets_pt[ptIndex[0]]")\
		.Define("AK8_etaJet0","AK8PuppiJets_eta[ptIndex[0]]")\
		.Define("AK8_sdmJet0","AK8PuppiJets_softdropmass[ptIndex[0]]")\
		.Define("AK8_ptJet1","AK8PuppiJets_pt[ptIndex[1]]")\
		.Define("AK8_etaJet1","AK8PuppiJets_eta[ptIndex[1]]")\
		.Define("AK8_sdmJet1","AK8PuppiJets_softdropmass[ptIndex[1]]")\
		.Define("PN_ZbbVsQCD","jet_pfParticleNetJetTags_probZbb.at(ptIndex[0])/(jet_pfParticleNetJetTags_probZbb.at(ptIndex[0])+jet_pfParticleNetJetTags_probQCDbb.at(ptIndex[0])+jet_pfParticleNetJetTags_probQCDcc.at(ptIndex[0])+jet_pfParticleNetJetTags_probQCDb.at(ptIndex[0])+jet_pfParticleNetJetTags_probQCDc.at(ptIndex[0])+jet_pfParticleNetJetTags_probQCDothers.at(ptIndex[0]))")\
		.Define("PNMD_XbbVsQC","jet_pfMassDecorrelatedParticleNetJetTags_probXbb.at(ptIndex[0])/(jet_pfMassDecorrelatedParticleNetJetTags_probXbb.at(ptIndex[0])+jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb.at(ptIndex[0])+jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc.at(ptIndex[0])+jet_pfMassDecorrelatedParticleNetJetTags_probQCDb.at(ptIndex[0])+jet_pfMassDecorrelatedParticleNetJetTags_probQCDc.at(ptIndex[0])+jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers.at(ptIndex[0]))")\
		.Define("DDX_XbbVsQCD","jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb.at(ptIndex[0])")\
		.Filter("AK8_ptJet0>450. && AK8_etaJet0<2.4","AK8Jet0 (pt>450 GeV && eta<2.4)")\
		.Filter("AK8_ptJet1>200. && AK8_etaJet1<2.4","AK8Jet1 (pt>200 GeV && eta<2.4)")

