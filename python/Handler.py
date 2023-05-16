"""
classes to handle ntuplizer output, skimming and adding new variables
"""

import ROOT
from ROOT import gROOT,RDataFrame,TFile,TChain,SetOwnership
from pathlib import Path
import logging,sys,os
logging.basicConfig(level=logging.INFO)


RDF = ROOT.ROOT.RDataFrame

class Handler(object):

	def __init__(self, name, nameTree, pathInput, cut='true', NEvents=-1,logLevel=logging.INFO):

		self.name = name
		self.nameTree = nameTree
		self.pathInput = pathInput
		self.cut = "true"
		self.HistoNEventsName = "Ana/nEvents"
		self.ListFile = []
		self.ListHisto = {}
		self.logger = logging.getLogger(name)
		self.logger.setLevel(logLevel)

		
		tc = TChain(self.nameTree)
		self.logger.info(f"Reading files: {self.pathInput}")
		for fName in Path(self.pathInput).rglob('*.root'): 
			if self.logger.level < logging.INFO: self.logger.debug(f" - {fName}")
			self.ListFile.append(str(fName))
			tc.Add(str(fName))

		SetOwnership(tc, False)
		rdSample = RDF(tc)
		SetOwnership(rdSample, False)
		rdSample = rdSample.Filter(self.cut,"Preselection")
		if not NEvents==-1: rdSample = rdSample.Range(NEvents)
		self.rdframe = rdSample 


	def AddLumiWeight(self,LumiTarget="1.",xsecName=""):

		from XSecDatabase import XSecDatabase
		xsecName= xsecName if not xsecName=="" else self.name
		if xsecName in XSecDatabase().database.keys():
			xsecDB = XSecDatabase().database[xsecName]

			TotEvents=0. 
			for fName in self.ListFile:
				f=TFile(fName)
				h = f.Get(self.HistoNEventsName)	
				TotEvents+=h.Integral()
				del h
				f.Close()
			self.logger.info(f"- LumiTarget={LumiTarget}")
			self.logger.info(f"- TotEvents={TotEvents}")
			self.logger.info(f"- xsecDB={xsecDB}")
			self.logger.info(f"- lumiweight={float(LumiTarget)/(TotEvents/float(xsecDB))}")
			self.rdframe = self.rdframe.Define("lumiWeight",f"double lumiweight={LumiTarget}/({TotEvents}/{xsecDB}); return lumiweight;")

		else:
			self.logger.error(f"Asked for lumiWeight but no xsec found in the database for {xsecName}")

	def Histo1D(self,hName,var,weight,Nbins,Xmin,Xmax,cut="true"):

		rdf = self.rdframe.Filter(cut)
		histo = rdf.Histo1D(ROOT.ROOT.RDF.TH1DModel(hName,hName,Nbins,Xmin,Xmax),var,weight).GetValue() 
		SetOwnership(histo, False)
		self.ListHisto[hName] = histo

	def Histo2D(self,hName,var,weight,NbinsX,Xmin,Xmax,NbinsY,Ymin,Ymax,cut="true"):

		histo = self.rdframe.Filter(cut).Histo1D(ROOT.ROOT.RDF.TH2DModel(hName,hName,NbinsX,Xmin,Xmax,NbinsY,Ymin,Ymax),var,weight).GetValue() 
		SetOwnership(histo, False)
		self.ListHisto[hName] = histo

	def WriteTree(self,outTreeName,outFileName,columns=None,fMode="RECREATE"):
		
		if columns==None: columns=self.rdframe.GetColumnNames()
		snapshotOptions = ROOT.RDF.RSnapshotOptions()
		snapshotOptions.fMode=fMode
		####
		#hack for ROOT 6.26 and previous - not supported conversion for RDF::Vec<string> to std::<string>
		self.rdframe = self.rdframe.Redefine("Trigger_hltname","std::vector<string> Trigger_hltname_v2; for(unsigned int i=0; i<Trigger_hltname.size(); i++){Trigger_hltname_v2.push_back(Trigger_hltname.at(i));} return Trigger_hltname_v2;")\
		.Redefine("Trigger_l1name","std::vector<string> Trigger_l1name_v2; for(unsigned int i=0; i<Trigger_l1name.size(); i++){Trigger_l1name_v2.push_back(Trigger_l1name.at(i));} return Trigger_l1name_v2;")

		####
		self.rdframe.Snapshot(outTreeName,outFileName,columns,snapshotOptions)
		cutflow_report = self.rdframe.Report()
		cutflow_report.Print()


	def WriteHisto(self,outFileName,fMode="RECREATE"):

		fOut = TFile(outFileName,fMode)
		fOut.cd()
		for k in self.ListHisto.keys():
			self.ListHisto[k].Write()
		fOut.Close()
		
	def __del__(self):
		del self.ListHisto
		del self.rdframe
		
