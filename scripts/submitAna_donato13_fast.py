#!/usr/bin/env python3
"""
basic code to submit a specific analysis case - i.e. AK8Analysis
"""
import argparse, sys, time, logging
import ROOT
from Analysis_Ak8Tagger13_fast import Analysis_Ak8Tagger
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F

def timeStr(seconds):

	hours=seconds//3600
	seconds=seconds-hours*3600
	minutes=seconds//60
	seconds=seconds-minutes*60

	return f"{hours}h:{minutes}m:{int(seconds*10+0.5)/10.}s"


if __name__ == "__main__":


	parser = argparse.ArgumentParser(description="submitAna options:")
	parser.add_argument("-l", "--logLevel", help="set log level", default="INFO", choices=["DEBUG", "INFO"])
	parser.add_argument("-MT","--MultiThread", help="enable multithread", action="store_true")
	parser.add_argument("-i","--InputPath", help="Input root files path",default="")
	parser.add_argument("-t","--TreeName", help="name of the tree to be processed, i.e. Ana/passedEvents",default="Ana/passedEvents")
	parser.add_argument("-H","--HistoName", help="name of the histo stoting the total number of events, i.e. Ana/nEvents",default="Ana/nEvents")
	parser.add_argument("-ot","--OutputTreeName", help="tree name to be used for the output",default="passedEvents")
	parser.add_argument("-otf","--OutputFileTreeName", help="file name of the output file to store trees",default="Output_tree.root")
	parser.add_argument("-ohf","--OutputFileHistoName", help="file name of the output file to store histos",default="Output_histo.root")
	parser.add_argument("-p","--PreselectionCut", help="preselection cut to be applied",default="")
	parser.add_argument("-L","--Lumi", help="luminosity to target in pb",default="-1")
	parser.add_argument("-ipEE","--ispreEE", help="is era C or era D",default="True")
	parser.add_argument("-s","--SampleName", help="sample name used in the xsec database",default="sample")
	parser.add_argument("-N","--NEvents", help="select NEvents to be processed", default="-1")

	inputs = parser.parse_args()

	StartTime = time.time()
	logLevel = logging.DEBUG if inputs.logLevel=="DEBUG" else logging.INFO
	if inputs.MultiThread and inputs.NEvents=="-1": ROOT.ROOT.EnableImplicitMT()
	elif inputs.MultiThread and not inputs.NEvents=="-1": print("Multithread can not be run with NEvents option")
	if inputs.InputPath=="": 
		print("no inputs provided... exiting")
		sys.exit(-1)
	
	AnaAK8 = Analysis_Ak8Tagger(inputs.SampleName,inputs.TreeName,inputs.InputPath,NEvents=int(inputs.NEvents))
	weight="1."
	if not inputs.Lumi=="-1": 
		AnaAK8.HistoNEventsName=inputs.HistoName
		AnaAK8.AddLumiWeight(inputs.Lumi,inputs.SampleName)
		weight="lumiWeight"
	if not inputs.PreselectionCut == "": 
		AnaAK8.cut = inputs.PreselectionCut
		
	h_vetomap = ROOT.TH2F("h_vetomap","hvetomap", 100, 0, 100, 100, 0, 100)

	if inputs.ispreEE == "True":
		#print("vero")
		inFile = ROOT.TFile.Open("/lustre/cms/store/user/dtroiano/Commissioning/Utilities/Summer22_23Sep2023_RunCD_v1.root","READ")
		h_vetomap = inFile.Get("jetvetomap")
	else :
		#print("falso")
		inFile = ROOT.TFile.Open("/lustre/cms/store/user/dtroiano/Commissioning/Utilities/Summer22EE_23Sep2023_RunEFG_v1.root","READ")
		h_vetomap = inFile.Get("jetvetomap")

	#print(h_vetomap.GetBinContent(38,39))

	AnaAK8.ApplyAK8Selection13(h_vetomap)

	#AnaAK8.Histo1D("hJet0_Fullsdm","AK8_sdmJet0",weight,30,40,160)
	#AnaAK8.Histo1D("hJet0_Windowsdm","AK8_sdmJet0",weight,30,40,160,"AK8_sdmJet0>80 && AK8_sdmJet0<110")
	#AnaAK8.WriteHisto(inputs.OutputFileHistoName)
	AnaAK8.WriteTree(inputs.OutputTreeName,inputs.OutputFileTreeName)
	print("All outputs produced in "+timeStr(time.time()-StartTime))
