#!/usr/bin/env python3
"""
basic code to submit a specific analysis case - i.e. AK8Analysis
"""
import argparse, sys, time, logging
import ROOT
from Analysis_Ak8Tagger13 import Analysis_Ak8Tagger
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


	AnaAK8.ApplyAK8Selection13()

	it0=-1
	cutted = AnaAK8.rdframe.GetColumnNames()
	for x in cutted:
                        it0=it0+1
                        #print(x)
                        #print(columns.at(it0))
                        if x=="AK8PuppiJets_Lpt_softdropmass_Down.first":
                                cutted.erase(cutted.begin()+it0)
                                cutted.erase(cutted.begin()+it0)
                                #print(columns.at(it0))
                                break
	it0=-1
	for x in cutted:
                        it0=it0+1
                        #print(x)
                        #print(columns.at(it0))
                        if x=="AK8PuppiJets_Lpt_softdropmass_Up.first":
                                cutted.erase(cutted.begin()+it0)
                                cutted.erase(cutted.begin()+it0)
                                #print(columns.at(it0))
                                break


	AnaAK8.WriteTree(inputs.OutputTreeName,inputs.OutputFileTreeName, cutted)
	print("All outputs produced in "+timeStr(time.time()-StartTime))
