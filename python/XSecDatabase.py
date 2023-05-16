"""
cross-section database: all the xsec are in pb
map sampleName-xsec, -1 if no xsec is known
"""


class XSecDatabase(object):

	def __init__(self):
		self.database = {}

		##ZJetsToQQ
		self.database["ZJetsToQQ_HT200to400_TuneCP5_13TeV-madgraphMLM-pythia8"] = 1029
		self.database["ZJetsToQQ_HT400to600_TuneCP5_13TeV-madgraphMLM-pythia8"] = 101.1
		self.database["ZJetsToQQ_HT600to800_TuneCP5_13TeV-madgraphMLM-pythia8"] = 25.81
		self.database["ZJetsToQQ_HT800toInf_TuneCP5_13TeV-madgraphMLM-pythia8"] = 13.20 	

		##QCD
		self.database["QCD_PT-120to170_TuneCP5_13p6TeV_pythia8"] = 442400
		self.database["QCD_PT-170to300_EMEnriched_TuneCP5_13p6TeV_pythia8"] = 17980
		self.database["QCD_PT-300to470_TuneCP5_13p6TeV_pythia8"] = 7616
		self.database["QCD_PT-470to600_TuneCP5_13p6TeV_pythia8"] = 625.6
		self.database["QCD_PT-600to800_TuneCP5_13p6TeV_pythia8"] = 180.1
		self.database["QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8"] = 30.63
		self.database["QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8"] = 8.92
		self.database["QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8"] = 0.8075
		self.database["QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8"] = 0.1149
		self.database["QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8"] = 0.007619
		self.database["QCD_PT-3200toInf_TuneCP5_13p6TeV_pythia8"] = -1.
