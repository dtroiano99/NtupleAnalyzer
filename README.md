NtupleAnalyzer
====================================
package to run different analysis configurations with the same class. Common classes defined in the package 
- `Handler.py`: main class to define RDataFrame object reading all the root files stored in a given directory; the following methods are defined:
	- `cut`: cut="preselection string" defined a preselection to be applied to all the events.
	- `AddLumiWeight(LumiTarget="1.",xsecName="")`: add a luminosity weight (labelled "lumiWeight") to the main tree.
	- `Histo1D(hName,var,weight,Nbins,Xmin,Xmax,cut="true")`: define a 1D histogram after having specified the weight to be applied (`weight`), interval and numbers of bins, additional cut to be applied.
	- `Histo2D(hName,var,weight,NbinsX,Xmin,Xmax,NbinsY,Ymin,Ymax,cut="true")`: same as `Histo1D` but for a 2D histo.
	- `WriteHisto(outFileName,fMode="RECREATE")`: write crated histos in a output root file.
	- `WriteTree(outTreeName,outFileName,columns=None,fMode="RECREATE")`: write a tree with the leaves specified in `columns` (if `None` is specified the same leaves of input tree plus the ones defined by the user selection are stored) in a new root file. Be careful that in this method there is a hack to let RDataFrame store the std::vector<string> objects in the output file; if needed modified this method.
- `XSecDataBase`: class storing a dictionary of the cross-section of the various samples used in the analysis.

Example of a user class:
- `Analysis_Ak8Tagger.py`: child class of `Handler.py` used to skim the initial tree, to add new variables and to produce 1D/2D histograms for AK8Analysis.

------------------------------------
Setup
------------------------------------

The code has been tested with `python3` and `root 6.26`. To run the code just simply source `setup.sh`.


------------------------------------
How to run
------------------------------------
Under scripts there is the basic example `submitAna_example.py` to process a MC dataset with `Analysis_Ak8Tagger.py`. 

```
submitAna.py -i /lustre/cms/store/user/azaza/QCD_PT-600to800_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-600to800_TuneCP5_13p6TeV_pythia8_Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2/230405_184254/ -L 8924. -s QCD_PT-600to800_TuneCP5_13p6TeV_pythia8 -N 100 |& tee output_test

INFO:QCD_PT-600to800_TuneCP5_13p6TeV_pythia8:Reading files: /lustre/cms/store/user/azaza/QCD_PT-600to800_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-600to800_TuneCP5_13p6TeV_pythia8_Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2/230405_184254/
INFO:QCD_PT-600to800_TuneCP5_13p6TeV_pythia8:- LumiTarget=8924.
INFO:QCD_PT-600to800_TuneCP5_13p6TeV_pythia8:- TotEvents=7840148.0
INFO:QCD_PT-600to800_TuneCP5_13p6TeV_pythia8:- xsecDB=180.1
INFO:QCD_PT-600to800_TuneCP5_13p6TeV_pythia8:- lumiweight=0.2049977117778899
Preselection: pass=100        all=100        -- eff=100.00 % cumulative eff=100.00 %
TriggerMatch: pass=100        all=100        -- eff=100.00 % cumulative eff=100.00 %
LeptonVeto: pass=100        all=100        -- eff=100.00 % cumulative eff=100.00 %
AK8nJet>2 : pass=100        all=100        -- eff=100.00 % cumulative eff=100.00 %
AK8Jet0 (pt>450 GeV && eta<2.4): pass=99         all=100        -- eff=99.00 % cumulative eff=99.00 %
AK8Jet1 (pt>200 GeV && eta<2.4): pass=98         all=99         -- eff=98.99 % cumulative eff=98.00 %
All outputs produced in 0.0h:0.0m:26.9s
```

The code can be run in multithread with the option `-MT` (this option can not be run together with `-N`). More option are available in `submit_example.py`:
```
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
parser.add_argument("-s","--SampleName", help="sample name used in the xsec database",default="sample")
parser.add_argument("-N","--NEvents", help="select NEvents to be processed", default="-1")
```
