import os
import sys
import array
from array import array
import ROOT
from ROOT import *
gROOT.SetBatch(True)
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
seaborn = Root.SeabornInterface()
seaborn.Initialize()
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_7_4_15/python/CMSDijet/StatisticalTools")
import limitPaths as limit_paths

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_7_4_15/python/MyTools/RootUtils")
import histogram_tools

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/CMSDIJET/QCDAnalysis")
import analysis_configuration_8TeV as analysis_config
import simulation_configuration_8TeV as simulation_config

import fits.py as fits

def bias_test(x, gen_pdf, test_pdf, seed, trials, output_file):
	RooRandom.generator().setSeed(seed)
	manager = RooMCStudy(gen_pdf, test_pdf, x, "", "")
	manager.generateAndFit(trials)
	output_file.cd()
	w = RooWorkspace("w", "workspace")
	getattr(w,'import')(manager)

def condor_bias_tests()


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Run bias tests")
	parser.add_argument("--run", action='store_true', help="Run bias test")
	parser.add_argument("--condor", action='store_true', help="Setup and run condor jobs")
	parser.add_argument("--njobs", type=int, default=50, help="Number of jobs")
	parser.add_argument("--ntoys", type=int, default=100, help="Number of toys per job")

	parser.add_argument("--gen_workspace", type=str, nargs=2 help="Path, PDF name of generating PDF")
	parser.add_argument("--test_workspace", type=str, nargs=2 help="Path, PDF name of generating PDF")

	parser.add_argument("--min_mjj", dest="min_mjj",
						default=0, type=int,
						help="Min m(jj) considered (default: %(default)s)",
						metavar="MASS_MIN")
	parser.add_argument("--max_mjj", dest="max_mjj",
						default=2000, type=int,
						help="Max m(jj) considered (default: %(default)s)",
						metavar="MASS_MAX")
	parser.add_argument("--seed", type=int, help="Set random seed ")
	args = parser.parse_args()

	if args.run:
		bias_test()

