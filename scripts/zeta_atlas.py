#!/usr/bin/env python
import os
import sys
import re
import math
import time
from argparse import ArgumentParser
from array import array
import numpy as np
import CMS_lumi
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/src/CMSDIJET/QCDAnalysis/python/")
import analysis_configuration_8TeV as analysis_config

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_7_4_15/python/MyTools/RootUtils")
import histogram_tools

import ROOT
from ROOT import *
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
seaborn = Root.SeabornInterface()
seaborn.Initialize()
Root.SetCanvasStyle()

# Get sigmatilde
import pickle
sigmatilde_file = open(os.path.expandvars("$CMSSW_BASE/src/CMSDIJET/StatisticalTools/data/sigmatildeATLAS.pkl"), 'r')
sigmatilde = pickle.load(sigmatilde_file)
print sigmatilde
sigmatilde_file.close()

# My values, taken from the tables in the paper
#acceptances = {
#	1500:0.085,
#	2000:0.072,
#	3000:0.035,
#	4000:0.025,
#	5000:0.021,
#}

# Kirtimaan's values
acceptances = {
	1500:0.0993114,
	2000:0.0742023,
	3000:0.0264772,
	4000:0.0173586,
	5000:0.00952033,
}
