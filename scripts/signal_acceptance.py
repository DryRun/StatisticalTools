import os
import sys
import pickle
from ROOT import *
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/CMSDIJET/QCDAnalysis")
import analysis_configuration_8TeV as analysis_config

models = ["Hbb", "RSG"]
analyses = ["trigbbl_CSVTM", "trigbbh_CSVTM", "trigbbl_CSVM", "trigbbh_CSVM"]
masses = [400, 500, 600, 750, 900, 1200]
#masses = {"trigbbl_CSVTM":[400, 500, 600, 750, 900], "trigbbh_CSVTM":[600, 750, 900, 1200]}

signal_acc_times_eff = {}

for model in models:
	signal_acc_times_eff[model] = {}
	for analysis in analyses:
		signal_acc_times_eff[model][analysis] = {}
		for mass in masses:
			f = TFile(analysis_config.get_b_histogram_filename(analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")))
			numerator = f.Get("BHistograms/h_pass_nevents").Integral()
			denominator = f.Get("BHistograms/h_sample_nevents").Integral()
			if denominator > 0:
				signal_acc_times_eff[model][analysis][mass] = numerator/denominator
			else:
				signal_acc_times_eff[model][analysis][mass] = 0.
print signal_acc_times_eff
pickle.dump(signal_acc_times_eff, open(analysis_config.simulation.get_signal_AE_filename(), "wb"))