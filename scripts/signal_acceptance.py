import os
import sys
import pickle
from ROOT import *
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/CMSDIJET/QCDAnalysis")
import analysis_configuration_8TeV as analysis_config

models = ["Hbb", "RSG", "ZPrime"]
analyses = ["trigbbl_CSVTM", "trigbbh_CSVTM"] # , "trigbbl_CSVM", "trigbbh_CSVM"
masses = [350, 400, 500, 600, 750, 900, 1200]
#masses = {"trigbbl_CSVTM":[400, 500, 600, 750, 900], "trigbbh_CSVTM":[600, 750, 900, 1200]}

signal_acc_times_eff = {}

for model in models:
	signal_acc_times_eff[model] = {}
	for analysis in analyses:
		signal_acc_times_eff[model][analysis] = {}
		if "bbl" in analysis:
			mjj_range = [296, 1246]
		elif "bbh" in analysis:
			mjj_range = [526, 1455]
		for mass in masses:
			f = TFile(analysis_config.get_b_histogram_filename(analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")))
			if not f.IsOpen():
				signal_acc_times_eff[model][analysis][mass] = 0.
				continue
			low_bin = f.Get("BHistograms/h_pfjet_mjj").GetXaxis().FindBin(mjj_range[0] + 1.e-5)
			high_bin = f.Get("BHistograms/h_pfjet_mjj").GetXaxis().FindBin(mjj_range[1] - 1.e-5)
			numerator = f.Get("BHistograms/h_pfjet_mjj").Integral(low_bin, high_bin)
			denominator = f.Get("BHistograms/h_sample_nevents").Integral()
			if denominator > 0:
				signal_acc_times_eff[model][analysis][mass] = numerator/denominator
			else:
				signal_acc_times_eff[model][analysis][mass] = 0.
			print "{} / {} / {} GeV : mjj acceptance = {}".format(model, analysis, mass, f.Get("BHistograms/h_pfjet_mjj").Integral(low_bin, high_bin) / f.Get("BHistograms/h_pfjet_mjj").Integral())
			print "\tsample nevents = {}".format(f.Get("BHistograms/h_sample_nevents").Integral())
			print "\tinput_nevents = {}".format(f.Get("BHistograms/h_input_nevents").Integral())
			print "\tinput_nevents_weighted = {}".format(f.Get("BHistograms/h_input_nevents_weighted").Integral())
			print "\tpass_nevents = {}".format(f.Get("BHistograms/h_pass_nevents").Integral())
			print "\tpass_nevents_weighted = {}".format(f.Get("BHistograms/h_pass_nevents_weighted").Integral())




print signal_acc_times_eff
pickle.dump(signal_acc_times_eff, open(analysis_config.simulation.get_signal_AE_filename(), "wb"))