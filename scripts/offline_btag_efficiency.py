import os
import sys
import ROOT
from ROOT import *

import CMSDIJET.StatisticalTools.limit_configuration as limit_config

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/src/CMSDIJET/QCDAnalysis/python/")
import analysis_configuration_8TeV as analysis_config

qcd_samples = ["QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6", "QCD_Pt-120to170_TuneZ2star_8TeV_pythia6", "QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6", "QCD_Pt-170to300_TuneZ2star_8TeV_pythia6", "QCD_Pt-1800_TuneZ2star_8TeV_pythia6", "QCD_Pt-300to470_TuneZ2star_8TeV_pythia6", "QCD_Pt-470to600_TuneZ2star_8TeV_pythia6", "QCD_Pt-600to800_TuneZ2star_8TeV_pythia6", "QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6", "QCD_Pt-80to120_TuneZ2star_8TeV_pythia6"]
def qcd_efficiency(wp, eta):
	h_num = None
	h_den = None

	numerator_analysis = "NoTrigger_{}_{}".format(eta, wp)
	denominator_analysis = "NoTrigger_{}".format(eta)

	for qcd_sample in qcd_samples:
		numerator_file = TFile(analysis_config.get_b_histogram_filename(numerator_analysis, qcd_sample))
		denominator_file = TFile(analysis_config.get_b_histogram_filename(denominator_analysis, qcd_sample))
		this_h_num = numerator_file.Get("BHistograms/h_pfjet_mjj")
		this_h_den = denominator_file.Get("BHistograms/h_pfjet_mjj")
		print "[debug] Scale factor for {} numerator = {}".format(qcd_sample, 19700 * analysis_config.simulation.background_cross_sections[qcd_sample] / numerator_file.Get("BHistograms/h_sample_nevents").Integral())
		print "[debug] \tdenominator = {}".format(19700 * analysis_config.simulation.background_cross_sections[qcd_sample] / numerator_file.Get("BHistograms/h_sample_nevents").Integral())
		this_h_num.Scale(19700 * analysis_config.simulation.background_cross_sections[qcd_sample] / numerator_file.Get("BHistograms/h_sample_nevents").Integral())
		this_h_den.Scale(19700 * analysis_config.simulation.background_cross_sections[qcd_sample] / denominator_file.Get("BHistograms/h_sample_nevents").Integral())
		if h_num:
			h_num.Add(this_h_num)
		else:
			h_num = this_h_num.Clone()
			h_num.SetDirectory(0)
		if h_den:
			h_den.Add(this_h_den)
		else:
			h_den = this_h_den.Clone()
			h_den.SetDirectory(0)
		numerator_file.Close()
		denominator_file.Close()
	f_out = TFile(analysis_config.get_offline_btag_file(wp, eta), "RECREATE")
	h_eff = h_num.Clone()
	h_eff.SetName("h_offline_btag_eff")
	h_eff.Divide(h_num, h_den, 1., 1., "B")
	h_num.SetName("h_num")
	h_num.Write()
	h_den.SetName("h_den")
	h_den.Write()
	f_out.cd()
	h_eff.Write()
	f_out.Close()

def signal_efficiency():
	pass

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Make offline b-tag efficiency histograms")
	parser.add_argument("--qcd", action='store_true', help="QCD MC")
	parser.add_argument("--analyses", type=str, default="trigbbl_CSVTM,trigbbh_CSVTM", help="Analyses to run (comaa-separated)")
	args = parser.parse_args()

	if args.qcd:
		for analysis in args.analyses.split(","):
			if "bbl" in analysis or "eta1p7" in analysis:
				eta = "eta1p7"
			elif "bbh" in analysis or "eta2p2" in analysis:
				eta = "eta2p2"
			else:
				print "ERROR : Can't figure out eta from analysis name"
				sys.exit(1)
			if "CSVTM" in analysis:
				wp = "CSVTM"
			else:
				print "ERROR : Only CSVTM implemented."
				sys.exit(1)
			qcd_efficiency(wp, eta)





