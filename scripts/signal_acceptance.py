import os
import sys
import pickle
from ROOT import *
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/CMSDIJET/QCDAnalysis")
import analysis_configuration_8TeV as analysis_config
import CMSDIJET.StatisticalTools.trigger_efficiency as trigger_efficiency

models = ["Hbb", "RSG", "ZPrime"]
#analyses = ["trigbbl_CSVTM", "trigbbh_CSVTM"] # , "trigbbl_CSVM", "trigbbh_CSVM"
analyses = []
for sr in ["trigbbl", "trigbbh"]:
	#for wp in ["CSVT", "CSVM", "CSVL", "CSVTL", "CSVML", "CSVTM"]:
	for wp in ["CSVTM"]:
		analysis = sr + "_" + wp
		analyses.append(analysis)

#masses = [325, 350, 400, 500, 600, 750, 900, 1200]
#masses = {"trigbbl_CSVTM":[400, 500, 600, 750, 900], "trigbbh_CSVTM":[600, 750, 900, 1200]}
masses = {"trigbbl_CSVTM":[325,350,400,500,600,750,900], "trigbbh_CSVTM":[600,750,900,1200]}

signal_acc_times_eff = {}

use_MC_trigger = True

for model in models:
	signal_acc_times_eff[model] = {}
	for analysis in analyses:
		if analysis == "trigbbl_CSVTM":
			notrig_analysis = "NoTrigger_eta1p7_CSVTM"
		elif analysis == "trigbbh_CSVTM":
			notrig_analysis = "NoTrigger_eta2p2_CSVTM"


        #signal_pdf_file = analysis_config.get_signal_fit_file(notrig_analysis, args.model, mass, "bukin", interpolated=(not mass in analysis_config.simulation.simulated_masses))
		signal_acc_times_eff[model][analysis] = {}
		if "bbl" in analysis:
			mjj_range = [296, 1058]
		elif "bbh" in analysis:
			mjj_range = [526, 1607]

		for mass in masses[analysis]:
			if use_MC_trigger:
				f = TFile(analysis_config.get_b_histogram_filename(analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")))
			else:
				f = TFile(analysis_config.get_b_histogram_filename(notrig_analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")))
			if not f.IsOpen():
				signal_acc_times_eff[model][analysis][mass] = 0.
				continue
			low_bin = f.Get("BHistograms/h_pfjet_mjj").GetXaxis().FindBin(mjj_range[0] + 1.e-5)
			high_bin = f.Get("BHistograms/h_pfjet_mjj").GetXaxis().FindBin(mjj_range[1] - 1.e-5)
			numerator = f.Get("BHistograms/h_pfjet_mjj").Integral(low_bin, high_bin)
			if "ZPrime" in model:
				denominator = f.Get("BHistograms/h_input_nevents").Integral()
				print "{} / {} den = {}".format(model, analysis, denominator)
			else:
				denominator = f.Get("BHistograms/h_sample_nevents").Integral()

			if denominator > 0:
				signal_acc_times_eff[model][analysis][mass] = numerator/denominator
			else:
				signal_acc_times_eff[model][analysis][mass] = 0.
			if not use_MC_trigger:
				signal_acc_times_eff[model][analysis][mass] *= trigger_efficiency.online_btag_eff[analysis][0]
			print "{} / {} / {} GeV : mjj acceptance = {}".format(model, analysis, mass, f.Get("BHistograms/h_pfjet_mjj").Integral(low_bin, high_bin) / f.Get("BHistograms/h_pfjet_mjj").Integral())
			print "\tsample nevents = {}".format(f.Get("BHistograms/h_sample_nevents").Integral())
			print "\tinput_nevents = {}".format(f.Get("BHistograms/h_input_nevents").Integral())
			print "\tinput_nevents_weighted = {}".format(f.Get("BHistograms/h_input_nevents_weighted").Integral())
			print "\tpass_nevents = {}".format(f.Get("BHistograms/h_pass_nevents").Integral())
			print "\tpass_nevents_weighted = {}".format(f.Get("BHistograms/h_pass_nevents_weighted").Integral())

print signal_acc_times_eff
print "Saving to {}".format(analysis_config.simulation.get_signal_AE_filename())
if use_MC_trigger:
	pickle.dump(signal_acc_times_eff, open(analysis_config.simulation.get_signal_AE_filename() + ".MCtrigger", "wb"))
else:
	pickle.dump(signal_acc_times_eff, open(analysis_config.simulation.get_signal_AE_filename(), "wb"))


# Print signal cutflow table
header_formatting = {
	"MaxMetOverSumEt":"$E_{\\mathrm{T}}^{\\mathrm(miss)}/E_{\\mathrm{T}}<0.5$",
	"GoodPFDijet":"$\\geq2$ good jets",
	"MinNCSVM":"$N_{CSVM}\\geq2$",
	"MinNCSVT":"$N_{CSVT}\\geq1$",
	"MinLeadingPFJetPt":"Leading $p_{\\mathrm{T}}$",
	"MinSubleadingPFJetPt":"Subleading $p_{\\mathrm{T}}$",
	"PFDijetMaxDeltaEta":"$\\Delta\\eta<1.3$"
}


for analysis in analyses:
	if analysis == "trigbbl_CSVTM":
		notrig_analysis = "NoTrigger_eta1p7_CSVTM"
	elif analysis == "trigbbh_CSVTM":
		notrig_analysis = "NoTrigger_eta2p2_CSVTM"
	headers = []
	rows = {}
	first = True
	for model in models:
		rows[model] = {}
		if "bbl" in analysis:
			mjj_range = [296, 1058]
		elif "bbh" in analysis:
			mjj_range = [526, 1607]

		for mass in masses[analysis]:
			if use_MC_trigger:
				f = TFile(analysis_config.get_b_histogram_filename(analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")))
			else:
				f = TFile(analysis_config.get_b_histogram_filename(notrig_analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")))
			if not f.IsOpen():
				print "ERROR : Failed to open {}".format(analysis_config.get_b_histogram_filename(notrig_analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")))
				sys.exit(1)
			cutflow_histogram = f.Get("BHistograms/CutFlowCounter_QCDEventSelector")
			if first:
				for bin in xrange(1, cutflow_histogram.GetNbinsX()):
					bin_label = cutflow_histogram.GetXaxis().GetBinLabel(bin)
					if bin_label in header_formatting:
						headers.append(header_formatting[bin_label])
					else:
						headers.append(bin_label)
				headers.append("b-tag #epsilon")
				headers.append("m_{jj} acceptance")
			rows[model][mass] = {}
			inclusive = f.Get("BHistograms/h_sample_nevents").Integral()
			# For ZPrime, h_sample_nevents includes charms decays, which you have decided to ignore.
			if model == "ZPrime":
				inclusive = inclusive / 2.
			if inclusive <= 0:
				print "ERROR : Cutflow histograms inclusive bin has zero entries"
				sys.exit(1)
			rows[model][mass]["Inclusive"] = 1.0
			btag_sf = 1.
			for bin in xrange(2, cutflow_histogram.GetNbinsX()):
				bin_label = cutflow_histogram.GetXaxis().GetBinLabel(bin)
				if bin_label == "GoodPFDijet":
					btag_sf = f.Get("BHistograms/h_pass_nevents_weighted").Integral() / f.Get("BHistograms/h_pass_nevents").Integral()
					print ""
				if bin_label in header_formatting:
					bin_label = header_formatting[bin_label]
				rows[model][mass][bin_label] = cutflow_histogram.GetBinContent(bin) / inclusive * btag_sf
			final_cutflow_bin = cutflow_histogram.GetBinContent(cutflow_histogram.GetNbinsX() - 1) * btag_sf
			mjj_hist_entries = f.Get("BHistograms/h_pfjet_mjj").Integral(0, f.Get("BHistograms/h_pfjet_mjj").GetNbinsX() + 1)
			print "Final cutflow bin = {}, mjj hist entries = {}".format(final_cutflow_bin, mjj_hist_entries)

			rows[model][mass]["b-tag #epsilon"] = final_cutflow_bin / inclusive * trigger_efficiency.online_btag_eff[analysis][0]
			low_bin = f.Get("BHistograms/h_pfjet_mjj").GetXaxis().FindBin(mjj_range[0] + 1.e-5)
			high_bin = f.Get("BHistograms/h_pfjet_mjj").GetXaxis().FindBin(mjj_range[1] - 1.e-5)
			numerator = f.Get("BHistograms/h_pfjet_mjj").Integral(low_bin, high_bin)
			denominator = f.Get("BHistograms/h_pfjet_mjj").Integral(0, f.Get("BHistograms/h_pfjet_mjj").GetNbinsX() + 1)
			rows[model][mass]["m_{jj} acceptance"] = final_cutflow_bin / inclusive * trigger_efficiency.online_btag_eff[analysis][0] * numerator / denominator

			first = False
	print "\\begin{table}"
	print "\t\\centering"
	print "\t\\resizebox{7in}{!}{"
	print "\t\\begin{tabular}{|c|",
	for header in headers:
		print "c|",
	print "}"
	print "\t\t\\hline"
	print "\t\tModel/Mass\t&",
	for i, header in enumerate(headers):
		print "\t" + header + "\t",
		if i != len(headers) - 1:
			print "&",
	print "\\\\"
	print "\t\t\\hline"
	for mass in masses[analysis]:
		for model in models:
			print "\t\t{}/{}\t".format(model, mass),
			for header in headers:
				print "&\t{:.2f}\t".format(rows[model][mass][header] * 100),
			print "\\\\"
			print "\t\t\\hline"
		print "\t\t\\hline"
	print "\t\\end{tabular}"
	print "\t}"
	print "\t\\caption{{Percent of signal events remaining after each cut for the {} signal region.}}".format(analysis)
	print "\t\\label{{table:cut-efficiency-{}}}".format(analysis)
	print "\\end{table}"

#print "\\begin{table}"
#print "\t\\centering"
#print "\t\\begin{tabular}{|c|",
#for mass in masses[analysis]:
#	print "c|",
#print "}"
#print "\t\t\\hline"
#print "\t\tModel / Signal Region "
#for mass in masses[analysis]:
#	print "\t&\t {} GeV".format(mass),
#print "\\\\"
#print "\t\t\\hline"
#for model in models:
#	for analysis in analyses:
#		print "\t\t{} / {}".format(model, analysis),
#		for mass in masses[analysis]:
#			print "\t&\t{:.2f}\\%".format(signal_acc_times_eff[model][analysis][mass]*100),
#		print "\\\\"
#		print "\t\t\\hline"
#print "\t\\end{tabular}"
#print "\\caption{Efficiency of the event selection on signal events.}"
#print "\\label{table:total-signal-efficiency}"
#print "\\end{table}"