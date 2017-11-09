import os
import sys
import ROOT
from array import array
from ROOT import *
import math
from math import sqrt

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/src/CMSDIJET/QCDAnalysis/python/")
import analysis_configuration_8TeV as analysis_config

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_7_4_15/python/MyTools/RootUtils")
import histogram_tools

gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
seaborn = Root.SeabornInterface()
seaborn.Initialize()

import CMSDIJET.StatisticalTools.limit_configuration as limit_config
import CMSDIJET.StatisticalTools.trigger_efficiency as trigger_efficiency
import MyTools.RootUtils.root_plots as root_plots
dijet_binning = array("d", [0, 1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5000]) #5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8000])

srs = ["lowmass", "highmass"]
notrig_analyses = {"lowmass":"NoTrigger_eta1p7_CSVTM", "highmass":"NoTrigger_eta2p2_CSVTM"}
trig_analyses = {"lowmass":"trigbbl_CSVTM", "highmass":"trigbbh_CSVTM"}

# QCD samples and stuff
qcd_samples = ["QCD_Pt-80to120_TuneZ2star_8TeV_pythia6",
"QCD_Pt-120to170_TuneZ2star_8TeV_pythia6",
"QCD_Pt-170to300_TuneZ2star_8TeV_pythia6",
"QCD_Pt-300to470_TuneZ2star_8TeV_pythia6",
"QCD_Pt-470to600_TuneZ2star_8TeV_pythia6",
"QCD_Pt-600to800_TuneZ2star_8TeV_pythia6",
"QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6",
"QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6",
"QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6",
"QCD_Pt-1800_TuneZ2star_8TeV_pythia6"]
qcd_cross_sections = {}
qcd_cross_sections["QCD_Pt-80to120_TuneZ2star_8TeV_pythia6"] = 1033680
qcd_cross_sections["QCD_Pt-120to170_TuneZ2star_8TeV_pythia6"] = 156293.3
qcd_cross_sections["QCD_Pt-170to300_TuneZ2star_8TeV_pythia6"] = 34138.15
qcd_cross_sections["QCD_Pt-300to470_TuneZ2star_8TeV_pythia6"] = 1759.549
qcd_cross_sections["QCD_Pt-470to600_TuneZ2star_8TeV_pythia6"] = 113.8791
qcd_cross_sections["QCD_Pt-600to800_TuneZ2star_8TeV_pythia6"] = 26.9921
qcd_cross_sections["QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6"] = 3.550036
qcd_cross_sections["QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6"] = 0.737844
qcd_cross_sections["QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6"] = 0.03352235
qcd_cross_sections["QCD_Pt-1800_TuneZ2star_8TeV_pythia6"] = 0.001829005
online_btag_eff = {
	"trigbbl_CSVTM":1.82719e-01,
	"trigbbh_CSVTM":4.89446e-01,
}
qcd_hists_num = {} # QCD with trigbb*
qcd_hists_den = {} # QCD with NoTrigger*
qcd_hists_eff = {}
for sr in srs:
	for qcd_sample in qcd_samples:
		f_num = TFile(analysis_config.get_b_histogram_filename(trig_analyses[sr], qcd_sample), "READ")
		f_den = TFile(analysis_config.get_b_histogram_filename(notrig_analyses[sr], qcd_sample), "READ")
		input_nevents_num = f_num.Get("BHistograms/h_input_nevents").Integral()
		this_num_hist = f_num.Get("BHistograms/h_pfjet_mjj")
		#this_num_hist.Scale(19700.*qcd_cross_sections[qcd_sample] / input_nevents_num)
		input_nevents_den = f_den.Get("BHistograms/h_input_nevents").Integral()
		this_den_hist = f_den.Get("BHistograms/h_pfjet_mjj")
		#this_den_hist.Scale(19700.*qcd_cross_sections[qcd_sample] / input_nevents_den)
		if not sr in qcd_hists_num:
			qcd_hists_num[sr] = this_num_hist.Clone()
			qcd_hists_num[sr].SetDirectory(0)
			qcd_hists_den[sr] = this_den_hist.Clone()
			qcd_hists_den[sr].SetDirectory(0)
		else:
			qcd_hists_num[sr].Add(this_num_hist)
			qcd_hists_den[sr].Add(this_den_hist)
		f_num.Close()
		f_den.Close()
	qcd_hists_num[sr] = histogram_tools.rebin_histogram(qcd_hists_num[sr], dijet_binning)
	qcd_hists_den[sr] = histogram_tools.rebin_histogram(qcd_hists_den[sr], dijet_binning)
	for bin in xrange(1, qcd_hists_num[sr].GetNbinsX()+1):
		bin_center = qcd_hists_num[sr].GetXaxis().GetBinCenter(bin)
		if sr == "lowmass":
			if bin_center < 296 or bin_center > 1058:
				qcd_hists_num[sr].SetBinContent(bin, 0)
				qcd_hists_num[sr].SetBinError(bin, 0)
				qcd_hists_den[sr].SetBinContent(bin, 0)
				qcd_hists_den[sr].SetBinError(bin, 0)
		elif sr == "highmass":
			if bin_center < 526 or bin_center > 1607:
				qcd_hists_num[sr].SetBinContent(bin, 0)
				qcd_hists_num[sr].SetBinError(bin, 0)
				qcd_hists_den[sr].SetBinContent(bin, 0)
				qcd_hists_den[sr].SetBinError(bin, 0)

	qcd_hists_eff[sr] = qcd_hists_num[sr].Clone()
	qcd_hists_eff[sr].Divide(qcd_hists_num[sr], qcd_hists_den[sr], 1, 1, "B")

# Signal samples and stuff
signals = ["Hbb", "ZPrime", "RSG"]
signal_masses = {"lowmass":[325,350,400,500,600,750], "highmass":[600,750,900,1200]}
signal_hists_num = {}
signal_hists_den = {}
signal_hists_eff = {}
for sr in srs:
	signal_hists_num[sr] = {}
	signal_hists_den[sr] = {}
	signal_hists_eff[sr] = {}
	for signal in signals:
		for mass in signal_masses[sr]:
			signal_tag = analysis_config.simulation.get_signal_tag(signal, mass, "FULLSIM")
			f_num = TFile(analysis_config.get_b_histogram_filename(trig_analyses[sr], signal_tag), "READ")
			f_den = TFile(analysis_config.get_b_histogram_filename(notrig_analyses[sr], signal_tag), "READ")
			this_num_hist = f_num.Get("BHistograms/h_pfjet_mjj")
			this_den_hist = f_den.Get("BHistograms/h_pfjet_mjj")
			if not signal in signal_hists_num[sr]:
				signal_hists_num[sr][signal] = this_num_hist.Clone()
				signal_hists_num[sr][signal].SetDirectory(0)
				signal_hists_den[sr][signal] = this_den_hist.Clone()
				signal_hists_den[sr][signal].SetDirectory(0)
			else:
				signal_hists_num[sr][signal].Add(this_num_hist)
				signal_hists_den[sr][signal].Add(this_den_hist)
			f_num.Close()
			f_den.Close()
		signal_hists_num[sr][signal] = histogram_tools.rebin_histogram(signal_hists_num[sr][signal], dijet_binning)
		signal_hists_den[sr][signal] = histogram_tools.rebin_histogram(signal_hists_den[sr][signal], dijet_binning)
		for bin in xrange(1, signal_hists_num[sr][signal].GetNbinsX()+1):
			bin_center = signal_hists_num[sr][signal].GetXaxis().GetBinCenter(bin)
			if sr == "lowmass":
				if bin_center < 296 or bin_center > 1058:
					signal_hists_num[sr][signal].SetBinContent(bin, 0)
					signal_hists_num[sr][signal].SetBinError(bin, 0)
					signal_hists_den[sr][signal].SetBinContent(bin, 0)
					signal_hists_den[sr][signal].SetBinError(bin, 0)
			elif sr == "highmass":
				if bin_center < 526 or bin_center > 1607:
					signal_hists_num[sr][signal].SetBinContent(bin, 0)
					signal_hists_num[sr][signal].SetBinError(bin, 0)
					signal_hists_den[sr][signal].SetBinContent(bin, 0)
					signal_hists_den[sr][signal].SetBinError(bin, 0)
		signal_hists_eff[sr][signal] = signal_hists_num[sr][signal].Clone()
		signal_hists_eff[sr][signal].Divide(signal_hists_num[sr][signal], signal_hists_den[sr][signal], 1, 1, "B")

# Data stuff
jetht_hists_num = {}
jetht_hists_den = {}
jetht_hists_eff = {}
ht_slices = {
	"lowmass":{
		#"trigjetht200_eta1p7_CSVTM":[220, 386],
		"trigjetht250_eta1p7_CSVTM":[386, 489],
		"trigjetht300_eta1p7_CSVTM":[489, 526],
		"trigjetht350_eta1p7_CSVTM":[526, 606],
		"trigjetht400_eta1p7_CSVTM":[606, 649],
		"trigjetht450_eta1p7_CSVTM":[649, 740],
		"trigjetht500_eta1p7_CSVTM":[740, 788],
		"trigjetht550_eta1p7_CSVTM":[788, 890],
		#"HT650":[800, 890],
		"trigjetht_eta1p7_CSVTM":[890, 2000]
	}, "highmass":{
		#"trigjetht200_CSVTM":[220, 386],
		"trigjetht250_CSVTM":[386, 489],
		"trigjetht300_CSVTM":[489, 526],
		"trigjetht350_CSVTM":[526, 606],
		"trigjetht400_CSVTM":[606, 649],
		"trigjetht450_CSVTM":[649, 740],
		"trigjetht500_CSVTM":[740, 788],
		"trigjetht550_CSVTM":[788, 890],
		#"HT650":[800, 890],
		"trigjetht_CSVTM":[890, 2000]
	}
}
for sr in srs:
	f_num = TFile(analysis_config.get_b_histogram_filename(trig_analyses[sr], "BJetPlusX_2012BCD"), "READ")
	jetht_hists_num[sr] = f_num.Get("BHistograms/h_pfjet_mjj")
	jetht_hists_num[sr].SetDirectory(0)
	for analysis, slice_range in ht_slices[sr].iteritems():
		f = TFile(analysis_config.get_b_histogram_filename(analysis, "JetHT_2012BCD"), "READ")
		this_hist = f.Get("BHistograms/h_pfjet_mjj")
		for bin in xrange(1, this_hist.GetNbinsX()+1):
			bin_center = this_hist.GetXaxis().GetBinCenter(bin)
			if bin_center < slice_range[0] or bin_center > slice_range[1]:
				this_hist.SetBinContent(bin, 0)
				this_hist.SetBinError(bin, 0)
		if not sr in jetht_hists_den:
			jetht_hists_den[sr] = this_hist.Clone()
			jetht_hists_den[sr].SetDirectory(0)
		else:
			jetht_hists_den[sr].Add(this_hist)
		f.Close()
	jetht_hists_num[sr] = histogram_tools.rebin_histogram(jetht_hists_num[sr], dijet_binning)
	jetht_hists_den[sr] = histogram_tools.rebin_histogram(jetht_hists_den[sr], dijet_binning)
	for bin in xrange(1, jetht_hists_num[sr].GetNbinsX()+1):
		bin_center = jetht_hists_num[sr].GetXaxis().GetBinCenter(bin)
		if sr == "lowmass":
			if bin_center < 296 or bin_center > 1058:
				jetht_hists_num[sr].SetBinContent(bin, 0)
				jetht_hists_num[sr].SetBinError(bin, 0)
				jetht_hists_den[sr].SetBinContent(bin, 0)
				jetht_hists_den[sr].SetBinError(bin, 0)
		elif sr == "highmass":
			if bin_center < 526 or bin_center > 1607:
				jetht_hists_num[sr].SetBinContent(bin, 0)
				jetht_hists_num[sr].SetBinError(bin, 0)
				jetht_hists_den[sr].SetBinContent(bin, 0)
				jetht_hists_den[sr].SetBinError(bin, 0)

	jetht_hists_eff[sr] = jetht_hists_num[sr].Clone()
	jetht_hists_eff[sr].Reset()
	# Efficiency computation is a little complication because of the denominator prescale
	for bin in xrange(1, jetht_hists_num[sr].GetNbinsX()+1):
		num = jetht_hists_num[sr].GetBinContent(bin)
		num_err = jetht_hists_num[sr].GetBinError(bin)
		den = jetht_hists_den[sr].GetBinContent(bin)
		den_err = jetht_hists_den[sr].GetBinError(bin)
		if den > 0.:
			#print "Bin {} eff = {}/{}={}".format(bin, num, den, num/den)
			eff = num / den
			den_prescale = den_err**2 / den
			eff_err = num / den * sqrt((den + num * (den_prescale - 2)) / (num * den))
			jetht_hists_eff[sr].SetBinContent(bin, eff)
			jetht_hists_eff[sr].SetBinError(bin, eff_err)
		else:
			#print "Bin {} den={}, num={}".format(bin, den, num)
			jetht_hists_eff[sr].SetBinContent(bin, 0.)
			jetht_hists_eff[sr].SetBinError(bin, 0.)

# SingleMu
singlemu_num_analyses = {
	"lowmass":"trigmu24ibbl_lowmass_CSVTM",
	"highmass":"trigmu24ibbh_highmass_CSVTM",
}
singlemu_den_analyses = {
	"lowmass":"trigmu24i_lowmass_CSVTM",
	"highmass":"trigmu24i_highmass_CSVTM",
}
singlemu_hists_num = {}
singlemu_hists_den = {}
singlemu_hists_eff = {}
for sr in srs:
	f_num = TFile(analysis_config.get_b_histogram_filename(singlemu_num_analyses[sr], "SingleMu_2012"), "READ")
	singlemu_hists_num[sr] = f_num.Get("BHistograms/h_pfjet_mjj")
	singlemu_hists_num[sr].SetDirectory(0)
	f_den = TFile(analysis_config.get_b_histogram_filename(singlemu_den_analyses[sr], "SingleMu_2012"), "READ")
	singlemu_hists_den[sr] = f_den.Get("BHistograms/h_pfjet_mjj")
	singlemu_hists_den[sr].SetDirectory(0)
	singlemu_hists_num[sr] = histogram_tools.rebin_histogram(singlemu_hists_num[sr], dijet_binning)
	singlemu_hists_den[sr] = histogram_tools.rebin_histogram(singlemu_hists_den[sr], dijet_binning)
	for bin in xrange(1, qcd_hists_num[sr].GetNbinsX()+1):
		bin_center = singlemu_hists_num[sr].GetXaxis().GetBinCenter(bin)
		if sr == "lowmass":
			if bin_center < 296 or bin_center > 1058:
				singlemu_hists_num[sr].SetBinContent(bin, 0)
				singlemu_hists_num[sr].SetBinError(bin, 0)
				singlemu_hists_den[sr].SetBinContent(bin, 0)
				singlemu_hists_den[sr].SetBinError(bin, 0)
		elif sr == "highmass":
			if bin_center < 526 or bin_center > 1607:
				singlemu_hists_num[sr].SetBinContent(bin, 0)
				singlemu_hists_num[sr].SetBinError(bin, 0)
				singlemu_hists_den[sr].SetBinContent(bin, 0)
				singlemu_hists_den[sr].SetBinError(bin, 0)
	singlemu_hists_eff[sr] = singlemu_hists_num[sr].Clone()
	singlemu_hists_eff[sr].SetDirectory(0)
	singlemu_hists_eff[sr].Divide(singlemu_hists_num[sr], singlemu_hists_den[sr], 1., 1., "B")
	f_num.Close()
	f_den.Close()


# Plot
for sr in srs:
	c = TCanvas("c_online_btag_eff_{}".format(sr), "c_online_btag_eff_{}".format(sr), 800, 1000)
	top = TPad("top", "top", 0., 0.5, 1., 1.)
	top.SetBottomMargin(0.03)
	top.Draw()
	top.cd()
	jetht_hists_eff[sr].SetMarkerStyle(20)
	jetht_hists_eff[sr].GetXaxis().SetTitleSize(0)
	jetht_hists_eff[sr].GetXaxis().SetLabelSize(0)
	jetht_hists_eff[sr].GetYaxis().SetTitle("Online b-tag efficiency")
	jetht_hists_eff[sr].SetMinimum(-0.1)
	if sr == "lowmass":
		jetht_hists_eff[sr].SetMaximum(0.5)
	else:
		jetht_hists_eff[sr].SetMaximum(1.0)
	if sr == "lowmass":
		jetht_hists_eff[sr].GetXaxis().SetRangeUser(200., 1500.)
	else:
		jetht_hists_eff[sr].GetXaxis().SetRangeUser(400., 2000.)
	jetht_hists_eff[sr].Draw("e1")
	signal_hists_eff[sr]["Hbb"].SetMarkerStyle(20)
	signal_hists_eff[sr]["Hbb"].SetMarkerSize(0)
	signal_hists_eff[sr]["Hbb"].SetLineStyle(1)
	signal_hists_eff[sr]["Hbb"].SetLineWidth(1)
	signal_hists_eff[sr]["Hbb"].SetLineColor(seaborn.GetColorRoot("default", 1))
	signal_hists_eff[sr]["Hbb"].Draw("hist same")

	singlemu_hists_eff[sr].SetMarkerStyle(24)
	singlemu_hists_eff[sr].SetMarkerColor(seaborn.GetColorRoot("default", 0))
	singlemu_hists_eff[sr].SetLineColor(seaborn.GetColorRoot("default", 0))
	singlemu_hists_eff[sr].Draw("e1 same")
	qcd_hists_eff[sr].SetMarkerStyle(25)
	qcd_hists_eff[sr].SetMarkerColor(seaborn.GetColorRoot("default", 2))
	qcd_hists_eff[sr].SetLineColor(seaborn.GetColorRoot("default", 2))
	qcd_hists_eff[sr].Draw("e1 same")
	jetht_hists_eff[sr].Draw("e1 same")
	l = TLegend(0.75, 0.3, 0.9, 0.7)
	l.SetFillColor(0)
	l.SetBorderSize(0)
	l.AddEntry(jetht_hists_eff[sr], "JetHT data", "pl")
	l.AddEntry(singlemu_hists_eff[sr], "mu24i data", "pl")
	l.AddEntry(qcd_hists_eff[sr], "QCD MC", "pl")
	l.AddEntry(signal_hists_eff[sr]["Hbb"], "Signal MC", "l")
	l.Draw()
	c.cd()
	bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
	bottom.SetTopMargin(0.02)
	bottom.Draw()
	bottom.cd()
	bottom.SetGrid()
	if sr == "lowmass":
		frame_ratio = TH1F("frame_ratio", "frame_ratio", 100, 200., 1500.)
	else:
		frame_ratio = TH1F("frame_ratio", "frame_ratio", 100, 400., 2000.)

	frame_ratio.SetMinimum(0.5)
	frame_ratio.SetMaximum(1.5)
	frame_ratio.GetXaxis().SetTitle("m_{jj} [GeV]")
	frame_ratio.GetYaxis().SetTitle("#frac{#epsilon_{data}}{#epsilon_{MC}}")
	frame_ratio.Draw()
	ratio_hist_jetht = jetht_hists_eff[sr].Clone()
	ratio_hist_jetht.Reset()
	for bin in xrange(1, ratio_hist_jetht.GetNbinsX()+1):
		num = jetht_hists_eff[sr].GetBinContent(bin)
		num_err = jetht_hists_eff[sr].GetBinError(bin)
		den = qcd_hists_eff[sr].GetBinContent(bin)
		den_err = qcd_hists_eff[sr].GetBinError(bin)
		if num > 0 and den > 0:
			r = num / den
			r_err = num/den * sqrt((num_err/num)**2 + (den_err/den)**2)
		else:
			r = 0.
			r_err = 0.
		ratio_hist_jetht.SetBinContent(bin, r)
		ratio_hist_jetht.SetBinError(bin, r_err)
	ratio_hist_jetht.Draw("e1 same")

	ratio_hist_singlemu = singlemu_hists_eff[sr].Clone()
	ratio_hist_singlemu.Reset()
	for bin in xrange(1, ratio_hist_singlemu.GetNbinsX()+1):
		num = singlemu_hists_eff[sr].GetBinContent(bin)
		num_err = singlemu_hists_eff[sr].GetBinError(bin)
		den = qcd_hists_eff[sr].GetBinContent(bin)
		den_err = qcd_hists_eff[sr].GetBinError(bin)
		if num > 0 and den > 0:
			r = num / den
			r_err = num/den * sqrt((num_err/num)**2 + (den_err/den)**2)
		else:
			r = 0.
			r_err = 0.
		ratio_hist_singlemu.SetBinContent(bin, r)
		ratio_hist_singlemu.SetBinError(bin, r_err)
	ratio_hist_singlemu.Draw("e1 same")
	ratio_hist_jetht.Draw("e1 same")
	# Average from fit
	#avg_eff_ratio = (qcd_hists_num[sr].Integral() / qcd_hists_den[sr].Integral()) / (jetht_hists_num[sr].Integral() / jetht_hists_den[sr].Integral())
	#if sr == "lowmass":
	#	avg_eff_ratio_line = TLine(296, avg_eff_ratio, 1058, avg_eff_ratio)
	#else:
	#	avg_eff_ratio_line = TLine(526, avg_eff_ratio, 1607, avg_eff_ratio)
	#avg_eff_ratio_line.SetLineStyle(2)
	#avg_eff_ratio_line.SetLineColor(seaborn.GetColorRoot("pastel", 2))
	#avg_eff_ratio_line.Draw("same")
	ratio_hist_combination = ratio_hist_jetht.Clone()
	for bin in xrange(1, ratio_hist_combination.GetNbinsX()+1):
		bin_center = ratio_hist_combination.GetXaxis().GetBinCenter(bin)
		if bin_center < 386:
			ratio_hist_combination.SetBinContent(bin, ratio_hist_singlemu.GetBinContent(bin))
			ratio_hist_combination.SetBinError(bin, ratio_hist_singlemu.GetBinError(bin))
	if sr == "lowmass":
		constant_fit = TF1("effratio_{}".format(sr), "[0]", 296, 1058)
		constant_fit.SetParameter(0, 1.0)
	else:
		constant_fit = TF1("effratio_{}".format(sr), "[0]", 526, 1607)
		constant_fit.SetParameter(0, 1.0)
	ratio_hist_combination.Fit(constant_fit, "R0")

	constant_fit.SetLineStyle(2)
	constant_fit.SetLineColor(seaborn.GetColorRoot("muted", 2))
	constant_fit.Draw("same")
	c.cd()
	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")
	ROOT.SetOwnership(c, False)
	ROOT.SetOwnership(top, False)
	ROOT.SetOwnership(bottom, False)
	print "{} avg eff ratio = {} +/- {}" .format(sr, constant_fit.GetParameter(0), constant_fit.GetParError(0))

	# Fit efficiencies with constants
	# Data
	eff_hist_combination = jetht_hists_eff[sr].Clone()
	for bin in xrange(1, eff_hist_combination.GetNbinsX()+1):
		bin_center = eff_hist_combination.GetXaxis().GetBinCenter(bin)
		if bin_center < 386:
			eff_hist_combination.SetBinContent(bin, singlemu_hists_eff[sr].GetBinContent(bin))
			eff_hist_combination.SetBinError(bin, singlemu_hists_eff[sr].GetBinError(bin))
	if sr == "lowmass":
		data_eff_fit = TF1("data_eff_{}".format(sr), "[0]", 296, 1058)
		data_eff_fit.SetParameter(0, 0.2)
	else:
		data_eff_fit = TF1("data_eff_{}".format(sr), "[0]", 526, 1607)
		data_eff_fit.SetParameter(0, 0.5)
	eff_hist_combination.Fit(data_eff_fit, "R0")
	print "{} data efficiency = {} +/- {}".format(sr, data_eff_fit.GetParameter(0), data_eff_fit.GetParError(0))

	# QCD MC
	if sr == "lowmass":
		qcd_eff_fit = TF1("qcd_eff_{}".format(sr), "[0]", 296, 1058)
		qcd_eff_fit.SetParameter(0, 0.2)
	else:
		qcd_eff_fit = TF1("qcd_eff_{}".format(sr), "[0]", 526, 1607)
		qcd_eff_fit.SetParameter(0, 0.5)
	qcd_hists_eff[sr].Fit(qcd_eff_fit, "R0")
	print "{} QCD MC efficiency = {} +/- {}".format(sr, qcd_eff_fit.GetParameter(0), qcd_eff_fit.GetParError(0))
	print "{} eff ratio data/MC = {}".format(sr, data_eff_fit.GetParameter(0) / qcd_eff_fit.GetParameter(0))


