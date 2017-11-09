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

analyses = ["NoTrigger_eta1p7_CSVTM", "NoTrigger_eta2p2_CSVTM"]
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

truth_categories = ["bb", "bc", "bg", "bl", "cc", "cg", "cl", "gg", "gl", "ll"]
truth_categories_incl = ["inclusive", "bb", "bc", "bg", "bl", "cc", "cg", "cl", "gg", "gl", "ll"]



for analysis in analyses:
	h_mjj = None
	h_mjj_truthbb = None
	for qcd_sample in qcd_samples:
		f = TFile(analysis_config.get_b_histogram_filename(analysis, qcd_sample), "READ")
		input_nevents = f.Get("BHistograms/h_input_nevents").Integral()
		this_h_mjj = f.Get("BHistograms/h_pfjet_mjj").Clone()
		this_h_mjj.Scale(19700.*qcd_cross_sections[qcd_sample] / input_nevents)
		this_h_mjj_truthbb = f.Get("BHistograms/h_pfjet_mjj_truthbb").Clone()
		this_h_mjj_truthbb.Scale(19700.*qcd_cross_sections[qcd_sample] / input_nevents)
		if not h_mjj:
			h_mjj = this_h_mjj.Clone()
			h_mjj.SetDirectory(0)
			h_mjj_truthbb = this_h_mjj_truthbb.Clone()
			h_mjj_truthbb.SetDirectory(0)
		else:
			h_mjj.Add(this_h_mjj)
			h_mjj_truthbb.Add(this_h_mjj_truthbb)
		f.Close()
	h_mjj = histogram_tools.rebin_histogram(h_mjj, dijet_binning, normalization_bin_width=1.)
	h_mjj_truthbb = histogram_tools.rebin_histogram(h_mjj_truthbb, dijet_binning, normalization_bin_width=1.)
	c = TCanvas("c_purity_{}".format(analysis), "c_purity_{}".format(analysis), 800, 1000)
	top = TPad("top", "top", 0., 0.5, 1., 1.)
	top.SetBottomMargin(0.02)
	top.SetLogy()
	top.Draw()
	top.cd()
	h_mjj.SetLineColor(1)
	h_mjj.SetLineWidth(2)
	h_mjj.GetXaxis().SetTitleSize(0)
	h_mjj.GetXaxis().SetLabelSize(0)
	h_mjj.GetYaxis().SetTitle("Events / GeV")
	h_mjj.GetXaxis().SetRangeUser(0., 2000.)
	h_mjj.Draw("hist")
	h_mjj_truthbb.SetLineColor(seaborn.GetColorRoot("cubehelixlarge", 2, 20))
	h_mjj_truthbb.SetLineWidth(1)
	h_mjj_truthbb.Draw("hist same")
	c.cd()
	bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
	bottom.SetTopMargin(0.02)
	bottom.SetBottomMargin(0.2)
	bottom.Draw()
	bottom.cd()
	h_ratio = h_mjj_truthbb.Clone()
	h_ratio.Divide(h_mjj_truthbb, h_mjj, 1., 1., "B")
	h_ratio.SetMarkerStyle(20)
	h_ratio.SetMarkerSize(1)
	h_ratio.SetMarkerColor(1)
	h_ratio.SetLineWidth(1)
	h_ratio.SetLineColor(1)
	h_ratio.GetYaxis().SetTitle("bb purity")
	h_ratio.GetXaxis().SetTitle("m_{jj} [GeV]")
	h_ratio.GetXaxis().SetRangeUser(0., 2000.)
	h_ratio.Draw()
	c.cd()
	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")


# Purity plot for truth-categorized events
for analysis in analyses:
	histograms = {}
	histograms["inclusive"] = None
	for truth_category in truth_categories:
		histograms[truth_category] = None
	for qcd_sample in qcd_samples:
		print "Opening " + analysis_config.get_b_histogram_filename(analysis, qcd_sample)
		f = TFile(analysis_config.get_b_histogram_filename(analysis, qcd_sample), "READ")
		if not f.IsOpen():
			print "ERROR : Unable to open file at " + analysis_config.get_b_histogram_filename(analysis, qcd_sample)
			sys.exit(1)
		input_nevents = f.Get("BHistograms/h_input_nevents").Integral()
		this_histograms = {}
		this_histograms["inclusive"] = f.Get("BHistograms/h_pfjet_mjj").Clone()
		this_histograms["inclusive"].Scale(19700.*qcd_cross_sections[qcd_sample] / input_nevents)
		if not histograms["inclusive"]:
			histograms["inclusive"] = this_histograms["inclusive"].Clone()
			histograms["inclusive"].SetDirectory(0)
		else:
			histograms["inclusive"].Add(this_histograms["inclusive"])
		for truth_category in truth_categories:
			this_histograms[truth_category] = f.Get("BHistograms/h_pfjet_mjj_truth{}".format(truth_category)).Clone()
			this_histograms[truth_category].Scale(19700.*qcd_cross_sections[qcd_sample] / input_nevents)
			if not histograms[truth_category]:
				histograms[truth_category] = this_histograms[truth_category].Clone()
				histograms[truth_category].SetDirectory(0)
			else:
				histograms[truth_category].Add(this_histograms[truth_category])
		f.Close()
	for cat, hist in histograms.iteritems():
		histograms[cat] = histogram_tools.rebin_histogram(hist, dijet_binning, normalization_bin_width=1.)

	c = TCanvas("c_purity_truthbreakdown_{}".format(analysis), "c_purity_truthbreakdown_{}", 800, 1000)
	l = TLegend(0.65, 0.42, 0.9, 0.85)
	l.SetFillColor(0)
	l.SetBorderSize(0)
	top = TPad("top", "top", 0., 0.5, 1., 1.)
	top.SetLogy()
	top.SetBottomMargin(0.02)
	top.Draw()
	top.cd()
	stack = THStack("truthstack", "truthstack")
	for style_counter, truth_category in enumerate(truth_categories):
		histograms[truth_category].SetFillColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, 20))
		histograms[truth_category].SetFillStyle(1001)
		stack.Add(histograms[truth_category])
	frame_top = TH1D("frame_top", "frame_top", 100, 0., 2000.)
	frame_top.SetMinimum(0.1)
	frame_top.SetMaximum(histograms["inclusive"].GetMaximum() * 10.)
	frame_top.GetXaxis().SetLabelSize(0)
	frame_top.GetXaxis().SetTitleSize(0)
	frame_top.GetYaxis().SetTitle("Events / GeV")
	frame_top.Draw()
	stack.Draw("hist same")
	histograms["inclusive"].SetMarkerStyle(20)
	histograms["inclusive"].Draw("p same")
	l.AddEntry(histograms["inclusive"], "Inclusive", "p")
	for style_counter, truth_category in enumerate(reversed(truth_categories)):
		l.AddEntry(histograms[truth_category], truth_category, "f")
	l.Draw()

	c.cd()
	bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
	bottom.SetTopMargin(0.02)
	bottom.SetBottomMargin(0.2)
	bottom.Draw()
	bottom.cd()
	frame_bottom = TH1D("frame_bottom", "frame_bottom", 100, 0., 2000.)
	frame_bottom.SetMinimum(-0.2)
	frame_bottom.SetMaximum(1.2)
	frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
	frame_bottom.GetYaxis().SetTitle("Fraction of Events")
	frame_bottom.Draw()
	purity_hists = {}
	purity_stack = THStack("puritystack", "puritystack")
	for style_counter, truth_category in enumerate(truth_categories):
		purity_hists[truth_category] = histograms[truth_category].Clone()
		purity_hists[truth_category].Divide(histograms[truth_category], histograms["inclusive"], 1, 1, "B")
		purity_hists[truth_category].SetFillColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, 20))
		purity_hists[truth_category].SetFillStyle(1001)
		purity_hists[truth_category].SetLineWidth(0)
		purity_stack.Add(purity_hists[truth_category])
		#purity_hists[truth_category].Draw("p same")
	purity_stack.Draw("hist same")
	c.cd()
	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")

# Online b-tag efficiency for truth categories
for num_an, den_an in [("trigbbl_CSVTM", "NoTrigger_eta1p7_CSVTM"), ("trigbbh_CSVTM", "NoTrigger_eta2p2_CSVTM")]:
	num_hists = {}
	den_hists = {}
	eff_hists = {}
	for truth_category in truth_categories_incl:
		for qcd_sample in qcd_samples:
			f = TFile(analysis_config.get_b_histogram_filename(num_an, qcd_sample), "READ")
			if not f.IsOpen():
				print "ERROR : Unable to open file at " + analysis_config.get_b_histogram_filename(num_an, qcd_sample)
				sys.exit(1)
			input_nevents = f.Get("BHistograms/h_input_nevents").Integral()
			if truth_category == "inclusive":
				this_histogram = f.Get("BHistograms/h_pfjet_mjj")
			else:
				this_histogram = f.Get("BHistograms/h_pfjet_mjj_truth{}".format(truth_category))
			this_histogram = histogram_tools.rebin_histogram(this_histogram, dijet_binning, normalization_bin_width=1.)
			#this_histogram.Scale(19700.*qcd_cross_sections[qcd_sample] / input_nevents)
			if truth_category in num_hists:
				num_hists[truth_category].Add(this_histogram)
			else:
				num_hists[truth_category] = this_histogram.Clone()
				num_hists[truth_category].SetDirectory(0)
			f.Close()
		for qcd_sample in qcd_samples:
			f = TFile(analysis_config.get_b_histogram_filename(den_an, qcd_sample), "READ")
			if not f.IsOpen():
				print "ERROR : Unable to open file at " + analysis_config.get_b_histogram_filename(den_an, qcd_sample)
				sys.exit(1)
			input_nevents = f.Get("BHistograms/h_input_nevents").Integral()
			if truth_category == "inclusive":
				this_histogram = f.Get("BHistograms/h_pfjet_mjj")
			else:
				this_histogram = f.Get("BHistograms/h_pfjet_mjj_truth{}".format(truth_category))
			this_histogram = histogram_tools.rebin_histogram(this_histogram, dijet_binning, normalization_bin_width=1.)
			#this_histogram.Scale(19700.*qcd_cross_sections[qcd_sample] / input_nevents)
			if truth_category in den_hists:
				den_hists[truth_category].Add(this_histogram)
			else:
				den_hists[truth_category] = this_histogram.Clone()
				den_hists[truth_category].SetDirectory(0)
			f.Close()

		eff_hists[truth_category] = num_hists[truth_category].Clone()
		eff_hists[truth_category].Divide(num_hists[truth_category], den_hists[truth_category], 1., 1., "B")

	if "bbl" in num_an:
		cname = "c_online_btag_eff_truthbreakdown_lowmass"
	else:
		cname = "c_online_btag_eff_truthbreakdown_highmass"
	c = TCanvas(cname, cname, 700, 500)
	c.SetRightMargin(0.2)
	l = TLegend(0.82, 0.3, 0.99, 0.8)
	l.SetFillColor(0)
	l.SetBorderSize(0)
	if "bbl" in num_an:
		frame = TH1F("frame", "frame", 100, 200., 1300.)
	else:
		frame = TH1F("frame", "frame", 100, 200., 1900.)
	frame.SetMinimum(0.)
	frame.SetMaximum(0.7)
	frame.GetXaxis().SetTitle("m_{jj} [GeV]")
	frame.GetYaxis().SetTitle("Online b-tag efficiency")
	frame.Draw()
	for style_counter, truth_category in enumerate(truth_categories_incl):
		if truth_category == "inclusive":
			eff_hists[truth_category].SetMarkerStyle(20)
			eff_hists[truth_category].Draw("p same")
			l.AddEntry(eff_hists[truth_category], "Inclusive", "p")
		else:
			eff_hists[truth_category].SetMarkerStyle(21+style_counter - 1)
			eff_hists[truth_category].SetMarkerColor(seaborn.GetColorRoot("cubehelixlarge", style_counter-1, 20))
			eff_hists[truth_category].SetLineColor(seaborn.GetColorRoot("cubehelixlarge", style_counter-1, 20))
			eff_hists[truth_category].Draw("p same")
			l.AddEntry(eff_hists[truth_category], truth_category, "p")
	l.Draw()
	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")
