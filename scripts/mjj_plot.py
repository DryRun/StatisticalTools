import os
import sys
import ROOT
from array import array
from ROOT import *

import CMSDIJET.StatisticalTools.limit_configuration as limit_config

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

import CMSDIJET.StatisticalTools.trigger_efficiency as trigger_efficiency


def BackgroundFit_f1(x, par):
	return par[0] * (1. - (x[0] / 8.e3))**par[1] / ((x[0] / 8.e3)**(par[2] + par[3] * TMath.Log((x[0] / 8.e3))))

def BackgroundFit_f2(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3))**par[2]

def BackgroundFit_f3(x, par):
	return par[0] / (1 + par[1] * (x[0] / 8.e3))**par[2]

def BackgroundFit_f4(x, par):
	if 1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2 <= 0:
		return 0
	elif ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2)**par[3]) < 1.e-15:
		return 0
	else:
		return par[0] / ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2)**par[3])

def BackgroundFit_f5(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3)**(1./3.))**par[2]

def BackgroundFit_f6(x, par):
	if 1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2 + par[3] * (x[0]/8.e3)**3 <= 0:
		return 0
	elif ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2  + par[3] * (x[0]/8.e3)**3)**par[4]) < 1.e-15:
		return 0
	else:
		return par[0] / ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2  + par[3] * (x[0]/8.e3)**3)**par[4])

def BackgroundFit_f1_trigcorr_bbl(x, par):
	return par[0] * (1. - (x[0] / 8.e3))**par[1] / ((x[0] / 8.e3)**(par[2] + par[3] * TMath.Log((x[0] / 8.e3)))) * trigger_efficiency.trigger_efficiency_bbl(x[0])

def BackgroundFit_f2_trigcorr_bbl(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3))**par[2] * trigger_efficiency.trigger_efficiency_bbl(x[0])

def BackgroundFit_f3_trigcorr_bbl(x, par):
	return par[0] / (1 + par[1] * (x[0] / 8.e3))**par[2] * trigger_efficiency.trigger_efficiency_bbl(x[0])

def BackgroundFit_f4_trigcorr_bbl(x, par):
	if 1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2 <= 0:
		return 0
	else:
		return par[0] / ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2)**par[3]) * trigger_efficiency.trigger_efficiency_bbl(x[0])

def BackgroundFit_f5_trigcorr_bbl(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3)**(1./3.))**par[2] * trigger_efficiency.trigger_efficiency_bbl(x[0])

def BackgroundFit_f6_trigcorr_bbl(x, par):
	if 1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2 + par[3] * (x[0]/8.e3)**3 <= 0:
		return 0
	elif ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2  + par[3] * (x[0]/8.e3)**3)**par[4]) < 1.e-15:
		return 0
	else:
		return (par[0] / ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2  + par[3] * (x[0]/8.e3)**3)**par[4]))  * trigger_efficiency.trigger_efficiency_bbl(x[0])

def BackgroundFit_f1_trigcorr_bbh(x, par):
	return par[0] * (1. - (x[0] / 8.e3))**par[1] / ((x[0] / 8.e3)**(par[2] + par[3] * TMath.Log((x[0] / 8.e3)))) * trigger_efficiency.trigger_efficiency_bbh(x[0])

def BackgroundFit_f2_trigcorr_bbh(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3))**par[2] * trigger_efficiency.trigger_efficiency_bbh(x[0])

def BackgroundFit_f3_trigcorr_bbh(x, par):
	return par[0] / (1 + par[1] * (x[0] / 8.e3))**par[2] * trigger_efficiency.trigger_efficiency_bbh(x[0])

def BackgroundFit_f4_trigcorr_bbh(x, par):
	if (1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2) <= 0.:
		return 0
	elif ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2)**par[3]) < 1.e-15:
		return 0
	else:
		return par[0] / ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2)**par[3]) * trigger_efficiency.trigger_efficiency_bbh(x[0])

def BackgroundFit_f5_trigcorr_bbh(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3)**(1./3.))**par[2] * trigger_efficiency.trigger_efficiency_bbh(x[0])

def BackgroundFit_f6_trigcorr_bbh(x, par):
	if 1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2 + par[3] * (x[0]/8.e3)**3 <= 0:
		return 0
	elif ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2  + par[3] * (x[0]/8.e3)**3)**par[4]) < 1.e-15:
		return 0
	else:
		return (par[0] / ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2  + par[3] * (x[0]/8.e3)**3)**par[4]))  * trigger_efficiency.trigger_efficiency_bbh(x[0])


def make_background_tf1_from_roofitresult(fit_function, roofitresult, mjj_range=None, trigger_correction=None):
	print "[make_background_tf1_from_roofitresult] INFO : make_background_tf1_from_roofitresult for fit function " + fit_function
	print mjj_range
	roofitresult.Print()
	function = make_background_tf1(fit_function, mjj_range, trigger_correction=trigger_correction)
	function.SetParameter(0, roofitresult.floatParsFinal().find("background_" + fit_function + "_norm").getVal())
	function.SetParError(0, roofitresult.floatParsFinal().find("background_" + fit_function + "_norm").getError())

	for i in xrange(1, function.GetNpar()):
		function.SetParameter(i, roofitresult.floatParsFinal().find(fit_function + "_p" + str(i)).getVal())
		function.SetParError(i, roofitresult.floatParsFinal().find(fit_function + "_p" + str(i)).getError())

	return function

def make_background_tf1(fit_function, mjj_range, trigger_correction=None):
	mjj_min = mjj_range[0]
	mjj_max = mjj_range[1]

	if not trigger_correction:
		print "Making trigger correction-less fit function"
		if fit_function == "f1":
			background_tf1 = TF1('background_tf1_f1',BackgroundFit_f1, mjj_min, mjj_max, 4)
		elif fit_function == "f2":
			background_tf1 = TF1('background_tf1_f2',BackgroundFit_f2, mjj_min, mjj_max, 3)
		elif fit_function == "f3":
			background_tf1 = TF1('background_tf1_f3',BackgroundFit_f3, mjj_min, mjj_max, 3)
		elif fit_function == "f4":
			background_tf1 = TF1('background_tf1_f4',BackgroundFit_f4, mjj_min, mjj_max, 4)
		elif fit_function == "f5":
			background_tf1 = TF1('background_tf1_f5',BackgroundFit_f5, mjj_min, mjj_max, 3)
		elif fit_function == "f6":
			background_tf1 = TF1('background_tf1_f6',BackgroundFit_f6, mjj_min, mjj_max, 5)
		else:
			print "[make_background_tf1] ERROR : Unrecognized fit function " + fit_function
			sys.exit(1)
	elif trigger_correction == "bbl":
		if fit_function == "f1":
			background_tf1 = TF1('background_tf1_f1',BackgroundFit_f1_trigcorr_bbl, mjj_min, mjj_max, 4)
		elif fit_function == "f2":
			background_tf1 = TF1('background_tf1_f2',BackgroundFit_f2_trigcorr_bbl, mjj_min, mjj_max, 3)
		elif fit_function == "f3":
			background_tf1 = TF1('background_tf1_f3',BackgroundFit_f3_trigcorr_bbl, mjj_min, mjj_max, 3)
		elif fit_function == "f4":
			background_tf1 = TF1('background_tf1_f4',BackgroundFit_f4_trigcorr_bbl, mjj_min, mjj_max, 4)
		elif fit_function == "f5":
			background_tf1 = TF1('background_tf1_f5',BackgroundFit_f5_trigcorr_bbl, mjj_min, mjj_max, 3)
		elif fit_function == "f6":
			background_tf1 = TF1('background_tf1_f6',BackgroundFit_f6_trigcorr_bbl, mjj_min, mjj_max, 5)
		else:
			print "[make_background_tf1] ERROR : Unrecognized fit function " + fit_function
			sys.exit(1)
	elif trigger_correction == "bbh":
		if fit_function == "f1":
			background_tf1 = TF1('background_tf1_f1',BackgroundFit_f1_trigcorr_bbh, mjj_min, mjj_max, 4)
		elif fit_function == "f2":
			background_tf1 = TF1('background_tf1_f2',BackgroundFit_f2_trigcorr_bbh, mjj_min, mjj_max, 3)
		elif fit_function == "f3":
			background_tf1 = TF1('background_tf1_f3',BackgroundFit_f3_trigcorr_bbh, mjj_min, mjj_max, 3)
		elif fit_function == "f4":
			background_tf1 = TF1('background_tf1_f4',BackgroundFit_f4_trigcorr_bbh, mjj_min, mjj_max, 4)
		elif fit_function == "f5":
			background_tf1 = TF1('background_tf1_f5',BackgroundFit_f5_trigcorr_bbh, mjj_min, mjj_max, 3)
		elif fit_function == "f6":
			background_tf1 = TF1('background_tf1_f6',BackgroundFit_f6_trigcorr_bbh, mjj_min, mjj_max, 5)
		else:
			print "[make_background_tf1] ERROR : Unrecognized fit function " + fit_function
			sys.exit(1)
	return background_tf1


def rooplot(save_tag, fit_functions, background_workspace, fitted_signal_workspaces=None, expected_signal_workspaces=None, log=False, x_range=None, data_binning=None, normalization_bin_width=1, draw_chi2ndf=False, draw_chi2prob=False, data_histogram=None, trigger_correction=None, draw_trigeff=False):
	print "Making plot " + save_tag
	c = TCanvas("c_" + save_tag, "c_" + save_tag, 800, 1200)
	l = TLegend(0.65, 0.55, 0.88, 0.8)
	l.SetFillColor(0)
	l.SetBorderSize(0)
	top = TPad("top", "top", 0., 0.5, 1., 1.)
	top.SetBottomMargin(0.03)
	top.SetLeftMargin(0.15)
	top.Draw()
	if log:
		top.SetLogy()
	c.cd()
	bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
	bottom.SetLeftMargin(0.15)
	bottom.SetTopMargin(0.02)
	bottom.SetBottomMargin(0.2)
	bottom.Draw()
	ROOT.SetOwnership(c, False)
	ROOT.SetOwnership(top, False)
	ROOT.SetOwnership(bottom, False)
	top.cd()

	background_histograms = {}
	background_histograms_fitted_range = {}
	pull_histograms = {}
	pull_histograms_fitted_range = {}
	chi2s = {}
	ndfs = {}
	chi2ndfs = {}
	style_counter = 0
	first = True
	print "[rooplot] INFO : Opening workspace file {}".format(background_workspace)
	f_workspace = TFile(background_workspace, "READ")
	workspace = f_workspace.Get("w")
	workspace.Print()
	for fit_function in fit_functions:
		mjj = workspace.var("mjj")
		
		# First time: draw frame and data histogram
		if first:
			first = False
			if not data_histogram:
				fine_bins_array = array('d', range(int(x_min), int(x_max + 1)))
				if data_binning:
					roobinning = RooBinning(len(data_binning) - 1, data_binning, "data_binning")
				else:
					roobinning = RooBinning(len(fine_bins_array) - 1, fine_bins_array, "data_binning")
				fine_roobinning = RooBinning(len(fine_bins_array) - 1, fine_bins_array, "data_fine_binning")
				data_histogram = workspace.data("data_obs").createHistogram("data", mjj, RooFit.Binning(roobinning))
				data_histogram.SetDirectory(0)
				data_histogram_unrebinned = workspace.data("data_obs").createHistogram("data_fine", mjj, RooFit.Binning(fine_roobinning))
				data_histogram_unrebinned.SetDirectory(0)
			else:
				data_histogram_unrebinned = data_histogram.Clone()
				data_histogram_unrebinned.SetDirectory(0)
				data_histogram = histogram_tools.rebin_histogram(data_histogram, mass_bins, normalization_bin_width=1)

			if x_range:
				x_min = x_range[0]
				x_max = x_range[1]
			else:
				x_min = mjj.GetMin()
				x_max = mjj.GetMax()
			frame_top = TH1D("frame_top", "frame_top", 100, x_min, x_max)
			frame_top.SetDirectory(0)
			if log:
				frame_top.SetMaximum(data_histogram.GetMaximum() * 50.)
				frame_top.SetMinimum(0.1)
			else:
				frame_top.SetMaximum(data_histogram.GetMaximum() * 1.3)
				frame_top.SetMinimum(0.)
			if normalization_bin_width == 1:
				frame_top.GetYaxis().SetTitle("Events / GeV")
			else:
				frame_top.GetYaxis().SetTitle("Events / " + str(int(normalization_bin_width)) + " GeV")
			frame_top.GetYaxis().SetTitleSize(0.06)
			frame_top.GetYaxis().SetTitleOffset(1)
			frame_top.GetYaxis().SetLabelSize(0.06)
			frame_top.GetXaxis().SetTitleSize(0)
			frame_top.GetXaxis().SetLabelSize(0)
			frame_top.Draw()

			data_histogram.SetMarkerStyle(20)
			data_histogram.SetMarkerColor(1)
			data_histogram.SetMarkerSize(1)
			data_histogram.Draw("p same")
			l.AddEntry(data_histogram, "Data", "pl")

		# Make histogram from fitted background
		fitresult_name = "fitresult_model_" + fit_function + "_rooDatahist"
		fitresult = workspace.genobj(fitresult_name)
		fitresult.Print()
		fit = make_background_tf1_from_roofitresult(fit_function, fitresult, mjj_range=x_range, trigger_correction=None)
		print "[debug] Printing fit function"
		fit.Print()

		# Normalize fit
		scale_factor = workspace.var("background_" + fit_function + "_norm").getVal() / fit.Integral(mjj.getMin(), mjj.getMax())
		fit.SetParameter(0, fit.GetParameter(0) * scale_factor)
		fit.SetParError(0, fit.GetParError(0) * scale_factor)

		background_histograms[fit_function] = data_histogram.Clone()
		background_histograms[fit_function].Reset()
		background_histograms[fit_function].SetDirectory(0)
		background_histograms[fit_function].SetName(fit_function)
		for bin in xrange(1, background_histograms[fit_function].GetNbinsX() + 1):
			low_edge = background_histograms[fit_function].GetXaxis().GetBinLowEdge(bin)
			up_edge = background_histograms[fit_function].GetXaxis().GetBinUpEdge(bin)
			background_histograms[fit_function].SetBinContent(bin, fit.Integral(low_edge, up_edge) / (up_edge - low_edge) * normalization_bin_width)

		#if draw_trigeff:
		#	fit_notrigcorr = make_background_tf1_from_roofitresult(fit_function, fitresult, mjj_range=x_range, trigger_correction=None)
		#	scale_factor = workspace.var("background_" + fit_function + "_norm").getVal() / fit_notrigcorr.Integral(mjj.getMin(), mjj.getMax())
		#	fit_notrigcorr.SetParameter(0, fit_notrigcorr.GetParameter(0) * scale_factor)
		#	fit_notrigcorr.SetParError(0, fit_notrigcorr.GetParError(0) * scale_factor)
		#	background_histograms[fit_function + "_notrigcorr"] = data_histogram.Clone()
		#	background_histograms[fit_function + "_notrigcorr"].Reset()
		#	background_histograms[fit_function + "_notrigcorr"].SetDirectory(0)
		#	background_histograms[fit_function + "_notrigcorr"].SetName(fit_function)
		#	for bin in xrange(1, background_histograms[fit_function + "_notrigcorr"].GetNbinsX() + 1):
		#		low_edge = background_histograms[fit_function + "_notrigcorr"].GetXaxis().GetBinLowEdge(bin)
		#		up_edge = background_histograms[fit_function + "_notrigcorr"].GetXaxis().GetBinUpEdge(bin)
		#		background_histograms[fit_function + "_notrigcorr"].SetBinContent(bin, fit_notrigcorr.Integral(low_edge, up_edge) / (up_edge - low_edge) * normalization_bin_width)

		chi2s[fit_function] = 0.
		chi2ndfs[fit_function] = 0.
		ndfs[fit_function] = -1 * fit.GetNpar()
		for bin in xrange(1, data_histogram_unrebinned.GetNbinsX() + 1):
			low_edge = data_histogram_unrebinned.GetXaxis().GetBinLowEdge(bin)
			up_edge = data_histogram_unrebinned.GetXaxis().GetBinUpEdge(bin)
			if up_edge <= mjj.getMin() or low_edge >= mjj.getMax():
				continue
			# Data errors
			#this_chi2 = ((data_histogram_unrebinned.GetBinContent(bin) - fit.Integral(low_edge, up_edge)) / data_histogram_unrebinned.GetBinError(bin))**2
			# Fit errors
			this_fit = fit.Integral(low_edge, up_edge)
			this_chi2 = ((data_histogram_unrebinned.GetBinContent(bin) - this_fit) / (this_fit**0.5))**2
			chi2s[fit_function] += this_chi2
			ndfs[fit_function] += 1
			#print "[debug] Bin {}, chi2 = (({} - {})/{})**2=\t{}".format(bin, data_histogram_unrebinned.GetBinContent(bin), fit.Integral(low_edge, up_edge), this_fit**0.5, this_chi2)

		chi2ndfs[fit_function] = chi2s[fit_function] / ndfs[fit_function]
		print "#chi^2/NDF(" + fit_function + ") = " + str(round(chi2s[fit_function], 3)) + "/" + str(ndfs[fit_function]) + " = " + str(round(chi2ndfs[fit_function], 3)) + " / p = " + str(TMath.Prob(chi2s[fit_function], ndfs[fit_function]))

		# Draw total background histogram as a colored line
		background_histograms[fit_function].SetLineColor(seaborn.GetColorRoot("default", style_counter))
		background_histograms[fit_function].SetLineWidth(1)
		background_histograms[fit_function].SetLineStyle(1)
		background_histograms[fit_function].Draw("hist same")
		if draw_chi2ndf:
			legend_entry = fit_function + " (#chi^{2}/NDF=" + str(round(chi2ndfs[fit_function], 2)) + ")"
		elif draw_chi2prob:
			legend_entry = fit_function + " (p=" + str(round(ROOT.TMath.Prob(chi2s[fit_function], ndfs[fit_function]), 2)) + ")"
		else:
			legend_entry = fit_function
		l.AddEntry(background_histograms[fit_function], legend_entry, "l")
		background_histograms[fit_function].SetDirectory(0)

		if draw_trigeff:
			background_histograms[fit_function + "_notrigcorr"].SetLineColor(seaborn.GetColorRoot("pastel", style_counter))
			background_histograms[fit_function + "_notrigcorr"].SetLineWidth(1)
			background_histograms[fit_function + "_notrigcorr"].SetLineStyle(3)
			background_histograms[fit_function + "_notrigcorr"].Draw("hist same")
			legend_entry = fit_function + ", no trig corr"
			l.AddEntry(background_histograms[fit_function + "_notrigcorr"], legend_entry, "l")
			background_histograms[fit_function + "_notrigcorr"].SetDirectory(0)

		pull_histograms[fit_function] = data_histogram.Clone()
		pull_histograms[fit_function].Reset()
		pull_histograms[fit_function].SetDirectory(0)
		for bin in xrange(1, data_histogram.GetNbinsX() + 1):
			if data_histogram.GetBinError(bin) > 0:
				pull = (data_histogram.GetBinContent(bin) - background_histograms[fit_function].GetBinContent(bin)) / (data_histogram.GetBinError(bin))
			else:
				pull = 0.
			#print "[debug] Pull = " + str(pull)
			pull_histograms[fit_function].SetBinContent(bin, pull)
			pull_histograms[fit_function].SetBinError(bin, 0.)

		if draw_trigeff:
			pull_histograms[fit_function + "_notrigcorr"] = data_histogram.Clone()
			pull_histograms[fit_function + "_notrigcorr"].Reset()
			pull_histograms[fit_function + "_notrigcorr"].SetDirectory(0)
			for bin in xrange(1, data_histogram.GetNbinsX() + 1):
				if data_histogram.GetBinError(bin) > 0:
					pull = (data_histogram.GetBinContent(bin) - background_histograms[fit_function + "_notrigcorr"].GetBinContent(bin)) / (data_histogram.GetBinError(bin))
				else:
					pull = 0.
				#print "[debug] Pull = " + str(pull)
				pull_histograms[fit_function + "_notrigcorr"].SetBinContent(bin, pull)
				pull_histograms[fit_function + "_notrigcorr"].SetBinError(bin, 0.)

		pull_histograms_fitted_range[fit_function] = pull_histograms[fit_function].Clone()
		pull_histograms_fitted_range[fit_function].SetDirectory(0)
		for bin in xrange(1, pull_histograms_fitted_range[fit_function].GetNbinsX() + 1):
			bin_center = pull_histograms_fitted_range[fit_function].GetXaxis().GetBinCenter(bin)
			if bin_center < mjj.getMin() or bin_center > mjj.getMax():
				pull_histograms_fitted_range[fit_function].SetBinContent(bin, 0)
				pull_histograms_fitted_range[fit_function].SetBinError(bin, 0)

		style_counter += 1
	f_workspace.Close()

	l.Draw()
	Root.CMSLabelTwoPane(0.65, 0.84, "Internal", 1, 0.6)
	# Pull histogram
	c.cd()
	bottom.cd()

	#pull_histogram = frame_top.pullHist("Data", "B Fit")
	frame_bottom = TH1D("frame_bottom", "frame_bottom", 100, x_min, x_max)
	frame_bottom.SetMinimum(-5.)
	frame_bottom.SetMaximum(5.)
	#pull_histogram.plotOn(frame_bottom, RooFit.Name(fit_pdf_name))
	frame_bottom.GetXaxis().SetNdivisions(505)
	frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
	frame_bottom.GetXaxis().SetLabelSize(0.06)
	frame_bottom.GetXaxis().SetTitleSize(0.06)
	frame_bottom.GetXaxis().SetTitleOffset(1.0)
	frame_bottom.GetYaxis().SetTitle("#frac{Data - Fit}{#sigma(Data)}")
	frame_bottom.GetYaxis().SetLabelSize(0.06)
	frame_bottom.GetYaxis().SetTitleSize(0.06)
	frame_bottom.GetYaxis().SetTitleOffset(1.0)
	frame_bottom.Draw()

	style_counter = 0
	for fit_function in fit_functions:
		pull_histograms[fit_function].SetLineColor(seaborn.GetColorRoot("default", style_counter))
		pull_histograms[fit_function].SetLineWidth(2)
		pull_histograms[fit_function].SetLineStyle(2)
		pull_histograms[fit_function].Draw("hist same")

		if draw_trigeff:
			pull_histograms[fit_function].SetLineColor(seaborn.GetColorRoot("pastel", style_counter))
			pull_histograms[fit_function].SetLineWidth(1)
			pull_histograms[fit_function].SetLineStyle(3)
			pull_histograms[fit_function].Draw("hist same")

		pull_histograms_fitted_range[fit_function].SetLineColor(seaborn.GetColorRoot("dark", style_counter))
		pull_histograms_fitted_range[fit_function].SetLineWidth(2)
		pull_histograms_fitted_range[fit_function].SetLineStyle(1)
		pull_histograms_fitted_range[fit_function].Draw("hist same")
		style_counter += 1
	#pull_histogram.GetXaxis().SetTitle("m_{jj} [GeV]")
	#pull_histogram.GetYaxis().SetTitle("#frac{Data - Fit}{#sigma(Fit)}")
	#pull_histogram.Draw("same")
	c.cd()

	c.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Results/figures/c_" + save_tag + ".pdf")

def PlotAllSB(save_tag, fit_function, sb_names, sb_workspace_files, b_workspace_file, log=False, x_range=None, data_binning=None, normalization_bin_width=1, draw_chi2ndf=False, draw_chi2prob=False, data_histogram=None, trigger_correction=None, draw_trigeff=False):
	print "Making plot " + save_tag
	c = TCanvas("c_" + save_tag, "c_" + save_tag, 800, 1200)
	l = TLegend(0.5, 0.55, 0.88, 0.88)
	l.SetFillColor(0)
	l.SetBorderSize(0)
	top = TPad("top", "top", 0., 0.5, 1., 1.)
	top.SetBottomMargin(0.03)
	top.Draw()
	if log:
		top.SetLogy()
	c.cd()
	bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
	bottom.SetTopMargin(0.02)
	bottom.SetBottomMargin(0.2)
	bottom.Draw()
	ROOT.SetOwnership(c, False)
	ROOT.SetOwnership(top, False)
	ROOT.SetOwnership(bottom, False)
	top.cd()

	fit_histograms = {}
	f = TFile(b_workspace_file, "READ")
	workspace_b = f.Get("w")

	mjj = workspace_b.var("mjj")
	fitresult_name = "fitresult_model_" + fit_function + "_rooDatahist"
	if x_range:
		x_min = x_range[0]
		x_max = x_range[1]
	else:
		x_min = mjj.GetMin()
		x_max = mjj.GetMax()

	# Data histogram
	fine_bins_array = array('d', range(int(x_min), int(x_max + 1)))
	fine_roobinning = RooBinning(len(fine_bins_array) - 1, fine_bins_array, "data_fine_binning")
	if data_binning:
		roobinning = RooBinning(len(data_binning) - 1, data_binning, "data_binning")
	else:
		roobinning = RooBinning(len(fine_bins_array) - 1, fine_bins_array, "data_binning")
	data_histogram = workspace_b.data("data_obs").createHistogram("data", mjj, RooFit.Binning(roobinning))
	data_histogram.SetDirectory(0)

	frame_top = TH1D("frame_top", "frame_top", 100, x_min, x_max)
	frame_top.SetDirectory(0)
	if log:
		frame_top.SetMaximum(data_histogram.GetMaximum() * 1000.)
		frame_top.SetMinimum(0.1)
	else:
		frame_top.SetMaximum(data_histogram.GetMaximum() * 1.3)
		frame_top.SetMinimum(0.)
	frame_top.GetYaxis().SetTitle("Events / " + str(normalization_bin_width) + " GeV")
	frame_top.GetXaxis().SetTitleSize(0)
	frame_top.GetXaxis().SetLabelSize(0)
	frame_top.Draw()

	data_histogram.SetMarkerStyle(20)
	data_histogram.SetMarkerColor(1)
	data_histogram.SetMarkerSize(1)
	data_histogram.Draw("p same")
	l.AddEntry(data_histogram, "Data", "pl")

	# Background-only fit histogram
	fitresult_b = workspace_b.genobj(fitresult_name)
	#fitresult.Print()
	fit_b = make_background_tf1_from_roofitresult(fit_function, fitresult_b, mjj_range=x_range, trigger_correction=trigger_correction)
	fit_b.SetName("fit_Bonly")
	scale_factor = workspace_b.var("background_" + fit_function + "_norm").getVal() / fit_b.Integral(mjj.getMin(), mjj.getMax())
	fit_b.SetParameter(0, fit_b.GetParameter(0) * scale_factor)
	fit_b.SetParError(0, fit_b.GetParError(0) * scale_factor)
	fit_histograms["Background"] = data_histogram.Clone()
	fit_histograms["Background"].Reset()
	fit_histograms["Background"].SetDirectory(0)
	fit_histograms["Background"].SetName(fit_function)
	for bin in xrange(1, fit_histograms["Background"].GetNbinsX() + 1):
		low_edge = fit_histograms["Background"].GetXaxis().GetBinLowEdge(bin)
		up_edge = fit_histograms["Background"].GetXaxis().GetBinUpEdge(bin)
		fit_histograms["Background"].SetBinContent(bin, fit_b.Integral(low_edge, up_edge) / (up_edge - low_edge) * normalization_bin_width)
	fit_histograms["Background"].SetLineStyle(1)
	fit_histograms["Background"].SetLineColor(seaborn.GetColorRoot("cubehelix", 1, len(sb_names) + 1))
	fit_histograms["Background"].Draw("l same")
	l.AddEntry(fit_histograms["Background"], "Background only", "l")

	#S+B fit histograms
	style_counter = 2
	for sb_name in sb_names:
		f = TFile(sb_workspace_files[sb_name], "READ")
		workspace_sb = f.Get("w")
		fitresult_sb = workspace_sb.genobj(fitresult_name)
		#mjj = workspaces[sb_name].var("mjj")
		fit = make_background_tf1_from_roofitresult(fit_function, fitresult_sb, mjj_range=x_range, trigger_correction=trigger_correction)
		fit.SetName("fit_" + sb_name)
		scale_factor = (workspace_sb.var("background_" + fit_function + "_norm").getVal() + workspace_sb.var("signal_norm_" + fit_function).getVal()) / fit.Integral(mjj.getMin(), mjj.getMax())
		fit.SetParameter(0, fit.GetParameter(0) * scale_factor)
		fit.SetParError(0, fit.GetParError(0) * scale_factor)
		fit_histograms[sb_name] = data_histogram.Clone()
		fit_histograms[sb_name].Reset()
		fit_histograms[sb_name].SetDirectory(0)
		fit_histograms[sb_name].SetName("fit_hist_" + fit_function + "_" + sb_name)
		for bin in xrange(1, fit_histograms[sb_name].GetNbinsX() + 1):
			low_edge = fit_histograms[sb_name].GetXaxis().GetBinLowEdge(bin)
			up_edge = fit_histograms[sb_name].GetXaxis().GetBinUpEdge(bin)
			fit_histograms[sb_name].SetBinContent(bin, fit.Integral(low_edge, up_edge) / (up_edge - low_edge) * normalization_bin_width)
		fit_histograms[sb_name].SetLineStyle(style_counter % 10)
		fit_histograms[sb_name].SetLineColor(seaborn.GetColorRoot("cubehelix", style_counter, len(sb_names) + 1))
		fit_histograms[sb_name].Draw("l same")
		l.AddEntry(fit_histograms[sb_name], sb_name, "l")
		f.Close()
		style_counter += 1
	l.Draw()
	Root.CMSLabel(0.2, 0.8, "Internal", 1, 0.5)

	# Pull histogram
	c.cd()
	bottom.cd()

	#pull_histogram = frame_top.pullHist("Data", "B Fit")
	frame_bottom = TH1D("frame_bottom", "frame_bottom", 100, x_min, x_max)
	frame_bottom.SetMinimum(-5.)
	frame_bottom.SetMaximum(5.)
	#pull_histogram.plotOn(frame_bottom, RooFit.Name(fit_pdf_name))
	frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
	frame_bottom.GetYaxis().SetTitle("#frac{Data - Fit}{#sigma(Data)}")
	frame_bottom.Draw()

	style_counter = 2
	pull_histograms = {}
	for sb_name in sb_names:
		pull_histograms[sb_name] = data_histogram.Clone()
		pull_histograms[sb_name].SetName("pulls_" + sb_name)
		pull_histograms[sb_name].Reset()
		pull_histograms[sb_name].SetDirectory(0)
		for bin in xrange(1, pull_histograms[sb_name].GetNbinsX() + 1):
			data_int = data_histogram.GetBinContent(bin)
			data_err = data_histogram.GetBinError(bin)
			fit_int = fit_histograms[sb_name].GetBinContent(bin)
			if data_err > 0:
				pull_histograms[sb_name].SetBinContent(bin, (data_int - fit_int) / data_err)
				pull_histograms[sb_name].SetBinError(bin, 0.)
			else:
				pull_histograms[sb_name].SetBinContent(bin, 0.)
				pull_histograms[sb_name].SetBinError(bin, 0.)

		pull_histograms[sb_name].SetLineColor(seaborn.GetColorRoot("default", style_counter))
		pull_histograms[sb_name].SetLineWidth(2)
		pull_histograms[sb_name].SetLineStyle(2)
		pull_histograms[sb_name].Draw("hist same")

		style_counter += 1
	#pull_histogram.GetXaxis().SetTitle("m_{jj} [GeV]")
	#pull_histogram.GetYaxis().SetTitle("#frac{Data - Fit}{#sigma(Fit)}")
	#pull_histogram.Draw("same")
	c.cd()

	c.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Results/figures/c_" + save_tag + ".pdf")

def CorrectTriggerEfficiency(hist, analysis):
    hist_corr = hist.Clone()
    hist_corr.SetName(hist.GetName() + "_trigcorr")
    if "bbl" in analysis:
        mjj_min = 296
        mjj_max = 5000
    elif "bbh" in analysis:
        mjj_min = 526
        mjj_max = 5000
    bin_min = hist_corr.GetXaxis().FindBin(mjj_min + 1.e-5)
    bin_max = hist_corr.GetXaxis().FindBin(mjj_max - 1.e-5)

    # Fit background*eff
    if "bbl" in analysis:
        background_fit_times_eff = ROOT.TF1("tmp_background_fit_times_eff", BackgroundFit_f4_trigcorr_bbl, mjj_min, mjj_max, 4)
        background_fit_times_eff.SetParameter(1, 41.)
        background_fit_times_eff.SetParameter(2, -45.)
        background_fit_times_eff.SetParameter(3, 10.)
    elif "bbh" in analysis:
        background_fit_times_eff = ROOT.TF1("tmp_background_fit_times_eff", BackgroundFit_f4_trigcorr_bbh, mjj_min, mjj_max, 4)
        background_fit_times_eff.SetParameter(1, 35.)
        background_fit_times_eff.SetParameter(2, -28.)
        background_fit_times_eff.SetParameter(3, 10.)
    background_fit_times_eff.SetParameter(0, hist_corr.Integral(bin_min, bin_max))
    background_fit_times_eff.SetParameter(0, hist_corr.Integral(bin_min, bin_max) / background_fit_times_eff.Integral(mjj_min, mjj_max))
    hist_corr.Fit(background_fit_times_eff, "QR0")

    # Make background TF1
    background_fit = ROOT.TF1("tmp_background_fit", BackgroundFit_f4, mjj_min, mjj_max, 4)
    for i in xrange(4):
        background_fit.SetParameter(i, background_fit_times_eff.GetParameter(i))

    # Correct histogram bins
    for bin in xrange(1, hist_corr.GetNbinsX() + 1):
        if bin < bin_min or bin > bin_max:
            hist_corr.SetBinContent(bin, 0)
            hist_corr.SetBinError(bin, 0)
        else:
            old_content = hist_corr.GetBinContent(bin)
            old_error = hist_corr.GetBinError(bin)
            this_bin_min = hist_corr.GetXaxis().GetBinLowEdge(bin)
            this_bin_max = hist_corr.GetXaxis().GetBinUpEdge(bin)
            num = background_fit.Integral(this_bin_min, this_bin_max)
            den = background_fit_times_eff.Integral(this_bin_min, this_bin_max)
            if ("bbl" in analysis and this_bin_min <= 350.) or ("bbh" in analysis and this_bin_min <= 575.):
                print "mjj={}: observed={}, corr={}".format(0.5*(this_bin_min+this_bin_max), den, num)
            if den > 0:
                correction = num / den
                if ("bbl" in analysis and this_bin_min <= 350.) or ("bbh" in analysis and this_bin_min <= 575.):
                    print "Trigger correction for mjj={}: {}".format(0.5*(this_bin_min+this_bin_max), correction)
            else:
                correction = 1.
            hist_corr.SetBinContent(bin, old_content * correction)
            hist_corr.SetBinError(bin, old_error * correction)
    return hist_corr



if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Run and plot fits")
	parser.add_argument("--what", type=str, default="mjj", help="mjj = normal fit+data plot, sb = all s+b fits plus b only")
	parser.add_argument("--analyses", type=str, default="trigbbh_CSVTM,trigbbl_CSVTM", help='Analysis name (see analysis_configuration_8TeV.py)')
	parser.add_argument("--models", type=str, default="Hbb", help='Model name')
	parser.add_argument("--fit_functions", type=str, default="f1,f2,f3,f4,f5", help="Fit functions")
	parser.add_argument("--x_range", type=int, nargs=2, help="Plot xrange")
	parser.add_argument("--draw_trigeff", action="store_true", help="Plot with and without trigger efficiency (assumes create_datacards with run with correctTrigger)")
	parser.add_argument("--hide_chi2prob", action="store_false", help="Hide fit chi2 probabilities in legend.")
	parser.add_argument("--central", action="store_true", help="Draw central value of fit only")
	parser.add_argument("--correctTrigger", action="store_true", help="")
	parser.add_argument("--fitTrigger", action="store_true", help="")
	parser.add_argument("--fitOffB", action="store_true", help="")
	# Fit options
	parser.add_argument("-l", "--lumi", dest="lumi",
						default=19700., type=float,
						help="Integrated luminosity in pb-1 (default: %(default).1f)",
						metavar="LUMI")
	args = parser.parse_args()

	analyses = args.analyses.split(",")
	models = args.models.split(",")
	fit_functions = args.fit_functions.split(",")

	print "Plotting"
	#fitted_signal_workspaces = []
	#expected_signal_workspaces = []
	#if args.signal:
	#	for signal_model in args.signal.split(","):
	#		fitted_signal_workspaces.append(limit_paths.get_workspace_filename(args.analysis_name, signal_model))
	#if args.fixed_signal:
	#	for signal_model in args.fixed_signal.split(","):
	#		expected_signal_workspaces.append(limit_paths.get_workspace_filename(args.analysis_name, signal_model))

	mass_bins = array("d", [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5000]) #5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8000])
	if args.x_range:
		x_range = args.x_range
	else:
		x_range = [0., 2000.]

	if args.what == "mjj":
		for analysis in analyses:
			histogram_file = TFile(analysis_config.get_b_histogram_filename(analysis, "BJetPlusX_2012"), "READ")
			data_histogram_raw = histogram_file.Get("BHistograms/h_pfjet_mjj")
			data_histogram = CorrectTriggerEfficiency(data_histogram_raw, analysis)
			print "Data integral = {}".format(data_histogram.Integral())
			#if "trigbbh" in analysis:
			#	trigger_correction = "bbh"
			#elif "trigbbl" in analysis:
			#	trigger_correction = "bbl"
			data_histogram.SetDirectory(0)
			for model in models:
				background_workspace = limit_config.get_workspace_filename(analysis, model, 750, fitBonly=True, correctTrigger=args.correctTrigger, fitTrigger=args.fitTrigger, fitOffB=args.fitOffB)
				save_tag = "mjj_combinefits_" + analysis + "_" + model
				if len(fit_functions) == 1:
					save_tag += "_" + fit_functions[0]
				rooplot(save_tag, fit_functions, background_workspace, log=True, x_range=x_range, data_binning=mass_bins, normalization_bin_width=1., data_histogram=data_histogram, draw_chi2prob=args.hide_chi2prob, trigger_correction=None, draw_trigeff=args.draw_trigeff) # fitted_signal_workspaces=fitted_signal_workspaces, expected_signal_workspaces=expected_signal_workspaces, 
	elif args.what == "sb":
		for analysis in analyses:
			for fit_function in fit_functions:
				histogram_file = TFile(analysis_config.get_b_histogram_filename(analysis, "BJetPlusX_2012"), "READ")
				data_histogram = histogram_file.Get("BHistograms/h_pfjet_mjj")
				print "Data integral = {}".format(data_histogram.Integral())
				if "trigbbh" in analysis:
					trigger_correction = "bbh"
				elif "trigbbl" in analysis:
					trigger_correction = "bbl"
				data_histogram.SetDirectory(0)
				for model in models:
					b_workspace_file = limit_config.get_workspace_filename(analysis, model, 750, fitBonly=True, correctTrigger=False, fitTrigger=("trigbbl" in analysis))
					sb_names = []
					sb_workspaces = {}
					if "bbl" in analysis:
						signal_masses = xrange(400, 650, 50)
					else:
						signal_masses = xrange(600, 1200, 50)
					for signal_mass in signal_masses:
						sb_name = "m_{X}=" + str(signal_mass) + " [GeV]"
						sb_names.append(sb_name)
						sb_workspaces[sb_name] = limit_config.get_workspace_filename(analysis, model, signal_mass, fitBonly=False, correctTrigger=False, fitTrigger=("trigbbl" in analysis))
					save_tag = "mjj_sbfits_" + analysis + "_" + model + "_" + fit_function
					PlotAllSB(save_tag=save_tag, fit_function=fit_function, sb_names=sb_names, sb_workspace_files=sb_workspaces, b_workspace_file=b_workspace_file, log=True, x_range=x_range, data_binning=mass_bins, normalization_bin_width=1, draw_chi2ndf=False, draw_chi2prob=False, data_histogram=data_histogram, trigger_correction=None, draw_trigeff=False)
