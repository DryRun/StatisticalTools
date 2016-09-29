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
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
seaborn = Root.SeabornInterface()
seaborn.Initialize()

def BackgroundFit_f1(x, par):
	return par[0] * (1. - (x[0] / 8.e3))**par[1] / ((x[0] / 8.e3)**(par[2] + par[3] * TMath.Log((x[0] / 8.e3))))

def BackgroundFit_f2(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3))**par[2]

def BackgroundFit_f3(x, par):
	return par[0] / (1 + par[1] * (x[0] / 8.e3))**par[2]

def BackgroundFit_f4(x, par):
	if 1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2 < 0:
		return 0
	else:
		return par[0] / ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2)**par[3])

def BackgroundFit_f5(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3)**(1./3.))**par[2]

def make_background_tf1_from_roofitresult(fit_function, roofitresult, mjj_range=None):
	print "[make_background_tf1_from_roofitresult] INFO : make_background_tf1_from_roofitresult for fit function " + fit_function
	print mjj_range
	roofitresult.Print()
	function = make_background_tf1(fit_function, mjj_range)
	function.SetParameter(0, roofitresult.floatParsFinal().find("background_" + fit_function + "_norm").getVal())
	function.SetParError(0, roofitresult.floatParsFinal().find("background_" + fit_function + "_norm").getError())

	for i in xrange(1, function.GetNpar()):
		function.SetParameter(i, roofitresult.floatParsFinal().find(fit_function + "_p" + str(i)).getVal())
		function.SetParError(i, roofitresult.floatParsFinal().find(fit_function + "_p" + str(i)).getError())

	return function

def make_background_tf1(fit_function, mjj_range):
	mjj_min = mjj_range[0]
	mjj_max = mjj_range[1]

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
	else:
		print "[make_background_tf1] ERROR : Unrecognized fit function " + fit_function
		sys.exit(1)
	return background_tf1


def rooplot(save_tag, fit_functions, background_workspace, fitted_signal_workspaces=None, expected_signal_workspaces=None, log=False, x_range=None, data_binning=None, normalization_bin_width=1, draw_chi2ndf=False, data_histogram=None):
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

	background_histograms = {}
	background_histograms_fitted_range = {}
	pull_histograms = {}
	pull_histograms_fitted_range = {}
	chi2s = {}
	chi2ndfs = {}
	style_counter = 0
	first = True
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
			frame_top.GetYaxis().SetTitle("Events / " + str(normalization_bin_width) + " GeV")
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
		fit = make_background_tf1_from_roofitresult(fit_function, fitresult, mjj_range=x_range)

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

		chi2s[fit_function] = 0.
		chi2ndfs[fit_function] = 0.
		ndf = -1 * fit.GetNpar()
		for bin in xrange(1, data_histogram_unrebinned.GetNbinsX() + 1):
			low_edge = data_histogram_unrebinned.GetXaxis().GetBinLowEdge(bin)
			up_edge = data_histogram_unrebinned.GetXaxis().GetBinUpEdge(bin)
			if up_edge < mjj.getMin() or low_edge > mjj.getMax():
				continue
			if data_histogram_unrebinned.GetBinError(bin):
				chi2s[fit_function] += ((data_histogram_unrebinned.GetBinContent(bin) - fit.Integral(low_edge, up_edge)) / data_histogram_unrebinned.GetBinError(bin))**2
				ndf += 1
		chi2ndfs[fit_function] = chi2s[fit_function] / ndf
		print "#chi^2/NDF(" + fit_function + ") = " + str(round(chi2s[fit_function], 3)) + "/" + str(ndf) + " = " + str(round(chi2ndfs[fit_function], 3)) + " / p = " + str(TMath.Prob(chi2s[fit_function], ndf))

		# Draw total background histogram as a colored line
		background_histograms[fit_function].SetLineColor(seaborn.GetColorRoot("default", style_counter))
		background_histograms[fit_function].SetLineWidth(1)
		background_histograms[fit_function].SetLineStyle(1)
		background_histograms[fit_function].Draw("hist same")
		if draw_chi2ndf:
			legend_entry = "Background " + fit_function + "(#chi^{2}/NDF=" + str(round(chi2ndfs[fit_function], 2)) + ")"
		else:
			legend_entry = "Background " + fit_function
		l.AddEntry(background_histograms[fit_function], legend_entry, "l")
		background_histograms[fit_function].SetDirectory(0)

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
	Root.CMSLabel()
	
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

	style_counter = 0
	for fit_function in fit_functions:
		pull_histograms[fit_function].SetLineColor(seaborn.GetColorRoot("default", style_counter))
		pull_histograms[fit_function].SetLineWidth(2)
		pull_histograms[fit_function].SetLineStyle(2)
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

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Run and plot fits")
	parser.add_argument("--analyses", type=str, default="trigbbh_CSVTM,trigbbl_CSVTM", help='Analysis name (see analysis_configuration_8TeV.py)')
	parser.add_argument("--models", type=str, default="Hbb,RSG", help='Model name')
	parser.add_argument("--plot", action="store_true", help="Plot mjj spectra and fits. Background fit is always plotted; signal fits are plotted if --signal is specified.")
	parser.add_argument("--x_range", type=int, nargs=2, help="Plot xrange")
	# Fit options
	parser.add_argument("-l", "--lumi", dest="lumi",
						default=19700., type=float,
						help="Integrated luminosity in pb-1 (default: %(default).1f)",
						metavar="LUMI")

	args = parser.parse_args()

	analyses = args.analyses.split(",")
	models = args.models.split(",")

	if args.plot:
		print "Plotting"
		#fitted_signal_workspaces = []
		#expected_signal_workspaces = []
		#if args.signal:
		#	for signal_model in args.signal.split(","):
		#		fitted_signal_workspaces.append(limit_paths.get_workspace_filename(args.analysis_name, signal_model))
		#if args.fixed_signal:
		#	for signal_model in args.fixed_signal.split(","):
		#		expected_signal_workspaces.append(limit_paths.get_workspace_filename(args.analysis_name, signal_model))

		fit_functions = ["f1", "f2", "f3", "f4", "f5"]
		mass_bins = array("d", [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8000])
		if args.x_range:
			x_range = args.x_range
		else:
			x_range = [0., 2000.]

		for analysis in analyses:
			histogram_file = TFile(analysis_config.get_b_histogram_filename(analysis, "BJetPlusX_2012"), "READ")
			data_histogram = histogram_file.Get("BHistograms/h_pfjet_mjj")
			data_histogram.SetDirectory(0)
			for model in models:
				background_workspace = limit_config.get_workspace_filename(analysis, model, 750, fitBonly=False, fitSignal=True)
				rooplot("mjj_combinefits_" + analysis + "_" + model, fit_functions, background_workspace, log=True, x_range=x_range, data_binning=mass_bins, normalization_bin_width=1., data_histogram=data_histogram, draw_chi2ndf=True) # fitted_signal_workspaces=fitted_signal_workspaces, expected_signal_workspaces=expected_signal_workspaces, 

