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

def MakeMjjComparisonPlot(hist_num, name_num, hist_den, name_den, x_range=None, y_range=None):
	hist_num.Scale(1. / hist_num.Integral())
	hist_den.Scale(1. / hist_den.Integral())
	for bin in xrange(1, hist_num.GetNbinsX() + 1):
		if hist_num.GetBinContent(bin) < 0:
			hist_num.SetBinContent(bin, 0)
			hist_num.SetBinError(bin, 0)
		if hist_den.GetBinContent(bin) < 0:
			hist_den.SetBinContent(bin, 0)
			hist_den.SetBinError(bin, 0)

	# Binning
	#mass_bins = array("d", [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8000])
	#hist_num = histogram_tools.rebin_histogram(hist_num, mass_bins, normalization_bin_width=1)
	#hist_den = histogram_tools.rebin_histogram(hist_den, mass_bins, normalization_bin_width=1)

	for bin in xrange(1, hist_num.GetNbinsX() + 1):
		if hist_num.GetBinContent(bin) < 0:
			hist_num.SetBinContent(bin, 0)
			hist_num.SetBinError(bin, 0)
		if hist_den.GetBinContent(bin) < 0:
			hist_den.SetBinContent(bin, 0)
			hist_den.SetBinError(bin, 0)


	c = TCanvas("c_mjj_comparison_" + name_num + "_" + name_den, "c_mjj_comparison_" + name_num + "_" + name_den, 800, 1200)
	c.SetLogy()
	l = TLegend(0.6, 0.65, 0.88, 0.88)
	l.SetFillColor(0)
	l.SetBorderSize(0)
	top = TPad("top", "top", 0., 0.5, 1., 1.)
	top.SetBottomMargin(0.03)
	top.Draw()
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

	# Frame
	if x_range:
		x_min = x_range[0]
		x_max = x_range[1]
	else:
		x_min = hist_num.GetXaxis().GetXmin()
		x_max = hist_den.GetXaxis().GetXmax()
	if y_range:
		y_min = y_range[0]
		y_max = y_range[1]
	else:
		y_min = 0.1
		y_max = max(hist_den.GetMaximum(), hist_den.GetMaximum())
		y_max = y_max * 10.
	frame_top = TH1D("frame_top", "frame_top", 100, x_min, x_max)
	frame_top.SetMinimum(y_min)
	frame_top.SetMaximum(y_max)
	print "[debug] Frame y range = {} {}".format(frame_top.GetMinimum(), frame_top.GetMaximum())
	#frame_top.GetXaxis().SetTitle("m_{jj} [GeV]")
	frame_top.GetXaxis().SetTitleSize(0)
	frame_top.GetXaxis().SetLabelSize(0)
	frame_top.GetYaxis().SetTitle("A.U.")
	frame_top.Draw("axis")

	hist_num.SetLineColor(seaborn.GetColorRoot("dark", 0))
	hist_num.Draw("hist same")
	l.AddEntry(hist_num, name_num, "l")
	hist_den.SetLineColor(seaborn.GetColorRoot("dark", 1))
	hist_den.Draw("hist same")
	l.AddEntry(hist_den, name_den, "l")
	l.Draw()

	# Denominator: derivative plots
	derivative_min = 1.e20
	derivative_max = -1.e20
	graph_dnum = TGraph(hist_num.GetNbinsX() - 1)
	for point in xrange(hist_num.GetNbinsX() - 1):
		x = (hist_num.GetXaxis().GetBinCenter(point) + hist_num.GetXaxis().GetBinCenter(point+1)) / 2.
		dx = (hist_num.GetXaxis().GetBinCenter(point+1) - hist_num.GetXaxis().GetBinCenter(point))
		dy = (hist_num.GetBinContent(point+1) - hist_num.GetBinContent(point))
		print "[debug] Point {}, x={}, dx={}, dy={}".format(point, x, dx, dy)
		if dx > 0:
			graph_dnum.SetPoint(point, x, dy/dx)
			if dy/dx < derivative_min:
				derivative_min = dy/dx
			if dy/dx > derivative_max:
				derivative_max = dy/dx
		else:
			graph_dnum.SetPoint(point, x, 0.)

	graph_dden = TGraph(hist_den.GetNbinsX() - 1)
	for point in xrange(hist_den.GetNbinsX() - 1):
		x = (hist_den.GetXaxis().GetBinCenter(point) + hist_den.GetXaxis().GetBinCenter(point+1)) / 2.
		dx = (hist_den.GetXaxis().GetBinCenter(point+1) - hist_den.GetXaxis().GetBinCenter(point))
		dy = (hist_den.GetBinContent(point+1) - hist_den.GetBinContent(point))
		if dx > 0:
			graph_dden.SetPoint(point, x, dy/dx)
			if dy/dx < derivative_min:
				derivative_min = dy/dx
			if dy/dx > derivative_max:
				derivative_max = dy/dx
		else:
			graph_dden.SetPoint(point, x, 0.)
	c.cd()
	bottom.cd()
	frame_bottom = TH1D("frame_bottom", "frame_bottom", 100, x_min, x_max)
	frame_bottom.SetMinimum(derivative_min - (derivative_max - derivative_min)/4.)
	frame_bottom.SetMaximum(derivative_max + (derivative_max - derivative_min)/4.)
	frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
	frame_bottom.GetYaxis().SetTitle("Derivative")
	frame_bottom.Draw()

	graph_dnum.SetLineColor(seaborn.GetColorRoot("dark", 0))
	graph_dnum.SetLineWidth(2)
	graph_dnum.Draw("l")
	graph_dden.SetLineColor(seaborn.GetColorRoot("dark", 1))
	graph_dden.SetLineWidth(2)
	graph_dden.Draw("l")
	##hist_ratio = hist_num.Clone()
	#hist_ratio = hist_num.Clone()
	#hist_ratio.Reset()
	#for bin in xrange(1, hist_ratio.GetNbinsX() + 1):
	#	num = hist_num.GetBinContent(bin)
	#	den = hist_den.GetBinContent(bin)
	#	if den > 0 and num > 0 and num <= den:
	#		ratio = num / den
	#		dratio = (ratio * (1. - ratio) / den)**0.5
	#	else:
	#		ratio = 0.
	#		dratio = 0.
	#	hist_ratio.SetBinContent(bin, ratio)
	#	hist_ratio.SetBinError(bin, dratio)
	#hist_ratio.SetLineColor(1)
	#hist_ratio.SetLineWidth(1)
	#hist_ratio.SetMarkerStyle(20)
	#hist_ratio.SetMarkerSize(1)
	#hist_ratio.SetMarkerColor(1)
	#hist_ratio.Draw("same")

	c.cd()
	c.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Results/figures/{}.pdf".format(c.GetName()))

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Run and plot fits")
	parser.add_argument("analysis_num", type=str, default="trigbbl_raw", help='Numerator analysis')
	parser.add_argument("analysis_den", type=str, default="trigbbl_CSVTM", help='Denominator analysis')
	args = parser.parse_args()

	f_num = TFile(analysis_config.get_b_histogram_filename(args.analysis_num, "BJetPlusX_2012"), "READ")
	h_num = f_num.Get("BHistograms/h_pfjet_mjj")
	h_num.SetDirectory(0)
	h_num.Rebin(5)
	f_den = TFile(analysis_config.get_b_histogram_filename(args.analysis_den, "BJetPlusX_2012"), "READ")
	h_den = f_den.Get("BHistograms/h_pfjet_mjj")
	h_den.SetDirectory(0)
	h_den.Rebin(5)
	print "[debug] h_num.Integral() = {}, h_den.Integral() = {}".format(h_num.Integral(), h_den.Integral())

	MakeMjjComparisonPlot(h_num, args.analysis_num, h_den, args.analysis_den, x_range=[0,1000], y_range=[1.e-5, 1.e-1])