import os
import sys
import ROOT
import copy
from ROOT import *
gROOT.SetBatch(True)
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/RooCBPlusVoigtian.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
seaborn = Root.SeabornInterface()
seaborn.Initialize()
from CMSDIJET.StatisticalTools.systematics import *
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
import CMSDIJET.StatisticalTools.trigger_efficiency as trigger_efficiency

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/CMSDIJET/QCDAnalysis")
import analysis_configuration_8TeV as analysis_config

from array import array
dijet_binning = array("d", [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8000])
dijet_roobinning = RooBinning(len(dijet_binning) - 1, dijet_binning, "dijet_binning")


def MakeSignalPlot(canvas_name, histogram_names, histograms, x_range=None, y_range=None, logx=False, logy=False, legend_coordinates=[0.6,0.6,0.88,0.88], x_title=None, y_title=None):
	c = TCanvas(canvas_name, canvas_name, 800, 600)
	l = TLegend(legend_coordinates[0], legend_coordinates[1], legend_coordinates[2], legend_coordinates[3])
	l.SetFillColor(0)
	l.SetBorderSize(0)

	# Frame
	if x_range:
		x_min = x_range[0]
		x_max = x_range[1]
	else:
		x_min = 1.e20
		x_max = -1.e20
		for histogram_name, histogram in histograms.iteritems():
			if histogram.GetXaxis().GetXmin() < x_min:
				x_min = histogram.GetXaxis().GetXmin()
			if histogram.GetXaxis().GetXmax() > x_max:
				x_max = histogram.GetXaxis().GetXmax()
	if y_range:
		y_min = y_range[0]
		y_max = y_range[1]
	else:
		y_min = 1.e20
		y_max = -1.e20
		for histogram_name, histogram in histograms.iteritems():
			if histogram.GetMinimum() < y_min:
				y_min = histogram.GetMinimum()
			if histogram.GetMaximum() > y_max:
				y_max = histogram.GetMaximum()
		if logy:
			y_min = y_min / 5.
			y_max = y_max * 5.
		else:
			y_min = y_min - (y_max - y_min)*0.1
			y_max = y_max + (y_max - y_min)*0.1
	frame = TH1F("frame", "frame", 100, x_min, x_max)
	frame.SetMinimum(y_min)
	frame.SetMaximum(y_max)
	if x_title:
		frame.GetXaxis().SetTitle(x_title)
	else:
		frame.GetXaxis().SetTitle(histograms[histogram_names[0]].GetXaxis().GetTitle())
	if y_title:
		frame.GetYaxis().SetTitle(y_title)
	else:
		frame.GetYaxis().SetTitle(histograms[histogram_names[0]].GetYaxis().GetTitle())
	frame.Draw("axis")

	style_counter = 0
	for histogram_name in histogram_names:
		histograms[histogram_name].SetMarkerStyle(20)
		histograms[histogram_name].SetMarkerSize(0)
		histograms[histogram_name].SetLineWidth(2)
		if len(histogram_names) < 6:
			histograms[histogram_name].SetLineColor(seaborn.GetColorRoot("default", style_counter))
			histograms[histogram_name].SetFillColor(seaborn.GetColorRoot("pastel", style_counter))
			histograms[histogram_name].SetFillStyle(3004)
			l.AddEntry(histograms[histogram_name], histogram_name, "lf")
			histograms[histogram_name].DrawCopy("e3 same")
			histograms[histogram_name].SetFillStyle(0)
			histograms[histogram_name].DrawCopy("hist same")
			histograms[histogram_name].SetFillStyle(3004)
		else:
			histograms[histogram_name].SetLineColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(histogram_names)))
			l.AddEntry(histograms[histogram_name], histogram_name, "l")
			histograms[histogram_name].Draw("hist same")
		style_counter += 1
	l.Draw()
	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")

def GetLegendEntry(analysis, model, mass):
	legend_entry = ""
	if "trigbbl" in analysis:
		legend_entry += "Low mass / "
	elif "trigbbh" in analysis: 
		legend_entry += "High mass / "
	else:
		legend_entry += analysis + " / "
	if model == "Hbb":
		legend_entry += "H / "
	elif model == "RSG":
		legend_entry += "G / "
	else:
		legend_entry += model + " / "
	legend_entry += str(mass) + " GeV"
	return legend_entry

if __name__ == "__main__":
	from argparse import ArgumentParser
	parser = ArgumentParser(description='Make various plots of signal distributions')
	parser.add_argument("--analyses", type=str, help="Analyses")
	parser.add_argument("--models", type=str, default="Hbb,RSG", help="Models (comma-separated list, options are Hbb, RSG, Zp")
	parser.add_argument("--masses", type=str, help="Masses")
	parser.add_argument("--prefix", type=str, default="BHistograms/h_pfjet_", help="Prefix of histograms in histogram file")
	parser.add_argument("--variables", type=str, help="1D variables to plot (suffix of h_pfjet_* in BHistograms")
	parser.add_argument("--x_range", type=float, nargs=2, help="Custom x range")
	parser.add_argument("--y_range", type=float, nargs=2, help="Custom y range")
	parser.add_argument("--logx", action="store_true", help="Log x axis")
	parser.add_argument("--logy", action="store_true", help="Log y axis")
	parser.add_argument("--normalize", action="store_true", help="Normalize histograms")
	parser.add_argument("--rebin", type=int, help="Rebin histograms")
	args = parser.parse_args()

	variables = args.variables.split(",")
	analyses = args.analyses.split(",")
	models = args.models.split(",")
	masses = [int(x) for x in args.masses.split(",")]

	for variable in variables:
		histograms = {}
		histogram_names = []
		for analysis in analyses:
			for model in models:
				for mass in masses:
					histogram_name = GetLegendEntry(analysis, model, mass)
					f = TFile(analysis_config.get_b_histogram_filename(analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")), "READ")
					histogram_names.append(histogram_name)
					histograms[histogram_name] = f.Get(args.prefix + variable)
					if not histograms[histogram_name]:
						print "[signal_plots] ERROR : Histogram {} not found in file {}".format(args.prefix + variable, analysis_config.get_b_histogram_filename(analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")))
					histograms[histogram_name].SetDirectory(0)
					f.Close()
					if args.normalize:
						if histograms[histogram_name].Integral() > 0:
							histograms[histogram_name].Scale(1. / histograms[histogram_name].Integral())
					if args.rebin:
						histograms[histogram_name].Rebin(args.rebin)
		canvas_name = "c_signal_distribution_{}_{}_{}_{}".format(variable, "_".join(analyses), "_".join(models), "_".join(args.masses.split(",")))
		MakeSignalPlot(canvas_name, histogram_names, histograms, x_range=args.x_range, y_range=args.y_range, logx=args.logx, logy=args.logy)
	