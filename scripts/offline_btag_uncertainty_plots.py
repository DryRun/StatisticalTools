import os
import sys
import ROOT
from array import array
from ROOT import *

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

def plot(analysis, model, mass):
	f = TFile(analysis_config.get_b_histogram_filename(analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")), "READ")
	h_central = f.Get("BHistograms/h_pfjet_mjj")
	h_central.SetDirectory(0)
	h_central = histogram_tools.rebin_histogram(h_central, dijet_binning, normalization_bin_width=1.)
	h_up = f.Get("BHistograms/h_pfjet_mjj_BTagOfflineSFUp")
	h_up.SetDirectory(0)
	h_up = histogram_tools.rebin_histogram(h_up, dijet_binning, normalization_bin_width=1.)
	h_down = f.Get("BHistograms/h_pfjet_mjj_BTagOfflineSFDown")
	h_down.SetDirectory(0)
	h_down = histogram_tools.rebin_histogram(h_down, dijet_binning, normalization_bin_width=1.)

	c = TCanvas("c_offb_unc_{}_{}{}".format(analysis, model, mass), "c_offb_unc_{}_{}{}".format(analysis, model, mass), 800, 1000)
	top = TPad("top", "top", 0., 0.5, 1., 1.)
	top.SetBottomMargin(0.02)
	#top.SetLogy()
	top.Draw()
	top.cd()
	frame_top = TH1D("frame_top", "frame_top", 100, 0., mass*2)
	frame_top.GetXaxis().SetTitleSize(0)
	frame_top.GetXaxis().SetLabelSize(0)
	frame_top.GetYaxis().SetTitle("Events / GeV")
	frame_top.SetMinimum(0.01)
	frame_top.SetMaximum(h_central.GetMaximum() * 1.3)
	frame_top.Draw()
	h_central.SetMarkerStyle(20)
	h_central.Draw("same")
	h_up.SetMarkerStyle(26)
	h_up.Draw("same")
	h_down.SetMarkerStyle(32)
	h_down.Draw("same")
	l = TLegend(0.6, 0.6, 0.88, 0.8)
	l.SetFillColor(0)
	l.SetBorderSize(0)
	l.AddEntry(h_central, "Central", "p")
	l.AddEntry(h_up, "SF unc up", "p")
	l.AddEntry(h_down, "SF unc down", "p")
	l.Draw()

	c.cd()
	bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
	bottom.SetTopMargin(0.02)
	bottom.SetBottomMargin(0.2)
	bottom.Draw()
	bottom.cd()
	frame_bottom = TH1D("frame_bottom", "frame_bottom", 100, 0., mass*2)
	frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
	frame_bottom.GetYaxis().SetTitle("Unc / Central")
	frame_bottom.SetMinimum(0.5)
	frame_bottom.SetMaximum(1.5)
	frame_bottom.Draw()
	ratio_up = h_up.Clone()
	ratio_up.Divide(h_central)
	ratio_up.SetMarkerStyle(26)
	ratio_up.Draw("p hist same")
	ratio_down = h_down.Clone()
	ratio_down.Divide(h_central)
	ratio_down.SetMarkerStyle(32)
	ratio_down.Draw("p hist same")

	avg_up = TLine(0., h_up.Integral() / h_central.Integral(), mass*2, h_up.Integral() / h_central.Integral())
	avg_up.SetLineStyle(2)
	avg_up.SetLineColor(seaborn.GetColorRoot("pastel", 2))
	avg_up.Draw("same")
	avg_down = TLine(0., h_down.Integral() / h_central.Integral(), mass*2, h_down.Integral() / h_central.Integral())
	avg_down.SetLineStyle(2)
	avg_down.SetLineColor(seaborn.GetColorRoot("pastel", 2))
	avg_down.Draw("same")


	c.cd()
	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")

	f.Close()

def plot_all_masses(analysis, model, masses):
	ratio_up = {}
	ratio_down = {}
	for mass in masses:
		f = TFile(analysis_config.get_b_histogram_filename(analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")), "READ")
		h_central = f.Get("BHistograms/h_pfjet_mjj")
		h_central.SetDirectory(0)
		h_central = histogram_tools.rebin_histogram(h_central, dijet_binning, normalization_bin_width=1.)
		h_up = f.Get("BHistograms/h_pfjet_mjj_BTagOfflineSFUp")
		h_up.SetDirectory(0)
		h_up = histogram_tools.rebin_histogram(h_up, dijet_binning, normalization_bin_width=1.)
		h_down = f.Get("BHistograms/h_pfjet_mjj_BTagOfflineSFDown")
		h_down.SetDirectory(0)
		h_down = histogram_tools.rebin_histogram(h_down, dijet_binning, normalization_bin_width=1.)
		ratio_up[mass] = h_up.Clone()
		ratio_up[mass].SetName("ratio_up_{}".format(mass))
		ratio_up[mass].SetDirectory(0)
		ratio_up[mass].Divide(h_central)
		ratio_up[mass].SetMarkerStyle(26)
		ratio_up[mass].Draw("p hist same")
		ratio_down[mass] = h_down.Clone()
		ratio_down[mass].SetName("ratio_down_{}".format(mass))
		ratio_down[mass].SetDirectory(0)
		ratio_down[mass].Divide(h_central)
		ratio_down[mass].SetMarkerStyle(32)
		ratio_down[mass].Draw("p hist same")
		norm = h_central.Integral()
		print "{}\t{}\t{}:{}".format(analysis, model, mass, h_up.Integral()/norm-1.)
		f.Close()

	c = TCanvas("c_offb_unc_{}_{}_all".format(analysis, model), "c_offb_unc_{}_{}_all".format(analysis, model), 800, 600)
	frame = TH1D("frame", "frame", 100, 0., 1300.)
	frame.SetMinimum(0.8)
	frame.SetMaximum(1.2)
	frame.GetXaxis().SetTitle("m_{jj}")
	frame.GetYaxis().SetTitle("SF Variation / Central")
	frame.Draw()
	for style_counter, mass in enumerate(masses):
		ratio_up[mass].SetMarkerStyle(20)
		ratio_up[mass].SetMarkerColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(masses)+2))
		ratio_up[mass].SetLineColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(masses)+2))
		ratio_up[mass].Draw("hist same")
		ratio_down[mass].SetMarkerStyle(20)
		ratio_down[mass].SetMarkerColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(masses)+2))
		ratio_down[mass].SetLineColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(masses)+2))
		ratio_down[mass].Draw("hist same")
	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")



if __name__ == "__main__":
	masses = {"NoTrigger_eta1p7_CSVTM":[325,350,400,500,600,750], "NoTrigger_eta2p2_CSVTM":[600,750,900,1200]}

	for analysis in ["NoTrigger_eta1p7_CSVTM", "NoTrigger_eta2p2_CSVTM"]:
		for model in ["Hbb", "ZPrime", "RSG"]:
			for mass in masses[analysis]:
				plot(analysis, model, mass)
			plot_all_masses(analysis, model, masses[analysis])