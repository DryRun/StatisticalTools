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

def data_jes_plot(new_sample, old_sample, analysis):
	print "Making data JES plot for {}".format(analysis)
	f_new = TFile(analysis_config.get_b_histogram_filename(analysis, new_sample), "READ")
	h_new = f_new.Get("BHistograms/h_pfjet_mjj")
	h_new.SetName("h_pfjet_mjj_newJEC")
	h_new.SetDirectory(0)
	h_new = histogram_tools.rebin_histogram(h_new, dijet_binning, normalization_bin_width=1.)
	print "New integral = {}".format(h_new.Integral())

	f_old = TFile(analysis_config.get_b_histogram_filename(analysis, old_sample).replace("EightTeeEeVeeBee/BHistograms", "EightTeeEeVeeBee/BHistograms/old/JEC2012"), "READ")
	h_old = f_old.Get("BHistograms/h_pfjet_mjj")
	h_old.SetName("h_pfjet_mjj_oldJEC")
	h_old.SetDirectory(0)
	h_old = histogram_tools.rebin_histogram(h_old, dijet_binning, normalization_bin_width=1.)
	print "Old integral = {}".format(h_old.Integral())

	c = TCanvas("c_jeccomparison_{}".format(analysis), "c_jeccomparison_{}".format(analysis), 800, 1200)
	top = TPad("top", "top", 0., 0.5, 1., 1.)
	top.SetBottomMargin(0.02)
	top.SetLogy()
	top.Draw()
	top.cd()
	h_new.SetMarkerStyle(20)
	h_new.GetXaxis().SetTitleSize(0)
	h_new.GetXaxis().SetLabelSize(0)
	h_new.GetYaxis().SetTitle("Events / GeV")
	h_new.Draw()
	h_old.SetMarkerStyle(25)
	h_old.Draw("same")

	c.cd()
	bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
	bottom.SetTopMargin(0.02)
	bottom.SetBottomMargin(0.2)
	bottom.Draw()
	bottom.cd()
	h_ratio = h_new.Clone()
	h_ratio.Divide(h_old)
	h_ratio.GetXaxis().SetTitle("m_{jj} [GeV]")
	h_ratio.GetYaxis().SetTitle("New JEC / Old JEC")
	h_ratio.SetMarkerStyle(21)
	h_ratio.Draw()

	c.cd()
	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")


	f_new.Close()
	f_old.Close()

def mc_jes_plot(analysis, model, mass):
	f_new = TFile(analysis_config.get_b_histogram_filename(analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")), "READ")
	h_new = f_new.Get("BHistograms/h_pfjet_mjj")
	h_new.SetName("h_pfjet_mjj_newJEC")
	h_new.SetDirectory(0)
	h_new = histogram_tools.rebin_histogram(h_new, dijet_binning, normalization_bin_width=1.)

	f_old = TFile(analysis_config.get_b_histogram_filename(analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")).replace("EightTeeEeVeeBee/BHistograms", "EightTeeEeVeeBee/BHistograms/old/JEC2012"), "READ")
	h_old = f_old.Get("BHistograms/h_pfjet_mjj")
	h_old.SetName("h_pfjet_mjj_oldJEC")
	h_old.SetDirectory(0)
	h_old = histogram_tools.rebin_histogram(h_old, dijet_binning, normalization_bin_width=1.)

	c = TCanvas("c_jeccomparison_{}_{}_{}".format(analysis, model, mass), "c_jeccomparison_{}_{}_{}".format(analysis, model, mass), 800, 1200)
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
	frame_top.SetMaximum(h_new.GetMaximum() * 1.3)
	frame_top.Draw()
	h_new.SetMarkerStyle(20)
	h_new.GetXaxis().SetTitleSize(0)
	h_new.GetXaxis().SetLabelSize(0)
	h_new.GetYaxis().SetTitle("Events / GeV")
	h_new.Draw("same")
	h_old.SetMarkerStyle(25)
	h_old.Draw("same")
	l = TLegend(0.6, 0.6, 0.88, 0.8)
	l.SetFillColor(0)
	l.SetBorderSize(0)
	l.AddEntry(h_new, "New JEC")
	l.AddEntry(h_old, "Old JEC")
	l.Draw()

	c.cd()
	bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
	bottom.SetTopMargin(0.02)
	bottom.SetBottomMargin(0.2)
	bottom.Draw()
	bottom.cd()
	frame_bottom = TH1D("frame_bottom", "frame_bottom", 100, 0., mass*2)
	frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
	frame_bottom.GetYaxis().SetTitle("New JEC / Old JEC")
	frame_bottom.SetMinimum(0.5)
	frame_bottom.SetMaximum(1.5)
	frame_bottom.Draw()
	h_ratio = h_new.Clone()
	h_ratio.Divide(h_old)
	h_ratio.SetMarkerStyle(21)
	h_ratio.Draw("same")

	c.cd()
	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")


	f_new.Close()
	f_old.Close()


if __name__ == "__main__":
	data_jes_plot("BJetPlusXJEC13_2012", "BJetPlusX_2012", "trigbbl_CSVTM")
	data_jes_plot("BJetPlusXJEC13_2012", "BJetPlusX_2012", "trigbbh_CSVTM")

	masses = {"NoTrigger_eta1p7_CSVTM":[350,400,500,600], "NoTrigger_eta2p2_CSVTM":[600,750,900,1200]}

	for analysis in ["NoTrigger_eta1p7_CSVTM", "NoTrigger_eta2p2_CSVTM"]:
		for model in ["Hbb", "ZPrime", "RSG"]:
			for mass in masses[analysis]:
				mc_jes_plot(analysis, model, mass)
