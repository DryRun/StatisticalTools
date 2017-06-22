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

	f_old = TFile(analysis_config.get_b_histogram_filename(analysis, old_sample), "READ")
	h_old = f_old.Get("BHistograms/h_pfjet_mjj")
	h_old.SetName("h_pfjet_mjj_oldJEC")
	h_old.SetDirectory(0)
	h_old = histogram_tools.rebin_histogram(h_old, dijet_binning, normalization_bin_width=1.)
	print "Old integral = {}".format(h_old.Integral())

	c = TCanvas("c_jeccomparison_{}".format(analysis), "c_jeccomparison_{}".format(analysis), 1000, 800)
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

if __name__ == "__main__":
	data_jes_plot("BJetPlusXJEC13_2012", "BJetPlusX_2012", "trigbbl_CSVTM")
	data_jes_plot("BJetPlusXJEC13_2012", "BJetPlusX_2012", "trigbbh_CSVTM")
