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
import MyTools.RootUtils.root_plots as root_plots

signal_regions = ["lowmass", "highmass"]
strategies = ["mctrigger", "datatrigger"]
mc_analyses = {
	"lowmass":{
		"mctrigger":"trigbbl_CSVTM",
		"datatrigger":"NoTrigger_eta1p7_CSVTM"
	}, 
	"highmass":{
		"mctrigger":"trigbbh_CSVTM",
		"datatrigger":"NoTrigger_eta2p2_CSVTM"
	}
}
data_analyses = {"lowmass":"trigbbl_CSVTM", "highmass":"trigbbh_CSVTM"}

dijet_binning = array("d", [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5000]) #5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8000])

for signal_region in signal_regions:
	f_data = TFile(analysis_config.get_b_histogram_filename(data_analyses[signal_region], "BJetPlusX_2012"), "READ")
	data_hist_raw = f_data.Get("BHistograms/h_pfjet_mjj")
	data_hist_raw.SetName("mjj_data_" + signal_region)
	data_hist_raw.SetDirectory(0)
	f_data.Close()

	for strategy in strategies:
		f_mc = TFile(analysis_config.get_b_histogram_filename(mc_analyses[signal_region][strategy], "QCD_TuneZ2star_8TeV_pythia6"), "READ")
		mc_hist = f_mc.Get("BHistograms/h_pfjet_mjj")
		mc_hist.SetName("mjj_mc_" + signal_region + "_" + strategy)
		mc_hist.SetDirectory(0)
		if strategy == "datatrigger":
			if signal_region == "lowmass":
				mc_hist.Scale(trigger_efficiency.online_btag_eff["trigbbl_CSVTM"][0])
			elif signal_region == "highmass":
				mc_hist.Scale(trigger_efficiency.online_btag_eff["trigbbh_CSVTM"][0])

			data_hist = data_hist_raw.Clone()
			for bin in xrange(1, data_hist.GetNbinsX() + 1):
				mjj = data_hist.GetXaxis().GetBinCenter(bin)
				if signal_region == "lowmass":
					data_hist.SetBinContent(bin, data_hist.GetBinContent(bin) / trigger_efficiency.trigger_efficiency_bbl(mjj))
					data_hist.SetBinError(bin, data_hist.GetBinError(bin) / trigger_efficiency.trigger_efficiency_bbl(mjj))
				elif signal_region == "highmass":
					data_hist.SetBinContent(bin, data_hist.GetBinContent(bin) / trigger_efficiency.trigger_efficiency_bbh(mjj))
					data_hist.SetBinError(bin, data_hist.GetBinError(bin) / trigger_efficiency.trigger_efficiency_bbh(mjj))
				# Zero histograms below where the trigger correction is known
				if signal_region == "lowmass":
					if mjj < 176.:
						data_hist.SetBinContent(bin, 0.)
						data_hist.SetBinError(bin, 0.)
						mc_hist.SetBinContent(bin, 0.)
						mc_hist.SetBinError(bin, 0.)
				elif signal_region == "highmass":
					if mjj < 325.:
						data_hist.SetBinContent(bin, 0.)
						data_hist.SetBinError(bin, 0.)
						mc_hist.SetBinContent(bin, 0.)
						mc_hist.SetBinError(bin, 0.)
		elif strategy == "mctrigger":
			data_hist = data_hist_raw.Clone()


		plot = root_plots.DataMCPlot()
		plot.add_data(data_hist, "Data 2012")
		plot.add_mc(mc_hist, "QCD MC", color=17)
		plot.draw(complex_rebinning=dijet_binning, logy=True, save_tag = "data_vs_MC_" + signal_region + "_" + strategy, draw_ratio=True, y_title="Events / GeV", cms_label="Internal", ratio_range=[0.,4.], x_range=[0.,2000.])
