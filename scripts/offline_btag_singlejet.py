import os
import sys
from array import array
import ROOT
from ROOT import *
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/src/CMSDIJET/QCDAnalysis/python/")
import analysis_configuration_8TeV as analysis_config
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/MyTools/RootUtils")
import histogram_tools
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
seaborn = Root.SeabornInterface()
seaborn.Initialize()
import time

# Make single jet b-tag efficiency histograms vs. pT vs. eta
# Events are taken from various JetHT histograms, with no prescale weighting. The HT slice is determined from m(jj).

ht_slices = ["HT200","HT250","HT300","HT350","HT400","HT450","HT500","HT550","HTUnprescaled"]
ht_ranges = {
	"HT200":[220, 386],
	"HT250":[386, 489],
	"HT300":[489, 526],
	"HT350":[526, 606],
	"HT400":[606, 649],
	"HT450":[649, 740],
	"HT500":[740, 788],
	"HT550":[788, 890],
	#"HT650":[890, 2000],
	"HTUnprescaled":[890,2000]
}
numerator_histograms = {}
denominator_histograms = {}
efficiency_histograms = {}
for sr in ["lowmass", "highmass"]:
	f_out = TFile("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/src/CMSDIJET/QCDAnalysis/data/single_jet_btag_eff_" + sr + ".root", "RECREATE")
	numerator_histograms[sr] = {}
	denominator_histograms[sr] = {}
	efficiency_histograms[sr] = {}
	if sr == "lowmass":
		analyses = {
			"HT200":"trigjetht200_eta1p7",
			"HT250":"trigjetht250_eta1p7",
			"HT300":"trigjetht300_eta1p7",
			"HT350":"trigjetht350_eta1p7",
			"HT400":"trigjetht400_eta1p7",
			"HT450":"trigjetht450_eta1p7",
			"HT500":"trigjetht500_eta1p7",
			"HT550":"trigjetht550_eta1p7",
			"HT650":"trigjetht650_eta1p7",
			"HTUnprescaled":"trigjetht_eta1p7"
		}
	elif sr == "highmass":
		analyses = {
			"HT200":"trigjetht200",
			"HT250":"trigjetht250",
			"HT300":"trigjetht300",
			"HT350":"trigjetht350",
			"HT400":"trigjetht400",
			"HT450":"trigjetht450",
			"HT500":"trigjetht500",
			"HT550":"trigjetht550",
			"HT650":"trigjetht650",
			"HTUnprescaled":"trigjetht"
		}
	for wp in ["CSVT", "CSVM", "CSVL"]:
		numerator_histograms[sr][wp] = {}
		denominator_histograms[sr][wp] = {}
		efficiency_histograms[sr][wp] = {}
		for which in [0, 1]:
			numerator_histograms[sr][wp][which] = None
			denominator_histograms[sr][wp][which] = None
			efficiency_histograms[sr][wp][which] = None
			# Load histograms
			first = True
			for ht_slice in ht_slices:
				f = TFile(analysis_config.get_b_histogram_filename(analyses[ht_slice], "JetHT_2012BCD"), "READ")
				this_numerator_histogram_3D = f.Get("BHistograms/h_pfjet_jet{}_pt_eta_{}".format(which, wp))
				this_denominator_histogram_3D = f.Get("BHistograms/h_pfjet_jet{}_pt_eta".format(which))
				zbin_min = this_numerator_histogram_3D.GetZaxis().FindBin(ht_ranges[ht_slice][0] + 1.e-5)
				zbin_max = this_numerator_histogram_3D.GetZaxis().FindBin(ht_ranges[ht_slice][0] - 1.e-5)
				if first:
					numerator_histograms[sr][wp][which] = TH2D(this_numerator_histogram_3D.GetName() + "_xy" + wp + str(which), this_numerator_histogram_3D.GetName() + "_xy", this_numerator_histogram_3D.GetXaxis().GetNbins(), this_numerator_histogram_3D.GetXaxis().GetXmin(), this_numerator_histogram_3D.GetXaxis().GetXmax(), this_numerator_histogram_3D.GetYaxis().GetNbins(), this_numerator_histogram_3D.GetYaxis().GetXmin(), this_numerator_histogram_3D.GetYaxis().GetXmax())
					numerator_histograms[sr][wp][which].SetDirectory(0)
					numerator_histograms[sr][wp][which].Sumw2()
					denominator_histograms[sr][wp][which] = TH2D(this_denominator_histogram_3D.GetName() + "_xy" + wp + str(which), this_denominator_histogram_3D.GetName() + "_xy", this_denominator_histogram_3D.GetXaxis().GetNbins(), this_denominator_histogram_3D.GetXaxis().GetXmin(), this_denominator_histogram_3D.GetXaxis().GetXmax(), this_denominator_histogram_3D.GetYaxis().GetNbins(), this_denominator_histogram_3D.GetYaxis().GetXmin(), this_denominator_histogram_3D.GetYaxis().GetXmax())
					denominator_histograms[sr][wp][which].SetDirectory(0)
					denominator_histograms[sr][wp][which].Sumw2()
					first = False
				for zbin in xrange(zbin_min, zbin_max+1):
					for xbin in xrange(1, this_numerator_histogram_3D.GetNbinsX()+1):
						for ybin in xrange(1, this_numerator_histogram_3D.GetNbinsY()+1):
							num_old_bin_content = numerator_histograms[sr][wp][which].GetBinContent(xbin, ybin)
							num_new_bin_content = num_old_bin_content + this_numerator_histogram_3D.GetBinContent(xbin, ybin, zbin)
							den_old_bin_content = denominator_histograms[sr][wp][which].GetBinContent(xbin, ybin)
							den_new_bin_content = den_old_bin_content + this_denominator_histogram_3D.GetBinContent(xbin, ybin, zbin)

							num_old_bin_error = numerator_histograms[sr][wp][which].GetBinError(xbin, ybin)
							num_new_bin_error = (num_old_bin_content**2 + this_numerator_histogram_3D.GetBinError(xbin, ybin, zbin)**2)**0.5
							den_old_bin_error = denominator_histograms[sr][wp][which].GetBinError(xbin, ybin)
							den_new_bin_error = (den_old_bin_content**2 + this_denominator_histogram_3D.GetBinError(xbin, ybin, zbin)**2)**0.5

							numerator_histograms[sr][wp][which].SetBinContent(xbin, ybin, num_new_bin_content)
							numerator_histograms[sr][wp][which].SetBinError(xbin, ybin, num_new_bin_error)
							denominator_histograms[sr][wp][which].SetBinContent(xbin, ybin, den_new_bin_content)
							denominator_histograms[sr][wp][which].SetBinError(xbin, ybin, den_new_bin_error)
				f.Close()
			efficiency_histograms[sr][wp][which] = numerator_histograms[sr][wp][which].Clone()
			efficiency_histograms[sr][wp][which].SetName("effb_" + "_" + wp + "_jet" + str(which))
			efficiency_histograms[sr][wp][which].Divide(numerator_histograms[sr][wp][which], denominator_histograms[sr][wp][which], 1, 1, "B")
		# Average
		numerator_histograms[sr][wp]["both"] = numerator_histograms[sr][wp][0].Clone()
		numerator_histograms[sr][wp]["both"].SetName("numerator_" + sr + wp + "both")
		numerator_histograms[sr][wp]["both"].Add(numerator_histograms[sr][wp][1])
		denominator_histograms[sr][wp]["both"] = denominator_histograms[sr][wp][0].Clone()
		denominator_histograms[sr][wp]["both"].SetName("denominator_" + sr + wp + "both")
		denominator_histograms[sr][wp]["both"].Add(denominator_histograms[sr][wp][1])
		efficiency_histograms[sr][wp]["both"] = numerator_histograms[sr][wp]["both"].Clone()
		efficiency_histograms[sr][wp]["both"].Divide(numerator_histograms[sr][wp]["both"], denominator_histograms[sr][wp]["both"], 1, 1, "B")
		for which in [0, 1, "both"]:
			c = TCanvas("c_single_btag_eff_{}_{}_jet{}".format(sr, wp, which), "c_single_btag_eff_{}_{}_jet{}".format(sr, wp, which), 800, 600)
			c.SetRightMargin(0.2)
			c.SetLogz()
			efficiency_histograms[sr][wp][which].GetXaxis().SetTitle("p_{T} [GeV]")
			efficiency_histograms[sr][wp][which].GetYaxis().SetTitle("#eta")
			efficiency_histograms[sr][wp][which].GetZaxis().SetTitle("b-tag #epsilon")
			efficiency_histograms[sr][wp][which].Draw("colz")
			c.SaveAs(analysis_config.figure_directory + "/OfflineBTag/" + c.GetName() + ".pdf")
		for which in [0,1]:
			f_out.cd()
			efficiency_histograms[sr][wp][which].Write()
	f_out.Close()