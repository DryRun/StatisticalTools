import os
import sys
import array
from math import sqrt
import ROOT
from ROOT import *
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/src/CMSDIJET/QCDAnalysis/python/")
import analysis_configuration_8TeV as analysis_config
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/MyTools/RootUtils")
import histogram_tools
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
seaborn = Root.SeabornInterface()
seaborn.Initialize()

gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)


class TrigEffPlotter():
	def __init__(self):
		self._analyses = []
		self._efficiency_combinations = {}
		self._mjj_histograms = {}
		self._mjj_histograms_csvorder = {}
		self._mjj_histograms_vetothirdjet = {}
		self._efficiency_histograms = {}
		self._efficiency_histograms_csvorder = {}
		self._efficiency_histograms_vetothirdjet = {}
		self._mjj_histograms_fine = {}
		self._efficiency_histograms_fine = {}
		self._mass_bins = array.array("d", [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8000])
		self.GetJetHTHistogram()
		self.GetSingleMuHistograms()
		self.GetBJetPlusXHistograms()
		self.MakeEfficiencyHistograms()
		self.FitEfficiencies()

	def GetJetHTHistogram(self):
		for sr_name in ["highmass", "lowmass"]:
			analyses = {}
			HT_slices = []
			for mass in xrange(200, 600, 50):
				HT_slices.append("HT" + str(mass))
				analyses["HT" + str(mass)] = "trigjetht" + str(mass)
				if sr_name == "lowmass":
					analyses["HT" + str(mass)] += "_eta1p7"
				analyses["HT" + str(mass)] += "_CSVTM"
			sample = "JetHT_2012BCD"
			HT_slice_histograms = {}
			for HT_slice in HT_slices:
				f = TFile(analysis_config.get_b_histogram_filename(analyses[HT_slice], sample), "READ")
				HT_slice_histograms[HT_slice] = f.Get("BHistograms/h_pfjet_mjj")
				print "On file " + analysis_config.get_b_histogram_filename(analyses[HT_slice], sample)
				HT_slice_histograms[HT_slice].SetName(HT_slice_histograms[HT_slice].GetName() + "_" + analyses[HT_slice])
				HT_slice_histograms[HT_slice].SetDirectory(0)
				f.Close()
			HT_slices.append("HTUnprescaled")
			unprescaled_analysis_name = "trigjetht"
			if sr_name == "lowmass":
				unprescaled_analysis_name += "_eta1p7"
			unprescaled_analysis_name += "_CSVTM"
			analyses["HTUnprescaled"] = unprescaled_analysis_name
			f_unprescaled = TFile(analysis_config.get_b_histogram_filename(unprescaled_analysis_name, sample), "READ")
			HT_slice_histograms["HTUnprescaled"] = f_unprescaled.Get("BHistograms/h_pfjet_mjj")
			HT_slice_histograms["HTUnprescaled"].SetName(HT_slice_histograms["HTUnprescaled"].GetName() + "_" + analyses["HTUnprescaled"])
			HT_slice_histograms["HTUnprescaled"].SetDirectory(0)
			f_unprescaled.Close()
			ranges = {
				"HT200":[220, 386],
				"HT250":[386, 489],
				"HT300":[489, 526],
				"HT350":[526, 606],
				"HT400":[606, 649],
				"HT450":[649, 740],
				"HT500":[740, 788],
				"HT550":[788, 890],
				#"HT650":[800, 890],
				"HTUnprescaled":[890, 2000]
			}

			self._analyses.append("JetHT")
			self._mjj_histograms_fine["JetHT"] = self.FrankenHist(HT_slices, HT_slice_histograms, ranges)
			self._mjj_histograms["JetHT"] = histogram_tools.rebin_histogram(self._mjj_histograms_fine["JetHT"], self._mass_bins, normalization_bin_width=1)
			self._mjj_histograms_fine["JetHT"].Rebin(5)

	# Cobble together m(jj) histogram from several JetHT triggers. 
	def FrankenHist(self, HT_slices, histograms, ranges):
		frankenhist = histograms[HT_slices[0]].Clone()
		frankenhist.SetDirectory(0)
		frankenhist.Sumw2()
		frankenhist.Reset()
		for bin in xrange(1, frankenhist.GetNbinsX() + 1):
			x = frankenhist.GetXaxis().GetBinCenter(bin)
			# Find range
			for HT_slice in HT_slices:
				if x > ranges[HT_slice][0] and x < ranges[HT_slice][1]:
					frankenhist.SetBinContent(bin, histograms[HT_slice].GetBinContent(bin))
					frankenhist.SetBinError(bin, histograms[HT_slice].GetBinError(bin))
					break
		return frankenhist

	def GetSingleMuHistograms(self):
		#this_analyses = ["trigmu_highmass_CSVTM", "trigmu_lowmass_CSVTM", "trigmubbh_highmass_CSVTM", "trigmubbl_lowmass_CSVTM", "trigmubbll_lowmass_CSVTM", "trigmu24i_lowmass_CSVTM", "trigmu24ibbl_lowmass_CSVTM", "trigmu24ibbll_lowmass_CSVTM", "trigmu40_lowmass_CSVTM", "trigmu40bbl_lowmass_CSVTM", "trigmu40bbll_lowmass_CSVTM"]
		this_analyses = ["trigmu24i_lowmass_CSVTM", "trigmu24ibbl_lowmass_CSVTM"]
		this_analyses.extend(["trigmu24i_highmass_CSVTM", "trigmu24ibbh_highmass_CSVTM"])
		for analysis in this_analyses:
			print "Opening " + analysis_config.get_b_histogram_filename(analysis, "SingleMu_2012")
			f = ROOT.TFile(analysis_config.get_b_histogram_filename(analysis, "SingleMu_2012"), "READ")
			self._mjj_histograms_fine[analysis] = f.Get("BHistograms/h_pfjet_mjj")
			self._mjj_histograms_fine[analysis].SetName("h_" + analysis + "_mjj_fine")
			self._mjj_histograms_fine[analysis].SetDirectory(0)
			self._mjj_histograms[analysis] = histogram_tools.rebin_histogram(self._mjj_histograms_fine[analysis], self._mass_bins, normalization_bin_width=1)
			self._mjj_histograms[analysis].SetName("h_" + analysis + "_mjj")
			self._mjj_histograms[analysis].SetDirectory(0)

			if not "highmass" in analysis:
				self._mjj_histograms_csvorder[analysis] = f.Get("BHistograms/h_pfjet_mjj_csvorder")
				self._mjj_histograms_csvorder[analysis].SetName("h_" + analysis + "_mjj_csvorder")
				self._mjj_histograms_csvorder[analysis].SetDirectory(0)
				self._mjj_histograms_csvorder[analysis] = histogram_tools.rebin_histogram(self._mjj_histograms_csvorder[analysis], self._mass_bins, normalization_bin_width=1)

				self._mjj_histograms_vetothirdjet[analysis] = f.Get("BHistograms/h_pfjet_mjj_vetothirdjet")
				self._mjj_histograms_vetothirdjet[analysis].SetName("h_" + analysis + "_mjj_vetothirdjet")
				self._mjj_histograms_vetothirdjet[analysis].SetDirectory(0)
				self._mjj_histograms_vetothirdjet[analysis] = histogram_tools.rebin_histogram(self._mjj_histograms_vetothirdjet[analysis], self._mass_bins, normalization_bin_width=1)

			if "bbll" in analysis:
				self._mjj_histograms_fine[analysis].Scale(1.7) # Prescale for singlemu + 60/53. The prescale was not computer for these analyses.
				self._mjj_histograms[analysis].Scale(1.7) # Prescale for singlemu + 60/53. The prescale was not computer for these analyses.
			self._mjj_histograms_fine[analysis].Rebin(5)
			self._analyses.append(analysis)
			f.Close()


	def GetBJetPlusXHistograms(self):
		this_analyses = ["trigbbh_CSVTM", "trigbbl_CSVTM", "trigbbll_CSVTM", "trigbbh_trigbbl_CSVTM"]
		this_samples = ["BJetPlusX_2012", "BJetPlusX_2012BCD"]
		for analysis in this_analyses:
			for sample in this_samples:
				name = analysis + "_" + sample
				f = ROOT.TFile(analysis_config.get_b_histogram_filename(analysis, sample), "READ")
				self._mjj_histograms_fine[name] = f.Get("BHistograms/h_pfjet_mjj")
				self._mjj_histograms_fine[name].SetName("h_" + name + "_mjj_fine")
				self._mjj_histograms_fine[name].SetDirectory(0)
				self._mjj_histograms[name] = histogram_tools.rebin_histogram(self._mjj_histograms_fine[name], self._mass_bins, normalization_bin_width=1)
				self._mjj_histograms[name].SetName("h_" + name + "_mjj")
				self._mjj_histograms[name].SetDirectory(0)
				self._mjj_histograms_csvorder[name] = f.Get("BHistograms/h_pfjet_mjj_csvorder")
				self._mjj_histograms_csvorder[name].SetName("h_" + name + "_mjj_csvorder")
				self._mjj_histograms_csvorder[name].SetDirectory(0)
				self._mjj_histograms_csvorder[name] = histogram_tools.rebin_histogram(self._mjj_histograms_csvorder[name], self._mass_bins, normalization_bin_width=1)
				self._mjj_histograms_vetothirdjet[name] = f.Get("BHistograms/h_pfjet_mjj_vetothirdjet")
				self._mjj_histograms_vetothirdjet[name].SetName("h_" + name + "_mjj_vetothirdjet")
				self._mjj_histograms_vetothirdjet[name].SetDirectory(0)
				self._mjj_histograms_vetothirdjet[name] = histogram_tools.rebin_histogram(self._mjj_histograms_vetothirdjet[name], self._mass_bins, normalization_bin_width=1)
				self._mjj_histograms_fine[name].Rebin(5)
				self._analyses.append(name)
				f.Close()

				if "bbll" in analysis:
					self._mjj_histograms[name].Scale(1.7)
					self._mjj_histograms_fine[name].Scale(1.7)

	def MakeEfficiencyHistograms(self):
		self._efficiency_combinations.update({
			"JetHT_highmass":["trigbbh_CSVTM_BJetPlusX_2012BCD", "JetHT"],
			"JetHT_lowmass":["trigbbl_CSVTM_BJetPlusX_2012BCD", "JetHT"], 
			#"JetHT_llowmass":["trigbbll_CSVTM_BJetPlusX_2012BCD", "JetHT"],
			#"SingleMu_highmass":["trigmubbh_highmass_CSVTM", "trigmu_highmass_CSVTM"], 
			#"SingleMu_lowmass":["trigmubbl_lowmass_CSVTM", "trigmu_lowmass_CSVTM"], 
			#"SingleMu_llowmass":["trigmubbll_lowmass_CSVTM", "trigmu_lowmass_CSVTM"], 
			"SingleMu24i_highmass":["trigmu24ibbh_highmass_CSVTM", "trigmu24i_highmass_CSVTM"], 
			"SingleMu24i_lowmass":["trigmu24ibbl_lowmass_CSVTM", "trigmu24i_lowmass_CSVTM"], 
			#"SingleMu24i_llowmass":["trigmu24ibbll_lowmass_CSVTM", "trigmu24i_lowmass_CSVTM"], 
			#"SingleMu40_highmass":["trigmu40bbh_highmass_CSVTM", "trigmu40_highmass_CSVTM"], 
			#"SingleMu40_lowmass":["trigmu40bbl_lowmass_CSVTM", "trigmu40_lowmass_CSVTM"], 
			#"SingleMu40_llowmass":["trigmu40bbll_lowmass_CSVTM", "trigmu40_lowmass_CSVTM"], 
			"BJet60_53_lowmass":["trigbbl_CSVTM_BJetPlusX_2012", "trigbbll_CSVTM_BJetPlusX_2012"], 
			"BJet60_53_highmass":["trigbbh_CSVTM_BJetPlusX_2012", "trigbbll_CSVTM_BJetPlusX_2012"], 
			"BJet80_70_highmass":["trigbbh_trigbbl_CSVTM_BJetPlusX_2012", "trigbbl_CSVTM_BJetPlusX_2012"], 
			})
		self._legend_entries = {
			"trigbbh_CSVTM":"160/120 + b-tag",
			"trigbbl_CSVTM":"80/70 + b-tag",
			"trigbbll_CSVTM":"60/53 + b-tag",
			"trigbbh_CSVTM_BJetPlusX_2012":"160/120 + b-tag",
			"trigbbh_trigbbl_CSVTM_BJetPlusX_2012":"((160/120) && (80/70))",
			"trigbbl_CSVTM_BJetPlusX_2012":"80/70 + b-tag",
			"trigbbll_CSVTM_BJetPlusX_2012":"60/53 + b-tag",
			"trigbbh_CSVTM_BJetPlusX_2012BCD":"160/120 + b-tag",
			"trigbbl_CSVTM_BJetPlusX_2012BCD":"80/70 + b-tag",
			"trigbbll_CSVTM_BJetPlusX_2012BCD":"60/53 + b-tag",
			"trigmu_highmass_CSVTM":"mu24i||mu40",
			"trigmu_lowmass_CSVTM":"mu24i||mu40",
			"trigmubbh_highmass_CSVTM":"(mu24i||mu40)&&(160/120+b-tag)",
			"trigmubbl_lowmass_CSVTM":"(mu24i||mu40)&&(80/70+b-tag)",
			"trigmubbll_lowmass_CSVTM":"(mu24i||mu40)&&(60/53+b-tag)", 
			"trigmu24ibbh_highmass_CSVTM":"(mu24i)&&(160/120+b-tag)",
			"trigmu24i_lowmass_CSVTM":"mu24i",
			"trigmu24ibbl_lowmass_CSVTM":"(mu24i)&&(80/70+b-tag)",
			"trigmu24ibbll_lowmass_CSVTM":"(mu24i)&&(60/53+b-tag)", 
			"trigmu40bbh_highmass_CSVTM":"(mu40)&&(160/120+b-tag)",
			"trigmu40bbl_lowmass_CSVTM":"(mu40)&&(80/70+b-tag)",
			"trigmu40bbll_lowmass_CSVTM":"(mu40)&&(60/53+b-tag)", 
			"JetHT":"JetHT"
		}
		self._efficiency_guesses = {
				"JetHT_highmass":0.5,
				"JetHT_lowmass":0.2,
				"JetHT_llowmass":0.1,
				"SingleMu_highmass":0.2,
				"SingleMu_lowmass":0.1,
				"SingleMu_llowmass":0.05,
				"SingleMu24i_highmass":0.2,
				"SingleMu24i_lowmass":0.18,
				"SingleMu24i_llowmass":0.05,
				"SingleMu40_highmass":0.2,
				"SingleMu40_lowmass":0.1,
				"SingleMu40_llowmass":0.05,
				"BJet60_53_lowmass":1,
				"BJet60_53_highmass":1,
				"BJet80_70_highmass":1,
		}
		self._colors = {
			"JetHT_highmass":seaborn.GetColorRoot("default", 0),
			"JetHT_lowmass":seaborn.GetColorRoot("default", 1),
			"JetHT_llowmass":seaborn.GetColorRoot("default", 2),
			"SingleMu_highmass":seaborn.GetColorRoot("default", 3),
			"SingleMu_lowmass":seaborn.GetColorRoot("default", 3),
			"SingleMu_llowmass":seaborn.GetColorRoot("default", 5),
			"SingleMu24i_highmass":seaborn.GetColorRoot("pastel", 3),
			"SingleMu24i_lowmass":seaborn.GetColorRoot("pastel", 3),
			"SingleMu24i_llowmass":seaborn.GetColorRoot("pastel", 5),
			"SingleMu40_highmass":seaborn.GetColorRoot("dark", 3),
			"SingleMu40_lowmass":seaborn.GetColorRoot("dark", 3),
			"SingleMu40_llowmass":seaborn.GetColorRoot("dark", 5),
			"BJet60_53_lowmass":seaborn.GetColorRoot("cubehelix", 0, 30),
			"BJet60_53_highmass":seaborn.GetColorRoot("cubehelix", 10, 30),
			"BJet80_70_highmass":seaborn.GetColorRoot("cubehelix", 20, 30),
		}
		self._line_styles = {
			"JetHT_highmass":1,
			"JetHT_lowmass":1,
			"JetHT_llowmass":1,
			"SingleMu_highmass":1,
			"SingleMu_lowmass":1,
			"SingleMu_llowmass":1,
			"SingleMu24i_highmass":2,
			"SingleMu24i_lowmass":2,
			"SingleMu24i_llowmass":2,
			"SingleMu40_highmass":3,
			"SingleMu40_lowmass":3,
			"SingleMu40_llowmass":3,
			"BJet60_53_lowmass":1,
			"BJet60_53_highmass":1,
			"BJet80_70_highmass":1,
		}
		for efficiency_name, hist_pair in self._efficiency_combinations.iteritems():
			print "[MakeEfficiencyHistograms] INFO : Making efficiency histogram " + efficiency_name
			print "[MakeEfficiencyHistograms] INFO : \tNumerator (" + hist_pair[0] + ") integral = " + str(self._mjj_histograms[hist_pair[0]].Integral())
			print "[MakeEfficiencyHistograms] INFO : \tDenominator (" + hist_pair[1] + ") integral = " + str(self._mjj_histograms[hist_pair[1]].Integral())
			self._efficiency_histograms[efficiency_name] = self._mjj_histograms[hist_pair[0]].Clone()
			self._efficiency_histograms[efficiency_name].SetName("h_efficiency_" + efficiency_name)
			self._efficiency_histograms[efficiency_name].SetDirectory(0)

			if "JetHT" in efficiency_name or efficiency_name == "BJet60_53_highmass":
				self._efficiency_histograms[efficiency_name].Divide(self._mjj_histograms[hist_pair[1]])
				# Set the histogram errors to N/n * sqrt((N+n(x-2))/(nN)), where x is the prescale
				for bin in xrange(1, self._efficiency_histograms[efficiency_name].GetNbinsX() + 1):
					bin_width = self._efficiency_histograms[efficiency_name].GetXaxis().GetBinWidth(bin)
					n = self._mjj_histograms[hist_pair[0]].GetBinContent(bin) * bin_width
					N = self._mjj_histograms[hist_pair[1]].GetBinContent(bin) * bin_width
					print "n = " + str(n) + " / N = " + str(N) + " / dN = " + str(self._mjj_histograms[hist_pair[1]].GetBinError(bin) * bin_width)
					if N == 0:
						continue
					x = 1. * (self._mjj_histograms[hist_pair[1]].GetBinError(bin) * bin_width)**2 / N
					print "Bin center = " + str(self._efficiency_histograms[efficiency_name].GetXaxis().GetBinCenter(bin)) + " / x=" + str(x)
					if x == 0 or n == 0:
						continue
					err = n / N * sqrt((N + n * (x - 2)) / (n * N))
					self._efficiency_histograms[efficiency_name].SetBinError(bin, err)
			else:
				self._efficiency_histograms[efficiency_name].Divide(self._mjj_histograms[hist_pair[0]], self._mjj_histograms[hist_pair[1]], 1, 1, "B")
			self._efficiency_histograms[efficiency_name].SetMarkerStyle(20)
			self._efficiency_histograms[efficiency_name].SetMarkerSize(1)
			self._efficiency_histograms[efficiency_name].SetMarkerColor(self._colors[efficiency_name])
			self._efficiency_histograms[efficiency_name].SetLineWidth(2)
			self._efficiency_histograms[efficiency_name].SetLineStyle(self._line_styles[efficiency_name])
			self._efficiency_histograms[efficiency_name].SetLineColor(self._colors[efficiency_name])

			self._mjj_histograms[hist_pair[0]].SetMarkerStyle(24)
			self._mjj_histograms[hist_pair[0]].SetMarkerColor(self._colors[efficiency_name])
			self._mjj_histograms[hist_pair[0]].SetLineColor(self._colors[efficiency_name])
			self._mjj_histograms[hist_pair[1]].SetMarkerStyle(20)
			self._mjj_histograms[hist_pair[1]].SetMarkerColor(self._colors[efficiency_name])
			self._mjj_histograms[hist_pair[1]].SetLineColor(self._colors[efficiency_name])

			# Finely binned
			self._efficiency_histograms_fine[efficiency_name] = self._mjj_histograms_fine[hist_pair[0]].Clone()
			self._efficiency_histograms_fine[efficiency_name].SetName("h_efficiency_" + efficiency_name)
			self._efficiency_histograms_fine[efficiency_name].SetDirectory(0)
			if "JetHT" in efficiency_name or efficiency_name == "BJet60_53_highmass":
				self._efficiency_histograms_fine[efficiency_name].Divide(self._mjj_histograms_fine[hist_pair[1]])
			else:
				self._efficiency_histograms_fine[efficiency_name].Divide(self._mjj_histograms_fine[hist_pair[0]], self._mjj_histograms_fine[hist_pair[1]], 1, 1, "B")
			self._efficiency_histograms_fine[efficiency_name].SetMarkerStyle(20)
			self._efficiency_histograms_fine[efficiency_name].SetMarkerSize(1)
			self._efficiency_histograms_fine[efficiency_name].SetMarkerColor(self._colors[efficiency_name])
			self._efficiency_histograms_fine[efficiency_name].SetLineWidth(2)
			self._efficiency_histograms_fine[efficiency_name].SetLineStyle(self._line_styles[efficiency_name])
			self._efficiency_histograms_fine[efficiency_name].SetLineColor(self._colors[efficiency_name])

			self._mjj_histograms_fine[hist_pair[0]].SetMarkerStyle(24)
			self._mjj_histograms_fine[hist_pair[0]].SetMarkerColor(self._colors[efficiency_name])
			self._mjj_histograms_fine[hist_pair[1]].SetMarkerStyle(20)
			self._mjj_histograms_fine[hist_pair[1]].SetMarkerColor(self._colors[efficiency_name])


			# Variations with no 3rd jet or CSV ordering
			if self._mjj_histograms_csvorder.has_key(hist_pair[0]) and self._mjj_histograms_csvorder.has_key(hist_pair[1]):
				self._efficiency_histograms_csvorder[efficiency_name] = self._mjj_histograms[hist_pair[0]].Clone()
				self._efficiency_histograms_csvorder[efficiency_name].SetName("h_efficiency_" + efficiency_name + "_csvorder")
				self._efficiency_histograms_csvorder[efficiency_name].SetDirectory(0)
				if "JetHT" in efficiency_name:
					self._efficiency_histograms_csvorder[efficiency_name].Divide(self._mjj_histograms[hist_pair[1]])
				else:
					self._efficiency_histograms_csvorder[efficiency_name].Divide(self._mjj_histograms[hist_pair[0]], self._mjj_histograms[hist_pair[1]], 1, 1, "B")
				self._efficiency_histograms_csvorder[efficiency_name].SetMarkerStyle(21)
				self._efficiency_histograms_csvorder[efficiency_name].SetMarkerSize(1)
				self._efficiency_histograms_csvorder[efficiency_name].SetMarkerColor(self._colors[efficiency_name])
				self._efficiency_histograms_csvorder[efficiency_name].SetLineWidth(2)
				self._efficiency_histograms_csvorder[efficiency_name].SetLineStyle(self._line_styles[efficiency_name])
				self._efficiency_histograms_csvorder[efficiency_name].SetLineColor(self._colors[efficiency_name])

				self._mjj_histograms_csvorder[hist_pair[0]].SetMarkerStyle(25)
				self._mjj_histograms_csvorder[hist_pair[0]].SetMarkerColor(self._colors[efficiency_name])
				self._mjj_histograms_csvorder[hist_pair[1]].SetMarkerStyle(20)
				self._mjj_histograms_csvorder[hist_pair[1]].SetMarkerColor(self._colors[efficiency_name])

			if self._mjj_histograms_vetothirdjet.has_key(hist_pair[0]) and self._mjj_histograms_vetothirdjet.has_key(hist_pair[1]):
				self._efficiency_histograms_vetothirdjet[efficiency_name] = self._mjj_histograms_vetothirdjet[hist_pair[0]].Clone()
				self._efficiency_histograms_vetothirdjet[efficiency_name].SetName("h_efficiency_" + efficiency_name + "_vetothirdjet")
				self._efficiency_histograms_vetothirdjet[efficiency_name].SetDirectory(0)
				if "JetHT" in efficiency_name:
					self._efficiency_histograms_vetothirdjet[efficiency_name].Divide(self._mjj_histograms_vetothirdjet[hist_pair[1]])
				else:
					self._efficiency_histograms_vetothirdjet[efficiency_name].Divide(self._mjj_histograms_vetothirdjet[hist_pair[0]], self._mjj_histograms_vetothirdjet[hist_pair[1]], 1, 1, "B")
				self._efficiency_histograms_vetothirdjet[efficiency_name].SetMarkerStyle(22)
				self._efficiency_histograms_vetothirdjet[efficiency_name].SetMarkerSize(1)
				self._efficiency_histograms_vetothirdjet[efficiency_name].SetMarkerColor(self._colors[efficiency_name])
				self._efficiency_histograms_vetothirdjet[efficiency_name].SetLineWidth(2)
				self._efficiency_histograms_vetothirdjet[efficiency_name].SetLineStyle(self._line_styles[efficiency_name])
				self._efficiency_histograms_vetothirdjet[efficiency_name].SetLineColor(self._colors[efficiency_name])

				self._mjj_histograms_vetothirdjet[hist_pair[0]].SetMarkerStyle(26)
				self._mjj_histograms_vetothirdjet[hist_pair[0]].SetMarkerColor(self._colors[efficiency_name])
				self._mjj_histograms_vetothirdjet[hist_pair[1]].SetMarkerStyle(20)
				self._mjj_histograms_vetothirdjet[hist_pair[1]].SetMarkerColor(self._colors[efficiency_name])

	def FitEfficiencies(self):
		self._efficiency_fits = {}
		self._efficiency_fits2 = {} # Extra fit for systematic
		self._fit_ranges = {}
		self._fit_chi2s = {}
		self._fit_ndfs = {}

		# Fit JetHT with constant first, to get online b-tag efficiencies
		# - Use coarse histograms.
		for name, hist in self._efficiency_histograms.iteritems():
			if "JetHT" in name:
				print "Fitting " + name					
				if "highmass" in name:
					self._fit_ranges[name] = [526, 1607]
				elif "lowmass" in name:
					self._fit_ranges[name] = [386, 1058]
				#self._efficiency_fits[name] = ROOT.TF1("linear_" + name, "[0]+[1]*x", self._fit_ranges[name][0], self._fit_ranges[name][1])
				self._efficiency_fits[name] = ROOT.TF1("constant_" + name, "[0]", self._fit_ranges[name][0], self._fit_ranges[name][1])
				#self._efficiency_fits[name].SetLineColor(self._colors[name])
				self._efficiency_fits[name].SetLineColor(kRed)
				if "lowmass" in name:
					self._efficiency_fits[name].SetParameter(0, 0.2)
					#self._efficiency_fits[name].SetParameter(1, 0.)
				elif "highmass" in name:
					self._efficiency_fits[name].SetParameter(0, 0.5)
					#self._efficiency_fits[name].SetParameter(1, 0.)
				hist.Fit(self._efficiency_fits[name], "R0")
				chi2 = self._efficiency_fits[name].GetChisquare()
				ndf = self._efficiency_fits[name].GetNDF()
				prob = TMath.Prob(chi2, ndf)
				print "chi2/ndf = {}/{}, p = {}".format(chi2, ndf, prob)
				# Calculate RSS as well
				bin_min = hist.GetXaxis().FindBin(self._fit_ranges[name][0]+1.e-5)
				bin_max = hist.GetXaxis().FindBin(self._fit_ranges[name][1]-1.e-5)
				rss_constant = 0
				for bin in xrange(bin_min, bin_max + 1):
					bin_center = hist.GetXaxis().GetBinCenter(bin)
					rss_constant += (hist.GetBinContent(bin) - self._efficiency_fits[name].Eval(bin_center))**2
				print "RSS = {}".format(rss_constant)

				print "\nDoing a linear fit for systematic"
				self._efficiency_fits2[name] = ROOT.TF1("linear_" + name, "[0]+[1]*x", self._fit_ranges[name][0], self._fit_ranges[name][1])
				if "lowmass" in name:
					self._efficiency_fits2[name].SetParameter(0, 0.2)
					self._efficiency_fits2[name].SetParameter(1, 0.)
				elif "highmass" in name:
					self._efficiency_fits2[name].SetParameter(0, 0.5)
					self._efficiency_fits2[name].SetParameter(1, 0.)
				hist.Fit(self._efficiency_fits2[name], "R0")
				chi2_linear = self._efficiency_fits2[name].GetChisquare()
				ndf_linear = self._efficiency_fits2[name].GetNDF()
				prob = TMath.Prob(chi2_linear, ndf_linear)
				print "chi2/ndf = {}/{}, p = {}".format(chi2_linear, ndf_linear, prob)
				# Calculate RSS as well
				rss_linear = 0
				for bin in xrange(bin_min, bin_max + 1):
					bin_center = hist.GetXaxis().GetBinCenter(bin)
					rss_linear += (hist.GetBinContent(bin) - self._efficiency_fits2[name].Eval(bin_center))**2
				print "RSS = {}".format(rss_linear)
				f21 = (rss_constant-rss_linear) / (rss_linear / (hist.GetNbinsX() - 2))
				cl21 = 1. - TMath.FDistI(f21, 2 - 1, hist.GetNbinsX() - 2)
				print "F21 = {}".format(f21)
				print "CS21 = {}".format(cl21)

		# Fit SingleMu and BJetPlusX with sigmoid functions
		for name, hist in self._efficiency_histograms_fine.iteritems():
			if "SingleMu" in name or "BJet" in name:
				print "Fitting " + name
				if "highmass" in name:
					self._fit_ranges[name] = [330, 605]
				elif "llowmass" in name:
					self._fit_ranges[name] = [175, 420]
				elif "lowmass" in name:
					self._fit_ranges[name] = [175, 420]

				self._efficiency_fits[name] = ROOT.TF1("sigmoid_" + name, "[2] * (1. / (1. + TMath::Exp(-1. * (x - [0]) / [1])))", self._fit_ranges[name][0], self._fit_ranges[name][1])
				#self._efficiency_fits[name].SetLineColor(self._colors[name])
				self._efficiency_fits[name].SetLineColor(kRed)
				self._efficiency_fits[name].SetParameter(0, 200.)
				if "lowmass" in name:
					self._efficiency_fits[name].SetParLimits(0, 0., 400.)
				elif "highmass" in name:
					self._efficiency_fits[name].SetParLimits(0, 0., 1000.)
				self._efficiency_fits[name].SetParameter(1, 50.)
				#self._efficiency_fits[name].SetParameter(2, 1.)
				if "SingleMu" in name:
					self._efficiency_fits[name].SetParameter(2, 0.17)
				elif name == "BJet60_53_lowmass":
					self._efficiency_fits[name].SetParameter(2, 1.)
				elif name == "BJet80_70_highmass" or name == "BJet60_53_highmass":
					self._efficiency_fits[name].SetParameter(2, 3.)
	
				fit_result = hist.Fit(self._efficiency_fits[name], "R0S")
				self._fit_chi2s[name] = self._efficiency_fits[name].GetChisquare()
				self._fit_ndfs[name] = self._efficiency_fits[name].GetNDF()
				fit_result.Print("V")

	def MakeSingleMuComparisons(self):
		for sr_name in ["lowmass", "llowmass"]:
			c = TCanvas("c_trigeff_SingleMu_comparison_" + sr_name)
			l = TLegend(0.7, 0.7, 0.88, 0.88)
			l.SetFillColor(0)
			l.SetBorderSize(0)
			frame = TH1F("frame_" + sr_name, "frame_" + sr_name, 100, 0., 1800.)
			frame.SetMinimum(0.)
			frame.SetMaximum(self._efficiency_guesses["SingleMu_" + sr_name] * 3.)
			frame.GetXaxis().SetTitle("m_{jj} [GeV]")
			frame.GetYaxis().SetTitle("Efficiency")
			frame.Draw()
			for mu_name in ["SingleMu", "SingleMu24i", "SingleMu40"]:
				self._efficiency_histograms[mu_name + "_" + sr_name].Draw("same")
				l.AddEntry(self._efficiency_histograms[mu_name + "_" + sr_name], mu_name, "lp")
				if self._efficiency_histograms_csvorder.has_key(mu_name + "_" + sr_name):
					self._efficiency_histograms_csvorder[mu_name + "_" + sr_name].Draw("same")
					l.AddEntry(self._efficiency_histograms_csvorder[mu_name + "_" + sr_name], mu_name + "/CSV", "lp")
				if self._efficiency_histograms_vetothirdjet.has_key(mu_name + "_" + sr_name):
					self._efficiency_histograms_vetothirdjet[mu_name + "_" + sr_name].Draw("same")
					l.AddEntry(self._efficiency_histograms_vetothirdjet[mu_name + "_" + sr_name], mu_name + "/<3j", "lp")
			l.Draw()
			c.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/TriggerEfficiency/figures/" + c.GetName() + ".pdf")

	def MakeSingleEfficiencyPlots(self):
		for efficiency_name, efficiency_histogram in self._efficiency_histograms.iteritems():
			print "[MakeSingleEfficiencyPlots] INFO : Plotting " + efficiency_name
			c = TCanvas("c_trigeff_" + efficiency_name, "c_trigeff_" + efficiency_name, 800, 1000)
			l = TLegend(0.5, 0.6, 0.88, 0.88)
			l.SetFillColor(0)
			l.SetBorderSize(0)
			l.SetHeader(efficiency_name)

			top = TPad("top", "top", 0., 0.5, 1., 1.)
			top.SetBottomMargin(0.02)
			top.Draw()
			top.SetLogy()
			top.cd()
			if "SingleMu" in efficiency_name:
				frame_top = TH1F("frame_top", "frame_top", 100, 0., 1000.)
			elif "BJet80_70_highmass" in efficiency_name:
				frame_top = TH1F("frame_top", "frame_top", 100, 100., 800.)
			else:
				frame_top = TH1F("frame_top", "frame_top", 100, 0., 1800.)
			frame_top.GetXaxis().SetTitleSize(0)
			frame_top.GetXaxis().SetLabelSize(0)
			frame_top.GetYaxis().SetTitle("Events / 1 GeV")
			frame_top.SetMinimum(0.01)
			frame_top.SetMaximum(10. * max(self._mjj_histograms[self._efficiency_combinations[efficiency_name][0]].GetMaximum(), self._mjj_histograms[self._efficiency_combinations[efficiency_name][1]].GetMaximum()))
			frame_top.Draw()
			self._mjj_histograms[self._efficiency_combinations[efficiency_name][0]].SetMarkerColor(seaborn.GetColorRoot("default", 0))
			self._mjj_histograms[self._efficiency_combinations[efficiency_name][0]].SetLineColor(seaborn.GetColorRoot("default", 0))
			self._mjj_histograms[self._efficiency_combinations[efficiency_name][0]].Draw("same")
			self._mjj_histograms[self._efficiency_combinations[efficiency_name][1]].SetMarkerColor(seaborn.GetColorRoot("default", 2))
			self._mjj_histograms[self._efficiency_combinations[efficiency_name][1]].SetLineColor(seaborn.GetColorRoot("default", 2))
			self._mjj_histograms[self._efficiency_combinations[efficiency_name][1]].Draw("same")
			if self._efficiency_combinations[efficiency_name][0] in self._legend_entries:
				l.AddEntry(self._mjj_histograms[self._efficiency_combinations[efficiency_name][0]], self._legend_entries[self._efficiency_combinations[efficiency_name][0]], "lp")
			else:
				l.AddEntry(self._mjj_histograms[self._efficiency_combinations[efficiency_name][0]], self._efficiency_combinations[efficiency_name][0], "lp")
			if self._efficiency_combinations[efficiency_name][1] in self._legend_entries:
				l.AddEntry(self._mjj_histograms[self._efficiency_combinations[efficiency_name][1]], self._legend_entries[self._efficiency_combinations[efficiency_name][1]], "lp")
			else:
				l.AddEntry(self._mjj_histograms[self._efficiency_combinations[efficiency_name][1]], self._efficiency_combinations[efficiency_name][1], "lp")
			l.Draw()

			c.cd()
			bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
			bottom.SetTopMargin(0.02)
			bottom.SetBottomMargin(0.2)
			bottom.Draw()
			bottom.cd()
			if "SingleMu" in efficiency_name:
				frame_bottom = TH1F("frame_bottom", "frame_bottom", 100, 0., 1000.)
			elif "BJet80_70_highmass" in efficiency_name:
				frame_bottom = TH1F("frame_bottom", "frame_bottom", 100, 100., 800.)
			else:
				frame_bottom = TH1F("frame_bottom", "frame_bottom", 100, 0., 1800.)
			frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
			frame_bottom.GetYaxis().SetTitle("Efficiency")
			frame_bottom.SetMinimum(0.)
			frame_bottom.SetMaximum(self._efficiency_guesses[efficiency_name] * 1.8)
			frame_bottom.Draw()
			efficiency_histogram.SetMarkerColor(kBlack)
			efficiency_histogram.SetLineColor(kBlack)
			efficiency_histogram.SetLineStyle(1)
			efficiency_histogram.SetLineWidth(1)
			efficiency_histogram.Draw("same")
			#if "SingleMu" in efficiency_name or "BJet" in efficiency_name:
			self._efficiency_fits[efficiency_name].Draw("same")
			if efficiency_name in self._efficiency_fits2:
				self._efficiency_fits2[efficiency_name].SetLineStyle(3)
				self._efficiency_fits2[efficiency_name].SetLineWidth(2)
				self._efficiency_fits2[efficiency_name].SetLineColor(seaborn.GetColorRoot("dark", 2))
				self._efficiency_fits2[efficiency_name].Draw("same")
			c.cd()
			c.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/TriggerEfficiency/figures/" + c.GetName() + ".pdf")
			ROOT.SetOwnership(c, False)
			ROOT.SetOwnership(top, False)
			ROOT.SetOwnership(bottom, False)

			# Efficiency only
			c_eff = TCanvas("c_trigeffonly_" + efficiency_name, "c_trigeffonly_" + efficiency_name, 800, 600)
			frame_bottom.Draw()
			efficiency_histogram.Draw("same")
			if "JetHT" in efficiency_name:
				self._efficiency_fits[efficiency_name].SetLineColor(kRed)
				self._efficiency_fits[efficiency_name].Draw("same")
				if efficiency_name in self._efficiency_fits2:
					self._efficiency_fits2[efficiency_name].Draw("same")
			c_eff.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/TriggerEfficiency/figures/" + c_eff.GetName() + ".pdf")


	def MakeSingleEfficiencyFinePlots(self):
		for efficiency_name, efficiency_histogram in self._efficiency_histograms_fine.iteritems():
			print "[MakeSingleEfficiencyPlots] INFO : Plotting " + efficiency_name
			c = TCanvas("c_trigeff_" + efficiency_name + "_fine", "c_trigeff_" + efficiency_name + "_fine", 800, 1200)
			l = TLegend(0.5, 0.6, 0.88, 0.88)
			l.SetFillColor(0)
			l.SetBorderSize(0)
			l.SetHeader(efficiency_name)

			top = TPad("top", "top", 0., 0.7, 1., 1.)
			top.SetBottomMargin(0.03)
			top.Draw()
			top.SetLogy()
			top.cd()
			if "SingleMu" in efficiency_name:
				frame_top = TH1F("frame_top", "frame_top", 100, 0., 800.)
			elif "BJet80_70_highmass" in efficiency_name:
				frame_top = TH1F("frame_top", "frame_top", 100, 100., 800.)
			else:
				frame_top = TH1F("frame_top", "frame_top", 100, 0., 1800.)
			frame_top.GetXaxis().SetTitleSize(0)
			frame_top.GetXaxis().SetLabelSize(0)
			frame_top.GetYaxis().SetLabelSize(0.06)
			frame_top.GetYaxis().SetTitleSize(0.06)
			frame_top.GetYaxis().SetTitleOffset(0.8)
			frame_top.GetYaxis().SetTitle("Events")
			frame_top.SetMinimum(0.5)
			frame_top.SetMaximum(10. * max(self._mjj_histograms_fine[self._efficiency_combinations[efficiency_name][0]].GetMaximum(), self._mjj_histograms_fine[self._efficiency_combinations[efficiency_name][1]].GetMaximum()))
			frame_top.Draw()
			self._mjj_histograms_fine[self._efficiency_combinations[efficiency_name][0]].SetMarkerColor(seaborn.GetColorRoot("default", 0))
			self._mjj_histograms_fine[self._efficiency_combinations[efficiency_name][0]].SetLineColor(seaborn.GetColorRoot("default", 0))
			self._mjj_histograms_fine[self._efficiency_combinations[efficiency_name][0]].Draw("same")
			self._mjj_histograms_fine[self._efficiency_combinations[efficiency_name][1]].SetMarkerColor(seaborn.GetColorRoot("default", 2))
			self._mjj_histograms_fine[self._efficiency_combinations[efficiency_name][1]].SetLineColor(seaborn.GetColorRoot("default", 2))
			self._mjj_histograms_fine[self._efficiency_combinations[efficiency_name][1]].Draw("same")
			if self._efficiency_combinations[efficiency_name][0] in self._legend_entries:
				l.AddEntry(self._mjj_histograms_fine[self._efficiency_combinations[efficiency_name][0]], self._legend_entries[self._efficiency_combinations[efficiency_name][0]], "lp")
			else:
				l.AddEntry(self._mjj_histograms_fine[self._efficiency_combinations[efficiency_name][0]], self._efficiency_combinations[efficiency_name][0], "lp")
			if self._efficiency_combinations[efficiency_name][1] in self._legend_entries:
				l.AddEntry(self._mjj_histograms_fine[self._efficiency_combinations[efficiency_name][1]], self._legend_entries[self._efficiency_combinations[efficiency_name][1]], "lp")
			else:
				l.AddEntry(self._mjj_histograms_fine[self._efficiency_combinations[efficiency_name][1]], self._efficiency_combinations[efficiency_name][1], "lp")
			l.Draw()

			c.cd()
			middle = TPad("middle", "middle", 0., 0.4, 1., 0.7)
			middle.SetTopMargin(0.03)
			middle.SetBottomMargin(0.03)
			middle.Draw()
			middle.cd()
			if "SingleMu" in efficiency_name:
				frame_middle = TH1F("frame_middle", "frame_middle", 100, 0., 800.)
			elif "BJet80_70_highmass" in efficiency_name:
				frame_middle = TH1F("frame_middle", "frame_middle", 100, 100., 800.)
			else:
				frame_middle = TH1F("frame_middle", "frame_middle", 100, 0., 1800.)
			#frame_middle.GetXaxis().SetTitle("m_{jj} [GeV]")
			frame_middle.GetXaxis().SetLabelSize(0.)
			frame_middle.GetXaxis().SetTitleSize(0.)
			frame_middle.GetYaxis().SetLabelSize(0.06)
			frame_middle.GetYaxis().SetTitleSize(0.06)
			frame_middle.GetYaxis().SetTitleOffset(0.8)
			frame_middle.GetYaxis().SetTitle("Efficiency")
			frame_middle.SetMinimum(0.)
			frame_middle.SetMaximum(self._efficiency_guesses[efficiency_name] * 1.5)
			frame_middle.Draw()
			efficiency_histogram.SetMarkerColor(kBlack)
			efficiency_histogram.SetLineColor(kBlack)
			efficiency_histogram.SetLineStyle(1)
			efficiency_histogram.SetLineWidth(1)
			efficiency_histogram.Draw("same")
			#if "SingleMu" in efficiency_name or "BJet" in efficiency_name:
			self._efficiency_fits[efficiency_name].Draw("same")

			c.cd()
			bottom = TPad("bottom", "bottom", 0., 0., 1., 0.4)
			bottom.SetTopMargin(0.02)
			bottom.SetBottomMargin(0.25)
			bottom.Draw()
			bottom.cd()
			if "SingleMu" in efficiency_name:
				frame_bottom = TH1F("frame_bottom", "frame_bottom", 100, 0., 800.)
			elif "BJet80_70_highmass" in efficiency_name:
				frame_bottom = TH1F("frame_bottom", "frame_bottom", 100, 100., 800.)
			else:
				frame_bottom = TH1F("frame_bottom", "frame_bottom", 100, 0., 1800.)
			frame_bottom.GetXaxis().SetTitleSize(0.05)
			frame_bottom.GetXaxis().SetLabelSize(0.05)
			frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
			frame_bottom.GetYaxis().SetLabelSize(0.06 * 0.3 / 0.4)
			frame_bottom.GetYaxis().SetTitleSize(0.06 * 0.3 / 0.4)
			frame_bottom.GetYaxis().SetTitleOffset(0.8)
			frame_bottom.GetYaxis().SetTitle("(eff - fit) / #sigma(eff)")
			frame_bottom.SetMinimum(-3.)
			frame_bottom.SetMaximum(3.)
			frame_bottom.Draw()
			fit_ratio_histogram = efficiency_histogram.Clone()
			for bin in xrange(1, fit_ratio_histogram.GetNbinsX() + 1):
				this_eff = efficiency_histogram.GetBinContent(bin)
				this_deff = efficiency_histogram.GetBinError(bin)
				this_mjj = efficiency_histogram.GetXaxis().GetBinCenter(bin)
				if this_mjj < self._efficiency_fits[efficiency_name].GetXmin() or this_mjj > self._efficiency_fits[efficiency_name].GetXmax() or this_deff == 0.:
					fit_ratio_histogram.SetBinContent(bin, 0.)
					fit_ratio_histogram.SetBinError(bin, 0.)
				else:
					this_fit = self._efficiency_fits[efficiency_name].Eval(this_mjj)
					fit_ratio_histogram.SetBinContent(bin, (this_eff - this_fit) / this_deff)
					fit_ratio_histogram.SetBinError(bin, 0.)
			fit_ratio_histogram.SetMarkerStyle(20)
			fit_ratio_histogram.SetMarkerColor(kBlack)
			fit_ratio_histogram.SetMarkerSize(1)
			fit_ratio_histogram.SetFillStyle(1001)
			fit_ratio_histogram.SetFillColor(seaborn.GetColorRoot("default", 2))
			fit_ratio_histogram.Draw("hist same")

			c.cd()
			c.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/TriggerEfficiency/figures/" + c.GetName() + ".pdf")
			ROOT.SetOwnership(c, False)
			ROOT.SetOwnership(top, False)
			ROOT.SetOwnership(middle, False)

	def MakeJetHTSingleMuComparisons(self):
		for sr_name in ["lowmass"]:
			c = TCanvas("c_trigeff_jetht_vs_singlemu_" + sr_name, "c_trigeff_jetht_vs_singlemu_" + sr_name, 800, 600)
			l = TLegend(0.7, 0.7, 0.88, 0.88)
			if sr_name == "highmass":
				l.SetHeader("HLT 160/120+btags")
			elif sr_name == "lowmass":
				l.SetHeader("Low mass trigger #varepsilon_{b-tag}")
			elif sr_name == "llowmass":
				l.SetHeader("HLT 60/53+btags")
			l.SetFillColor(0)
			l.SetBorderSize(0)
			if sr_name == "lowmass":
				frame = TH1F("frame_" + sr_name, "frame_" + sr_name, 100, 0., 1200.)
			else:
				frame = TH1F("frame_" + sr_name, "frame_" + sr_name, 100, 0., 1800.)
			frame.SetMinimum(0.)
			frame.SetMaximum(0.5)
			#frame.SetMaximum(self._efficiency_guesses["SingleMu24i_" + sr_name] * 10.)
			frame.GetXaxis().SetTitle("m_{jj} [GeV]")
			frame.GetYaxis().SetTitle("Efficiency")
			frame.Draw()

			# SingleMu-based efficiency
			if sr_name == "lowmass":
				singlemu_fit = TF1("singlemu_eff_" + sr_name, "[0]", 296, 1058)
				self._efficiency_histograms["SingleMu24i_lowmass"].Fit(singlemu_fit, "R0")
				self._efficiency_histograms["SingleMu24i_lowmass"].SetMarkerColor(seaborn.GetColorRoot("default", 2))
				self._efficiency_histograms["SingleMu24i_lowmass"].SetLineColor(seaborn.GetColorRoot("default", 2))
				self._efficiency_histograms["SingleMu24i_lowmass"].SetMarkerStyle(20)
				singlemu_fit.SetLineColor(seaborn.GetColorRoot("default", 2))
				self._efficiency_histograms["SingleMu24i_lowmass"].Draw("same")
			elif sr_name == "highmass":
				singlemu_fit = TF1("singlemu_eff_" + sr_name, "[0]", 526, 1607)
			singlemu_fit.Draw("same")
			l.AddEntry(self._efficiency_histograms["SingleMu24i_" + sr_name], "Mu24i", "lp")
			l.AddEntry(singlemu_fit, "Mu24i fit", "lp")

			# JetHT_based efficiency
			if sr_name == "lowmass":
				jetht_fit = TF1("jetht_eff_" + sr_name, "[0]", 386, 1058)
			elif sr_name == "highmass":
				jetht_fit = TF1("jetht_eff_" + sr_name, "[0]", 526, 1607)
			self._efficiency_histograms["JetHT_" + sr_name].Fit(jetht_fit, "R0")
			self._efficiency_histograms["JetHT_" + sr_name].SetLineColor(seaborn.GetColorRoot("default", 0))
			self._efficiency_histograms["JetHT_" + sr_name].SetMarkerStyle(24)
			self._efficiency_histograms["JetHT_" + sr_name].SetMarkerColor(seaborn.GetColorRoot("default", 0))
			self._efficiency_histograms["JetHT_" + sr_name].Draw("same")
			jetht_fit.SetLineColor(seaborn.GetColorRoot("default", 0))
			jetht_fit.Draw("same")
			l.AddEntry(self._efficiency_histograms["JetHT_" + sr_name], "JetHT", "lp")
			l.AddEntry(jetht_fit, "JetHT fit", "lp")
			l.Draw()
			c.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/TriggerEfficiency/figures/" + c.GetName() + ".pdf")

	def MakeSingleMu6053Comparison(self):
		c = TCanvas("c_trigeff_80_70_over_60_53", "c_trigeff_80_70_over_60_53", 800, 600)
		l = TLegend(0.7, 0.7, 0.88, 0.88)
		l.SetHeader("80/70 vs 60/53")
		l.SetFillColor(0)
		l.SetBorderSize(0)
		frame = TH1F("frame_SingleMu6053Comparison", "frame_SingleMu6053Comparison", 100, 0., 1800.)
		frame.SetMinimum(0.)
		frame.SetMaximum(1.5)
		frame.GetXaxis().SetTitle("m_{jj} [GeV]")
		frame.GetYaxis().SetTitle("Efficiency")
		frame.Draw()

		# SingleMu-based efficiency
		self._efficiency_histograms["BJet60_53_lowmass"].Draw("same")
		l.AddEntry(self._efficiency_histograms["BJet60_53_lowmass"], "60/53 + b-tag", "lp")

		# JetHT_based efficiency
		self._efficiency_histograms["SingleMu_lowmass"].Draw("same")
		l.AddEntry(self._efficiency_histograms["SingleMu_lowmass"], "SingleMu", "lp")
		l.Draw()
		c.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/TriggerEfficiency/figures/" + c.GetName() + ".pdf")

if __name__ == "__main__":
	trig_eff_plotter = TrigEffPlotter()
	#trig_eff_plotter.MakeSingleMuComparisons()
	trig_eff_plotter.MakeSingleEfficiencyPlots()
	trig_eff_plotter.MakeSingleEfficiencyFinePlots()
	trig_eff_plotter.MakeJetHTSingleMuComparisons()
	#trig_eff_plotter.MakeSingleMu6053Comparison()
