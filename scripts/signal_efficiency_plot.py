import os
import sys
import ROOT
from ROOT import *
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/src/CMSDIJET/QCDAnalysis/python/")
import analysis_configuration_8TeV as analysis_config
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
seaborn = Root.SeabornInterface()
seaborn.Initialize()
import CMSDIJET.StatisticalTools.trigger_efficiency as trigger_efficiency

masses = {"trigbbl_CSVTM":[325,350,400,500,600,750],"trigbbh_CSVTM":[600,750,900,1200]}
tg_eff = {}
for analysis_style, analysis in enumerate(["trigbbl_CSVTM", "trigbbh_CSVTM"]):
	tg_eff[analysis] = {}
	max_eff = -1.
	for model in ["Hbb", "ZPrime", "RSG"]:
		#if analysis == "trigbbl_CSVTM":
		#	notrig_analysis = "NoTrigger_eta1p7_CSVTM"
		#elif analysis == "trigbbh_CSVTM":
		#	notrig_analysis = "NoTrigger_eta2p2_CSVTM"

        #signal_pdf_file = analysis_config.get_signal_fit_file(notrig_analysis, args.model, mass, "bukin", interpolated=(not mass in analysis_config.simulation.simulated_masses))
		if "bbl" in analysis:
			mjj_range = [296, 1058]
		elif "bbh" in analysis:
			mjj_range = [526, 1607]

		tg_eff[analysis][model] = TGraphErrors(len(masses[analysis]))
		for i, mass in enumerate(masses[analysis]):
			f = TFile(analysis_config.get_b_histogram_filename(analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")))
			if not f.IsOpen():
				continue
			low_bin = f.Get("BHistograms/h_pfjet_mjj").GetXaxis().FindBin(mjj_range[0] + 1.e-5)
			high_bin = f.Get("BHistograms/h_pfjet_mjj").GetXaxis().FindBin(mjj_range[1] - 1.e-5)
			numerator = f.Get("BHistograms/h_pfjet_mjj").Integral(low_bin, high_bin)
			if model == "ZPrime":
				denominator = f.Get("BHistograms/h_input_nevents").Integral()
				print "{} / {} den = {}".format(model, analysis, denominator)
			else:
				denominator = f.Get("BHistograms/h_sample_nevents").Integral()
			if denominator > 0:
				eff = numerator/denominator
				deff = (eff*(1.-eff)/denominator)**0.5
				eff *= 100
				deff *= 100
				if eff + deff > max_eff:
					max_eff = eff + deff
				tg_eff[analysis][model].SetPoint(i, mass, eff)
				tg_eff[analysis][model].SetPointError(i, 0., deff)
				print "Model = " + model + " / mass = " + str(mass) + " / eff = " + str(eff) + " +/- " + str(deff)
	tg_eff[analysis]["Hbb"].SetMarkerStyle(24)
	tg_eff[analysis]["Hbb"].SetMarkerColor(seaborn.GetColorRoot("dark", 0))
	tg_eff[analysis]["Hbb"].SetMarkerSize(1)
	tg_eff[analysis]["Hbb"].SetLineStyle(1 + analysis_style)
	tg_eff[analysis]["Hbb"].SetLineColor(seaborn.GetColorRoot("dark", 0))
	tg_eff[analysis]["Hbb"].SetLineWidth(1)
	tg_eff[analysis]["Hbb"].SetFillStyle(3002)
	tg_eff[analysis]["Hbb"].SetFillColor(seaborn.GetColorRoot("pastel", 0))

	tg_eff[analysis]["ZPrime"].SetMarkerStyle(25)
	tg_eff[analysis]["ZPrime"].SetMarkerColor(seaborn.GetColorRoot("dark", 1))
	tg_eff[analysis]["ZPrime"].SetMarkerSize(1)
	tg_eff[analysis]["ZPrime"].SetLineStyle(1 + analysis_style)
	tg_eff[analysis]["ZPrime"].SetLineColor(seaborn.GetColorRoot("dark", 1))
	tg_eff[analysis]["ZPrime"].SetLineWidth(1)
	tg_eff[analysis]["ZPrime"].SetFillStyle(3002)
	tg_eff[analysis]["ZPrime"].SetFillColor(seaborn.GetColorRoot("pastel", 1))

	tg_eff[analysis]["RSG"].SetMarkerStyle(26)
	tg_eff[analysis]["RSG"].SetMarkerColor(seaborn.GetColorRoot("dark", 2))
	tg_eff[analysis]["RSG"].SetMarkerSize(1)
	tg_eff[analysis]["RSG"].SetLineStyle(1 + analysis_style)
	tg_eff[analysis]["RSG"].SetLineColor(seaborn.GetColorRoot("dark", 2))
	tg_eff[analysis]["RSG"].SetLineWidth(1)
	tg_eff[analysis]["RSG"].SetFillStyle(3002)
	tg_eff[analysis]["RSG"].SetFillColor(seaborn.GetColorRoot("pastel", 2))

c = TCanvas("c_signal_efficiency_combined", "c_signal_efficiency_combined", 800, 600)
c.SetBottomMargin(0.12)
c.SetRightMargin(0.075)
frame = TH1D("frame", "frame", 100, min(masses["trigbbl_CSVTM"]) - 100, max(masses["trigbbh_CSVTM"]) + 150)
frame.SetMinimum(0.)
frame.SetMaximum(max_eff * 1.5)
frame.GetXaxis().SetTitle("Resonance Mass [GeV]")
frame.GetYaxis().SetTitle("Acceptance #times efficiency [%]")
frame.GetXaxis().SetTitleSize(0.05)
frame.GetYaxis().SetTitleSize(0.05)
frame.GetXaxis().SetTitleOffset(1.06)
frame.GetYaxis().SetTitleOffset(1.04)
frame.GetXaxis().SetLabelSize(0.05)
frame.GetYaxis().SetLabelSize(0.05)
frame.GetXaxis().SetLabelOffset(0.012)
frame.GetYaxis().SetLabelOffset(0.007)
frame.Draw("axis")

tg_eff["trigbbl_CSVTM"]["Hbb"].Draw("lp3")
tg_eff["trigbbl_CSVTM"]["ZPrime"].Draw("lp3")
tg_eff["trigbbl_CSVTM"]["RSG"].Draw("lp3")
tg_eff["trigbbh_CSVTM"]["Hbb"].Draw("lp3")
tg_eff["trigbbh_CSVTM"]["ZPrime"].Draw("lp3")
tg_eff["trigbbh_CSVTM"]["RSG"].Draw("lp3")

l = TLegend(0.68, 0.5, 0.88, 0.88)
l.SetFillColor(0)
l.SetBorderSize(0)
l.AddEntry(tg_eff["trigbbl_CSVTM"]["Hbb"], "Scalar (SR1)", "lpf3")
l.AddEntry(tg_eff["trigbbh_CSVTM"]["Hbb"], "Scalar (SR2)", "lpf3")
l.AddEntry(tg_eff["trigbbl_CSVTM"]["ZPrime"], "Z' (SR1)", "lpf")
l.AddEntry(tg_eff["trigbbh_CSVTM"]["ZPrime"], "Z' (SR2)", "lpf")
l.AddEntry(tg_eff["trigbbl_CSVTM"]["RSG"], "RS Graviton (SR1)", "lpf")
l.AddEntry(tg_eff["trigbbh_CSVTM"]["RSG"], "RS Graviton (SR2)", "lpf")
l.Draw()

Root.CMSLabel(0.15, 0.83, "Simulation", 1, 0.6); 

c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")
c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".eps")
c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".C")

