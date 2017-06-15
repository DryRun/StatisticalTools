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

masses = {"trigbbl_CSVTM":[350,400,500,600,750],"trigbbh_CSVTM":[600,750,900,1200]}
for analysis in ["trigbbl_CSVTM", "trigbbh_CSVTM"]:
	tg_eff = {}
	max_eff = -1.
	for model in ["Hbb", "ZPrime", "RSG"]:
		if analysis == "trigbbl_CSVTM":
			notrig_analysis = "NoTrigger_eta1p7_CSVTM"
		elif analysis == "trigbbh_CSVTM":
			notrig_analysis = "NoTrigger_eta2p2_CSVTM"

        #signal_pdf_file = analysis_config.get_signal_fit_file(notrig_analysis, args.model, mass, "bukin", interpolated=(not mass in analysis_config.simulation.simulated_masses))
		if "bbl" in analysis:
			mjj_range = [296, 1058]
		elif "bbh" in analysis:
			mjj_range = [526, 1607]

		tg_eff[model] = TGraphErrors(len(masses[analysis]))
		for i, mass in enumerate(masses[analysis]):
			f = TFile(analysis_config.get_b_histogram_filename(notrig_analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")))
			if not f.IsOpen():
				continue
			low_bin = f.Get("BHistograms/h_pfjet_mjj").GetXaxis().FindBin(mjj_range[0] + 1.e-5)
			high_bin = f.Get("BHistograms/h_pfjet_mjj").GetXaxis().FindBin(mjj_range[1] - 1.e-5)
			numerator = f.Get("BHistograms/h_pfjet_mjj").Integral(low_bin, high_bin)
			denominator = f.Get("BHistograms/h_sample_nevents").Integral()
			if denominator > 0:
				eff = numerator/denominator * trigger_efficiency.online_btag_eff[analysis][0]
				deff = (eff*(1.-eff)/denominator)**0.5 * trigger_efficiency.online_btag_eff[analysis][0]
				eff *= 100
				deff *= 100
				if eff + deff > max_eff:
					max_eff = eff + deff
				tg_eff[model].SetPoint(i, mass, eff)
				tg_eff[model].SetPointError(i, 0., deff)
				print "Mass = " + str(mass) + " / eff = " + str(eff) + " +/- " + str(deff)
	tg_eff["Hbb"].SetMarkerStyle(24)
	tg_eff["Hbb"].SetMarkerColor(seaborn.GetColorRoot("dark", 0))
	tg_eff["Hbb"].SetMarkerSize(1)
	tg_eff["Hbb"].SetLineStyle(1)
	tg_eff["Hbb"].SetLineColor(seaborn.GetColorRoot("dark", 0))
	tg_eff["Hbb"].SetLineWidth(1)
	tg_eff["Hbb"].SetFillStyle(3002)
	tg_eff["Hbb"].SetFillColor(seaborn.GetColorRoot("pastel", 0))

	tg_eff["ZPrime"].SetMarkerStyle(25)
	tg_eff["ZPrime"].SetMarkerColor(seaborn.GetColorRoot("dark", 1))
	tg_eff["ZPrime"].SetMarkerSize(1)
	tg_eff["ZPrime"].SetLineStyle(1)
	tg_eff["ZPrime"].SetLineColor(seaborn.GetColorRoot("dark", 1))
	tg_eff["ZPrime"].SetLineWidth(1)
	tg_eff["ZPrime"].SetFillStyle(3002)
	tg_eff["ZPrime"].SetFillColor(seaborn.GetColorRoot("pastel", 1))

	tg_eff["RSG"].SetMarkerStyle(26)
	tg_eff["RSG"].SetMarkerColor(seaborn.GetColorRoot("dark", 2))
	tg_eff["RSG"].SetMarkerSize(1)
	tg_eff["RSG"].SetLineStyle(1)
	tg_eff["RSG"].SetLineColor(seaborn.GetColorRoot("dark", 2))
	tg_eff["RSG"].SetLineWidth(1)
	tg_eff["RSG"].SetFillStyle(3002)
	tg_eff["RSG"].SetFillColor(seaborn.GetColorRoot("pastel", 2))

	c = TCanvas("c_signal_efficiency_{}".format(analysis), "c_signal_efficiency_{}".format(analysis), 800, 600)
	frame = TH1D("frame", "frame", 100, min(masses[analysis]) - 100, max(masses[analysis]) + 200)
	frame.SetMinimum(0.)
	frame.SetMaximum(max_eff * 1.2)
	frame.GetXaxis().SetTitle("Resonance Mass [GeV]")
	frame.GetYaxis().SetTitle("Selection Efficiency [%]")
	frame.Draw("axis")

	tg_eff["Hbb"].Draw("lp3")
	tg_eff["ZPrime"].Draw("lp3")
	tg_eff["RSG"].Draw("lp3")

	l = TLegend(0.68, 0.68, 0.88, 0.88)
	l.SetFillColor(0)
	l.SetBorderSize(0)
	l.AddEntry(tg_eff["Hbb"], "Scalar", "lpf3")
	l.AddEntry(tg_eff["ZPrime"], "Z'", "lpf")
	l.AddEntry(tg_eff["RSG"], "RS Graviton", "lpf")
	l.Draw()

	Root.CMSLabel(0.15, 0.92, "Internal", 1, 0.5); 

	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")
	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".eps")
	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".C")

