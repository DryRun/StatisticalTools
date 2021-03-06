import os
import sys
import ROOT
from ROOT import *
from math import sqrt
from array import array

from CMSDIJET.StatisticalTools.roofit_functions_huge import *
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
seaborn = Root.SeabornInterface()
seaborn.Initialize()

ftest_dir = "/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/FTest"
dijet_binning = array("d", [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5000]) #5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8000])

def mlfit(fit_function, analysis, condor=True):
	directory = "{}/{}/{}".format(ftest_dir, analysis, fit_function)
	os.system("mkdir -pv {}".format(directory))
	os.chdir(directory)

	# Load some stuff from default workspace
	f_w_in = TFile(limit_config.get_workspace_filename(analysis, "Hbb", 600, fitBonly=False, fitTrigger=False, correctTrigger=True, useMCTrigger=False, qcd=False, fitOffB=False), "READ")
	w_in = f_w_in.Get("w")
	data_rdh = w_in.data("data_obs")
	data_integral = data_rdh.sum(False)
	signal_pdf = w_in.pdf("signal")
	signal_pdf.SetName("signal_rbp")
	mjj = w_in.var("mjj")

	# Convert signal to histogram PDF for performance (without systematics, no need for full Bukin)
	signal_hist = signal_pdf.createHistogram("mjj", int(mjj.getMax()-mjj.getMin()))
	signal_rdh = RooDataHist("signal_rdh", "signal_rdh", RooArgList(mjj), signal_hist)
	signal_rhp = RooHistPdf("signal", "signal", RooArgSet(mjj), signal_rdh)

	# Make background pdf
	background_pdf, background_parameters = make_background_pdf(fit_function, mjj, collision_energy=8000.)
	set_initial_values(background_parameters, fit_function, analysis)
	background_norm = RooRealVar('background_' + fit_function + '_norm', 'background_' + fit_function + '_norm', data_integral, 0., 1.e8)
	background_epdf = RooExtendPdf('ebackground_' + fit_function, 'ebackground_' + fit_function, background_pdf, background_norm)

	# Make S+B model
	signal_norm = RooRealVar('signal_norm', 'signal_norm', 0., 0., 1e+08)
	signal_epdf = RooExtendPdf("esignal", "esignal", signal_rhp, signal_norm)
	#model = RooAddPdf("model_" + fit_function, "s+b", RooArgList(background_epdf, signal_epdf))
	model = background_epdf

	# Note: in the first iteration, you did an initial fit to set the initial values.
	# Now that you have a first iteration, you are using those results as the starting point.
	# This way, the whole family starts from the same place.
	signal_norm.setConstant()
	fit_result = model.fitTo(data_rdh, RooFit.Save(True), RooFit.Extended(True), RooFit.Verbose(0))
	signal_norm.setConstant(False)

	# Make new workspace
	w_out = RooWorkspace('w','workspace')
	getattr(w_out,'import')(data_rdh, ROOT.RooCmdArg(), RooFit.Rename("data_obs"))
	getattr(w_out,'import')(background_pdf, ROOT.RooCmdArg(), RooFit.Rename("background_" + fit_function), RooFit.RecycleConflictNodes())
	getattr(w_out,'import')(background_norm,ROOT.RooCmdArg(),RooFit.Rename("background_" + fit_function + "_norm"))
	getattr(w_out,'import')(signal_rhp, ROOT.RooCmdArg(), RooFit.Rename("signal"), RooFit.RecycleConflictNodes())
	getattr(w_out,'import')(signal_norm, ROOT.RooCmdArg())
	w_out.writeToFile(directory + "/workspace.root")


	# Make datacard
	datacard = open(directory + "/datacard.txt", "w")
	datacard.write('imax 1\n')
	datacard.write('jmax 1\n')
	datacard.write('kmax *\n')
	datacard.write('---------------\n')
	datacard.write('shapes * * workspace.root w:$PROCESS\n')
	datacard.write('---------------\n')
	datacard.write('bin 1\n')
	datacard.write('observation -1\n')
	datacard.write('------------------------------\n')
	datacard.write('bin          1          1\n')
	datacard.write('process      signal     background_' + fit_function + '\n')
	datacard.write('process      0          1\n')
	datacard.write('rate         1         1\n')
	datacard.write('------------------------------\n')
	# Background fit parameters --- flat prior
	datacard.write('background_' + fit_function + '_norm  flatParam\n')
	for par_name, par in background_parameters.iteritems():
		datacard.write(fit_function + "_" + par_name + '  flatParam\n')
	datacard.close()

	# Make script for running combine
	combine_script = open("run_mlfit.sh", "w")
	combine_script.write("#!/bin/bash\n")
	combine_script.write("mkdir -pv plots\n")
	combine_script.write("combine datacard.txt -M MaxLikelihoodFit --saveWithUncertainties --saveShapes --plots -v 2 -t -1 2>&1 | tee combine_log.txt\n")
	combine_script.close()

	if condor:
		print "Running mlfit on condor"
		os.chdir(directory)
		condor_script = open(directory + "/condor_submit.sh", "w")
		condor_script.write("#!/bin/bash\n")
		condor_script.write("csub run_mlfit.sh --cmssw --no_retar -F {}/workspace.root,{}/datacard.txt\n".format(directory, directory))
		condor_script.close()
		os.system("source {}/condor_submit.sh".format(directory))
	else:
		print "Running mlfit locally"
		print "JK I dont' want to"
		#os.system("source run_mlfit.sh")

# npars = list of #parameters in family (e.g. [4,5,6])
## npar_function = {npar:[functions], npar:[functions], ...}
def ftest(npars, npar_functions, analysis):
	print "Ftest for " + analysis
	top_directory = "{}/{}".format(ftest_dir, analysis)

	# Loop over all functions, and plot/calculate stuff for each
	rss = {}
	chi2 = {}
	ndf = {}
	ad = {}
	ks = {}
	nll = {}
	first = True
	for npar, functions in npar_functions.iteritems():
		for function in functions:
			if first:
				data_file = TFile(top_directory + "/" + function + "/workspace.root", "READ")
				data_rdh = data_file.Get("w").data("data_obs")
				mjj = data_file.Get("w").var("mjj")
				data_hist = data_rdh.createHistogram("mjj", int(mjj.getMax()-mjj.getMin()))
				data_hist.SetDirectory(0)
				data_file.Close()
				first = False
			fit_file = TFile(top_directory + "/" + function + "/mlfit.root", "READ")
			fit_hist = fit_file.Get("shapes_fit_b/bin1/background_{}".format(function))
			if not fit_hist:
				print "[ftest] ERROR : Couldn't find histogram {} in file {}. Abandoning.".format("shapes_fit_b/bin1/background_{}".format(function), top_directory + "/" + function + "/mlfit.root")
				rss[function] = None
				chi2[function] = None
				ndf[function] = None
				continue
			fit_hist.SetDirectory(0)
			for bin in xrange(1, fit_hist.GetNbinsX() + 1):
				fit_hist.SetBinError(bin, 0)

			rss[function] = calculate_rss(fit_hist, data_hist)
			chi2[function] = calculate_chi2(fit_hist, data_hist)
			ndf[function] = data_hist.GetNbinsX() - npar
			ad[function] = calculate_andersondarling(fit_hist, data_hist)
			ks[function] = calculate_ks(fit_hist, data_hist)

			nll[function] = fit_file.Get("fit_b").minNll()

			fit_file.Close()

	print "Printing F-test results for this family:"
	print "\\begin{table}"
	print "\t\\begin{tabular}{|c|c|c|c|c|}"
	print "\t\t\\hline\n"
	print "\t\t$f_1$ & $f_2$ & $F_{21}$ & CL$_{21}$ & Choose $f_2$ \\\\\n\t\t\\hline"

	for i in xrange(len(npars) - 1):
		n1 = npars[i]
		n2 = npars[i+1]
		if not n2 > n1:
			print "[ftest] ERROR : npar must be in ascending order."
			sys.exit(1)
		for f1 in npar_functions[n1]:
			if rss[f1] == None:
				continue
			for f2 in npar_functions[n2]:
				if rss[f2] == None:
					continue
				f21 = ((rss[f1]-rss[f2]) / (n2 - n1)) / (rss[f2] / (data_hist.GetNbinsX() - n2))
				cl21 = 1. - TMath.FDistI(f21, n2 - n1, data_hist.GetNbinsX() - n2)
				if cl21 < 0.05:
					sigsymbol = "$\\bullet$"
				else:
					sigsymbol = "$\\circ$"
				print "\t\t{}\t&\t{}\t&\t{:.2f}\t&\t{:.2e}\t&\t{}\t\\\\\n\t\t\\hline".format(f1, f2, f21, cl21, sigsymbol)
	print "\t\\end{tabular}"
	print "\t\\caption{}"
	print "\t\\label{table:ftest-family}"
	print "\\end{table}\n"


def qof_table(fit_functions, analysis):
	top_directory = "{}/{}".format(ftest_dir, analysis)

	# Loop over all functions, and plot/calculate stuff for each
	rss = {}
	chi2 = {}
	ndf = {}
	ad = {}
	ks = {}
	nll = {}
	first = True
	for f in fit_functions:
		if first:
			data_file = TFile(top_directory + "/" + f + "/workspace.root", "READ")
			data_rdh = data_file.Get("w").data("data_obs")
			mjj = data_file.Get("w").var("mjj")
			data_hist = data_rdh.createHistogram("mjj", int(mjj.getMax()-mjj.getMin()))
			data_hist.SetDirectory(0)
			data_file.Close()
			first = False
		fit_file = TFile(top_directory + "/" + f + "/mlfit.root", "READ")
		fit_hist = fit_file.Get("shapes_fit_b/bin1/background_{}".format(f))
		if not fit_hist:
			print "[ftest] ERROR : Couldn't find histogram {} in file {}. Abandoning.".format("shapes_fit_b/bin1/background_{}".format(f), top_directory + "/" + f + "/mlfit.root")
			rss[f] = None
			chi2[f] = None
			ndf[f] = None
			continue
		fit_hist.SetDirectory(0)
		for bin in xrange(1, fit_hist.GetNbinsX() + 1):
			fit_hist.SetBinError(bin, 0)

		rss[f] = calculate_rss(fit_hist, data_hist)
		chi2[f] = calculate_chi2(fit_hist, data_hist)
		ndf[f] = data_hist.GetNbinsX() - npar
		ad[f] = calculate_andersondarling(fit_hist, data_hist)
		ks[f] = calculate_ks(fit_hist, data_hist)

		nll[f] = fit_file.Get("fit_b").minNll()

		fit_file.Close()

	print "\nPrinting quality of fit results:"
	print "\\begin{table}"
	print "\t\\begin{tabular}{|c|c|c|c|}" # c|
	print "\t\t\\hline"
	print "\t\tFunction & Min NLL & $\\chi^2/$NDF ($p$) & KS test statistic ($p$) \\\\\n\t\t\\hline" # & AD test statistic ($p$)
	for f in fit_functions:
			print "\t\t{} & {:.2e} & {:.2f}/{} ({:.2f}) & {:.2e} ({:.2f}) \\\\\n\t\t\\hline\n".format(f, nll[f], chi2[f], ndf[f], TMath.Prob(chi2[f], ndf[f]), ks[f][1], ks[f][0]) # & {:.2f} ({:.2f}) || , ad[f][1], ad[f][0]
	print "\t\\end{tabular}"
	print "\t\\caption{}"
	print "\t\\label{table:fit-quality-family}"
	print "\\end{table}\n"

def make_data_fit_plots(family_fit_functions, analysis):
	top_directory = "{}/{}".format(ftest_dir, analysis)

	# Loop over all functions, and plot/calculate stuff for each
	rss = {}
	chi2 = {}
	ndf = {}
	for family, fit_functions in family_fit_functions.iteritems():
		first = True
		fit_hists = {}
		for fit_function in fit_functions:
			if first:
				data_file = TFile(top_directory + "/" + fit_function + "/workspace.root", "READ")
				data_rdh = data_file.Get("w").data("data_obs")
				mjj = data_file.Get("w").var("mjj")
				data_hist = data_rdh.createHistogram("mjj", int(mjj.getMax()-mjj.getMin()))
				data_hist.SetDirectory(0)
				data_file.Close()
				first = False
			fit_file = TFile(top_directory + "/" + fit_function + "/mlfit.root", "READ")
			if not fit_file.IsOpen():
				print "[make_data_fit_plots] WARNING : Couldn't open file {}. Skipping.".format(top_directory + "/" + fit_function + "/mlfit.root")
				continue
			fit_hists[fit_function] = fit_file.Get("shapes_fit_b/bin1/background_{}".format(fit_function))
			if not fit_hists[fit_function]:
				print "[make_data_fit_plots] WARNING : Couldn't find histogram {} in file {}. Skipping.".format("shapes_fit_b/bin1/background_{}".format(function), top_directory + "/" + function + "/mlfit.root")
				continue
			fit_hists[fit_function].SetDirectory(0)
			for bin in xrange(1, fit_hists[fit_function].GetNbinsX() + 1):
				fit_hists[fit_function].SetBinError(bin, 0)
			fit_file.Close()

			chi2[fit_function] = calculate_chi2(fit_hists[fit_function], data_hist)
			ndf[fit_function] = int(data_hist.GetNbinsX() - npar)

		make_data_fit_plot(data_hist, fit_hists, family, analysis)
		this_binning = array("d", [])
		for bin_value in dijet_binning:
			if bin_value >= data_hist.GetXaxis().GetXmin() and bin_value <= data_hist.GetXaxis().GetXmax():
				this_binning.append(bin_value)
		make_data_fit_plot(data_hist, fit_hists, family, analysis + "_coarse", binning=this_binning)

def make_data_fit_plot(data_hist, fit_hists, fit_family, analysis, binning=None):
	save_tag = fit_family + "_" + analysis
	if binning:
		data_hist = data_hist.Rebin(len(binning) - 1, data_hist.GetName() + "_coarse", binning)
		for bin in xrange(1, data_hist.GetNbinsX() + 1):
			bin_width = data_hist.GetXaxis().GetBinWidth(bin)
			data_hist.SetBinContent(bin, data_hist.GetBinContent(bin) / bin_width)
			data_hist.SetBinError(bin, data_hist.GetBinError(bin) / bin_width)
		for fit_function in fit_hists.keys():
			fit_hists[fit_function] = fit_hists[fit_function].Rebin(len(binning) - 1, fit_hists[fit_function].GetName() + "_coarse", binning)
			for bin in xrange(1, data_hist.GetNbinsX() + 1):
				bin_width = data_hist.GetXaxis().GetBinWidth(bin)
				fit_hists[fit_function].SetBinContent(bin, fit_hists[fit_function].GetBinContent(bin) / bin_width)
				fit_hists[fit_function].SetBinError(bin, fit_hists[fit_function].GetBinError(bin) / bin_width)

	c = TCanvas("c_data_vs_fit_{}".format(save_tag), "c_{}".format(save_tag), 800, 1000)
	top = TPad("top", "top", 0., 0.5, 1., 1.)
	top.SetBottomMargin(0.02)
	top.SetLogy()
	top.Draw()
	top.cd()

	x_min = data_hist.GetXaxis().GetXmin()
	x_max = data_hist.GetXaxis().GetXmax()
	initial_xrange = x_max - x_min
	x_max += initial_xrange * 0.15
	x_min -= initial_xrange * 0.15
	frame_top = TH1D("frame_top", "frame_top", 100, x_min, x_max)
	frame_top.SetMinimum(0.1)
	frame_top.SetMaximum(max(data_hist.GetMaximum(), max([x.GetMaximum() for x in fit_hists.values()])) * 10.)
	frame_top.GetXaxis().SetLabelSize(0)
	frame_top.GetXaxis().SetTitleSize(0)
	frame_top.GetYaxis().SetTitle("Events / GeV")
	frame_top.Draw("axis")

	for style_counter, fit_function in enumerate(sorted(fit_hists.keys())):
		fit_hists[fit_function].SetMarkerStyle(20)
		fit_hists[fit_function].SetMarkerSize(0)
		fit_hists[fit_function].SetLineStyle(style_counter % 4 + 1)
		fit_hists[fit_function].SetLineWidth(2)
		fit_hists[fit_function].SetLineColor(seaborn.GetColorRoot("default", style_counter))
		fit_hists[fit_function].Draw("hist same")

	data_hist.SetMarkerStyle(20)
	data_hist.SetMarkerSize(1)
	data_hist.SetMarkerColor(1)
	data_hist.Draw("p same")

	frame_top.Draw("axis same")

	l = TLegend(0.7, 0.7, 0.88, 0.88)
	l.SetFillColor(0)
	l.SetBorderSize(0)
	l.AddEntry(data_hist, "Data", "p")
	for style_counter, fit_function in enumerate(sorted(fit_hists.keys())):
		l.AddEntry(fit_hists[fit_function], fit_function, "l")
	l.Draw()

	print analysis
	if "trigbbl_CSVTM" in analysis:
		sr_text = "Low mass SR"
	else:
		sr_text = "High mass SR"
	Root.myText(0.3, 0.24, 1, sr_text, 0.5)
	Root.myText(0.3, 0.20, 1, fit_family, 0.5)
	#chi2_text = "#chi^{{2}}/NDF = {:.2f}/{}, p={:.2f}".format(chi2, ndf, TMath.Prob(chi2, ndf))
	#Root.myText(0.3, 0.16, 1, chi2_text, 0.5)

	c.cd()
	bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
	bottom.SetTopMargin(0.02)
	bottom.SetBottomMargin(0.2)
	bottom.Draw()
	bottom.cd()

	frame_bottom = TH1D("frame_bottom", "frame_bottom", 100, x_min, x_max)
	frame_bottom.SetMinimum(-3.)
	frame_bottom.SetMaximum(3.)
	frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
	frame_bottom.GetYaxis().SetTitle("#frac{data - fit}{#sqrt{fit}}")
	frame_bottom.Draw("axis")

	zero = TLine(x_min, 0., x_max, 0.)

	pull_hists = {}
	for style_counter, fit_function in enumerate(sorted(fit_hists.keys())):
		pull_hists[fit_function] = data_hist.Clone()
		pull_hists[fit_function].Reset()
		for bin in xrange(1, data_hist.GetNbinsX() + 1):
			if fit_hists[fit_function].GetBinContent(bin) > 0:
				if binning:
					bin_width = data_hist.GetXaxis().GetBinWidth(bin)
					pull_hists[fit_function].SetBinContent(bin, sqrt(bin_width) * (data_hist.GetBinContent(bin) - fit_hists[fit_function].GetBinContent(bin)) / sqrt(fit_hists[fit_function].GetBinContent(bin)))

				else:
					pull_hists[fit_function].SetBinContent(bin, (data_hist.GetBinContent(bin) - fit_hists[fit_function].GetBinContent(bin)) / sqrt(fit_hists[fit_function].GetBinContent(bin)))
		pull_hists[fit_function].SetMarkerStyle(20)
		pull_hists[fit_function].SetMarkerSize(0)
		pull_hists[fit_function].SetLineColor(seaborn.GetColorRoot("default", style_counter))
		pull_hists[fit_function].SetLineWidth(2)
		pull_hists[fit_function].SetLineStyle(style_counter % 4 + 1)
		#pull_hists[fit_function].SetFillStyle(1001)
		#pull_hists[fit_function].SetFillColor(seaborn.GetColorRoot("default", style_counter))
		pull_hists[fit_function].Draw("hist same ][")
	zero = TLine(x_min, 0., x_max, 0.)
	zero.SetLineStyle(1)
	zero.SetLineWidth(2)
	zero.SetLineColor(1)
	zero.Draw("same")

	c.cd()
	c.SaveAs("{}/figures/{}.pdf".format(ftest_dir, c.GetName()))

	ROOT.SetOwnership(c, False)
	ROOT.SetOwnership(top, False)
	ROOT.SetOwnership(bottom, False)

# Helper functions
def calculate_rss(hist1, hist2):
	if hist1.GetNbinsX() != hist2.GetNbinsX():
		print "[calculate_rss] ERROR : Histogram bin mismatch"
		print "{}=>{}, {}=>{}".format(hist1.GetName(), hist1.GetNbinsX(), hist2.GetName(), hist2.GetNbinsX())
		sys.exit(1)
	rss = 0
	for bin in xrange(1, hist1.GetNbinsX() + 1):
		rss += (hist1.GetBinContent(bin) - hist2.GetBinContent(bin))**2
	return rss

def calculate_chi2(model_hist, data_hist):
	if model_hist.GetNbinsX() != data_hist.GetNbinsX():
		print "[calculate_rss] ERROR : Histogram bin mismatch"
		print "{}=>{}, {}=>{}".format(model_hist.GetName(), model_hist.GetNbinsX(), data_hist.GetName(), data_hist.GetNbinsX())
		sys.exit(1)
	chi2 = 0
	for bin in xrange(1, model_hist.GetNbinsX() + 1):
		if model_hist.GetBinContent(bin) > 0:
			fit = model_hist.GetBinContent(bin)
			data = data_hist.GetBinContent(bin)
			chi2 += ((fit - data) / sqrt(fit))**2
	return chi2

def calculate_andersondarling(model_hist, data_hist):
	ad_prob = data_hist.AndersonDarlingTest(model_hist)
	ad_ts = data_hist.AndersonDarlingTest(model_hist, "T")
	return (ad_prob, ad_ts)

def calculate_ks(model_hist, data_hist):
	ks_prob = data_hist.KolmogorovTest(model_hist)
	ks_ts = data_hist.KolmogorovTest(model_hist, "M")
	return (ks_prob, ks_ts)

if __name__ == "__main__":
	from argparse import ArgumentParser
	parser = ArgumentParser()
	parser.add_argument("--analyses", type=str, default="trigbbl_CSVTM,trigbbh_CSVTM", help="Analysis name")
	functions_group = parser.add_mutually_exclusive_group() 
	functions_group.add_argument("--all", action="store_true", help="Run all fit functions")
	functions_group.add_argument("--family", type=str, help="Run single fit family (see roofit_functions_huge)")
	functions_group.add_argument("--functions", type=str, help="Specify specific functions (mlfit only, ftest requires a family)")
	action_group = parser.add_mutually_exclusive_group()
	action_group.add_argument("--mlfit", action="store_true", help="Run fits")
	action_group.add_argument("--ftest", action="store_true", help="Run f-test (have to run mlfit first)")
	action_group.add_argument("--qof", action="store_true", help="Make quality of fit table")
	action_group.add_argument("--plots", action="store_true", help="Make data-fit plots")
	args = parser.parse_args()

	analyses = args.analyses.split(",")

	family_npars = {
		"dijet":[3,4,5,6],
		"dijet4_OffB":[4,5,6],
		"rational":[3,4],
		"modexp":[3,4,5],
		"polyx":[5,6,7],
		"atlas":[4,5,6],
		"polypower":[4,5],
		"onethird":[3],
	}

	# Maps for keeping track of the number of parameters in each function
	fit_functions_npar = {}
	npar_fit_functions = {}
	family_npar_fit_functions = {}
	if args.all:
		fit_families = ["dijet", "polyx", "atlas", "polypower", "modexp", "dijet4_OffB"]
		fit_functions = []
		family_fit_functions = {}
		for fit_family in fit_families:
			family_fit_functions[fit_family] = []
			family_npar_fit_functions[fit_family] = {}
			for npar in family_npars[fit_family]:
				fit_functions.append(fit_family + str(npar))
				fit_functions_npar[fit_family + str(npar)] = npar
				if not npar in npar_fit_functions:
					npar_fit_functions[npar] = []
				if not npar in family_npar_fit_functions[fit_family]:
					family_npar_fit_functions[fit_family][npar] = []
				npar_fit_functions[npar].append(fit_family + str(npar))
				family_npar_fit_functions[fit_family][npar].append(fit_family + str(npar))
				family_fit_functions[fit_family].append(fit_family + str(npar))
	elif args.family:
		fit_families = [args.family]
		fit_functions = []
		family_fit_functions = {}
		for fit_family in fit_families:
			family_fit_functions[fit_family] = []
			family_npar_fit_functions[fit_family] = {}
			for npar in family_npars[fit_family]:
				fit_functions.append(fit_family + str(npar))
				fit_functions_npar[fit_family + str(npar)] = npar
				if not npar in npar_fit_functions:
					npar_fit_functions[npar] = []
				if not npar in family_npar_fit_functions[fit_family]:
					family_npar_fit_functions[fit_family][npar] = []
				npar_fit_functions[npar].append(fit_family + str(npar))
				family_npar_fit_functions[fit_family][npar].append(fit_family + str(npar))
				family_fit_functions[fit_family].append(fit_family + str(npar))
	elif args.functions:
		if args.ftest or args.plots:
			print "[do_f_tests] ERROR : Cannot specify args.functions in ftest or plots mode"
			sys.exit(1)
		fit_families = []
		fit_functions = args.functions.split(",")

	if args.mlfit:
		for analysis in analyses:
			for fit_function in fit_functions:
				mlfit(fit_function, analysis, condor=True)

	if args.ftest:
		for analysis in analyses:
			for fit_family in fit_families:
				ftest(family_npars[fit_family], family_npar_fit_functions[fit_family], analysis)

	if args.qof:
		for analysis in analyses:
			qof_table(fit_functions, analysis)


	if args.plots:
		for analysis in analyses:
			make_data_fit_plots(family_fit_functions, analysis)
