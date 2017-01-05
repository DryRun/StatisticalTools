#!/usr/bin/env python

import sys, os, re
from joblib import Parallel, delayed
from argparse import ArgumentParser
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/CMSDIJET/QCDAnalysis")
import analysis_configuration_8TeV as analysis_config
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
from argparse import ArgumentParser
import subprocess

import ROOT
from ROOT import *
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_7_4_15/python/MyTools/RootUtils")
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
seaborn = Root.SeabornInterface()
seaborn.Initialize()


def RunAsimovCheck(datacard, expected_signal, dry_run=False):
	name = os.path.basename(datacard).replace("datacard", "").replace(".txt", "") + "_" + str(expected_signal)
	command = "combine {} -M MaxLikelihoodFit -t -1 --expectSignal {} --name {}".format(datacard, expected_signal, name)
	if dry_run:
		print command
	else:
		command = "combine {} -M MaxLikelihoodFit -t -1 --expectSignal {} --name {}".format(datacard, expected_signal, name)
		logfile = "/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/AsimovFits/asimovcheck" + name + ".log"
		os.system(command + " >& " + logfile)

		#with open("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/Logs/asimovcheck" + name + ".log", "w") as out:
		#	subprocess.call("combine {} -M MaxLikelihoodFit -t -1 --expectSignal {}".format(datacard, expected_signal), stdout=out)
		
def AsimovFitTables(model, analysis, fit_functions, expected_signal, masses):
	import ROOT
	from ROOT import TFile, TTree, TLeaf
	table_file = open("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/AsimovFits/table_{}_{}_{}.tex".format(model, analysis, expected_signal), 'w')
	table_file.write("\\begin{table}\n")
	table_file.write("\t\\centering\n")
	table_file.write("\t\\begin{tabular}{|c|")
	for i in xrange(len(fit_functions)):
		table_file.write("c|")
	table_file.write("}\n")
	table_file.write("\t\t\\hline\n")
	table_file.write("\t\tMass\t")
	for fit_function in fit_functions:
		table_file.write("\t&\t$\\mu_{" + fit_function + "}$")
	table_file.write("\\\\\n\t\t\\hline\n")
	for mass in masses:
		table_file.write("\t\t{} GeV ".format(mass))
		for fit_function in fit_functions:
			tree_file = "/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/AsimovFits/mlfit_{}_{}_{}_{}_fitSignal_correctTrigger_{}.root".format(analysis, model, mass, fit_function, expected_signal)
			f = TFile(tree_file, "READ")
			if not f.IsOpen():
				print "[AsimovFitTables] WARNING : File {} was not found. Skipping.".format(tree_file)
			t = f.Get("tree_fit_sb")
			if not t:
				print "[AsimovFitTables] ERROR : Didn't find tree_fit_sb in file {}. Printing contents of file.".format(tree_file)
				f.Print()
				sys.exit(1)
			t.GetEntry(0)
			mu = t.GetLeaf("mu").GetValue()
			mu_err = t.GetLeaf("muErr").GetValue()
			mu_sig = (mu / mu_err if mu_err > 0 else -1.)
			#table_file.write("\t\t{} GeV & {} & ${:.2}$ & ${:.2}$ & ${:.2}$ \n".format(mass, t.GetLeaf("fit_status").GetValue(), mu, mu_err, mu_sig))
			table_file.write("\t&\t ${:.2} \pm {:.2}$".format(mu, mu_err))
			#" & ${:.2}$ & ${:.2}$ & ${:.2}$ \n".format(mass, t.GetLeaf("fit_status").GetValue(), mu, mu_err, mu_sig))
			f.Close()
		table_file.write("\t\\\\\n\t\t\\hline\n")
	table_file.write("\t\\end{tabular}\n")
	table_file.write("\\end{table}\n")
	table_file.close()

def AsimovFitPlots(model, analysis, fit_functions, expected_signal, masses):
	graph_mu = {}
	graph_pull = {}
	canvas_mu = TCanvas("c_asimov_fit_mu_{}_{}_{}".format(model, analysis, expected_signal), "c_asimov_fit_mu_{}_{}_{}".format(model, analysis, expected_signal), 800, 600)
	legend_mu = TLegend(0.83, 0.7, 0.89, 0.88)
	legend_mu.SetFillColor(0)
	legend_mu.SetBorderSize(0)
	frame_mu = TH1D("frame_mu", "frame_mu", 100, min(masses) - 50., max(masses) + 150.)
	frame_mu.GetXaxis().SetTitle("Signal Mass [GeV]")
	frame_mu.GetYaxis().SetTitle("#mu")
	frame_mu.SetMinimum(expected_signal -1.)
	frame_mu.SetMaximum(expected_signal + 1.)
	frame_mu.Draw("axis")

	canvas_pull = TCanvas("c_asimov_fit_pull_{}_{}_{}".format(model, analysis, expected_signal), "c_asimov_fit_pull_{}_{}_{}".format(model, analysis, expected_signal), 800, 600)
	legend_pull = TLegend(0.83, 0.7, 0.89, 0.88)
	legend_pull.SetFillColor(0)
	legend_pull.SetBorderSize(0)
	frame_pull = TH1D("frame_pull", "frame_pull", 100, min(masses) - 50., max(masses) + 150.)
	frame_pull.GetXaxis().SetTitle("Signal Mass [GeV]")
	frame_pull.GetYaxis().SetTitle("#mu/#sigma_{#mu}")
	frame_pull.SetMinimum(-3.)
	frame_pull.SetMaximum(3.)
	frame_pull.Draw("axis")

	style_counter = 0
	for fit_function in fit_functions:
		graph_mu[fit_function] = TGraphErrors(len(masses))
		graph_pull[fit_function] = TGraph(len(masses))
		for i in xrange(len(masses)):
			tree_file = "/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/AsimovFits/mlfit_{}_{}_{}_{}_fitSignal_correctTrigger_{}.root".format(analysis, model, masses[i], fit_function, expected_signal)
			f = TFile(tree_file, "READ")
			if not f.IsOpen():
				print "[AsimovFitTables] WARNING : File {} was not found. Skipping.".format(tree_file)
			t = f.Get("tree_fit_sb")
			if not t:
				print "[AsimovFitTables] ERROR : Didn't find tree_fit_sb in file {}. Printing contents of file.".format(tree_file)
				f.Print()
				sys.exit(1)
			t.GetEntry(0)
			mu = t.GetLeaf("mu").GetValue()
			mu_err = t.GetLeaf("muErr").GetValue()
			mu_sig = ((mu - expected_signal) / mu_err if mu_err > 0 else -3.)
			graph_mu[fit_function].SetPoint(i, masses[i], mu)
			graph_mu[fit_function].SetPointError(i, 0., mu_err)
			graph_pull[fit_function].SetPoint(i, masses[i], mu_sig)
		canvas_mu.cd()
		graph_mu[fit_function].SetMarkerStyle(20 + style_counter)
		graph_mu[fit_function].SetMarkerColor(seaborn.GetColorRoot("default", style_counter))
		graph_mu[fit_function].Draw("p")
		legend_mu.AddEntry(graph_mu[fit_function], fit_function, "p")

		canvas_pull.cd()
		graph_pull[fit_function].SetMarkerStyle(20 + style_counter)
		graph_pull[fit_function].SetMarkerColor(seaborn.GetColorRoot("default", style_counter))
		graph_pull[fit_function].Draw("p")
		legend_pull.AddEntry(graph_pull[fit_function], fit_function, "p")

		style_counter += 1
	canvas_mu.cd()
	line_expected = TLine(min(masses) - 50., expected_signal, max(masses) + 150., expected_signal)
	line_expected.SetLineColor(kGray)
	line_expected.SetLineStyle(2)
	line_expected.Draw("same")
	legend_mu.Draw()
	canvas_mu.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/AsimovFits/{}.pdf".format(canvas_mu.GetName()))

	canvas_pull.cd()
	legend_pull.Draw()
	line_zero = TLine(min(masses) - 50., 0., max(masses) + 150., 0.)
	line_zero.SetLineColor(kGray)
	line_zero.SetLineStyle(2)
	line_zero.Draw("same")
	canvas_pull.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/AsimovFits/{}.pdf".format(canvas_pull.GetName()))



if __name__ == "__main__":
	parser = ArgumentParser(description='Script that runs limit calculation for specified mass points')
	parser.add_argument("--models", type=str, default="Hbb,RSG", help="Models")
	parser.add_argument("--analyses", type=str, default="trigbbl_CSVTM,trigbbh_CSVTM", help="Analyses")
	parser.add_argument("--masses", type=str, default=None, help="Masses")
	parser.add_argument("--fit_functions", type=str, default="f3,f4", help="Fit functions")
	parser.add_argument("--expected_signals", type=str, default="0,1", help="Injected signal strength")
	parser.add_argument("--fit", action="store_true", help="Run fits")
	parser.add_argument("--table", action="store_true", help="Make table from results")
	parser.add_argument("--plot", action="store_true", help="Make plot from results")
	args = parser.parse_args()

	models = args.models.split(",")
	analyses = args.analyses.split(",")
	fit_functions = args.fit_functions.split(",")
	expected_signals = [float(x) for x in args.expected_signals.split(",")]
	masses = {}
	if args.masses:
		for analysis in analyses:
			masses[analysis] = [int(x) for x in args.masses.split(",")]
	else:
		for analysis in analyses:
			if "bbl" in analysis:
				masses[analysis] = xrange(400, 650, 50)
			elif "bbh" in analysis:
				masses[analysis] = xrange(600, 1250, 50)

	if args.fit:
		cwd = os.getcwd()
		os.chdir("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/AsimovFits")
		for analysis in analyses:
			#for model in models:
			#	for mass in masses[analysis]:
			#		for fit_function in fit_functions:
			#			for expected_signal in expected_signals:
			#				datacard = limit_config.get_datacard_filename(analysis, model, mass, fit_function, #fitSignal=True, correctTrigger=True)
			#				name = os.path.basename(datacard).replace("datacard", "").replace(".txt", "") + "_" + #str(expected_signal)
			#				command = "combine {} -M MaxLikelihoodFit -t -1 --expectSignal {} --name {}".format(#datacard, expected_signal, name)
			#				log = "/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/Logs/asimovcheck" + name + "#.log"
			#				os.system(command + " >& " + log)
			Parallel(n_jobs=4)(delayed(RunAsimovCheck)(limit_config.get_datacard_filename(analysis, model, mass, fit_function, fitSignal=True, correctTrigger=True), expected_signal) for mass in masses[analysis] for model in models for fit_function in fit_functions for expected_signal in expected_signals)

		os.chdir(cwd)
	if args.table:
		for analysis in analyses:
			Parallel(n_jobs=4)(delayed(AsimovFitTables)(model, analysis, fit_functions, expected_signal, masses[analysis]) for model in models for expected_signal in expected_signals)
	if args.plot:
		for analysis in analyses:
			Parallel(n_jobs=4)(delayed(AsimovFitPlots)(model, analysis, fit_functions, expected_signal, masses[analysis]) for model in models for expected_signal in expected_signals)
