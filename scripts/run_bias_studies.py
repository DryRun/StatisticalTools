import os
import sys
import array
import ROOT
from ROOT import *
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/src/CMSDIJET/QCDAnalysis/python")
import analysis_configuration_8TeV as analysis_config
import CMSDIJET.StatisticalTools.limit_configuration as limit_config

gROOT.SetBatch(True)
#gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPalette(1)
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
seaborn = Root.SeabornInterface()
seaborn.Initialize()

# Light temperature palette
stops = array.array('d', [0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000])
red = array.array('d', [  31./255.,  71./255., 123./255., 160./255., 210./255., 222./255., 214./255., 199./255., 183./255.])
green = array.array('d', [  40./255., 117./255., 171./255., 211./255., 231./255., 220./255., 190./255., 132./255.,  65./255.])
blue = array.array('d', [ 234./255., 214./255., 228./255., 222./255., 210./255., 160./255., 105./255.,  60./255.,  34./255.])
palette_index = TColor.CreateGradientColorTable(9, stops, red, green, blue, 255, 1.)
gStyle.SetNumberContours(255)
gStyle.SetPaintTextFormat('1.2f')

datacard_folders = {
	"Hbb":{
		"trigbbl_CSVTM":"/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/BiasStudies/davidlGcards",
		"trigbbh_CSVTM":"/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/BiasStudies/davidhGcards",
	},
	"RSG":{
		"trigbbl_CSVTM":"/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/BiasStudies/davidlRcards",
		"trigbbh_CSVTM":"/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/BiasStudies/davidhRcards",
	},
}
#functions = ["f" + str(x) for x in xrange(1, 6)]

top_directory = "/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/BiasStudies/"

def run_bias_study(name, gen_datacard=None, gen_workspace=None, fit_datacard=None, fit_workspace=None, n_toys=1000, mu=0, dry_run=True):
	print "Running bias study " + name
	working_directory = top_directory + name
	print "\tDirectory = " + working_directory
	os.system("mkdir -pv " + working_directory)
	script_path = working_directory + "/run.sh"
	script = open(script_path, 'w')
	script.write("#!/bin/bash\n")
	script.write("mkdir higgsCombine" + name)
	script.write("combine -M GenerateOnly --expectSignal " + str(mu) + " -n " + name + " -t " + str(n_toys) + " --saveToys --rMin -1000 --rMax 1000 " + os.path.basename(gen_datacard) + "\n")
	script.write("combine -M MaxLikelihoodFit -n " + name + " -t " + str(n_toys) + " --out test_fit --saveNormalizations --toysFile higgsCombine" + name + ".GenerateOnly.mH120.123456.root --rMin -1000 --rMax 1000 --robustFit=0 --verbose 2 " + os.path.basename(fit_datacard) + "\n")
	script.close()
	if not dry_run:
		os.system("csub " + script_path + " --cmssw --no_retar -d " + working_directory + " -F " + gen_datacard + "," + fit_datacard + "," + gen_workspace + "," + fit_workspace)

def run_many_bias_studies(top_name, job_names, gen_datacards, gen_workspaces, fit_datacards, fit_workspaces, n_toys, mu, dry_run=True):
	working_directory = top_directory + top_name
	print "\tDirectory = " + working_directory
	os.system("mkdir -pv " + working_directory)
	script_paths = []
	files_to_transfer = []
	for name in job_names:
		script_path = working_directory + "/run_" + name + ".sh"
		script_paths.append(script_path)
		script = open(script_path, 'w')
		script.write("#!/bin/bash\n")

		script.write("echo $PWD\n")
		script.write("ls -lrth\n")
		script.write("mkdir -pv higgsCombine" + name + "\n")
		script.write("combine -M GenerateOnly --expectSignal " + str(mu) + " -n " + name + " -t " + str(n_toys) + " --saveToys --rMin -1000 --rMax 1000 " + os.path.basename(gen_datacards[name]) + "\n")
		script.write("combine -M MaxLikelihoodFit -n " + name + " -t " + str(n_toys) + " --saveNormalizations --toysFile higgsCombine" + name + ".GenerateOnly.mH120.123456.root --rMin -1000 --rMax 1000 --robustFit=0 --verbose 2 " + os.path.basename(fit_datacards[name]) + "\n") #  --out test_fit
		script.write("echo $PWD\n")
		script.write("ls -lrth\n")
		script.write("rm -f higgsCombine" + name + ".GenerateOnly.mH120.123456.root")
		files_to_transfer.append(gen_datacards[name])
		files_to_transfer.append(fit_datacards[name])
		files_to_transfer.append(gen_workspaces[name])
		files_to_transfer.append(fit_workspaces[name])
		files_to_transfer.append(script_path)
		script.close()

	files_to_transfer = list(set(files_to_transfer))

	# Master script
	master_script_path = working_directory + "/run_master.sh"
	master_script = open(master_script_path, "w")
	master_script.write("#!/bin/bash\n")
	master_script.write("job_scripts=( " + " ".join([os.path.basename(x) for x in script_paths]) + " )\n")
	master_script.write("source ${job_scripts[$1]}\n")
	master_script.close()

	if not dry_run:
		print "csub " + master_script_path + " --cmssw --no_retar -d " + working_directory + " -F " + ",".join(files_to_transfer) + " -n " + str(len(script_paths))
		os.system("csub " + master_script_path + " --cmssw --no_retar -d " + working_directory + " -F " + ",".join(files_to_transfer) + " -n " + str(len(script_paths)))
	else:
		print "csub " + master_script_path + " --cmssw --no_retar -d " + working_directory + " -F " + ",".join(files_to_transfer) + " -n " + str(len(script_paths))

# 2D histogram plot
# x-axis = functions (25 combinations)
# y-axis = mass
# One per model
def plot_all_averages(model, analysis, gen_functions, fit_functions, mu=0, twosig_mu_values=None):
	#results_tree = TTree("bias_study_results", "bias_study_results")
	#containers = {}
	##containers["f_gen"] = array.array('i', [0])
	##containers["f_fit"] = array.array('i', [0])
	#containers["mass"] = array.array('d', [0])
	#containers["average_mu"] = array.array('d', [0])
	#containers["average_pull"] = array.array('d', [0])
	#containers["average_pull_centered"] = array.array('d', [0])
	#containers["rms_pull"] = array.array('d', [0])
	#for container_name, container in containers.iteritems():
	#	print "Making branch for " + container_name
	#	results_tree.Branch(container_name, container)

	print "[debug] model = {}, analysis = {}, mu = {}".format(model, analysis, mu)
	print "[debug] twosig_mu_values = ",
	print twosig_mu_values

	h_avg_mu = TH2D("h_avg_mu", "h_avg_mu", len(fit_functions)*len(gen_functions), 0.5, len(fit_functions)*len(gen_functions)+0.5, len(analysis_config.simulation.limit_signal_masses[analysis]), analysis_config.simulation.limit_signal_masses[analysis][0] - 25., analysis_config.simulation.limit_signal_masses[analysis][-1] + 25.)
	h_avg_pull = TH2D("h_avg_pull", "h_avg_pull", len(fit_functions)*len(gen_functions), 0.5, len(fit_functions)*len(gen_functions)+0.5, len(analysis_config.simulation.limit_signal_masses[analysis]), analysis_config.simulation.limit_signal_masses[analysis][0] - 25., analysis_config.simulation.limit_signal_masses[analysis][-1] + 25.)
	h_avg_centered_pull = TH2D("h_avg_centered_pull", "h_avg_pull", len(fit_functions)*len(gen_functions), 0.5, len(fit_functions)*len(gen_functions)+0.5, len(analysis_config.simulation.limit_signal_masses[analysis]), analysis_config.simulation.limit_signal_masses[analysis][0] - 25., analysis_config.simulation.limit_signal_masses[analysis][-1] + 25.)
	h_rms_pull = TH2D("h_rms_pull", "h_rms_pull", len(fit_functions)*len(gen_functions), 0.5, len(fit_functions)*len(gen_functions)+0.5, len(analysis_config.simulation.limit_signal_masses[analysis]), analysis_config.simulation.limit_signal_masses[analysis][0] - 25., analysis_config.simulation.limit_signal_masses[analysis][-1] + 25.)
	h_avg_centered_pull_single_fit = {}
	h_avg_centered_pull_single_fit_single_gen = {}
	x_bin = 1
	for f_fit in fit_functions:
		h_avg_centered_pull_single_fit_single_gen[f_fit] = {}
		h_avg_centered_pull_single_fit[f_fit] = TH2D("h_avg_centered_pull_" + f_fit, "h_avg_pull_" + f_fit, len(gen_functions), -0.5, len(gen_functions) + 0.5, len(analysis_config.simulation.limit_signal_masses[analysis]), analysis_config.simulation.limit_signal_masses[analysis][0] - 25., analysis_config.simulation.limit_signal_masses[analysis][-1] + 25.)
		h_avg_centered_pull_single_fit[f_fit].SetDirectory(0)
		x_bin_single_fit = 1
		#containers["f_fit"][0] = int(f_fit[1:])
		for f_gen in gen_functions:
			h_avg_centered_pull_single_fit_single_gen[f_fit][f_gen] = TH1D("h_avg_centered_pull_fit{}_gen{}".format(f_fit, f_gen), "h_avg_centered_pull_fit{}_gen{}".format(f_fit, f_gen), len(analysis_config.simulation.limit_signal_masses[analysis]), analysis_config.simulation.limit_signal_masses[analysis][0] - 25., analysis_config.simulation.limit_signal_masses[analysis][-1] + 25.)
			h_avg_centered_pull_single_fit_single_gen[f_fit][f_gen].SetDirectory(0)
			for binm1, this_mass in enumerate(analysis_config.simulation.limit_signal_masses[analysis]):
				h_avg_centered_pull_single_fit_single_gen[f_fit][f_gen].GetXaxis().SetBinLabel(binm1+1, "{} GeV".format(this_mass))

			#containers["f_gen"][0] = int(f_gen[1:])
			h_avg_mu.GetXaxis().SetBinLabel(x_bin, "Gen " + f_gen + " / Fit " + f_fit)
			h_avg_pull.GetXaxis().SetBinLabel(x_bin, "Gen " + f_gen + " / Fit " + f_fit)
			h_avg_centered_pull.GetXaxis().SetBinLabel(x_bin, "Gen " + f_gen + " / Fit " + f_fit)
			h_rms_pull.GetXaxis().SetBinLabel(x_bin, "Gen " + f_gen + " / Fit " + f_fit)
			h_avg_centered_pull_single_fit[f_fit].GetXaxis().SetBinLabel(x_bin_single_fit, "Gen " + f_gen + " / Fit " + f_fit)
			y_bin = 1
			for binm1, mass in enumerate(analysis_config.simulation.limit_signal_masses[analysis]):
				#containers["mass"][0] = mass
				#print "[plot_all_averages] INFO : Opening " + analysis_config.get_bias_study_results(model, analysis, mass, mu, f_gen, f_fit)
				results_file = TFile(analysis_config.get_bias_study_results(model, analysis, mass, mu, f_gen, f_fit), "READ")
				if not results_file.IsOpen():
					print "[plot_all_averages] WARNING : Results file not found for model={}, analysis={}, mass={}, mu={}, f_gen={}, f_fit={}".format(model, analysis, mass, mu, f_gen, f_fit)
					continue
				t = results_file.Get("tree_fit_sb")
				if not t:
					print "[plot_all_averages] WARNING : File found but tree not available for model={}, analysis={}, mass={}, mu={}, f_gen={}, f_fit={}".format(model, analysis, mass, mu, f_gen, f_fit)
					continue
				weights = 0
				mu_sum = 0.
				mu2_sum = 0.
				pull_sum = 0.
				pull2_sum = 0.
				centered_pull_sum = 0.
				centered_pull2_sum = 0.
				t.SetBranchStatus("*", 0)
				containers = {}
				containers["mu"] = array.array('d', [0.])
				containers["muErr"] = array.array('d', [0.])
				containers["fit_status"] = array.array('i', [0])
				for branch_name, branch_container in containers.iteritems():
					t.SetBranchStatus(branch_name, 1)
					t.SetBranchAddress(branch_name, branch_container)
				for entry in xrange(t.GetEntriesFast()):
					t.GetEntry(entry)
					if containers["fit_status"][0] == -1:
						continue
					if containers["muErr"][0] <= 0:
						continue
					if abs(containers["mu"][0] / containers["muErr"][0]) > 10:
						continue
					mu_sum += containers["mu"][0]
					mu2_sum += containers["mu"][0]**2
					pull_sum += containers["mu"][0] / containers["muErr"][0]
					pull2_sum += (containers["mu"][0] / containers["muErr"][0])**2
					if mu != 0:
						injected_mu = mu * twosig_mu_values[mass]
					else:
						injected_mu = 0.
					centered_pull_sum += (containers["mu"][0] - injected_mu) / containers["muErr"][0]
					centered_pull2_sum += ((containers["mu"][0] - injected_mu) / containers["muErr"][0])**2
					weights += 1
				if weights > 0:
					mu_avg = mu_sum / weights
					mu2_avg = mu2_sum / weights
					pull_avg = pull_sum / weights
					pull2_avg = pull2_sum / weights
					centered_pull_avg = centered_pull_sum / weights
					centered_pull2_avg = centered_pull2_sum / weights
					mu_rms = (mu2_avg - (mu_avg**2))**0.5
					pull_rms = (pull2_avg - (pull_avg**2))**0.5
				else:
					mu_avg = 1.e20
					mu2_avg = 1.e20
					pull_avg = 1.e20
					pull2_avg = 1.e20
					centered_pull_avg = 1.e20
					centered_pull2_avg = 1.e20
					mu_rms = 1.e20
					pull_rms = 1.e20

				h_avg_mu.SetBinContent(x_bin, y_bin, mu_avg)
				h_avg_pull.SetBinContent(x_bin, y_bin, pull_avg)
				h_avg_centered_pull.SetBinContent(x_bin, h_avg_centered_pull.GetYaxis().FindBin(mass), centered_pull_avg)
				h_rms_pull.SetBinContent(x_bin, y_bin, pull_rms)
				h_avg_centered_pull_single_fit[f_fit].SetBinContent(x_bin_single_fit, h_avg_centered_pull_single_fit[f_fit].GetYaxis().FindBin(mass), centered_pull_avg)
				h_avg_centered_pull_single_fit_single_gen[f_fit][f_gen].SetBinContent(binm1+1, centered_pull_avg)
				#containers["average_mu"][0] = mu_avg
				#containers["average_pull"][0] = pull_avg
				#containers["average_pull_centered"][0] = centered_pull_avg
				#containers["rms_pull"][0] = pull_rms
				#results_tree.Fill()
				y_bin += 1
			x_bin += 1
			x_bin_single_fit += 1
	c_avg_mu = TCanvas("c_avg_mu_f_vs_mass_" + model + "_" + analysis + "_" + str(mu), "c_avg_mu_f_vs_mass", 800, 600)
	c_avg_mu.SetRightMargin(0.2)
	c_avg_mu.SetBottomMargin(0.15)
	h_avg_mu.GetXaxis().SetTitle("Gen/fit function")
	h_avg_mu.GetYaxis().SetTitle("Signal mass")
	h_avg_mu.GetZaxis().SetTitle("#LT#mu#GT")
	h_avg_mu.GetXaxis().SetTitleOffset(1.3)
	h_avg_mu.GetYaxis().SetTitleOffset(1.1)
	h_avg_mu.SetMinimum(mu - 0.1)
	h_avg_mu.SetMaximum(mu + 0.1)

	h_avg_mu.Draw("colz")
	c_avg_mu.SaveAs(analysis_config.figure_directory + "/" + c_avg_mu.GetName() + ".pdf")

	c_avg_pull = TCanvas("c_avg_pull_f_vs_mass_" + model + "_" + analysis + "_" + str(mu), "c_avg_pull_f_vs_mass", 800, 600)
	c_avg_pull.SetRightMargin(0.2)
	c_avg_pull.SetBottomMargin(0.15)
	h_avg_pull.GetXaxis().SetTitle("Gen/fit function")
	h_avg_pull.GetYaxis().SetTitle("Signal mass")
	h_avg_pull.GetZaxis().SetTitle("#LT#mu/#sigma_{#mu}#GT")
	h_avg_pull.GetXaxis().SetTitleOffset(1.3)
	h_avg_pull.GetYaxis().SetTitleOffset(1.1)
	if mu == 0:
		h_avg_pull.SetMinimum(-5.)
		h_avg_pull.SetMaximum(5.)
	h_avg_pull.Draw("colz")
	c_avg_pull.SaveAs(analysis_config.figure_directory + "/" + c_avg_pull.GetName() + ".pdf")

	c_avg_centered_pull = TCanvas("c_avg_centered_pull_f_vs_mass_" + model + "_" + analysis + "_" + str(mu), "c_avg_centered_pull_f_vs_mass", 800, 600)
	c_avg_centered_pull.SetRightMargin(0.2)
	c_avg_centered_pull.SetBottomMargin(0.15)
	h_avg_centered_pull.GetXaxis().SetTitle("Gen/fit function")
	h_avg_centered_pull.GetYaxis().SetTitle("Signal mass")
	h_avg_centered_pull.GetZaxis().SetTitle("#LT(#mu-#mu_{inj})/#sigma_{#mu}#GT")
	h_avg_centered_pull.GetXaxis().SetTitleOffset(1.3)
	h_avg_centered_pull.GetYaxis().SetTitleOffset(1.1)
	h_avg_centered_pull.SetMinimum(-3)
	h_avg_centered_pull.SetMaximum(3)
	h_avg_centered_pull.Draw("colz")
	c_avg_centered_pull.SaveAs(analysis_config.figure_directory + "/" + c_avg_centered_pull.GetName() + ".pdf")

	c_rms_pull = TCanvas("c_rms_pull_f_vs_mass_" + model + "_" + analysis + "_" + str(mu), "c_rms_pull_f_vs_mass", 800, 600)
	c_rms_pull.SetRightMargin(0.2)
	c_rms_pull.SetBottomMargin(0.15)
	h_rms_pull.GetXaxis().SetTitle("Gen/fit function")
	h_rms_pull.GetYaxis().SetTitle("Signal mass")
	h_rms_pull.GetZaxis().SetTitle("#sigma(#mu/#sigma_{#mu})")
	h_rms_pull.GetXaxis().SetTitleOffset(1.3)
	h_rms_pull.GetYaxis().SetTitleOffset(1.1)
	h_rms_pull.SetMinimum(-1.)
	h_rms_pull.SetMaximum(1.)
	h_rms_pull.Draw("colz")
	c_rms_pull.SaveAs(analysis_config.figure_directory + "/" + c_rms_pull.GetName() + ".pdf")

	for f_fit in fit_functions:
		c_avg_centered_pull_single_fit = TCanvas("c_avg_centered_pull_f_vs_mass_" + model + "_" + analysis + "_" + str(mu) + "_" + f_fit, "c_avg_centered_pull_f_vs_mass", 800, 600)
		c_avg_centered_pull_single_fit.SetRightMargin(0.2)
		c_avg_centered_pull_single_fit.SetBottomMargin(0.15)
		h_avg_centered_pull_single_fit[f_fit].GetXaxis().SetTitle("Gen/fit function")
		h_avg_centered_pull_single_fit[f_fit].GetYaxis().SetTitle("Signal mass")
		h_avg_centered_pull_single_fit[f_fit].GetZaxis().SetTitle("#LT(#mu-#mu_{inj})/#sigma_{#mu}#GT")
		h_avg_centered_pull_single_fit[f_fit].GetXaxis().SetTitleOffset(1.3)
		h_avg_centered_pull_single_fit[f_fit].GetYaxis().SetTitleOffset(1.1)
		h_avg_centered_pull_single_fit[f_fit].SetMinimum(-3)
		h_avg_centered_pull_single_fit[f_fit].SetMaximum(3)
		h_avg_centered_pull_single_fit[f_fit].Draw("colz text")
		c_avg_centered_pull_single_fit.SaveAs(analysis_config.figure_directory + "/" + c_avg_centered_pull_single_fit.GetName() + ".pdf")
		for f_gen in gen_functions:
			c_avg_centered_pull_single_fit_single_gen = TCanvas("c_avg_centered_pull_f_vs_mass_" + model + "_" + analysis + "_" + str(mu) + "_fit" + f_fit + "_gen" + f_gen, "c_avg_centered_pull_f_vs_mass", 800, 600)
			c_avg_centered_pull_single_fit_single_gen.SetGrid()
			h_avg_centered_pull_single_fit_single_gen[f_fit][f_gen].GetXaxis().SetTitle("Signal mass [GeV]")
			h_avg_centered_pull_single_fit_single_gen[f_fit][f_gen].GetYaxis().SetTitle("#LT(#mu-#mu_{inj})/#sigma_{#mu}#GT")
			h_avg_centered_pull_single_fit_single_gen[f_fit][f_gen].SetMinimum(-1.)
			h_avg_centered_pull_single_fit_single_gen[f_fit][f_gen].SetMaximum(1.)
			h_avg_centered_pull_single_fit_single_gen[f_fit][f_gen].SetLineColor(seaborn.GetColorRoot("default", 2))
			h_avg_centered_pull_single_fit_single_gen[f_fit][f_gen].SetLineWidth(2)
			h_avg_centered_pull_single_fit_single_gen[f_fit][f_gen].Draw("hist")
			c_avg_centered_pull_single_fit_single_gen.SaveAs(analysis_config.figure_directory + "/" + c_avg_centered_pull_single_fit_single_gen.GetName() + ".pdf")
	results_file = TFile(analysis_config.bias_study_directory + "/summary_" + model + "_" + analysis + "_mu" + str(mu) + ".root", "RECREATE")
	h_avg_mu.Write()
	h_avg_pull.Write()
	h_avg_centered_pull.Write()
	h_rms_pull.Write()
	#results_tree.Write()
	results_file.Close()

# Plot: x-axis = gen function, y-axis = fit function
# One per mass point and model

# Plot: pull histogram
def plot_pull_distributions(model, analysis, mass, f_gen, f_fit, mu=0, twosig_mu_values=None):
	h = TH1D("h", "h", 40, -5., 5.)
	results_file = TFile(analysis_config.get_bias_study_results(model, analysis, mass, mu, f_gen, f_fit), "READ")
	if not results_file.IsOpen():
		print "[plot_all_averages] WARNING : Results file not found for model={}, analysis={}, mass={}, mu={}, f_gen={}, f_fit={}".format(model, analysis, mass, mu, f_gen, f_fit)
		return
	t = results_file.Get("tree_fit_sb")
	if not t:
		print "[plot_all_averages] WARNING : File found but tree not available for model={}, analysis={}, mass={}, mu={}, f_gen={}, f_fit={}".format(model, analysis, mass, mu, f_gen, f_fit)
		return

	containers = {}
	containers["fit_status"] = array.array('i', [0])
	containers["mu"] = array.array('d', [0.])
	containers["muErr"] = array.array('d', [0.])
	for branch_name, branch_container in containers.iteritems():
		t.SetBranchStatus(branch_name, 1)
		t.SetBranchAddress(branch_name, branch_container)
	for entry in xrange(t.GetEntriesFast()):
		t.GetEntry(entry)
		if containers["fit_status"][0] == -1:
			continue
		if containers["muErr"][0] <= 0:
			continue
		if abs(containers["mu"][0] / containers["muErr"][0]) > 10:
			continue

		if mu != 0:
			injected_mu = mu * twosig_mu_values[mass]
		else:
			injected_mu = 0.
		centered_pull = (containers["mu"][0] - injected_mu) / containers["muErr"][0]
		h.Fill(centered_pull)
	c = TCanvas("c_pulldist_{}_{}_{}_gen{}_fit{}_mu{}".format(model, analysis, mass, f_gen, f_fit, mu), "Pull", 700, 500)
	h.GetXaxis().SetTitle("#frac{#hat{#mu}-#mu_{inj}}{#sigma_{#hat{#mu}}}")
	h.Draw("hist")
	if "trigbbl_CSVTM" in analysis:
		sr = "Low mass SR"
	else:
		sr = "High mass SR"
	Root.myText(0.15, 0.85, 1, "{} / {}".format(model, sr), 0.5)
	Root.myText(0.15, 0.8, 1, "m_{{X}}={} GeV".format(mass), 0.5)

	Root.myText(0.15, 0.75, 1, "Gen {} / Fit {}".format(f_gen, f_fit), 0.5)
	if mu == 0:
		Root.myText(0.15, 0.7, 1, "#mu_{inj}=0", 0.5)
	elif mu == 1:
		Root.myText(0.15, 0.7, 1, "#mu_{inj}#approx +2#sigma_{exp}", 0.5)
	Root.myText(0.7, 0.8, 1, "Mean={:.2f}".format(h.GetMean()), 0.5)
	Root.myText(0.7, 0.75, 1, "RMS={:.2f}".format(h.GetRMS()), 0.5)

	c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")



if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description = 'Run bias studies and plot output')
	parser.add_argument('--model', type=str, default="Hbb,RSG", help='Model name')
	parser.add_argument('--analysis', type=str, default="trigbbl_CSVTM,trigbbh_CSVTM", help='Analysis name')
	parser.add_argument('--run', action='store_true', help='Run bias studies')
	parser.add_argument('--retry', action='store_true', help='Rerun failed bias studies')
	parser.add_argument('--dry_run', action='store_true', help='Setup but do not run bias studies')
	parser.add_argument('--n_toys', type=int, default=1000, help='Number of toys to run')
	parser.add_argument('--mu', type=int, default=0, help='Signal mu (1 = inject signal at +2sigma expected limit)')
	parser.add_argument('--test', action='store_true', help='Small test job')
	parser.add_argument('--plots', action='store_true', help='Plot results')
	args = parser.parse_args()

	models = args.model.split(",")
	analyses = args.analysis.split(",")
	masses = {"trigbbl_CSVTM":range(350, 850, 50), "trigbbh_CSVTM":range(600, 1250, 50)}
	masses["trigbbl_CSVTM"].append(325)
	#masses = {"trigbbl_CSVTM":[325], "trigbbh_CSVTM":[]}

	job_mu_values = {}
	for model in models:
		job_mu_values[model] = {}
		for analysis in analyses:
			job_mu_values[model][analysis] = {}
			for mass in masses[analysis]:
				if args.mu == 0:
					job_mu_values[model][analysis][mass] = 0.
				else:
					workspace_file = TFile(limit_config.get_workspace_filename(analysis, model, mass, correctTrigger=False, useMCTrigger=False, qcd=False, fitTrigger=True, fitBonly=False), "READ")
					#workspace_file = TFile(datacard_folders[model][analysis] + "/workspace_qq_m" + str(mass) + ".root", "READ")
					workspace = workspace_file.Get("w")
					signal_norm = workspace.var("signal_norm").getVal()
					job_mu_values[model][analysis][mass] = 19700. * limit_config.limit_p2sigma_estimates[analysis][model][mass] / signal_norm

	#if args.run:
	#	if args.test:
	#		run_bias_study("test", gen_datacard=datacard_folders["Hbb"]["trigbbh_CSVTM"] + "/datacard_qq_m" + str(750) + "_f3.txt", fit_datacard=datacard_folders["Hbb"]["trigbbh_CSVTM"] + "/datacard_qq_m" + str(750) + "_f1.txt", gen_workspace=datacard_folders["Hbb"]["trigbbh_CSVTM"] + "/workspace_qq_m" + str(750) + ".root", fit_workspace=datacard_folders["Hbb"]["trigbbh_CSVTM"] + "/workspace_qq_m" + str(750) + ".root", n_toys=args.n_toys, mu=args.mu)
	#	else:
	#		for model in models:
	#			for analysis in analyses:
	#				for mass in masses[analysis]:
	#					for f_gen in functions:
	#						gen_datacard = datacard_folders[model][analysis] + "/datacard_qq_m" + str(mass) + "_" + f_gen + ".txt"
	#						gen_workspace = datacard_folders[model][analysis] + "/workspace_qq_m" + str(mass) + ".root"
	#						for f_fit in functions:
	#							fit_datacard = datacard_folders[model][analysis] + "/datacard_qq_m" + str(mass) + "_" + f_fit + ".txt"
	#							fit_workspace = datacard_folders[model][analysis] + "/workspace_qq_m" + str(mass) + ".root"
	#							name = model + "_" + analysis + "/m" + str(mass) + "_gen_" + f_gen + "_fit_" + f_fit + "_mu" + str(args.mu)
	#							run_bias_study(name, gen_datacard=gen_datacard, fit_datacard=fit_datacard, gen_workspace=gen_workspace, fit_workspace=fit_workspace, n_toys=args.n_toys, mu=args.mu, dry_run=args.dry_run)
	if args.run:
		for model in models:
			for analysis in analyses:
				if analysis == "trigbbl_CSVTM":
					gen_functions = ["dijet4", "polypower4", "polyx6"]
					fit_functions = ["dijet4"]
					#gen_functions = ["f1", "f3", "f4"]
					#fit_functions = ["f1"]
				else:
					gen_functions = ["dijet4", "polypower4"]
					fit_functions = ["dijet4"]
					#gen_functions = ["f1", "f2", "f4", "f5"]
					#fit_functions = ["f1"]
				for mass in masses[analysis]:
					names = []
					gen_datacards = {}
					fit_datacards = {}
					gen_workspaces = {}
					fit_workspaces = {}
					top_name = model + "_" + analysis + "/m" + str(mass) + "_mu" + str(args.mu)
					job_names = []
					for f_gen in gen_functions:
						for f_fit in fit_functions:
							name = model + "_" + analysis + "_m" + str(mass) + "_gen_" + f_gen + "_fit_" + f_fit + "_mu" + str(args.mu)
							job_names.append(name)
							gen_datacards[name] = limit_config.get_datacard_filename(analysis, model, mass, f_gen, correctTrigger=False, useMCTrigger=False, qcd=False, fitTrigger=True, fitBonly=False)
							gen_workspaces[name] = limit_config.get_workspace_filename(analysis, model, mass, correctTrigger=False, useMCTrigger=False, qcd=False, fitTrigger=True, fitBonly=False)
							fit_datacards[name] = limit_config.get_datacard_filename(analysis, model, mass, f_fit, correctTrigger=False, useMCTrigger=False, qcd=False, fitTrigger=True, fitBonly=False)
							fit_workspaces[name] = limit_config.get_workspace_filename(analysis, model, mass, correctTrigger=False, useMCTrigger=False, qcd=False, fitTrigger=True, fitBonly=False)
							#gen_datacards[name] = datacard_folders[model][analysis] + "/datacard_qq_m" + str(mass) + "_" + f_gen + ".txt"
							#gen_workspaces[name] = datacard_folders[model][analysis] + "/workspace_qq_m" + str(mass) + ".root"
							#fit_datacards[name] = datacard_folders[model][analysis] + "/datacard_qq_m" + str(mass) + "_" + f_fit + ".txt"
							#fit_workspaces[name] = datacard_folders[model][analysis] + "/workspace_qq_m" + str(mass) + ".root"
					run_many_bias_studies(top_name, job_names, gen_datacards=gen_datacards, fit_datacards=fit_datacards, gen_workspaces=gen_workspaces, fit_workspaces=fit_workspaces, n_toys=args.n_toys, mu=job_mu_values[model][analysis][mass], dry_run=args.dry_run)

	if args.retry:
		for model in models:
			for analysis in analyses:
				for mass in masses[analysis]:
					names = []
					gen_datacards = {}
					fit_datacards = {}
					gen_workspaces = {}
					fit_workspaces = {}
					top_name = model + "_" + analysis + "/m" + str(mass) + "_mu" + str(args.mu)
					job_names = []
					for f_gen in functions:
						for f_fit in functions:
							# Check if the output file exists, and the tree is readable
							results_file = TFile(analysis_config.get_bias_study_results(model, analysis, mass, args.mu, f_gen, f_fit), "READ")
							retry_this_job = False
							if not results_file.IsOpen():
								retry_this_job = True
							else:
								t = results_file.Get("tree_fit_sb")
								if not t:
									retry_this_job = True
							if retry_this_job:
								name = model + "_" + analysis + "_m" + str(mass) + "_gen_" + f_gen + "_fit_" + f_fit + "_mu" + str(args.mu)
								job_names.append(name)
								gen_datacards[name] = limit_config.get_datacard_filename(analysis, model, mass, f_gen, correctTrigger=False)
								gen_workspaces[name] = limit_config.get_workspace_filename(analysis, model, mass, correctTrigger=False)
								fit_datacards[name] = limit_config.get_datacard_filename(analysis, model, mass, f_fit, correctTrigger=False)
								fit_workspaces[name] = limit_config.get_workspace_filename(analysis, model, mass, correctTrigger=False)
								#gen_datacards[name] = datacard_folders[model][analysis] + "/datacard_qq_m" + str(mass) + "_" + f_gen + ".txt"
								#gen_workspaces[name] = datacard_folders[model][analysis] + "/workspace_qq_m" + str(mass) + ".root"
								#fit_datacards[name] = datacard_folders[model][analysis] + "/datacard_qq_m" + str(mass) + "_" + f_fit + ".txt"
								#fit_workspaces[name] = datacard_folders[model][analysis] + "/workspace_qq_m" + str(mass) + ".root"
					if len(job_names) > 0:
						run_many_bias_studies(top_name, job_names, gen_datacards=gen_datacards, fit_datacards=fit_datacards, gen_workspaces=gen_workspaces, fit_workspaces=fit_workspaces, n_toys=args.n_toys, mu=args.mu, dry_run=args.dry_run)

	if args.plots:
		print job_mu_values
		for model in models:
			for analysis in analyses:
				if analysis == "trigbbl_CSVTM":
					gen_functions = ["dijet4", "polypower4", "polyx6"]
					fit_functions = ["dijet4"]
					#gen_functions = ["f1", "f3", "f4"]
					#fit_functions = ["f1"]
				else:
					gen_functions = ["dijet4", "polypower4"]
					fit_functions = ["dijet4"]
					#gen_functions = ["f1", "f2", "f4", "f5"]
					#fit_functions = ["f1"]
				plot_all_averages(model, analysis, gen_functions, fit_functions, mu=args.mu, twosig_mu_values=job_mu_values[model][analysis])
				for mass in masses[analysis]:
					for f_gen in gen_functions:
						for f_fit in fit_functions:
							plot_pull_distributions(model, analysis, mass, f_gen, f_fit, mu=args.mu, twosig_mu_values=job_mu_values[model][analysis])
