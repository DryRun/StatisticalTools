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

top_directory = "/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/GranularityStudies/"

def run_many_granularity_studies(top_name, job_names, gen_datacard, gen_workspace, fit_datacards, fit_workspaces, n_toys, mu, dry_run=False):
	working_directory = top_directory + top_name
	print "\tDirectory = " + working_directory
	os.system("mkdir -pv " + working_directory)
	script_paths = []
	files_to_transfer = []
	for name in job_names:
		print "Adding job " + name
		script_path = working_directory + "/run_" + name + ".sh"
		script_paths.append(script_path)
		script = open(script_path, 'w')
		script.write("#!/bin/bash\n")

		script.write("echo $PWD\n")
		script.write("ls -lrth\n")
		script.write("mkdir -pv higgsCombine" + name + "\n")
		script.write("combine -M GenerateOnly --expectSignal " + str(mu) + " -n " + name + " -t " + str(n_toys) + " --saveToys --rMin -1 --rMax 1000 " + os.path.basename(gen_datacard) + "\n")
		script.write("combine -M Asymptotic -n " + name + " -t " + str(n_toys) + " --toysFile higgsCombine" + name + ".GenerateOnly.mH120.123456.root --rMin -1 --rMax 300 --verbose 10 " + os.path.basename(fit_datacards[name]) + "\n") #  --out test_fit
		#script.write("combine -M MaxLikelihoodFit -n " + name + " -t " + str(n_toys) + " --saveNormalizations --toysFile higgsCombine" + name + ".GenerateOnly.mH120.123456.root --rMin -1 --rMax 300 --robustFit=0 --verbose 10 " + os.path.basename(fit_datacards[name]) + "\n") #  --out test_fit
		script.write("echo $PWD\n")
		script.write("ls -lrth\n")
		script.write("rm -f higgsCombine" + name + ".GenerateOnly.mH120.123456.root")
		files_to_transfer.append(gen_datacard)
		files_to_transfer.append(fit_datacards[name])
		files_to_transfer.append(gen_workspace)
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

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description = 'Run granularity studies and plot output')
	parser.add_argument('--model', type=str, default="Hbb,RSG", help='Model name')
	parser.add_argument('--analysis', type=str, default="trigbbl_CSVTM", help='Analysis name')
	parser.add_argument('--run', action='store_true', help='Run bias studies')
	parser.add_argument('--n_toys', type=int, default=100, help='Number of toys to run')
	parser.add_argument('--mu', type=int, default=1, help='Signal mu (1 = inject signal at +2sigma expected limit)')
	parser.add_argument('--genmass', type=str, default="375", help='Masses for generation, i.e. signal injection mass values')
	parser.add_argument('--fitmass', type=str, default="350,400", help='Masses for fitting')
	args = parser.parse_args()

	models = args.model.split(",")
	analyses = args.analysis.split(",")
	gen_masses = [int(x) for x in args.genmass.split(",")]
	fit_masses = [int(x) for x in args.fitmass.split(",")]
	job_mu_values = {}

	for model in models:
		job_mu_values[model] = {}
		for analysis in analyses:
			job_mu_values[model][analysis] = {}
			for gen_mass in gen_masses:
				if args.mu == 0:
					job_mu_values[model][analysis][gen_mass] = 0.
				else:
					workspace_file = TFile(limit_config.get_workspace_filename(analysis, model, gen_mass, correctTrigger=False, useMCTrigger=True, qcd=False, fitTrigger=True, fitBonly=False), "READ")
					#workspace_file = TFile(datacard_folders[model][analysis] + "/workspace_qq_m" + str(gen_mass) + ".root", "READ")
					workspace = workspace_file.Get("w")
					signal_norm = workspace.var("signal_norm").getVal()
					job_mu_values[model][analysis][gen_mass] = 19700. * limit_config.limit_p2sigma_estimates[analysis][model][gen_mass] / signal_norm
	if args.run:
		for model in models:
			for analysis in analyses:
				for gen_mass in gen_masses:
					names = []
					fit_datacards = {}
					fit_workspaces = {}
					gen_name = model + "_" + analysis + "/m" + str(gen_mass) + "_mu" + str(args.mu)
					gen_datacard = limit_config.get_datacard_filename(analysis, model, gen_mass, "dijet4", correctTrigger=False, useMCTrigger=True, qcd=False, fitTrigger=True, fitBonly=False)
					gen_workspace = limit_config.get_workspace_filename(analysis, model, gen_mass, correctTrigger=False, useMCTrigger=True, qcd=False, fitTrigger=True, fitBonly=False)
					top_name = model + "_" + analysis + "/genmass_" + str(gen_mass) + "_mu" + str(args.mu)
					job_names = []
					for fit_mass in fit_masses:
						fit_name = model + "_" + analysis + "/m" + str(fit_mass) + "_mu" + str(args.mu)
						name = model + "_" + analysis + "_genmass" + str(gen_mass) + "_fitmass_" + str(fit_mass) + "_mu" + str(args.mu)
						job_names.append(name)
						fit_datacards[name] = limit_config.get_datacard_filename(analysis, model, fit_mass, "dijet4", correctTrigger=False, useMCTrigger=True, qcd=False, fitTrigger=True, fitBonly=False)
						fit_workspaces[name] = limit_config.get_workspace_filename(analysis, model, fit_mass, correctTrigger=False, useMCTrigger=True, qcd=False, fitTrigger=True, fitBonly=False)
					run_many_granularity_studies(top_name, job_names, gen_datacard=gen_datacard, fit_datacards=fit_datacards, gen_workspace=gen_workspace, fit_workspaces=fit_workspaces, n_toys=args.n_toys, mu=job_mu_values[model][analysis][gen_mass])
