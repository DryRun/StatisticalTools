import os
import sys

models = ["Hbb", "RSG"]
analyses = ["trigbbl_CSVTM", "trigbbh_CSVTM"]
masses = {"trigbbl_CSVTM":range(400, 950, 50), "trigbbh_CSVTM":range(550, 1250, 50)}
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
functions = ["f" + str(x) for x in xrange(1, 6)]

top_directory = "/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/BiasStudies/"

def run_bias_study(name, gen_datacard=None, gen_workspace=None, fit_datacard=None, fit_workspace=None, n_toys=1000, mu=0, dry_run=False):
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

def run_many_bias_studies(top_name, job_names, gen_datacards, gen_workspaces, fit_datacards, fit_workspaces, n_toys, mu, dry_run=False):
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

		script.close()

	files_to_transfer = list(set(files_to_transfer))

	if not dry_run:
		os.system("csub " + ",".join(script_paths) + " --cmssw --no_retar -d " + working_directory + " -F " + ",".join(files_to_transfer))
	else:
		print "csub " + ",".join(script_paths) + " --cmssw --no_retar -d " + working_directory + " -F " + ",".join(files_to_transfer)

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description = 'Run bias studies and plot output')
	parser.add_argument('--run', action='store_true', help='Run bias studies')
	parser.add_argument('--run2', action='store_true', help='Run bias studies')
	parser.add_argument('--dry_run', action='store_true', help='Setup but do not run bias studies')
	parser.add_argument('--n_toys', type=int, default=1000, help='Number of toys to run')
	parser.add_argument('--mu', type=int, default=0, help='Signal mu')
	parser.add_argument('--test', action='store_true', help='Small test job')
	args = parser.parse_args()

	if args.run:
		if args.test:
			run_bias_study("test", gen_datacard=datacard_folders["Hbb"]["trigbbh_CSVTM"] + "/datacard_qq_m" + str(750) + "_f3.txt", fit_datacard=datacard_folders["Hbb"]["trigbbh_CSVTM"] + "/datacard_qq_m" + str(750) + "_f1.txt", gen_workspace=datacard_folders["Hbb"]["trigbbh_CSVTM"] + "/workspace_qq_m" + str(750) + ".root", fit_workspace=datacard_folders["Hbb"]["trigbbh_CSVTM"] + "/workspace_qq_m" + str(750) + ".root", n_toys=args.n_toys, mu=args.mu)
		else:
			for model in models:
				for analysis in analyses:
					for mass in masses[analysis]:
						for f_gen in functions:
							gen_datacard = datacard_folders[model][analysis] + "/datacard_qq_m" + str(mass) + "_" + f_gen + ".txt"
							gen_workspace = datacard_folders[model][analysis] + "/workspace_qq_m" + str(mass) + ".root"
							for f_fit in functions:
								fit_datacard = datacard_folders[model][analysis] + "/datacard_qq_m" + str(mass) + "_" + f_fit + ".txt"
								fit_workspace = datacard_folders[model][analysis] + "/workspace_qq_m" + str(mass) + ".root"
								name = model + "_" + analysis + "/m" + str(mass) + "_gen_" + f_gen + "_fit_" + f_fit + "_mu" + str(args.mu)
								run_bias_study(name, gen_datacard=gen_datacard, fit_datacard=fit_datacard, gen_workspace=gen_workspace, fit_workspace=fit_workspace, n_toys=args.n_toys, mu=args.mu, dry_run=args.dry_run)
	if args.run2:
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
							name = model + "_" + analysis + "_m" + str(mass) + "_gen_" + f_gen + "_fit_" + f_fit + "_mu" + str(args.mu)
							job_names.append(name)
							gen_datacards[name] = datacard_folders[model][analysis] + "/datacard_qq_m" + str(mass) + "_" + f_gen + ".txt"
							gen_workspaces[name] = datacard_folders[model][analysis] + "/workspace_qq_m" + str(mass) + ".root"
							fit_datacards[name] = datacard_folders[model][analysis] + "/datacard_qq_m" + str(mass) + "_" + f_fit + ".txt"
							fit_workspaces[name] = datacard_folders[model][analysis] + "/workspace_qq_m" + str(mass) + ".root"
					run_many_bias_studies(top_name, job_names, gen_datacards=gen_datacards, fit_datacards=fit_datacards, gen_workspaces=gen_workspaces, fit_workspaces=fit_workspaces, n_toys=args.n_toys, mu=args.mu, dry_run=args.dry_run)
