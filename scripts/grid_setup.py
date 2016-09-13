#!/usr/bin/env python

import sys, os, re
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
from CMSDIJET.StatisticalTools.systematics import *
from CMSDIJET.StatisticalTools.roofit_functions import *
from math import log

def setup_grid(name, datacard, workspace, r_values, run=False):
	start_dir = os.getcwd()
	working_dir = limit_config.paths["r_grid"] + "/" + name + "/"
	os.system("mkdir -pv " + working_dir)
	os.chdir(working_dir)
	bash_script_path = working_dir + "/setup_" + name + ".sh"
	bash_script = open(bash_script_path, 'w')
	bash_script.write("#!/bin/bash\n")
	bash_script.write("r_values=( " + " ".join([str(x) for x in r_values]) + " )\n")
	bash_script.write("R_VALUE=${r_values[$1]}\n")
	bash_script.write("echo \"Running combine with r=$R_VALUE\"\n")
	bash_script.write("combine {} -M HybridNew --frequentist -s 123456 --singlePoint $R_VALUE --saveToys --saveHybridResult --name {}_i$1 \n".format(os.path.basename(datacard), name))
	bash_script.close()

	submit_command = "csub {} --cmssw --no_retar -F {} -n {}".format(bash_script_path, ",".join([datacard, workspace]), len(r_values))
	if run:
		os.system(submit_command)
	else:
		print "To run:"
		print submit_command
		print "All done."
	os.chdir(start_dir)
	return submit_command

def merge(analysis, model, mass, fit_function):
	merge_command = "hadd " + limit_config.get_hn_grid(analysis, model, mass, args.fit_function) + " " + " ".join(limit_config.get_hn_grid_subfiles(analysis, model, mass, args.fit_function))
	os.system(merge_command)


if __name__ == "__main__":
	from argparse import ArgumentParser
	parser = ArgumentParser(description='Setup HybridNew r-grid for computing expected limits')
	parser.add_argument("--run", action='store_true', help="Run grid setup jobs")
	parser.add_argument("--merge", action='store_true', help="Merge outputs of grid setup jobs")

	parser.add_argument('--analyses', type=str, default="trigbbh_CSVTM,trigbbl_CSVTM", help="Analysis names (comma-separated)")
	parser.add_argument('--models', type=str, default="Hbb,RSG", help="Model names (comma-separated)")
	parser.add_argument('--fit_function', type=str, default="f3", help="Name of central fit function")
	parser.add_argument('--n_points', type=int, default=200, help="Number of r values in the grid")
	parser.add_argument('--r_min', type=float, default=0.005, help="Minimum r value")
	parser.add_argument('--r_max', type=float, default=0.5, help="Maximum r value")

	# Combine options to pass through
	parser.add_argument("--rRelAcc", dest="rRelAcc", type=float, help="rRelAcc")
	parser.add_argument("--rAbsAcc", dest="rAbsAcc", type=float, help="rAbsAcc")
	parser.add_argument("--toysH", dest="toysH", type=int, default=500, help="Number of Toy MC extractions for HybridNew")
	parser.add_argument("--noSyst", dest="noSyst", default=False, action="store_true", help="Run without systematic uncertainties")
	parser.add_argument("--fork", type=int, default=4, help="Fork N processes")
	args = parser.parse_args()

	models = args.models.split(",")
	analyses = args.analyses.split(",")

	if args.run:
		# Make list of r values in grid
		min_logr = log(args.r_min, 10)
		max_logr = log(args.r_max, 10)
		logr_values = [min_logr + ((max_logr - min_logr)*i/(args.n_points - 1)) for i in range(args.n_points)]
		r_values = [10**x for x in logr_values]

		postfix = args.fit_function
		if args.noSyst:
			postfix += "_noSyst"

		for analysis in analyses:
			for model in models:
				masses = limit_config.limit_signal_masses[analysis]
				for mass in masses:
					name = "grid_{}_{}_{}_{}".format(analysis, model, mass, postfix)
					datacard = limit_config.get_datacard_filename(analysis, model, mass, args.fit_function, fitSignal=True)
					workspace = limit_config.get_workspace_filename(analysis, model, mass, fitSignal=True)
					setup_grid(name, datacard, workspace, r_values, run=True)

	if args.merge:
		from joblib import Parallel, delayed
		for analysis in analyses:
			Parallel(n_jobs=4)(delayed(merge)(analysis, model, mass, args.fit_function) for mass in limit_config.limit_signal_masses[analysis] for model in models)

