#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
from CMSDIJET.StatisticalTools.systematics import *
from CMSDIJET.StatisticalTools.roofit_functions import *

if __name__ == "__main__":
	# input parameters
	parser = ArgumentParser(description='Script that runs limit calculation for specified mass points')
	parser.add_argument('--analyses', type=str, default="trigbbh_CSVTM,trigbbl_CSVTM", help="Analysis names")
	parser.add_argument('--models', type=str, default="Hbb", help="Model names")
	parser.add_argument('--qcd', action='store_true', help="Use QCD instead of data (assumes no trigger emulation)")
	parser.add_argument('--correctTrigger', action='store_true', help="Use model with trigger correction (has to have been specified in create_datacards.py)")
	parser.add_argument('--fitTrigger', action='store_true', help="Use model with trigger fit (has to have been specified in create_datacards.py)")
	parser.add_argument('--fit_function', type=str, default="f4", help="Name of central fit function")
	parser.add_argument('--algo', type=str, default="AD", help="GoF statistic: AD, saturated, KS")
	parser.add_argument('--toys', type=int, help="Number of toys")
	parser.add_argument('--seed', type=int, default=123456, help="Random seed")
	parser.add_argument('--plots', action="store_true", help="Make plots")
	parser.add_argument('--signal', type=float, help="Fixed signal strength")
	parser.add_argument('-m', '--masses', type=str, help='Manually specify masses (comma-separated list). Otherwise, taken from limit_configuration.')
	parser.add_argument('--no_retar', action='store_true', help='Don\'t retar CMSSW directory before submitting')
	args = parser.parse_args()

	analyses = args.analyses.split(",")
	models = args.models.split(",")
	masses = {}
	if args.masses:
		for analysis in analyses:
			masses[analysis] = [int(x) for x in args.masses.split(",")]
	else:
		masses = limit_config.limit_signal_masses

	datacards_path = limit_config.paths["datacards"]
	condor_path = limit_config.paths["gof"] + "/condor/"
	# change to the appropriate directory
	os.chdir(condor_path)
	first = True

	prefix = 'gof'

	for analysis in analyses:
		for model in models:
			for mass in masses[analysis]:

				job_name = '{}_{}_{}_m{}_{}_{}'.format(prefix, analysis, model, int(mass), args.fit_function, args.algo)
				log_name = job_name + ".log"

				combine_options =  "-M GoodnessOfFit -v3 --algo {} -s {} --name {} --mass {} ".format(
					args.algo,
					args.seed,
					job_name,
					mass
				)
				if args.plots:
					combine_options += " --plots "
				if args.signal != None:
					combine_options += " --fixedSignalStrength {} ".format(args.signal)

				cmd = "combine {} {} 2>&1 | tee {}".format(
					combine_options, 
					os.path.basename(limit_config.get_datacard_filename(analysis, model, mass, args.fit_function, correctTrigger=args.correctTrigger, qcd=args.qcd, fitTrigger=args.fitTrigger)),
					log_name)

				if args.toys:
					toy_cmd = "combine {} {} -t {} 2>&1 | tee {}".format(
						combine_options, 
						os.path.basename(limit_config.get_datacard_filename(analysis, model, mass, args.fit_function, correctTrigger=args.correctTrigger, qcd=args.qcd, fitTrigger=args.fitTrigger)),
						args.toys,
						log_name)

				submission_dir = condor_path
				start_dir = os.getcwd()
				os.chdir(submission_dir)

				# create the executable script
				bash_script_path = os.path.join(submission_dir,'run_{}.sh'.format(job_name))
				#bash_content = bash_template
				#bash_content = re.sub('DUMMY_CMD',cmd,bash_content)

				bash_script = open(bash_script_path,'w')
				bash_script.write("echo '" + cmd + "'\n")
				bash_script.write(cmd + "\n")
				if args.toys:
					bash_script.write("echo '" + toy_cmd + "'\n")
					bash_script.write(toy_cmd + "\n")
				bash_script.close()

				# Create the condor command
				condor_command = "csub " + bash_script_path
				files_to_transfer = []
				files_to_transfer.append(bash_script_path)
				files_to_transfer.append(limit_config.get_datacard_filename(analysis, model, mass, args.fit_function, correctTrigger=args.correctTrigger, qcd=args.qcd, fitTrigger=args.fitTrigger))
				files_to_transfer.append(limit_config.get_workspace_filename(analysis, model, mass, correctTrigger=args.correctTrigger, qcd=args.qcd, fitTrigger=args.fitTrigger))
				condor_command += " -F " + ",".join(files_to_transfer)
				condor_command += " -l combine_{}_\$\(Cluster\)_\$\(Process\)".format(job_name)
				condor_command += " -s submit_combine_{}.jdl".format(job_name)
				condor_command += " -d " + submission_dir
				condor_command += " --cmssw"
				if not first or args.no_retar:
					condor_command += " --no_retar "
				print ">> Submitting job {}".format(job_name)
				print "Submission command: "
				print condor_command
				os.system(condor_command)

				print "---------------------------------------------------------------------------"
				os.chdir(start_dir)
				first = False
