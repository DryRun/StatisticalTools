#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
from CMSDIJET.StatisticalTools.systematics import *
from CMSDIJET.StatisticalTools.roofit_functions import *

from ROOT import *
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
seaborn = Root.SeabornInterface()
seaborn.Initialize()

titles = {"KS":"d_{KS}", "AD":"A^{2}", "saturated":""}

if __name__ == "__main__":
	# input parameters
	parser = ArgumentParser(description='Run goodness of fit tests using combine')
	parser.add_argument('--analyses', type=str, default="trigbbh_CSVTM,trigbbl_CSVTM", help="Analysis names")
	parser.add_argument('--models', type=str, default="Hbb", help="Model names")
	parser.add_argument('--qcd', action='store_true', help="Use QCD instead of data (assumes no trigger emulation)")
	parser.add_argument('--correctTrigger', action='store_true', help="Use model with trigger correction (has to have been specified in create_datacards.py)")
	parser.add_argument('--fitTrigger', action='store_true', help="Use model with trigger fit (has to have been specified in create_datacards.py)")
	parser.add_argument('--fit_function', type=str, default="dijet4", help="Name of central fit function")
	parser.add_argument('--algo', type=str, default="AD", help="GoF statistic: AD, saturated, KS")
	parser.add_argument('--toys', type=int, help="Number of toys")
	parser.add_argument('--seed', type=int, default=123456, help="Random seed")
	parser.add_argument('--cplots', action="store_true", help="Make combine plots (not sure if this does anything)")
	parser.add_argument('--no_signal', action='store_true', help="Fix signal strength to zero")
	parser.add_argument('--signal', type=float, help="Fixed signal strength")
	parser.add_argument('-m', '--masses', type=str, help='Manually specify masses (comma-separated list). Otherwise, taken from limit_configuration.')
	parser.add_argument('--no_retar', action='store_true', help='Don\'t retar CMSSW directory before submitting')
	parser.add_argument('--fit', action='store_true', help='Run GoF tests')
	parser.add_argument('--plot', action='store_true', help='GoF plots')
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
				if args.no_signal:
					job_name += "_mu0"

				log_name = job_name + ".log"
				
				if args.fit:
					combine_options =  "-M GoodnessOfFit -v3 --algo {} -s {} --name {} --mass {} ".format(
						args.algo,
						args.seed,
						job_name,
						mass
					)
					if args.cplots:
						combine_options += " --plots "
					if args.no_signal:
						combine_options += " --fixedSignalStrength 0"
					elif args.signal != None:
						combine_options += " --fixedSignalStrength {} ".format(args.signal)

					cmd = "combine {} {} 2>&1 | tee {}".format(
						combine_options, 
						os.path.basename(limit_config.get_datacard_filename(analysis, model, mass, args.fit_function, correctTrigger=args.correctTrigger, qcd=args.qcd, fitTrigger=args.fitTrigger)),
						log_name)

					if args.toys:
						combine_toy_options =  "-M GoodnessOfFit -v3 --algo {} -s {} --name {} --mass {} ".format(
							args.algo,
							args.seed,
							job_name + "_toys",
							mass
						)
						if args.no_signal:
							combine_toy_options += " --fixedSignalStrength 0"
						elif args.signal != None:
							combine_toy_options += " --fixedSignalStrength {} ".format(args.signal)

						toy_cmd = "combine {} {} -t {} 2>&1 | tee {}".format(
							combine_toy_options, 
							os.path.basename(limit_config.get_datacard_filename(analysis, model, mass, args.fit_function, correctTrigger=args.correctTrigger, qcd=args.qcd, fitTrigger=args.fitTrigger)),
							args.toys,
							job_name + "_toys.log")

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

				if args.plot:
					# higgsCombinegof_trigbbh_CSVTM_Hbb_m750_f1_AD.GoodnessOfFit.mH750.root
					# higgsCombinegof_trigbbh_CSVTM_Hbb_m750_f1_AD_toys.GoodnessOfFit.mH750.123456.root
					data_file = TFile(limit_config.paths["gof"] + "/higgsCombine{}.GoodnessOfFit.mH{}.root".format(job_name, mass), "READ")
					toy_file = TFile(limit_config.paths["gof"] + "/higgsCombine{}_toys.GoodnessOfFit.mH{}.{}.root".format(job_name, mass, args.seed), "READ")
					data_tree = data_file.Get("limit")
					data_tree.GetEntry(0)
					data_gof = data_tree.GetLeaf("limit").GetValue(0)
					toy_tree = toy_file.Get("limit")
					toy_gofs = []
					toy_gof_hist = TH1D("toy_gof_" + job_name, "toy_gof_" + job_name, 50, 0., 10.)
					total_toys = 0
					gt_toys = 0
					for entry in xrange(toy_tree.GetEntriesFast()):
						toy_tree.GetEntry(entry)
						this_gof = toy_tree.GetLeaf("limit").GetValue(0)
						toy_gofs.append(this_gof)
						toy_gof_hist.Fill(this_gof)
						if this_gof >= data_gof:
							gt_toys += 1
						total_toys += 1

					c = TCanvas("c_gof_" + job_name, "c_gof_" + job_name, 800, 600)
					x_min = 0.
					x_max = max(max(toy_gofs), data_gof) * 1.2
					frame =  TH1D("frame_" + job_name, "frame_" + job_name, 100, x_min, x_max)
					frame.GetXaxis().SetTitle(titles[args.algo])
					frame.SetMinimum(0.)
					frame.SetMaximum(toy_gof_hist.GetMaximum()*1.3)
					frame.Draw("axis")
					toy_gof_hist.Draw("hist same")
					data_gof_line = TLine(data_gof, frame.GetMinimum(), data_gof, frame.GetMaximum())
					data_gof_line.SetLineStyle(1)
					data_gof_line.SetLineWidth(3)
					data_gof_line.SetLineColor(seaborn.GetColorRoot("default", 2))
					data_gof_line.Draw("same")

					Root.myText(0.6, 0.7, kBlack, "p={}".format(round(1.*gt_toys/total_toys, 2)), 0.5)
					Root.myText(0.6, 0.8, kBlack, "Data GOF=".format(round(data_gof, 2)), 0.5)

					c.SaveAs(limit_config.paths["gof"] + "/figures/" + c.GetName() + ".pdf")
