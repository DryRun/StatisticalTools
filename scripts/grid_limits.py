import sys, os, re
from argparse import ArgumentParser
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
from CMSDIJET.StatisticalTools.systematics import *
from CMSDIJET.StatisticalTools.roofit_functions import *

def limits(analysis, model, mass, fit_function, noSyst=False, freezeNuisances=None):
	print "\nCalculating limits for {}, {}, {}, {}".format(analysis, model, mass, fit_function)
	prefix = 'limits'

	postfix = ""
	if noSyst:
		postfix += "_noSyst"
	if freezeNuisances:
		postfix += "_" + args.freezeNuisances.replace(",", "_")

	datacard = limit_config.get_datacard_filename(analysis, model, mass, fit_function, fitSignal=True)
	log_base = limit_config.get_combine_log_path_grid(analysis, model, mass, fit_function, "HybridNewGrid", systematics=(not noSyst), frozen_nps=freezeNuisances)
	command_base = "combine {} -M HybridNew --frequentist --grid={}".format(datacard, limit_config.get_hn_grid(analysis, model, mass, fit_function))

	# Observed
	#log_observed = log_base.replace(".log", "_obs.log")
	print "Observed:"
	log_observed = limit_config.get_combine_log_path_grid(analysis, model, mass, fit_function, "obs", method="HybridNewGrid", systematics=(not noSyst), frozen_nps=freezeNuisances)

	command_observed = command_base + " 2>&1 | tee {}".format(log_observed)
	os.system(command_observed)

	# Expected
	expected_r = {-2:0.025, -1:0.16, 0:0.5, 1:0.84, 2:0.975}
	for exp in [-2, -1, 0, 1, 2]:
		print "Expected, " + str(exp)
		#log_expected = log_base.replace(".log", "_exp{}.log".format(exp))
		log_expected = limit_config.get_combine_log_path_grid(analysis, model, mass, fit_function, "exp" + str(exp), method="HybridNewGrid", systematics=(not noSyst), frozen_nps=freezeNuisances)
		print "[debug] Writing log file to " + log_expected
		command_expected = command_base + " --expectedFromGrid {} 2>&1 | tee {}".format(expected_r[exp], log_expected)
		os.system(command_expected)

if __name__ == "__main__":
	from joblib import Parallel, delayed
	from argparse import ArgumentParser

	# input parameters
	parser = ArgumentParser(description='Run combine over r-grid files to calculate the HybridNew limits')
	parser.add_argument('--analyses', type=str, default="trigbbh_CSVTM,trigbbl_CSVTM", help="Analysis names")
	parser.add_argument('--models', type=str, default="Hbb,RSG", help="Model names")
	parser.add_argument('--fit_function', type=str, default="f3", help="Name of central fit function")
	parser.add_argument('--mass', type=str, help="Run a single mass")
	args = parser.parse_args()

	analyses = args.analyses.split(",")
	models = args.models.split(",")

	for analysis in analyses:
		if args.mass:
			Parallel(n_jobs=4)(delayed(limits)(analysis, model, args.mass, args.fit_function) for model in models)
		else:
			Parallel(n_jobs=4)(delayed(limits)(analysis, model, mass, args.fit_function) for mass in limit_config.limit_signal_masses[analysis] for model in models)
