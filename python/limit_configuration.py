import glob

paths = {}
paths["limits"] = "/uscms_data/d1/dryu/Dijets/EightTeeEeVeeBee/Fits/"
paths["datacards"] = paths["limits"] + "/Datacards/"
paths["condor"] = paths["limits"] + "/condor/"
paths["combine_logs"] = paths["limits"] + "/Logs/"
paths["resonance_shapes"] = "/uscms_data/d1/dryu/Dijets/EightTeeEeVeeBee/Fits/ResonanceShapes/"
paths["limit_plots"] = "/uscms_data/d1/dryu/Dijets/EightTeeEeVeeBee/Results/figures/"
paths["r_grid"] = "/uscms_data/d1/dryu/Dijets/EightTeeEeVeeBee/Fits/Grid/"

limit_signal_masses = {
	"trigbbl_CSVTM":range(400, 850, 50),
	"trigbbh_CSVTM":range(600, 1250, 50)
}

# Estimates of 2-sigma expected limits, in pb
limit_p2sigma_estimates = {
	"trigbbh_CSVTM":{
		"Hbb":{1200.0:0.018899, 1150.0:0.037124, 1100.0:0.039879, 1050.0:0.023405, 1000.0:0.025091, 950.0:0.025012, 900.0:0.02878, 850.0:0.029659, 800.0:0.036289, 750.0:0.042728, 700.0:0.057456, 650.0:0.102621, 600.0:0.303954},
		"RSG":{1200.0:0.019388, 1150.0:0.023826, 1100.0:0.035423, 1050.0:0.021946, 1000.0:0.024242, 950.0:0.02489, 900.0:0.022151, 850.0:0.030239, 800.0:0.037759, 750.0:0.042609, 700.0:0.058826, 650.0:0.09928, 600.0:0.295702, }
	},
	"trigbbl_CSVTM":{
		"Hbb":{800.0:0.044405, 750.0:0.04907, 700.0:0.055848, 650.0:0.065613, 600.0:0.081398, 550.0:0.114254, 500.0:0.135334, 450.0:(0.135334+0.627745)/2., 400.0:0.627745},
		"RSG":{800.0:0.044812, 750.0:0.049796, 700.0:0.056204, 650.0:0.065271, 600.0:0.079637, 550.0:0.117042, 500.0:0.14398, 450.0:(0.14398+0.598384)/2., 400.0:0.598384}
	}
}
limit_m2sigma_estimates = {
	"trigbbh_CSVTM":{
		"Hbb":{600.0:0.116381, 650.0:0.025939, 700.0:0.006439, 750.0:0.004601, 800.0:0.004001, 850.0:0.003405, 900.0:0.003926, 950.0:0.003535, 1000.0:0.004532, 1050.0:0.004533, 1100.0:0.006448, 1150.0:(0.006448+0.004828)/2., 1200.0:0.004828},
		"RSG":{600.0:0.091513, 650.0:0.013531, 700.0:0.006568, 750.0:0.004622, 800.0:0.004485, 850.0:0.003476, 900.0:0.002364, 950.0:0.003509, 1000.0:0.004443, 1050.0:0.004149, 1100.0:(0.004149+0.004126)/2., 1150.0:0.004126, 1200.0:0.005845}
	},
	"trigbbl_CSVTM":{
		"Hbb":{400.0:0.175381, 450.0:(0.175381+0.015603)/2., 500.0:0.015603, 550.0:0.01208, 600.0:0.010857, 650.0:0.010387, 700.0:0.010275, 750.0:0.010221, 800.0:0.01033},
		"RSG":{400.0:0.166479, 450.0:(0.166479+0.016075)/2., 500.0:0.016075, 550.0:0.012205, 600.0:0.010532, 650.0:0.010202, 700.0:0.010257, 750.0:0.010426, 800.0:0.010378}
	}
}

# Get the path to a workspace. 
def get_workspace_filename(analysis_name, model, mass, fitBonly=False, fitSignal=False):
	path = paths["datacards"] + "/workspace_" + analysis_name + "_" + model + "_" + str(mass)
	if fitBonly:
		path += "_fitBonly"
	if fitSignal:
		path += "_fitSignal"
	path += ".root"
	return path
	
def get_datacard_filename(analysis_name, model, mass, fit_function, fitSignal=False):
	path = paths["datacards"] + "/datacard_" + analysis_name + "_" + model + "_" + str(mass) + "_" + fit_function
	if fitSignal:
		path += "_fitSignal"
	path += ".txt"
	return path

def get_combine_log_path(analysis_name, model, mass, fit_function, method, systematics=True, frozen_nps=None, what="limits"):
	postfix = ""
	if not systematics:
		postfix += "_noSyst"
	if frozen_nps:
		postfix += "_" + frozen_nps.replace(",", "_")
	path = paths["combine_logs"] + "/" + what + "_" + method + "_" + analysis_name + "_" + model + "_m" + str(mass) + postfix + "_" + fit_function + ".log"
	return path

# Example: limits_HybridNewGrid_trigbbh_CSVTM_Hbb_m1200_f3_exp2.log
def get_combine_log_path_grid(analysis_name, model, mass, fit_function, what, method="HybridNewGrid", systematics=True, frozen_nps=None):
	postfix = ""
	if not systematics:
		postfix += "_noSyst"
	if frozen_nps:
		postfix += "_" + frozen_nps.replace(",", "_")
	path = paths["combine_logs"] + "/limits_" + method + "_" + analysis_name + "_" + model + "_m" + str(mass) + postfix + "_" + fit_function + "_" + what + ".log"
	return path

def get_data_input(analysis):
	path = paths["resonance_shapes"] + "/tyler/databb"
	if analysis == "trigbbh_CSVTM":
		path += "hTM.root"
	elif analysis == "trigbbl_CSVTM":
		path += "lTM.root"
	else:
		raise ValueError("[get_data_input] ERROR : Unknown analysis " + analysis)
	return path

def get_resonance_shapes(analysis, model):
	path = paths["resonance_shapes"] + "/tyler/ResonanceShapes_"
	if analysis == "trigbbh_CSVTM":
		path += "trigbbh"
	elif analysis == "trigbbl_CSVTM":
		path += "trigbbl"
	else:
		raise ValueError("[get_resonance_shapes] ERROR : Unknown analysis " + analysis)
	path += "_aug1"
	if model == "Hbb":
		path += "GG"
	elif model == "RSG":
		path += "RS"
	else:
		raise ValueError("[get_resonance_shapes] ERROR : Unknown model " + model)
	path += ".root"
	return path

def get_hn_grid(analysis, model, mass, fit_function):
	path = paths["r_grid"] + "/grid_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function + ".root" 
	return path

def get_hn_grid_subfiles(analysis, model, mass, fit_function):
	pattern = paths["r_grid"] + "/grid_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function + "/higgsCombinegrid_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function + "_i*root"
	return glob.glob(pattern)

