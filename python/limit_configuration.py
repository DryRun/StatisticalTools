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

def get_combine_log_path(analysis_name, model, mass, fit_function, method, systematics=True, frozen_nps=None):
	postfix = ""
	if not systematics:
		postfix += "_noSyst"
	if frozen_nps:
		postfix += "_" + frozen_nps.replace(",", "_")
	path = paths["combine_logs"] + "/limits_" + method + "_" + analysis_name + "_" + model + "_m" + str(mass) + postfix + "_" + fit_function + ".log"
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

