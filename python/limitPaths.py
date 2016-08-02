paths = {}
paths["limits"] = "/uscms_data/d1/dryu/Dijets/EightTeeEeVeeBee/Limits/"
paths["fits"] = "/uscms_data/d1/dryu/Dijets/EightTeeEeVeeBee/Fits/"
paths["datacards"] = paths["limits"] + "/datacards/"
paths["workspaces"] = paths["fits"] + "/workspaces/"
#paths["condor"] = paths["fits"] + "/condor/"
paths["resonance_shapes"] = "/uscms_data/d1/dryu/Dijets/EightTeeEeVeeBee/Fits/ResonanceShapes/"

# Get the path to a workspace. 
def get_workspace_filename(analysis_name, model):
	return paths["workspaces"] + "/workspace_" + analysis_name + "_" + model + ".root"