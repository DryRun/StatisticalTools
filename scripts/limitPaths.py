paths = {}
paths["limits"] = "/uscms_data/d1/dryu/Dijets/EightTeeEeVeeBee/Results/Limits/"
paths["datacards"] = paths["limits"] + "/datacards/"
paths["workspaces"] = paths["limits"] + "/workspaces/"
paths["condor"] = paths["limits"] + "/condor/"
paths["resonance_shapes"] = "/uscms_data/d1/dryu/Dijets/EightTeeEeVeeBee/Results/ResonanceShapes/"

# Get the path to a workspace. 
def GetWorkspacePath(analysis_name, model):
	return paths["workspaces"] + "/workspace_" + analysis_name + "_" + model + ".root"