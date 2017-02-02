from ROOT import *
import math

# trigbbl from SingleMu (sigmoid) and JetHT (online b-tag efficiency)
# trigbbh from 60/53 (doesn't really matter where from)

# sigmoid**N
#sigmoid_parameters = {
#	"trigbbl_CSVTM":{
# 		"trigeff_p0":[2.31467e+02, 9.65421e+00],
# 		"trigeff_p1":[9.72373e+00, 7.29044e+00],
# 		"trigeff_p2":[1.22108e-01, 1.12468e-01],
# 		"trigeff_p3":[1.69494e-01, 4.45973e-03],
#	},
#	"trigbbh_CSVTM":{
#   		"trigeff_p0":[3.55647e+02, 2.84631e+01],
#   		"trigeff_p1":[2.31348e+01, 3.88590e+00],
#   		"trigeff_p2":[8.37131e-01, 7.61813e-01],
#   		"trigeff_p3":[2.76534e+00, 1.54388e-02],
#	}, 
#}
# Sigmoid
sigmoid_parameters = {
	"trigbbl_CSVTM":{
 		"trigeff_p0":[1.90191e+02, 4.68756e+00],
 		"trigeff_p1":[1.71505e+01, 4.52246e+00],
	},
	"trigbbh_CSVTM":{
		"trigeff_p0":[3.50048e+02, 5.33442e-01],
		"trigeff_p1":[2.38416e+01, 7.49631e-01],
	}, 
}
online_btag_eff = {
	"trigbbl_CSVTM":[1.82719e-01, 3.99115e-03],
	"trigbbh_CSVTM":[4.89446e-01, 5.91335e-03],
}

# From 80/70 / 60/53
#sigmoid_parameters = {
#	"trigbbl_CSVTM":{
#		"p0":1.10510e+02,
#		"p1":3.15726e+01,
#		"p2":6.72139e+00,
#	}, 
#	"trigbbh_CSVTM":{
#		"p0":4.03206e+02,
#		"p1":2.74656e+01,
#		"p2":1.53865e-01,
#	}, 
#}
## Old values
#sigmoid_parameters = {
#	"trigbbl_CSVTM":{
#		"p0":1.82469e+02,
#		"p1":2.87768e+01,
#		"p2":9.11659e-01,
#	}, 
#	"trigbbh_CSVTM":{
#		"p0":3.61785e+02,
#		"p1":3.16523e+01,
#		"p2":4.84357e-01,
#	}, 
#}


def get_pdf(analysis, mjj_var, name_tag=""):
	if analysis == "trigbbl_CSVTM":
		return RooGenericPdf("trigger_efficiency_trigbbl_CSVTM" + name_tag, 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbl_CSVTM"]["p0"][0], sigmoid_parameters["trigbbl_CSVTM"]["p1"][0], sigmoid_parameters["trigbbl_CSVTM"]["p2"][0]), 
		RooArgList(mjj_var))
	elif analysis == "trigbbh_CSVTM":
		return RooGenericPdf("trigger_efficiency_trigbbh_CSVTM" + name_tag, 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbh_CSVTM"]["p0"][0], sigmoid_parameters["trigbbh_CSVTM"]["p1"][0],sigmoid_parameters["trigbbh_CSVTM"]["p2"][0]), 
		RooArgList(mjj_var))
	elif analysis == "trigbbl_CSVM":
		return RooGenericPdf("trigger_efficiency_trigbbl_CSVM" + name_tag, 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbl_CSVTM"]["p0"][0], sigmoid_parameters["trigbbl_CSVTM"]["p1"][0], sigmoid_parameters["trigbbl_CSVTM"]["p2"][0]), 
		RooArgList(mjj_var))
	elif analysis == "trigbbh_CSVM":
		return RooGenericPdf("trigger_efficiency_trigbbh_CSVM" + name_tag, 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbh_CSVTM"]["p0"][0], sigmoid_parameters["trigbbh_CSVTM"]["p1"][0],sigmoid_parameters["trigbbh_CSVTM"]["p2"][0]), 
		RooArgList(mjj_var))
	elif analysis == "trigbbll_CSVTM":
		# This is just the trigbbl trigger efficiency. Consider calculating this properly; hopefully this is close enough.
		return RooGenericPdf("trigger_efficiency_trigbbll_CSVTM" + name_tag, 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbl_CSVTM"]["p0"][0], sigmoid_parameters["trigbbl_CSVTM"]["p1"][0], sigmoid_parameters["trigbbl_CSVTM"]["p2"][0]), 
		RooArgList(mjj_var))
	elif analysis == "trigbbhl_CSVTM":
		return RooGenericPdf("trigger_efficiency_trigbbhl_CSVTM" + name_tag, 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbh_CSVTM"]["p0"][0], sigmoid_parameters["trigbbh_CSVTM"]["p1"][0],sigmoid_parameters["trigbbh_CSVTM"]["p2"][0]), 
		RooArgList(mjj_var))

def get_formula(analysis, mjj_var):
	if "bbl" in analysis:
		return RooFormulaVar("trigger_efficiency_bbl", 
			"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbl_CSVTM"]["p0"][0], sigmoid_parameters["trigbbl_CSVTM"]["p1"][0], sigmoid_parameters["trigbbl_CSVTM"]["p2"][0]), 
			RooArgList(mjj_var))
	elif "bbh" in analysis:
		return RooFormulaVar("trigger_efficiency_trigbbh_CSVTM", 
			"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbh_CSVTM"]["p0"][0], sigmoid_parameters["trigbbh_CSVTM"]["p1"][0],sigmoid_parameters["trigbbh_CSVTM"]["p2"][0]), 
			RooArgList(mjj_var))

def get_var_formula(analysis, mjj_var):
	if "bbl" in analysis:
		eff_analysis = "trigbbl_CSVTM"
	elif "bbh" in analysis:
		eff_analysis = "trigbbh_CSVTM"

	eff_vars = {}
	for varname in ["trigeff_p0", "trigeff_p1"]:
		eff_vars[varname] = RooRealVar(varname, varname, sigmoid_parameters[eff_analysis][varname][0], sigmoid_parameters[eff_analysis][varname][0]-10.*sigmoid_parameters[eff_analysis][varname][1], sigmoid_parameters[eff_analysis][varname][0]+10.*sigmoid_parameters[eff_analysis][varname][1])
	return [
				RooFormulaVar("trigger_efficiency_bb", 
					"1. / (1. + exp(-1. * (@0 - @1) / @2))", 
					RooArgList(mjj_var, eff_vars["trigeff_p0"], eff_vars["trigeff_p1"])), 
				eff_vars
			]


def get_formula_with_btag(analysis, mjj_var):
	if "bbl" in analysis:
		return RooFormulaVar("trigger_efficiency_bbl", 
			"pow(%.1f / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(online_btag_eff["trigbbl_CSVTM"], sigmoid_parameters["trigbbl_CSVTM"]["p0"][0], sigmoid_parameters["trigbbl_CSVTM"]["p1"][0], sigmoid_parameters["trigbbl_CSVTM"]["p2"][0]), 
			RooArgList(mjj_var))
	elif "bbh" in analysis:
		return RooFormulaVar("trigger_efficiency_trigbbh_CSVTM", 
			"pow($.1f / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(online_btag_eff["trigbbh_CSVTM"], sigmoid_parameters["trigbbh_CSVTM"]["p0"][0], sigmoid_parameters["trigbbh_CSVTM"]["p1"][0],sigmoid_parameters["trigbbh_CSVTM"]["p2"][0]), 
			RooArgList(mjj_var))


def get_trivial_pdf(mjj_var):
	return RooGenericPdf("trigger_efficiency_dummy", "1.", RooArgList(mjj_var))

def get_trivial_formula(mjj_var):
	return RooFormulaVar("trigger_efficiency_dummy", "1.", RooArgList(mjj_var))

def trigger_efficiency_bbl(mjj):
	return math.pow(1. / (1. + math.exp(-1. * (mjj - sigmoid_parameters["trigbbl_CSVTM"]["p0"][0]) / sigmoid_parameters["trigbbl_CSVTM"]["p1"][0])), sigmoid_parameters["trigbbl_CSVTM"]["p2"][0]) 

def trigger_efficiency_bbh(mjj):
	return math.pow(1. / (1. + math.exp(-1. * (mjj - sigmoid_parameters["trigbbh_CSVTM"]["p0"][0]) / sigmoid_parameters["trigbbh_CSVTM"]["p1"][0])), sigmoid_parameters["trigbbh_CSVTM"]["p2"][0])

