from ROOT import *
import math

sigmoid_parameters = {
	"trigbbl_CSVTM":{
		"p0":1.10510e+02,
		"p1":3.15726e+01,
		"p2":6.72139e+00,
	}, 
	"trigbbh_CSVTM":{
		"p0":4.03206e+02,
		"p1":2.74656e+01,
		"p2":1.53865e-01,
	}, 
}
# Old values
sigmoid_parameters = {
	"trigbbl_CSVTM":{
		"p0":1.82469e+02,
		"p1":2.87768e+01,
		"p2":9.11659e-01,
	}, 
	"trigbbh_CSVTM":{
		"p0":3.61785e+02,
		"p1":3.16523e+01,
		"p2":4.84357e-01,
	}, 
}


def get_pdf(analysis, mjj_var):
	if analysis == "trigbbl_CSVTM":
		return RooGenericPdf("trigger_efficiency_trigbbl_CSVTM", 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbl_CSVTM"]["p0"], sigmoid_parameters["trigbbl_CSVTM"]["p1"], sigmoid_parameters["trigbbl_CSVTM"]["p2"]), 
		RooArgList(mjj_var))
	elif analysis == "trigbbh_CSVTM":
		return RooGenericPdf("trigger_efficiency_trigbbh_CSVTM", 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbh_CSVTM"]["p0"], sigmoid_parameters["trigbbh_CSVTM"]["p1"],sigmoid_parameters["trigbbh_CSVTM"]["p2"]), 
		RooArgList(mjj_var))
	elif analysis == "trigbbl_CSVM":
		return RooGenericPdf("trigger_efficiency_trigbbl_CSVM", 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbl_CSVTM"]["p0"], sigmoid_parameters["trigbbl_CSVTM"]["p1"], sigmoid_parameters["trigbbl_CSVTM"]["p2"]), 
		RooArgList(mjj_var))
	elif analysis == "trigbbh_CSVM":
		return RooGenericPdf("trigger_efficiency_trigbbh_CSVM", 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbh_CSVTM"]["p0"], sigmoid_parameters["trigbbh_CSVTM"]["p1"],sigmoid_parameters["trigbbh_CSVTM"]["p2"]), 
		RooArgList(mjj_var))
	elif analysis == "trigbbll_CSVTM":
		# This is just the trigbbl trigger efficiency. Consider calculating this properly; hopefully this is close enough.
		return RooGenericPdf("trigger_efficiency_trigbbll_CSVTM", 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbl_CSVTM"]["p0"], sigmoid_parameters["trigbbl_CSVTM"]["p1"], sigmoid_parameters["trigbbl_CSVTM"]["p2"]), 
		RooArgList(mjj_var))
	elif analysis == "trigbbhl_CSVTM":
		return RooGenericPdf("trigger_efficiency_trigbbhl_CSVTM", 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbh_CSVTM"]["p0"], sigmoid_parameters["trigbbh_CSVTM"]["p1"],sigmoid_parameters["trigbbh_CSVTM"]["p2"]), 
		RooArgList(mjj_var))

def get_formula(analysis, mjj_var):
	if "bbl" in analysis:
		return RooFormulaVar("trigger_efficiency_bbl", 
			"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbl_CSVTM"]["p0"], sigmoid_parameters["trigbbl_CSVTM"]["p1"], sigmoid_parameters["trigbbl_CSVTM"]["p2"]), 
			RooArgList(mjj_var))
	elif "bbh" in analysis:
		return RooFormulaVar("trigger_efficiency_trigbbhl_CSVTM", 
			"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(sigmoid_parameters["trigbbh_CSVTM"]["p0"], sigmoid_parameters["trigbbh_CSVTM"]["p1"],sigmoid_parameters["trigbbh_CSVTM"]["p2"]), 
			RooArgList(mjj_var))


def get_trivial_pdf(mjj_var):
	return RooGenericPdf("trigger_efficiency_dummy", "1.", RooArgList(mjj_var))

def get_trivial_formula(mjj_var):
	return RooFormulaVar("trigger_efficiency_dummy", "1.", RooArgList(mjj_var))

def trigger_efficiency_bbl(mjj):
	return math.pow(1. / (1. + math.exp(-1. * (mjj - sigmoid_parameters["trigbbl_CSVTM"]["p0"]) / sigmoid_parameters["trigbbl_CSVTM"]["p1"])), sigmoid_parameters["trigbbl_CSVTM"]["p2"]) 

def trigger_efficiency_bbh(mjj):
	return math.pow(1. / (1. + math.exp(-1. * (mjj - sigmoid_parameters["trigbbh_CSVTM"]["p0"]) / sigmoid_parameters["trigbbh_CSVTM"]["p1"])), sigmoid_parameters["trigbbh_CSVTM"]["p2"])

