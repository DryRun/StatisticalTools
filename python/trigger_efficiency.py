from ROOT import *
import math

def get_pdf(analysis, mjj_var):
	if analysis == "trigbbl_CSVTM":
		return RooGenericPdf("trigger_efficiency_trigbbl_CSVTM", 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(1.82469e+02, 2.87768e+01, 9.11659e-01), 
		RooArgList(mjj_var))
	elif analysis == "trigbbh_CSVTM":
		return RooGenericPdf("trigger_efficiency_trigbbh_CSVTM", 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(3.61785e+02, 3.16523e+01,4.84357e-01), 
		RooArgList(mjj_var))
	elif analysis == "trigbbl_CSVM":
		return RooGenericPdf("trigger_efficiency_trigbbl_CSVM", 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(1.82469e+02, 2.87768e+01, 9.11659e-01), 
		RooArgList(mjj_var))
	elif analysis == "trigbbh_CSVM":
		return RooGenericPdf("trigger_efficiency_trigbbh_CSVM", 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(3.61785e+02, 3.16523e+01,4.84357e-01), 
		RooArgList(mjj_var))
	elif analysis == "trigbbll_CSVTM":
		# This is just the trigbbl trigger efficiency. Consider calculating this properly; hopefully this is close enough.
		return RooGenericPdf("trigger_efficiency_trigbbll_CSVTM", 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(1.82469e+02, 2.87768e+01, 9.11659e-01), 
		RooArgList(mjj_var))
	elif analysis == "trigbbhl_CSVTM":
		return RooGenericPdf("trigger_efficiency_trigbbhl_CSVTM", 
		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(3.61785e+02, 3.16523e+01,4.84357e-01), 
		RooArgList(mjj_var))

def get_formula(analysis, mjj_var):
	if "bbl" in analysis:
		return RooFormulaVar("trigger_efficiency_bbl", 
			"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(1.82469e+02, 2.87768e+01, 9.11659e-01), 
			RooArgList(mjj_var))
	elif "bbh" in analysis:
		return RooFormulaVar("trigger_efficiency_trigbbhl_CSVTM", 
			"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(3.61785e+02, 3.16523e+01,4.84357e-01), 
			RooArgList(mjj_var))


def get_trivial_pdf(mjj_var):
	return RooGenericPdf("trigger_efficiency_dummy", "1.", RooArgList(mjj_var))

def get_trivial_formula(mjj_var):
	return RooFormulaVar("trigger_efficiency_dummy", "1.", RooArgList(mjj_var))

def trigger_efficiency_bbl(mjj):
	return math.pow(1. / (1. + math.exp(-1. * (mjj - 1.82469e+02) / 2.87768e+01)), 9.11659e-01) 

def trigger_efficiency_bbh(mjj):
	return math.pow(1. / (1. + math.exp(-1. * (mjj - 3.61785e+02) / 3.16523e+01)), 4.84357e-01)

