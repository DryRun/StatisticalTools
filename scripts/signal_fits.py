import os
import sys
sys.stdout.flush()
import ROOT
import copy
from ROOT import *
gROOT.SetBatch(True)
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/RooCBPlusVoigtian.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so"))
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
seaborn = Root.SeabornInterface()
seaborn.Initialize()
from CMSDIJET.StatisticalTools.systematics import *
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
import CMSDIJET.StatisticalTools.trigger_efficiency as trigger_efficiency

sys.path.append(os.path.expandvars("$CMSSW_BASE/python/CMSDIJET/StatisticalTools"))
import analysis_configuration_8TeV as analysis_config

from array import array
dijet_binning = array("d", [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8000])
dijet_roobinning = RooBinning(len(dijet_binning) - 1, dijet_binning, "dijet_binning")

def make_signal_pdf(fit_function, mjj, tag=None, mass=None):
	if fit_function == "voigt":
		signal_vars = {
			"mean":RooRealVar("mean", "mean", mass, 0., 1530),
			"width":RooRealVar("width", "width", 100., 0., 1000.),
			"sigma":RooRealVar("sigma", "sigma", 100., 0., 1000.),
		}
		signal_pdf = ROOT.RooVoigtian("voigt", "voigt", mjj, signal_vars["mean"], signal_vars["width"], signal_vars["sigma"])
	elif fit_function == "bukin":
		# Bukin
		signal_vars = {
			"xp":RooRealVar("xp", "xp", mass, mass - 200., mass + 200.),
			"sigp":RooRealVar("sigp", "sigp", 100., 15., 200.),
			"xi":RooRealVar("xi", "xi", -0.25, -1.5, -1.e-3),
			"rho1":RooRealVar("rho1", "rho1", 0.5, 0., 10.),
			"rho2":RooRealVar("rho2", "rho2", 0.5, 0.01, 10.),
		}
		signal_vars["rho1"].setVal(0.)
		signal_vars["rho1"].setConstant()
		signal_pdf = ROOT.RooBukinPdf("bukin", "bukin", mjj, signal_vars["xp"], signal_vars["sigp"], signal_vars["xi"], signal_vars["rho1"], signal_vars["rho2"])
	elif fit_function == "cb":
		# Crystal Ball
		signal_vars = {
			"mean":RooRealVar("mean", "mean", mass, 0., 1530.),
			"sigma":RooRealVar("sigma", "sigma", 100., 0., 1000.),
			"alpha":RooRealVar("alpha", "alpha", 1., -20., 20.),
			"n":RooRealVar("n", "n", 1., 0., 100.),
		}
		signal_pdf = RooCBShape("cb", "cb", mjj, signal_vars["mean"], signal_vars["sigma"], signal_vars["alpha"], signal_vars["n"])
	elif fit_function == "cb_voigt":
		signal_vars = {
			"mean":RooRealVar("mean", "mean", mass, 0., 1530),
			"voigt_width":RooRealVar("voigt_width", "voigt_width", 100., 0., 1000.),
			"voigt_sigma":RooRealVar("voigt_sigma", "voigt_sigma", 100., 0., 1000.),
			"cb_sigma":RooRealVar("cb_sigma", "cb_sigma", 100., 0., 1000.),
			"cb_alpha":RooRealVar("cb_alpha", "cb_alpha", 1., -20., 20.),
			"cb_n":RooRealVar("cb_n", "cb_n", 1., 0., 100.),
			"f_cb":RooRealVar("f_cb", "f_cb", 0.5, 0., 1.),
		}
		signal_pdf = ROOT.RooCBPlusVoigtian("cb_voigt", "cb_voigt", mjj, signal_vars["mean"], signal_vars["voigt_width"], signal_vars["voigt_sigma"], signal_vars["cb_sigma"], signal_vars["cb_alpha"], signal_vars["cb_n"], signal_vars["f_cb"])
	else:
		print "[make_signal_pdf] ERROR : Fit function " + fit_function + " not known."
		sys.exit(1)

	if tag:
		signal_pdf.SetName(signal_pdf.GetName() + tag)
		for var_name, var in signal_vars.iteritems():
			var.SetName(var.GetName() + tag)
	return signal_pdf, signal_vars

# Make signal PDFs for JER/JES/BOff variations, with nuisance parameter built in. 
def make_signal_pdf_systematic(fit_function, mjj, tag=None, mass=None):
	if fit_function == "bukin":
		# Bukin
		signal_vars = {
			"xp_0":RooRealVar("xp_0", "xp_0", mass, mass - 200., mass + 200.),
			"sigp_0":RooRealVar("sigp_0", "sigp_0", 100., 15., 200.),
			"xi_0":RooRealVar("xi_0", "xi_0", 0., -1.5, -1.e-3),
			"rho1_0":RooRealVar("rho1_0", "rho1_0", 0.5, 0., 10.),
			"rho2_0":RooRealVar("rho2_0", "rho2_0", 0.5, 1.e-3, 10.),
			"alpha_jes":RooRealVar("alpha_jes", "alpha_jes", 0., -4., 4.),
			"alpha_jer":RooRealVar("alpha_jer", "alpha_jer", 0., -4., 4.),
			"dsigp":RooRealVar("dsigp", "dsigp", 0., -1.e10, 1.e10),
			"dxp":RooRealVar("dxp", "dxp", 0., -1.e10, 1.e10),
		}
		signal_vars["rho1_0"].setVal(0.)
		signal_vars["rho1_0"].setConstant()
		signal_vars["xp"] = ROOT.RooFormulaVar("xp", "xp", "xp_0 + (alpha_jes * dxp)", RooArgList(signal_vars["xp_0"], signal_vars["alpha_jes"], signal_vars["dxp"]))
		signal_vars["sigp"] = ROOT.RooFormulaVar("sigp", "sigp", "sigp_0 + (alpha_jer * dsigp)", RooArgList(signal_vars["sigp_0"], signal_vars["alpha_jer"], signal_vars["dsigp"]))
		signal_pdf = ROOT.RooBukinPdf("bukin", "bukin", mjj, signal_vars["xp"], signal_vars["sigp"], signal_vars["xi_0"], signal_vars["rho1_0"], signal_vars["rho2_0"])
	else:
		print "[make_signal_pdf_systematic] ERROR : Fit function " + fit_function + " not known."
		sys.exit(1)

	if tag:
		signal_pdf.SetName(signal_pdf.GetName() + tag)
		for var_name, var in signal_vars.iteritems():
			var.SetName(var.GetName() + tag)
	return signal_pdf, signal_vars


#systematics_floating_parameters = {}
#systematics_floating_parameters["BOffUp"] = {
#	"voigt":["mean"],
#	"bukin":["xp"],
#	"cb":["mean", "alpha"],
#	"cb_voigt":["mean", "cb_alpha"]
#}
#systematics_floating_parameters["BOffDown"] = copy.deepcopy(systematics_floating_parameters["BOffUp"])
#
#systematics_floating_parameters["JERUp"] = {
#	"voigt":["width", "sigma"],
#	"bukin":["sigp"],
#	"cb":["sigma", "alpha", "n"],
#	"cb_voigt":["voigt_width", "voigt_sigma", "cb_sigma", "cb_alpha", "cb_n"]
#}
#systematics_floating_parameters["JERDown"] = copy.deepcopy(systematics_floating_parameters["JERUp"])

def setfix_systematic_parameters(fit_function, systematic_name, central_parameters, systematic_parameters):
	systematics_floating_parameters = {}
	systematics_floating_parameters["JESUp"] = {
		"voigt":["mean"],
		"bukin":["xp"],
		"cb":["mean", "alpha"],
		"cb_voigt":["mean", "cb_alpha"]
	}
	systematics_floating_parameters["JESDown"] = copy.deepcopy(systematics_floating_parameters["JESUp"])

	systematics_floating_parameters["JERUp"] = {
		"voigt":["width", "sigma"],
		#"bukin":["sigp", "xi", "rho1", "rho2"],
		"bukin":["sigp"],
		"cb":["sigma", "alpha", "n"],
		"cb_voigt":["voigt_width", "voigt_sigma", "cb_sigma", "cb_alpha", "cb_n"]
	}
	systematics_floating_parameters["JERDown"] = copy.deepcopy(systematics_floating_parameters["JERUp"])

	systematics_floating_parameters["BOffUp"] = {
		"voigt":["width", "sigma"],
		#"bukin":["sigp", "xi", "rho1", "rho2"],
		"bukin":["sigp", "xp"],
		"cb":["sigma", "alpha", "n"],
		"cb_voigt":["voigt_width", "voigt_sigma", "cb_sigma", "cb_alpha", "cb_n"]
	}
	systematics_floating_parameters["BOffDown"] = copy.deepcopy(systematics_floating_parameters["JERUp"])

	# Add norm to all lists of floating parameters
	for name1, map1 in systematics_floating_parameters.iteritems():
		for name2, list2 in map1.iteritems():
			list2.append("norm")

	if not systematic_name in systematics_floating_parameters:
		raise ValueError("[setfix_systematic_parameters] ERROR : Systematic " + systematic_name + " not known.")

	modification_factors = {
		"JESUp":1.01,
		"JESDown":0.99,
		"JERUp":1.1,
		"JERDown":0.95,
		"BOffUp":1.,
		"BOffDown":1.
	}

	for parameter_name in central_parameters:
		systematic_parameters[parameter_name].setVal(central_parameters[parameter_name].getVal())
		if not parameter_name in systematics_floating_parameters[systematic_name][fit_function]:
			systematic_parameters[parameter_name].setConstant()
		else:
			systematic_parameters[parameter_name].setVal(central_parameters[parameter_name].getVal() * modification_factors[systematic_name])

def fix_parameters(pdf):
	parameter_argset = pdf.getVariables()
	iterator = parameter_argset.createIterator()
	parameter = iterator.Next()
	while parameter:
		parameter.setConstant()
		parameter = iterator.Next()

def get_parameters(pdf):
	parameters = {}
	parameter_argset = pdf.getVariables()
	iterator = parameter_argset.createIterator()
	parameter = iterator.Next()
	while parameter:
		name = parameter.GetName()
		value = parameter.getVal()
		error = parameter.getError()
		parameters[name] = [value, error]
		parameter = iterator.Next()
	return parameters

def copy_signal_pdf(fit_function, input_function, mjj, tag=None, include_systematics=False, mass=None):
	if include_systematics:
		copy_pdf, copy_vars = make_signal_pdf_systematic(fit_function, mjj, tag)
	else:
		copy_pdf, copy_vars = make_signal_pdf(fit_function, mjj, tag, mass=mass)
	input_vars = input_function.getVariables()
	input_var_values = {}
	iterator = input_vars.createIterator()
	input_var = iterator.Next()
	while input_var:
		name = input_var.GetName()
		if tag:
			name = name.replace(tag, "")
		if not "mjj" in name:
			input_var_values[name] = input_var.getVal()
		input_var = iterator.Next()

	for input_var_name, input_var_value in input_var_values.iteritems():
		#print copy_vars
		copy_vars[input_var_name].setVal(input_var_value)

	#if fit_function == "voigt":
	#	copy_vars["mean"].setVal(input_var_values["mean"])
	#	copy_vars["width"].setVal(input_var_values["width"])
	#	copy_vars["sigma"].setVal(input_var_values["sigma"])
	#elif fit_function == "bukin":
	#	copy_vars["xp"].setVal(input_var_values["xp"])
	#	copy_vars["sigp"].setVal(input_var_values["sigp"])
	#	copy_vars["xi"].setVal(input_var_values["xi"])
	#	copy_vars["rho1"].setVal(input_var_values["rho1"])
	#	copy_vars["rho2"].setVal(input_var_values["rho2"])
	#elif fit_function == "cb":
	#	copy_vars["mean"].setVal(input_var_values["mean"])
	#	copy_vars["sigma"].setVal(input_var_values["sigma"])
	#	copy_vars["alpha"].setVal(input_var_values["alpha"])
	#	copy_vars["n"].setVal(input_var_values["n"])
	#elif fit_function == "cb_voigt":
	#	print input_var_values
	#	copy_vars["mean"].setVal(input_var_values["mean"])
	#	copy_vars["voigt_width"].setVal(input_var_values["voigt_width"])
	#	copy_vars["voigt_sigma"].setVal(input_var_values["voigt_sigma"])
	#	copy_vars["cb_sigma"].setVal(input_var_values["cb_sigma"])
	#	copy_vars["cb_alpha"].setVal(input_var_values["cb_alpha"])
	#	copy_vars["cb_n"].setVal(input_var_values["cb_n"])
	#	copy_vars["f_cb"].setVal(input_var_values["f_cb"])
	#else:
	#	print "[copy_signal_pdf] ERROR : Fit function " + fit_function + " not known."
	#	sys.exit(1)
	return copy_pdf, copy_vars

def signal_fit(analysis, model, mass, fit_functions, systematics=None):
	print "[signal_fit] INFO : Welcome to signal_fit({}, {}, {}, {})".format(analysis, model, mass, ",".join(fit_functions))
	fit_trigger = not ("notrig" in analysis or "NoTrigger" in analysis)
	if systematics:
		print "\t",
		print systematics
	histogram_file = TFile(analysis_config.get_b_histogram_filename(analysis, analysis_config.simulation.get_signal_tag(model, mass, "FULLSIM")), "READ")
	histogram = histogram_file.Get("BHistograms/h_pfjet_mjj")
	#if trigger_strategy="CorrectMC":
	#	if "bbl" in analysis:
	#		for bin in xrange(1, histogram.GetNbinsX() + 1):
	#			bin_center = histogram.GetXaxis().GetBinCenter(bin)
	#			if bin_center > 240.:
	#				correction = 1. / trigger_efficiency.trigger_efficiency_bbl(bin_center)
	#			else:
	#				correction = 0.
	#			histogram.SetBinContent(bin, histogram.GetBinContent(bin) * correction)
	#			histogram.SetBinError(bin, histogram.GetBinError(bin) * correction)
	#	elif "bbh" in analysis:
	#		for bin in xrange(1, histogram.GetNbinsX() + 1):
	#			bin_center = histogram.GetXaxis().GetBinCenter(bin)
	#			if bin_center > 400.:
	#				correction = 1. / trigger_efficiency.trigger_efficiency_bbl(bin_center)
	#			else:
	#				correction = 0.
	#			histogram.SetBinContent(bin, histogram.GetBinContent(bin) * correction)
	#			histogram.SetBinError(bin, histogram.GetBinError(bin) * correction)

	# Make systematic variation histograms
	systematic_variations = []
	if systematics:
		print "[signal_fit] INFO : Fitting systematic variations: ",
		print systematics
		signal_cdf = TGraph(histogram.GetNbinsX()+1)
		histograms_syst = {}
		if "jer" in systematics or "jes" in systematics:
			signal_cdf.SetPoint(0,0.,0.)
			integral = 0.
			for i in range(1, histogram.GetNbinsX()+1):
				x = histogram.GetXaxis().GetBinLowEdge(i+1)
				integral = integral + histogram.GetBinContent(i)
				signal_cdf.SetPoint(i,x,integral)

		# Initialize histograms
		if "jes" in systematics:
			systematic_variations.append("JESUp")
			systematic_variations.append("JESDown")
			histograms_syst['JESUp'] = copy.deepcopy(histogram)
			histograms_syst['JESDown'] = copy.deepcopy(histogram)
		if "jer" in systematics:
			systematic_variations.append("JERUp")
			systematic_variations.append("JERDown")
			histograms_syst['JERUp'] = copy.deepcopy(histogram)
			histograms_syst['JERDown'] = copy.deepcopy(histogram)
		if "boff" in systematics:
			systematic_variations.append("BOffUp")
			systematic_variations.append("BOffDown")
			histograms_syst['BOffUp'] = copy.deepcopy(histogram)
			histograms_syst['BOffDown'] = copy.deepcopy(histogram)
		for key in histograms_syst.keys():
			histograms_syst[key].Reset()
			histograms_syst[key].SetName(histograms_syst[key].GetName() + '_' + key)

		# produce JES signal shapes
		if "jes" in systematics:
			for i in range(1, histogram.GetNbinsX()+1):
				xLow = histogram.GetXaxis().GetBinLowEdge(i)
				xUp = histogram.GetXaxis().GetBinLowEdge(i+1)
				jes = 1. - systematics["jes"]
				xLowPrime = jes*xLow
				xUpPrime = jes*xUp
				histograms_syst['JESUp'].SetBinContent(i, signal_cdf.Eval(xUpPrime) - signal_cdf.Eval(xLowPrime))
				jes = 1. + systematics["jes"]
				xLowPrime = jes*xLow
				xUpPrime = jes*xUp
				histograms_syst['JESDown'].SetBinContent(i, signal_cdf.Eval(xUpPrime) - signal_cdf.Eval(xLowPrime))
			
		# produce JER signal shapes
		if "jer" in systematics:
			for i in range(1, histogram.GetNbinsX()+1):
				xLow = histogram.GetXaxis().GetBinLowEdge(i)
				xUp = histogram.GetXaxis().GetBinLowEdge(i+1)
				jer = 1. - systematics["jer"]
				xLowPrime = jer*(xLow-float(mass))+float(mass)
				xUpPrime = jer*(xUp-float(mass))+float(mass)
				histograms_syst['JERUp'].SetBinContent(i, signal_cdf.Eval(xUpPrime) - signal_cdf.Eval(xLowPrime))
				jer = 1. + systematics["jer"]
				xLowPrime = jer*(xLow-float(mass))+float(mass)
				xUpPrime = jer*(xUp-float(mass))+float(mass)
				histograms_syst['JERDown'].SetBinContent(i, signal_cdf.Eval(xUpPrime) - signal_cdf.Eval(xLowPrime))

		# B tag SFs: histograms are already made
		if "boff" in systematics:
			systematic_variations.append("BOffUp")
			systematic_variations.append("BOffDown")
			histograms_syst['BOffUp'] = histogram_file.Get("BHistograms/h_pfjet_mjj_BTagOfflineSFUp")
			histograms_syst['BOffDown'] = histogram_file.Get("BHistograms/h_pfjet_mjj_BTagOfflineSFDown")


	# End making systematic variation histograms

	for fit_function in fit_functions:
		if "trigbbh" in analysis or "NoTrigger_eta2p2" in analysis:
			#mjj = RooRealVar('mjj','mjj',354, 1530)
			mjj = RooRealVar('mjj','mjj',489, 1945)
		elif "trigbbl" in analysis or "NoTrigger_eta1p7" in analysis:
			#mjj = RooRealVar('mjj','mjj',296, 1530)
			mjj = RooRealVar('mjj','mjj',270, 1945)
		elif "trigbbll" in analysis:
			#mjj = RooRealVar('mjj','mjj',296, 1530)
			mjj = RooRealVar('mjj','mjj',296, 1945)
		elif "trigbbhl" in analysis:
			#mjj = RooRealVar('mjj','mjj',296, 1530)
			mjj = RooRealVar('mjj','mjj',354, 1945)

		# Edit 10/23/2016: moved this to an independent module, to enforce consistency across scripts.
		#trigger_efficiency_pdfs = {}
		#trigger_efficiency_pdfs["trigbbl_CSVTM"] = RooGenericPdf("trigger_efficiency_trigbbl_CSVTM", 
		#		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(1.82469e+02, 2.87768e+01, 9.11659e-01), 
		#		RooArgList(mjj))
		#trigger_efficiency_pdfs["trigbbh_CSVTM"] = RooGenericPdf("trigger_efficiency_trigbbh_CSVTM", 
		#		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(3.61785e+02, 3.16523e+01,4.84357e-01), 
		#		RooArgList(mjj))
		#trigger_efficiency_pdfs["trigbbl_CSVM"] = RooGenericPdf("trigger_efficiency_trigbbl_CSVM", 
		#		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(1.82469e+02, 2.87768e+01, 9.11659e-01), 
		#		RooArgList(mjj))
		#trigger_efficiency_pdfs["trigbbh_CSVM"] = RooGenericPdf("trigger_efficiency_trigbbh_CSVM", 
		#		"pow(1. / (1. + exp(-1. * (@0 - %.1f) / %.1f)), %.1f)"%(3.61785e+02, 3.16523e+01,4.84357e-01), 
		#		RooArgList(mjj))

		signal_rdh = RooDataHist(histogram.GetName() + "_rdh", histogram.GetName() + "_rdh", ROOT.RooArgList(mjj), histogram)
		signal_pdf_raw, signal_vars = make_signal_pdf(fit_function, mjj, mass=mass)
		signal_pdf_raw.SetName("signal_bukin_pdf")
		if fit_trigger:
			#trigger_efficiency_pdf = trigger_efficiency.get_pdf(analysis, mjj)
			#signal_pdf = RooProdPdf("signal_pdf", "signal_pdf", signal_pdf_raw, trigger_efficiency_pdf)
			trigger_efficiency_formula = trigger_efficiency.get_formula(analysis, mjj)
			signal_pdf = RooEffProd("signal_pdf", "signal_pdf", signal_pdf_raw, trigger_efficiency_formula)
		else:
			signal_pdf = signal_pdf_raw
			signal_pdf.SetName("signal_pdf")

		# Some fits are special and need help
		if "trigbbl" in analysis or "NoTrigger_eta1p7" in analysis:
			#if (analysis == "trigbbl_CSVTM" or analysis == "trigbbll_CSVTM" or analysis == "trigbbl_notrig_CSVTM"):
			xi_initial = -0.1 + (-0.45 + 0.1) / (900 - 400) * (mass - 400)
			signal_vars["xi"].setMax(-0.05)
			signal_vars["xi"].setVal(xi_initial)
			#signal_vars["rho2"].setMin(0.04)
			#signal_vars["rho2"].setVal(0.08)
			rho2_initial = 0.03 + (0.135 - 0.03) / (900 - 400) * (mass - 400)
			signal_vars["rho2"].setMin(0.01)
			signal_vars["rho2"].setMax(0.2)
			signal_vars["rho2"].setVal(rho2_initial)

			sigp_initial = 35 + (100 -35) / (900 - 400) * (mass - 400)
			signal_vars["sigp"].setVal(sigp_initial)
			if mass == 325 and model == "Hbb":
				signal_vars["xp"].setVal(305)
				signal_vars["sigp"].setVal(34)
				signal_vars["xi"].setVal(-0.05)
				signal_vars["rho2"].setVal(0.03)
			elif mass == 350 and model == "ZPrime":
				signal_vars["xp"].setVal(325)
				signal_vars["sigp"].setVal(37)
				signal_vars["xi"].setVal(-0.05)
				signal_vars["rho2"].setVal(0.033)
			elif mass == 400 and model == "ZPrime":
				signal_vars["xp"].setVal(372)
				signal_vars["sigp"].setVal(46)
				signal_vars["xi"].setVal(-0.15)
				signal_vars["rho2"].setVal(0.042)
			elif mass == 500 and model == "RSG":
				signal_vars["xp"].setVal(466)
				signal_vars["sigp"].setVal(58)
				signal_vars["xi"].setVal(-0.24)
				signal_vars["rho2"].setVal(0.056)
			elif mass == 600:
				signal_vars["xp"].setVal(5.5398e+02)
				signal_vars["rho2"].setVal(7.4604e-02)
				signal_vars["xi"].setVal(-1.9649e-01)
				signal_vars["sigp"].setVal(7.2282e+01)
			elif mass == 400:
				#signal_vars["xi"].setMax(10.)
				signal_vars["xp"].setMin(300.)
			elif mass == 500:
				signal_vars["sigp"].setMin(45.)
				signal_vars["xi"].setMin(-0.5)
			elif mass == 750:
				signal_vars["sigp"].setVal(82.)
				signal_vars["sigp"].setMin(60.)
				signal_vars["sigp"].setMax(110.)
				signal_vars["xi"].setVal(-0.3)
				signal_vars["xi"].setMax(-0.1)
				signal_vars["xi"].setMin(-0.5)
				signal_vars["xp"].setVal(710.)
				signal_vars["xp"].setMin(600.)
				signal_vars["xp"].setMax(800.)
				signal_vars["rho2"].setVal(0.09)
				signal_vars["rho2"].setMin(0.01)
				signal_vars["rho2"].setMax(0.2)
			elif mass == 900:
				signal_vars["sigp"].setMin(70.)
			if model == "ZPrime" and mass == 400 and "CSVTM" in analysis:
				signal_vars["rho2"].setMin(0.02)
				signal_vars["sigp"].setMax(80.)
				signal_vars["xi"].setMin(-0.5)

		elif "trigbbh" in analysis or "NoTrigger_eta2p2" in analysis:
			#elif analysis == "trigbbh_CSVTM" or analysis == "trigbbh_notrig_CSVTM":
			sigp_initial = 55. + (160.-55.) / (1200 - 600) * (mass - 600)
			signal_vars["sigp"].setVal(sigp_initial)
			signal_vars["sigp"].setMin(sigp_initial - 30.)
			#signal_vars["sigp"].setMa(sigp_initial)
			rho2_initial = 0.06 + (0.22 - 0.06) / (1200 - 600) * (mass - 600)
			signal_vars["rho2"].setMin(0.04)
			signal_vars["rho2"].setVal(rho2_initial)
			xi_initial = -0.17 + (-0.53 + 0.17) / (1200 - 600) * (mass - 600)
			signal_vars["xi"].setMax(-0.05)
			signal_vars["xi"].setVal(xi_initial)
			if model == "RSG":
				# Take initial values from RSG 
				#if mass == 750:
				#	signal_vars["xi"].setMax(-0.1)
				#	signal_vars["xi"].setVal(-0.3)
				#	signal_vars["rho2"].setMin(0.07)
				#	signal_vars["rho2"].setVal(0.15)
				#	signal_vars["sigp"].setMin(50.)
				if mass == 1200:
					signal_vars["xp"].setVal(1.1395e+03)
					signal_vars["sigp"].setMin(75.)
					signal_vars["sigp"].setVal(1.6098e+02)
					signal_vars["xi"].setVal(-5.8503e-01)
					signal_vars["rho2"].setVal(2.5920e-01)
			elif model == "Hbb":
				# Take initial values from RSG 
				rho2_initial = 0.06 + (0.22 - 0.06) / (1200 - 600) * (mass - 600)
				signal_vars["rho2"].setMin(0.04)
				signal_vars["rho2"].setVal(rho2_initial)
				xi_initial = -0.17 + (-0.53 + 0.17) / (1200 - 600) * (mass - 600)
				signal_vars["xi"].setMax(-0.05)
				signal_vars["xi"].setVal(xi_initial)


		signal_vars["norm"] = RooRealVar(signal_pdf.GetName() + "_norm", signal_pdf.GetName() + "_norm", histogram.Integral(), histogram.Integral() / 50., histogram.Integral() * 50.)
		signal_epdf = ROOT.RooExtendPdf("signal_epdf", "signal_epdf", signal_pdf, signal_vars["norm"])

		if len(systematic_variations) > 0:
			signal_rdhs_syst = {}
			signal_pdfs_raw_syst = {}
			signal_pdfs_syst = {}
			signal_epdfs_syst = {}
			signal_vars_syst = {}
			for variation_name in systematic_variations:
				signal_rdhs_syst[variation_name] = RooDataHist(histograms_syst[variation_name].GetName() + "_rdh", histograms_syst[variation_name].GetName() + "_rdh", ROOT.RooArgList(mjj), histograms_syst[variation_name])
				signal_pdfs_raw_syst[variation_name], signal_vars_syst[variation_name] = make_signal_pdf(fit_function, mjj, tag=variation_name, mass=mass)
				signal_pdfs_raw_syst[variation_name].SetName("signal_pdf__" + variation_name + "_raw")
				if fit_trigger:
					signal_pdfs_syst[variation_name] = RooEffProd("signal_pdf__" + variation_name, "signal_pdf__" + variation_name, signal_pdfs_raw_syst[variation_name], trigger_efficiency_formula)
				else:
					signal_pdfs_syst[variation_name] = signal_pdfs_raw_syst[variation_name]
					signal_pdfs_syst[variation_name].SetName("signal_pdf__" + variation_name)
				signal_vars_syst[variation_name]["norm"] = RooRealVar(signal_pdfs_syst[variation_name].GetName() + "_norm", signal_pdfs_syst[variation_name].GetName() + "_norm", histograms_syst[variation_name].Integral(), histograms_syst[variation_name].Integral() / 50., histograms_syst[variation_name].Integral() * 50.)
				signal_epdfs_syst[variation_name] = ROOT.RooExtendPdf("signal_epdf__" + variation_name, "signal_epdf__" + variation_name, signal_pdfs_syst[variation_name], signal_vars_syst[variation_name]["norm"])

		print "[signal_fits] INFO : Printing starting parameters"
		for varname, var in signal_vars.iteritems():
			print "[signal_fits] INFO : \t{}".format(varname),
			var.Print()

		print "[signal_fits] INFO : Fitting " + fit_function + " nominal"
		fit_results = signal_epdf.fitTo(signal_rdh, RooFit.Save(kTRUE))
		print fit_results
		fit_results.Print()
		if len(systematic_variations) > 0:
			fit_results_syst = {}
			for variation_name in systematic_variations:
				print "[signal_fits] INFO : Fitting " + fit_function + " systematic " + variation_name
				setfix_systematic_parameters(fit_function, variation_name, signal_vars, signal_vars_syst[variation_name])
				fit_results_syst[variation_name] = signal_epdfs_syst[variation_name].fitTo(signal_rdhs_syst[variation_name], RooFit.Save(kTRUE))
				fit_results_syst[variation_name].Print()

		# Save
		w = RooWorkspace('w_signal','w_signal')
		getattr(w,'import')(signal_rdh, ROOT.RooCmdArg(), RooFit.Rename("signal_hist"))
		if fit_trigger:
			signal_pdf_raw.SetName("signal_bukin")
			signal_pdf.SetName("signal")
			getattr(w,'import')(signal_pdf, ROOT.RooCmdArg())
		else:
			signal_pdf.SetName("signal")
			getattr(w,'import')(signal_pdf, ROOT.RooCmdArg())
		signal_vars["norm"].SetName(signal_pdf.GetName() + "_norm")
		getattr(w,'import')(signal_vars["norm"], ROOT.RooCmdArg())
		fit_results.SetName("fit_results_signal")
		getattr(w,'import')(fit_results)
		#for name, var in signal_vars.iteritems():
		#	print "Import: " + var.GetName()
		#	getattr(w,'import')(var, ROOT.RooCmdArg(), RooFit.Rename(var.GetName()))
		for variation_name in systematic_variations:
			getattr(w,'import')(signal_rdhs_syst[variation_name], ROOT.RooCmdArg(), RooFit.Rename("signal_hist__" + variation_name))
			if fit_trigger:
				signal_pdfs_raw_syst[variation_name].SetName("signal_bukin__" + variation_name)
				signal_pdfs_syst[variation_name].SetName("signal__" + variation_name)
				getattr(w,'import')(signal_pdfs_syst[variation_name], ROOT.RooCmdArg(), RooFit.RecycleConflictNodes())
			else:
				signal_pdfs_syst[variation_name].SetName("signal__" + variation_name)
				getattr(w,'import')(signal_pdfs_syst[variation_name], ROOT.RooCmdArg(), RooFit.RecycleConflictNodes())
			signal_vars_syst[variation_name]["norm"].SetName(signal_pdfs_syst[variation_name].GetName() + "_norm")
			getattr(w,'import')(signal_vars_syst[variation_name]["norm"], ROOT.RooCmdArg())
			fit_results_syst[variation_name].SetName("fit_results_signal_" + variation_name)
			getattr(w,'import')(fit_results_syst[variation_name])
			#for name, var in signal_vars_syst[variation_name].iteritems():
			#	print "Import: " + var.GetName() + "__" + variation_name
			#	getattr(w,'import')(var, ROOT.RooCmdArg(), RooFit.Rename(var.GetName() + "__" + variation_name))
		w.Print()
		print "Saving fit results to " + analysis_config.get_signal_fit_file(analysis, model, mass, fit_function)
		f = TFile(analysis_config.get_signal_fit_file(analysis, model, mass, fit_function), "RECREATE")
		w.writeToFile(analysis_config.get_signal_fit_file(analysis, model, mass, fit_function))
		f.Close()
	print "[signal_fits] INFO : Done with this point."

def fit_parameter_table(analysis, model, masses):
	print "Saving table to " + analysis_config.figure_directory + "/table_bukin_parameters_{}_{}.tex".format(analysis, model)
	table = open(analysis_config.figure_directory + "/table_bukin_parameters_{}_{}.tex".format(analysis, model), 'w')
	variables = ["xp", "sigp", "rho2", "xi"]
	variables_pretty = {"xp":"$x_p$", "sigp":"$\\sigma_p$", "rho2":"$\\rho_2$", "xi":"$\\xi$"}
	table.write("\\begin{table}[h]\n")
	table.write("\t\\begin{tabular}{|c|c|c|c|c|}\n")
	table.write("\t\t\\hline\n")
	table.write("\t\tMass\t")
	for variable in variables:
		table.write("&\t{}\t".format(variables_pretty[variable]))
	table.write("\\\\\n\t\t\\hline\n")
	for mass in masses:
		table.write("\t\t{} GeV\t".format(mass))
		f = TFile(analysis_config.get_signal_fit_file(analysis, model, mass, "bukin"), "READ")
		w = f.Get("w_signal")
		#w.Print()
		for variable in variables:
			var = w.var(variable)
			table.write("&\t${}\pm{}$\t".format(round(var.getVal(), 2), round(var.getError(), 2)))
		table.write("\\\\\n\t\t\\hline\n")
	table.write("\t\\end{tabular}\n")
	table.write("\\end{table}\n")
	table.close()

def plot_fits(analysis, model, mass, fit_functions, systematics=[]):
	print "Welcome to plot_fits"
	print "analysis=",
	print analysis
	print "model=",
	print model
	print "mass=",
	print mass
	print "fit_functions=",
	print fit_functions
	print "systematics=",
	print systematics

	fit_trigger = not ("notrig" in analysis or "NoTrigger" in analysis)

	systematic_variations = []
	if "jer" in systematics:
		systematic_variations.append("JERUp")
		systematic_variations.append("JERDown")
	if "jes" in systematics:
		systematic_variations.append("JESUp")
		systematic_variations.append("JESDown")
	if "boff" in systematics:
		systematic_variations.append("BOffUp")
		systematic_variations.append("BOffDown")

	for fit_function in fit_functions:
		# Load objects
		print "[plot_fits] INFO : Opening " + analysis_config.get_signal_fit_file(analysis, model, mass, fit_function)
		f = TFile(analysis_config.get_signal_fit_file(analysis, model, mass, fit_function), "READ")
		w = f.Get("w_signal")
		mjj = w.var("mjj")
		fit_results = w.genobj("fit_results_signal")
		signal_rdh = w.data("signal_hist")
		signal_pdf = w.pdf("signal")
		signal_rdhs_syst = {}
		signal_pdfs_syst = {}
		fit_results_syst = {}
		for variation in systematic_variations:
			signal_rdhs_syst[variation] = w.data("signal_hist__" + variation)
			signal_pdfs_syst[variation] = w.pdf("signal__" + variation)
			fit_results_syst[variation] = w.genobj("fit_results_signal_" + variation)

		# Make extended PDFs
		signal_vars = {}
		signal_vars["norm"] = RooRealVar(signal_pdf.GetName() + "_norm", signal_pdf.GetName() + "_norm", signal_rdh.sum(False), signal_rdh.sum(False) / 50., signal_rdh.sum(False) * 50.)
		signal_epdf = ROOT.RooExtendPdf("signal_" + fit_function + "_epdf", "signal_" + fit_function + "_epdf", signal_pdf, signal_vars["norm"])
		signal_epdfs_syst = {}
		signal_vars_syst = {}
		for variation in systematic_variations:
			signal_vars_syst[variation] = {}
			signal_vars_syst[variation]["norm"] = RooRealVar(
				signal_pdfs_syst[variation].GetName() + "_norm", signal_pdfs_syst[variation].GetName() + "_norm", 
				signal_rdhs_syst[variation].sum(False), 
				signal_rdhs_syst[variation].sum(False) / 50., 
				signal_rdhs_syst[variation].sum(False) * 50.)
			signal_epdfs_syst[variation] = ROOT.RooExtendPdf("signal_" + fit_function + "_epdf__" + variation, "signal_" + fit_function + "_epdf__" + variation, signal_pdfs_syst[variation], signal_vars_syst[variation]["norm"])

		# Draw central value fit
		l = TLegend(0.6, 0.3, 0.88, 0.88)
		l.SetFillColor(0)
		l.SetBorderSize(0)
		frame = mjj.frame()
		frame_fine = mjj.frame()
		signal_rdh.plotOn(frame, RooFit.Name("Data"), RooFit.Binning(dijet_roobinning))
		signal_rdh.plotOn(frame_fine, RooFit.Name("Data"), RooFit.Binning(dijet_roobinning))
		#style_counter = 0
		signal_epdf.plotOn(frame, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(2), RooFit.Name(fit_function))
		signal_epdf.plotOn(frame_fine, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(2), RooFit.Name(fit_function))
		legend_entry = fit_function
		chi2ndf = frame_fine.chiSquare(fit_function, "Data", fit_results.floatParsFinal().getSize())
		legend_entry += " (#chi^{2}/NDF = " + str(round(chi2ndf, 2)) + ")"
		l.AddEntry(frame.findObject(fit_function), legend_entry, "l")
		#style_counter += 1
		c = TCanvas("c_signal_fits_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, "c_signal_fits_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, 800, 600)
		frame.Draw()
		l.Draw()

		# For Bukin, draw x1 and x2
		if fit_function == "bukin":
			consts = 2*(2*TMath.Log(2.))**0.5
			hp = w.var("sigp").getVal() * consts
			r4 = ((w.var("xi").getVal())**2+1)**0.5
			r1 = w.var("xi").getVal()/r4; 
			x1 = w.var("xp").getVal() + (hp / 2) * (r1-1)
			x2 = w.var("xp").getVal() + (hp / 2) * (r1+1)
			x1_line = TLine(x1, frame.GetMinimum(), x1, frame.GetMaximum())
			x1_line.SetLineStyle(1)
			x1_line.SetLineColor(kGray)
			x1_line.Draw()
			x2_line = TLine(x2, frame.GetMinimum(), x2, frame.GetMaximum())
			x2_line.SetLineStyle(1)
			x2_line.SetLineColor(kGray)
			x2_line.Draw()
		c.SaveAs(analysis_config.figure_directory + "/" + c.GetName() + ".pdf")

		# Systematic variations
		if "jer" in systematics:
			frame_jer = mjj.frame()
			c_jer = TCanvas("c_signal_hist_jer_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, "c_signal_hist_jer_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, 800, 600)
			signal_rdh.plotOn(frame_jer, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.Name(fit_function), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
			signal_rdhs_syst["JERUp"].plotOn(frame_jer, RooFit.LineColor(seaborn.GetColorRoot("default", 1)), RooFit.LineStyle(2), RooFit.Name(fit_function + ", JER+"), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
			signal_rdhs_syst["JERDown"].plotOn(frame_jer, RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(3), RooFit.Name(fit_function + ", JER-"), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
			frame_jer.Draw()
			l_jer = TLegend(0.6, 0.3, 0.88, 0.88)
			l_jer.SetFillColor(0)
			l_jer.SetBorderSize(0)
			l_jer.AddEntry(frame_jer.findObject(fit_function), fit_function, "l")
			l_jer.AddEntry(frame_jer.findObject(fit_function + ", JER+"), fit_function + ", JER+", "l")
			l_jer.AddEntry(frame_jer.findObject(fit_function + ", JER-"), fit_function + ", JER-", "l")
			l_jer.Draw()
			c_jer.SaveAs(analysis_config.figure_directory + "/" + c_jer.GetName() + ".pdf")

			frame_jer_fit = mjj.frame()
			c_jer_fit = TCanvas("c_signal_fits_jer_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, "c_signal_fits_jer_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, 800, 600)
			signal_epdf.plotOn(frame_jer_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.Name(fit_function))
			signal_epdfs_syst["JERUp"].plotOn(frame_jer_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 1)), RooFit.LineStyle(2), RooFit.Name(fit_function + ", JER+"))
			signal_epdfs_syst["JERDown"].plotOn(frame_jer_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(3), RooFit.Name(fit_function + ", JER-"))
			frame_jer_fit.Draw()
			l_jer_fit = TLegend(0.6, 0.3, 0.88, 0.88)
			l_jer_fit.SetFillColor(0)
			l_jer_fit.SetBorderSize(0)
			l_jer_fit.AddEntry(frame_jer_fit.findObject(fit_function), fit_function, "l")
			l_jer_fit.AddEntry(frame_jer_fit.findObject(fit_function + ", JER+"), fit_function + ", JER+", "l")
			l_jer_fit.AddEntry(frame_jer_fit.findObject(fit_function + ", JER-"), fit_function + ", JER-", "l")
			l_jer_fit.Draw()
			c_jer_fit.SaveAs(analysis_config.figure_directory + "/" + c_jer_fit.GetName() + ".pdf")

		if "jes" in systematics:
			frame_jes = mjj.frame()
			c_jes = TCanvas("c_signal_hist_jes_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, "c_signal_hist_jes_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, 800, 600)
			signal_rdh.plotOn(frame_jes, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.Name(fit_function), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
			signal_rdhs_syst["JESUp"].plotOn(frame_jes, RooFit.LineColor(seaborn.GetColorRoot("default", 1)), RooFit.LineStyle(2), RooFit.Name(fit_function + ", JES+"), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
			signal_rdhs_syst["JESDown"].plotOn(frame_jes, RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(3), RooFit.Name(fit_function + ", JES-"), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
			frame_jes.Draw()
			l_jes = TLegend(0.6, 0.3, 0.88, 0.88)
			l_jes.SetFillColor(0)
			l_jes.SetBorderSize(0)
			l_jes.AddEntry(frame_jes.findObject(fit_function), fit_function, "l")
			l_jes.AddEntry(frame_jes.findObject(fit_function + ", JES+"), fit_function + ", JES+", "l")
			l_jes.AddEntry(frame_jes.findObject(fit_function + ", JES-"), fit_function + ", JES-", "l")
			l_jes.Draw()
			c_jes.SaveAs(analysis_config.figure_directory + "/" + c_jes.GetName() + ".pdf")

			frame_jes_fit = mjj.frame()
			c_jes_fit = TCanvas("c_signal_fits_jes_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, "c_signal_fits_jes_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, 800, 600)
			signal_epdf.plotOn(frame_jes_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.Name(fit_function))
			signal_epdfs_syst["JESUp"].plotOn(frame_jes_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 1)), RooFit.LineStyle(2), RooFit.Name(fit_function + ", JES+"))
			signal_epdfs_syst["JESDown"].plotOn(frame_jes_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(3), RooFit.Name(fit_function + ", JES-"))
			frame_jes_fit.Draw()
			l_jes_fit = TLegend(0.6, 0.3, 0.88, 0.88)
			l_jes_fit.SetFillColor(0)
			l_jes_fit.SetBorderSize(0)
			l_jes_fit.AddEntry(frame_jes_fit.findObject(fit_function), fit_function, "l")
			l_jes_fit.AddEntry(frame_jes_fit.findObject(fit_function + ", JES+"), fit_function + ", JES+", "l")
			l_jes_fit.AddEntry(frame_jes_fit.findObject(fit_function + ", JES-"), fit_function + ", JES-", "l")
			l_jes_fit.Draw()
			c_jes_fit.SaveAs(analysis_config.figure_directory + "/" + c_jes_fit.GetName() + ".pdf")

		if "boff" in systematics:
			frame_boff = mjj.frame()
			c_boff = TCanvas("c_signal_hist_boff_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, "c_signal_hist_boff_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, 800, 600)
			signal_rdh.plotOn(frame_boff, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.Name(fit_function), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
			signal_rdhs_syst["BOffUp"].plotOn(frame_boff, RooFit.LineColor(seaborn.GetColorRoot("default", 1)), RooFit.LineStyle(2), RooFit.Name(fit_function + ", BOff+"), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
			signal_rdhs_syst["BOffDown"].plotOn(frame_boff, RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(3), RooFit.Name(fit_function + ", BOff-"), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
			frame_boff.Draw()
			l_boff = TLegend(0.6, 0.3, 0.88, 0.88)
			l_boff.SetFillColor(0)
			l_boff.SetBorderSize(0)
			l_boff.AddEntry(frame_boff.findObject(fit_function), fit_function, "l")
			l_boff.AddEntry(frame_boff.findObject(fit_function + ", BOff+"), fit_function + ", BOff+", "l")
			l_boff.AddEntry(frame_boff.findObject(fit_function + ", BOff-"), fit_function + ", BOff-", "l")
			l_boff.Draw()
			c_boff.SaveAs(analysis_config.figure_directory + "/" + c_boff.GetName() + ".pdf")

			frame_boff_fit = mjj.frame()
			c_boff_fit = TCanvas("c_signal_fits_boff_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, "c_signal_fits_boff_" + analysis + "_" + model + "_" + str(mass) + "_" + fit_function, 800, 600)
			signal_epdf.plotOn(frame_boff_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.Name(fit_function))
			signal_epdfs_syst["BOffUp"].plotOn(frame_boff_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 1)), RooFit.LineStyle(2), RooFit.Name(fit_function + ", BOff+"))
			signal_epdfs_syst["BOffDown"].plotOn(frame_boff_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(3), RooFit.Name(fit_function + ", BOff-"))
			frame_boff_fit.Draw()
			l_boff_fit = TLegend(0.6, 0.3, 0.88, 0.88)
			l_boff_fit.SetFillColor(0)
			l_boff_fit.SetBorderSize(0)
			l_boff_fit.AddEntry(frame_boff_fit.findObject(fit_function), fit_function, "l")
			l_boff_fit.AddEntry(frame_boff_fit.findObject(fit_function + ", BOff+"), fit_function + ", BOff+", "l")
			l_boff_fit.AddEntry(frame_boff_fit.findObject(fit_function + ", BOff-"), fit_function + ", BOff-", "l")
			l_boff_fit.Draw()
			c_boff_fit.SaveAs(analysis_config.figure_directory + "/" + c_boff_fit.GetName() + ".pdf")

		f.Close()

def plot_multiple_fits(analysis, model, masses, fit_functions, systematics=[]):
	print "Welcome to plot_multiple_fits"

	systematic_variations = []
	if "jer" in systematics:
		systematic_variations.append("JERUp")
		systematic_variations.append("JERDown")
	if "jes" in systematics:
		systematic_variations.append("JESUp")
		systematic_variations.append("JESDown")
	if "boff" in systematics:
		systematic_variations.append("BOffUp")
		systematic_variations.append("BOffDown")

	for fit_function in fit_functions:
		signal_rdhs = {}
		signal_pdfs = {}
		signal_epdfs = {}
		signal_rdhs_syst = {}
		signal_pdfs_syst = {}
		fit_results_syst = {}
		signal_vars = {}
		signal_epdfs_syst = {}
		signal_vars_syst = {}
		for mass in masses:
			# Load objects
			print "Loading objects from {}".format(analysis_config.get_signal_fit_file(analysis, model, mass, fit_function))
			f = TFile(analysis_config.get_signal_fit_file(analysis, model, mass, fit_function), "READ")
			w = f.Get("w_signal")
			mjj = w.var("mjj")
			fit_results = w.genobj("fit_results_signal")
			signal_rdhs[mass] = w.data("signal_hist")
			signal_rdhs[mass].SetName(signal_rdhs[mass].GetName() + str(mass))
			signal_pdfs[mass] = w.pdf("signal")
			signal_pdfs[mass].SetName(signal_pdfs[mass].GetName() + str(mass))
			signal_rdhs_syst[mass] = {}
			signal_pdfs_syst[mass] = {}
			fit_results_syst[mass] = {}
			for variation in systematic_variations:
				signal_rdhs_syst[mass][variation] = w.data("signal_hist__" + variation)
				signal_rdhs_syst[mass][variation].SetName(signal_rdhs_syst[mass][variation].GetName() + str(mass))
				signal_pdfs_syst[mass][variation] = w.pdf("signal__" + variation)
				signal_pdfs_syst[mass][variation].SetName(signal_pdfs_syst[mass][variation].GetName() + str(mass))
				fit_results_syst[mass][variation] = w.genobj("fit_results_signal_" + variation)
				fit_results_syst[mass][variation].SetName(fit_results_syst[mass][variation].GetName() + str(mass))

			# Make extended PDFs
			signal_vars[mass] = {}
			signal_vars[mass]["norm"] = RooRealVar(signal_pdfs[mass].GetName() + "_norm" + str(mass), signal_pdfs[mass].GetName() + "_norm", signal_rdhs[mass].sum(False), signal_rdhs[mass].sum(False) / 50., signal_rdhs[mass].sum(False) * 50.)
			signal_epdfs[mass] = ROOT.RooExtendPdf("signal_" + fit_function + "_epdf" + str(mass), "signal_" + fit_function + "_epdf", signal_pdfs[mass], signal_vars[mass]["norm"])
			signal_epdfs_syst[mass] = {}
			signal_vars_syst[mass] = {}
			for variation in systematic_variations:
				signal_vars_syst[mass][variation] = {}
				signal_vars_syst[mass][variation]["norm"] = RooRealVar(signal_pdfs_syst[mass][variation].GetName() + "_norm" + str(mass), signal_pdfs_syst[mass][variation].GetName() + "_norm", signal_rdhs_syst[mass][variation].sum(False), signal_rdhs_syst[mass][variation].sum(False) / 50., signal_rdhs_syst[mass][variation].sum(False) * 50.)
				signal_epdfs_syst[mass][variation] = ROOT.RooExtendPdf("signal_" + fit_function + "_epdf__" + variation + str(mass), "signal_" + fit_function + "_epdf__" + variation, signal_pdfs_syst[mass][variation], signal_vars_syst[mass][variation]["norm"])


		# Systematic variations
		if "jer" in systematics:
			frame_jer = mjj.frame()
			c_jer = TCanvas("c_signal_hist_jer_" + analysis + "_" + model + "_multiplemasses_" + fit_function, "c_signal_hist_jer_" + analysis + "_" + model + "_multiplemasses_" + fit_function, 800, 600)
			l_jer = TLegend(0.55, 0.6, 0.88, 0.88)
			l_jer.SetFillColor(0)
			l_jer.SetBorderSize(0)
			for mass in masses:
				signal_rdhs[mass].plotOn(frame_jer, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.Name("Central, {} GeV".format(mass)), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
				signal_rdhs_syst[mass]["JERUp"].plotOn(frame_jer, RooFit.LineColor(seaborn.GetColorRoot("default", 1)), RooFit.LineStyle(2), RooFit.Name("JER+, {} GeV".format(mass)), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
				signal_rdhs_syst[mass]["JERDown"].plotOn(frame_jer, RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(3), RooFit.Name("JER-, {} GeV".format(mass)), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
			l_jer.AddEntry(frame_jer.findObject("Central, {} GeV".format(masses[0])), "Central", "l")
			l_jer.AddEntry(frame_jer.findObject("JER+, {} GeV".format(masses[0])), "JER+", "l")
			l_jer.AddEntry(frame_jer.findObject("JER-, {} GeV".format(masses[0])), "JER-", "l")
			frame_jer.Draw()
			l_jer.Draw()
			c_jer.SaveAs(analysis_config.figure_directory + "/" + c_jer.GetName() + ".pdf")

			frame_jer_fit = mjj.frame()
			c_jer_fit = TCanvas("c_signal_fits_jer_" + analysis + "_" + model + "_multiplemasses_" + fit_function, "c_signal_fits_jer_" + analysis + "_" + model + "_multiplemasses_" + fit_function, 800, 600)
			l_jer_fit = TLegend(0.55, 0.6, 0.88, 0.88)
			l_jer_fit.SetFillColor(0)
			l_jer_fit.SetBorderSize(0)
			for mass in masses:
				signal_epdfs[mass].plotOn(frame_jer_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.Name("Central, {} GeV".format(mass)))
				signal_epdfs_syst[mass]["JERUp"].plotOn(frame_jer_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 1)), RooFit.LineStyle(2), RooFit.Name("JER+, {} GeV".format(mass)))
				signal_epdfs_syst[mass]["JERDown"].plotOn(frame_jer_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(3), RooFit.Name("JER-, {} GeV".format(mass)))
			l_jer_fit.AddEntry(frame_jer_fit.findObject("Central, {} GeV".format(masses[0])), "Central", "l")
			l_jer_fit.AddEntry(frame_jer_fit.findObject("JER+, {} GeV".format(masses[0])), "JER+", "l")
			l_jer_fit.AddEntry(frame_jer_fit.findObject("JER-, {} GeV".format(masses[0])), "JER-", "l")
			frame_jer_fit.Draw()
			l_jer_fit.Draw()
			c_jer_fit.SaveAs(analysis_config.figure_directory + "/" + c_jer_fit.GetName() + ".pdf")

		if "jes" in systematics:
			frame_jes = mjj.frame()
			c_jes = TCanvas("c_signal_hist_jes_" + analysis + "_" + model + "_multiplemasses_" + fit_function, "c_signal_hist_jes_" + analysis + "_" + model + "_multiplemasses_" + fit_function, 800, 600)
			l_jes = TLegend(0.55, 0.6, 0.88, 0.88)
			l_jes.SetFillColor(0)
			l_jes.SetBorderSize(0)
			for mass in masses:
				signal_rdhs[mass].plotOn(frame_jes, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.Name("Central, {} GeV".format(mass)), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
				signal_rdhs_syst[mass]["JESUp"].plotOn(frame_jes, RooFit.LineColor(seaborn.GetColorRoot("default", 1)), RooFit.LineStyle(2), RooFit.Name("JES+, {} GeV".format(mass)), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
				signal_rdhs_syst[mass]["JESDown"].plotOn(frame_jes, RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(3), RooFit.Name("JES-, {} GeV".format(mass)), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
			l_jes.AddEntry(frame_jes.findObject("Central, {} GeV".format(masses[0])), "Central", "l")
			l_jes.AddEntry(frame_jes.findObject("JES+, {} GeV".format(masses[0])), "JES+", "l")
			l_jes.AddEntry(frame_jes.findObject("JES-, {} GeV".format(masses[0])), "JES-", "l")
			frame_jes.Draw()
			l_jes.Draw()
			c_jes.SaveAs(analysis_config.figure_directory + "/" + c_jes.GetName() + ".pdf")

			frame_jes_fit = mjj.frame()
			c_jes_fit = TCanvas("c_signal_fits_jes_" + analysis + "_" + model + "_multiplemasses_" + fit_function, "c_signal_fits_jes_" + analysis + "_" + model + "_multiplemasses_" + fit_function, 800, 600)
			l_jes_fit = TLegend(0.55, 0.6, 0.88, 0.88)
			l_jes_fit.SetFillColor(0)
			l_jes_fit.SetBorderSize(0)
			for mass in masses:
				signal_epdfs[mass].plotOn(frame_jes_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.Name("Central, {} GeV".format(mass)))
				signal_epdfs_syst[mass]["JESUp"].plotOn(frame_jes_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 1)), RooFit.LineStyle(2), RooFit.Name("JES+, {} GeV".format(mass)))
				signal_epdfs_syst[mass]["JESDown"].plotOn(frame_jes_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(3), RooFit.Name("JES-, {} GeV".format(mass)))
			l_jes_fit.AddEntry(frame_jes_fit.findObject("Central, {} GeV".format(masses[0])), "Central", "l")
			l_jes_fit.AddEntry(frame_jes_fit.findObject("JES+, {} GeV".format(masses[0])), "JES+", "l")
			l_jes_fit.AddEntry(frame_jes_fit.findObject("JES-, {} GeV".format(masses[0])), "JES-", "l")
			frame_jes_fit.Draw()
			l_jes_fit.Draw()
			c_jes_fit.SaveAs(analysis_config.figure_directory + "/" + c_jes_fit.GetName() + ".pdf")

		if "boff" in systematics:
			frame_boff = mjj.frame()
			c_boff = TCanvas("c_signal_hist_boff_" + analysis + "_" + model + "_multiplemasses_" + fit_function, "c_signal_hist_boff_" + analysis + "_" + model + "_multiplemasses_" + fit_function, 800, 600)
			l_boff = TLegend(0.55, 0.6, 0.88, 0.88)
			l_boff.SetFillColor(0)
			l_boff.SetBorderSize(0)
			for mass in masses:
				signal_rdhs[mass].plotOn(frame_boff, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.Name("Central, {} GeV".format(mass)), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
				signal_rdhs_syst[mass]["BOffUp"].plotOn(frame_boff, RooFit.LineColor(seaborn.GetColorRoot("default", 1)), RooFit.LineStyle(2), RooFit.Name("BOff+, {} GeV".format(mass)), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
				signal_rdhs_syst[mass]["BOffDown"].plotOn(frame_boff, RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(3), RooFit.Name("BOff-, {} GeV".format(mass)), RooFit.Binning(dijet_roobinning), RooFit.MarkerStyle(20), RooFit.MarkerSize(0))
			l_boff.AddEntry(frame_boff.findObject("Central, {} GeV".format(masses[0])), "Central", "l")
			l_boff.AddEntry(frame_boff.findObject("BOff+, {} GeV".format(masses[0])), "BOff+", "l")
			l_boff.AddEntry(frame_boff.findObject("BOff-, {} GeV".format(masses[0])), "BOff-", "l")
			frame_boff.Draw()
			l_boff.Draw()
			c_boff.SaveAs(analysis_config.figure_directory + "/" + c_boff.GetName() + ".pdf")

			frame_boff_fit = mjj.frame()
			c_boff_fit = TCanvas("c_signal_fits_boff_" + analysis + "_" + model + "_multiplemasses_" + fit_function, "c_signal_fits_boff_" + analysis + "_" + model + "_multiplemasses_" + fit_function, 800, 600)
			l_boff_fit = TLegend(0.55, 0.6, 0.88, 0.88)
			l_boff_fit.SetFillColor(0)
			l_boff_fit.SetBorderSize(0)
			for mass in masses:
				signal_epdfs[mass].plotOn(frame_boff_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.Name("Central, {} GeV".format(mass)))
				signal_epdfs_syst[mass]["BOffUp"].plotOn(frame_boff_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 1)), RooFit.LineStyle(2), RooFit.Name("BOff+, {} GeV".format(mass)))
				signal_epdfs_syst[mass]["BOffDown"].plotOn(frame_boff_fit, RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(3), RooFit.Name("BOff-, {} GeV".format(mass)))
			l_boff_fit.AddEntry(frame_boff_fit.findObject("Central, {} GeV".format(masses[0])), "Central", "l")
			l_boff_fit.AddEntry(frame_boff_fit.findObject("BOff+, {} GeV".format(masses[0])), "BOff+", "l")
			l_boff_fit.AddEntry(frame_boff_fit.findObject("BOff-, {} GeV".format(masses[0])), "BOff-", "l")
			frame_boff_fit.Draw()
			l_boff_fit.Draw()
			c_boff_fit.SaveAs(analysis_config.figure_directory + "/" + c_boff_fit.GetName() + ".pdf")
		f.Close()


# Get values of fit parameters as a dict, excluding the independent variable.
def get_signal_parameters(pdf, exclude=["mjj"], tag=None):
	parameters = pdf.getVariables()
	parameter_values = {}
	parameter_uncertainties = {}
	iterator = parameters.createIterator()
	parameter = iterator.Next()
	while parameter:
		name = parameter.GetName()
		if tag:
			name = name.replace(tag, "")
		if not name in exclude:
			parameter_values[name] = parameter.getVal()
			parameter_uncertainties[name] = parameter.getError()
		parameter = iterator.Next()
	return parameter_values, parameter_uncertainties



def signal_interpolations(analysis, model, input_masses, output_masses, fit_function, systematic_variations=[], spline_interpolation=False):
	print "Welcome to signal_interpolations"
	print "analysis = {}".format(analysis)
	print "model = {}".format(model)
	print "input_masses = ",
	print input_masses
	print "output_masses = ",
	print output_masses
	print "fit_function = {}".format(fit_function)
	print "systematic_variations = ",
	print systematic_variations

	input_masses.sort()
	# Load input shapes
	input_shapes = {}
	input_parameters = {}
	input_parameter_uncertainties = {}
	input_parameter_graphs = {}
	input_parameter_splines = {}
	parameter_names = {}

	fit_trigger = not ("notrig" in analysis or "NoTrigger" in analysis)

	for input_mass in input_masses:
		index = input_masses.index(input_mass)
		f = TFile(analysis_config.get_signal_fit_file(analysis, model, input_mass, fit_function), "READ")
		w = f.Get("w_signal")
		if index == 0:
			mjj = w.var("mjj")

		input_shapes[input_mass] = {}
		input_parameters[input_mass] = {}
		input_parameter_uncertainties[input_mass] = {}
		# Nominal
		if fit_trigger:
			input_shapes[input_mass]["nominal"] = w.pdf("signal_bukin")
		else:
			input_shapes[input_mass]["nominal"] = w.pdf("signal")
		input_parameters[input_mass]["nominal"], input_parameter_uncertainties[input_mass]["nominal"] = get_signal_parameters(input_shapes[input_mass]["nominal"])

		for systematic_variation in systematic_variations:
			if fit_trigger:
				input_shapes[input_mass][systematic_variation] = w.pdf("signal_bukin__" + systematic_variation)
			else:
				input_shapes[input_mass][systematic_variation] = w.pdf("signal__" + systematic_variation)
			input_parameters[input_mass][systematic_variation], input_parameter_uncertainties[input_mass][systematic_variation] = get_signal_parameters(input_shapes[input_mass][systematic_variation], tag=systematic_variation)

		# Make TGraphs
		if index == 0:
			input_parameter_graphs["nominal"] = {}
			for systematic_variation in systematic_variations:
				input_parameter_graphs[systematic_variation] = {}

			for variation, parameter_name_values in input_parameters[input_mass].iteritems():
				parameter_names[variation] = []
				for parameter_name, parameter_value in parameter_name_values.iteritems():
					parameter_names[variation].append(parameter_name)
					input_parameter_graphs[variation][parameter_name] = TGraphErrors(len(input_masses))

		for variation, parameter_name_value in input_parameters[input_mass].iteritems():
			for parameter_name, parameter_value in parameter_name_value.iteritems():
				input_parameter_graphs[variation][parameter_name].SetPoint(index, input_mass, parameter_value)
				input_parameter_graphs[variation][parameter_name].SetPointError(index, 0., input_parameter_uncertainties[input_mass][variation][parameter_name])
		f.Close()

	# Convert graphs to splines
	if spline_interpolation:
		for systematic_variation, parameter_graphs in input_parameter_graphs.iteritems():
			input_parameter_splines[systematic_variation] = {}
			for parameter, graph in parameter_graphs.iteritems():
				input_parameter_splines[systematic_variation][parameter] = TSpline3("spline_" + systematic_variation + "_" + parameter, graph)

	if fit_trigger:
		trigger_efficiency_formula = trigger_efficiency.get_formula(analysis, mjj)


	output_shapes = {}
	output_shapes_raw = {}
	output_parameters = {}
	for output_mass in output_masses:
		# Warn if output_mass is outside of input_mass range
		if output_mass < min(input_masses):
			print "[signal_interpolations] WARNING : output mass {} is below lowest input mass {}. Value will be extrapolated, which may be unreliable.".format(output_mass, min(input_masses))
		if output_mass > max(input_masses):
			print "[signal_interpolations] WARNING : output mass {} is above highest input mass {}. Value will be extrapolated, which may be unreliable.".format(output_mass, max(input_masses))

		output_shapes_raw[output_mass] = {}
		output_shapes[output_mass] = {}
		output_parameters[output_mass] = {}
		output_shapes_raw[output_mass]["nominal"], output_parameters[output_mass]["nominal"] = copy_signal_pdf(fit_function, input_shapes[input_masses[0]]["nominal"], mjj, mass=output_mass)
		if "xp" in output_parameters[output_mass]["nominal"]:
			output_parameters[output_mass]["nominal"]["xp"].setMin(output_mass - 200.)
			output_parameters[output_mass]["nominal"]["xp"].setMax(output_mass + 200.)
		for parameter_name in parameter_names["nominal"]:
			if spline_interpolation:	
				output_parameters[output_mass]["nominal"][parameter_name].setVal(input_parameter_splines["nominal"][parameter_name].Eval(output_mass))
			else:
				output_parameters[output_mass]["nominal"][parameter_name].setVal(input_parameter_graphs["nominal"][parameter_name].Eval(output_mass))
			print "[debug] Mass {}, set par {} to {}".format(output_mass, parameter_name, output_parameters[output_mass]["nominal"][parameter_name].getVal())
		output_shapes_raw[output_mass]["nominal"].SetName("signal_bukin")
		if fit_trigger:
			output_shapes[output_mass]["nominal"] = RooEffProd("signal", "signal", output_shapes_raw[output_mass]["nominal"], trigger_efficiency_formula)
		else:
			output_shapes[output_mass]["nominal"] = output_shapes_raw[output_mass]["nominal"]
			output_shapes[output_mass]["nominal"].SetName("signal")

		for systematic_variation in systematic_variations:
			output_shapes_raw[output_mass][systematic_variation], output_parameters[output_mass][systematic_variation] = copy_signal_pdf(fit_function, input_shapes[input_masses[0]][systematic_variation], mjj, tag=systematic_variation, mass=output_mass)
			for parameter_name in parameter_names[systematic_variation]:
				if spline_interpolation:
					output_parameters[output_mass][systematic_variation][parameter_name].setVal(input_parameter_splines[systematic_variation][parameter_name].Eval(output_mass))
				else:
					output_parameters[output_mass][systematic_variation][parameter_name].setVal(input_parameter_graphs[systematic_variation][parameter_name].Eval(output_mass))
			output_shapes_raw[output_mass][systematic_variation].SetName("signal_bukin__" + systematic_variation)
			if fit_trigger:
				output_shapes[output_mass][systematic_variation] = RooEffProd("signal__" + variation, "signal__" + variation, output_shapes_raw[output_mass][systematic_variation], trigger_efficiency_formula)
			else:
				output_shapes[output_mass][systematic_variation] = output_shapes_raw[output_mass][systematic_variation]
				output_shapes[output_mass][systematic_variation].SetName("signal__" + variation)

		# Save
		w = RooWorkspace('w_signal','w_signal')
		output_shapes[output_mass]["nominal"].SetName("signal")
		getattr(w,'import')(output_shapes[output_mass]["nominal"], ROOT.RooCmdArg(), RooFit.RecycleConflictNodes())
		for systematic_variation in systematic_variations:
			output_shapes[output_mass][systematic_variation].SetName("signal__" + systematic_variation)
			getattr(w,'import')(output_shapes[output_mass][systematic_variation], ROOT.RooCmdArg(), RooFit.RecycleConflictNodes())
			#for name, var in signal_vars_syst[variation_name].iteritems():
			#	print "Import: " + var.GetName() + "__" + variation_name
			#	getattr(w,'import')(var, ROOT.RooCmdArg(), RooFit.Rename(var.GetName() + "__" + variation_name))
		w.Print()
		f = TFile(analysis_config.get_signal_fit_file(analysis, model, output_mass, fit_function, interpolated=True), "RECREATE")
		w.writeToFile(analysis_config.get_signal_fit_file(analysis, model, output_mass, fit_function, interpolated=True))
		f.Close()

	# Draw parameter interpolations
	for systematic_variation, parameter_graphs in input_parameter_graphs.iteritems():
		for parameter in parameter_graphs:
			c = TCanvas("c_signal_interpolation_" + systematic_variation + "_" + parameter + "_" + analysis + "_" + model + "_" + fit_function, "c_signal_interpolation_" + systematic_variation + "_" + parameter + "_" + analysis + "_" + model + "_" + fit_function, 800, 600)
			input_parameter_graphs[systematic_variation][parameter].SetMarkerStyle(21)
			input_parameter_graphs[systematic_variation][parameter].SetMarkerColor(1)
			input_parameter_graphs[systematic_variation][parameter].SetMarkerSize(1)
			input_parameter_graphs[systematic_variation][parameter].GetXaxis().SetTitle(parameter)
			input_parameter_graphs[systematic_variation][parameter].Draw("ap")
			if spline_interpolation:
				input_parameter_splines[systematic_variation][parameter].SetLineStyle(1)
				input_parameter_splines[systematic_variation][parameter].SetLineWidth(1)
				input_parameter_splines[systematic_variation][parameter].SetLineColor(2)
				input_parameter_splines[systematic_variation][parameter].Draw("same")
			c.SaveAs(analysis_config.figure_directory + "/"+ c.GetName() + ".pdf")

	# Draw all signal fits
	frame_shapes = mjj.frame()
	for input_mass in input_masses:
		input_shapes[input_mass]["nominal"].plotOn(frame_shapes, RooFit.Name("{} GeV".format(input_mass)), RooFit.LineStyle(1), RooFit.LineWidth(2), RooFit.LineColor(seaborn.GetColorRoot("cubehelix", int(input_mass / 1200. * 20), 20)))
	for output_mass in output_masses:
		output_shapes[output_mass]["nominal"].plotOn(frame_shapes, RooFit.Name("{} GeV".format(output_mass)), RooFit.LineStyle(2), RooFit.LineWidth(2), RooFit.LineColor(seaborn.GetColorRoot("cubehelix", int(output_mass / 1200. * 20), 20)))
	c_frame_shapes = TCanvas("c_signal_fits_all_{}_{}_{}".format(analysis, model, fit_function), "c_signal_fits_{}_{}_{}".format(analysis, model, fit_function), 800, 600)
	frame_shapes.Draw()
	c_frame_shapes.SaveAs(analysis_config.figure_directory + "/{}.pdf".format(c_frame_shapes.GetName()))

if __name__ == "__main__":
	from argparse import ArgumentParser
	parser = ArgumentParser(description='Run signal roofits and save to workspace')
	parser.add_argument("--fit", action="store_true", help="Run signal fits")
	parser.add_argument("--plots", action="store_true", help="Plot signal fits")
	parser.add_argument("--table", action="store_true", help="Table of fit parameters")
	parser.add_argument("--interpolate", action="store_true", help="Run signal interpolation")
	parser.add_argument("--interpolate375", action="store_true", help="Run signal interpolation at 375 GeV (for granularity study)")
	parser.add_argument("--validate_interpolation", action="store_true", help="Compare actual fits to interpolations")
	#parser.add_argument("--analyses", type=str, default="NoTrigger_eta1p7_CSVTM,NoTrigger_eta2p2_CSVTM", help="List of analyses to run (comma-separated)")
	parser.add_argument("--analyses", type=str, default="trigbbl_CSVTM,trigbbh_CSVTM", help="List of analyses to run (comma-separated)")
	parser.add_argument("--models", type=str, default="Hbb,ZPrime,RSG", help="List of models to run")
	parser.add_argument("--mass", type=int, help="Manually specify mass (otherwise, script runs all mass points)")
	parser.add_argument("--fit_functions", type=str, default="bukin", help="List of fit functions")
	parser.add_argument("--no_systematics", action="store_true", help="Run without systematics")
	args = parser.parse_args()

	analyses = args.analyses.split(",")
	models = args.models.split(",")
	masses = {}
	if args.mass:
		for analysis in analyses:
			masses[analysis] = [args.mass]
	else:
		# No need for limited masses, Tyler's files have all the interpolations
		for analysis in analyses:
			if "trigbbl" in analysis or "eta1p7" in analysis:
				masses[analysis] = [325,350,400,500,600,750,900]
			elif "trigbbh" in analysis or "eta2p2" in analysis:
				masses[analysis] = [600, 750, 900, 1200]
	fit_functions = args.fit_functions.split(",")
	if args.no_systematics:
		print "[signal_fits] INFO : Running without systematics"
		systematics_arg = []
	else:
		systematics_arg = systematics

	if args.fit:
		from joblib import Parallel, delayed
		for analysis in analyses:
			for model in models:
				Parallel(n_jobs=4)(delayed(signal_fit)(analysis, model, mass, fit_functions, systematics=systematics_arg) for mass in masses[analysis])
				#for mass in masses[analysis]:
				#	#signal_fit(analysis, model, mass, fit_functions, systematics=systematics_arg, fitTrigger=True)
				#	#signal_fit(analysis, model, mass, fit_functions, systematics=systematics_arg, correctTrigger=True)
				#	signal_fit(analysis, model, mass, fit_functions, systematics=systematics_arg)

	if args.plots:
		for analysis in analyses:
			for model in models:
				for mass in masses[analysis]:
					plot_fits(analysis, model, mass, fit_functions, systematics=systematics_arg)
				plot_multiple_fits(analysis, model, masses[analysis], fit_functions, systematics=systematics_arg)
	if args.table:
		for analysis in analyses:
			for model in models:
				fit_parameter_table(analysis, model, masses[analysis])


	if args.interpolate:
		output_masses = {}
		if args.mass:
			print "[signal_fits] ERROR : Argument mass not supported for interpolation. The input/output masses are too complicated for the command line, and are written into the code."
			sys.exit(1)
		for analysis in analyses:
			if "trigbbl" in analysis or "eta1p7" in analysis:
				output_masses[analysis] = list(set(range(350, 950, 50)) - set(masses[analysis]))
			elif "trigbbh" in analysis or "eta2p2" in analysis:
				output_masses[analysis] = list(set(range(600, 1250, 50)) - set(masses[analysis]))
		#output_masses = {
		#	"trigbbll_CSVTM":list(set(range(400, 950, 50)) - set(masses["trigbbll_CSVTM"])),
		#	"trigbbl_CSVTM":list(set(range(400, 950, 50)) - set(masses["trigbbl_CSVTM"])),
		#	"trigbbhl_CSVTM":list(set(range(600, 1250, 50)) - set(masses["trigbbhl_CSVTM"])),
		#	"trigbbh_CSVTM":list(set(range(600, 1250, 50)) - set(masses["trigbbh_CSVTM"])),
		#	"trigbbl_CSVM":list(set(range(400, 950, 50)) - set(masses["trigbbl_CSVM"])),
		#	"trigbbh_CSVM":list(set(range(600, 1250, 50)) - set(masses["trigbbh_CSVM"])),
		#	"trigbbl_notrig_CSVTM":list(set(range(400, 950, 50)) - set(masses["trigbbl_notrig_CSVTM"])),
		#	"trigbbh_notrig_CSVTM":list(set(range(600, 1250, 50)) - set(masses["trigbbh_notrig_CSVTM"])),
		#}
		for analysis in analyses:
			for model in models:
				for fit_function in fit_functions:
					signal_interpolations(analysis, model, masses[analysis], output_masses[analysis], fit_function, systematic_variations=["JERUp", "JERDown", "JESUp", "JESDown", "BOffUp", "BOffDown"])

	if args.interpolate375:
		output_masses = {}
		if args.mass:
			print "[signal_fits] ERROR : Argument mass not supported for interpolation. The input/output masses are too complicated for the command line, and are written into the code."
			sys.exit(1)
		for analysis in analyses:
			if "trigbbl" in analysis or "eta1p7" in analysis:
				output_masses[analysis] = [375]
			elif "trigbbh" in analysis or "eta2p2" in analysis:
				output_masses[analysis] = [375]
		#output_masses = {
		#	"trigbbll_CSVTM":list(set(range(400, 950, 50)) - set(masses["trigbbll_CSVTM"])),
		#	"trigbbl_CSVTM":list(set(range(400, 950, 50)) - set(masses["trigbbl_CSVTM"])),
		#	"trigbbhl_CSVTM":list(set(range(600, 1250, 50)) - set(masses["trigbbhl_CSVTM"])),
		#	"trigbbh_CSVTM":list(set(range(600, 1250, 50)) - set(masses["trigbbh_CSVTM"])),
		#	"trigbbl_CSVM":list(set(range(400, 950, 50)) - set(masses["trigbbl_CSVM"])),
		#	"trigbbh_CSVM":list(set(range(600, 1250, 50)) - set(masses["trigbbh_CSVM"])),
		#	"trigbbl_notrig_CSVTM":list(set(range(400, 950, 50)) - set(masses["trigbbl_notrig_CSVTM"])),
		#	"trigbbh_notrig_CSVTM":list(set(range(600, 1250, 50)) - set(masses["trigbbh_notrig_CSVTM"])),
		#}
		for analysis in analyses:
			for model in models:
				for fit_function in fit_functions:
					signal_interpolations(analysis, model, masses[analysis], output_masses[analysis], fit_function, systematic_variations=["JERUp", "JERDown", "JESUp", "JESDown", "BOffUp", "BOffDown"])					