import os
import sys
import ROOT
from ROOT import *

def make_background_pdf(function_name, mjj, collision_energy=8000.):
	fit_parameters = {}
	if function_name == "f1":
		fit_parameters["p1"] = RooRealVar('f1_p1','f1_p1',5.e-3,0.,100.)
		fit_parameters["p2"] = RooRealVar('f1_p2','f1_p2',9.1,0.,60.)
		fit_parameters["p3"] = RooRealVar('f1_p3','f1_p3',0.5,-10.,10.)
		background_pdf = RooGenericPdf('background_' + function_name,'(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)))'%(collision_energy,collision_energy,collision_energy),RooArgList(mjj,fit_parameters["p1"],fit_parameters["p2"],fit_parameters["p3"]))
	elif function_name == "f2":
		fit_parameters["p1"] = RooRealVar('f2_p1','f2_p1',5.6,0.,100.)
		fit_parameters["p2"] = RooRealVar('f2_p2','f2_p2',10.,0.,60.)
		background_pdf = RooGenericPdf('background_' + function_name,'(pow(@0/%.1f,-@1)*pow(1-@0/%.1f,@2))'%(collision_energy, collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"]))
	elif function_name == "f3":
		fit_parameters["p1"] = RooRealVar('f3_p1','f3_p1',80, 0., 100.)
		fit_parameters["p2"] = RooRealVar('f3_p2','f3_p2',8., 0., 100.)
		background_pdf = RooGenericPdf('background_' + function_name,'(1/pow(1+(@1*@0/%.1f),@2))'%(collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"]))
	elif function_name == "f4":
		fit_parameters["p1"] = RooRealVar('f4_p1', 'f4_p1', 80, 0., 100.)
		fit_parameters["p2"] = RooRealVar('f4_p2', 'f4_p2', 0., -100., 100.)
		fit_parameters["p3"] = RooRealVar('f4_p3', 'f4_p3', 8.0, 0., 100.)
		background_pdf = RooGenericPdf('background_' + function_name,'(1/pow(1+(@1*@0/%.1f)+(@2*pow(@0/%.1f,2)),@3))'%(collision_energy,collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"], fit_parameters["p3"]))
	elif function_name == "f5":
		fit_parameters["p1"] = RooRealVar('f5_p1','f5_p1',4.8,-100., 100.)
		fit_parameters["p2"] = RooRealVar('f5_p2','f5_p2',7., -100., 100.)
		background_pdf = RooGenericPdf('background_' + function_name,'(pow(@0/%.1f,-@1)*pow(1-pow(@0/%.1f,1/3),@2))'%(collision_energy, collision_energy),RooArgList(mjj,fit_parameters["p1"],fit_parameters["p2"]))
	else:
		print "[roofit_functions.make_background_pdf] ERROR : Unrecognized fit function " + function_name
		sys.exit(1)
	return (background_pdf, fit_parameters)
