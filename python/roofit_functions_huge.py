# This is similar to module roofit_functions, but has a bazillion more options for parameters.
# Intended for F-tests, not for use in the analysis itself!

import os
import sys
import ROOT
from ROOT import *

def make_background_pdf(function_name, mjj, collision_energy=8000.):
	fit_parameters = {}
	if function_name == "dijet3":
		fit_parameters["p1"] = RooRealVar('dijet4_p1','dijet4_p1',-1.,-50.,50.)
		fit_parameters["p2"] = RooRealVar('dijet4_p2','dijet4_p2',9.1,0.,60.)
		background_pdf = RooGenericPdf('background_' + function_name,'(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2))'%(collision_energy,collision_energy),RooArgList(mjj,fit_parameters["p1"],fit_parameters["p2"]))
	elif function_name == "dijet4":
		fit_parameters["p1"] = RooRealVar('dijet4_p1','dijet4_p1',-1.,-50.,50.)
		fit_parameters["p2"] = RooRealVar('dijet4_p2','dijet4_p2',9.1,0.,60.)
		fit_parameters["p3"] = RooRealVar('dijet4_p3','dijet4_p3',0.5,-10.,10.)
		background_pdf = RooGenericPdf('background_' + function_name,'(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)))'%(collision_energy,collision_energy,collision_energy),RooArgList(mjj,fit_parameters["p1"],fit_parameters["p2"],fit_parameters["p3"]))
	elif function_name == "dijet5":
		fit_parameters["p1"] = RooRealVar('dijet5_p1','dijet5_p1',-1.,-50.,50.)
		fit_parameters["p2"] = RooRealVar('dijet5_p2','dijet5_p2',9.1,0.,60.)
		fit_parameters["p3"] = RooRealVar('dijet5_p3','dijet5_p3',0.5,-10.,10.)
		fit_parameters["p4"] = RooRealVar('dijet5_p4','dijet5_p4',0.5,-10.,10.)
		background_pdf = RooGenericPdf('background_' + function_name,'(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)+@4*pow(log(@0/%.1f), 2)))'%(collision_energy,collision_energy,collision_energy,collision_energy),RooArgList(mjj,fit_parameters["p1"],fit_parameters["p2"],fit_parameters["p3"],fit_parameters["p4"]))
	elif function_name == "dijet6":
		fit_parameters["p1"] = RooRealVar('dijet6_p1','dijet6_p1',-1.,-50.,50.)
		fit_parameters["p2"] = RooRealVar('dijet6_p2','dijet6_p2',9.1,0.,60.)
		fit_parameters["p3"] = RooRealVar('dijet6_p3','dijet6_p3',0.5,-10.,10.)
		fit_parameters["p4"] = RooRealVar('dijet6_p4','dijet6_p4',0.5,-10.,10.)
		fit_parameters["p5"] = RooRealVar('dijet6_p5','dijet6_p5',0.5,-10.,10.)
		background_pdf = RooGenericPdf('background_' + function_name,'(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)+@4*pow(log(@0/%.1f), 2)+@5*pow(log(@0/%.1f), 3)))'%(collision_energy,collision_energy,collision_energy,collision_energy,collision_energy),RooArgList(mjj,fit_parameters["p1"],fit_parameters["p2"],fit_parameters["p3"],fit_parameters["p4"],fit_parameters["p5"]))
	elif function_name == "dijet4_OffB4":
		fit_parameters["p1"] = RooRealVar('dijet4_OffB4_p1','dijet4_OffB4_p1',-1.,-50.,50.)
		fit_parameters["p2"] = RooRealVar('dijet4_OffB4_p2','dijet4_OffB4_p2',9.1,0.,60.)
		fit_parameters["p3"] = RooRealVar('dijet4_OffB4_p3','dijet4_OffB4_p3',0.5,-10.,10.)
		background_pdf = RooGenericPdf('background_' + function_name,'(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)))'%(collision_energy,collision_energy,collision_energy),RooArgList(mjj,fit_parameters["p1"],fit_parameters["p2"],fit_parameters["p3"]))
	elif function_name == "dijet4_OffB5":
		fit_parameters["p1"] = RooRealVar('dijet4_OffB5_p1','dijet4_OffB5_p1',-1.,-50.,50.)
		fit_parameters["p2"] = RooRealVar('dijet4_OffB5_p2','dijet4_OffB5_p2',9.1,0.,60.)
		fit_parameters["p3"] = RooRealVar('dijet4_OffB5_p3','dijet4_OffB5_p3',0.5,-10.,10.)
		fit_parameters["p4"] = RooRealVar('dijet4_OffB5_p4','dijet4_OffB5_p4',-11.043684776942095, -20., 0.)
		background_pdf = RooGenericPdf('background_' + function_name,'(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f))) * (1+@4*@0/%.1f)'%(collision_energy,collision_energy,collision_energy, collision_energy),RooArgList(mjj,fit_parameters["p1"],fit_parameters["p2"],fit_parameters["p3"], fit_parameters["p4"]))
	elif function_name == "dijet4_OffB6":
		fit_parameters["p1"] = RooRealVar('dijet4_OffB6_p1','dijet4_OffB6_p1',-1.,-50.,50.)
		fit_parameters["p2"] = RooRealVar('dijet4_OffB6_p2','dijet4_OffB6_p2',9.1,0.,60.)
		fit_parameters["p3"] = RooRealVar('dijet4_OffB6_p3','dijet4_OffB6_p3',0.5,-10.,10.)
		fit_parameters["p4"] = RooRealVar('dijet4_OffB6_p4','dijet4_OffB6_p4',-11.043684776942095, -20., 0.)
		fit_parameters["p5"] = RooRealVar('dijet4_OffB6_p5','dijet4_OffB6_p5',38.32538741315025, 0., 60.)
		background_pdf = RooGenericPdf('background_' + function_name,'(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f))) * (1+@4*@0/%.1f+@5*pow(@0/%.1f,2))'%(collision_energy,collision_energy,collision_energy, collision_energy, collision_energy),RooArgList(mjj,fit_parameters["p1"],fit_parameters["p2"],fit_parameters["p3"], fit_parameters["p4"], fit_parameters["p5"]))

	elif function_name == "rational2":
		fit_parameters["p1"] = RooRealVar('rational2_p1','rational2_p1',8., 0., 100.)
		fit_parameters["p2"] = RooRealVar('rational2_p2','rational2_p2',80, 0., 100.)

		background_pdf = RooGenericPdf('background_' + function_name,'(1/pow(1+(@2*@0/%.1f),@1))'%(collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"]))
	elif function_name == "rational3":
		fit_parameters["p1"] = RooRealVar('rational3_p3', 'rational3_p1', 8.0, 0., 100.)
		fit_parameters["p2"] = RooRealVar('rational3_p2', 'rational3_p2', 80, 0., 100.)
		fit_parameters["p3"] = RooRealVar('rational3_p3', 'rational3_p3', 0., -100., 100.)

		background_pdf = RooGenericPdf('background_' + function_name,'(1/pow(1+(@2*@0/%.1f)+(@3*pow(@0/%.1f,2)),@1))'%(collision_energy,collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"], fit_parameters["p3"]))
	elif function_name == "rational4":
		fit_parameters["p1"] = RooRealVar('rational4_p4', 'rational4_p1', 8.0, 0., 100.)
		fit_parameters["p2"] = RooRealVar('rational4_p2', 'rational4_p2', 80, 0., 100.)
		fit_parameters["p3"] = RooRealVar('rational4_p3', 'rational4_p3', 0., -100., 100.)
		fit_parameters["p4"] = RooRealVar('rational4_p4', 'rational4_p4', 0., -100., 100.)

		background_pdf = RooGenericPdf('background_' + function_name,'(1/pow(1+(@2*@0/%.1f)+(@3*pow(@0/%.1f,2))+(@4*pow(@0/%.1f,3)),@1))'%(collision_energy,collision_energy, collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"], fit_parameters["p3"], fit_parameters["p4"]))

	elif function_name == "modexp3":
		fit_parameters["p1"] = RooRealVar("modexp3_p1", "modexp3_p1", -1., -100., 100.)
		fit_parameters["p2"] = RooRealVar("modexp3_p2", "modexp3_p2", 1., -10., 10.)
		background_pdf = RooGenericPdf('background_' + function_name,'exp(@1*pow(@0/%.1f, @2))'%(collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"]))
	elif function_name == "modexp4":
		fit_parameters["p1"] = RooRealVar("modexp4_p1", "modexp4_p1", -1., -100., 100.)
		fit_parameters["p2"] = RooRealVar("modexp4_p2", "modexp4_p2", 1., -10., 10.)
		fit_parameters["p3"] = RooRealVar("modexp4_p3", "modexp4_p3", 1., -10., 10.)
		background_pdf = RooGenericPdf('background_' + function_name,'exp(@1*pow(@0/%.1f, @2) + @1*pow((1.-@0/%.1f), @3))'%(collision_energy, collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"], fit_parameters["p3"]))
	elif function_name == "modexp5":
		fit_parameters["p1"] = RooRealVar("modexp5_p1", "modexp5_p1", -1., -100., 100.)
		fit_parameters["p2"] = RooRealVar("modexp5_p2", "modexp5_p2", 1., -10., 10.)
		fit_parameters["p3"] = RooRealVar("modexp5_p3", "modexp5_p3", 0., -10., 10.)
		fit_parameters["p4"] = RooRealVar("modexp5_p4", "modexp5_p4", 0., -10., 10.)
		background_pdf = RooGenericPdf('background_' + function_name,'exp(@1*pow(@0/%.1f, @2) + @4*pow((1.-@0/%.1f), @3))'%(collision_energy, collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"], fit_parameters["p3"], fit_parameters["p4"]))

	elif function_name == "polyx5":
		fit_parameters["p1"] = RooRealVar("polyx5_p1", "polyx5_p1", 1., -100., 100.)
		fit_parameters["p2"] = RooRealVar("polyx5_p2", "polyx5_p2", 1., -10., 10.)
		fit_parameters["p3"] = RooRealVar("polyx5_p3", "polyx5_p3", 0., -10., 10.)
		fit_parameters["p4"] = RooRealVar("polyx5_p4", "polyx5_p4", 0., -10., 10.)
		background_pdf = RooGenericPdf('background_' + function_name,'pow(1-(@0/%.1f), @1)*(1+@4*(@0/%.1f))/pow((@0/%.1f), @2+@3*log((@0/%.1f)))'%(collision_energy, collision_energy, collision_energy, collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"], fit_parameters["p3"], fit_parameters["p4"]))
	elif function_name == "polyx6":
		fit_parameters["p1"] = RooRealVar("polyx6_p1", "polyx6_p1", 1., -100., 100.)
		fit_parameters["p2"] = RooRealVar("polyx6_p2", "polyx6_p2", 1., -10., 10.)
		fit_parameters["p3"] = RooRealVar("polyx6_p3", "polyx6_p3", 0., -10., 10.)
		fit_parameters["p4"] = RooRealVar("polyx6_p4", "polyx6_p4", 0., -10., 10.)
		fit_parameters["p5"] = RooRealVar("polyx6_p5", "polyx6_p5", 0., -10., 10.)
		background_pdf = RooGenericPdf('background_' + function_name,'pow(1-(@0/%.1f), @1)*(1+@4*(@0/%.1f)+@5*pow(@0/%.1f,2))/pow((@0/%.1f), @2+@3*log((@0/%.1f)))'%(collision_energy, collision_energy, collision_energy, collision_energy, collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"], fit_parameters["p3"], fit_parameters["p4"], fit_parameters["p5"]))
	elif function_name == "polyx7":
		fit_parameters["p1"] = RooRealVar("polyx7_p1", "polyx7_p1", 1., -100., 100.)
		fit_parameters["p2"] = RooRealVar("polyx7_p2", "polyx7_p2", 1., -10., 10.)
		fit_parameters["p3"] = RooRealVar("polyx7_p3", "polyx7_p3", 0., -10., 10.)
		fit_parameters["p4"] = RooRealVar("polyx7_p4", "polyx7_p4", 0., -10., 10.)
		fit_parameters["p5"] = RooRealVar("polyx7_p5", "polyx7_p5", 0., -10., 10.)
		fit_parameters["p6"] = RooRealVar("polyx7_p6", "polyx7_p6", 0., -10., 10.)
		background_pdf = RooGenericPdf('background_' + function_name,'pow(1-(@0/%.1f), @1)*(1+@4*(@0/%.1f)+@5*pow(@0/%.1f,2)+@6*pow(@0/%.1f,3))/pow((@0/%.1f), @2+@3*log((@0/%.1f)))'%(collision_energy, collision_energy, collision_energy, collision_energy, collision_energy, collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"], fit_parameters["p3"], fit_parameters["p4"], fit_parameters["p5"], fit_parameters["p6"]))

	elif function_name == "atlas4":
		fit_parameters["p1"] = RooRealVar("atlas4_p1", "atlas4_p1", 1., -100., 100.)
		fit_parameters["p2"] = RooRealVar("atlas4_p2", "atlas4_p2", 1., -10., 10.)
		fit_parameters["p3"] = RooRealVar("atlas4_p3", "atlas4_p3", 0., -10., 10.)
		background_pdf = RooGenericPdf('background_' + function_name, 
			'exp(-1.*(@2*(@0/%.1f)+@3*pow((@0/%.1f),2)))/(pow((@0/%.1f),@1))'
				%(collision_energy,collision_energy,collision_energy), 
			RooArgList(mjj, fit_parameters["p1"],fit_parameters["p2"],fit_parameters["p3"]))
	elif function_name == "atlas5":
		fit_parameters["p1"] = RooRealVar("atlas5_p1", "atlas5_p1", 1., -100., 100.)
		fit_parameters["p2"] = RooRealVar("atlas5_p2", "atlas5_p2", 1., -10., 10.)
		fit_parameters["p3"] = RooRealVar("atlas5_p3", "atlas5_p3", 0., -10., 10.)
		fit_parameters["p4"] = RooRealVar("atlas5_p4", "atlas5_p4", 0., -10., 10.)
		background_pdf = RooGenericPdf('background_' + function_name, 
			'exp(-1.*(@2*(@0/%.1f)+@3*pow((@0/%.1f),2)+@4*pow((@0/%.1f),3)))/(pow((@0/%.1f),@1))'
				%(collision_energy,collision_energy,collision_energy,collision_energy), 
			RooArgList(mjj, fit_parameters["p1"],fit_parameters["p2"],fit_parameters["p3"],fit_parameters["p4"]))
	elif function_name == "atlas6":
		fit_parameters["p1"] = RooRealVar("atlas6_p1", "atlas6_p1", 1., -100., 100.)
		fit_parameters["p2"] = RooRealVar("atlas6_p2", "atlas6_p2", 1., -10., 10.)
		fit_parameters["p3"] = RooRealVar("atlas6_p3", "atlas6_p3", 0., -10., 10.)
		fit_parameters["p4"] = RooRealVar("atlas6_p4", "atlas6_p4", 0., -10., 10.)
		fit_parameters["p5"] = RooRealVar("atlas6_p5", "atlas6_p5", 0., -10., 10.)
		background_pdf = RooGenericPdf('background_' + function_name, 
			'exp(-1.*(@2*(@0/%.1f)+@3*pow((@0/%.1f),2)+@4*pow((@0/%.1f),3)+@5*pow((@0/%.1f),4)))/(pow((@0/%.1f),@1))'
				%(collision_energy,collision_energy,collision_energy,collision_energy,collision_energy), 
			RooArgList(mjj, fit_parameters["p1"],fit_parameters["p2"],fit_parameters["p3"],fit_parameters["p4"],fit_parameters["p5"]))

	elif function_name == "polypower3":
		fit_parameters["p1"] = RooRealVar('polypower3_p1','polypower3_p1',8., 0., 40.)
		fit_parameters["p2"] = RooRealVar('polypower3_p2','polypower3_p2',80, -100., 100.)
		background_pdf = RooGenericPdf('background_' + function_name,'(1/pow(max(1+(@2*@0/%.1f), 1.e-5),@1))'%(collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"]))
	elif function_name == "polypower4":
		fit_parameters["p1"] = RooRealVar('polypower4_p1', 'polypower4_p1', 8.0, 0., 40.)
		fit_parameters["p2"] = RooRealVar('polypower4_p2', 'polypower4_p2', 80, -100., 100.)
		fit_parameters["p3"] = RooRealVar('polypower4_p3', 'polypower4_p3', 0., -100., 100.)
		background_pdf = RooGenericPdf('background_' + function_name,'(1/pow(max(1+(@2*@0/%.1f)+(@3*pow(@0/%.1f,2)), 1.e-5),@1))'%(collision_energy,collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"], fit_parameters["p3"]))
	elif function_name == "polypower5":
		fit_parameters["p1"] = RooRealVar('polypower5_p1', 'polypower5_p1', 8.0, 0., 40.)
		fit_parameters["p2"] = RooRealVar('polypower5_p2', 'polypower5_p2', 80, -100., 100.)
		fit_parameters["p3"] = RooRealVar('polypower5_p3', 'polypower5_p3', 0., -100., 100.)
		fit_parameters["p4"] = RooRealVar('polypower5_p4', 'polypower5_p4', 0., -100, 100.)
		background_pdf = RooGenericPdf('background_' + function_name,'(1/pow(max(1+(@2*@0/%.1f)+(@3*pow(@0/%.1f,2))+(@4*pow(@0/%.1f,3)), 1.e-5),@1))'%(collision_energy,collision_energy,collision_energy), RooArgList(mjj, fit_parameters["p1"], fit_parameters["p2"], fit_parameters["p3"],fit_parameters["p4"]))
	elif function_name == "onethird3":
		fit_parameters["p1"] = RooRealVar('onethird3_p1','onethird3_p1',4.8,-100., 100.)
		fit_parameters["p2"] = RooRealVar('onethird3_p2','onethird3_p2',7., -100., 100.)
		background_pdf = RooGenericPdf('background_' + function_name,'(pow(@0/%.1f,-@1)*pow(1-pow(@0/%.1f,1/3),@2))'%(collision_energy, collision_energy),RooArgList(mjj,fit_parameters["p1"],fit_parameters["p2"]))

	else:
		print "[roofit_functions_huge.make_background_pdf] ERROR : Unrecognized fit function " + function_name
		sys.exit(1)
	return (background_pdf, fit_parameters)

initial_values = {
	"dijet":{
		"trigbbl_CSVTM":{
			"p1":-2.7601e+01,
			"p2": 1.7256e+01,
			"p3": 1.6284e+00,
		}, "trigbbh_CSVTM":{
			"p1":-1.3545e+01,
			"p2": 1.4047e+01,
			"p3": 1.2425e+00,
		}
	}, 
	"rational":{
		"trigbbl_CSVTM":{
			"p1": 2.7487e+01,
			"p2":-4.5282e+01,
			"p3": 1.1914e+01,
		}, 
		"trigbbh_CSVTM":{
			"p1": 3.1388e+01,
			"p2":-2.2446e+01,
			"p3": 9.9301e+00,
		}
	},
	"modexp":{
		"trigbbl_CSVTM":{
			"p1":-4.8711e+01,
			"p2": 3.9427e-01,
		},
		"trigbbh_CSVTM":{
			"p1":1.7639e+01,
			"p2":7.6084e+00,
		}
	},
	"polyx":{
		"trigbbl_CSVTM":{
			"p1":-8.8328e-01,
			"p2": 9.9998e+00,
			"p3": 6.8157e-01,
			"p4":-3.5396e+00,
			"p5": 9.9897e+00,
			"p6": 9.8764e+00,
		}, 
		"trigbbh_CSVTM":{
			"p1":3.9829e+00,
			"p2":8.0021e+00,
			"p3":2.5962e-01,
			"p4":9.9818e+00,
		}
	}, 
	"polypower":{
		"trigbbl_CSVTM":{
			"p1": 1.1914e+01,
			"p2": 2.7487e+01,
			"p3":-4.5282e+01,
		},
		"trigbbh_CSVTM":{
			"p1": 9.9301e+00,
			"p2": 3.1388e+01,
			"p3":-2.2446e+01,
		}
	},
	"atlas":{
		"trigbbl_CSVTM":{
			"p1":5.4787e+00,
			"p2":9.9980e+00,
			"p3":1.0000e+01,
		}, 
		"trigbbh_CSVTM":{
			"p1": 5.8752e+00,
			"p2": 9.9995e+00,
			"p3":-6.2029e+00,
		}
	}
}

def set_initial_values(parameters, function_name, analysis):
	# Keep track of whether we have found an initial value for this parameter.
	# If not, it will be set to zero later.
	found_initial_value = {}
	for par_name, par_rv in parameters.iteritems():
		found_initial_value[par_name] = False

	for dict_function_name in initial_values:
		if dict_function_name in function_name:
			if analysis in initial_values[dict_function_name]:
				for par_name, par_rv in parameters.iteritems():
					print initial_values[dict_function_name][analysis]
					if par_name in initial_values[dict_function_name][analysis]:
						par_rv.setVal(initial_values[dict_function_name][analysis][par_name])
						found_initial_value[par_name] = True

	for par_name, par_rv in parameters.iteritems():
		if not found_initial_value[par_name]:
			par_rv.setVal(0.)
