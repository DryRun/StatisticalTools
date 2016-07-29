#! /usr/bin/env python

import os
import sys
import array
from array import array
import ROOT
from ROOT import *
gROOT.SetBatch(True)
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
seaborn = Root.SeabornInterface()
seaborn.Initialize()

if not os.path.exists(os.path.join(os.environ["CMSSW_BASE"],"src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer_cc.so")):
	ROOT.gROOT.ProcessLine(".L " + os.path.join(os.environ["CMSSW_BASE"],"src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer.cc")+"+")
else:
	ROOT.gSystem.Load(os.path.join(os.environ["CMSSW_BASE"],"src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer_cc.so"))
from ROOT import PdfDiagonalizer
RooMsgService.instance().setSilentMode(kTRUE)
RooMsgService.instance().setStreamStatus(0,kFALSE)
RooMsgService.instance().setStreamStatus(1,kFALSE)

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_7_4_15/python/CMSDijet/StatisticalTools")
import limitPaths as limit_paths

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/CMSDIJET/QCDAnalysis")
import analysis_configuration_8TeV as analysis_config
import simulation_configuration_8TeV as simulation_config

def BackgroundFit_f1(x, par):
	return par[0] * (1. - (x[0] / 8.e3))**par[1] / ((x[0] / 8.e3)**(par[2] + par[3] * TMath.Log((x[0] / 8.e3))))

def BackgroundFit_f2(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3))**par[2]

def BackgroundFit_f3(x, par):
	return par[0] / (par[1] + (x[0] / 8.e3))**par[2]

def BackgroundFit_f4(x, par):
	return par[0] / ((par[1] + par[2]*x[0]/8.e3 + (x[0]/8.e3)**2)**par[3])

def BackgroundFit_f5(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3)**(1./3.))**par[2]

def CrystalBallFit(x, par):
	#Crystal ball function for signal, parameters are 0:alpha,1:n,2:mean,3:sigma,4:normalization;
	t = (x[0]-par[2])/par[3]
	if par[0] < 0:
		t = -t
	absAlpha = TMath.Abs(par[0])
	if (t >= -1. * absAlpha):
		return par[4]*TMath.Exp(-0.5*t*t)
	else:
		a =  TMath.Power(par[1]/absAlpha,par[1])*TMath.Exp(-0.5*absAlpha*absAlpha);
		b= par[1]/absAlpha - absAlpha; 
		return par[4]*(a/TMath.Power(b - t, par[1]))

class MjjFit:
	def __init__(self, mass_range=[0., 2000.]):
		self.data_histogram_ = None
		self.data_roohistogram_ = None

		self.backgrounds_ = []
		self.background_histograms_ = {}
		self.background_initial_normalizations_ = {}
		self.background_dnorm_ = {}

		self.signals_ = []
		self.signal_histograms_ = {}
		self.signal_roohistograms_ = {}
		#self.signal_names_ = []
		self.signal_initial_normalizations_ = {}

		self.luminosity_ = 0.
		self.collision_energy_ = 8000.

		# Fit storage
		self.simple_fit_ = None
		self.mjj_ = RooRealVar('mjj','mjj',float(mass_range[0]),float(mass_range[1]))
		self.fit_range_ = mass_range
		#self.mjj_.setRange("fit_range", float(mass_range[0]), float(mass_range[1]))
		#self.workspace_ = None
		self.saved_workspaces_ = {}
		self.fit_parameters = {}

		# Trigger correction
		self.trigger_correction_applied_ = False
		self.trigbbl_efficiency_ = TF1("trigbbl_efficiency", "(1. / (1. + TMath::Exp(-1. * (x - [0]) / [1])))**[2]", 175, 400)
		self.trigbbl_efficiency_.SetParameter(0, 1.82469e+02)
		self.trigbbl_efficiency_.SetParameter(1,  2.87768e+01)
		self.trigbbl_efficiency_.SetParameter(2,  9.11659e-01)

	def set_fit_range(self, fit_range):
		self.fit_range_ = fit_range

	def add_data(self, data_histogram, luminosity):
		print "[MjjFit.add_data] INFO : Adding data histogram"
		# Add a data histogram
		self.data_histogram_ = data_histogram.Clone()
		self.data_histogram_.SetDirectory(0)
		self.luminosity_ = luminosity

	def add_background(self, background_histogram, background, normalization_uncertainty = 0):
		if self.trigger_correction_applied_:
			raise Exception("Attempt to subtract background after correcting for trigger. For consistency, all backgrounds must be subtracted before correcting.")
		self.backgrounds_.append(background)
		self.background_histograms_[background] = background_histogram.Clone()
		self.background_histograms_[background].SetDirectory(0)
		self.background_initial_normalizations_[background] = self.background_histograms_[background].Integral()
		self.background_dnorm_[background] = normalization_uncertainty

	def correct_trigger(self):
		self.trigger_correction_applied_ = True
		print "[MjjFit.correct_trigger] INFO : Correcting trigger efficiency using Jet60_Jet53"
		if not self.data_histogram_:
			print "[MjjFit.correct_trigger] ERROR : Call to correct_trigger before add_data"
			sys.exit(1)
		for bin in xrange(1, self.data_histogram_.GetNbinsX() + 1):
			bin_center = self.data_histogram_.GetXaxis().GetBinCenter(bin)
			if bin_center > self.trigbbl_efficiency_.GetXmin() and bin_center < self.trigbbl_efficiency_.GetXmax():
				self.data_histogram_.SetBinContent(bin, self.data_histogram_.GetBinContent(bin) / self.trigbbl_efficiency_.Eval(bin_center))
				self.data_histogram_.SetBinError(bin, self.data_histogram_.GetBinError(bin) / self.trigbbl_efficiency_.Eval(bin_center))

				for background in self.backgrounds_:
					self.background_histograms_[background].SetBinContent(bin, self.background_histograms_[background].GetBinContent(bin) / self.trigbbl_efficiency_.Eval(bin_center))
					self.background_histograms_[background].SetBinError(bin, self.background_histograms_[background].GetBinError(bin) / self.trigbbl_efficiency_.Eval(bin_center))
					# Recompute background normalization!
					self.background_initial_normalizations_[background] = self.background_histograms_[background].Integral()

				for signal in self.signals_:
					self.signal_histograms_[background].SetBinContent(bin, self.signal_histograms_[background].GetBinContent(bin) / self.trigbbl_efficiency_.Eval(bin_center))
					self.signal_histograms_[background].SetBinError(bin, self.signal_histograms_[background].GetBinError(bin) / self.trigbbl_efficiency_.Eval(bin_center))

	def add_signal(self, signal_name, signal_histogram):
		print "[MjjFit.add_signal] INFO : Adding signal histogram " + signal_name
		if self.trigger_correction_applied_:
			raise Exception("Attempt to subtract background after correcting for trigger. For consistency, all backgrounds must be subtracted before correcting.")
		# Add a signal histogram. 
		# Scale to sigma=1 (in whatever units the data luminosity is given in), so the 'r' parameter corresponds to the limit on the cross section.
		if self.luminosity_ == 0:
			print "[MjjFit.add_signal] ERROR : Please set luminosity first (MjjFit.set_luminosity(###))."
			sys.exit(1)
		self.signals_.append(signal_name)
		self.signal_histograms_[signal_name] = signal_histogram.Clone()
		self.signal_histograms_[signal_name].SetDirectory(0)
		self.signal_initial_normalizations_[signal_name] = self.signal_histograms_[signal_name].Integral()
		self.signal_histograms_[signal_name].Scale(1. * self.luminosity_ / self.signal_histograms_[signal_name].Integral())
		self.signal_roohistograms_[signal_name] = RooDataHist(signal_histogram.GetName() + "_rdh", signal_histogram.GetName() + "_rdh", RooArgList(self.mjj_), self.signal_histograms_[signal_name])
		self.signal_roohistograms_[signal_name].Print()

	def simple_fit_B(self, fit_range=None):
		# Run a simple ROOT fit (B only)
		# - No S+B equivalent. This is done with RooFit.
		if fit_range:
			fit_min = fit_range[0]
			fit_max = fit_range[1]
		else:
			fit_min = data_histogram.GetXaxis().GetXmin()
			fit_max = data_histogram.GetXaxis().GetXmax()
		fit = TF1(data_histogram.GetName() + "_fit_B", BackgroundFit_f0, fit_min, fit_max, 4)
		fit.SetParameter(0, 2.e-4)
		fit.SetParameter(1, 3)
		fit.SetParameter(2, 10)
		fit.SetParameter(3, 1)
		fit.SetParLimits(0, 1.e-6, 1.e2)
		fit.SetParLimits(1, -25., 25.)
		fit.SetParLimits(2, -25., 25.)
		fit.SetParLimits(3, -5., 5.)
		data_histogram.Fit(fit, "ER0I")
		fit_ratio = self.make_fit_pull_histogram(data_histogram, fit)
		print "Fit chi2/ndf = " + str(fit.GetChisquare()) + " / " + str(fit.GetNDF()) + " = " + str(fit.GetChisquare() / fit.GetNDF())
		self.simple_fit_ = {"fit":fit, "fit_ratio":fit_ratio}
		return self.simple_fit_

	def make_fit_pull_histogram(self, hist, fit):
		#print "Fit xmin = " + str(fit.GetXmin())
		hist_ratio = hist.Clone()
		hist_ratio.SetName(hist.GetName() + "_fit_ratio")
		for bin in xrange(1, hist_ratio.GetNbinsX() + 1):
			xmin = hist_ratio.GetXaxis().GetBinLowEdge(bin)
			xmax = hist_ratio.GetXaxis().GetBinUpEdge(bin)
			if xmax < fit.GetXmin() or xmin > fit.GetXmax():
				hist_ratio.SetBinContent(bin, 0.)
				hist_ratio.SetBinError(bin, 0.)
				continue
			fit_integral = fit.Integral(xmin, xmax)
			if hist.GetBinError(bin) > 0:
				hist_ratio.SetBinContent(bin, (hist.GetBinContent(bin) * hist.GetBinWidth(bin) - fit_integral) / (hist.GetBinError(bin) * hist.GetBinWidth(bin)))
				hist_ratio.SetBinError(bin, 0.)
			else:
				hist_ratio.SetBinContent(bin, 0.)
				hist_ratio.SetBinError(bin, 0.)
		return hist_ratio


	def make_background_pdf(self, fit_function):
		self.fit_parameters[fit_function] = {}
		if fit_function == "f1":
			self.fit_parameters[fit_function]["p1"] = RooRealVar('p1','p1',5.e-3,0.,100.)
			self.fit_parameters[fit_function]["p2"] = RooRealVar('p2','p2',9.1,0.,60.)
			self.fit_parameters[fit_function]["p3"] = RooRealVar('p3','p3',0.5,-10.,10.)
			background_pdf = RooGenericPdf('background_pdf','(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)))'%(self.collision_energy_,self.collision_energy_,self.collision_energy_),RooArgList(self.mjj_,self.fit_parameters[fit_function]["p1"],self.fit_parameters[fit_function]["p2"],self.fit_parameters[fit_function]["p3"]))
		elif fit_function == "f2":
			self.fit_parameters[fit_function]["p1"] = RooRealVar('p1','p1',5.6,0.,100.)
			self.fit_parameters[fit_function]["p2"] = RooRealVar('p2','p2',10.,0.,60.)
			background_pdf = RooGenericPdf('background_pdf','(pow(@0/%.1f,-@1)*pow(1-@0/%.1f,@2))'%(self.collision_energy_, self.collision_energy_), RooArgList(self.mjj_, self.fit_parameters[fit_function]["p1"], self.fit_parameters[fit_function]["p2"]))
		elif fit_function == "f3":
			self.fit_parameters[fit_function]["p1"] = RooRealVar('p1','p1',0.016, -50., 50.)
			self.fit_parameters[fit_function]["p2"] = RooRealVar('p2','p2',8., -50., 50.)
			background_pdf = RooGenericPdf('background_pdf','(1/pow(@1+@0/%.1f,@2))'%(self.collision_energy_), RooArgList(self.mjj_, self.fit_parameters[fit_function]["p1"], self.fit_parameters[fit_function]["p2"]))
		elif fit_function == "f4":
			self.fit_parameters[fit_function]["p1"] = RooRealVar('p1', 'p1', 22, -50., 50.)
			self.fit_parameters[fit_function]["p2"] = RooRealVar('p2', 'p2', 14.1, -50., 50.)
			self.fit_parameters[fit_function]["p3"] = RooRealVar('p3', 'p3', 8.0,-50., 50.)
			background_pdf = RooGenericPdf('background_pdf','(1/pow(@1+@2*@0/%.1f+pow(@0/%.1f,2),@3))'%(self.collision_energy_,self.collision_energy_), RooArgList(self.mjj_, self.fit_parameters[fit_function]["p1"], self.fit_parameters[fit_function]["p2"], self.fit_parameters[fit_function]["p3"]))
		elif fit_function == "f5":
			self.fit_parameters[fit_function]["p1"] = RooRealVar('p1','p1',4.8,-50., 50.)
			self.fit_parameters[fit_function]["p2"] = RooRealVar('p2','p2',7., -50., 50.)
			background_pdf = RooGenericPdf('background_pdf','(pow(@0/%.1f,-@1)*pow(1-pow(@0/%.1f,1/3),@2))'%(self.collision_energy_, self.collision_energy_),RooArgList(self.mjj_,self.fit_parameters[fit_function]["p1"],self.fit_parameters[fit_function]["p2"]))
		else:
			print "[MjjFit::make_background_pdf] ERROR : Unrecognized fit function " + fit_function
			sys.exit(1)
		return background_pdf

	def make_background_tf1(self, fit_function, mjj_range=None):
		if mjj_range:
			mjj_min = mjj_range[0]
			mjj_max = mjj_range[1]
		else:
			mjj_min = self.mjj_.getMin()
			mjj_max = self.mjj_.getMax()

		if fit_function == "f1":
			background_tf1 = TF1('background_tf1_f1',BackgroundFit_f1, mjj_min, mjj_max, 4)
		elif fit_function == "f2":
			background_tf1 = TF1('background_tf1_f2',BackgroundFit_f2, mjj_min, mjj_max, 3)
		elif fit_function == "f3":
			background_tf1 = TF1('background_tf1_f3',BackgroundFit_f3, mjj_min, mjj_max, 3)
		elif fit_function == "f4":
			background_tf1 = TF1('background_tf1_f4',BackgroundFit_f4, mjj_min, mjj_max, 4)
		elif fit_function == "f5":
			background_tf1 = TF1('background_tf1_f5',BackgroundFit_f5, mjj_min, mjj_max, 3)
		else:
			print "[MjjFit::make_background_tf1] ERROR : Unrecognized fit function " + fit_function
			sys.exit(1)
		return background_tf1

	def make_background_tf1_from_roofitresult(self, fit_function, roofitresult, mjj_range=None):
		print "[MjJFit::make_background_tf1_from_roofitresult] INFO : make_background_tf1_from_roofitresult for fit function " + fit_function
		print mjj_range
		function = self.make_background_tf1(fit_function, mjj_range)
		for i in xrange(function.GetNpar()):
			function.SetParameter(i, roofitresult.floatParsFinal()[i].getVal())
			function.SetParError(i, roofitresult.floatParsFinal()[i].getError())
		return function


	def fit(self, fit_function, save_to, fit_options, signal_name=None, fit_strategy=1):
		print "[MjjFit::fit] INFO : Fitting with function " + fit_function + ", save_to = " + save_to
		# Run a RooFit fit

		# Make data roohistogram
		self.data_roohistogram_ = RooDataHist('data_roohistogram','data_roohistogram',RooArgList(self.mjj_),self.data_histogram_)
		self.data_roohistogram_.Print()

		# Create background PDF
		background_pdfs = {}
		background_epdfs = {}
		background_rdhs = {}
		background_normalizations = {}
		background_constraints = {}
		background_epdf_list = RooArgList("background_epdfs")
		background_constraint_list = RooArgSet("background_constraints")
		for background in self.backgrounds_:
			print "[debug] Creating background RooHistPdf for background " + background
			print "[debug] \tInput = histogram ",
			print self.background_histograms_[background],
			print " with integral " + str(self.background_histograms_[background].Integral())
			background_rdhs[background] = RooDataHist(self.background_histograms_[background].GetName() + "_rdf", self.background_histograms_[background].GetName() + "_rdh", RooArgList(self.mjj_), self.background_histograms_[background])
			background_pdfs[background]           = RooHistPdf(background + "_pdf", background + "_pdf", RooArgSet(self.mjj_), background_rdhs[background])
			background_normalizations[background] = RooRealVar(background + "_norm", background + "_norm", self.background_initial_normalizations_[background], 0., self.background_initial_normalizations_[background] * 2.)
			background_epdfs[background]          = RooExtendPdf(background + "_epdf", background + "_epdf", background_pdfs[background], background_normalizations[background])
			background_epdf_list.add(background_epdfs[background])
			if self.background_dnorm_[background] > 0:
				background_constraints[background] = RooGaussian(background + "_constraint", background + "_constraint", background_normalizations[background], RooFit.RooConst(self.background_initial_normalizations_[background]), RooFit.RooConst(self.background_initial_normalizations_[background] * self.background_dnorm_[background]))
				background_constraint_list.add(background_constraints[background])
			else:
				background_normalizations[background].setConstant(True)
		background_pdfs["fit"] = self.make_background_pdf(fit_function)
		data_integral = self.data_histogram_.Integral(self.data_histogram_.GetXaxis().FindBin(float(fit_options["min_mjj"])),self.data_histogram_.GetXaxis().FindBin(float(fit_options["max_mjj"])))
		background_normalizations["fit"] = RooRealVar('fit_norm','fit_norm',data_integral,0.,data_integral * 1.e4)
		background_epdfs["fit"] = RooExtendPdf("fit_epdf", "fit_epdf", background_pdfs["fit"], background_normalizations["fit"])
		background_epdf_list.add(background_epdfs["fit"])

		background_model = RooAddPdf("total_background", "b", background_epdf_list)

		# Create signal PDF and fit model
		if signal_name:
			signal_pdf = RooHistPdf('signal_pdf', 'signal_pdf', RooArgSet(self.mjj_), self.signal_roohistograms_[signal_name])
			#signal_pdf.Print()
			signal_norm = RooRealVar('signal_norm','signal_norm',0,-1e+05,1e+05)
			signal_epdf = RooExtendPdf("signal_epdf", "signal_epdf", signal_pdf, signal_norm)
			#signal_norm.Print()
			signal_initial_norm = RooRealVar('signal_initial_norm', 'signal_initial_norm', self.signal_initial_normalizations_[signal_name])
			model = RooAddPdf("model_" + fit_function,"s+b",RooArgList(background_model,signal_epdf))
		else:
			model = background_model
			#model = RooAddPdf("model_" + fit_function,"b",RooArgList(background_model))

		# Run fit
		res = model.fitTo(self.data_roohistogram_, RooFit.Save(kTRUE), RooFit.Strategy(fit_strategy), RooFit.Range(float(self.fit_range_[0]), float(self.fit_range_[1])), RooFit.ExternalConstraints(background_constraint_list), RooFit.Verbose(True))

		# Print some results
		print "Background normalizations:"
		for background, var in background_normalizations.iteritems():
			var.Print()
		print "Constraints:"
		for background, var in background_constraints.iteritems():
			var.Print()


		# Save to workspace
		self.workspace_ = RooWorkspace('w','workspace')
		#getattr(w,'import')(background,ROOT.RooCmdArg())
		print model
		model.Print()
		getattr(self.workspace_, 'import')(model, RooFit.Rename("total_background_model"))
		for background in background_epdfs.keys():
			getattr(self.workspace_,'import')(background_epdfs[background],RooFit.Rename("background_epdf_" + background))
		for background in background_normalizations.keys():
			getattr(self.workspace_,'import')(background_normalizations[background],ROOT.RooCmdArg())
		for background in self.backgrounds_:
			getattr(self.workspace_,'import')(RooRealVar(background + '_initial_norm', 'background_' + background + '_initial_norm', self.background_initial_normalizations_[background], 0., self.background_initial_normalizations_[background]*1.e3), ROOT.RooCmdArg())
		for background in background_constraints.keys():
			getattr(self.workspace_,'import')(background_constraints[background],ROOT.RooCmdArg())
		getattr(self.workspace_,'import')(self.data_roohistogram_,RooFit.Rename("data_obs"))
		getattr(self.workspace_, 'import')(res)
		if signal_name:
			getattr(self.workspace_,'import')(self.signal_roohistograms_[signal_name],RooFit.Rename("signal"))
			getattr(self.workspace_,'import')(signal_pdf,RooFit.Rename("signal_pdf"))
			getattr(self.workspace_,'import')(signal_norm,ROOT.RooCmdArg())
			getattr(self.workspace_,'import')(signal_initial_norm,RooFit.Rename("signal_initial_norm"))
	
		self.workspace_.Print()
		self.workspace_.writeToFile(save_to)
		if signal_name:
			self.saved_workspaces_[signal_name] = save_to
		else:
			self.saved_workspaces_["background"] = save_to

		print "[MjjFit::fit] INFO : Done with fit."

	def rooplot(self, save_tag, fit_functions, background_workspaces, fitted_signal_workspaces=None, expected_signal_workspaces=None, log=False, x_range=None, data_binning=None, normalization_bin_width=1):
		print "RooPlotting " + save_tag
		c = TCanvas("c_" + save_tag, "c_" + save_tag, 800, 1200)
		l = TLegend(0.55, 0.6, 0.88, 0.88)
		l.SetFillColor(0)
		l.SetBorderSize(0)
		top = TPad("top", "top", 0., 0.5, 1., 1.)
		top.SetBottomMargin(0.03)
		top.Draw()
		if log:
			top.SetLogy()
		c.cd()
		bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
		bottom.SetTopMargin(0.02)
		bottom.SetBottomMargin(0.2)
		bottom.Draw()
		ROOT.SetOwnership(c, False)
		ROOT.SetOwnership(top, False)
		ROOT.SetOwnership(bottom, False)
		top.cd()

		# Make data histogram
		data_hist_raw = self.data_histogram_.Clone()

		if data_binning:
			data_hist = data_hist_raw.Rebin(len(data_binning) - 1, "data_hist_rebinned", data_binning)
			for bin in xrange(1, data_hist.GetNbinsX() + 1):
				data_hist.SetBinContent(bin, data_hist.GetBinContent(bin) / data_hist.GetXaxis().GetBinWidth(bin) * normalization_bin_width)
				data_hist.SetBinError(bin, data_hist.GetBinError(bin) / data_hist.GetXaxis().GetBinWidth(bin) * normalization_bin_width)
			roobinning = RooBinning(len(data_binning) - 1, data_binning, "data_binning")
		else:
			data_hist = data_hist_raw
			bins_array = array('d', [])
			for bin in xrange(1, data_hist.GetNbinsX() + 1):
				bins_array.append(data_hist.GetXaxis().GetBinLowEdge(bin))
			bins_array.append(data_hist.GetXaxis().GetBinUpEdge(data_hist.GetNbinsX()))
			roobinning = RooBinning(len(bins_array) - 1, bins_array, "data_binning")

		if x_range:
			x_min = x_range[0]
			x_max = x_range[1]
		else:
			x_min = self.mjj_.GetMin()
			x_max = self.mjj_.GetMax()
		frame_top = TH1D("frame_top", "frame_top", 100, x_min, x_max)
		if log:
			frame_top.SetMaximum(data_hist.GetMaximum() * 50.)
			frame_top.SetMinimum(0.1)
		else:
			frame_top.SetMaximum(data_hist.GetMaximum() * 1.3)
			frame_top.SetMinimum(0.)
		frame_top.GetYaxis().SetTitle("Events / " + str(normalization_bin_width) + " GeV")
		frame_top.GetXaxis().SetTitleSize(0)
		frame_top.GetXaxis().SetLabelSize(0)
		#if log:
		#	frame_top.SetMaximum(y_max * 10.)
		#	frame_top.SetMinimum(y_min / 100.)
		#else:
		#	frame_top.SetMaximum(y_max * 1.3)
		#	frame_top.SetMinimum(0.)
		frame_top.Draw()
		#data_rdh.plotOn(frame_top, RooFit.Name("Data"))
		data_hist.SetMarkerStyle(20)
		data_hist.SetMarkerColor(1)
		data_hist.SetMarkerSize(1)
		data_hist.Draw("p same")
		l.AddEntry(data_hist, "Data", "pl")



		background_histograms = {}
		background_histograms_fitted_range = {}
		mc_stack = THStack("mc_stack", "")
		style_counter = 0
		for fit_function in fit_functions:
			f_workspace = TFile(background_workspaces[fit_function], "READ")
			workspace = f_workspace.Get("w")

			# First time: scale the MC histograms
			if style_counter == 0:
				mc_stack = THStack("mc_stack", "")
				mc_counter = 0
				for mc_background in self.backgrounds_:
					mc_background_pdf = workspace.pdf(mc_background + "_epdf")
					background_histograms[background] = mc_background_pdf.createHistogram(mc_background, self.mjj_, RooFit.Binning(roobinning), RooFit.Extended(True))
					background_histograms[background].SetDirectory(0)
					background_histograms[background].SetFillColor(seaborn.GetColorRoot("default", mc_counter))
					mc_stack.Add(background_histograms[background])
					l.AddEntry(background_histograms[background], background, "f")
					mc_counter += 1
				mc_stack.Draw("hist same")

			# Make TF1 from fit background
			fitresult_name = "fitresult_total_background_data_roohistogram"
			fitresult = workspace.genobj(fitresult_name)
			fitresult.Print()
			background_tf1 = self.make_background_tf1_from_roofitresult(fit_function, fitresult, mjj_range=x_range)
			fit_normalization = workspace.var("fit_norm").getVal()
			scale_factor = fit_normalization / background_tf1.Integral(self.mjj_.getMin(), self.mjj_.getMax())
			background_tf1.SetParameter(0, background_tf1.GetParameter(0) * scale_factor)

			print "[MjjFit::plot] INFO : Normalizing fit to data over range [" + str(self.mjj_.getMin()) + ", " + str(self.mjj_.getMax()) + "]"
			data_integral = data_hist.Integral(data_hist.GetXaxis().FindBin(self.mjj_.getMin()), data_hist.GetXaxis().FindBin(self.mjj_.getMax()), "width")
			integration_xmin = data_hist.GetXaxis().GetBinLowEdge(data_hist.GetXaxis().FindBin(self.mjj_.getMin()))
			integration_xmax = data_hist.GetXaxis().GetBinUpEdge(data_hist.GetXaxis().FindBin(self.mjj_.getMax()))
			scale_factor = data_integral / background_tf1s[fit_function].Integral(integration_xmin, integration_xmax)
			background_tf1s[fit_function].SetParameter(0, background_tf1s[fit_function].GetParameter(0) * scale_factor)
			background_tf1s[fit_function].SetParError(0, background_tf1s[fit_function].GetParError(0) * scale_factor)
			background_tf1s[fit_function].SetLineColor(seaborn.GetColorRoot("default", style_counter))
			background_tf1s[fit_function].SetLineWidth(1)
			background_tf1s[fit_function].SetLineStyle(2)
			background_tf1s[fit_function].Draw("l same")

			background_tf1s_fitted_range[fit_function] = background_tf1s[fit_function].Clone()
			background_tf1s_fitted_range[fit_function].SetRange(self.mjj_.getMin(), self.mjj_.getMax())
			background_tf1s_fitted_range[fit_function].SetLineColor(seaborn.GetColorRoot("dark", style_counter))
			background_tf1s_fitted_range[fit_function].SetLineWidth(2)
			background_tf1s_fitted_range[fit_function].SetLineStyle(1)
			background_tf1s_fitted_range[fit_function].Draw("l same")




			background_pdf = workspace.pdf("total_background")
			background_histograms[fit_function] = background_pdf.createHistogram("total_background_" + fit_function, self.mjj_, RooFit.Binning(roobinning), RooFit.Extended(True))
			background_histograms[fit_function].SetDirectory(0)

			background_histograms[fit_function].SetLineColor(seaborn.GetColorRoot("default", style_counter))
			background_histograms[fit_function].SetLineWidth(1)
			background_histograms[fit_function].SetLineStyle(2)
			background_histograms[fit_function].Draw("hist same")

			background_histograms_fitted_range[fit_function] = background_histograms[fit_function].Clone()
			background_histograms_fitted_range[fit_function].SetDirectory(0)
			low_bin = background_histograms[fit_function].GetXaxis().FindBin(self.mjj_.getMin())
			high_bin = background_histograms[fit_function].GetXaxis().FindBin(self.mjj_.getMax())
			for bin in xrange(1, background_histograms[fit_function].GetNbinsX() + 1):
				if bin < low_bin or bin > high_bin:
					background_histograms_fitted_range[fit_function].SetBinContent(bin, 0)
					background_histograms_fitted_range[fit_function].SetBinError(bin, 0)

			background_histograms_fitted_range[fit_function].SetLineColor(seaborn.GetColorRoot("dark", style_counter))
			background_histograms_fitted_range[fit_function].SetLineWidth(2)
			background_histograms_fitted_range[fit_function].SetLineStyle(1)
			background_histograms_fitted_range[fit_function].Draw("hist same")


			l.AddEntry(background_histograms[fit_function], "Background " + fit_function, "l")

			f_workspace.Close()
			style_counter += 1

		#background_pdf = w_background.pdf("background_pdf")
		#background_hist_raw = background_pdf.createHistogram("mjj", 1000)
		#background_hist_raw.Scale(data_integral / background_hist_raw.Integral())
		#if data_binning:
		#	background_hist = background_hist_raw.Rebin(len(data_binning) - 1, "background_hist_rebinned", data_binning)
		#else:
		#	background_hist = background_hist_raw
		#for bin in xrange(1, data_hist.GetNbinsX() + 1):
		#	background_hist.SetBinContent(bin, background_hist.GetBinContent(bin) / background_hist.GetXaxis().GetBinWidth(bin) * normalization_bin_width)
		#	background_hist.SetBinError(bin, background_hist.GetBinError(bin) / background_hist.GetXaxis().GetBinWidth(bin) * normalization_bin_width)

		#background_hist.SetLineColor(seaborn.GetColorRoot("default", 2))
		#background_hist.SetLineStyle(1)
		#background_hist.SetLineWidth(2)
		#background_hist.Draw("hist same")
		#background_pdf.plotOn(frame_top, RooFit.Name("B Fit"), RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(1), RooFit.LineWidth(2))

		#l.AddEntry(background_hist, "B Fit", "l")

		l.Draw()
		
		# Pull histogram
		c.cd()
		bottom.cd()
		pull_histograms = {}
		pull_histograms_fitted_range = {}
		for fit_function in fit_functions:
			pull_histograms[fit_function] = data_hist.Clone()
			pull_histograms[fit_function].Reset()
			for bin in xrange(1, data_hist.GetNbinsX() + 1):
				if data_hist.GetBinError(bin) > 0:
					pull = (data_hist.GetBinContent(bin) - background_histograms[fit_function].GetBinContent(bin)) / (data_hist.GetBinError(bin))
				else:
					pull = 0.
				#print "[debug] Pull = " + str(pull)
				pull_histograms[fit_function].SetBinContent(bin, pull)
				pull_histograms[fit_function].SetBinError(bin, 0.)

			pull_histograms_fitted_range[fit_function] = pull_histograms[fit_function].Clone()
			for bin in xrange(1, pull_histograms_fitted_range[fit_function].GetNbinsX() + 1):
				bin_center = pull_histograms_fitted_range[fit_function].GetXaxis().GetBinCenter(bin)
				if bin_center < self.mjj_.getMin() or bin_center > self.mjj_.getMax():
					pull_histograms_fitted_range[fit_function].SetBinContent(bin, 0)
					pull_histograms_fitted_range[fit_function].SetBinError(bin, 0)

			#pull_histogram = frame_top.pullHist("Data", "B Fit")
		frame_bottom = TH1D("frame_bottom", "frame_bottom", 100, x_min, x_max)
		frame_bottom.SetMinimum(-5.)
		frame_bottom.SetMaximum(5.)
		#pull_histogram.plotOn(frame_bottom, RooFit.Name(fit_pdf_name))
		frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
		frame_bottom.GetYaxis().SetTitle("#frac{Data - Fit}{#sigma(Data)}")
		frame_bottom.Draw()

		style_counter = 0
		for fit_function in fit_functions:
			pull_histograms[fit_function].SetLineColor(seaborn.GetColorRoot("default", style_counter))
			pull_histograms[fit_function].SetLineWidth(2)
			pull_histograms[fit_function].SetLineStyle(2)
			pull_histograms[fit_function].Draw("hist same")

			pull_histograms_fitted_range[fit_function].SetLineColor(seaborn.GetColorRoot("dark", style_counter))
			pull_histograms_fitted_range[fit_function].SetLineWidth(2)
			pull_histograms_fitted_range[fit_function].SetLineStyle(1)
			pull_histograms_fitted_range[fit_function].Draw("hist same")
			style_counter += 1
		#pull_histogram.GetXaxis().SetTitle("m_{jj} [GeV]")
		#pull_histogram.GetYaxis().SetTitle("#frac{Data - Fit}{#sigma(Fit)}")
		#pull_histogram.Draw("same")
		c.cd()

		c.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Results/figures/c_" + save_tag + ".pdf")

	# fitted_signal_shapes = list of fitted S(+B) shapes to plot.
	# expected_signal_shapes = list of S shapes to plot, scaled to cross sections taken from configuration file
	def plot(self, save_tag, fit_functions, background_workspaces, fitted_signal_workspaces=None, expected_signal_workspaces=None, log=False, x_range=None, data_binning=None, normalization_bin_width=1):
		print "Plotting " + save_tag
		c = TCanvas("c_" + save_tag, "c_" + save_tag, 800, 1200)
		l = TLegend(0.55, 0.6, 0.88, 0.88)
		l.SetFillColor(0)
		l.SetBorderSize(0)
		top = TPad("top", "top", 0., 0.5, 1., 1.)
		top.SetBottomMargin(0.03)
		top.Draw()
		if log:
			top.SetLogy()
		c.cd()
		bottom = TPad("bottom", "bottom", 0., 0., 1., 0.5)
		bottom.SetTopMargin(0.02)
		bottom.SetBottomMargin(0.2)
		bottom.Draw()
		ROOT.SetOwnership(c, False)
		ROOT.SetOwnership(top, False)
		ROOT.SetOwnership(bottom, False)
		top.cd()

		# Make data histogram
		data_hist_raw = self.data_histogram_.Clone()

		if data_binning:
			data_hist = data_hist_raw.Rebin(len(data_binning) - 1, "data_hist_rebinned", data_binning)
			for bin in xrange(1, data_hist.GetNbinsX() + 1):
				data_hist.SetBinContent(bin, data_hist.GetBinContent(bin) / data_hist.GetXaxis().GetBinWidth(bin) * normalization_bin_width)
				data_hist.SetBinError(bin, data_hist.GetBinError(bin) / data_hist.GetXaxis().GetBinWidth(bin) * normalization_bin_width)
		else:
			data_hist = data_hist_raw

		if x_range:
			x_min = x_range[0]
			x_max = x_range[1]
		else:
			x_min = self.mjj_.GetMin()
			x_max = self.mjj_.GetMax()
		frame_top = TH1D("frame_top", "frame_top", 100, x_min, x_max)
		if log:
			frame_top.SetMaximum(data_hist.GetMaximum() * 50.)
			frame_top.SetMinimum(0.1)
		else:
			frame_top.SetMaximum(data_hist.GetMaximum() * 1.3)
			frame_top.SetMinimum(0.)
		frame_top.GetYaxis().SetTitle("Events / " + str(normalization_bin_width) + " GeV")
		frame_top.GetXaxis().SetTitleSize(0)
		frame_top.GetXaxis().SetLabelSize(0)
		#if log:
		#	frame_top.SetMaximum(y_max * 10.)
		#	frame_top.SetMinimum(y_min / 100.)
		#else:
		#	frame_top.SetMaximum(y_max * 1.3)
		#	frame_top.SetMinimum(0.)
		frame_top.Draw()
		#data_rdh.plotOn(frame_top, RooFit.Name("Data"))
		data_hist.SetMarkerStyle(20)
		data_hist.SetMarkerColor(1)
		data_hist.SetMarkerSize(1)
		data_hist.Draw("p same")
		l.AddEntry(data_hist, "Data", "pl")

		#if fitted_signal_workspaces:
		#	for fitted_signal_workspace in fitted_signal_workspaces:
		#		# Load from workspace
		#		print "Opening signal workspace from " + fitted_signal_workspace
		#		f_in = TFile(fitted_signal_workspace, "READ")
		#		w = f_in.Get("w")
		#		w.Print()
		#		fit_pdf = w.pdf("signal_pdf")
		#		fit_pdf_name = "S+B Fit"
		#		fit_pdf.SetName(fit_pdf_name)
		#		fit_pdf.plotOn(frame_top, RooFit.Name(fit_pdf_name), RooFit.LineColor(seaborn.GetColorRoot("default", style_counter)), RooFit.LineStyle(1), RooFit.LineWidth(2))
		#		l.AddEntry(frame_top.findObject(fit_pdf_name), fit_pdf_name, "l")
		#		style_counter += 1
		#		f_in.Close()
		#if expected_signal_workspaces:
		#	for expected_signal_workspace in expected_signal_workspaces:
		#		f_in = TFile(expected_signal_workspace, "READ")
		#		w = f_in.Get("w")
		#		fit_pdf = w.pdf("model")

	#	#		self.signal_roohistograms_[signal_name].plotOn(frame_top, RooFit.Name(signal_name), RooFit.Rescale(cross_section * self.luminosity_ / self.signal_roohistograms_[signal_name].sum()), RooFit.LineColor(seaborn.GetColorRoot("pastel", style_counter)), RooFit.LineStyle(2), RooFit.LineWidth(2))
		#		l.AddEntry(frame_top.findObject(signal_name), signal_name, "l")
		#		style_counter += 1

		background_tf1s = {}
		background_tf1s_fitted_range = {}

		background_histograms = {}
		style_counter = 0
		for fit_function in fit_functions:
			f_workspace = TFile(background_workspaces[fit_function], "READ")
			workspace = f_workspace.Get("w")
			fitresult_name = "fitresult_model_" + fit_function + "_data_roohistogram"
			fitresult = workspace.genobj(fitresult_name)
			fitresult.Print()
			background_tf1s[fit_function] = self.make_background_tf1_from_roofitresult(fit_function, fitresult, mjj_range=x_range)
			# Force normalization?
			print "[MjjFit::plot] INFO : Normalizing fit to data over range [" + str(self.mjj_.getMin()) + ", " + str(self.mjj_.getMax()) + "]"
			data_integral = data_hist.Integral(data_hist.GetXaxis().FindBin(self.mjj_.getMin()), data_hist.GetXaxis().FindBin(self.mjj_.getMax()), "width")
			integration_xmin = data_hist.GetXaxis().GetBinLowEdge(data_hist.GetXaxis().FindBin(self.mjj_.getMin()))
			integration_xmax = data_hist.GetXaxis().GetBinUpEdge(data_hist.GetXaxis().FindBin(self.mjj_.getMax()))
			scale_factor = data_integral / background_tf1s[fit_function].Integral(integration_xmin, integration_xmax)
			background_tf1s[fit_function].SetParameter(0, background_tf1s[fit_function].GetParameter(0) * scale_factor)
			background_tf1s[fit_function].SetParError(0, background_tf1s[fit_function].GetParError(0) * scale_factor)
			background_tf1s[fit_function].SetLineColor(seaborn.GetColorRoot("default", style_counter))
			background_tf1s[fit_function].SetLineWidth(1)
			background_tf1s[fit_function].SetLineStyle(2)
			background_tf1s[fit_function].Draw("l same")

			background_tf1s_fitted_range[fit_function] = background_tf1s[fit_function].Clone()
			background_tf1s_fitted_range[fit_function].SetRange(self.mjj_.getMin(), self.mjj_.getMax())
			background_tf1s_fitted_range[fit_function].SetLineColor(seaborn.GetColorRoot("dark", style_counter))
			background_tf1s_fitted_range[fit_function].SetLineWidth(2)
			background_tf1s_fitted_range[fit_function].SetLineStyle(1)
			background_tf1s_fitted_range[fit_function].Draw("l same")


			l.AddEntry(background_tf1s[fit_function], "Background " + fit_function, "l")

			f_workspace.Close()
			style_counter += 1

		#background_pdf = w_background.pdf("background_pdf")
		#background_hist_raw = background_pdf.createHistogram("mjj", 1000)
		#background_hist_raw.Scale(data_integral / background_hist_raw.Integral())
		#if data_binning:
		#	background_hist = background_hist_raw.Rebin(len(data_binning) - 1, "background_hist_rebinned", data_binning)
		#else:
		#	background_hist = background_hist_raw
		#for bin in xrange(1, data_hist.GetNbinsX() + 1):
		#	background_hist.SetBinContent(bin, background_hist.GetBinContent(bin) / background_hist.GetXaxis().GetBinWidth(bin) * normalization_bin_width)
		#	background_hist.SetBinError(bin, background_hist.GetBinError(bin) / background_hist.GetXaxis().GetBinWidth(bin) * normalization_bin_width)

		#background_hist.SetLineColor(seaborn.GetColorRoot("default", 2))
		#background_hist.SetLineStyle(1)
		#background_hist.SetLineWidth(2)
		#background_hist.Draw("hist same")
		#background_pdf.plotOn(frame_top, RooFit.Name("B Fit"), RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(1), RooFit.LineWidth(2))

		#l.AddEntry(background_hist, "B Fit", "l")

		l.Draw()
		
		# Pull histogram
		c.cd()
		bottom.cd()
		pull_histograms = {}
		pull_histograms_fitted_range = {}
		for fit_function in fit_functions:
			pull_histograms[fit_function] = data_hist.Clone()
			pull_histograms[fit_function].Reset()
			for bin in xrange(1, data_hist.GetNbinsX() + 1):
				if data_hist.GetBinError(bin) > 0:
					background_integral = background_tf1s[fit_function].Integral(data_hist.GetXaxis().GetBinLowEdge(bin), data_hist.GetXaxis().GetBinUpEdge(bin))
					bin_width = data_hist.GetXaxis().GetBinWidth(bin)
					pull = (data_hist.GetBinContent(bin)*bin_width - background_integral) / (data_hist.GetBinError(bin)*bin_width)
				else:
					pull = 0.
				#print "[debug] Pull = " + str(pull)
				pull_histograms[fit_function].SetBinContent(bin, pull)
				pull_histograms[fit_function].SetBinError(bin, 0.)

			pull_histograms_fitted_range[fit_function] = pull_histograms[fit_function].Clone()
			for bin in xrange(1, pull_histograms_fitted_range[fit_function].GetNbinsX() + 1):
				bin_center = pull_histograms_fitted_range[fit_function].GetXaxis().GetBinCenter(bin)
				if bin_center < self.mjj_.getMin() or bin_center > self.mjj_.getMax():
					pull_histograms_fitted_range[fit_function].SetBinContent(bin, 0)
					pull_histograms_fitted_range[fit_function].SetBinError(bin, 0)

			#pull_histogram = frame_top.pullHist("Data", "B Fit")
		frame_bottom = TH1D("frame_bottom", "frame_bottom", 100, x_min, x_max)
		frame_bottom.SetMinimum(-5.)
		frame_bottom.SetMaximum(5.)
		#pull_histogram.plotOn(frame_bottom, RooFit.Name(fit_pdf_name))
		frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
		frame_bottom.GetYaxis().SetTitle("#frac{Data - Fit}{#sigma(Data)}")
		frame_bottom.Draw()

		style_counter = 0
		for fit_function in fit_functions:
			pull_histograms[fit_function].SetLineColor(seaborn.GetColorRoot("default", style_counter))
			pull_histograms[fit_function].SetLineWidth(2)
			pull_histograms[fit_function].SetLineStyle(2)
			pull_histograms[fit_function].Draw("hist same")

			pull_histograms_fitted_range[fit_function].SetLineColor(seaborn.GetColorRoot("dark", style_counter))
			pull_histograms_fitted_range[fit_function].SetLineWidth(2)
			pull_histograms_fitted_range[fit_function].SetLineStyle(1)
			pull_histograms_fitted_range[fit_function].Draw("hist same")
			style_counter += 1
		#pull_histogram.GetXaxis().SetTitle("m_{jj} [GeV]")
		#pull_histogram.GetYaxis().SetTitle("#frac{Data - Fit}{#sigma(Fit)}")
		#pull_histogram.Draw("same")
		c.cd()

		c.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Results/figures/c_" + save_tag + ".pdf")


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Run and plot fits")
	parser.add_argument("analysis_name", type=str, help='Analysis name (see analysis_configuration_8TeV.py)')
	parser.add_argument("data_sample", type=str, default="BJetPlusX_2012", help='Data sample name (see analysis_configuration_8TeV.py)')
	parser.add_argument("--fit", action="store_true", help="Run mjj fits. Background fit is always run; signal fit is run if --signal is specified.")
	parser.add_argument("--signal", type=str, help='Signal model(s) to fit. Use comma-separated list to specify multiple.')
	parser.add_argument("--fixed_signal", type=str, help='Signal model(s) to plot, normalized to expected cross sections')
	parser.add_argument("--plot", action="store_true", help="Plot mjj spectra and fits. Background fit is always plotted; signal fits are plotted if --signal is specified.")
	parser.add_argument("--x_range", type=int, nargs=2, help="Plot xrange")
	parser.add_argument('--backgrounds', type=str, help='Backgrounds to subtract (comma-separated list)')

	# Fit options
	parser.add_argument("-l", "--lumi", dest="lumi",
						default=19700., type=float,
						help="Integrated luminosity in pb-1 (default: %(default).1f)",
						metavar="LUMI")

	parser.add_argument("--min_mjj", dest="min_mjj",
						default=0., type=int,
						help="Min m(jj) considered (default: %(default)s)",
						metavar="MASS_MIN")

	parser.add_argument("--max_mjj", dest="max_mjj",
						default=2000., type=int,
						help="Max m(jj) considered (default: %(default)s)",
						metavar="MASS_MAX")

	parser.add_argument("--fit_range", dest="fit_range",
						type=int, nargs=2,
						help="m(jj) range used in fit (default: %(default)s)",
						metavar="MASS_MIN")

	#parser.add_argument("--fit_max_mjj", dest="fit_max_mjj",
	#					default=1500., type=float,
	#					help="Upper bound of the mass range used for fitting (default: %(default)s)",
	#					metavar="MASS_MAX")

	args = parser.parse_args()

	# Pack up fit options
	fit_options = {}
	fit_options["min_mjj"] = args.min_mjj
	fit_options["max_mjj"] = args.max_mjj
	#fit_options["fit_range"] = [args.fit_min_mjj, args.fit_max_mjj]
	mjj_fit = MjjFit([args.min_mjj, args.max_mjj])
	if args.fit_range:
		mjj_fit.set_fit_range(args.fit_range)
	data_file = TFile(analysis_config.get_b_histogram_filename(args.analysis_name, args.data_sample), "READ")
	mjj_fit.add_data(data_file.Get("BHistograms/h_pfjet_mjj"), args.lumi)
	data_file.Close()

	if args.backgrounds:
		backgrounds = args.backgrounds.split(",")
		for background in backgrounds:
			if analysis_config.simulation.background_supersamples.has_key(background):
				first = True
				for subbackground in analysis_config.simulation.background_supersamples[background]:
					background_file = TFile(analysis_config.get_b_histogram_filename(args.analysis_name, subbackground), "READ")
					input_nevents = background_file.Get("BHistograms/h_input_nevents").Integral()
					if first:
						background_histogram = background_file.Get("BHistograms/h_pfjet_mjj").Clone()
						background_histogram.Scale(args.lumi * analysis_config.simulation.background_cross_sections[subbackground] / input_nevents)
						background_histogram.SetDirectory(0)
					else:
						background_histogram.Add(background_file.Get("BHistograms/h_pfjet_mjj").Scale(args.lumi * analysis_config.simulation.background_cross_sections[subbackground] / input_nevents))
					background_file.Close()
			else:
				background_file = TFile(analysis_config.get_b_histogram_filename(args.analysis_name, background), "READ")
				input_nevents = background_file.Get("BHistograms/h_input_nevents").Integral()
				background_histogram = background_file.Get("BHistograms/h_pfjet_mjj")
				background_histogram.Scale(args.lumi * analysis_config.simulation.background_cross_sections[background] / input_nevents)
			mjj_fit.add_background(background_histogram, background, analysis_config.simulation.background_cross_section_uncertainties[background])

	signal_models = []
	if args.signal:
		signal_models = args.signal.split(",")
		for signal_model in signal_models:
			signal_file = TFile(analysis_config.get_b_histogram_filename(args.analysis_name, args.data_sample), "READ")
			mjj_fit.add_signal(signal_model, signal_file.Get("BHistograms/h_pfjet_mjj"))
			signal_file.Close()

	# For trigbbl, apply a trigger correction. This has to be done after loading all the histograms.
	if "trigbbl" in args.analysis_name:
		mjj_fit.correct_trigger()

	if args.fit:
		#for fit_function in ["f1", "f2", "f3", "f4", "f5"]:
		for fit_function in ["f1"]:
			mjj_fit.fit(fit_function, limit_paths.get_workspace_filename(args.analysis_name, "background_" + fit_function), fit_options)
			for signal_model in signal_models:
				mjj_fit.fit(fit_function, limit_paths.get_workspace_filename(args.analysis_name, signal_model), fit_options, signal_name=signal_model)

	mass_bins = array("d", [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8000])

	if args.plot:
		print "Plotting"
		fitted_signal_workspaces = []
		expected_signal_workspaces = []
		if args.signal:
			for signal_model in args.signal.split(","):
				fitted_signal_workspaces.append(limit_paths.get_workspace_filename(args.analysis_name, signal_model))
		if args.fixed_signal:
			for signal_model in args.fixed_signal.split(","):
				expected_signal_workspaces.append(limit_paths.get_workspace_filename(args.analysis_name, signal_model))

		fit_functions = ["f1", "f2", "f3", "f4", "f5"]
		background_workspaces = {}
		for fit_function in fit_functions:
			background_workspaces[fit_function] = limit_paths.get_workspace_filename(args.analysis_name, "background_" + fit_function)
		if args.x_range:
			x_range = args.x_range
		else:
			x_range = [0., 2000.]
		mjj_fit.rooplot("mjj_fits_" + args.analysis_name, fit_functions, background_workspaces, fitted_signal_workspaces=fitted_signal_workspaces, expected_signal_workspaces=expected_signal_workspaces, log=True, x_range=x_range, data_binning=mass_bins, normalization_bin_width=1.)

