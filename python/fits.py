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

def BackgroundFit(x, par):
	return par[0] * (1. - (x[0] / 8.e3))**par[1] / ((x[0] / 8.e3)**(par[2] + par[3] * TMath.Log((x[0] / 8.e3))))

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
		self.luminosity_ = 0.
		self.collision_energy_ = 8000.
		self.signal_histograms_ = {}
		self.signal_roohistograms_ = {}
		self.signal_names_ = []
		self.signal_initial_normalizations_ = {}

		# Fit storage
		self.simple_fit_ = None
		self.mjj_ = RooRealVar('mjj','mjj',float(mass_range[0]),float(mass_range[1]))
		#self.workspace_ = None
		self.saved_workspaces_ = {}

	def add_data(self, data_histogram, luminosity):
		print "[MjFit.add_data] INFO : Adding data histogram"
		# Add a data histogram
		self.data_histogram_ = data_histogram.Clone()
		self.data_histogram_.SetDirectory(0)		
		self.data_roohistogram_ = RooDataHist('data_roohistogram','data_roohistogram',RooArgList(self.mjj_),self.data_histogram_)
		self.data_roohistogram_.Print()
		self.luminosity_ = luminosity

	def add_signal(self, signal_name, signal_histogram):
		print "[MjjFit.add_signal] INFO : Adding signal histogram " + signal_name
		# Add a signal histogram. 
		# Scale to sigma=1 (in whatever units the data luminosity is given in), so the 'r' parameter corresponds to the limit on the cross section.
		if self.luminosity_ == 0:
			print "[MjjFit.add_signal] ERROR : Please set luminosity first (MjjFit.set_luminosity(###))."
			sys.exit(1)
		self.signal_names_.append(signal_name)
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
		fit = TF1(data_histogram.GetName() + "_fit_B", BackgroundFit, fit_min, fit_max, 4)
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
		 
	def fit(self, save_to, fit_options, signal_name=None, fix_p3=False, fit_strategy=1):
		# Run a RooFit fit

		# Create background PDF
		p1 = RooRealVar('p1','p1',fit_options["p1"],0.,100.)
		p2 = RooRealVar('p2','p2',fit_options["p2"],0.,60.)
		p3 = RooRealVar('p3','p3',fit_options["p3"],-10.,10.)
		if fix_p3:
			p3.setConstant()
		background_pdf = RooGenericPdf('background_pdf','(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)))'%(self.collision_energy_,self.collision_energy_,self.collision_energy_),RooArgList(self.mjj_,p1,p2,p3))
		background_pdf.Print()
		data_integral = self.data_histogram_.Integral(self.data_histogram_.GetXaxis().FindBin(float(fit_options["min_mjj"])),self.data_histogram_.GetXaxis().FindBin(float(fit_options["max_mjj"])))
		background_norm = RooRealVar('background_norm','background_norm',data_integral,0.,1e+08)
		background_norm.Print()

		# Create signal PDF and fit model
		if signal_name:
			signal_pdf = RooHistPdf('signal_pdf', 'signal_pdf', RooArgSet(self.mjj_), self.signal_roohistograms_[signal_name])
			signal_pdf.Print()
			signal_norm = RooRealVar('signal_norm','signal_norm',0,-1e+05,1e+05)
			signal_norm.Print()
			signal_initial_norm = RooRealVar('signal_initial_norm', 'signal_initial_norm', self.signal_initial_normalizations_[signal_name])
			model = RooAddPdf("model","s+b",RooArgList(background_pdf,signal_pdf),RooArgList(background_norm,signal_norm))
		else:
			model = RooAddPdf("model","b",RooArgList(background_pdf),RooArgList(background_norm))

		# Run fit
		res = model.fitTo(self.data_roohistogram_, RooFit.Save(kTRUE), RooFit.Strategy(fit_strategy))

		# Save to workspace
		self.workspace_ = RooWorkspace('w','workspace')
		#getattr(w,'import')(background,ROOT.RooCmdArg())
		getattr(self.workspace_,'import')(background_pdf,RooFit.Rename("background"))
		getattr(self.workspace_,'import')(background_norm,ROOT.RooCmdArg())
		getattr(self.workspace_,'import')(self.data_roohistogram_,RooFit.Rename("data_obs"))
		getattr(self.workspace_, 'import')(model, RooFit.Rename("model"))
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

	# fitted_signal_shapes = list of fitted S(+B) shapes to plot.
	# expected_signal_shapes = list of S shapes to plot, scaled to cross sections taken from configuration file
	def plot(self, save_tag, background_workspace, fitted_signal_workspaces=None, expected_signal_workspaces=None, log=False, x_range=None, data_binning=None, normalization_bin_width=20):
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

		# Frame from background fit
		print "Opening background workspace from " + background_workspace
		f_background = TFile(background_workspace, "READ")
		w_background = f_background.Get("w")
		w_background.Print()
		data_rdh = w_background.data("data_obs")
		data_hist_raw = data_rdh.createHistogram("mjj", 1000)
		data_integral = data_hist_raw.Integral()
		if data_binning:
			data_hist = data_hist_raw.Rebin(len(data_binning) - 1, "data_hist_rebinned", data_binning)
		else:
			data_hist = data_hist_raw
		for bin in xrange(1, data_hist.GetNbinsX() + 1):
			data_hist.SetBinContent(bin, data_hist.GetBinContent(bin) / data_hist.GetXaxis().GetBinWidth(bin) * normalization_bin_width)
			data_hist.SetBinError(bin, data_hist.GetBinError(bin) / data_hist.GetXaxis().GetBinWidth(bin) * normalization_bin_width)
		for bin in xrange(1, data_hist.GetNbinsX() + 1):
			print "Bin " + str(bin) + " = " + str(data_hist.GetBinContent(bin))
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

		style_counter = 0
		if fitted_signal_workspaces:
			for fitted_signal_workspace in fitted_signal_workspaces:
				# Load from workspace
				print "Opening signal workspace from " + fitted_signal_workspace
				f_in = TFile(fitted_signal_workspace, "READ")
				w = f_in.Get("w")
				w.Print()
				fit_pdf = w.pdf("signal_pdf")
				fit_pdf_name = "S+B Fit"
				fit_pdf.SetName(fit_pdf_name)
				fit_pdf.plotOn(frame_top, RooFit.Name(fit_pdf_name), RooFit.LineColor(seaborn.GetColorRoot("default", style_counter)), RooFit.LineStyle(1), RooFit.LineWidth(2))
				l.AddEntry(frame_top.findObject(fit_pdf_name), fit_pdf_name, "l")
				style_counter += 1
				f_in.Close()
		if expected_signal_workspaces:
			for expected_signal_workspace in expected_signal_workspaces:
				f_in = TFile(expected_signal_workspace, "READ")
				w = f_in.Get("w")
				fit_pdf = w.pdf("model")

				self.signal_roohistograms_[signal_name].plotOn(frame_top, RooFit.Name(signal_name), RooFit.Rescale(cross_section * self.luminosity_ / self.signal_roohistograms_[signal_name].sum()), RooFit.LineColor(seaborn.GetColorRoot("pastel", style_counter)), RooFit.LineStyle(2), RooFit.LineWidth(2))
				l.AddEntry(frame_top.findObject(signal_name), signal_name, "l")
				style_counter += 1
		background_pdf = w_background.pdf("background_pdf")
		background_hist_raw = background_pdf.createHistogram("mjj", 1000)
		background_hist_raw.Scale(data_integral / background_hist_raw.Integral())
		if data_binning:
			background_hist = background_hist_raw.Rebin(len(data_binning) - 1, "background_hist_rebinned", data_binning)
		else:
			background_hist = background_hist_raw
		for bin in xrange(1, data_hist.GetNbinsX() + 1):
			background_hist.SetBinContent(bin, background_hist.GetBinContent(bin) / background_hist.GetXaxis().GetBinWidth(bin) * normalization_bin_width)
			background_hist.SetBinError(bin, background_hist.GetBinError(bin) / background_hist.GetXaxis().GetBinWidth(bin) * normalization_bin_width)

		background_hist.SetLineColor(seaborn.GetColorRoot("default", 2))
		background_hist.SetLineStyle(1)
		background_hist.SetLineWidth(2)
		background_hist.Draw("hist same")
		#background_pdf.plotOn(frame_top, RooFit.Name("B Fit"), RooFit.LineColor(seaborn.GetColorRoot("default", 2)), RooFit.LineStyle(1), RooFit.LineWidth(2))
		l.AddEntry(background_hist, "B Fit", "l")
		l.Draw()
		
		# Pull histogram
		c.cd()
		bottom.cd()
		pull_histogram = data_hist.Clone()
		for bin in xrange(1, data_hist.GetNbinsX() + 1):
			if data_hist.GetBinError(bin) > 0:
				pull = (data_hist.GetBinContent(bin) - background_hist.GetBinContent(bin)) / data_hist.GetBinError(bin)
			else:
				pull = 0.
			pull_histogram.SetBinContent(bin, pull)
			pull_histogram.SetBinError(bin, 0.)
		#pull_histogram = frame_top.pullHist("Data", "B Fit")
		frame_bottom = TH1D("frame_bottom", "frame_bottom", 100, x_min, x_max)
		frame_bottom.SetMinimum(-5.)
		frame_bottom.SetMaximum(5.)
		#pull_histogram.plotOn(frame_bottom, RooFit.Name(fit_pdf_name))
		frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
		frame_bottom.GetYaxis().SetTitle("#frac{Data - Fit}{#sigma(Data)}")
		frame_bottom.Draw()
		pull_histogram.Draw("hist same")
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

	# Fit options
	parser.add_argument("-l", "--lumi", dest="lumi",
						default=19700., type=float,
						help="Integrated luminosity in pb-1 (default: %(default).1f)",
						metavar="LUMI")

	parser.add_argument("--fit_min_mjj", dest="fit_min_mjj",
						default=500, type=int,
						help="Lower bound of the mass range used for fitting (default: %(default)s)",
						metavar="MASS_MIN")

	parser.add_argument("--fit_max_mjj", dest="fit_max_mjj",
						default=1500, type=int,
						help="Upper bound of the mass range used for fitting (default: %(default)s)",
						metavar="MASS_MAX")

	parser.add_argument("--fit_p1", dest="fit_p1",
						default=1.0003e+01, type=float,
						help="Fit function p1 parameter (default: %(default)e)",
						metavar="P1")

	parser.add_argument("--fit_p2", dest="fit_p2",
						default=5.2648e+00, type=float,
						help="Fit function p2 parameter (default: %(default)e)",
						metavar="P2")

	parser.add_argument("--fit_p3", dest="fit_p3",
						default=0.0000e+00, type=float,
						help="Fit function p3 parameter (default: %(default)e)",
						metavar="P3")

	args = parser.parse_args()

	# Pack up fit options
	fit_options = {}
	fit_options["min_mjj"] = args.fit_min_mjj
	fit_options["max_mjj"] = args.fit_max_mjj
	fit_options["p1"]      = args.fit_p1
	fit_options["p2"]      = args.fit_p2
	fit_options["p3"]      = args.fit_p3
	mjj_fit = MjjFit([fit_options["min_mjj"], fit_options["max_mjj"]])
	data_file = TFile(analysis_config.get_b_histogram_filename(args.analysis_name, args.data_sample), "READ")
	mjj_fit.add_data(data_file.Get("BHistograms/h_pfjet_mjj"), args.lumi)
	data_file.Close()

	signal_models = None
	if args.signal:
		signal_models = args.signal.split(",")
		for signal_model in signal_models:
			signal_file = TFile(analysis_config.get_b_histogram_filename(args.analysis_name, args.data_sample), "READ")
			mjj_fit.add_signal(signal_model, signal_file.Get("BHistograms/h_pfjet_mjj"))
			signal_file.Close()

	if args.fit:
		mjj_fit.fit(limit_paths.get_workspace_filename(args.analysis_name, "background"), fit_options)
		if args.signal:
			for signal_model in signal_models:
				mjj_fit.fit(limit_paths.get_workspace_filename(args.analysis_name, signal_model), fit_options, signal_name=signal_model)

	mass_bins = array("d", [500, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1500])


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
		mjj_fit.plot("mjj_fits_" + args.analysis_name, limit_paths.get_workspace_filename(args.analysis_name, "background"), fitted_signal_workspaces=fitted_signal_workspaces, expected_signal_workspaces=expected_signal_workspaces, log=True, x_range=[300., 2000.], data_binning=mass_bins)

