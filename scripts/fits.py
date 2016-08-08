#! /usr/bin/env python

import os
import sys
from ROOT import TFile
from ROOT import RooRealVar
from ROOT import RooDataHist
from ROOT import RooArgList
from ROOT import RooArgSet
from ROOT import RooAddPdf
from ROOT import RooFit
from ROOT import RooGenericPdf
from ROOT import RooWorkspace
from ROOT import RooMsgService
from ROOT import RooHistPdf
from ROOT import PdfDiagonalizer
RooMsgService.instance().setSilentMode(kTRUE)
RooMsgService.instance().setStreamStatus(0,kFALSE)
RooMsgService.instance().setStreamStatus(1,kFALSE)

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python")
import CMSDIJET.QCDAnalysis.analysis_configuration_8TeV as config

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

		# Fit storage
		self.simple_fit_ = None
		self.mjj_ = RooRealVar('mjj','mjj',float(mass_range[0]),float(mass_range[1]))
		self.workspace_ = None

	def add_data(self, data_histogram):
		print "[MjFit.add_data] INFO : Adding data histogram"
		# Add a data histogram
		self.data_histogram_ = data_histogram.Clone()
		self.data_histogram_.SetDirectory(0)		
		self.data_roohistogram_ = RooDataHist('data_roohistogram','data_roohistogram',RooArgList(self.mjj_),self.data_histogram_)
		self.data_roohistogram_.Print()

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
		self.signal_histograms_[signal_name].Scale(1. * self.luminosity_ / self.signal_histograms_[signal_name].Integral())
		self.signal_roohistograms_[signal_name] = RooDataHist(signal_histogram.GetName() + "_rdh", signal_histogram.GetName() + "_rdh", RooArgList(self.mjj_), signal_histograms[signal_name])
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
		 
	def fit(self, save_to, signal_name=None, fix_p3=False, fit_range=[300., 1200.], fit_strategy=1):
		# Run a RooFit fit

		# Create background PDF
		p1 = RooRealVar('p1','p1',args.p1,0.,100.)
		p2 = RooRealVar('p2','p2',args.p2,0.,60.)
		p3 = RooRealVar('p3','p3',args.p3,-10.,10.)
		if args.fix_p3:
			p3.setConstant()
		background_pdf = RooGenericPdf('background_pdf','(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)))'%(self.collision_energy,self.collision_energy,self.collision_energy),RooArgList(self.mjj_,p1,p2,p3))
		background_pdf.Print()
		data_integral = data_histogram.Integral(data_histogram.GetXaxis().FindBin(float(fit_range[0])),data_histogram.GetXaxis().FindBin(float(fit_range[1])))
		background_norm = RooRealVar('background_norm','background_norm',data_integral,0.,1e+08)
		background_norm.Print()

		# Create signal PDF and fit model
		if signal_name:
			signal_pdf = RooHistPdf('signal_pdf', 'signal_pdf', RooArgSet(self.mjj_), self.signal_roohistograms_[signal_name])
			signal_pdf.Print()
			signal_norm = RooRealVar('signal_norm','signal_norm',0,-1e+05,1e+05)
			signal_norm.Print()
			model = RooAddPdf("model","s+b",RooArgList(background_pdf,signal_pdf),RooArgList(background_norm,signal_norm))
		else:
			model = RooAddPdf("model","b",RooArgList(background_pdf),RooArgList(background_norm))

		# Run fit
		res = model.fitTo(data_, RooFit.Save(kTRUE), RooFit.Strategy(fit_strategy))

		# Save to workspace
		self.workspace_ = RooWorkspace('w','workspace')
		#getattr(w,'import')(background,ROOT.RooCmdArg())
		getattr(self.workspace_,'import')(background_pdf,RooFit.Rename("background"))
		getattr(self.workspace_,'import')(background_norm,ROOT.RooCmdArg())
		getattr(self.workspace_,'import')(self.data_roohistogram_,RooFit.Rename("data_obs"))
		getattr(self.workspace_, 'import')(model, RooFit.Rename("model"))
		if signal_name:
			getattr(self.workspace_,'import')(signal_roohistogram,RooFit.Rename("signal"))
			getattr(self.workspace_,'import')(signal_pdf,RooFit.Rename("signal_pdf"))
			getattr(self.workspace_,'import')(signal_norm,ROOT.RooCmdArg())
	
		self.workspace_.Print()
		self.workspace_.writeToFile(save_to)
		if signal_name:
			roofit_results[signal_name] = save_to
		else:
			roofit_results["background"] = save_to

	# fitted_signal_shapes = list of fitted S(+B) shapes to plot.
	# expected_signal_shapes = list of S shapes to plot, scaled to cross sections taken from configuration file
	def plot(self, save_tag, fitted_signal_shapes=None, expected_signal_shapes=None, log=False, x_range=None):
		c = TCanvas("c_" + save_tag, "c_" + save_tag, 800, 1600)
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
		if not roofit_results.has_key("background"):
			print "[MjjFit.plot] ERROR : Please run background-only fit before plotting (e.g. MjjFit.plot(\"bkgd_file.root\"))"
		f_background = TFile(roofit_results["background"], "READ")
		w_background = f_background.Get("w")
		data_rdh = w_background.data("data_obs")
		if x_range:
			x_min = x_range[0]
			x_max = x_range[1]
		else:
			x_min = self.mjj_.GetMin()
			x_max = self.mjj_.GetMax()
		frame_top = self.mjj_.frame(x_min, x_max)
		frame_top.GetYaxis().SetTitle("Events / " + str(int(data_histogram.GetXaxis().GetBinWidth(1))) + " GeV")
		frame_top.GetXaxis().SetTitleSize(0)
		frame_top.GetXaxis().SetLabelSize(0)
		#if log:
		#	frame_top.SetMaximum(y_max * 10.)
		#	frame_top.SetMinimum(y_min / 100.)
		#else:
		#	frame_top.SetMaximum(y_max * 1.3)
		#	frame_top.SetMinimum(0.)
		frame_top.Draw()
		data_rdh.plotOn(frame_top, RooFit.Name("Data"))
		l.AddEntry(frame_top.findObject("Data"), "Data", "pl")
		background_pdf = w_background.pdf("background")
		background_pdf.plotOn(frame_top, RooFit.Name("Background Fit"), RooFit.LineColor(seaborn.GetColorRoot("default", 0)), RooFit.LineStyle(1), RooFit.LineWidth(2))

		style_counter = 1
		if fitted_signal_shapes:
			for signal_name in fitted_signal_shapes:
				# Load from workspace
				if not roofit_results.has_key(signal_name):
					print "[MjjFit.plot] ERROR : Signal name " + signal_name + " does not exist in roofit_results. Did you run the fit first?"
					sys.exit(1)
				f_in = TFile(roofit_results[signal_name], "READ")
				w = f_in.Get("workspace")
				fit_pdf = w.pdf("model")
				fit_pdf_name = "Fit (" + signal_name + ")"
				fit_pdf.SetName(fit_pdf_name)
				fit_pdf.plotOn(frame_top, RooFit.Name(fit_pdf_name), RooFit.LineColor(seaborn.GetColorRoot("default", style_counter)), RooFit.LineStyle(1), RooFit.LineWidth(2))
				l.AddEntry(frame_top.findObject(signal_name), signal_name, "l")
				style_counter += 1
				f_in.Close()
		if expected_signal_shapes:
			for signal_name in expected_signal_shapes:
				self.signal_roohistograms_[signal_name].plotOn(frame_top, RooFit.Name(signal_name), RooFit.Rescale(cross_section * self.luminosity_ / self.signal_roohistograms_[signal_name].sum()), RooFit.LineColor(seaborn.GetColorRoot("pastel", style_counter)), RooFit.LineStyle(2), RooFit.LineWidth(2))
				l.AddEntry(frame_top.findObject(signal_name), signal_name, "l")
				style_counter += 1
		frame_top.Draw()
		l.Draw()
		
		# Pull histogram
		c.cd()
		bottom.cd()
		pull_histogram = frame_top.pullHist("Data", "Background Fit")
		frame_bottom = self.mjj_.frame(x_min, x_max)
		pull_histogram.plotOn(frame_bottom, RooFit.Name(fit_pdf_name))
		frame_bottom.GetXaxis().SetTitle("m_{jj} [GeV]")
		frame_bottom.GetYaxis().SetTitle("#frac{Data - Fit}{#sigma(Fit)}")
		frame_bottom.Draw()
		c.cd()

		c.SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Results/figures/c_" + save_tag + ".pdf")


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Run and plot fits")
	parser.add_argument("analysis_name", type=str, help='Analysis name (see analysis_configuration_8TeV.py)')
	parser.add_argument("--fit", action="store_true", help="Run mjj fits. Background fit is always run; signal fit is run if --signal is specified.")
	parser.add_argument("--signal", type=str, help='Signal model(s) to fit. Use comma-separated list to specify multiple.')
	parser.add_argument("--plot", action="store_true", help="Plot mjj spectra and fits. Background fit is always plotted; signal fits are plotted if --signal is specified.")


usage = "python submit_fits_batch.py -q <queue> -i <input list> -o <output dir> -t <num of toys> --run"

parser = argparse.ArgumentParser(description='Process options.')

parser.add_argument('-q', '--queue', type=str, action='store',     dest='queue',       
		help='run in batch in queue specified as option (default -q cmslong)', 
		default='cmslong',
		metavar="QUEUE"
		)
parser.add_argument("--inputListGen", type=str, dest="inputListGen", default="../lists/list_Qstar_genOnly.txt",
		help="list of datacards for the generation",
		)
parser.add_argument("--inputListFit", type=str, dest="inputListFit", default="../lists/list_Qstar.txt",
		help="list of datacards for the fit with flatParameter option",
		)
parser.add_argument("--inputFileLimits", type=str, dest="inputFileLimits", default="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/dijet_limits_interpolated_pseudodatasetDinko/higgsCombine'+sample+'_limit.Asymptotic.mH120.root",
		help="inputfile containing the limits (obtained with combine -m Asymptotic command)",
		)
parser.add_argument("-o", "--output", type=str, dest="output", default="./",
		help="the directory OUTDIR contains the output of the program",
		metavar="OUTDIR"
		)
parser.add_argument("-t", "--toys", type=int, dest="toys", default=10,
		help="number of toys",
		metavar="TOYS"
		)
parser.add_argument("--mu", action="store", type=float, dest="mu", default=-999.,
		help="MU",
		)
parser.add_argument("--tag", action="store", type=str, dest="tag", default="",
		help="some useful tag",
		)


parser.add_argument("--run", action="store_true", dest="run", default=False,
		help="run in batch",
		)



args = parser.parse_args()
print args
#####################################

os.system("mkdir -p batch")
pwd = os.environ['PWD']

#print(args.queue)
ins = open(args.inputListGen,"r")
insFit = open(args.inputListFit,"r")

############# read datacard for generation #####################
ii=0
k=0
datacards_fit = []

for line in  insFit:
	datacards_fit.append(line)

for line in  ins:
	sample = os.path.basename(line)
	#sample = os.path.splitext(line)[0]
	sample = sample.split("_")[0]
	print ("process %s" % sample)
	line = line.rstrip('\n')
	sample = sample.rstrip('\n')
	mass = sample.split("Qstar")[1]
	print mass

	if not(os.path.isdir("args.output")):
		os.system("mkdir -p "+args.output)
	
	
	### when you want to inject 0 (or "fixed" signal strenght) signal you pass the signal strenght as an argument with --mu option 
	if not(args.mu==-999.):
		command = "combine -M GenerateOnly --expectSignal "+str(args.mu)+"  -n "+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+"  -t "+str(args.toys)+" --saveToys --rMin -1000 --rMax 1000 "+line
		command_fit = "combine -M MaxLikelihoodFit  -n "+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+" -t "+str(args.toys)+" --out "+args.output+" --saveNormalizations  --toysFile "+args.output+"/higgsCombine"+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+".GenerateOnly.mH120.123456.root --rMin -1000 --rMax 1000  --robustFit=0 --verbose=2 "+datacards_fit[ii]
	
	### when you want to inject the 2-sigma signal you take the corresponding signal strenght from the limits output file 
	else:
		#filename = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/dijet_limits_interpolated_pseudodatasetDinko/higgsCombine'+sample+'_limit.Asymptotic.mH120.root'
		filename = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_5/src/StatisticalTools/scripts/output_limits_biasStudy_2fb-1/higgsCombine'+sample+'_limit.Asymptotic.mH120.root'
		inf = TFile.Open(filename)
		tr = inf.Get('limit')
		y = []
		N = tr.GetEntriesFast()
		for i in xrange(N):
			tr.GetEntry(i)
			print tr.limit
			y.append(tr.limit)
			k+=1
	
		if (N>1): 
			twosigmalimit = y[4]
		else: 
			twosigmalimit = 0
		print "mu lim hi band : "+str(twosigmalimit)
		command = "combine -M GenerateOnly --expectSignal "+str(twosigmalimit)+"  -n "+sample+args.tag+"_MLfit_mu_limit  -t "+str(args.toys)+"  --saveToys   --rMin "+str(twosigmalimit-10.*twosigmalimit)+" --rMax "+str(twosigmalimit+10.*twosigmalimit)+"   "+line
		command_fit = "combine -M MaxLikelihoodFit  -n "+sample+args.tag+"_MLfit_mu_limit  -t "+str(args.toys)+" --out "+args.output+" --saveNormalizations --toysFile "+args.output+"/higgsCombine"+sample+args.tag+"_MLfit_mu_limit.GenerateOnly.mH120.123456.root --rMin "+str(twosigmalimit-10.*twosigmalimit)+" --rMax "+str(twosigmalimit+10.*twosigmalimit)+" --robustFit=0 --verbose=2 "+datacards_fit[ii]

	ii += 1  
	
	print command
	print command_fit
	print ""

	if not args.mu==-999:
		logfile = "logfile_"+sample+args.tag+"_mu"+str(int(args.mu))+".log"
		outputname = "batch/submit_"+sample+args.tag+"_mu"+str(int(args.mu))+".src"
	else:
		logfile = "logfile_"+sample+args.tag+"_mu_limit.log"
		outputname = "batch/submit_"+sample+args.tag+"_mu_limit.src"

	outputfile = open(outputname,'w')
	outputfile.write('#!/bin/bash\n')
	outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc481\n')
	outputfile.write('cd '+pwd+' \n')
	outputfile.write('eval `scramv1 runtime -sh`\n')
	outputfile.write(command+"\n")
	if not(args.mu==-999.):
		outputfile.write("mv higgsCombine"+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+".GenerateOnly.mH120.123456.root  "+args.output+"\n")  
		outputfile.write(command_fit+"\n")
		outputfile.write("mv higgsCombine"+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+".MaxLikelihoodFit.mH120.123456.root  "+args.output+"\n")  
		
	else:
		outputfile.write("mv higgsCombine"+sample+args.tag+"_MLfit_mu_limit.GenerateOnly.mH120.123456.root  "+args.output+"/\n")  
		outputfile.write(command_fit+"\n")
		outputfile.write("mv higgsCombine"+sample+args.tag+"_MLfit_mu_limit.MaxLikelihoodFit.mH120.123456.root  "+args.output+"/\n")  

	print outputname 
	if args.run:
		os.system("bsub -q "+args.queue+" -o batch/"+logfile+" source "+pwd+"/"+outputname)

