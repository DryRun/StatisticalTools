#!/usr/bin/env python
import os
import sys
import re
import math
import time
from argparse import ArgumentParser
from array import array
import numpy as np
import CMS_lumi
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/src/CMSDIJET/QCDAnalysis/python/")
import analysis_configuration_8TeV as analysis_config
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_7_4_15/python/MyTools/RootUtils")
import histogram_tools

import ROOT
from ROOT import *
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
seaborn = Root.SeabornInterface()
seaborn.Initialize()
Root.SetCanvasStyle()

# Get sigmatilde
import pickle
sigmatilde_file = open(os.path.expandvars("$CMSSW_BASE/src/CMSDIJET/StatisticalTools/data/sigmatilde2.pkl"), 'r')
sigmatilde = pickle.load(sigmatilde_file)
print sigmatilde
sigmatilde_file.close()

# Analytical function for zeta
# Equal coupling, 
def get_zeta(x):
	gB = 0.25 * 6
	tpart = (1 + 2*172.**2/x[0]**2) * math.sqrt(max(0., 1.-4.*172.**2/x[0]**2))
	BR_initial = 2. / (5. + tpart)
	BR_final = 1. / (5. + tpart)
	width = (gB**2 * x[0] * (5. + tpart))/(144.*math.pi)
	zeta = BR_initial * BR_final * width / x[0]
	print "zeta({}) = {} * {} * {} / {} = {}".format(x[0], BR_initial, BR_final, width, x[0], zeta)
	return zeta

class ZetaPlot():
	def __init__(self):
		self._zeta_graphs = {}
		self._zeta_tf1s = {}
		self._legend_entries = {}
		self._theory_legend_entries = {}
		self._shadings = {}
		self._xrange = [1.e20, -1.e20]
		self._yrange = [1.e20, -1.e20]
		self._names = []
		self._theory_names = []

	# initial_states = u, d, or ud
	# final_states = b
	def add_xsbr_graph(self, name, legend, xsbr_limit, initial_states, marker_style=None, marker_size=None, marker_color=None, line_color=None, line_style=None, line_width=None):
		print "add_xsbr_graph({})".format(name)
		self._names.append(name)
		self._zeta_graphs[name] = TGraph(xsbr_limit.GetN())
		for i in xrange(xsbr_limit.GetN()):
			mass = xsbr_limit.GetX()[i]
			xsbr = xsbr_limit.GetY()[i]
			zeta = xsbr / sigmatilde[mass][initial_states]
			print "zeta = {} / {} = {}".format(xsbr, sigmatilde[mass][initial_states], xsbr / sigmatilde[mass][initial_states])
			self._zeta_graphs[name].SetPoint(i, mass, zeta)

			self._legend_entries[name] = legend

			if marker_style:
				self._zeta_graphs[name].SetMarkerStyle(marker_style)
			if marker_size:
				self._zeta_graphs[name].SetMarkerSize(marker_size)
			self._zeta_graphs[name].SetMarkerSize(0)
			if marker_color:
				self._zeta_graphs[name].SetMarkerColor(marker_color)
			if line_color:
				self._zeta_graphs[name].SetLineColor(line_color)
			if line_style:
				self._zeta_graphs[name].SetLineStyle(line_style)

			# Keep track of plot boundaries
			if mass < self._xrange[0]:
				self._xrange[0] = mass
			if mass > self._xrange[1]:
				self._xrange[1] = mass
			if zeta < self._yrange[0]:
				self._yrange[0] = zeta
			if zeta > self._yrange[0]:
				self._yrange[1] = zeta

	def add_zeta_tf1(self, name, legend, zeta_tf1, line_style=None, line_width=None, line_color=None):
		print "add_zeta_tf1({})".format(name)
		self._theory_names.append(name)
		self._zeta_tf1s[name] = zeta_tf1
		if line_style:
			self._zeta_tf1s[name].SetLineStyle(line_style)
		if line_width:
			self._zeta_tf1s[name].SetLineWidth(line_width)
		if line_color:
			self._zeta_tf1s[name].SetLineColor(line_color)
		self._theory_legend_entries[name] = legend

	def add_shading(self, name1, name2, line_color=None, fill_color=kGray, fill_style=3001):
		print "add_shading({}, {})".format(name1, name2)
		if not name1 in self._zeta_graphs:
			print "ERROR : In add_shading, {} doesn't exist in self._zeta_graphs".format(name1)
			sys.exit(1)
		if not name2 in self._zeta_graphs:
			print "ERROR : In add_shading, {} doesn't exist in self._zeta_graphs".format(name2)
			sys.exit(1)
		self._shadings[name1 + name2] = TGraph(self._zeta_graphs[name1].GetN() + self._zeta_graphs[name2].GetN())
		for i in xrange(self._zeta_graphs[name1].GetN()):
			self._shadings[name1 + name2].SetPoint(i, self._zeta_graphs[name1].GetX()[i], self._zeta_graphs[name1].GetY()[i])
		for i in xrange(self._zeta_graphs[name2].GetN()):
			j = self._zeta_graphs[name2].GetN() - i - 1
			self._shadings[name1 + name2].SetPoint(self._zeta_graphs[name1].GetN() + i, self._zeta_graphs[name2].GetX()[j], self._zeta_graphs[name2].GetY()[j])
		self._shadings[name1 + name2].SetFillColor(fill_color)
		if line_color:
			self._shadings[name1 + name2].SetLineColor(line_color)
		self._shadings[name1 + name2].SetFillStyle(fill_style)
		print "Done with add_shading"

	def draw(self, canvas_name=None, logy=True, xtitle="m_{Z'} [GeV]", ytitle="#zeta", yrange=None):
		print "draw()"
		print "x range = ",
		print self._xrange
		print "y range = ",
		print self._yrange
		if not canvas_name:
			canvas_name = str(time.time())
		self._canvas = TCanvas(canvas_name, canvas_name, 700, 500)
		self._legend = TLegend(0.22, 0.6, 0.44, 0.85)
		self._legend.SetFillColor(0)
		self._legend.SetBorderSize(0)
		self._legend.SetHeader("95% CL upper limits")
		if logy:
			self._canvas.SetLogy()
		self._frame = TH1D("frame", "frame", 100, self._xrange[0] - 100, self._xrange[1] + 100)
		if yrange:
			self._frame.SetMinimum(yrange[0])
			self._frame.SetMaximum(yrange[1])
		else:
			if logy:
				self._frame.SetMinimum(self._yrange[0] / 2.)
				self._frame.SetMaximum(self._yrange[1] * 5.)
			else:
				self._frame.SetMinimum(self._yrange[0] * 0.8)
				self._frame.SetMaximum(self._yrange[1] * 1.2)
		self._frame.GetXaxis().SetTitle(xtitle)
		self._frame.GetYaxis().SetTitle(ytitle)
		self._frame.Draw()

		for name1name2, shading in self._shadings.iteritems():
			shading.Draw("F")

		for name in self._names:
			if name in self._zeta_graphs:
				self._zeta_graphs[name].Draw("pl")
			if name in self._legend_entries:
				if name in self._zeta_graphs:
					self._legend.AddEntry(self._zeta_graphs[name], self._legend_entries[name], "pl")
		self._legend.Draw()

		self._theory_legend = TLegend(0.65, 0.27, 0.91, 0.32)
		self._theory_legend.SetBorderSize(0)
		self._theory_legend.SetFillColor(0)
		for name in self._theory_names:
			self._zeta_tf1s[name].Draw("l same")
			self._theory_legend.AddEntry(self._zeta_tf1s[name], self._theory_legend_entries[name], "l")
		self._theory_legend.Draw()
		Root.CMSLabel(0.22, 0.87, "", 1, 0.5); 
		Root.myText(0.7, 0.96, 1, "19.7 fb^{-1} (8 TeV)", 0.5)

	def save(self, save_path):
		print "save()"
		self._canvas.SaveAs(save_path)

def JoinTGraphs(graph1, graph2):
	graph = TGraph(graph1.GetN() + graph2.GetN())
	for i in xrange(graph1.GetN()):
		graph.SetPoint(i, graph1.GetX()[i], graph1.GetY()[i])
	for i in xrange(graph2.GetN()):
		j = i + graph1.GetN()
		graph.SetPoint(j, graph2.GetX()[i], graph2.GetY()[i])
	return graph

if __name__ == "__main__":
	zeta_plot = ZetaPlot()

	f_Xbb_8TeVl = TFile("~/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_trigbbl_CSVTM_ZPrime_dijet4.root", "READ")
	f_Xbb_8TeVh = TFile("~/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_trigbbh_CSVTM_ZPrime_dijet4.root", "READ")

	gr_Xbb_8TeVl_exp = f_Xbb_8TeVl.Get("graph_exp")
	gr_Xbb_8TeVh_exp = f_Xbb_8TeVh.Get("graph_exp")
	gr_Xbb_8TeV_exp = JoinTGraphs(gr_Xbb_8TeVl_exp, gr_Xbb_8TeVh_exp)
	print "Printing xs*br limits, sigmatilde, and zeta:"
	for i in xrange(gr_Xbb_8TeV_exp.GetN()):
		mass = gr_Xbb_8TeV_exp.GetX()[i]
		print "{}\t:\t{}\t/\t{}\t=\t{}".format(mass , gr_Xbb_8TeV_exp.GetY()[i], sigmatilde[mass]["u"], gr_Xbb_8TeV_exp.GetY()[i] / sigmatilde[mass]["u"])
	zeta_plot.add_xsbr_graph("exp_uu", "Exp. u#bar{u}#rightarrowZ'#rightarrowb#bar{b}", gr_Xbb_8TeV_exp, "u", marker_style=25, marker_size=0, line_color=seaborn.GetColorRoot("default", 2), line_style=2, line_width=1)
	zeta_plot.add_xsbr_graph("exp_dd", "Exp. d#bar{d}#rightarrowZ'#rightarrowb#bar{b}", gr_Xbb_8TeV_exp, "d", marker_style=25, marker_size=0, line_color=seaborn.GetColorRoot("default", 3), line_style=2, line_width=1)

	gr_Xbb_8TeVl_obs = f_Xbb_8TeVl.Get("graph_obs")
	gr_Xbb_8TeVh_obs = f_Xbb_8TeVh.Get("graph_obs")
	gr_Xbb_8TeV_obs = JoinTGraphs(gr_Xbb_8TeVl_obs, gr_Xbb_8TeVh_obs)
	zeta_plot.add_xsbr_graph("obs_uu", "Obs. u#bar{u}#rightarrowZ'#rightarrowb#bar{b}", gr_Xbb_8TeV_obs, "u", marker_style=24, marker_size=0, line_color=seaborn.GetColorRoot("default", 2), line_style=1, line_width=1)
	zeta_plot.add_xsbr_graph("obs_dd", "Obs. d#bar{d}#rightarrowZ'#rightarrowb#bar{b}", gr_Xbb_8TeV_obs, "d", marker_style=24, marker_size=0, line_color=seaborn.GetColorRoot("default", 3), line_style=1, line_width=1)

	zeta_plot.add_shading("obs_uu", "obs_dd", fill_color=seaborn.GetColorRoot("pastel", 2), fill_style=3004)

	zp_zeta = TF1("zp_zeta", get_zeta, 325., 1200.)
	zp_zeta.SetParameter(0, 0.25 * 6)
	print zp_zeta
	zeta_plot.add_zeta_tf1("zp", "pp#rightarrowZ'#rightarrowb#bar{b}, g_{q}=0.25", zp_zeta, line_style=3, line_color=seaborn.GetColorRoot("default", 1), line_width=2)

	zeta_plot.draw(logy=True, yrange=[2e-4, 2e-2])
	zeta_plot.save("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/zeta/zeta.pdf")
