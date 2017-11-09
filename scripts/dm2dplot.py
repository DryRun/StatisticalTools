#!/usr/bin/env python
import os
import sys
import re
import math
from argparse import ArgumentParser
from array import array
import numpy as np
import CMS_lumi
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/src/CMSDIJET/QCDAnalysis/python/")
import analysis_configuration_8TeV as analysis_config

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
gStyle.SetPalette(1)

#import MyTools.RootUtils.tdrstyle as tdrstyle
#tdrstyle.setTDRStyle()

# Computed cross section times BR(jj) of dark matter model. These are the basis for the whole conversion.
# 13 TeV values
#dm_x = array('d', [528.7733914858716,631.2905390292997,795.7485275479014,980.9839334973528,1166.2193394468052,1331.0361216730032,1578.3102214269736,1825.7147916200695,2073.0541265936026,2217.3544322672033,2454.5170729888905,2701.98687840155,2928.9075896518298,3104.1946246179077,3310.4357712666806,3475.513494371131,3568.343211809438])
#dm_y = array('d', [57.269239965262344,26.146418729332762,10.943588903119819,4.998039360562294,2.2826513012206147,1.3151916426227461,0.6008150098368491,0.30829460763931776,0.1492637583468037,0.09943936422165509,0.0556698333169542,0.03027483510044679,0.018491825891702274,0.011967992892333293,0.007309405483957529,0.005313471607526083,0.004337281648532702])

# 8 TeV values from Bogdan, gB=1
# gq = 1./6.
#dm_x = array('d', [400, 500, 700, 1000, 1300, 1500, 2000, 2300, 2400, 3000, 4000, 5000])
#dm_y = array('d', [122.59814203237992,50.135970631685996,12.276050098847596,2.38572459898666,0.616459792663593,0.2752307540653187,0.043874144927211766,0.015493830332160183,0.011048705099045492,0.001478458283748929,0.00005101771602514974,0.0000019151463694130615,])
#dm_sigma_times_BRjj = TGraph(len(dm_x), dm_x, dm_y)

dm_x = array('d', [325,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200])
dm_y = array('d', [631.01208, 489.916584, 305.702112, 197.110728, 131.071932, 89.81076, 62.855808, 45.185256, 32.940768, 24.4725504, 18.4185048, 14.04279, 10.8025608, 8.3939928, 6.58188, 5.20744332, 4.13820204, 3.31039356, 2.66449692])
dm_sigma = TGraph(len(dm_x), dm_x, dm_y)
gq=0.25

from limit_combination import *

# Convert xs (total Z' production) to gq limit
def xs2gq_vector(xs, mZp, gdm, mdm):
	nodm_width = total_width_vector(gq, 0., mZp, 1.e20)
	yesdm_width = total_width_vector(gq, gdm, mZp, mdm)
	this_gq = gq * math.sqrt(xs / dm_sigma.Eval(mZp)) * math.sqrt(0.5 + 0.5*math.sqrt(1. + 4. * (yesdm_width - nodm_width) / nodm_width))
	return this_gq

def make_mdm_mZp_graph(xs_graph, gdm):
	mdm_array = array('d', [])
	mZp_array = array('d', [])
	gq_array  = array('d', [])
	
	max_mdm = 1500.
	n_mdm = 150
	for i_mdm in xrange(150):
		this_mdm = max_mdm / n_mdm * i_mdm
		for i_mZp in xrange(xs_graph.GetN()):
			this_mZp = xs_graph.GetX()[i_mZp]
			this_xs = xs_graph.GetY()[i_mZp]
			this_gq = xs2gq_vector(this_xs, this_mZp, gdm, this_mdm)
			
			mdm_array.append(this_mdm)
			mZp_array.append(this_mZp)
			gq_array.append(this_gq)
	return TGraph2D(len(mdm_array), mZp_array, mdm_array, gq_array)

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description = 'Dijet mass spectrum fits')
	parser.add_argument('--contour', type=float, help='Draw contours at specified gq')
	parser.add_argument('--gdm', type=float, default=1.0, help='g_DM value')
	args = parser.parse_args()

	# Load xs
	limit_collection = {}
	# 8 TeV bb
	f_Xbb_8TeVl = TFile("~/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_trigbbl_CSVTM_ZPrime_dijet4.root", "READ")
	gr_Xbb_8TeVl_exp = f_Xbb_8TeVl.Get("graph_exp")
	BR_bbl_8TeV = []
	for i in xrange(gr_Xbb_8TeVl_exp.GetN()):
		mass = gr_Xbb_8TeVl_exp.GetX()[i]
		BR_bbl_8TeV.append(ZpBranchingRatio(mass, selected_decays=["b"]))
		#print "[debug] BR(Z-->bb, M={} GeV) = {}".format(mass, BR_bb_8TeVl[-1])
	print "Loading 2012bbl expected limits"
	limit_collection["CMS bbl, 2012 exp"] = LimitGraph()
	limit_collection["CMS bbl, 2012 exp"].LoadXSBR(gr_Xbb_8TeVl_exp, BR_bbl_8TeV)

	gr_Xbb_8TeVl_exp1sig = f_Xbb_8TeVl.Get("graph_exp_1sigma")
	BR_Xbb_8TeVl_exp1sig = []
	for i in xrange(gr_Xbb_8TeVl_exp1sig.GetN()):
		mass = gr_Xbb_8TeVl_exp1sig.GetX()[i]
		BR_Xbb_8TeVl_exp1sig.append(ZpBranchingRatio(mass, selected_decays=["b"]))
	print "Loading 2012bbl expected1sig limits"
	limit_collection["CMS bbl, 2012 exp1sig"] = LimitGraph()
	limit_collection["CMS bbl, 2012 exp1sig"].LoadXSBR(gr_Xbb_8TeVl_exp1sig, BR_Xbb_8TeVl_exp1sig)

	gr_Xbb_8TeVl_exp2sig = f_Xbb_8TeVl.Get("graph_exp_2sigma")
	BR_Xbb_8TeVl_exp2sig = []
	for i in xrange(gr_Xbb_8TeVl_exp2sig.GetN()):
		mass = gr_Xbb_8TeVl_exp2sig.GetX()[i]
		BR_Xbb_8TeVl_exp2sig.append(ZpBranchingRatio(mass, selected_decays=["b"]))
	print "Loading 2012bbl expected2sig limits"
	limit_collection["CMS bbl, 2012 exp2sig"] = LimitGraph()
	limit_collection["CMS bbl, 2012 exp2sig"].LoadXSBR(gr_Xbb_8TeVl_exp2sig, BR_Xbb_8TeVl_exp2sig)

	gr_Xbb_8TeVl_obs = f_Xbb_8TeVl.Get("graph_obs")
	print "Loading 2012bbl observed limits"
	limit_collection["CMS bbl, 2012 obs"] = LimitGraph()
	limit_collection["CMS bbl, 2012 obs"].LoadXSBR(gr_Xbb_8TeVl_obs, BR_bbl_8TeV)

	f_Xbb_8TeVh = TFile("~/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_trigbbh_CSVTM_ZPrime_dijet4.root", "READ")
	gr_Xbb_8TeVh_exp = f_Xbb_8TeVh.Get("graph_exp")
	BR_bbh_8TeV = []
	for i in xrange(gr_Xbb_8TeVh_exp.GetN()):
		mass = gr_Xbb_8TeVh_exp.GetX()[i]
		BR_bbh_8TeV.append(ZpBranchingRatio(mass, selected_decays=["b"]))
	print "Loading 2012bbh expected limits"
	limit_collection["CMS bbh, 2012 exp"] = LimitGraph()
	limit_collection["CMS bbh, 2012 exp"].LoadXSBR(gr_Xbb_8TeVh_exp, BR_bbh_8TeV)

	gr_Xbb_8TeVh_exp1sig = f_Xbb_8TeVh.Get("graph_exp_1sigma")
	BR_Xbb_8TeVh_exp1sig = []
	for i in xrange(gr_Xbb_8TeVh_exp1sig.GetN()):
		mass = gr_Xbb_8TeVh_exp1sig.GetX()[i]
		BR_Xbb_8TeVh_exp1sig.append(ZpBranchingRatio(mass, selected_decays=["b"]))
	print "Loading 2012bbh expected1sig limits"
	limit_collection["CMS bbh, 2012 exp1sig"] = LimitGraph()
	limit_collection["CMS bbh, 2012 exp1sig"].LoadXSBR(gr_Xbb_8TeVh_exp1sig, BR_Xbb_8TeVh_exp1sig)

	gr_Xbb_8TeVh_exp2sig = f_Xbb_8TeVh.Get("graph_exp_2sigma")
	BR_Xbb_8TeVh_exp2sig = []
	for i in xrange(gr_Xbb_8TeVh_exp2sig.GetN()):
		mass = gr_Xbb_8TeVh_exp2sig.GetX()[i]
		BR_Xbb_8TeVh_exp2sig.append(ZpBranchingRatio(mass, selected_decays=["b"]))
	print "Loading 2012bbh expected2sig limits"
	limit_collection["CMS bbh, 2012 exp2sig"] = LimitGraph()
	limit_collection["CMS bbh, 2012 exp2sig"].LoadXSBR(gr_Xbb_8TeVh_exp2sig, BR_Xbb_8TeVh_exp2sig)

	gr_Xbb_8TeVh_obs = f_Xbb_8TeVh.Get("graph_obs")
	print "Loading 2012bbh observed limits"
	limit_collection["CMS bbh, 2012 obs"] = LimitGraph()
	limit_collection["CMS bbh, 2012 obs"].LoadXSBR(gr_Xbb_8TeVh_obs, BR_bbh_8TeV)

	gr2D_bbl_obs = make_mdm_mZp_graph(limit_collection["CMS bbl, 2012 obs"].GetXS(), args.gdm)
	gr2D_bbl_obs.SetName("gr2D_bbl_obs")
	gr2D_bbl_exp = make_mdm_mZp_graph(limit_collection["CMS bbl, 2012 exp"].GetXS(), args.gdm)
	gr2D_bbl_exp.SetName("gr2D_bbl_exp")
	gr2D_bbh_obs = make_mdm_mZp_graph(limit_collection["CMS bbh, 2012 obs"].GetXS(), args.gdm)
	gr2D_bbh_obs.SetName("gr2D_bbh_obs")
	gr2D_bbh_exp = make_mdm_mZp_graph(limit_collection["CMS bbh, 2012 exp"].GetXS(), args.gdm)
	gr2D_bbh_exp.SetName("gr2D_bbh_exp")


	# Contours
	if args.contour:
		gr2D_bbl_obs_copy = gr2D_bbl_obs.Clone()	
		gr2D_bbl_obs_copy.GetHistogram().SetContour(1, array('d', [args.contour]))
		contours_bbl_obs = gr2D_bbl_obs_copy.GetContourList(args.contour)
		if contours_bbl_obs:
			print contours_bbl_obs.GetSize()
			for i in xrange(contours_bbl_obs.GetSize()):
				contour_graph = contours_bbl_obs.At(i)
				contour_graph.SetLineColor(1)
				contour_graph.SetLineStyle(1)
				contour_graph.SetLineWidth(2)

		gr2D_bbl_exp_copy = gr2D_bbl_exp.Clone()			
		gr2D_bbl_exp_copy.GetHistogram().SetContour(1, array('d', [args.contour]))
		contours_bbl_exp = gr2D_bbl_exp_copy.GetContourList(args.contour)
		if contours_bbl_exp:
			print contours_bbl_exp.GetSize()
			for i in xrange(contours_bbl_exp.GetSize()):
				contour_graph = contours_bbl_exp.At(i)
				contour_graph.SetLineColor(kGray)
				contour_graph.SetLineStyle(2)
				contour_graph.SetLineWidth(2)

		gr2D_bbh_obs_copy = gr2D_bbh_obs.Clone()			
		gr2D_bbh_obs_copy.GetHistogram().SetContour(1, array('d', [args.contour]))
		contours_bbh_obs = gr2D_bbh_obs_copy.GetContourList(args.contour)
		if contours_bbh_obs:
			print contours_bbh_obs.GetSize()
			for i in xrange(contours_bbh_obs.GetSize()):
				contour_graph = contours_bbh_obs.At(i)
				contour_graph.SetLineColor(1)
				contour_graph.SetLineStyle(1)
				contour_graph.SetLineWidth(2)

		gr2D_bbh_exp_copy = gr2D_bbh_exp.Clone()			
		gr2D_bbh_exp_copy.GetHistogram().SetContour(1, array('d', [args.contour]))
		contours_bbh_exp = gr2D_bbh_exp_copy.GetContourList(args.contour)
		if contours_bbh_exp:
			print contours_bbh_exp.GetSize()
			for i in xrange(contours_bbh_exp.GetSize()):
				contour_graph = contours_bbh_exp.At(i)
				contour_graph.SetLineColor(kGray)
				contour_graph.SetLineStyle(2)
				contour_graph.SetLineWidth(2)

	# Drawing
	c = TCanvas("c_gq_vs_mdm_vs_mZp", "c_gq_vs_mdm_vs_mZp", 800, 600)
	c.SetRightMargin(0.25)
	c.SetLogz()
	frame = TH2D("frame", "frame", 100, 0., 1500., 100, 0., 1500.)
	frame.GetXaxis().SetTitle("m_{Z'} [GeV]")
	frame.GetYaxis().SetTitle("m_{#chi} [GeV]")
	frame.GetZaxis().SetTitle("g_{q}, 95% CL upper limit")
	frame.SetMinimum(0.05)
	frame.SetMaximum(5.)
	frame.Draw("axis colz")

	gr2D_bbl_obs.SetMinimum(0.05)
	gr2D_bbl_obs.SetMaximum(5.)
	gr2D_bbl_obs.GetZaxis().SetTitle("g_{q}, 95% CL upper limit")
	gr2D_bbl_obs.Draw("colz same")
	gr2D_bbh_obs.SetMinimum(0.05)
	gr2D_bbh_obs.SetMaximum(5.)
	gr2D_bbh_obs.Draw("col same")
	if args.contour:
		print contours_bbl_obs
		print contours_bbl_exp
		print contours_bbh_obs
		print contours_bbh_exp
		l = TLegend(0.2, 0.65, 0.4, 0.8)
		l.SetBorderSize(1)
		l.SetFillColor(0)
		first_obs = True
		first_exp = True		
		if contours_bbl_obs:
			for i in xrange(contours_bbl_obs.GetSize()):
				contour_graph = contours_bbl_obs.At(i)
				contour_graph.Draw("l")
				if first_obs:
					l.AddEntry(contour_graph, "g_{{q}}={}, obs ".format(args.contour), "l")
					first_obs = False
		if contours_bbl_exp:
			for i in xrange(contours_bbl_exp.GetSize()):
				contour_graph = contours_bbl_exp.At(i)
				contour_graph.Draw("l")
				if first_exp: 
					l.AddEntry(contour_graph, "g_{{q}}={}, exp ".format(args.contour), "l")
					first_exp = False
		if contours_bbh_obs:
			for i in xrange(contours_bbh_obs.GetSize()):
				contour_graph = contours_bbh_obs.At(i)
				contour_graph.Draw("l")
				if first_obs:
					l.AddEntry(contour_graph, "g_{{q}}={}, obs ".format(args.contour), "l")
					first_obs = False
		if contours_bbh_exp:
			for i in xrange(contours_bbh_exp.GetSize()):
				contour_graph = contours_bbh_exp.At(i)
				contour_graph.Draw("l")
				if first_exp: 
					l.AddEntry(contour_graph, "g_{{q}}={}, exp ".format(args.contour), "l")
					first_exp = False
		l.Draw()
	frame.Draw("axis same")

	c.SaveAs("{}/{}.pdf".format(analysis_config.figure_directory, c.GetName()))
