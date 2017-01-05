#!/usr/bin/env python

import sys
import os
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

class LimitGraph:
	def __init__(self):
		self._xs_limit = None
		self._xsBR_limit = None
		self._gB_limit = None

	def LoadXS(self, xs_limit, br=None):
		self._xs_limit = xs_limit
		self._gB_limit = ConvertXSToGB(self._xs_limit)
		if br:
			if len(br) != xs_limit.GetN():
				print "[LimitGraph] ERROR : len(br)={} != xs_limit.GetN()={}".format(len(br), xs_limit.GetN())
				sys.exit(1)
			self._xsBR_limit = self._xs_limit.Clone()
			self._xsBR_limit.SetName(self._xsBR_limit.GetName() + "_timesBR")
			for i in xrange(self._xsBR_limit.GetN()):
				self._xsBR_limit.SetPoint(i, self._xsBR_limit.GetX()[i], self._xsBR_limit.GetY()[i] * br[i])

	def LoadGB(self, gB_limit, br=None):
		self._gB_limit = gB_limit
		self._xs_limit = ConvertGBToXS(self._gB_limit)
		if br:
			if len(br) != gB_limit.GetN():
				print "[LimitGraph] ERROR : len(br)={} != gB_limit.GetN()={}".format(len(br), gB_limit.GetN())
				sys.exit(1)
			self._xsBR_limit = self._xs_limit.Clone()
			self._xsBR_limit.SetName(self._xsBR_limit.GetName() + "_timesBR")
			for i in xrange(self._xsBR_limit.GetN()):
				self._xsBR_limit.SetPoint(i, self._xsBR_limit.GetX()[i], self._xsBR_limit.GetY()[i] * br[i])

	def LoadXSBR(self, xsBR_limit, br):
		if len(br) != xsBR_limit.GetN():
			print "[LimitGraph] ERROR : len(br)={} != xsBR_limit.GetN()={}".format(len(br), xsBR_limit.GetN())
			sys.exit(1)
		self._xsBR_limit = xsBR_limit
		self._xs_limit = self._xsBR_limit.Clone()
		self._xs_limit.SetName(self._xs_limit.GetName() + "_noBR")
		for i in xrange(self._xs_limit.GetN()):
			if br[i] > 0:
				self._xs_limit.SetPoint(i, self._xs_limit.GetX()[i], self._xs_limit.GetY()[i] / br[i])
			else:
				self._xs_limit.SetPoint(i, self._xs_limit.GetX()[i], 0.)
		self._gB_limit = ConvertXSToGB(self._xs_limit)

	def GetXS(self):
		return self._xs_limit

	def GetXSBR(self):
		return self._xsBR_limit

	def GetGB(self):
		return self._gB_limit


# Computed cross section times BR(jj) of dark matter model
dm_x = array('d', [528.7733914858716,631.2905390292997,795.7485275479014,980.9839334973528,1166.2193394468052,1331.0361216730032,1578.3102214269736,1825.7147916200695,2073.0541265936026,2217.3544322672033,2454.5170729888905,2701.98687840155,2928.9075896518298,3104.1946246179077,3310.4357712666806,3475.513494371131,3568.343211809438])
dm_y = array('d', [57.269239965262344,26.146418729332762,10.943588903119819,4.998039360562294,2.2826513012206147,1.3151916426227461,0.6008150098368491,0.30829460763931776,0.1492637583468037,0.09943936422165509,0.0556698333169542,0.03027483510044679,0.018491825891702274,0.011967992892333293,0.007309405483957529,0.005313471607526083,0.004337281648532702])
dm_sigma_times_BRjj = TGraph(len(dm_x), dm_x, dm_y)

# Parton luminosity ratio 8/13
parton_lumi_qq_x = array('d', [212.44114380327838,256.61539097505755,189.6759157071346,166.63007851997236,148.77397175424056,129.29446527108541,231.60316570875963,293.6875978919667,281.275999744092,356.6759221822652,370.4099179447112,395.1955101394775,408.20347685923906,330.71672989214755,318.45447138645915,459.6711606703674,495.7524682412154,561.2799896042493,621.896007068797,621.896007068797,707.9075738492571,731.2085424206072,775.9371322362275,810.176171052523,869.066435119888,897.6719910415933,994.6170133018367,1078.4937378238508,1126.0834038140092,1207.936402435694,1268.06383234535,1360.2371356456406,1405.0097075435822,1507.1373628002914,1556.745195264016,1678.9399444796177,1743.5884184831468,1840.2852206152222])
parton_lumi_qq_y = array('d', [1. / x for x in [1.9084845217716158,1.9701099766065029,1.8782550849899042,1.8483835132927056,1.8191060463151,1.7900970823130031,1.938989007634503,2.017576452165467,2.0016431192325372,2.09898276099542,2.115644783157695,2.149272912557035,2.1662868937909714,2.0660513546505883,2.049779911668458,2.235602848865901,2.2889166804078145,2.380596053538979,2.47573150572491,2.47573150572491,2.615193824659969,2.656414604730997,2.7407557944954624,2.805749111310745,2.9174930594233572,2.98654717428537,3.203741702660958,3.409888265666115,3.54530646023206,3.773266592268887,3.9537405526141987,4.273729358185959,4.443259768096302,4.840253677824623,5.071429377205673,5.6546923942954574,5.971016709192964,6.554715672402531]])
parton_lumi_qq = TGraph(len(parton_lumi_qq_x), parton_lumi_qq_x, parton_lumi_qq_y)


def PartonLuminosityRatio(x, initial_state):
	ratio = -1.
	if initial_state == "gg":
		ratio = (1.43368e-21*(x**6)+-2.037e-17*(x**5)+1.14943e-13*(x**4)+-3.31131e-10*(x**3)+5.43977e-07*x*x+-0.000607347*x+0.483842)
	elif initial_state == "qq":
		ratio = parton_lumi_qq.Eval(x)
	else:
		print "[PartonLuminosityRatio] ERROR : initial_state={} not known.".format(initial_state)
		sys.exit(1)
	return ratio

# Branching ratios of vector/axial vector Z' to quarks, assuming equal coupling to quarks and no other decays
quark_masses = {
	"u":2.3e-3,
	"d":4.8e-3,
	"c":1.275,
	"s":0.095,
	"t":173.210,
	"b":4.180
}
def ZpBranchingRatio(m_med, selected_decays=["b"], zp_type="vector"):
	total_width = 0.
	selected_width = 0.
	for quark in ["u", "d", "s", "c", "b", "t"]:
		z = (quark_masses[quark] / m_med)**2
		if z > 0.25:
			continue
		if zp_type == "vector":
			this_width = (1 - 4*z)**0.5 * (1 + 2*z)
		elif zp_type == "avector":
			this_width = (1 - 4*z)**1.5
		else:
			print "[ZpBranchingRatio] ERROR : zp_type must be vector or avector"
			sys.exit(1)
		total_width += this_width
		if quark in selected_decays:
			selected_width += this_width
	return selected_width / total_width


# Convert a cross section TGraph to a g_B TGraph
# xs_limit = TGraph with cross section (*BR(jj)) limits
# gDM = Z'-DM coupling
# gQ = Z'-quark coupling
# mDM = dark matter mass
def ConvertXSToGB(xs_limit, gDM=1, gQ=0.25, mDM=1.):
	gb_x = array('d', [])
	gb_y = array('d', [])
	for i in xrange(xs_limit.GetN()):
		this_mass = xs_limit.GetX()[i]
		this_xs = xs_limit.GetY()[i]
		width_DM_total = avtotwidth(2, gDM,gQ,this_mass,mDM)
		width_qq   = avtotwidth(2, 0., gQ,this_mass,mDM)
		this_gB = math.sqrt(this_xs / (dm_sigma_times_BRjj.Eval(this_mass)*ZpBranchingRatio(this_mass, selected_decays=["u","d","s","c","b"])) * width_qq / width_DM_total) * 0.25 * 6.
		gb_x.append(this_mass)
		gb_y.append(this_gB)
	return TGraph(len(gb_x), gb_x, gb_y)

def ConvertGBToXS(gb_limit, gDM=1, gQ=0.25, mDM=1.):
	xs_x = array('d', [])
	xs_y = array('d', [])
	for i in xrange(gb_limit.GetN()):
		this_mass = gb_limit.GetX()[i]
		this_gb = gb_limit.GetY()[i]
		width_DM_total = avtotwidth(2, gDM,gQ,this_mass,mDM)
		width_qq   = avtotwidth(2, 0., gQ,this_mass,mDM)
		this_xs = (this_gb / (0.25 * 6.))**2 * width_DM_total / width_qq * (dm_sigma_times_BRjj.Eval(this_mass)*ZpBranchingRatio(this_mass, selected_decays=["u","d","s","c","b"]))
		xs_x.append(this_mass)
		xs_y.append(this_xs)
	return TGraph(len(xs_x), xs_x, xs_y)

def csvToGraph(fn, linecolor=1, addFactor=False):

	factor = 1.;
	if addFactor: factor = 6.;

	a_m = array('d', []);
	a_g = array('d', []);

	ifile = open(fn,'r');
	npoints = 0;
	for line in ifile: 
		lline = line.strip().split(',');
		a_m.append(float(lline[0]))
		a_g.append(float(lline[1])*factor)
		npoints += 1;

	gr = ROOT.TGraph(npoints,a_m,a_g);
	gr.SetLineColor(linecolor);
	gr.SetLineWidth(2);

	return gr

def avwidth(iType,g,med,mfm):
    front=g*g*(med*med+2*mfm*mfm)/(12*med*3.14159265)
    if abs(iType) == 2:
        front=g*g*(med*med-4*mfm*mfm)/(12*med*3.14159265)
    if 2.*mfm > med:
        return 0.001
    sqrtV=math.sqrt(1-4*(mfm/med)*(mfm/med))
    return front*sqrtV

def avtotwidth(iType,gdm,gsm,med,mdm):
    u=avwidth(iType,gsm,med,0.001)
    d=u
    s=avwidth(iType,gsm,med,0.135)
    c=avwidth(iType,gsm,med,1.5)
    b=avwidth(iType,gsm,med,5.1)
    t=0
    if med > 2.*172.5:
        t=avwidth(iType,gsm,med,172.5)
    quarks=3*(u+d+s+c+b+t)
    dm=avwidth(iType,gdm,med,mdm)
    #print u,d,s,c,b,t,dm,quarks
    return dm+quarks


def MakeLimitPlot(limit_names, limit_graphs, save_tag, what="gB", logx=False, logy=False, line_styles={}, line_colors={}, x_range=None, x_title=None, y_title=None):
	if not what in ["gB", "xs", "xsBR"]:
		print "[MakeLimitPlot] ERROR : what must be gB, xs, or xsBR"
		sys.exit(1)

	c = TCanvas("c_{}".format(save_tag), "c_{}".format(save_tag), 800, 600)
	if logx:
		c.SetLogx()
	if logy:
		c.SetLogy()
	#if what == "gB":
	#	l = TLegend(0.17, 0.6, 0.37, 0.9)
	#else:
	l = TLegend(0.64, 0.6, 0.9, 0.9)
	l.SetFillColor(0)
	l.SetBorderSize(0)

	if x_range:
		x_min = x_range[0]
		x_max = x_range[1]
	else:
		x_min = 70.
		x_max = 2500.
	y_min = 1.e20
	y_max = -1.e20
	for name in limit_names:
		for i in xrange(limit_graphs[name].GetGB().GetN()):
			this_x = limit_graphs[name].GetGB().GetX()[i]
			if what == "gB":
				this_y = limit_graphs[name].GetGB().GetY()[i]
			elif what == "xs":
				this_y = limit_graphs[name].GetXS().GetY()[i]
			elif what == "xsBR":
				this_y = limit_graphs[name].GetXSBR().GetY()[i]
			if this_y > y_max:
				y_max = this_y
			if this_y < y_min:
				y_min = this_y

	if what == "gB":
		if logy:
			y_min = 0.1
			y_max = 10.
		else:
			y_min = 0.
			y_max = 5.
	else:
		if logy:
			if y_min < 0.:
				y_min = 0.01
			else:
				y_min = y_min / 5.
			y_max = y_max * 10.
		else:
			y_max = y_max + 0.125 * (y_max - y_min)
			y_min = y_min - 0.125 * (y_max - y_min)

	frame = TH1D("frame", "frame", 100, x_min, x_max)
	frame.SetMinimum(y_min)
	frame.SetMaximum(y_max)
	if x_title:
		frame.GetXaxis().SetTitle(x_title)
	else:
		frame.GetXaxis().SetTitle("Mediator Mass [GeV]")
	
	if y_title:
		frame.GetYaxis().SetTitle(y_title)
	else:
		if what == "gB":
			frame.GetYaxis().SetTitle("g_{B}")
		elif what == "xs":
			frame.GetYaxis().SetTitle("#sigma [pb]")
		elif what == "xsBR":
			frame.GetYaxis().SetTItle("#sigma #times BR [pb]")
	print "[debug] frame min/max = {}/{}".format(frame.GetXaxis().GetXmin(), frame.GetXaxis().GetXmax())
	frame.Draw("axis")

	style_counter = 0
	for name in limit_names:
		if what == "gB":
			if name in line_colors:
				limit_graphs[name].GetGB().SetLineColor(line_colors[name])
			else:
				limit_graphs[name].GetGB().SetLineColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(limit_names)))
			if name in line_widths:
				limit_graphs[name].GetGB().SetLineWidth(line_widths[name])
			else:
				limit_graphs[name].GetGB().SetLineWidth(2)
			if name in line_styles:
				limit_graphs[name].GetGB().SetLineStyle(line_styles[name])
			limit_graphs[name].GetGB().Draw("l")
			l.AddEntry(limit_graphs[name].GetGB(), name, "l")
		elif what == "xs":
			if name in line_colors:
				limit_graphs[name].GetXS().SetLineColor(line_colors[name])
			else:
				limit_graphs[name].GetXS().SetLineColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(limit_names)))
			limit_graphs[name].GetXS().SetLineWidth(2)
			if name in line_styles:
				limit_graphs[name].GetXS().SetLineStyle(line_styles[name])
			limit_graphs[name].GetXS().Draw("l")
			l.AddEntry(limit_graphs[name].GetXS(), name, "l")
		elif what == "xsBR":
			if name in line_colors:
				limit_graphs[name].GetXSBR().SetLineColor(line_colors[name])
			else:
				limit_graphs[name].GetXSBR().SetLineColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(limit_names)))
			limit_graphs[name].GetXSBR().SetLineWidth(2)
			if name in line_styles:
				limit_graphs[name].GetXSBR().SetLineStyle(line_styles[name])
			limit_graphs[name].GetXSBR().Draw("l")
			l.AddEntry(limit_graphs[name].GetXSBR(), name, "l")
		style_counter += 1

	l.Draw()
	c.SaveAs("{}/{}.pdf".format(analysis_config.figure_directory, c.GetName()))

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description = 'Dijet mass spectrum fits')
	parser.add_argument('--bb_comparison', action='store_true', help='Compare 8 TeV and 13 TeV bb analyses')
	parser.add_argument('--gB', action='store_true', help='Compare all analyses with gB')

	args = parser.parse_args()


	limit_names = []
	limit_collection = {}

	# Load external data
	gr_UA2 = csvToGraph(analysis_config.external_limits + "/UA2.csv",4,True);
	limit_collection["UA2"] = LimitGraph()
	limit_collection["UA2"].LoadGB(gr_UA2)
	gr_CDFRun1 = csvToGraph(analysis_config.external_limits + "/CDF_Run1.csv",2,True );
	limit_collection["CDFRun1"] = LimitGraph()
	limit_collection["CDFRun1"].LoadGB(gr_CDFRun1)
	gr_CDFRun2 = csvToGraph(analysis_config.external_limits + "/CDF_Run2.csv",6,True );
	limit_collection["CDFRun2"] = LimitGraph()
	limit_collection["CDFRun2"].LoadGB(gr_CDFRun2)
	gr_ATLAS = csvToGraph(analysis_config.external_limits + "/gBMZB_ATLAS_all_fbinv.csv",7,False );
	limit_collection["ATLAS jj"] = LimitGraph()
	limit_collection["ATLAS jj"].LoadGB(gr_ATLAS)
	gr_CMS = csvToGraph(analysis_config.external_limits + "/CMS_Scouting.csv",8,False );
	limit_collection["CMS jj"] = LimitGraph()
	limit_collection["CMS jj"].LoadGB(gr_CMS)
	gr_CMSBoostedDijet = csvToGraph(analysis_config.external_limits + "/BoostedDijet.csv")
	limit_collection["CMS boosted jj"] = LimitGraph()
	limit_collection["CMS boosted jj"].LoadGB(gr_CMSBoostedDijet)

	# HIG-16-025 (2015/13 TeV X->bb)
	# Raw limit values are sigma*BR(bb). Assuming equal BR to all quarks, this should be multiplied by 5-6 depending on ttbar threshold effects. 
	x_Xbb_13TeV_obs = array('d', [550,600,650,700,750,800,850,900,1000,1100,1200])
	y_Xbb_13TeV_obs = array('d', [11.1814,5.15109,5.68285,5.6275,5.09532,4.8379,5.24706,4.38625,2.11611,2.06881,4.6223])
	BR_bb_13TeV = []
	for mass in x_Xbb_13TeV_obs:
		BR_bb_13TeV.append(ZpBranchingRatio(mass, selected_decays=["b"]))
	for i in xrange(len(y_Xbb_13TeV_obs)):
		mass = x_Xbb_13TeV_obs[i]
		y_Xbb_13TeV_obs[i] = y_Xbb_13TeV_obs[i]# / ZpBranchingRatio(mass, selected_decays=["b"])
	gr_Xbb_13TeV_obs = TGraph(len(x_Xbb_13TeV_obs), x_Xbb_13TeV_obs, y_Xbb_13TeV_obs)
	limit_collection["CMS bb, 2015 obs"] = LimitGraph()
	limit_collection["CMS bb, 2015 obs"].LoadXSBR(gr_Xbb_13TeV_obs, BR_bb_13TeV)

	x_Xbb_13TeV_exp = array('d', [550,600,650,700,750,800,850,900,1000,1100,1200])
	y_Xbb_13TeV_exp = array('d', [12.5625,8.34375,6.96875,5.79688,4.79688,4.26562,3.75,3.54688,3.20312,3.02344,2.96094])
	for i in xrange(len(y_Xbb_13TeV_exp)):
		mass = x_Xbb_13TeV_exp[i]
		y_Xbb_13TeV_exp[i] = y_Xbb_13TeV_exp[i] #/ ZpBranchingRatio(mass, selected_decays=["b"])
	gr_Xbb_13TeV_exp = TGraph(len(x_Xbb_13TeV_exp), x_Xbb_13TeV_exp, y_Xbb_13TeV_exp)
	limit_collection["CMS bb, 2015 exp"] = LimitGraph()
	limit_collection["CMS bb, 2015 exp"].LoadXSBR(gr_Xbb_13TeV_exp, BR_bb_13TeV)

	line_styles = {
		"UA2":2,
		"CDFRun1":2,
		"CDFRun2":2,
		"ATLAS jj":2,
		"CMS boosted jj":1,
		"CMS jj":1,
		"CMS bb, 2015 exp":3,
		"CMS bb, 2012 low exp":3,
		"CMS bb, 2012 high exp":3,
		"CMS bb, 2015 obs":1,
		"CMS bb, 2012 low obs":1,
		"CMS bb, 2012 high obs":1,
	}
	line_widths = {
		"CMS bb, 2012 low exp":4,
		"CMS bb, 2012 high exp":4,
		"CMS bb, 2012 low obs":4,
		"CMS bb, 2012 high obs":4,
	}
	line_colors = {
		"UA2":seaborn.GetColorRoot("cubehelixlarge", 1, 6),
		"CDFRun1":seaborn.GetColorRoot("cubehelixlarge", 2, 6),
		"CDFRun2":seaborn.GetColorRoot("cubehelixlarge", 3, 6),
		"ATLAS jj":seaborn.GetColorRoot("cubehelixlarge", 4, 6),
		"CMS jj":seaborn.GetColorRoot("cubehelixlarge", 5, 6),
		"CMS boosted jj":seaborn.GetColorRoot("default", 5, 6),
		"CMS bb, 2015 exp":seaborn.GetColorRoot("default", 1),
		"CMS bb, 2015 obs":seaborn.GetColorRoot("default", 1),
		"CMS bb, 2012 low exp":seaborn.GetColorRoot("default", 2),
		"CMS bb, 2012 low obs":seaborn.GetColorRoot("default", 2),
		"CMS bb, 2012 high exp":seaborn.GetColorRoot("default", 3),
		"CMS bb, 2012 high obs":seaborn.GetColorRoot("default", 3),
	}

	if args.bb_comparison:
		# 8 TeV bb
		f_Xbb_8TeVl = TFile("~/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_trigbbl_CSVTM_RSG_f4.root", "READ")
		gr_Xbb_8TeVl_exp = f_Xbb_8TeVl.Get("graph_exp")
		BR_bb_8TeVl = []
		for i in xrange(gr_Xbb_8TeVl_exp.GetN()):
			mass = gr_Xbb_8TeVl_exp.GetX()[i]
			gr_Xbb_8TeVl_exp.SetPoint(i, mass, gr_Xbb_8TeVl_exp.GetY()[i] / PartonLuminosityRatio(mass, "gg")) # / ZpBranchingRatio(mass, selected_decays=["b"]) 
			BR_bb_8TeVl.append(ZpBranchingRatio(mass, selected_decays=["b"]) )
		limit_collection["CMS bb, 2012 low exp"] = LimitGraph()
		limit_collection["CMS bb, 2012 low exp"].LoadXSBR(gr_Xbb_8TeVl_exp, BR_bb_8TeVl)

		gr_Xbb_8TeVl_obs = f_Xbb_8TeVl.Get("graph_obs")
		for i in xrange(gr_Xbb_8TeVl_obs.GetN()):
			mass = gr_Xbb_8TeVl_obs.GetX()[i]
			gr_Xbb_8TeVl_obs.SetPoint(i, mass, gr_Xbb_8TeVl_obs.GetY()[i] / PartonLuminosityRatio(mass, "gg"))
		limit_collection["CMS bb, 2012 low obs"] = LimitGraph()
		limit_collection["CMS bb, 2012 low obs"].LoadXSBR(gr_Xbb_8TeVl_obs, BR_bb_8TeVl)

		f_Xbb_8TeVh = TFile("~/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_trigbbh_CSVTM_Hbb_f1.root", "READ")
		gr_Xbb_8TeVh_exp = f_Xbb_8TeVh.Get("graph_exp")
		BR_bb_8TeVh = []
		for i in xrange(gr_Xbb_8TeVh_exp.GetN()):
			mass = gr_Xbb_8TeVh_exp.GetX()[i]
			gr_Xbb_8TeVh_exp.SetPoint(i, mass, gr_Xbb_8TeVh_exp.GetY()[i] / PartonLuminosityRatio(mass, "gg"))
			BR_bb_8TeVh.append(ZpBranchingRatio(mass, selected_decays=["b"]) )
		limit_collection["CMS bb, 2012 high exp"] = LimitGraph()
		limit_collection["CMS bb, 2012 high exp"].LoadXSBR(gr_Xbb_8TeVh_exp, BR_bb_8TeVh)

		gr_Xbb_8TeVh_obs = f_Xbb_8TeVh.Get("graph_obs")
		for i in xrange(gr_Xbb_8TeVh_obs.GetN()):
			mass = gr_Xbb_8TeVh_obs.GetX()[i]
			gr_Xbb_8TeVh_obs.SetPoint(i, mass, gr_Xbb_8TeVh_obs.GetY()[i] / PartonLuminosityRatio(mass, "gg"))
			print "[debug] Set point {} {}".format(mass, gr_Xbb_8TeVh_obs.GetY()[i])
		limit_collection["CMS bb, 2012 high obs"] = LimitGraph()
		limit_collection["CMS bb, 2012 high obs"].LoadXSBR(gr_Xbb_8TeVh_obs, BR_bb_8TeVh)

		limit_names = ["CMS bb, 2015 exp", "CMS bb, 2015 obs","CMS bb, 2012 low exp", "CMS bb, 2012 low obs","CMS bb, 2012 high exp", "CMS bb, 2012 high obs"]
		MakeLimitPlot(limit_names, limit_collection, "xs_combination_log", what="xsBR", x_title="m_{X} [GeV]", y_title="#sigma_{13 TeV} #times BR(b#bar{b}) [pb]", logy=True, line_styles=line_styles, line_colors=line_colors, x_range=[0,1500])

	elif args.gB:
		# 8 TeV bb
		f_Xbb_8TeVl = TFile("~/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_trigbbl_CSVTM_RSG_f4.root", "READ")
		gr_Xbb_8TeVl_exp = f_Xbb_8TeVl.Get("graph_exp")
		BR_bb_8TeVl = []
		for i in xrange(gr_Xbb_8TeVl_exp.GetN()):
			mass = gr_Xbb_8TeVl_exp.GetX()[i]
			gr_Xbb_8TeVl_exp.SetPoint(i, mass, gr_Xbb_8TeVl_exp.GetY()[i] / PartonLuminosityRatio(mass, "qq")) # / ZpBranchingRatio(mass, selected_decays=["b"]) 
			BR_bb_8TeVl.append(ZpBranchingRatio(mass, selected_decays=["b"]))
			#print "[debug] BR(Z-->bb, M={} GeV) = {}".format(mass, BR_bb_8TeVl[-1])
		limit_collection["CMS bb, 2012 low exp"] = LimitGraph()
		limit_collection["CMS bb, 2012 low exp"].LoadXSBR(gr_Xbb_8TeVl_exp, BR_bb_8TeVl)

		gr_Xbb_8TeVl_obs = f_Xbb_8TeVl.Get("graph_obs")
		for i in xrange(gr_Xbb_8TeVl_obs.GetN()):
			mass = gr_Xbb_8TeVl_obs.GetX()[i]
			gr_Xbb_8TeVl_obs.SetPoint(i, mass, gr_Xbb_8TeVl_obs.GetY()[i] / PartonLuminosityRatio(mass, "qq"))
		limit_collection["CMS bb, 2012 low obs"] = LimitGraph()
		limit_collection["CMS bb, 2012 low obs"].LoadXSBR(gr_Xbb_8TeVl_obs, BR_bb_8TeVl)

		f_Xbb_8TeVh = TFile("~/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_trigbbh_CSVTM_Hbb_f1.root", "READ")
		gr_Xbb_8TeVh_exp = f_Xbb_8TeVh.Get("graph_exp")
		BR_bb_8TeVh = []
		for i in xrange(gr_Xbb_8TeVh_exp.GetN()):
			mass = gr_Xbb_8TeVh_exp.GetX()[i]
			gr_Xbb_8TeVh_exp.SetPoint(i, mass, gr_Xbb_8TeVh_exp.GetY()[i] / PartonLuminosityRatio(mass, "qq"))
			BR_bb_8TeVh.append(ZpBranchingRatio(mass, selected_decays=["b"]) )
		limit_collection["CMS bb, 2012 high exp"] = LimitGraph()
		limit_collection["CMS bb, 2012 high exp"].LoadXSBR(gr_Xbb_8TeVh_exp, BR_bb_8TeVh)

		gr_Xbb_8TeVh_obs = f_Xbb_8TeVh.Get("graph_obs")
		for i in xrange(gr_Xbb_8TeVh_obs.GetN()):
			mass = gr_Xbb_8TeVh_obs.GetX()[i]
			gr_Xbb_8TeVh_obs.SetPoint(i, mass, gr_Xbb_8TeVh_obs.GetY()[i] / PartonLuminosityRatio(mass, "qq"))
		limit_collection["CMS bb, 2012 high obs"] = LimitGraph()
		limit_collection["CMS bb, 2012 high obs"].LoadXSBR(gr_Xbb_8TeVh_obs, BR_bb_8TeVh)

		limit_names = ["UA2","CDFRun1","CDFRun2","ATLAS jj","CMS jj", "CMS boosted jj", "CMS bb, 2015 obs","CMS bb, 2012 low obs", "CMS bb, 2012 low exp","CMS bb, 2012 high obs", "CMS bb, 2012 high exp"]

		#MakeLimitPlot(limit_names, limit_collection, "gB_combination_log", what="gB", logy=True, line_styles=line_styles)
		#MakeLimitPlot(limit_names, limit_collection, "xs_combination", what="xs", line_styles=line_styles, line_colors=line_colors, x_range=[0,1500])
		MakeLimitPlot(limit_names, limit_collection, "gB_combination", what="gB", line_styles=line_styles, line_colors=line_colors, x_title="m_{Z'} [GeV]", y_title="g_{B}", x_range=[0,1000])
		MakeLimitPlot(limit_names, limit_collection, "gB_combination_log", what="gB", line_styles=line_styles, line_colors=line_colors, x_title="m_{Z'} [GeV]", y_title="g_{B}", x_range=[0,1000], logy=True)
