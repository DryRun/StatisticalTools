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



class LimitGraph:
	def __init__(self):
		self._xs_limit = None
		self._xsBR_limit = None
		self._gB_limit = None
		self._gq_limit = None

	def LoadXS(self, xs_limit, br=None):
		self._xs_limit = xs_limit
		self._gB_limit = ConvertXSToGB(self._xs_limit)
		self._gq_limit = ConvertXSToGQ(self._xs_limit)
		if br:
			if len(br) != xs_limit.GetN():
				print "[LimitGraph] ERROR : len(br)={} != xs_limit.GetN()={}".format(len(br), xs_limit.GetN())
				sys.exit(1)
			self._xsBR_limit = self._xs_limit.Clone()
			self._xsBR_limit.SetName(self._xsBR_limit.GetName() + "_timesBR")
			for i in xrange(self._xsBR_limit.GetN()):
				self._xsBR_limit.SetPoint(i, self._xsBR_limit.GetX()[i], self._xsBR_limit.GetY()[i] * br[i])

	def LoadGB(self, gB_limit, br=None, initial_states=["uu","dd"], final_states=["uu","dd","cc","ss","bb"]):
		self._gB_limit = gB_limit
		self._gq_limit = gB_limit.Clone()
		for i in xrange(self._gq_limit.GetN()):
			self._gq_limit.SetPoint(i, self._gq_limit.GetX()[i], self._gq_limit.GetY()[i] / 6.)
		self._xs_limit = ConvertGBToXS(self._gB_limit)
		if br:
			if len(br) != gB_limit.GetN():
				print "[LimitGraph] ERROR : len(br)={} != gB_limit.GetN()={}".format(len(br), gB_limit.GetN())
				sys.exit(1)
			self._xsBR_limit = self._xs_limit.Clone()
			self._xsBR_limit.SetName(self._xsBR_limit.GetName() + "_timesBR")
			for i in xrange(self._xsBR_limit.GetN()):
				self._xsBR_limit.SetPoint(i, self._xsBR_limit.GetX()[i], self._xsBR_limit.GetY()[i] * br[i])

	def LoadGQ(self, gq_limit, br=None, initial_states=["uu","dd"], final_states=["uu","dd","cc","ss","bb"]):
		self._gq_limit = gq_limit
		self._gB_limit = gq_limit.Clone()
		for i in xrange(self._gB_limit.GetN()):
			self._gB_limit.SetPoint(i, self._gq_limit.GetX()[i], self._gq_limit.GetY()[i] * 6.)
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
				self._xs_limit.SetPoint(i, self._xs_limit.GetX()[i], self._xsBR_limit.GetY()[i] / br[i])
				print "xs({}) = xsbr/br  = {}/{} = {}".format(self._xs_limit.GetX()[i], self._xsBR_limit.GetY()[i], br[i], self._xsBR_limit.GetY()[i] / br[i])
			else:
				self._xs_limit.SetPoint(i, self._xs_limit.GetX()[i], 0.)
		self._gB_limit = ConvertXSToGB(self._xs_limit)
		self._gq_limit = ConvertXSToGQ(self._xs_limit)

	def GetXS(self):
		return self._xs_limit

	def GetXSBR(self):
		return self._xsBR_limit

	def GetGB(self):
		return self._gB_limit

	def GetGQ(self):
		return self._gq_limit


# Parton luminosity ratio 8/13
#parton_lumi_qq_x = array('d', [212.44114380327838,256.61539097505755,189.6759157071346,166.63007851997236,148.77397175424056,129.29446527108541,231.60316570875963,293.6875978919667,281.275999744092,356.6759221822652,370.4099179447112,395.1955101394775,408.20347685923906,330.71672989214755,318.45447138645915,459.6711606703674,495.7524682412154,561.2799896042493,621.896007068797,621.896007068797,707.9075738492571,731.2085424206072,775.9371322362275,810.176171052523,869.066435119888,897.6719910415933,994.6170133018367,1078.4937378238508,1126.0834038140092,1207.936402435694,1268.06383234535,1360.2371356456406,1405.0097075435822,1507.1373628002914,1556.745195264016,1678.9399444796177,1743.5884184831468,1840.2852206152222])
#parton_lumi_qq_y = array('d', [1. / x for x in [1.9084845217716158,1.9701099766065029,1.8782550849899042,1.8483835132927056,1.8191060463151,1.7900970823130031,1.938989007634503,2.017576452165467,2.0016431192325372,2.09898276099542,2.115644783157695,2.149272912557035,2.1662868937909714,2.0660513546505883,2.049779911668458,2.235602848865901,2.2889166804078145,2.380596053538979,2.47573150572491,2.47573150572491,2.615193824659969,2.656414604730997,2.7407557944954624,2.805749111310745,2.9174930594233572,2.98654717428537,3.203741702660958,3.409888265666115,3.54530646023206,3.773266592268887,3.9537405526141987,4.273729358185959,4.443259768096302,4.840253677824623,5.071429377205673,5.6546923942954574,5.971016709192964,6.554715672402531]])
#parton_lumi_qq = TGraph(len(parton_lumi_qq_x), parton_lumi_qq_x, parton_lumi_qq_y)


#def PartonLuminosityRatio(x, initial_state):
#	ratio = -1.
#	if initial_state == "gg":
#		ratio = (1.43368e-21*(x**6)+-2.037e-17*(x**5)+1.14943e-13*(x**4)+-3.31131e-10*(x**3)+5.43977e-07*x*x+-0.000607347*x+0.483842)
#	elif initial_state == "qq":
#		ratio = parton_lumi_qq.Eval(x)
#	else:
#		print "[PartonLuminosityRatio] ERROR : initial_state={} not known.".format(initial_state)
#		sys.exit(1)
#	return ratio

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
# xs_limit = TGraph with cross section (total Z' production) limits
def ConvertXSToGB(xs_limit):
	gb_x = array('d', [])
	gb_y = array('d', [])
	for i in xrange(xs_limit.GetN()):
		this_mass = xs_limit.GetX()[i]
		this_xs = xs_limit.GetY()[i]
		this_gB = gq * 6. * math.sqrt(this_xs / dm_sigma.Eval(this_mass)) 
		gb_x.append(this_mass)
		gb_y.append(this_gB)
	return TGraph(len(gb_x), gb_x, gb_y)

def ConvertXSToGQ(xs_limit):
	gq_x = array('d', [])
	gq_y = array('d', [])
	for i in xrange(xs_limit.GetN()):
		this_mass = xs_limit.GetX()[i]
		this_xs = xs_limit.GetY()[i]
		this_gq = gq * math.sqrt(this_xs / dm_sigma.Eval(this_mass))
		gq_x.append(this_mass)
		gq_y.append(this_gq)
		print "Mass {}: gq = {} * sqrt({} / {}) = {}".format(this_mass, gq, this_xs, dm_sigma.Eval(this_mass), this_gq)
		#gq_nodm = 1./6. * math.sqrt(this_xs / (dm_sigma_times_BRjj.Eval(this_mass)/ZpBranchingRatio(this_mass, selected_decays=["u","d","s","c","b"]))) * gq * 6.
		#print "[debug] Mass {}, gq = {}, gq without DM = {}".format(mass, this_gq, gq_nodm)
	return TGraph(len(gq_x), gq_x, gq_y)

def ConvertGBToXS(gb_limit):
	xs_x = array('d', [])
	xs_y = array('d', [])
	for i in xrange(gb_limit.GetN()):
		this_mass = gb_limit.GetX()[i]
		this_gb = gb_limit.GetY()[i]
		this_xs = (this_gb / (gq * 6.))**2 * dm_sigma.Eval(this_mass)
		xs_x.append(this_mass)
		xs_y.append(this_xs)
	return TGraph(len(xs_x), xs_x, xs_y)

def ConvertGQToXS(gq_limit):
	xs_x = array('d', [])
	xs_y = array('d', [])
	for i in xrange(gq_limit.GetN()):
		this_mass = gq_limit.GetX()[i]
		this_gb = 6. * gq_limit.GetY()[i]
		this_xs = (this_gb / (gq * 6.))**2 * dm_sigma.Eval(this_mass)
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

def width_vector(gq, mZp, mf):
	if 2. * mf > mZp:
		return 0.
	else:
	    return gq**2 / (12 * math.pi * mZp) * (mZp**2 + (2 * mf**2)) * math.sqrt(1. - (4. * mf**2 / mZp**2))

def total_width_vector(gq, gdm, mZp, mdm):
	width = 0.
	for q in ["u", "d", "s", "c", "b", "t"]:
		width += 3 * width_vector(gq, mZp, quark_masses[q])

	if gdm > 0:
		width += width_vector(gdm, mZp, mdm)

	return width

def width_axialvector(gq, mZp, mf):
	if 2. * mf > mZp:
		return 0.
	else:
	    return gq**2 / (12 * math.pi * mZp) * (mZp**2 + (2 * mf**2))**1.5

def total_width_axialvector(gq, gdm, mZp, mdm):
	width = 0.
	for q in ["u", "d", "s", "c", "b", "t"]:
		width += 3 * width_axialvector(gq, mZp, quark_masses[q])

	if gdm > 0:
		width += width_axialvector(gdm, mZp, mdm)

	return width

#def avtotwidth(iType,gsm,med):
#	u=avwidth(iType,gsm,med,0.001)
#	d=u
#	s=avwidth(iType,gsm,med,0.135)
#	c=avwidth(iType,gsm,med,1.5)
#	b=avwidth(iType,gsm,med,5.1)
#	t=0
#	if med > 2.*172.5:
#		t=avwidth(iType,gsm,med,172.5)
#	quarks=3*(u+d+s+c+b+t)
#	#print u,d,s,c,b,t,dm,quarks
#	return quarks

def Zconstraints_gq(x, what="gB"):
	
	# Zwidth = 2.4952;
	# ZwidthError = 0.0023*3; # times 3 to give 3 sigma
	# relZwidthUnc = ZwidthError/Zwidth;
	# sin2thetaW = 0.2312;
	# sinthetaW = math.sqrt(sin2thetaW);
	# costhetaW = math.sqrt(1 - sin2thetaW);
	# mW = 80.385;
	# vev = 246.;
	# g = mW*2./vev;
	# Vu = 0.25 - (4. * sin2thetaW / 6.);
	# Vd = -0.25 - (2. * sin2thetaW / 6.);
	# mZ = 91.18;
	# mZp = x[0];

	# # y = gZ
	# ynum = relZwidthUnc * 3 * g * -1. * math.fabs(1-(mZp*mZp/(mZ*mZ))) * (2*Vu*Vu + 3*Vd*Vd + 5/16.)
	# yden = 2*0.01*costhetaW*sinthetaW*(2*Vu+3*Vd);
	# # print ynum,yden,x[0],math.sqrt(ynum/yden)
	# y = math.sqrt(ynum/yden);
	# y *= 1.5;

	mZ = 91.18;
	mZp = x[0];
	ynum = 4. * math.sqrt( 4. * math.pi ) * 1.96 * 1.1e-3 * ( 1-(mZp*mZp/(mZ*mZ)) );
	yden = 1.193 * 0.02;
	if ynum < 0: ynum *= -1.;
	y = math.sqrt(ynum/yden);

	return y / 6.

def Zconstraints_gB(x, what="gB"):
	
	# Zwidth = 2.4952;
	# ZwidthError = 0.0023*3; # times 3 to give 3 sigma
	# relZwidthUnc = ZwidthError/Zwidth;
	# sin2thetaW = 0.2312;
	# sinthetaW = math.sqrt(sin2thetaW);
	# costhetaW = math.sqrt(1 - sin2thetaW);
	# mW = 80.385;
	# vev = 246.;
	# g = mW*2./vev;
	# Vu = 0.25 - (4. * sin2thetaW / 6.);
	# Vd = -0.25 - (2. * sin2thetaW / 6.);
	# mZ = 91.18;
	# mZp = x[0];

	# # y = gZ
	# ynum = relZwidthUnc * 3 * g * -1. * math.fabs(1-(mZp*mZp/(mZ*mZ))) * (2*Vu*Vu + 3*Vd*Vd + 5/16.)
	# yden = 2*0.01*costhetaW*sinthetaW*(2*Vu+3*Vd);
	# # print ynum,yden,x[0],math.sqrt(ynum/yden)
	# y = math.sqrt(ynum/yden);
	# y *= 1.5;

	mZ = 91.18;
	mZp = x[0];
	ynum = 4. * math.sqrt( 4. * math.pi ) * 1.96 * 1.1e-3 * ( 1-(mZp*mZp/(mZ*mZ)) );
	yden = 1.193 * 0.02;
	if ynum < 0: ynum *= -1.;
	y = math.sqrt(ynum/yden);

	return y

# main_limit_name: to emulate the style of EXO-17-001, the results from this analysis have their own legend, and the external results a separate legend.
#def MakeLimitPlot(limit_names, limit_graphs, save_tag, brazil_graphs=None, main_limit_name=None, what="gB", logx=False, logy=False, line_styles={}, line_colors={}, x_range=None, x_title=None, y_title=None, cms_label=None, z_constraint=False):

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description = 'Dijet mass spectrum fits')
	parser.add_argument('--bb_comparison', action='store_true', help='Compare 8 TeV and 13 TeV bb analyses')
	parser.add_argument('--gB', action='store_true', help='Compare all analyses with gB')
	parser.add_argument('--gq', action='store_true', help='Compare all analyses with gq')

	# Plot options
	parser.add_argument('--logx', action='store_true', help='')
	parser.add_argument('--logy', action='store_true', help='')
	parser.add_argument('--x_range', type=float, nargs=2, help='')
	parser.add_argument('--y_range', type=float, nargs=2, help='')
	parser.add_argument('--cms', action="store_true", help='Draw CMS label')
	parser.add_argument('--cms_txt', type=str, help='Draw CMS label with text')
	parser.add_argument('--theory', type=float, help='Draw a theory curve, i.e. flat line at gq or gB value')

	args = parser.parse_args()

	limit_names = []
	limit_collection = {}

	# Load external data
	gr_UA2 = csvToGraph(analysis_config.external_limits + "/UA2.csv",4,True);
	limit_collection["UA2"] = LimitGraph()
	limit_collection["UA2"].LoadGB(gr_UA2)
	gr_CDFRun1 = csvToGraph(analysis_config.external_limits + "/CDF_Run1.csv",2,True );
	limit_collection["CDF Run1"] = LimitGraph()
	limit_collection["CDF Run1"].LoadGB(gr_CDFRun1)
	gr_CDFRun2 = csvToGraph(analysis_config.external_limits + "/CDF_Run2.csv",6,True );
	limit_collection["CDF Run2"] = LimitGraph()
	limit_collection["CDF Run2"].LoadGB(gr_CDFRun2)
	gr_ATLAS = csvToGraph(analysis_config.external_limits + "/gBMZB_ATLAS_all_fbinv.csv",7,False );
	limit_collection["ATLAS dijet, 8 TeV"] = LimitGraph()
	limit_collection["ATLAS dijet, 8 TeV"].LoadGB(gr_ATLAS)
	gr_CMS_2012 = csvToGraph(analysis_config.external_limits + "/CMS_Scouting.csv",8,False );
	limit_collection["CMS dijet, 8 TeV"] = LimitGraph()
	limit_collection["CMS dijet, 8 TeV"].LoadGB(gr_CMS_2012)
	gr_CMS_2015 = csvToGraph(analysis_config.external_limits + "/CMS_2016combined.csv",8,False );
	limit_collection["CMS dijet, 13 TeV"] = LimitGraph()
	limit_collection["CMS dijet, 13 TeV"].LoadGB(gr_CMS_2015)
	gr_CMSBoostedDijet = csvToGraph(analysis_config.external_limits + "/BoostedDijet_final.csv")
	limit_collection["CMS boosted dijet, 13 TeV"] = LimitGraph()
	limit_collection["CMS boosted dijet, 13 TeV"].LoadGB(gr_CMSBoostedDijet)

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
		"CDF Run1":2,
		"CDF Run2":2,
		"ATLAS dijet, 8 TeV":6,
		"CMS boosted dijet, 13 TeV":7,
		"CMS dijet, 8 TeV":5,
		"CMS dijet, 13 TeV":8,
		"CMS bb, 2015 exp":3,
		"CMS bbl, 2012 exp":3,
		"CMS bbh, 2012 exp":3,
		"CMS bb, 2015 obs":1,
		"CMS bbh, 2012 obs":1,
	}
	line_widths = {
		"CMS bb, 2012 exp":2,
		"CMS bb, 2012 obs":2,
	}
	line_colors = {
		"UA2":seaborn.GetColorRoot("cubehelixlarge", 1, 15),
		"CDF Run1":seaborn.GetColorRoot("cubehelixlarge", 6, 15),
		"CDF Run2":seaborn.GetColorRoot("cubehelixlarge", 11, 15),
		"ATLAS dijet, 8 TeV":seaborn.GetColorRoot("cubehelixlarge", 2, 15),
		"CMS dijet, 8 TeV":seaborn.GetColorRoot("cubehelixlarge", 8, 15),
		"CMS dijet, 13 TeV":seaborn.GetColorRoot("cubehelixlarge", 13, 15),
		"CMS boosted dijet, 13 TeV":seaborn.GetColorRoot("cubehelixlarge", 7, 15),
		"CMS bb, 2015 exp":seaborn.GetColorRoot("default", 1),
		"CMS bbl, 2015 obs":seaborn.GetColorRoot("default", 1),
		"CMS bbl, 2012 exp":1,
		"CMS bbl, 2012 obs":1,
		"CMS bbh, 2012 exp":1,
		"CMS bbh, 2012 obs":1,
	}

	if args.bb_comparison:
		# 8 TeV bb
		f_Xbb_8TeVl = TFile("~/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_trigbbl_CSVTM_ZPrime_dijet4.root", "READ")
		gr_Xbb_8TeVl_exp = f_Xbb_8TeVl.Get("graph_exp")
		BR_bb_8TeVl = []
		for i in xrange(gr_Xbb_8TeVl_exp.GetN()):
			mass = gr_Xbb_8TeVl_exp.GetX()[i]
			gr_Xbb_8TeVl_exp.SetPoint(i, mass, gr_Xbb_8TeVl_exp.GetY()[i]) # / ZpBranchingRatio(mass, selected_decays=["b"]) 
			BR_bb_8TeVl.append(ZpBranchingRatio(mass, selected_decays=["c", "b"]) )
		limit_collection["CMS bb, 2012 low exp"] = LimitGraph()
		limit_collection["CMS bb, 2012 low exp"].LoadXSBR(gr_Xbb_8TeVl_exp, BR_bb_8TeVl)

		gr_Xbb_8TeVl_obs = f_Xbb_8TeVl.Get("graph_obs")
		for i in xrange(gr_Xbb_8TeVl_obs.GetN()):
			mass = gr_Xbb_8TeVl_obs.GetX()[i]
			gr_Xbb_8TeVl_obs.SetPoint(i, mass, gr_Xbb_8TeVl_obs.GetY()[i])
		limit_collection["CMS bb, 2012 low obs"] = LimitGraph()
		limit_collection["CMS bb, 2012 low obs"].LoadXSBR(gr_Xbb_8TeVl_obs, BR_bb_8TeVl)

		f_Xbb_8TeVh = TFile("~/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_trigbbh_CSVTM_ZPrime_dijet4.root", "READ")
		gr_Xbb_8TeVh_exp = f_Xbb_8TeVh.Get("graph_exp")
		BR_bb_8TeVh = []
		for i in xrange(gr_Xbb_8TeVh_exp.GetN()):
			mass = gr_Xbb_8TeVh_exp.GetX()[i]
			gr_Xbb_8TeVh_exp.SetPoint(i, mass, gr_Xbb_8TeVh_exp.GetY()[i])
			BR_bb_8TeVh.append(ZpBranchingRatio(mass, selected_decays=["c", "b"]) )
		limit_collection["CMS bb, 2012 high exp"] = LimitGraph()
		limit_collection["CMS bb, 2012 high exp"].LoadXSBR(gr_Xbb_8TeVh_exp, BR_bb_8TeVh)

		gr_Xbb_8TeVh_obs = f_Xbb_8TeVh.Get("graph_obs")
		for i in xrange(gr_Xbb_8TeVh_obs.GetN()):
			mass = gr_Xbb_8TeVh_obs.GetX()[i]
			gr_Xbb_8TeVh_obs.SetPoint(i, mass, gr_Xbb_8TeVh_obs.GetY()[i])
			print "[debug] Set point {} {}".format(mass, gr_Xbb_8TeVh_obs.GetY()[i])
		limit_collection["CMS bb, 2012 high obs"] = LimitGraph()
		limit_collection["CMS bb, 2012 high obs"].LoadXSBR(gr_Xbb_8TeVh_obs, BR_bb_8TeVh)

		limit_names = ["CMS bb, 2015 exp", "CMS bb, 2015 obs","CMS bb, 2012 low exp", "CMS bb, 2012 low obs","CMS bb, 2012 high exp", "CMS bb, 2012 high obs"]
		#MakeLimitPlot(limit_names, limit_collection, "xs_combination_log", what="xsBR", x_title="m_{X} [GeV]", y_title="#sigma_{13 TeV} #times BR(b#bar{b}) [pb]", logy=True, line_styles=line_styles, line_colors=line_colors, x_range=[0,1500])

	elif args.gB:
		# 8 TeV bb
		f_Xbb_8TeVl = TFile("~/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_trigbbl_CSVTM_ZPrime_dijet4.root", "READ")
		gr_Xbb_8TeVl_exp = f_Xbb_8TeVl.Get("graph_exp")
		BR_bb_8TeVl = []
		for i in xrange(gr_Xbb_8TeVl_exp.GetN()):
			mass = gr_Xbb_8TeVl_exp.GetX()[i]
			gr_Xbb_8TeVl_exp.SetPoint(i, mass, gr_Xbb_8TeVl_exp.GetY()[i]) # / ZpBranchingRatio(mass, selected_decays=["b"]) 
			BR_bb_8TeVl.append(ZpBranchingRatio(mass, selected_decays=["b"]))
			#print "[debug] BR(Z-->bb, M={} GeV) = {}".format(mass, BR_bb_8TeVl[-1])
		limit_collection["CMS bb, 2012 low exp"] = LimitGraph()
		limit_collection["CMS bb, 2012 low exp"].LoadXSBR(gr_Xbb_8TeVl_exp, BR_bb_8TeVl)

		gr_Xbb_8TeVl_obs = f_Xbb_8TeVl.Get("graph_obs")
		for i in xrange(gr_Xbb_8TeVl_obs.GetN()):
			mass = gr_Xbb_8TeVl_obs.GetX()[i]
			gr_Xbb_8TeVl_obs.SetPoint(i, mass, gr_Xbb_8TeVl_obs.GetY()[i])
		limit_collection["CMS bb, 2012 low obs"] = LimitGraph()
		limit_collection["CMS bb, 2012 low obs"].LoadXSBR(gr_Xbb_8TeVl_obs, BR_bb_8TeVl)

		f_Xbb_8TeVh = TFile("~/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_trigbbh_CSVTM_ZPrime_dijet4.root", "READ")
		gr_Xbb_8TeVh_exp = f_Xbb_8TeVh.Get("graph_exp")
		BR_bb_8TeVh = []
		for i in xrange(gr_Xbb_8TeVh_exp.GetN()):
			mass = gr_Xbb_8TeVh_exp.GetX()[i]
			gr_Xbb_8TeVh_exp.SetPoint(i, mass, gr_Xbb_8TeVh_exp.GetY()[i])
			BR_bb_8TeVh.append(ZpBranchingRatio(mass, selected_decays=["b"]) )
		limit_collection["CMS bb, 2012 high exp"] = LimitGraph()
		limit_collection["CMS bb, 2012 high exp"].LoadXSBR(gr_Xbb_8TeVh_exp, BR_bb_8TeVh)

		gr_Xbb_8TeVh_obs = f_Xbb_8TeVh.Get("graph_obs")
		for i in xrange(gr_Xbb_8TeVh_obs.GetN()):
			mass = gr_Xbb_8TeVh_obs.GetX()[i]
			gr_Xbb_8TeVh_obs.SetPoint(i, mass, gr_Xbb_8TeVh_obs.GetY()[i])
		limit_collection["CMS bb, 2012 high obs"] = LimitGraph()
		limit_collection["CMS bb, 2012 high obs"].LoadXSBR(gr_Xbb_8TeVh_obs, BR_bb_8TeVh)

		limit_names = ["UA2","CDF Run1","CDF Run2","ATLAS dijet, 8 TeV","CMS dijet, 8 TeV", "CMS boosted dijet, 13 TeV", "CMS bb, 2015 obs","CMS bb, 2012 low obs", "CMS bb, 2012 low exp","CMS bb, 2012 high obs", "CMS bb, 2012 high exp"]

		#MakeLimitPlot(limit_names, limit_collection, "gB_combination_log", what="gB", logy=True, line_styles=line_styles)
		#MakeLimitPlot(limit_names, limit_collection, "xs_combination", what="xs", line_styles=line_styles, line_colors=line_colors, x_range=[0,1500])
		#MakeLimitPlot(limit_names, limit_collection, "gB_combination", what="gB", line_styles=line_styles, line_colors=line_colors, x_title="m_{Z'} [GeV]", y_title="g_{B}", x_range=[0,1000])
		#MakeLimitPlot(limit_names, limit_collection, "gB_combination_log", what="gB", line_styles=line_styles, line_colors=line_colors, x_title="m_{Z'} [GeV]", y_title="g_{B}", x_range=[0,1000], logy=True)

	elif args.gq:
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

		limit_names = ["UA2", "CDF Run1", "CDF Run2", "ATLAS dijet, 8 TeV", "CMS dijet, 8 TeV", "CMS dijet, 13 TeV", "CMS boosted dijet, 13 TeV", "CMS bbl, 2012 obs", "CMS bbl, 2012 exp", "CMS bbh, 2012 obs", "CMS bbh, 2012 exp"]
	what = "gq"
	if not what in ["gB", "gq", "xs", "xsBR"]:
		print "[MakeLimitPlot] ERROR : what must be gB, gq, xs, or xsBR"
		sys.exit(1)

	cname = "c_{}_combination".format(what)
	if args.logx:
		cname += "_logx"
	if args.logy:
		cname += "_logy"
	if args.theory:
		cname += "_plustheory"
	c = TCanvas(cname, cname, 1000,800)
	if args.logx:
		c.SetLogx()
	if args.logy:
		c.SetLogy()

	if args.x_range:
		x_min = args.x_range[0]
		x_max = args.x_range[1]
	else:
		x_min = 0.
		x_max = 1300.
	y_min = 1.e20
	y_max = -1.e20
	for name in limit_names:
		for i in xrange(limit_collection[name].GetGB().GetN()):
			this_x = limit_collection[name].GetGB().GetX()[i]
			if what == "gB":
				this_y = limit_collection[name].GetGB().GetY()[i]
			elif what == "gq":
				this_y = limit_collection[name].GetGQ().GetY()[i]
			elif what == "xs":
				this_y = limit_collection[name].GetXS().GetY()[i]
			elif what == "xsBR":
				this_y = limit_collection[name].GetXSBR().GetY()[i]
			if this_y > y_max:
				y_max = this_y
			if this_y < y_min:
				y_min = this_y

	# Start plotting
	if what == "gB":
		if args.logy:
			y_min = 0.1
			y_max = 10.
		else:
			y_min = 0.
			y_max = 5.
	elif what == "gq":
		if args.logy:
			y_min = 0.02
			y_max = 1.
		else:
			y_min = 0.
			y_max = 1.
	else:
		if args.logy:
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
	frame.GetXaxis().SetTitle("m_{Z'} [GeV]")
	
	if what == "gB":
		frame.GetYaxis().SetTitle("g_{B}")
	elif what == "gq":
		frame.GetYaxis().SetTitle("g_{q}")
	elif what == "xs":
		frame.GetYaxis().SetTitle("#sigma [pb]")
	elif what == "xsBR":
		frame.GetYaxis().SetTItle("#sigma #times BR [pb]")
	print "[debug] frame min/max = {}/{}".format(frame.GetXaxis().GetXmin(), frame.GetXaxis().GetXmax())
	frame.Draw("axis")

	# Draw brazil plots first
	if what == "gq":
		limit_collection["CMS bbl, 2012 exp2sig"].GetGQ().SetFillColor(ROOT.kYellow)
		limit_collection["CMS bbl, 2012 exp2sig"].GetGQ().SetMarkerStyle(20)
		limit_collection["CMS bbl, 2012 exp2sig"].GetGQ().SetMarkerSize(0)
		limit_collection["CMS bbl, 2012 exp2sig"].GetGQ().SetLineWidth(0)
		limit_collection["CMS bbl, 2012 exp2sig"].GetGQ().Draw("f")

		limit_collection["CMS bbl, 2012 exp1sig"].GetGQ().SetFillColor(ROOT.kGreen)
		limit_collection["CMS bbl, 2012 exp1sig"].GetGQ().SetMarkerStyle(20)
		limit_collection["CMS bbl, 2012 exp1sig"].GetGQ().SetMarkerSize(0)
		limit_collection["CMS bbl, 2012 exp1sig"].GetGQ().SetLineWidth(0)
		limit_collection["CMS bbl, 2012 exp1sig"].GetGQ().Draw("f")

		limit_collection["CMS bbh, 2012 exp2sig"].GetGQ().SetFillColor(ROOT.kYellow)
		limit_collection["CMS bbh, 2012 exp2sig"].GetGQ().SetMarkerStyle(20)
		limit_collection["CMS bbh, 2012 exp2sig"].GetGQ().SetMarkerSize(0)
		limit_collection["CMS bbh, 2012 exp2sig"].GetGQ().SetLineWidth(0)
		limit_collection["CMS bbh, 2012 exp2sig"].GetGQ().Draw("f")

		print "[debug] Printing 2012 bbh exp2sig graph"
		limit_collection["CMS bbh, 2012 exp2sig"].GetGB().Print("all")

		limit_collection["CMS bbh, 2012 exp1sig"].GetGQ().SetFillColor(ROOT.kGreen)
		limit_collection["CMS bbh, 2012 exp1sig"].GetGQ().SetMarkerStyle(20)
		limit_collection["CMS bbh, 2012 exp1sig"].GetGQ().SetMarkerSize(0)
		limit_collection["CMS bbh, 2012 exp1sig"].GetGQ().SetLineWidth(0)
		limit_collection["CMS bbh, 2012 exp1sig"].GetGQ().Draw("f")
	elif what == "gB":
		limit_collection["CMS bbl, 2012 exp2sig"].GetGB().SetFillColor(ROOT.kYellow)
		limit_collection["CMS bbl, 2012 exp2sig"].GetGB().SetMarkerStyle(20)
		limit_collection["CMS bbl, 2012 exp2sig"].GetGB().SetMarkerSize(0)
		limit_collection["CMS bbl, 2012 exp2sig"].GetGB().SetLineWidth(0)
		limit_collection["CMS bbl, 2012 exp2sig"].GetGB().Draw("f")

		limit_collection["CMS bbl, 2012 exp1sig"].GetGB().SetFillColor(ROOT.kGreen)
		limit_collection["CMS bbl, 2012 exp1sig"].GetGB().SetMarkerStyle(20)
		limit_collection["CMS bbl, 2012 exp1sig"].GetGB().SetMarkerSize(0)
		limit_collection["CMS bbl, 2012 exp1sig"].GetGB().SetLineWidth(0)
		limit_collection["CMS bbl, 2012 exp1sig"].GetGB().Draw("f")

		print "[debug] Printing 2012 bbh exp2sig graph"
		limit_collection["CMS bbh, 2012 exp2sig"].GetGB().Print("all")
		limit_collection["CMS bbh, 2012 exp2sig"].GetGB().SetFillColor(ROOT.kYellow)
		limit_collection["CMS bbh, 2012 exp2sig"].GetGB().SetMarkerStyle(20)
		limit_collection["CMS bbh, 2012 exp2sig"].GetGB().SetMarfx_rangekerSize(0)
		limit_collection["CMS bbh, 2012 exp2sig"].GetGB().SetLineWidth(0)
		limit_collection["CMS bbh, 2012 exp2sig"].GetGB().Draw("f")

		limit_collection["CMS bbh, 2012 exp1sig"].GetGB().SetFillColor(ROOT.kGreen)
		limit_collection["CMS bbh, 2012 exp1sig"].GetGB().SetMarkerStyle(20)
		limit_collection["CMS bbh, 2012 exp1sig"].GetGB().SetMarkerSize(0)
		limit_collection["CMS bbh, 2012 exp1sig"].GetGB().SetLineWidth(0)
		limit_collection["CMS bbh, 2012 exp1sig"].GetGB().Draw("f")

	# Z constraint next
	if what in ["gB", "gq"]:
		if what == "gB":
			z_constraint = TF1('myfunc',Zconstraints_gB,0,1000,0)
		elif what == "gq":
			z_constraint = TF1('myfunc',Zconstraints_gq,0,1000,0)
		z_constraint.SetNpx(1000)
		z_constraint.SetLineColor(15)
		z_constraint.SetLineStyle(9)
		z_constraint.SetLineWidth(2)
		z_constraint.Draw("same")

	# Limit curves
	style_counter = 0
	for name in limit_names:
		if "CMS bb" in name and "2012" in name:
			draw_style = "l"
		else:
			draw_style = "c"
		if what == "gB":
			if name in line_colors:
				limit_collection[name].GetGB().SetLineColor(line_colors[name])
			else:
				limit_collection[name].GetGB().SetLineColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(limit_names)))
			if name in line_widths:
				limit_collection[name].GetGB().SetLineWidth(line_widths[name])
			else:
				limit_collection[name].GetGB().SetLineWidth(2)
			if name in line_styles:
				limit_collection[name].GetGB().SetLineStyle(line_styles[name])
			limit_collection[name].GetGB().Draw(draw_style)
		elif what == "gq":
			if name in line_colors:
				limit_collection[name].GetGQ().SetLineColor(line_colors[name])
			else:
				limit_collection[name].GetGQ().SetLineColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(limit_names)))
			if name in line_widths:
				limit_collection[name].GetGQ().SetLineWidth(line_widths[name])
			else:
				limit_collection[name].GetGQ().SetLineWidth(2)
			if name in line_styles:
				limit_collection[name].GetGQ().SetLineStyle(line_styles[name])
			limit_collection[name].GetGQ().Draw(draw_style)
			#limit_collection[name].GetGQ().Print("all")
		elif what == "xs":
			if name in line_colors:
				limit_collection[name].GetXS().SetLineColor(line_colors[name])
			else:
				limit_collection[name].GetXS().SetLineColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(limit_names)))
			limit_collection[name].GetXS().SetLineWidth(2)
			if name in line_styles:
				limit_collection[name].GetXS().SetLineStyle(line_styles[name])
			limit_collection[name].GetXS().Draw("l")
		elif what == "xsBR":
			if name in line_colors:
				limit_collection[name].GetXSBR().SetLineColor(line_colors[name])
			else:
				limit_collection[name].GetXSBR().SetLineColor(seaborn.GetColorRoot("cubehelixlarge", style_counter, len(limit_names)))
			limit_collection[name].GetXSBR().SetLineWidth(2)
			if name in line_styles:
				limit_collection[name].GetXSBR().SetLineStyle(line_styles[name])
			limit_collection[name].GetXSBR().Draw("l")
		style_counter += 1

	# Main legend
	if args.logy:
		l1 = TLegend(0.49, 0.71, 0.67, 0.87)
		l1.SetFillColor(0)
		l1.SetBorderSize(0)
		l1.SetHeader("95% CL upper limits")
		if what == "gB":
			l1.AddEntry(limit_collection["CMS bbl, 2012 obs"].GetGB(), "Observed", "l")
			l1.AddEntry(limit_collection["CMS bbl, 2012 exp"].GetGB(), "Expected", "l")
			l1.AddEntry(limit_collection["CMS bbl, 2012 exp1sig"].GetGB(), "Exp. #pm 1#sigma", "f")
			l1.AddEntry(limit_collection["CMS bbl, 2012 exp2sig"].GetGB(), "Exp. #pm 2#sigma", "f")
		elif what == "gq":
			l1.AddEntry(limit_collection["CMS bbl, 2012 obs"].GetGQ(), "Observed", "l")
			l1.AddEntry(limit_collection["CMS bbl, 2012 exp"].GetGQ(), "Expected", "l")
			l1.AddEntry(limit_collection["CMS bbl, 2012 exp1sig"].GetGQ(), "Exp. #pm 1#sigma", "f")
			l1.AddEntry(limit_collection["CMS bbl, 2012 exp2sig"].GetGQ(), "Exp. #pm 2#sigma", "f")
		elif what == "xs":
			l1.AddEntry(limit_collection["CMS bbl, 2012 obs"].GetXS(), "Observed", "l")
			l1.AddEntry(limit_collection["CMS bbl, 2012 exp"].GetXS(), "Expected", "l")
			l1.AddEntry(limit_collection["CMS bbl, 2012 exp1sig"].GetXS(), "Exp. #pm 1#sigma", "f")
			l1.AddEntry(limit_collection["CMS bbl, 2012 exp2sig"].GetXS(), "Exp. #pm 2#sigma", "f")
		elif what == "xsBR":
			l1.AddEntry(limit_collection["CMS bbl, 2012 obs"].GetXSBR(), "Observed", "l")
			l1.AddEntry(limit_collection["CMS bbl, 2012 exp"].GetXSBR(), "Expected", "l")
			l1.AddEntry(limit_collection["CMS bbl, 2012 exp1sig"].GetXSBR(), "Exp. #pm 1#sigma", "f")
			l1.AddEntry(limit_collection["CMS bbl, 2012 exp2sig"].GetXSBR(), "Exp. #pm 2#sigma", "f")
		l1.Draw()

		l2 = TLegend(0.4, 0.2, 0.92, 0.36)
		l2.SetNColumns(2)
		l2.SetFillColor(0)
		l2.SetBorderSize(0)
		external_legend_entries = {
			"UA2":"UA2 [2]",
			"CMS dijet, 8 TeV":"CMS dijet, 8 TeV [33]",
			"CDF Run1":"CDF Run1 [8]",
			"CMS dijet, 13 TeV":"CMS dijet, 13 TeV [34]",
			"CDF Run2":"CDF Run2 [9]",
			"CMS boosted dijet, 13 TeV":"CMS boosted dijet, 13 TeV [36]",
		}		
		for name in ["UA2", "CMS dijet, 8 TeV", "CDF Run1", "CMS dijet, 13 TeV", "CDF Run2", "CMS boosted dijet, 13 TeV"]:
			if what == "gq":
				l2.AddEntry(limit_collection[name].GetGQ(), external_legend_entries[name], "l")
			elif what == "gq":
				l2.AddEntry(limit_collection[name].GetGB(), external_legend_entries[name], "l")
		l2.AddEntry(z_constraint, "Z width (indirect) [70]", "l")
		if what == "gq":
			l2.AddEntry(limit_collection["ATLAS dijet, 8 TeV"].GetGQ(), "ATLAS dijet, 8 TeV [18]", "l")
		elif what == "gq":
			l2.AddEntry(limit_collection["ATLAS dijet, 8 TeV"].GetGB(), "ATLAS dijet, 8 TeV [18]", "l")
		l2.Draw()

	if args.cms:
		Root.CMSLabel(0.49, 0.88, "", 1, 0.5);
	elif args.cms_txt:
		Root.CMSLabel(0.49, 0.88, args.cms_txt, 1, 0.5);
	Root.myText(0.7, 0.96, 1, "19.7 fb^{-1} (8 TeV)", 0.5)

	if args.theory:
		theory_curve = TLine(frame.GetXaxis().GetXmin(), args.theory, frame.GetXaxis().GetXmax(), args.theory)
		theory_curve.SetLineStyle(2)
		theory_curve.SetLineWidth(3)
		theory_curve.SetLineColor(kRed)
		theory_curve.Draw("same")

	c.SaveAs("{}/{}.pdf".format(analysis_config.figure_directory, c.GetName()))
