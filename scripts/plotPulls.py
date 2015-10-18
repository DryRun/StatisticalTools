#! /usr/bin/env python

import os
import sys
import argparse
import math
from ROOT import *
from setTDRStyle import setTDRStyle
from array import array



usage = "usage: python plotPulls.py -i <inputlist> -o <outputdir> --mu <mu> (-999. if you wat to inject the signal corresponding to the 2sigma band of the upper limit)"
print usage

parser = argparse.ArgumentParser(description='Process options.')

parser.add_argument("--inputList", type=str, dest="inputList", default="",
    help="list of datacard processed",
    )
parser.add_argument("-o", "--output", type=str, dest="output", default="",
    help="the directory OUTDIR contains the output of the program",
    metavar="OUTDIR"
    )
parser.add_argument("--mu", type=float, dest="mu", default=0,
    help="injected signal strenght",
    metavar="MU"
    )
parser.add_argument("--inputFitRes", type=str, dest="inputFitRes", default="",
    help="directory containing files with results of fit",
    )
#parser.add_argument("--inputLimits", type=str, dest="inputLimits", default="",
#    help="directory containing results of limits",
#    )
parser.add_argument("--tag", type=str, dest="tag", default="",
    help="tag contained in the name of files",
    )

args = parser.parse_args()
print args
#####################################

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')


gStyle.SetOptFit(1111)

pull_mu_list = []
err_pull_mu_list = []
err_mass = []
#giulia --- to change 
mass = []
#mass = [2,3,4,5,6]

ins = open(args.inputList,"r")


list_median = []
list_sigmaD68 = []
list_sigmaU68 = []
list_sigmaD95 = []
list_sigmaU95 = []

k=0
for line in  ins:
  sample = os.path.basename(line)
  #sample = os.path.splitext(line)[0]
  sample = sample.split("_")[0]
  print ("process %s" % sample)
  line = line.rstrip('\n')
  sample = sample.rstrip('\n')
  mass_string = sample.split("Qstar")[1]
  print mass_string
  #if not( mass_string=="2500"): 
  #  continue
  mass.append(int(mass_string))
 
 #extract mu limit 2 sigma band 
  #filename = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/dijet_limits_interpolated_pseudodatasetDinko//higgsCombine'+sample+'_limit.Asymptotic.mH120.root'
  filename = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_5/src/StatisticalTools/scripts/output_limits_dataRunD_830pb-1/higgsCombine'+sample+'_limit.Asymptotic.mH120.root'
  inf = TFile.Open(filename)
  tr = inf.Get('limit')
  y = []
  N = tr.GetEntriesFast()
  for i in xrange(N):
    tr.GetEntry(i)
    y.append(tr.limit)
    k+=1
 

  if not(args.output):
    os.system("mkdir "+args.output)

  if not(args.mu==-999.):
    inputFile = TFile(args.inputFitRes+"/mlfit"+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+".root")
    #inputFile = TFile("output_toys_mu"+str(int(args.mu))+"_flatParam/mlfit"+sample+"_MLfit_mu"+str(int(args.mu))+".root","read") 
    #inputFile = TFile("test/mlfit"+sample+"_MLfit_mu"+str(int(args.mu))+".root","read") 
    #inputFile = TFile("output_toys_mu"+str(int(args.mu))+"_flatParam_genFromHist/mlfit"+sample+"_MLfit_mu"+str(int(args.mu))+".root","read") 
    print inputFile
  else:
    inputFile = TFile(args.inputFitRes+"/mlfit"+sample+args.tag+"_MLfit_mu_limit.root")
    #inputFile = TFile("output_toys_mu_limit_2Sband_flatParam/mlfit"+sample+"_MLfit_mu_limit.root","read") 
    #inputFile = TFile("output_toys_mu_limit_flatParam_genFromHist/mlfit"+sample+"_MLfit_mu_limit.root","read") 
   
  
  tree_fit_sb = inputFile.Get("tree_fit_sb")
  h_pull_mu = TH1F("h_pull_mu", "S_{true} - S_{fit} / err", 100, -3., 3.)
  if not (args.mu == -999):
    mu = args.mu
  else:
    mu = y[4]
      
  #giulia test --> use always high error
  print("WARNING : using muHiErr for both low and high band")

  for event in tree_fit_sb:
    if tree_fit_sb.fit_status==0:
      if tree_fit_sb.mu > mu:
        #err_mu = tree_fit_sb.muLoErr
        err_mu = tree_fit_sb.muHiErr
      else:
        err_mu = tree_fit_sb.muHiErr
  
      if err_mu ==0:
        pull_mu = 0
      else:  
        pull_mu = (tree_fit_sb.mu - mu)/err_mu
      
      if(abs(pull_mu)>10 ):
        print "mu = " + str(mu) 
        print "mu_fit - mu_exp  = "+str(tree_fit_sb.mu - mu )
        print "err_mu = "+str(err_mu)
        print "pull_mu = "+str(pull_mu)
      if(abs(pull_mu) < 0.01):
        print "mu = " + str(mu) 
        print "mu_fit - mu_exp  = "+str(tree_fit_sb.mu - mu )
        print "err_mu = "+str(err_mu)
        print "pull_mu = "+str(pull_mu)
         

      h_pull_mu.Fill( pull_mu )
  
  h_pull_mu.Print() 

  ##################################
  #### fit pulls ####
  mygaus = TF1("mygaus",  "[0]*TMath::Exp(-0.5*TMath::Power(((x-[1])/[2]),2))", -5., 5.)
  mygaus.SetParameter(0,10)
  mygaus.SetParameter(1,0.)
  mygaus.SetParameter(2,1.)
  h_pull_mu.Fit("mygaus")

  p0 = mygaus.GetParameter(0)
  p1 = mygaus.GetParameter(1)
  p2 = mygaus.GetParameter(2)
  e0 = mygaus.GetParError(0)
  e1 = mygaus.GetParError(1)
  e2 = mygaus.GetParError(2)
 
  ## append pulls and error to a list
  pull_mu_list.append(p1)
  err_pull_mu_list.append(p2)
  err_mass.append(0)

  ##### Draw and save fit to pulls #####
  c = TCanvas("c","",800,600)
  h_pull_mu.Draw()
  if not(args.mu==-999.):
    c.SaveAs(args.output+"/pull_"+sample+"_mu_"+str(int(args.mu))+".png")
    c.SaveAs(args.output+"/pull_"+sample+"_mu_"+str(int(args.mu))+".pdf")
  else:
    c.SaveAs(args.output+"/pull_"+sample+"_mu_limit_2Sband.png")
    c.SaveAs(args.output+"/pull_"+sample+"_mu_limit_2Sband.pdf")

  c.Clear()

  ###########################################
  ### compute quantiles
  nqm = 1
  xqm_= [0.5]
  yqm_= [0]
  xqm=array("d",xqm_)
  yqm=array("d",yqm_)
  h_pull_mu.GetQuantiles(nqm,yqm,xqm)
  median =yqm[0];
  print "**** median **** " + str(median)
  list_median.append(median)

  #68% bands
  nq68=3
  xq68_= []
  yq68_= [0,0,0]
  xq68_.append(0.5-0.68/2.)
  xq68_.append(0.5)
  xq68_.append(0.5+0.68/2.)
  xq68=array("d",xq68_)
  yq68=array("d",yq68_)
  h_pull_mu.GetQuantiles(nq68,yq68,xq68)
  mean68=yq68[1]
  sigmaD68=TMath.Abs(yq68[1]-yq68[0])
  sigmaU68=TMath.Abs(yq68[2]-yq68[1])
  list_sigmaD68.append(sigmaD68)
  list_sigmaU68.append(sigmaU68)
  
  #95% bands
  nq95=3
  xq95_ =[]
  yq95_ =[0,0,0]
  xq95_.append(0.5-0.95/2.)
  xq95_.append(0.5)
  xq95_.append(0.5+0.95/2)
  xq95=array("d",xq95_)
  yq95=array("d",yq95_)
  h_pull_mu.GetQuantiles(nq95,yq95,xq95)
  mean95=yq95[1]
  sigmaD95=TMath.Abs(yq95[1]-yq95[0])
  sigmaU95=TMath.Abs(yq95[2]-yq95[1])
  list_sigmaD95.append(sigmaD95)
  list_sigmaU95.append(sigmaU95)

  #LOOP OVER MASSES finished


### Draw bias plot ####
# with fit
v_mass = array("d",mass)
v_err_mass = array("d",err_mass)
v_pull_mu = array("d",pull_mu_list)
v_pull_mu_err = array("d", err_pull_mu_list)

g_pull_mu = TGraphErrors(len(v_mass),v_mass,v_pull_mu,v_err_mass, v_pull_mu_err)
g_pull_mu_simple = TGraph(len(v_mass),v_mass,v_pull_mu )

#with quantiles
######
v_median = array("d",list_median)
v_sigmaD68 = array("d",list_sigmaD68)
v_sigmaU68 = array("d",list_sigmaU68) 
v_sigmaD95 = array("d",list_sigmaD95) 
v_sigmaU95 = array("d",list_sigmaU95) 

g_pull_mu_median = TGraph(len(v_mass),v_mass,v_median )
g_pull_mu_68 = TGraphAsymmErrors(len(v_mass),v_mass,v_median,v_err_mass,v_err_mass,v_sigmaD68,v_sigmaU68)
g_pull_mu_95 = TGraphAsymmErrors(len(v_mass),v_mass,v_median,v_err_mass,v_err_mass,v_sigmaD95,v_sigmaU95)

#style
g_pull_mu.GetXaxis().SetTitle(' qg resonance Mass (TeV)')
g_pull_mu.GetYaxis().SetTitle('Pull in Mean (Standard Deviations)')
g_pull_mu.GetYaxis().SetNdivisions(510)
g_pull_mu.GetYaxis().SetRangeUser(-1.5,1.5)
g_pull_mu.SetFillStyle(3001)
g_pull_mu.SetFillColor(kAzure+1)
g_pull_mu_simple.SetMarkerStyle(20)
g_pull_mu_simple.SetMarkerSize(0.8)
g_pull_mu_simple.SetMarkerColor(kBlack)

g_pull_mu_95.GetXaxis().SetTitle(' qg resonance Mass (TeV)')
g_pull_mu_95.GetYaxis().SetTitle('Pull median')
g_pull_mu_95.GetYaxis().SetNdivisions(510)
g_pull_mu_95.GetYaxis().SetRangeUser(-3,3)
g_pull_mu_median.SetMarkerStyle(20)
g_pull_mu_median.SetMarkerSize(0.8)
g_pull_mu_median.SetMarkerColor(kBlack)
g_pull_mu_68.SetFillColor(kGreen)
g_pull_mu_68.SetMarkerSize(0)
g_pull_mu_95.SetFillColor(kYellow)
g_pull_mu_95.SetMarkerSize(0)


#xmin = g_pull_mu.GetXaxis().GetXmin()
#xmax = g_pull_mu.GetXaxis().GetXmax()
xmin = mass[0]
xmax = mass[len(mass)-1]

line_0  = TLine(xmin,0.,xmax,0.);
line_up = TLine(xmin,0.5,xmax,0.5);
line_do = TLine(xmin,-0.5,xmax,-0.5);

line_up.SetLineStyle(2)
line_do.SetLineStyle(2)
line_0.SetLineColor(kBlack)
line_up.SetLineColor(kBlack)
line_do.SetLineColor(kBlack)

leg = TLegend(0.5,0.8,0.9,0.9)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.SetLineWidth(0)

leg.AddEntry(g_pull_mu_simple,"mean gaussian fit","p")
leg.AddEntry(g_pull_mu,"1 #sigma","f")

leg2 = TLegend(0.5,0.7,0.9,0.9)
leg2.SetLineColor(0)
leg2.SetFillColor(0)
leg2.SetLineWidth(0)
leg2.AddEntry(g_pull_mu_median,"median","p")
leg2.AddEntry(g_pull_mu_68,"68%","f")
leg2.AddEntry(g_pull_mu_95,"95%","f")



### Draw canvas 
canvasName = 'Bias'
can = TCanvas(canvasName,canvasName,900,600)
if (args.mu ==0 ):
  g_pull_mu.SetTitle("Bias VS resonance mass, no signal")
if (args.mu == -999):
  g_pull_mu.SetTitle("Bias VS resonance mass, 2#sigma signal")

can.cd()
g_pull_mu.Draw("ae3")
g_pull_mu_simple.Draw("lp")
line_0.Draw()
line_up.Draw()
line_do.Draw()
leg.Draw()

if not(args.mu==-999.):
  can.SaveAs(args.output+"/bias_standard_param_mu"+str(int(args.mu))+".png")
  can.SaveAs(args.output+"/bias_standard_param_mu"+str(int(args.mu))+".pdf")
else:
  can.SaveAs(args.output+"/bias_standard_param_mu_limit_2Sband.png")
  can.SaveAs(args.output+"/bias_standard_param_mu_limit_2Sband.pdf")

canvasName = 'Bias median'
can2 = TCanvas(canvasName,canvasName,900,600)
if (args.mu ==0 ):
  g_pull_mu_95.SetTitle("Bias VS resonance mass, no signal")
if (args.mu == -999):
  g_pull_mu_95.SetTitle("Bias VS resonance mass, 2#sigma signal")

can2.cd()
g_pull_mu_95.Draw("ae3")
g_pull_mu_68.Draw("e3")
g_pull_mu_median.Draw("lp")
line_0.Draw()
line_up.Draw()
line_do.Draw()
leg2.Draw()

if not(args.mu==-999.):
  can2.SaveAs(args.output+"/bias_median_standard_param_mu"+str(int(args.mu))+".png")
  can2.SaveAs(args.output+"/bias_median_standard_param_mu"+str(int(args.mu))+".pdf")
else:
  can2.SaveAs(args.output+"/bias_median_standard_param_mu_limit_2Sband.png")
  can2.SaveAs(args.output+"/bias_median_standard_param_mu_limit_2Sband.pdf")


