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

parser.add_argument("-i", "--inputList", type=str, dest="inputList", default="",
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



args = parser.parse_args()
print args
#####################################

gStyle.SetOptFit(1111)

pull_mu_list = []
err_pull_mu_list = []
err_mass = []
#giulia --- to change 
mass = [2,3,4,5,6]

ins = open(args.inputList,"r")

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

  #extract mu limit 2 sigma band 
  filename = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/dijet_limits/higgsCombineNormal.Asymptotic.mH'+mass_string+'.root'
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
    inputFile = TFile("output_toys_mu"+str(int(args.mu))+"/mlfit"+sample+"_MLfit.root","read") 
  else:
    inputFile = TFile("output_toys_mu_limit_2Sband/mlfit"+sample+"_MLfit.root","read") 
   
  
  tree_fit_sb = inputFile.Get("tree_fit_sb")
  h_pull_mu = TH1F("h_pull_mu", "S_{true} - S_{fit} / err", 100, -3., 3.)
  if not (args.mu == -999):
    mu = args.mu
  else:
    mu = y[4]

  for event in tree_fit_sb:
    if tree_fit_sb.mu > mu:
      err_mu = tree_fit_sb.muLoErr
      pull_mu = (tree_fit_sb.mu - mu)/tree_fit_sb.muLoErr
    else:
      err_mu = tree_fit_sb.muHiErr
      pull_mu = (tree_fit_sb.mu - mu)/tree_fit_sb.muHiErr
    #print "mu_fit - mu_exp  = "+str(tree_fit_sb.mu)
    #print "err_mu = "+str(err_mu)
    #print "pull_mu = "+str(pull_mu)
    h_pull_mu.Fill( pull_mu )
  
  h_pull_mu.Print() 

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

  #LOOP finished

### Draw bias plot ####
v_mass = array("d",mass)
v_err_mass = array("d",err_mass)
v_pull_mu = array("d",pull_mu_list)
v_pull_mu_err = array("d", err_pull_mu_list)

g_pull_mu = TGraphErrors(len(v_mass),v_mass,v_pull_mu,v_err_mass, v_pull_mu_err)
g_pull_mu_simple = TGraph(len(v_mass),v_mass,v_pull_mu )

canvasName = 'Bias'
can = TCanvas(canvasName,canvasName,900,600)
g_pull_mu.GetXaxis().SetTitle(' qg resonance Mass (TeV)')
g_pull_mu.GetYaxis().SetTitle('#mu pull')
#g_pull_mu.GetYaxis().CenterTitle(ROOT.kTRUE)
g_pull_mu.GetYaxis().SetNdivisions(510)
g_pull_mu.GetYaxis().SetRangeUser(-3,3)
g_pull_mu.SetFillStyle(3001)
g_pull_mu.SetFillColor(kAzure+1)
g_pull_mu_simple.SetMarkerStyle(20)
g_pull_mu_simple.SetMarkerSize(0.5)
g_pull_mu_simple.SetMarkerColor(kBlack)


#xmin = g_pull_mu.GetXaxis().GetXmin()
#xmax = g_pull_mu.GetXaxis().GetXmax()
xmin = mass[0]
xmax = mass[len(mass)-1]

line_0  = TLine(xmin,0.,xmax,0.);
line_up = TLine(xmin,0.15,xmax,0.15);
line_do = TLine(xmin,-0.15,xmax,-0.15);

line_up.SetLineStyle(2)
line_do.SetLineStyle(2)
line_0.SetLineColor(kBlack)
line_up.SetLineColor(kBlack)
line_do.SetLineColor(kBlack)
can.cd()
g_pull_mu.Draw("ae3")
g_pull_mu_simple.Draw("p")
line_0.Draw()
line_up.Draw()
line_do.Draw()

if not(args.mu==-999.):
  can.SaveAs(args.output+"/bias_standard_param_mu"+str(int(args.mu))+".png")
  can.SaveAs(args.output+"/bias_standard_param_mu"+str(int(args.mu))+".pdf")
else:
  can.SaveAs(args.output+"/bias_standard_param_mu_limit_2Sband.png")
  can.SaveAs(args.output+"/bias_standard_param_mu_limit_2Sband.pdf")


