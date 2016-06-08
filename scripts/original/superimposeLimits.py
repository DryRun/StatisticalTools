#!usr/bin/python
import sys, getopt
import ROOT
from ROOT import TFile, TTree, TH1F, TCanvas, TLegend, TGraph, TGraphAsymmErrors, gROOT, gPad
from array import array
from setTDRStyle import setTDRStyle
import optparse

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--useSub",action="store_true",default=False,dest="useSub")
(options, args) = parser.parse_args()useSub = options.useSub

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

mass = [2000,3000,4000,5000,6000]
xsection = [21.49,1.625,0.176,0.0182,0.00186]
acceptance = [0.5887, 0.6160, 0.6161, 0.6126, 0.6117]
x = []
exl = []
exh = []
y1 = []
y2 = []
yExp = []
yObs = []
ey1l = []
ey1h = []
ey2l = []
ey2h = []  

##### Dinko's plot #####
masses_list =[1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6100.0, 6200.0, 6300.0, 6400.0, 6500.0, 6600.0, 6700.0, 6800.0, 6900.0, 7000.0]  
masses = array('d', [1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6100.0, 6200.0, 6300.0, 6400.0, 6500.0, 6600.0, 6700.0, 6800.0, 6900.0, 7000.0])

massesTeV = []
for i in range(0,len(masses_list)):
  massesTeV.append(0.001*masses_list[i])

massesTeV_v = array('d',massesTeV)

xs_obs_limits = array('d', [20.9545, 2.45626, 0.416225, 1.01931, 2.82771, 2.19572, 1.16475, 0.465517, 0.312011, 0.257445, 0.232343, 0.186241, 0.14177, 0.161448, 0.223433, 0.235894, 0.22298, 0.192951, 0.164016, 0.141445, 0.132672, 0.144567, 0.154099, 0.139979, 0.109849, 0.0806352, 0.0624284, 0.0528164, 0.0509088, 0.0567577, 0.0680257, 0.0716774, 0.0677652, 0.0582306, 0.0484861, 0.0427263, 0.0406525, 0.0400296, 0.0383613, 0.0356998, 0.0308406, 0.0253015, 0.0198666, 0.0152415, 0.0123197, 0.0110247, 0.00929554, 0.00821889, 0.00775016, 0.00722283, 0.00679991, 0.00656296, 0.00640062, 0.00628973, 0.0062155, 0.00616987, 0.00614263, 0.00613129, 0.00613597])

xs_exp_limits = array('d', [2.63206, 2.14471, 1.63841, 1.48217, 1.05316, 0.81135, 0.67681, 0.62097, 0.513012, 0.431082, 0.38607, 0.341819, 0.298624, 0.270375, 0.245451, 0.199897, 0.178578, 0.172584, 0.150738, 0.130131, 0.118085, 0.100697, 0.0934993, 0.0868947, 0.0863419, 0.0747488, 0.0630101, 0.0636891, 0.051635, 0.0489897, 0.0445535, 0.0413454, 0.0356299, 0.0323483, 0.0317702, 0.028163, 0.0258878, 0.0248127, 0.0219089, 0.0212312, 0.0193777, 0.0184784, 0.0187839, 0.0172933, 0.0165523, 0.0149275, 0.0144313, 0.0136047, 0.013115, 0.0127536, 0.0121895, 0.0116958, 0.0112944, 0.010989, 0.010528, 0.0103265, 0.00995773, 0.00966355, 0.00945783])
############

k = 0

for im in mass:
  filename = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/dijet_limits/higgsCombineQstar'+str(im)+'_limit.Asymptotic.mH120.root' 
  if useSub:
    filename = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/dijet_limits/higgsCombineSub.Asymptotic.mH'+str(im)+'.root' 
  print filename
  inf = TFile.Open(filename)
  tr = inf.Get('limit')
  y = []
  N = tr.GetEntriesFast()
  for i in xrange(N): 
    tr.GetEntry(i)
    y.append(tr.limit)
  
  x.append(im*0.001)
  exl.append(0.0)
  exh.append(0.0)
  y1.append(y[2]*xsection[k]*acceptance[k])
  y2.append(y[2]*xsection[k]*acceptance[k])
  yExp.append(y[2]*xsection[k]*acceptance[k])
  yObs.append(y[5]*xsection[k]*acceptance[k])

  ey1l.append((y[2]-y[1])*xsection[k]*acceptance[k])
  ey2l.append((y[2]-y[0])*xsection[k]*acceptance[k])
  ey1h.append((y[3]-y[2])*xsection[k]*acceptance[k])
  ey2h.append((y[4]-y[2])*xsection[k]*acceptance[k])
  
  k+=1
 

xsecTimesAcc = []
for i in range(0,len(xsection)):
  xsecTimesAcc.append(xsection[i]*acceptance[i])

#vxs = array('d',xsection)
vxs = array('d',xsecTimesAcc)
vx = array('d',x)
vyObs = array('d',yObs)
vyExp = array('d',yExp)
vy1 = array('d',y1)
vy2 = array('d',y2)
vexl = array('d',exl)
vexh = array('d',exh)
vey1l = array('d',ey1l)
vey1h = array('d',ey1h)
vey2l = array('d',ey2l)
vey2h = array('d',ey2h)

g1   = TGraphAsymmErrors(len(vx),vx,vy1,vexl,vexh,vey1l,vey1h)
g2   = TGraphAsymmErrors(len(vx),vx,vy2,vexl,vexh,vey2l,vey2h)
gObs = TGraph(len(vx),vx,vyObs)
gExp = TGraph(len(vx),vx,vyExp)
gxs  = TGraph(len(vx),vx,vxs)

gDinko = TGraph(len(masses), massesTeV_v, xs_exp_limits)

g1.SetFillColor(ROOT.kYellow)
g1.SetLineColor(ROOT.kYellow)
g2.SetFillColor(ROOT.kGreen)
g2.SetLineColor(ROOT.kGreen)
gExp.SetLineWidth(2)
gExp.SetLineStyle(9)
gObs.SetLineWidth(2)
gObs.SetLineStyle(1)
gObs.SetLineColor(ROOT.kBlue+1)
gObs.SetMarkerColor(ROOT.kBlue+1)
gObs.SetMarkerStyle(21)
gObs.SetMarkerSize(1.5)
gxs.SetLineStyle(5)
gxs.SetLineWidth(3)
gxs.SetLineColor(ROOT.kRed)
gDinko.SetLineColor(ROOT.kBlack)
gDinko.SetLineStyle(2)
gDinko.SetLineWidth(2)


canvasName = 'Limits'
if useSub:
  canvasName = 'Limits_Sub'
can = TCanvas(canvasName,canvasName,900,600)
gPad.SetLogy()
g2.GetXaxis().SetTitle(' qg resonance Mass (TeV)')
g2.GetYaxis().SetTitle('#sigma #times A #times BR (q*#rightarrow jj) (pb)')
#g2.GetYaxis().SetTitle('#sigma #times BR (q*#rightarrow jj) (pb)')
g2.GetYaxis().CenterTitle(ROOT.kTRUE)
g2.GetYaxis().SetNdivisions(510)
g2.GetYaxis().SetRangeUser(1e-3,2)
g2.Draw('AE3')
g1.Draw('sameE3')
gExp.Draw('sameL')
gObs.Draw('sameLP')
gxs.Draw('sameL')
gDinko.Draw('sameL')

leg = TLegend(0.6,0.65,0.9,0.9)
if useSub:
  leg.SetHeader('Substructure Selection')
else:
  leg.SetHeader('Simple Dijet Selection')
leg.AddEntry(gObs,'Observed','LP')
leg.AddEntry(gExp,'Expected (CLs Asymptotic combine)','L')
leg.AddEntry(gDinko,'Expected (Bayesian dijet code)','L')
leg.AddEntry(g1,'Expected #pm 1 #sigma','F')
leg.AddEntry(g2,'Expected #pm 2 #sigma','F')
leg.AddEntry(gxs,'q*#rightarrow jj','L')
leg.SetFillColor(0)
leg.SetTextFont(42)
leg.SetTextSize(0.02525)
leg.Draw()

can.SaveAs("limits_dijet13TeV_test.png")
can.SaveAs("limits_dijet13TeV_test.pdf")

#----- keep the GUI alive ------------
if __name__ == '__main__':
  rep = ''
  while not rep in ['q','Q']:
    rep = raw_input('enter "q" to quit: ')
    if 1 < len(rep):
      rep = rep[0]
