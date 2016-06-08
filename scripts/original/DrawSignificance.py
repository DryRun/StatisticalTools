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
parser.add_option("--toy",action="store", type="int", default=1,dest="toy")
(options, args) = parser.parse_args()
useSub = options.useSub
toy = options.toy

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

#mass = [2000,3000,4000,5000,6000]
#xsection = [21.49,1.625,0.176,0.0182,0.00186]
#acceptance = [0.5887, 0.6160, 0.6161, 0.6126, 0.6117]

sigXStimesA = {
    1000.0:  0.4101E+03,
    1100.0:  0.2620E+03,
    1200.0:  0.1721E+03,
    1300.0:  0.1157E+03,
    1400.0:  0.7934E+02,
    1500.0:  0.5540E+02,
    1600.0:  0.3928E+02,
    1700.0:  0.2823E+02,
    1800.0:  0.2054E+02,
    1900.0:  0.1510E+02,
    2000.0:  0.1121E+02,
    2100.0:  0.8390E+01,
    2200.0:  0.6328E+01,
    2300.0:  0.4807E+01,
    2400.0:  0.3674E+01,
    2500.0:  0.2824E+01,
    2600.0:  0.2182E+01,
    2700.0:  0.1694E+01,
    2800.0:  0.1320E+01,
    2900.0:  0.1033E+01,
    3000.0:  0.8116E+00,
    3100.0:  0.6395E+00,
    3200.0:  0.5054E+00,
    3300.0:  0.4006E+00,
    3400.0:  0.3182E+00,
    3500.0:  0.2534E+00,
    3600.0:  0.2022E+00,
    3700.0:  0.1616E+00,
    3800.0:  0.1294E+00,
    3900.0:  0.1038E+00,
    4000.0:  0.8333E-01,
    4100.0:  0.6700E-01,
    4200.0:  0.5392E-01,
    4300.0:  0.4344E-01,
    4400.0:  0.3503E-01,
    4500.0:  0.2827E-01,
    4600.0:  0.2283E-01,
    4700.0:  0.1844E-01,
    4800.0:  0.1490E-01,
    4900.0:  0.1205E-01,
    5000.0:  0.9743E-02,
    5100.0:  0.7880E-02,
    5200.0:  0.6373E-02,
    5300.0:  0.5155E-02,
    5400.0:  0.4169E-02,
    5500.0:  0.3371E-02,
    5600.0:  0.2725E-02,
    5700.0:  0.2202E-02,
    5800.0:  0.1779E-02,
    5900.0:  0.1437E-02,
    6000.0:  0.1159E-02,
    6100.0:  0.9353E-03,
    6200.0:  0.7541E-03,
    6300.0:  0.6076E-03,
    6400.0:  0.4891E-03,
    6500.0:  0.3935E-03,
    6600.0:  0.3164E-03,
    6700.0:  0.2541E-03,
    6800.0:  0.2039E-03,
    6900.0:  0.1635E-03,
    7000.0:  0.1310E-03,
    7100.0:  0.1049E-03,
    7200.0:  0.8385E-04,
    7300.0:  0.6699E-04,
    7400.0:  0.5347E-04,
    7500.0:  0.4264E-04,
    7600.0:  0.3397E-04,
    7700.0:  0.2704E-04,
    7800.0:  0.2151E-04,
    7900.0:  0.1709E-04,
    8000.0:  0.1357E-04,
    8100.0:  0.1077E-04,
    8200.0:  0.8544E-05,
    8300.0:  0.6773E-05,
    8400.0:  0.5367E-05,
    8500.0:  0.4251E-05,
    8600.0:  0.3367E-05,
    8700.0:  0.2666E-05,
    8800.0:  0.2112E-05,
    8900.0:  0.1673E-05,
    9000.0:  0.1326E-05
                       
}


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

k = 0

#for im in mass:
for im in range(1200,6100,100):
  filename = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/significance_pseudodataDinko/higgsCombineQstar'+str(int(im))+'_toy'+str(toy)+'.ProfileLikelihood.mH120.root' 
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
  #y1.append(y[2]*xsection[k]*acceptance[k])
  #y2.append(y[2]*xsection[k]*acceptance[k])
  #yExp.append(y[2]*xsection[k]*acceptance[k])
  yObs.append(y[0])


  #ey1l.append((y[2]-y[1])*xsection[k]*acceptance[k])
  #ey2l.append((y[2]-y[0])*xsection[k]*acceptance[k])
  #ey1h.append((y[3]-y[2])*xsection[k]*acceptance[k])
  #ey2h.append((y[4]-y[2])*xsection[k]*acceptance[k])
 
  print "observed significance :  m="+str(im)+"  "+str(y[0])  
  
  k+=1
 

#xsecTimesAcc = []
#for i in range(0,len(xsection)):
#  xsecTimesAcc.append(xsection[i]*acceptance[i])

#vxs = array('d',xsecTimesAcc)
vx = array('d',x)
vyObs = array('d',yObs)
#vyExp = array('d',yExp)
#vy1 = array('d',y1)
#vy2 = array('d',y2)
#vexl = array('d',exl)
#vexh = array('d',exh)
#vey1l = array('d',ey1l)
#vey1h = array('d',ey1h)
#vey2l = array('d',ey2l)
#vey2h = array('d',ey2h)

#g1   = TGraphAsymmErrors(len(vx),vx,vy1,vexl,vexh,vey1l,vey1h)
#g2   = TGraphAsymmErrors(len(vx),vx,vy2,vexl,vexh,vey2l,vey2h)
gObs = TGraph(len(vx),vx,vyObs)
#gExp = TGraph(len(vx),vx,vyExp)
#gxs  = TGraph(len(vx),vx,vxs)

#g1.SetFillColor(ROOT.kGreen)
#g1.SetLineColor(ROOT.kGreen)
#g2.SetFillColor(ROOT.kYellow)
#g2.SetLineColor(ROOT.kYellow)
#gExp.SetLineWidth(2)
#gExp.SetLineStyle(9)
gObs.SetLineWidth(2)
gObs.SetLineStyle(1)
gObs.SetLineColor(ROOT.kBlue+1)
gObs.SetMarkerColor(ROOT.kBlue+1)
gObs.SetMarkerStyle(21)
gObs.SetMarkerSize(1.5)
#gxs.SetLineStyle(5)
#gxs.SetLineWidth(3)
#gxs.SetLineColor(ROOT.kRed)
    
canvasName = 'Significance'
can = TCanvas(canvasName,canvasName,900,600)
#gPad.SetLogy()
gObs.GetXaxis().SetTitle(' qg resonance Mass (TeV)')
gObs.GetYaxis().SetTitle('Observed significance')
gObs.GetYaxis().CenterTitle(ROOT.kTRUE)
gObs.GetYaxis().SetNdivisions(510)
gObs.GetYaxis().SetRangeUser(0,2)
#g2.Draw('AE3')
#g1.Draw('sameE3')
#gExp.Draw('sameL')
gObs.Draw('ALP')
#gxs.Draw('sameL')

leg = TLegend(0.65,0.65,0.9,0.9)
leg.SetHeader('Simple Dijet Selection')
leg.AddEntry(gObs,'Observed','LP')
#leg.AddEntry(gExp,'Expected','L')
#leg.AddEntry(g1,'Expected #pm 1 #sigma','F')
#leg.AddEntry(g2,'Expected #pm 2 #sigma','F')
#leg.AddEntry(gxs,'q*#rightarrow jj','L')
leg.SetFillColor(0)
leg.SetTextFont(42)
leg.SetTextSize(0.04)
leg.Draw()

can.SaveAs("significance_dijet13TeV_test_toy"+str(toy)+".png")
can.SaveAs("significance_dijet13TeV_test_toy"+str(toy)+".pdf")
#can.SaveAs("significance_dijet13TeV_test_fitDataBonly.png")
#can.SaveAs("significance_dijet13TeV_test_fitDataBonly.pdf")


#----- keep the GUI alive ------------
if __name__ == '__main__':
  rep = ''
  while not rep in ['q','Q']:
    rep = raw_input('enter "q" to quit: ')
    if 1 < len(rep):
      rep = rep[0]
