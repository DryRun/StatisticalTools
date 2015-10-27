import os
import sys
import argparse
import math
from ROOT import *
from array import array
import CMS_lumi, tdrstyle

usage = "usage: python plotFits.py -i <inputList> -o <outputdir>"
print usage

parser = argparse.ArgumentParser(description='Process options.')

parser.add_argument("--inputFileData", type=str, dest="inputFileData", default="",
    help="the input file with the spectrum to fit"
    )
parser.add_argument("--inputFileRes", type=str, dest="inputFileRes", default="mlfitgoldenDataset.root",
    help="the input file with results of the fit"
    )
parser.add_argument("--inputWS", type=str, dest="inputWS", default="",
    help="the input workspace"
    )
parser.add_argument("--lumi", type=float, dest="lumi", default=1000,
    help="luminosity"
    )
parser.add_argument("--xsec", type=float, dest="xsec", default=0.6700E-01,
    help="luminosity"
    )
parser.add_argument("-o", "--output", type=str, dest="output", default="./",
    help="the directory OUTDIR contains the output of the program",
    metavar="OUTDIR"
    )
parser.add_argument("--tag", type=str, dest="suffix", default="test",
    help="tag to add at the end of the name",
    )
    
    

args = parser.parse_args()
print args
############
ROOT.gStyle.SetOptFit(0000)
gROOT.Reset()
#set the tdr style
tdrstyle.setTDRStyle()

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
CMS_lumi.lumi_13TeV = str(int(args.lumi))+" pb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = ""
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

os.system("mkdir -p "+args.output  )
minX_mass = 1181
maxX_mass = 7589 
#maxX_mass = 6564 


massBins_list = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]
massBins_list_limited = [ 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564]#, 6808, 7060, 7320, 7589]

massBins = array("d",massBins_list)
massBins_limited = array("d",massBins_list_limited)


#inputFileData = TFile("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/test_fit/dijetFitResults_FuncType0_nParFit4_MC_1fb-1_Dinko.root")
inputFileData = TFile(args.inputFileData)
#inputFitRes = TFile("mlfitgoldenDataset.root")
inputFitRes = TFile(args.inputFileRes)

hist_1GeV = inputFileData.Get("h_dat")
#hist_binned = inputFileData.Get("hist_binned")
hist_binned = hist_1GeV.Rebin(len(massBins)-1,"hist_binned",massBins)
hist_mass = TH1F("hist_mass","",len(massBins)-1,massBins) 
for i in range(1,len(massBins)+1):
  hist_mass.SetBinContent(i,hist_binned.GetBinContent(i)/(hist_binned.GetBinWidth(i)*args.lumi))
  hist_mass.SetBinError(i,hist_binned.GetBinError(i)/(hist_binned.GetBinWidth(i)*args.lumi))

tree_fit_b = inputFitRes.Get("tree_fit_b")
tree_fit_sb = inputFitRes.Get("tree_fit_sb")

#Set error to 1.8 in empty bins
for i in range(1,hist_mass.GetNbinsX()+1):
  if (hist_mass.GetBinContent(i) == 0):
    hist_mass.SetBinError(i,1.8/hist_mass.GetBinWidth(i))

input_workspace = TFile(args.inputWS)
print input_workspace
workspace = input_workspace.Get("w")
workspace.Print()
signal=workspace.pdf("signal")
datahist_sig = signal.dataHist()
set_sig=datahist_sig.get()
mjj_sig=set_sig.find("mjj")
h_sig_1GeV= signal.createHistogram("h_sig_1GeV",mjj_sig)
h_sig_1GeV.Print()
 
 
######################################################
#data in TGraph format (hist binned)
alpha = 1 - 0.6827;

x=[]
y=[]
exl=[]
exh=[]
eyl=[]
eyh=[]
x_mc=[]
y_mc=[]

for i in range(0,len(massBins)):
  n    = hist_binned.GetBinContent(i+1)
  dm   = hist_binned.GetBinWidth(i+1)
  mass = hist_binned.GetBinCenter(i+1)
  xl   = hist_binned.GetBinLowEdge(i+1)
  xh   = xl+dm
  x.append( (xl+xh)/2.)
  exl.append( dm/2.)
  exh.append( dm/2.)
  y.append( n / (dm*args.lumi))
  l = 0.5*TMath.ChisquareQuantile(alpha/2,2*n)
  h = 0.5*TMath.ChisquareQuantile(1-alpha/2,2*(n+1))
  eyl.append( (n-l)/(args.lumi*dm) )
  eyh.append( (h-n)/(args.lumi*dm) )
  #print "%f   %f    %f    %f    %f     %f" % (x[i],y[i],exl[i],exh[i],eyl[i],eyh[i])

vx = array("f",x)
vy = array("f",y)
vexl = array("f",exl)
vexh = array("f",exh)
veyl = array("f",eyl)
veyh = array("f",eyh)

#data in TGraph format
g = TGraphAsymmErrors(len(massBins)-1,vx,vy,vexl,vexh,veyl,veyh)
g.SetName("g_data")
g.Print()
###########################


tree_fit_sb.GetEntry(0)
p1_val = tree_fit_sb.p1
p2_val = tree_fit_sb.p2
p3_val = tree_fit_sb.p3
mu_val = tree_fit_sb.mu
integral_val = tree_fit_sb.n_exp_final_binbin1_proc_background  
tree_fit_b.GetEntry(0)
p1_val_b = tree_fit_b.p1
p2_val_b = tree_fit_b.p2
p3_val_b = tree_fit_b.p3
mu_val_b = tree_fit_b.mu
integral_val_b = tree_fit_b.n_exp_final_binbin1_proc_background  

background_noNorm = TF1("background","( TMath::Power(1-x/13000,[0]) ) / ( TMath::Power(x/13000,[1]+[2]*log(x/13000)) )",minX_mass,maxX_mass)
background_noNorm.SetParameter(0,p1_val_b)
background_noNorm.SetParameter(1,p2_val_b)
background_noNorm.SetParameter(2,p3_val_b)
norm_b = background_noNorm.Integral(minX_mass,maxX_mass)
p0_val_b = integral_val_b/norm_b
#print "integral_val: "+str(integral_val_b)+"   norm : "+str(norm_b)  

background_noNorm_SplusB = TF1("background_SplusB","( TMath::Power(1-x/13000,[0]) ) / ( TMath::Power(x/13000,[1]+[2]*log(x/13000)) )",minX_mass,maxX_mass)
background_noNorm_SplusB.SetParameter(0,p1_val)
background_noNorm_SplusB.SetParameter(1,p2_val)
background_noNorm_SplusB.SetParameter(2,p3_val)
norm = background_noNorm_SplusB.Integral(minX_mass,maxX_mass)
p0_val = integral_val/norm


#h_sig_1GeV.Scale(args.xsec*mu_val*args.lumi/h_sig_1GeV.Integral())
h_sig_1GeV.Scale(-1*mu_val/h_sig_1GeV.Integral())
int_sig = h_sig_1GeV.Integral()
#print "int_sig = %f" % int_sig
h_sig_rebin=h_sig_1GeV.Rebin(len(massBins)-1,"h_sig_rebin",massBins)
h_sig=TH1F("h_sig","",len(massBins)-1,massBins)
for i in range(1,len(massBins)):
  h_sig.SetBinContent(i, h_sig_rebin.GetBinContent(i)/(h_sig_rebin.GetBinWidth(i)*args.lumi))


#tree_fit_b.GetEntry(0)
#p1_val = tree_fit_b.p1
#p2_val = tree_fit_b.p2
#p3_val = tree_fit_b.p3
#integral_val = tree_fit_b.n_exp_final_binbin1_proc_background

#fit_b = inputFitRes.Get("fit_b")
#list = fit_b.floatParsFinal()
#p1 = list.find("p1")
#p2 = list.find("p2")
#p3 = list.find("p3")
#integral = list.find("shapeBkg_background_bin1__norm")

#p1_val = p1.getVal()
#p2_val = p2.getVal()
#p3_val = p3.getVal()
#integral_val = integral.getVal() #* hist_binned.Integral(hist_binned.FindBin(minX_mass),hist_binned.FindBin(maxX_mass))
#p1_error = p1.getError()
#p2_error = p2.getError()
#p3_error = p3.getError()
#integral_error = integral.getError() #* hist_binned.Integral(hist_binned.FindBin(minX_mass),hist_binned.FindBin(maxX_mass))


#print str(p1_val)+"  "+str(p2_val)+"  "+str(p3_val)+"  "+str(integral_val)
#print str(p1_error)+"  "+str(p2_error)+"  "+str(p3_error)+"  "+str(integral_error)
#print str(p1_val_b)+"  "+str(p2_val_b)+"  "+str(p3_val_b)+"  "+str(integral_val_b)

#background_noNorm = TF1("background","( TMath::Power(1-x/13000,[0]) ) / ( TMath::Power(x/13000,[1]+[2]*log(x/13000)) )",minX_mass,maxX_mass)
#background_noNorm.SetParameter(0,p1_val)
#background_noNorm.SetParameter(1,p2_val)
#background_noNorm.SetParameter(2,p3_val)
#norm = background_noNorm.Integral(minX_mass,maxX_mass)
#data_integral = hist_binned.Integral(1,len(massBins))
#print "data : "+str(data_integral)
#print "norm : "+str(norm)  
#p0_val = data_integral * (integral_val / norm)
#p0_error = data_integral * (integral_error / norm)
#p0_val = integral_val
#p0_error = integral_error
#print ("p0 : %.2e  +/- %.2e" %(p0_val,p0_error))  

background = TF1("background","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
background.SetParameter(0,p0_val_b/args.lumi)
background.SetParameter(1,p1_val_b)
background.SetParameter(2,p2_val_b)
background.SetParameter(3,p3_val_b)
background.SetLineColor(kRed)
##background component of S+B
background_SplusB = TF1("background_SplusB","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
background_SplusB.SetParameter(0,p0_val/args.lumi)
background_SplusB.SetParameter(1,p1_val)
background_SplusB.SetParameter(2,p2_val)
background_SplusB.SetParameter(3,p3_val)

#hist_mass.Draw("hist p")
#background.Draw("same")
## signal plus background ##
h_s_plus_b = TH1F("h_s_plus_b","",len(massBins_limited)-1,massBins_limited) 
for ibin in range(1,len(massBins_limited)):
  h_s_plus_b.SetBinContent(ibin, 0)
  bin_center=h_s_plus_b.GetBinCenter(ibin)
  if(h_s_plus_b.GetBinCenter(ibin) > minX_mass and h_s_plus_b.GetBinCenter(ibin) < maxX_mass):
    print "b = "+str(background_SplusB.Eval(bin_center))+"   s = "+str(h_sig.GetBinContent( h_sig.FindBin(bin_center)))+"   s+b = "+str(background_SplusB.Eval(bin_center) + h_sig.GetBinContent( h_sig.FindBin(bin_center)))
    #h_s_plus_b.SetBinContent(ibin, background_SplusB.Eval(h_s_plus_b.GetBinCenter(ibin)) + h_sig.GetBinContent( h_sig.FindBin(bin_center)))
    b_comp = background_SplusB.Integral(h_s_plus_b.GetXaxis().GetBinLowEdge(ibin),h_s_plus_b.GetXaxis().GetBinUpEdge(ibin)) / h_s_plus_b.GetBinWidth(ibin)
    h_s_plus_b.SetBinContent(ibin, b_comp + h_sig.GetBinContent( h_sig.FindBin(bin_center)))

h_s_plus_b.Print()


###########################################
# fit residuals and chi2 variable binning
###########################################
hist_fit_residual_vsMass = TH1D("hist_fit_residual_vsMass","hist_fit_residual_vsMass",len(massBins)-1,massBins)
hist_fit_residual_vsMass_bkg = TH1D("hist_fit_residual_vsMass_bkg","hist_fit_residual_vsMass_bkg",len(massBins)-1,massBins)
NumberOfObservations_VarBin = 0
chi2_VarBin_B = 0.
chi2_VarBin_SB = 0.

for bin in range(1,len(massBins)):
  
  if( hist_mass.GetXaxis().GetBinLowEdge(bin)>=minX_mass and hist_mass.GetXaxis().GetBinUpEdge(bin)<=maxX_mass ):
    NumberOfObservations_VarBin += 1
    #print "bin content = " + str(hist_mass.GetBinContent(bin)) + "   graph y = " + str(vy[bin-1]) + "  error y low = " + str(g.GetErrorYlow(bin-1))
    data = hist_mass.GetBinContent(bin)
    err_data_low = g.GetErrorYlow(bin-1) 
    err_data_high = g.GetErrorYhigh(bin-1)
    
    fit_SB =  h_s_plus_b.GetBinContent(h_s_plus_b.FindBin(hist_mass.GetXaxis().GetBinLowEdge(bin)))
    if(fit_SB > data): err_tot_SB = err_data_high
    else: err_tot_SB = err_data_low
    if err_tot_SB==0:
      print "*** data = %f  fit = %f  err = 0 ***" %  (data, fit_SB)
      fit_residual_SB = 0
    else:  
      fit_residual_SB = (data - fit_SB) / err_tot_SB
    chi2_VarBin_SB += pow( fit_residual_SB,2 )
    
    fit_B = background.Integral(hist_mass.GetXaxis().GetBinLowEdge(bin) , hist_mass.GetXaxis().GetBinUpEdge(bin) )
    fit_B = fit_B / ( hist_mass.GetBinWidth(bin) )
    if(fit_B > data): err_tot_B = err_data_high
    else: err_tot_B = err_data_low
    fit_residual_B = (data - fit_B) / err_tot_B
    chi2_VarBin_B += pow( fit_residual_B , 2 )

    hist_fit_residual_vsMass.SetBinContent(bin,fit_residual_SB)
    hist_fit_residual_vsMass_bkg.SetBinContent(bin,fit_residual_B)


ndf_VarBin = NumberOfObservations_VarBin #-1 
print "============ CHI2 ==============" 
print "NumberOfObservations_VarBin: %d" %  NumberOfObservations_VarBin
print "ndf_VarBin_SB: %d" % (NumberOfObservations_VarBin-5) 
print "chi2_VarBin_SB: %f" % chi2_VarBin_SB
print "ndf_VarBin_B: %d" % (NumberOfObservations_VarBin-4) 
print "chi2_VarBin_B: %f" % chi2_VarBin_B
print "============================"   
  
    
#  if( hist_mass.GetXaxis().GetBinLowEdge(bin)>=minX_mass  and hist_mass.GetXaxis().GetBinUpEdge(bin)<=maxX_mass ):
#    NumberOfObservations_VarBin+=1
#    data = hist_mass.GetBinContent(bin)
#    err_data = hist_mass.GetBinError(bin)
#    if( data == 0 ):
#      err_data = 1.8 / hist_mass.GetBinWidth(bin)
#      print "err_data %f" % err_data
#    #background only
#    fit_B = background.Integral(hist_mass.GetXaxis().GetBinLowEdge(bin), hist_mass.GetXaxis().GetBinUpEdge(bin) ) 
#    fit_B = fit_B / ( hist_mass.GetBinWidth(bin) )
#    #S+B
#    fit = h_s_plus_b.GetBinContent(h_s_plus_b.FindBin(hist_mass.GetBinCenter(bin) )) 
#    err_tot = err_data	  
#    fit_residual = (data - fit) / err_tot
#    fit_residual_B = (data - fit_B) / err_tot
#	    	  
#    chi2_VarBin += pow( (data - fit) , 2 ) / pow( err_data , 2 )
#    print "data, err_data, fit: "+str( data)+", "+str(err_data) +", " +str(fit)
#    print "bin, fit residual : " + str(bin) + ", " +str(fit_residual)	  
#
#   
#
#ndf_VarBin = NumberOfObservations_VarBin - 4 
#
#print "============================"
#print "NumberOfObservations_VarBin: " + str(NumberOfObservations_VarBin)
#print "ndf_VarBin: " + str(ndf_VarBin)
#print "chi2_VarBin: " +str(chi2_VarBin)
#print "============================"   

###############################
#### Draw ##########
W = 600
H = 700
H_ref = 700 
W_ref = 600 
T = 0.08*H_ref
B = 0.12*H_ref
L = 0.12*W_ref
R = 0.04*W_ref

c = TCanvas("c","DijetMass cross section with Fit and QCD MC",W,H)
#c.GetWindowHeight() 
#c.GetWindowWidth() 
c.SetLogy() 
c.Divide(1,2,0,0,0) 
c.cd(1) 


#------------ pad 1  ----------------
c.cd(1)
p11_1 = c.GetPad(1)
p11_1.SetPad(0.01,0.26,0.99,0.98)
p11_1.SetLogy()
p11_1.SetRightMargin(0.05)
p11_1.SetTopMargin(0.05)
p11_1.SetFillColor(0)
p11_1.SetBorderMode(0)
p11_1.SetFrameFillStyle(0)
p11_1.SetFrameBorderMode(0)

##Pave text for fit results
#pave_fit1 = TPaveText(0.56,0.55,0.9,0.85,"NDC") 
#pave_fit1.AddText("#chi^{2} / ndf = %.2f / %d" % (chi2_VarBin,ndf_VarBin)) 
#pave_fit1.AddText("p0 =  %.2e #pm %.2e" % (p0_val,p0_error))
#pave_fit1.AddText("p1 =  %.2e #pm %.2e" % (p1_val,p1_error))
#pave_fit1.AddText("p2 =  %.2e #pm %.2e" % (p2_val,p2_error))
#pave_fit1.AddText("p3 =  %.2e #pm %.2e" % (p3_val,p3_error))
#pave_fit1.SetFillColor(0) 
#pave_fit1.SetFillStyle(0) 
#pave_fit1.SetBorderSize(1) 
#pave_fit1.SetTextFont(42) 
#pave_fit1.SetTextSize(0.03) 
#pave_fit1.SetTextAlign(12)  

#Pave text
pave_fit = TPaveText(0.2358691,0.04035043,0.5050171,0.1870085,"NDC")
  
#pave_fit.AddText("AK4 Jets")
pave_fit.AddText("Wide Jets")
pave_fit.AddText("|#eta| < 2.5, |#Delta#eta_{jj}| < 1.3")
#pave_fit.AddText("M_{jj} > 1.2 TeV")
pave_fit.SetFillColor(0)
pave_fit.SetLineColor(0)
pave_fit.SetFillStyle(0)
pave_fit.SetBorderSize(0)
pave_fit.SetTextFont(42)
pave_fit.SetTextSize(0.040)
pave_fit.SetTextAlign(12) 

vFrame = p11_1.DrawFrame(minX_mass,0.0000005,maxX_mass,5.0)

vFrame.SetTitle("")
vFrame.SetXTitle("Dijet Mass [GeV]")
vFrame.SetYTitle("d#sigma / dm_{jj}   [pb / GeV]")
vFrame.GetXaxis().SetTitleSize(0.06)
vFrame.GetXaxis().SetTitleOffset(0.95)
vFrame.GetXaxis().SetLabelSize(0.05)
vFrame.GetYaxis().SetTitleSize(0.06)
#vFrame.GetYaxis().SetTitleOffset(1.0)
vFrame.GetYaxis().SetLabelSize(0.05)

#hist_mass.GetXaxis().SetRangeUser(minX_mass,maxX_mass) 
#hist_mass.SetTitle("") 
#hist_mass.SetLineColor(1) 
#hist_mass.SetFillColor(1) 
#hist_mass.SetLineColor(1) 
#hist_mass.SetMarkerColor(1) 
#hist_mass.SetMarkerStyle(20) 
#hist_mass.SetMinimum(0.001) 
#
#hist_mass.Draw("HIST P0E0") 
g.GetXaxis().SetNdivisions(405)
g.SetMarkerSize(0.9)
g.SetMarkerStyle(20)
#g.Draw("pe0")
background.SetLineWidth(2) 
background.SetLineStyle(2) 
background.SetLineColor(2) 
background.Draw("same") 
h_s_plus_b.SetLineWidth(2) 
h_s_plus_b.SetLineStyle(1) 
h_s_plus_b.SetLineColor(kBlue) 
h_s_plus_b.Draw("c same")
g.Draw("pe0 same")

leg =  TLegend(0.4564991,0.62,0.9303575,0.90)
leg.SetTextSize(0.038)
leg.SetLineColor(0)
leg.SetLineStyle(1)
leg.SetLineWidth(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetMargin(0.35)
leg.AddEntry(hist_mass,"Data" ,"PL") 
leg.AddEntry(background,"B Fit","L") 
leg.AddEntry(h_s_plus_b,"S+B Fit","L") 

pave_fit.Draw("same") 
#pave_fit1.Draw("same") 
leg.Draw("samw")

#redraw axis
p11_1.RedrawAxis() 
p11_1.Update() 
#cout << "MIN: " << p11_1.GetUxmin() << endl 
#cout << "MAX: " << p11_1.GetUxmax() << endl 
#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(p11_1, iPeriod, iPos)

#
#---- Next PAD

c.cd(2) 
p11_2 = c.GetPad(2) 
p11_2.SetPad(0.01,0.02,0.99,0.27)
p11_2.SetBottomMargin(0.35)
p11_2.SetRightMargin(0.05)
p11_2.SetGridx()
p11_2.SetGridy()

vFrame2 = p11_2.DrawFrame(p11_1.GetUxmin(), -3.5, p11_1.GetUxmax(), 3.5) 

vFrame2.SetTitle("")
vFrame2.SetXTitle("Dijet Mass [GeV]")
vFrame2.GetXaxis().SetTitleSize(0.06)
vFrame2.SetYTitle("#frac{Data-Fit}{#sigma_{Data}}")
vFrame2.GetYaxis().SetTitleSize(0.15)
vFrame2.GetYaxis().SetTitleOffset(0.40)
vFrame2.GetYaxis().SetLabelSize(0.09)
vFrame2.GetXaxis().SetTitleSize(0.18)
vFrame2.GetXaxis().SetTitleOffset(0.90)
vFrame2.GetXaxis().SetLabelSize(0.15)

hist_fit_residual_vsMass.GetXaxis().SetNdivisions(405)
hist_fit_residual_vsMass.GetXaxis().SetRangeUser(minX_mass,maxX_mass) 
hist_fit_residual_vsMass.GetYaxis().SetRangeUser(-3.5,3.5) 
hist_fit_residual_vsMass.SetLineWidth(0) 
hist_fit_residual_vsMass.SetFillColor(kBlue) 
hist_fit_residual_vsMass.SetLineColor(1) 
hist_fit_residual_vsMass.Draw("samehist") 
hist_fit_residual_vsMass_bkg.SetLineWidth(2) 
hist_fit_residual_vsMass_bkg.SetLineColor(2) 
hist_fit_residual_vsMass_bkg.SetLineStyle(2) 
hist_fit_residual_vsMass_bkg.Draw("samehist")

line =  TLine(minX_mass,0,maxX_mass,0) 
line.Draw("") 
p11_2.RedrawAxis()
line2=TLine()
line2.DrawLine(p11_2.GetUxmin(), p11_2.GetUymax(), p11_2.GetUxmax(), p11_2.GetUymax())
line2.DrawLine(p11_2.GetUxmax(), p11_2.GetUymin(), p11_2.GetUxmax(), p11_2.GetUymax())
#c.Close() 

#c.SaveAs(args.output+"/fit_Data.C")
c.SaveAs(args.output+"/fit_Data_"+args.suffix+".png")
c.SaveAs(args.output+"/fit_Data_"+args.suffix+".pdf")
c.Clear()


