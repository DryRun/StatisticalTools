import os
import sys
import argparse
import math
from ROOT import *
from setTDRStyle import setTDRStyle
from array import array

usage = "usage: python plotFits.py -i <inputList> -o <outputdir>"
print usage

parser = argparse.ArgumentParser(description='Process options.')

parser.add_argument("-o", "--output", type=str, dest="output", default="./",
    help="the directory OUTDIR contains the output of the program",
    metavar="OUTDIR"
    )
    

args = parser.parse_args()
print args
############
ROOT.gStyle.SetOptFit(0000)
gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

minX_mass = 1118
maxX_mass = 6099 

massBins_list = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]

massBins = array("d",massBins_list)

inputFileData = TFile("/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/test_fit/dijetFitResults_FuncType0_nParFit4_MC_1fb-1_Dinko.root")
inputFitRes = TFile("mlfitgoldenDataset.root")
hist_mass = inputFileData.Get("hist_mass")
hist_binned = inputFileData.Get("hist_binned")
tree_fit_b = inputFitRes.Get("tree_fit_b")

#Set error to 1.8 in empty bins
for i in range(1,hist_mass.GetNbinsX()+1):
  if (hist_mass.GetBinContent(i) == 0):
    hist_mass.SetBinError(i,1.8/hist_mass.GetBinWidth(i))

c = TCanvas("c", "",339,117,695,841)
tree_fit_b.GetEntry(0)
#p1_val = tree_fit_b.p1
#p2_val = tree_fit_b.p2
#p3_val = tree_fit_b.p3
#integral_val = tree_fit_b.n_exp_final_binbin1_proc_background

fit_b = inputFitRes.Get("fit_b")
list = fit_b.floatParsFinal()
p1 = list.find("p1")
p2 = list.find("p2")
p3 = list.find("p3")
integral = list.find("shapeBkg_background_bin1__norm")

p1_val = p1.getVal()
p2_val = p2.getVal()
p3_val = p3.getVal()
integral_val = integral.getVal() * hist_binned.Integral(hist_binned.FindBin(minX_mass),hist_binned.FindBin(maxX_mass))
p1_error = p1.getError()
p2_error = p2.getError()
p3_error = p3.getError()
integral_error = integral.getError() * hist_binned.Integral(hist_binned.FindBin(minX_mass),hist_binned.FindBin(maxX_mass))


print str(p1_val)+"  "+str(p2_val)+"  "+str(p3_val)+"  "+str(integral_val)
print str(p1_error)+"  "+str(p2_error)+"  "+str(p3_error)+"  "+str(integral_error)

background_noNorm = TF1("background","( TMath::Power(1-x/13000,[0]) ) / ( TMath::Power(x/13000,[1]+[2]*log(x/13000)) )",minX_mass,maxX_mass)
background_noNorm.SetParameter(0,p1_val)
background_noNorm.SetParameter(1,p2_val)
background_noNorm.SetParameter(2,p3_val)
norm = background_noNorm.Integral(minX_mass,maxX_mass)
print "norm : "+str(norm)  
p0_val = integral_val/norm
p0_error = integral_error/norm
print ("p0 : %.2e  +/- %.2e" %(p0_val,p0_error))  

background = TF1("background","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
background.SetParameter(0,p0_val)
background.SetParameter(1,p1_val)
background.SetParameter(2,p2_val)
background.SetParameter(3,p3_val)
background.SetLineColor(kRed)

#hist_mass.Draw("hist p")
#background.Draw("same")

###########################################
# fit residuals and chi2
###########################################
hist_fit_residual_vsMass = TH1D("hist_fit_residual_vsMass","hist_fit_residual_vsMass",len(massBins)-1,massBins)
NumberOfObservations_VarBin = 0
chi2_VarBin = 0.

for bin in range(1,len(massBins)):
    
  if( hist_mass.GetXaxis().GetBinLowEdge(bin)>=minX_mass  and hist_mass.GetXaxis().GetBinUpEdge(bin)<=maxX_mass ):
    NumberOfObservations_VarBin+=1
    data = hist_mass.GetBinContent(bin)
    err_data = hist_mass.GetBinError(bin)
    if( data == 0 ):
      err_data = 1.8 / hist_mass.GetBinWidth(bin)
      print "err_data %f" % err_data
    fit = background.Integral(hist_mass.GetXaxis().GetBinLowEdge(bin), hist_mass.GetXaxis().GetBinUpEdge(bin) ) 
    fit = fit / ( hist_mass.GetBinWidth(bin) )
    err_tot = err_data	  
    fit_residual = (data - fit) / err_tot
    err_fit_residual = 1
	    	  
    chi2_VarBin += pow( (data - fit) , 2 ) / pow( err_data , 2 )
    print "data, err_data, fit: "+str( data)+", "+str(err_data) +", " +str(fit)
    print "bin, fit residual : " + str(bin) + ", " +str(fit_residual)	  
    hist_fit_residual_vsMass.SetBinContent(bin,fit_residual)
    hist_fit_residual_vsMass.SetBinError(bin,err_fit_residual)
   
ndf_VarBin = NumberOfObservations_VarBin - 4 
print "============================"
print "NumberOfObservations_VarBin: " + str(NumberOfObservations_VarBin)
print "ndf_VarBin: " + str(ndf_VarBin)
print "chi2_VarBin: " +str(chi2_VarBin)
print "============================"   

###############################
c.GetWindowHeight() 
c.GetWindowWidth() 
c.SetLogy() 
c.Divide(1,2,0,0,0) 
c.cd(1) 
p11_1 = c.GetPad(1) 
p11_1.SetPad(0.01,0.23,0.99,0.98) 
p11_1.SetLogy() 
p11_1.SetRightMargin(0.05) 
p11_1.SetTopMargin(0.05) 

#Pave text for fit results
pave_fit1 = TPaveText(0.56,0.55,0.9,0.85,"NDC") 
pave_fit1.AddText("#chi^{2} / ndf = %.2f / %d" % (chi2_VarBin,ndf_VarBin)) 
pave_fit1.AddText("p0 =  %.2e #pm %.2e" % (p0_val,p0_error))
pave_fit1.AddText("p1 =  %.2e #pm %.2e" % (p1_val,p1_error))
pave_fit1.AddText("p2 =  %.2e #pm %.2e" % (p2_val,p2_error))
pave_fit1.AddText("p3 =  %.2e #pm %.2e" % (p3_val,p3_error))
pave_fit1.SetFillColor(0) 
pave_fit1.SetFillStyle(0) 
pave_fit1.SetBorderSize(1) 
pave_fit1.SetTextFont(42) 
pave_fit1.SetTextSize(0.03) 
pave_fit1.SetTextAlign(12)  

#Pave text
pave_fit = TPaveText(0.18,0.15,0.40,0.27,"NDC") 
pave_fit.AddText(" #sqrt{s} = 13 TeV") 
pave_fit.AddText("|#eta| < 2.5, |#Delta#eta| < 1.3") 
pave_fit.AddText("M_{jj} > 1100 GeV") 
pave_fit.AddText("Wide Jets") 
pave_fit.SetFillColor(0) 
pave_fit.SetLineColor(0) 
pave_fit.SetFillStyle(0) 
pave_fit.SetBorderSize(0) 
pave_fit.SetTextFont(42) 
pave_fit.SetTextSize(0.03) 
pave_fit.SetTextAlign(12)  


pt1 =  TPaveText(0.1284756,0.9602144,0.3887139,0.9902251,"brNDC") 
pt1.SetBorderSize(0) 
pt1.SetFillColor(0) 
pt1.SetFillStyle(0) 
pt1.SetLineColor(0) 
pt1.SetTextAlign(12) 
pt1.SetTextSize(0.035) 
text = pt1.AddText("CMS") 

pt2 =  TPaveText(0.45,0.96,0.65,0.99,"brNDC") 
pt2.SetBorderSize(0) 
pt2.SetFillColor(0) 
pt2.SetFillStyle(0) 
pt2.SetLineColor(0) 
pt2.SetTextAlign(12) 
pt2.SetTextSize(0.035) 
text2 = pt2.AddText("#sqrt{s} = 13 TeV") 

pt3 = TPaveText(0.7687988,0.9602144,0.9297357,0.9902251,"brNDC") 
pt3.SetBorderSize(0) 
pt3.SetFillColor(0) 
pt3.SetFillStyle(0) 
pt3.SetLineColor(0) 
pt3.SetTextAlign(12) 
pt3.SetTextSize(0.035) 
text3 = pt3.AddText("L= 1 fb^{-1}") 
#text3 = pt3.AddText("L= 10 fb^{-1}") 

vFrame = p11_1.DrawFrame(minX_mass,0.000000001,maxX_mass,10.0) 

vFrame.SetTitle("") 
vFrame.SetXTitle("Dijet Mass (GeV)") 
vFrame.GetXaxis().SetTitleSize(0.06) 
vFrame.SetYTitle("events / bin width") 

vFrame.GetYaxis().SetTitleSize(0.12) 
vFrame.GetYaxis().SetLabelSize(0.07) 
vFrame.GetYaxis().SetTitleOffset(0.50) 
vFrame.GetXaxis().SetTitleOffset(0.90) 
vFrame.GetXaxis().SetTitleSize(0.18) 
vFrame.GetXaxis().SetLabelSize(0.1) 

hist_mass.GetXaxis().SetRangeUser(minX_mass,maxX_mass) 
hist_mass.SetTitle("") 
hist_mass.SetLineColor(1) 
hist_mass.SetFillColor(1) 
hist_mass.SetLineColor(1) 
hist_mass.SetMarkerColor(1) 
hist_mass.SetMarkerStyle(20) 
hist_mass.SetMinimum(0.001) 

hist_mass.Draw("HIST P0E0") 
background.SetLineWidth(2) 
background.SetLineStyle(2) 
background.SetLineColor(2) 
background.Draw("same") 

leg = TLegend(0.5564991,0.4,0.8903575,0.575812) 
leg.SetTextSize(0.03146853) 
leg.SetLineColor(1) 
leg.SetLineStyle(1) 
leg.SetLineWidth(1) 
leg.SetFillColor(0) 
leg.SetMargin(0.35) 
leg.AddEntry(hist_mass,"pseudo data" ,"PL") 
leg.AddEntry(background,"fit to pseudodata","L") 
leg.Draw("same") 

pt1.Draw("same")  
pt2.Draw("same") 
pt3.Draw("same") 
pave_fit.Draw("same") 
pave_fit1.Draw("same") 

#redraw axis
#p11_1.RedrawAxis() 
#p11_1.Update() 
#cout << "MIN: " << p11_1.GetUxmin() << endl 
#cout << "MAX: " << p11_1.GetUxmax() << endl 
#
#---- Next PAD

c.cd(2) 
p11_2 = c.GetPad(2) 
p11_2.SetPad(0.01,0.02,0.99,0.24) 
p11_2.SetBottomMargin(0.35) 
p11_2.SetRightMargin(0.05) 
p11_2.SetGridx() 
p11_2.SetGridy() 
#c3_2.SetTickx(50) 

vFrame2 = p11_2.DrawFrame(p11_1.GetUxmin(), -4.5, p11_1.GetUxmax(), 4.5) 

vFrame2.SetTitle("") 
vFrame2.SetXTitle("Dijet Mass (GeV)") 
vFrame2.GetXaxis().SetTitleSize(0.06) 
vFrame2.SetYTitle("(pseudoData-Fit)/#sigma") 
vFrame2.GetYaxis().SetTitleSize(0.12) 
vFrame2.GetYaxis().SetLabelSize(0.07) 
vFrame2.GetYaxis().SetTitleOffset(0.50) 
vFrame2.GetXaxis().SetTitleOffset(0.90) 
vFrame2.GetXaxis().SetTitleSize(0.18) 
vFrame2.GetXaxis().SetLabelSize(0.1) 

hist_fit_residual_vsMass.GetXaxis().SetRangeUser(minX_mass,maxX_mass) 
hist_fit_residual_vsMass.GetYaxis().SetRangeUser(-4.,4.) 
hist_fit_residual_vsMass.SetLineWidth(0) 
hist_fit_residual_vsMass.SetFillColor(2) 
hist_fit_residual_vsMass.SetLineColor(1) 
hist_fit_residual_vsMass.Draw("SAMEHIST") 

line =  TLine(minX_mass,0,maxX_mass,0) 
line.Draw("") 
#c.Close() 

c.SaveAs(args.output+"/fit_goldenDataset.C")
c.SaveAs(args.output+"/fit_goldenDataset.png")
c.SaveAs(args.output+"/fit_goldenDataset.pdf")
c.Clear()


