#!/usr/bin/env python

import sys, os, copy, re, imp, math, multiprocessing
sys.argv.append( '-b' ) #for batch mode
from ROOT import * 
from array import array
import CMS_lumi, setTDRStyle


#==================
# Input files
#==================

#inputFileWorkspace = TFile("../datacards/workspace_gg_m1400.root")    
#inputFileFit = TFile("../datacards/mlfit_gg_m1400.root")
inputFileWorkspace = TFile("../datacards/workspace_qg_m4000.root")    
inputFileFit = TFile("../datacards/mlfit_qg_m4000.root")
#inputFileFit = TFile("../datacards/mlfit.root")

#==================
# Options
#==================

showCrossSection = 1 #1=cross section [pb] , 0=number of events/GeV 
drawSignalShapeAlsoAlone = 1
lumiValue = 1769 #[pb]        #---> take it from the workspace FIXME
fixedRange = 1 #1=YES , 0=NO  (the option works only if showCrossSection=1; otherwise=0)
minY = 0.0000001
maxY = 20
if showCrossSection==1:
    lumi = lumiValue 
else:
    lumi = 1 

minX_mass_plot = 1181
maxX_mass_plot = 7320
range_residual = 3
MinNumEvents = 10.
nParFit = 4
xaxisTitle = "Dijet Mass [GeV]"
if showCrossSection==1:
    yaxisTitle_main = "Cross section [pb]"
else:
    yaxisTitle_main = "Number of events / GeV"
yaxisTitle_secondary = "#frac{(Data-Fit)}{#sigma_{Data}}   "
outputLabel = "canvas_dataPlot_m4000"

xVariableWsName = "mjj"
dataWsName="data_obs"
signalWsName="signal"
parameterOfInterest="r"

massBins_list = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]

for ibb,bb in enumerate(massBins_list):
    newbin = int(bb + bb * 0.0) 
    massBins_list[ibb] = newbin
print massBins_list

#==================
# CMS style and lumi
#==================

#set tdr style
setTDRStyle.setTDRStyle()
tdrStyle.SetNdivisions(505, "XYZ")

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "%1.f fb^{-1}" % lumiValue
CMS_lumi.lumi_8TeV = "%1.f fb^{-1}" % lumiValue
CMS_lumi.lumi_13TeV = "%1.f pb^{-1}" % lumiValue
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.cmsTextSize = 1.0

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4


W = 600
H = 700
H_ref = 700 
W_ref = 600 
T = 0.08*H_ref
B = 0.12*H_ref
L = 0.12*W_ref
R = 0.04*W_ref    


#==============================================================================

def main():

    global minX_mass_plot
    global maxX_mass_plot

    #==================
    # Mass binning
    #==================

    workspace = inputFileWorkspace.Get("w")
    print ""
    print ""
    print "======== Workspace Content ========"            
    print ""
    workspace.Print()

    mjj = workspace.var(xVariableWsName)
    minX_mass = mjj.getMin()
    maxX_mass = mjj.getMax()
    print ""
    print ""
    print "======== Fitting Ranges ========"            
    print ""
    print "The fit was performed in this range:"
    print "minX_mass: "+str(minX_mass)
    print "maxX_mass: "+str(maxX_mass)
    if minX_mass_plot<minX_mass: 
        minX_mass_plot=minX_mass
    if maxX_mass_plot>maxX_mass:
        maxX_mass_plot=maxX_mass 

    massBins_list_actual = []
    firstBin = -1
    lastBin = -1
    index_bin_boundary=-1
    for bin_boundary in massBins_list:
        index_bin_boundary = index_bin_boundary+1
        #print bin_boundary, minX_mass_plot, maxX_mass_plot
        if (bin_boundary>=minX_mass_plot and firstBin==-1 ):
            minX_mass_plot = bin_boundary
            firstBin=1
            #print "FIRST BIN is "+str(minX_mass_plot)
        if (bin_boundary>(maxX_mass_plot+0.0000000001) and lastBin==-1 ):
            maxX_mass_plot = massBins_list[index_bin_boundary-1]
            lastBin=1
            #print "LAST BIN is "+str(maxX_mass_plot)
        if (bin_boundary>=minX_mass_plot and bin_boundary<=maxX_mass_plot):
            massBins_list_actual.append(bin_boundary)
    #print massBins_list_actual
    massBins = array("d",massBins_list_actual)
    N_massBins = len(massBins)-1

    print "This is the plotting range:"
    print "minX_mass_plot: "+str(minX_mass_plot)
    print "maxX_mass_plot: "+str(maxX_mass_plot)

    #==================
    # Data spectrum
    #==================

    data_obs = workspace.data(dataWsName)    
    data_obs_TH1_fineBinning = data_obs.createHistogram("data_obs_TH1_fineBinning",mjj)
    data_obs_TH1_fineBinning.Rebin(N_massBins,"data_obs_TH1",massBins)    
    data_obs_TGraph_Events = TGraphAsymmErrors(data_obs_TH1)
    data_obs_TGraph = TGraphAsymmErrors(data_obs_TH1)
    #print "%d %d" % (N_massBins , data_obs_TGraph.GetN())

    print ""
    print ""
    print "======== Number of events ========"            
    print ""
    alpha = 1-0.6827
    for i in range(0,data_obs_TGraph.GetN()):
        N = data_obs_TGraph.GetY()[i]
        binWidth = data_obs_TGraph.GetEXlow()[i] + data_obs_TGraph.GetEXhigh()[i]
        #print str(data_obs_TGraph.GetX()[i])+" "+
        print str(data_obs_TGraph.GetX()[i]-data_obs_TGraph.GetEXlow()[i])+" "+str(data_obs_TGraph.GetX()[i]+data_obs_TGraph.GetEXhigh()[i])+" "+str(N) 

        L = 0
        if N!=0:
            L = ROOT.Math.gamma_quantile(alpha/2,N,1.)
        U = ROOT.Math.gamma_quantile_c(alpha/2,N+1,1)

        data_obs_TGraph_Events.SetPointEYlow(i, (N-L));
        data_obs_TGraph_Events.SetPointEYhigh(i, (U-N));
        data_obs_TGraph_Events.SetPoint(i, data_obs_TGraph.GetX()[i], N)

        data_obs_TGraph.SetPointEYlow(i, (N-L)/(binWidth*lumi));
        data_obs_TGraph.SetPointEYhigh(i, (U-N)/(binWidth*lumi));
        data_obs_TGraph.SetPoint(i, data_obs_TGraph.GetX()[i], N/(binWidth*lumi))

    #==================
    # Background function
    #==================

    fitResults_BFit = inputFileFit.Get("fit_b")
    print ""
    print ""
    print "======== Fit results (B only) ========"            
    print ""
    fitResults_BFit.Print()

    ##-- Start edit --##

    p1_b = fitResults_BFit.floatParsFinal().find("p1") 
    p2_b = fitResults_BFit.floatParsFinal().find("p2") 
    p3_b = fitResults_BFit.floatParsFinal().find("p3") 
    bkgNorm_b = fitResults_BFit.floatParsFinal().find("shapeBkg_background_bin1__norm") 
    #fitResults_BFit.floatParsFinal().Print("s")
    #print "p1_b = " , p1_b.getVal() , " +" , p1_b.getAsymErrorHi() , " -" , p1_b.getAsymErrorLo()
    #print "p2_b = " , p2_b.getVal() , " +" , p2_b.getAsymErrorHi() , " -" , p2_b.getAsymErrorLo()
    #print "p3_b = " , p3_b.getVal() , " +" , p3_b.getAsymErrorHi() , " -" , p3_b.getAsymErrorLo()
    #print "bkgNorm_b = " , bkgNorm_b.getVal() , " +" , bkgNorm_b.getAsymErrorHi() , " -" , bkgNorm_b.getAsymErrorLo()

    background_noNorm = TF1("background_noNorm","( TMath::Power(1-x/13000,[0]) ) / ( TMath::Power(x/13000,[1]+[2]*log(x/13000)) )",minX_mass,maxX_mass)
    background_noNorm.SetParameter(0,p1_b.getVal())
    background_noNorm.SetParameter(1,p2_b.getVal())
    background_noNorm.SetParameter(2,p3_b.getVal())
    norm = background_noNorm.Integral(minX_mass,maxX_mass)
    print "norm : "+str(norm)  
    p0_b = bkgNorm_b.getVal()/(norm*lumi) #FIXME
    print "p0_b = " , bkgNorm_b.getVal()/(norm) , " +" , bkgNorm_b.getAsymErrorHi()/(norm) , " -" , bkgNorm_b.getAsymErrorLo()/(norm)
    print "p1_b = " , p1_b.getVal() , " +" , p1_b.getAsymErrorHi() , " -" , p1_b.getAsymErrorLo()
    print "p2_b = " , p2_b.getVal() , " +" , p2_b.getAsymErrorHi() , " -" , p2_b.getAsymErrorLo()
    print "p3_b = " , p3_b.getVal() , " +" , p3_b.getAsymErrorHi() , " -" , p3_b.getAsymErrorLo()

    background = TF1("background","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
    background.SetParameter(0,p0_b)
    background.SetParameter(1,p1_b.getVal())
    background.SetParameter(2,p2_b.getVal())
    background.SetParameter(3,p3_b.getVal())
    #print "norm (after scaling) : "+str(background.Integral(minX_mass,maxX_mass))

    ##-- End edit --##

    #==================
    # Background function (in S+B fits)
    #==================
    fitResults_SBFit = inputFileFit.Get("fit_s")
    print ""
    print ""
    print "======== Fit results (S+B) ========"            
    print ""
    fitResults_SBFit.Print()

    ##-- Start edit --##

    p1_b_SBFit = fitResults_SBFit.floatParsFinal().find("p1") 
    p2_b_SBFit = fitResults_SBFit.floatParsFinal().find("p2") 
    p3_b_SBFit = fitResults_SBFit.floatParsFinal().find("p3") 
    bkgNorm_b_SBFit = fitResults_SBFit.floatParsFinal().find("shapeBkg_background_bin1__norm") 

    background_noNorm_SBFit = TF1("background_noNorm","( TMath::Power(1-x/13000,[0]) ) / ( TMath::Power(x/13000,[1]+[2]*log(x/13000)) )",minX_mass,maxX_mass)
    background_noNorm_SBFit.SetParameter(0,p1_b_SBFit.getVal())
    background_noNorm_SBFit.SetParameter(1,p2_b_SBFit.getVal())
    background_noNorm_SBFit.SetParameter(2,p3_b_SBFit.getVal())
    norm_SBFit = background_noNorm_SBFit.Integral(minX_mass,maxX_mass)
    #print "norm_SBFit : "+str(norm_SBFit)  
    p0_b_SBFit = bkgNorm_b_SBFit.getVal()/(norm_SBFit*lumi) #FIXME

    background_SBFit = TF1("background_SBFit","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
    background_SBFit.SetParameter(0,p0_b_SBFit)
    background_SBFit.SetParameter(1,p1_b_SBFit.getVal())
    background_SBFit.SetParameter(2,p2_b_SBFit.getVal())
    background_SBFit.SetParameter(3,p3_b_SBFit.getVal())
    #print "norm_SBFit (after scaling) : "+str(background_SBFit.Integral(minX_mass,maxX_mass))

    ##-- End edit --##

    #==================
    # Signal shape (get)
    #==================
    signal_shape = workspace.data(signalWsName)    
    signal_shape_TH1_fineBinning = signal_shape.createHistogram("signal_shape_TH1_fineBinning",mjj)
    signal_shape_TH1_fineBinning.Rebin(N_massBins,"signal_shape_TH1",massBins)    
    integral_signal_shape_TH1 = signal_shape_TH1.Integral()

    #==================
    # Signal yield and errors
    #==================
    print ""
    print ""
    print "======== POI and Signal Yield ========"            
    print ""
    POI = fitResults_SBFit.floatParsFinal().find(parameterOfInterest)
    POI.Print()
    print ("POI : POI.getVal() %f POI.getAsymErrorHi() +%f POI.getAsymErrorLo() -%f" % (POI.getVal(), POI.getAsymErrorHi(), POI.getAsymErrorLo()) )
    print "signal_shape_TH1.Integral() : "+str(integral_signal_shape_TH1)
    Nsig = POI.getVal() * integral_signal_shape_TH1
    eh_Nsig = POI.getAsymErrorHi() * integral_signal_shape_TH1
    el_Nsig = POI.getAsymErrorLo() * integral_signal_shape_TH1
    if el_Nsig==0:
        el_Nsig = -eh_Nsig
    sigmaXaccSig = Nsig/lumiValue
    eh_sigmaXaccSig = eh_Nsig/lumiValue
    el_sigmaXaccSig = el_Nsig/lumiValue
    print ("Number of signal events : %f +%f %f" % (Nsig, eh_Nsig, el_Nsig) ) 
    print ("Cross section: %f +%f %f pb" % (sigmaXaccSig, eh_sigmaXaccSig, el_sigmaXaccSig) )

    #==================
    # Signal shape (rescale)
    #==================    
    for ii in range (0,N_massBins):                        
        if (showCrossSection == 1):#cross section
            value = (signal_shape_TH1.GetBinContent(ii+1)/signal_shape_TH1.GetBinWidth(ii+1))*(sigmaXaccSig/integral_signal_shape_TH1) 
        if (showCrossSection == 0):#numer of events
            value = (signal_shape_TH1.GetBinContent(ii+1)/signal_shape_TH1.GetBinWidth(ii+1))*(Nsig/integral_signal_shape_TH1) 
        signal_shape_TH1.SetBinContent(ii+1,value)
        #print "%d %f" % (ii+1, value)

    #==================
    # Fit Residuals, Chi2, and Make Plots
    #==================
    print ""
    print ""
    print "======== Background Fits ========"            
    background_TH1 = convertFunctionToHisto(background,"background_TH1",N_massBins,massBins)
    ##
    hist_fit_residual_vsMass = TH1D("hist_fit_residual_vsMass","hist_fit_residual_vsMass",N_massBins,massBins)
    list_chi2AndNdf_background = calculateChi2AndFillResiduals(data_obs_TGraph,background_TH1,hist_fit_residual_vsMass,0)
    ##
    drawAndSavePlot_background(data_obs_TGraph,background_TH1,hist_fit_residual_vsMass,outputLabel)

    print ""
    print ""
    print "======== Signal+Background Fits ========"            
    background_TH1_SBFit = convertFunctionToHisto(background_SBFit,"background_TH1_SBFit",N_massBins,massBins)
    signalPlusbackground_TH1_SBFit = convertFunctionToHisto(background_SBFit,"signalPlusbackground_TH1_SBFit",N_massBins,massBins)
    signalPlusbackground_TH1_SBFit.Add(signal_shape_TH1)
    ##
    hist_fit_residual_vsMass_SBFit =  TH1D("hist_fit_residual_vsMass_SBFit","hist_fit_residual_vsMass_SBFit",N_massBins,massBins)
    hist_fit_residual_vsMass_SBFit_withSig =  TH1D("hist_fit_residual_vsMass_SBFit_withSig","hist_fit_residual_vsMass_SBFit_withSig",N_massBins,massBins)
    list_chi2AndNdf_background_SBFit = calculateChi2AndFillResiduals(data_obs_TGraph,background_TH1_SBFit,hist_fit_residual_vsMass_SBFit,0)    
    list_chi2AndNdf_background_SBFit_withSig = calculateChi2AndFillResiduals(data_obs_TGraph,signalPlusbackground_TH1_SBFit,hist_fit_residual_vsMass_SBFit_withSig,0)    
    ##
    drawAndSavePlot_signalPlusBackground(data_obs_TGraph,background_TH1_SBFit,hist_fit_residual_vsMass_SBFit,hist_fit_residual_vsMass_SBFit_withSig,signal_shape_TH1,signalPlusbackground_TH1_SBFit,Nsig,eh_Nsig,el_Nsig,outputLabel)


#============
#end of  main
#============

#==============================================================================
        

def convertFunctionToHisto(background_,name_,N_massBins_,massBins_):
    
    background_hist_ = TH1D(name_,name_,N_massBins_,massBins_)
    
    for bin in range (0,N_massBins_):
        xbinLow = massBins_[bin]
        xbinHigh = massBins_[bin+1]
        binWidth_current = xbinHigh - xbinLow
        value = background_.Integral(xbinLow , xbinHigh) / binWidth_current
        background_hist_.SetBinContent(bin+1,value)

    return background_hist_

def calculateChi2AndFillResiduals(data_obs_TGraph_,background_hist_,hist_fit_residual_vsMass_,prinToScreen_=0):
    
    print "-- "+str(background_hist_.GetName())

    N_massBins_ = data_obs_TGraph_.GetN()

    chi2_FullRangeAll = 0
    chi2_PlotRangeAll = 0
    chi2_PlotRangeNonZero = 0
    chi2_PlotRangeMinNumEvents = 0 

    N_FullRangeAll = 0
    N_PlotRangeAll = 0
    N_PlotRangeNonZero = 0
    N_PlotRangeMinNumEvents = 0 

    if(prinToScreen_):
        print ""
        print ""
        print "======== Number of events / GeV (data, errors, fit, residuals) ========"            
        print ""

    for bin in range (0,N_massBins_):

        ## Values and errors

        value_data = data_obs_TGraph_.GetY()[bin]
        err_low_data = data_obs_TGraph_.GetEYlow()[bin]
        err_high_data = data_obs_TGraph_.GetEYhigh()[bin]
        xbinCenter = data_obs_TGraph_.GetX()[bin] 
        xbinLow = data_obs_TGraph_.GetX()[bin]-data_obs_TGraph_.GetEXlow()[bin] 
        xbinHigh = data_obs_TGraph_.GetX()[bin]+data_obs_TGraph_.GetEXhigh()[bin]
        binWidth_current = xbinHigh - xbinLow
        #value_fit = background_.Integral(xbinLow , xbinHigh) / binWidth_current
        value_fit = background_hist_.GetBinContent(bin+1)
        
        ## Fit residuals

        err_tot_data = 0
        if (value_fit >= value_data):
            err_tot_data = err_high_data  
        else:
            err_tot_data = err_low_data  
        fit_residual = (value_data - value_fit) / err_tot_data
        err_fit_residual = 1

        ## Fill histo with residuals

        hist_fit_residual_vsMass_.SetBinContent(bin+1,fit_residual)
        hist_fit_residual_vsMass_.SetBinError(bin+1,err_fit_residual)

        ## Chi2

        chi2_FullRangeAll += pow(fit_residual,2)
        N_FullRangeAll += 1
        if (xbinLow >= minX_mass_plot and xbinHigh<=maxX_mass_plot):
            chi2_PlotRangeAll += pow(fit_residual,2)
            N_PlotRangeAll += 1
            if (value_data > 0):
                chi2_PlotRangeNonZero += pow(fit_residual,2)
                N_PlotRangeNonZero += 1
                if(value_data * binWidth_current * lumi > MinNumEvents):
                    chi2_PlotRangeMinNumEvents += pow(fit_residual,2)
                    N_PlotRangeMinNumEvents += 1
    
        if(prinToScreen_):
            print str(xbinLow)+" "+str(xbinHigh)+" "+str(binWidth_current)+" : "+str(value_data)+" "+str(value_data * binWidth_current * lumi)+" - "+str(err_low_data)+" + "+str(err_high_data)+" fit: "+str(value_fit)+" fit residual: "+str(fit_residual) 

    #==================
    # Calculate chi2/ndf
    #==================

    # ndf
    ndf_FullRangeAll = N_FullRangeAll - nParFit    
    ndf_PlotRangeAll = N_PlotRangeAll - nParFit    
    ndf_PlotRangeNonZero = N_PlotRangeNonZero - nParFit    
    ndf_PlotRangeMinNumEvents = N_PlotRangeMinNumEvents - nParFit    

    chi2_ndf_FullRangeAll = chi2_FullRangeAll / ndf_FullRangeAll
    chi2_ndf_PlotRangeAll = chi2_PlotRangeAll / ndf_PlotRangeAll
    chi2_ndf_PlotRangeNonZero = chi2_PlotRangeNonZero / ndf_PlotRangeNonZero
    chi2_ndf_PlotRangeMinNumEvents = chi2_PlotRangeMinNumEvents / ndf_PlotRangeMinNumEvents

    print "chi2/ndf FullRangeAll : %.1f / %d = %.2f" % ( chi2_FullRangeAll , ndf_FullRangeAll , chi2_ndf_FullRangeAll ) 
    print "chi2/ndf PlotRangeAll : %.1f / %d = %.2f" % ( chi2_PlotRangeAll , ndf_PlotRangeAll , chi2_ndf_PlotRangeAll ) 
    print "chi2/ndf PlotRangeNonZero : %.1f / %d = %.2f" % ( chi2_PlotRangeNonZero , ndf_PlotRangeNonZero , chi2_ndf_PlotRangeNonZero ) 
    print "chi2/ndf PlotRangeMinNumEvents : %.1f / %d = %.2f" % ( chi2_PlotRangeMinNumEvents , ndf_PlotRangeMinNumEvents , chi2_ndf_PlotRangeMinNumEvents ) 

    return [chi2_FullRangeAll, ndf_FullRangeAll, chi2_PlotRangeAll, ndf_PlotRangeAll, chi2_PlotRangeNonZero, ndf_PlotRangeNonZero, chi2_PlotRangeMinNumEvents, ndf_PlotRangeMinNumEvents]


#==============================================================================
    

def drawAndSavePlot_background(data_obs_TGraph_,background_TH1_,hist_fit_residual_vsMass_,outputLabel_):

    global minY, maxY

    canvas = TCanvas("canvas","canvas",W,H)
    canvas.GetWindowHeight()
    canvas.GetWindowWidth()
    #canvas.SetLogy()
    canvas.SetTitle("")
    canvas.Divide(1,2,0,0,0)

    #pad1 - data spectrum
    canvas.cd(1)
    pad_1 = canvas.GetPad(1)
    #pad_1.SetPad(0.01,0.26,0.99,0.98) #FIXME
    pad_1.SetPad(0.01,0.36,0.99,0.98)
    pad_1.SetLogy()
    pad_1.SetRightMargin(0.05)
    pad_1.SetTopMargin(0.05)
    pad_1.SetFillColor(0)
    pad_1.SetBorderMode(0)
    pad_1.SetFrameFillStyle(0)
    pad_1.SetFrameBorderMode(0)

    if ( (fixedRange==1 and showCrossSection==1)==0 ):
        minY = 0.0001/lumi
        maxY = data_obs_TGraph_.GetY()[0]*10
        
    #vFrame = pad_1.DrawFrame(minX_mass_plot,0.0001/lumi,maxX_mass_plot,data_obs_TGraph_.GetY()[0]*10)
    vFrame = pad_1.DrawFrame(minX_mass_plot,minY,maxX_mass_plot,maxY)
    vFrame.SetTitle("")
    #vFrame.SetXTitle(xaxisTitle)
    #vFrame.SetYTitle(yaxisTitle_main)
    vFrame.GetXaxis().SetTitleSize(0.06)
    vFrame.GetXaxis().SetTitleOffset(0.95)
    vFrame.GetXaxis().SetLabelSize(0.05)
    vFrame.GetYaxis().SetTitleSize(0.06)
    #vFrame.GetYaxis().SetTitleOffset(1.0)
    vFrame.GetYaxis().SetLabelSize(0.05)

    #style data spectrum    
    gStyle.SetErrorX(1)

    data_obs_TGraph_.SetMarkerStyle(20)
    data_obs_TGraph_.SetMarkerColor(1)
    data_obs_TGraph_.SetLineColor(1)
    data_obs_TGraph_.SetTitle("")
    data_obs_TGraph_.GetXaxis().SetTitle(xaxisTitle)
    data_obs_TGraph_.GetYaxis().SetTitle(yaxisTitle_main)
    data_obs_TGraph_.GetXaxis().SetLimits(minX_mass_plot,maxX_mass_plot)
    #data_obs_TGraph_.GetYaxis().SetRangeUser(0.0001/lumi,data_obs_TGraph_.GetY()[0]*10)
    data_obs_TGraph_.GetYaxis().SetRangeUser(minY,maxY)

    #style background function
    background_TH1_.SetLineColor(2)
    background_TH1_.SetLineWidth(2)

    #draw objects
    data_obs_TGraph_.Draw("A P E0")
    background_TH1_.Draw("C SAME")

    #draw text
    pave_general = TPaveText(0.566772,0.794229,0.83557,0.940972,"NDC")    
    pave_general.AddText("background fit")
    pave_general.SetFillColor(0)
    pave_general.SetLineColor(1)
    pave_general.SetFillStyle(0)
    pave_general.SetBorderSize(0)
    pave_general.SetTextFont(42)
    pave_general.SetTextSize(0.040)
    pave_general.SetTextAlign(12) 
    pave_general.SetTextColor(1) 
    pave_general.Draw("SAME")

    #draw text
    pave_sel = TPaveText(0.229489,0.0817972,0.464046,0.254608,"NDC")    
    pave_sel.SetFillColor(0)
    pave_sel.SetBorderSize(0)
    pave_sel.SetFillStyle(0)
    pave_sel.AddText(0.5,1.2,"Wide Jets")
    pave_sel.AddText(0.5,0.5,"m_{jj} > 1.2 TeV")
    pave_sel.AddText(0.5,0.,"|#eta| < 2.5, |#Delta#eta| < 1.3")
    pave_sel.Draw("SAME")

    #draw legend
    leg = TLegend(0.5564991,0.58,0.9203575,0.835812)
    leg.SetTextSize(0.03546853)
    leg.SetLineColor(0)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetMargin(0.35)
    leg.AddEntry(data_obs_TGraph_,"data" ,"EPL")
    leg.AddEntry(background_TH1_,"background","L")
    leg.Draw("SAME")

    #draw pad
    pad_1.RedrawAxis()
    pad_1.Update()
    pad_1.GetFrame().Draw()
    CMS_lumi.CMS_lumi(pad_1, iPeriod, iPos)


    #pad2 - residuals
    canvas.cd(2)
    pad_2 = canvas.GetPad(2)
    #pad_2.SetPad(0.01,0.02,0.99,0.27) #FIXME
    pad_2.SetPad(0.01,0.02,0.99,0.37) 
    pad_2.SetBottomMargin(0.35)
    pad_2.SetRightMargin(0.05)
    pad_2.SetGridx()
    pad_2.SetGridy()

    vFrame2 = pad_2.DrawFrame(pad_1.GetUxmin(), -range_residual, pad_1.GetUxmax(), +range_residual)    
    vFrame2.SetTitle("")
    vFrame2.SetXTitle(xaxisTitle)
    vFrame2.GetXaxis().SetTitleSize(0.06)
    vFrame2.SetYTitle(yaxisTitle_secondary)
    vFrame2.GetYaxis().SetTitleSize(0.12)
    vFrame2.GetYaxis().SetTitleOffset(0.60)
    vFrame2.GetYaxis().SetLabelSize(0.09)
    vFrame2.GetXaxis().SetTitleSize(0.15)
    vFrame2.GetXaxis().SetTitleOffset(0.90)
    vFrame2.GetXaxis().SetLabelSize(0.12)

    #style residuals
    hist_fit_residual_vsMass_.GetXaxis().SetRangeUser(minX_mass_plot,maxX_mass_plot)
    hist_fit_residual_vsMass_.GetYaxis().SetRangeUser(-range_residual,+range_residual)
    hist_fit_residual_vsMass_.SetLineWidth(0)
    hist_fit_residual_vsMass_.SetFillColor(2)
    hist_fit_residual_vsMass_.SetLineColor(1)
    hist_fit_residual_vsMass_.Draw("SAME HIST")

    line = TLine(minX_mass_plot,0,maxX_mass_plot,0)
    line.Draw("")

    #draw pad
    pad_2.RedrawAxis()

    #============

    #write canvas
    canvas.SaveAs(outputLabel_+"_B"+".root")
    canvas.SaveAs(outputLabel_+"_B"+".png")


#==============================================================================


def drawAndSavePlot_signalPlusBackground(data_obs_TGraph_,background_TH1_SBFit_,hist_fit_residual_vsMass_SBFit_,hist_fit_residual_vsMass_SBFit_withSig_,signal_shape_TH1_,signalPlusbackground_TH1_SBFit_,Nsig_,eh_Nsig_,el_Nsig_,outputLabel_):

    global minY, maxY

    canvas = TCanvas("canvas","canvas",W,H)
    canvas.GetWindowHeight()
    canvas.GetWindowWidth()
    #canvas.SetLogy()
    canvas.SetTitle("")
    canvas.Divide(1,2,0,0,0)

    #pad1 - data spectrum
    canvas.cd(1)
    pad_1 = canvas.GetPad(1)
    pad_1.SetPad(0.01,0.36,0.99,0.98)
    pad_1.SetLogy()
    pad_1.SetRightMargin(0.05)
    pad_1.SetTopMargin(0.05)
    pad_1.SetFillColor(0)
    pad_1.SetBorderMode(0)
    pad_1.SetFrameFillStyle(0)
    pad_1.SetFrameBorderMode(0)

    if ( (fixedRange==1 and showCrossSection==1)==0 ):
        minY = 0.0001/lumi
        maxY = data_obs_TGraph_.GetY()[0]*10

    #vFrame = pad_1.DrawFrame(minX_mass_plot,0.0001/lumi,maxX_mass_plot,data_obs_TGraph_.GetY()[0]*10)
    vFrame = pad_1.DrawFrame(minX_mass_plot,minY,maxX_mass_plot,maxY)
    vFrame.SetTitle("")
    #vFrame.SetXTitle(xaxisTitle)
    #vFrame.SetYTitle(yaxisTitle_main)
    vFrame.GetXaxis().SetTitleSize(0.06)
    vFrame.GetXaxis().SetTitleOffset(0.95)
    vFrame.GetXaxis().SetLabelSize(0.05)
    vFrame.GetYaxis().SetTitleSize(0.06)
    #vFrame.GetYaxis().SetTitleOffset(1.0)
    vFrame.GetYaxis().SetLabelSize(0.05)

    #style data spectrum    
    gStyle.SetErrorX(1)

    data_obs_TGraph_.SetMarkerStyle(20)
    data_obs_TGraph_.SetMarkerColor(1)
    data_obs_TGraph_.SetLineColor(1)
    data_obs_TGraph_.SetTitle("")
    data_obs_TGraph_.GetXaxis().SetTitle(xaxisTitle)
    data_obs_TGraph_.GetYaxis().SetTitle(yaxisTitle_main)
    data_obs_TGraph_.GetXaxis().SetLimits(minX_mass_plot,maxX_mass_plot)
    #data_obs_TGraph_.GetYaxis().SetRangeUser(0.0001/lumi,data_obs_TGraph_.GetY()[0]*10)
    data_obs_TGraph_.GetYaxis().SetRangeUser(minY,maxY)

    #style background function
    background_TH1_SBFit_.SetLineColor(2)
    background_TH1_SBFit_.SetLineWidth(2)

    #create TGraph for signal
    signal_shape_TGraph_ = TGraphAsymmErrors(signal_shape_TH1_)
    signal_shape_TGraph_cut_ = TGraph()    
    for jj in range (0,signal_shape_TGraph_.GetN()):
        ## Values and errors
        value_signal = signal_shape_TGraph_.GetY()[jj]
        if value_signal>0:
            signal_shape_TGraph_cut_.SetPoint(jj,signal_shape_TGraph_.GetX()[jj],value_signal)

    #style signal
    signal_shape_TGraph_cut_.SetLineColor(4)
    signal_shape_TGraph_cut_.SetLineWidth(2)

    #style signal+background function
    signalPlusbackground_TH1_SBFit_.SetLineColor(4)
    signalPlusbackground_TH1_SBFit_.SetLineWidth(2)

    #draw objects
    data_obs_TGraph_.Draw("A P E0")
    signalPlusbackground_TH1_SBFit_.Draw("C SAME HIST")
    background_TH1_SBFit_.Draw("C SAME")
    if (drawSignalShapeAlsoAlone==1):
        signal_shape_TGraph_cut_.Draw("C SAME")

    #draw text
    pave_general = TPaveText(0.566772,0.794229,0.83557,0.940972,"NDC")    
    pave_general.AddText("signal+background fit")
    pave_general.SetFillColor(0)
    pave_general.SetLineColor(1)
    pave_general.SetFillStyle(0)
    pave_general.SetBorderSize(0)
    pave_general.SetTextFont(42)
    pave_general.SetTextSize(0.040)
    pave_general.SetTextAlign(12) 
    pave_general.SetTextColor(1) 
    pave_general.Draw("SAME")

    #draw text
    pave_sel = TPaveText(0.674634,0.361816,0.909191,0.533825,"NDC")    
    pave_sel.SetFillColor(0)
    pave_sel.SetBorderSize(0)
    pave_sel.SetFillStyle(0)
    pave_sel.AddText(0.5,1.2,"Wide Jets")
    pave_sel.AddText(0.5,0.5,"m_{jj} > 1.2 TeV")
    pave_sel.AddText(0.5,0.,"|#eta| < 2.5, |#Delta#eta| < 1.3")
    pave_sel.Draw("SAME")

    #draw legend
    leg = TLegend(0.5564991,0.58,0.9203575,0.835812)
    leg.SetTextSize(0.03546853)
    leg.SetLineColor(0)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetMargin(0.35)
    leg.AddEntry(data_obs_TGraph_,"data" ,"EPL")
    leg.AddEntry(background_TH1_SBFit_,"background","L")
    text_signal = "signal (N_{signal} = %2.f^{+%2.f}_{-%2.f})" % (Nsig_ , eh_Nsig_ , el_Nsig_) 
    leg.AddEntry(signal_shape_TGraph_cut_,text_signal,"L")
    leg.Draw("SAME")

    #draw pad
    pad_1.RedrawAxis()
    pad_1.Update()
    pad_1.GetFrame().Draw()
    CMS_lumi.CMS_lumi(pad_1, iPeriod, iPos)


    #pad2 - residuals
    canvas.cd(2)
    pad_2 = canvas.GetPad(2)
    pad_2.SetPad(0.01,0.02,0.99,0.37)
    pad_2.SetBottomMargin(0.35)
    pad_2.SetRightMargin(0.05)
    pad_2.SetGridx()
    pad_2.SetGridy()

    vFrame2 = pad_2.DrawFrame(pad_1.GetUxmin(), -range_residual, pad_1.GetUxmax(), +range_residual)    
    vFrame2.SetTitle("")
    vFrame2.SetXTitle(xaxisTitle)
    vFrame2.GetXaxis().SetTitleSize(0.06)
    vFrame2.SetYTitle(yaxisTitle_secondary)
    vFrame2.GetYaxis().SetTitleSize(0.12)
    vFrame2.GetYaxis().SetTitleOffset(0.60)
    vFrame2.GetYaxis().SetLabelSize(0.09)
    vFrame2.GetXaxis().SetTitleSize(0.15)
    vFrame2.GetXaxis().SetTitleOffset(0.90)
    vFrame2.GetXaxis().SetLabelSize(0.12)

    #style residuals
    hist_fit_residual_vsMass_SBFit_.GetXaxis().SetRangeUser(minX_mass_plot,maxX_mass_plot)
    hist_fit_residual_vsMass_SBFit_.GetYaxis().SetRangeUser(-range_residual,+range_residual)
    hist_fit_residual_vsMass_SBFit_.SetLineWidth(0)
    hist_fit_residual_vsMass_SBFit_.SetFillColor(2)
    hist_fit_residual_vsMass_SBFit_.SetLineColor(1)
    hist_fit_residual_vsMass_SBFit_.Draw("SAME HIST")
    ##
    hist_fit_residual_vsMass_SBFit_withSig_.SetLineWidth(0)
    hist_fit_residual_vsMass_SBFit_withSig_.SetLineColor(4)
    hist_fit_residual_vsMass_SBFit_withSig_.SetLineStyle(2)
    hist_fit_residual_vsMass_SBFit_withSig_.Draw("SAME HIST")

    line = TLine(minX_mass_plot,0,maxX_mass_plot,0)
    line.Draw("")

    #draw pad
    pad_2.RedrawAxis()

    #============

    #write canvas
    canvas.SaveAs(outputLabel_+"_SB"+".root")
    canvas.SaveAs(outputLabel_+"_SB"+".png")

#==============================================================================


if __name__ == '__main__':
    main()


