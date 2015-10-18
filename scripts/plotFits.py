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

parser.add_argument("-i", "--inputList", type=str, dest="inputList", default="",
    help="list of datacard processed to extract name of the samples",
    )
parser.add_argument("--inputFitRes", type=str, dest="inputFitRes", default="",
    help="directory containing files with results of fit",
    )
parser.add_argument("--inputToys", type=str, dest="inputToys", default="",
    help="directory containing files with toys",
    )
parser.add_argument("-o", "--output", type=str, dest="output", default="",
    help="the directory OUTDIR contains the output of the program",
    metavar="OUTDIR"
    )
parser.add_argument("--tag", type=str, dest="tag", default="",
    help="tag contained in the name of files",
    )
parser.add_argument("--mu", type=float, dest="mu", default=0,
    help="injected signal strenght",
    metavar="MU"
    )
parser.add_argument("-t", "--toys", type=int, dest="toys", default=1,
    help="toys to plot",
    metavar="TOYS"
    )

  
lumi=1000.
#lumi=830.

args = parser.parse_args()
print args
############
ROOT.gStyle.SetOptFit(0000)
gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

gStyle.SetOptStat(0)

ins = open(args.inputList,"r")
  #giulia debug for significance experiment

#inputToys = TFile("blindExercise_MLfit.MaxLikelihoodFit.root")
#inputToys = TFile("higgsCombineQstar5000_MLfit_mu_limit.MaxLikelihoodFit.mH120.123456.root")


massBins_list = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]
#massBins_list = [ 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099]
massBins_list_limited = [ 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564]

massBins = array("d",massBins_list)
massBins_limited = array("d",massBins_list_limited)



#giulia : TO FIX
dict = {
    #2000:1.26511629999999986e+01,
    #3000:1.00099999999999989e+00,
    #4000:1.08433599999999991e-01,
    #5000:1.11493200000000008e-02,
    #6000:1.86000000000000011e-03
    #1200.0:  0.1721E+03,
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
    6000.0:  0.1159E-02
    }

bins = 50

k=0
##this loop is useless, to remove
for line in  ins:
  if not k==0:
    continue
  print "====================="
  print "    %d     " % k
  print "====================="
  sample = os.path.basename(line)
  dir = os.path.dirname(line) 
  #sample = os.path.splitext(line)[0]
  sample = sample.split("_")[0]
  print ("process %s" % sample)
  line = line.rstrip('\n')
  sample = sample.rstrip('\n')
  mass_string = sample.split("Qstar")[1]
  mass = int(mass_string)
  #print mass_string
  #print line
  input_workspace = TFile(dir +"/"+ sample +"_workspace_toFit.root")
  print input_workspace
  workspace = input_workspace.Get("w")
  workspace.Print()
  signal=workspace.pdf("signal")
  datahist_sig = signal.dataHist()
  set_sig=datahist_sig.get()
  mjj_sig=set_sig.find("mjj")
  h_sig_1GeV= signal.createHistogram("h_sig_1GeV",mjj_sig)
  h_sig_1GeV.Print()

  #h_sig.Rebin(100)
  #datahist_sig_rebin=RooDataHist("datahist_sig_rebin","datahist_signal_rebin",RooArgList(mjj_sig),h_sig)
  #signal_rebin = RooHistPdf("signal_rebin","signal_rebin",RooArgSet(mjj_sig),datahist_sig_rebin)

  if not (args.mu == -999): 
    inputToys = TFile(args.inputToys+"/higgsCombine"+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+".GenerateOnly.mH120.123456.root") 
    inputFitRes = TFile(args.inputFitRes+"/mlfit"+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+".root")
    #inputToys = TFile("output_toys_mu"+str(int(args.mu))+"_flatParam/higgsCombine"+sample+"_MLfit_mu"+str(int(args.mu))+".GenerateOnly.mH120.123456.root") 
    #inputFitRes = TFile("output_toys_mu"+str(int(args.mu))+"_flatParam/mlfit"+sample+"_MLfit_mu"+str(int(args.mu))+".root")
    #inputToys = TFile("output_toys_mu"+str(int(args.mu))+"_flatParam_genFromHist/higgsCombine"+sample+"_MLfit_mu"+str(int(args.mu))+".GenerateOnly.mH120.123456.root") 
    #inputFitRes = TFile("output_toys_mu"+str(int(args.mu))+"_flatParam_genFromHist/mlfit"+sample+"_MLfit_mu"+str(int(args.mu))+".root")
  else:
    inputToys = TFile(args.inputToys+"/higgsCombine"+sample+args.tag+"_MLfit_mu_limit.GenerateOnly.mH120.123456.root") 
    inputFitRes = TFile(args.inputFitRes+"/mlfit"+sample+args.tag+"_MLfit_mu_limit.root")
    #inputToys = TFile("output_toys_mu_limit_2Sband_flatParam/higgsCombine"+sample+"_MLfit_mu_limit.GenerateOnly.mH120.123456.root") 
    #inputFitRes = TFile("output_toys_mu_limit_2Sband_flatParam/mlfit"+sample+"_MLfit_mu_limit.root")
    #inputToys = TFile("output_toys_mu_limit_flatParam_genFromHist/higgsCombine"+sample+"_MLfit_mu_limit.GenerateOnly.mH120.123456.root") 
    #inputFitRes = TFile("output_toys_mu_limit_flatParam_genFromHist/mlfit"+sample+"_MLfit_mu_limit.root")

  tree_fit_sb = inputFitRes.Get("tree_fit_sb")
  tree_fit_b = inputFitRes.Get("tree_fit_b")
  #tree_fit_sb.Print()
  #tree_fit_b.Print()

  #print inputToys
  j = 0
  list_toys = []
  chi2 = []
  h_chi2 = TH1F("h_chi2","",80,0,80 )
  #inputToys.cd()
  c = TCanvas("c","", 1)
  c_pulls = TCanvas("c_pulls","",1)
  print "now loop on toys"
  for j in range(1,args.toys+1):
  
    name_toy = "toy_"+str(j)
    toy =  inputToys.Get("toys/"+name_toy)
    print "toys/"+name_toy
    list_toys.append(toy)
    print toy
    set_ = toy.get()
    print set_
    set = RooArgSet(set_)
    x = set.find("mjj")
    #x.setRange(minX_mass,maxX_mass)

    minX_mass=x.getMin() 
    maxX_mass =x.getMax() 

    ###
    tree_fit_sb.GetEntry(j-1)
    p1_val = tree_fit_sb.p1
    p2_val = tree_fit_sb.p2
    p3_val = tree_fit_sb.p3
    mu_val = tree_fit_sb.mu
    integral_val = tree_fit_sb.n_exp_binbin1_proc_background  
    tree_fit_b.GetEntry(j-1)
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
    print "integral_val: "+str(integral_val_b)+"   norm : "+str(norm_b)  
    
    background_noNorm_SplusB = TF1("background_SplusB","( TMath::Power(1-x/13000,[0]) ) / ( TMath::Power(x/13000,[1]+[2]*log(x/13000)) )",minX_mass,maxX_mass)
    background_noNorm_SplusB.SetParameter(0,p1_val)
    background_noNorm_SplusB.SetParameter(1,p2_val)
    background_noNorm_SplusB.SetParameter(2,p3_val)
    norm = background_noNorm_SplusB.Integral(minX_mass,maxX_mass)
    p0_val = integral_val/norm
    
    h_sig_1GeV.Scale(dict[mass]*mu_val*lumi/h_sig_1GeV.Integral())
    h_sig_rebin=h_sig_1GeV.Rebin(len(massBins)-1,"h_sig_rebin",massBins)
    h_sig=TH1F("h_sig","",len(massBins)-1,massBins)
    for i in range(1,len(massBins)):
      h_sig.SetBinContent(i, h_sig_rebin.GetBinContent(i)/h_sig_rebin.GetBinWidth(i))
 
    
    Nobs_val = toy.sumEntries()
    rateSig_val = lumi*dict[mass]
    print str(p0_val)+"   "+ str(p1_val)+"  "+str(p2_val)+"  "+str(p3_val)+"  "+str(mu_val)
    print "rateSig_val = "+str(Nobs_val)+" * "+str(dict[mass])+" =  "+str(rateSig_val) 

    p1 = RooRealVar('p1','p1',p1_val,-1000,1000.)
    p2 = RooRealVar('p2','p2',p2_val,-1000.,1000.)
    p3 = RooRealVar('p3','p3',p3_val,-1000.,1000)
    r = RooRealVar('r','r',mu_val,-2,100)
    Nobs = RooRealVar('Nobs','Nobs',Nobs_val,0,10000000)
    rateSig = RooRealVar('rateSig','rateSig',rateSig_val,-100000,10000000)
    p1.Print()
    p2.Print()
    p3.Print()
    r.Print()
    Nobs.Print() 
    rateSig.Print()

#### Draw plots using Roofit ######
    #background = RooGenericPdf("background","background","(pow(1-@0/13000,@1)/pow(@0/13000,@2+@3*log(@0/13000)))",RooArgList(x,p1,p2,p3))
    #f_sb = RooFormulaVar('f_sb','f_sb','@0/@1*@2',RooArgList(rateSig,Nobs,r))
    #f_sb.Print()
    #bkg_sig = RooAddPdf('bkg_sig','bkg_sig',signal_rebin,background,f_sb)
  
#    frame = x.frame()
#    print "toy "+str(j)+" entries : "+str(toy.sumEntries())
#    toy.plotOn(frame,RooFit.Binning(bins,minX_mass,maxX_mass))
#    background.plotOn(frame,RooFit.LineColor(kRed),RooFit.LineStyle(1),RooFit.LineWidth(2))#,RooFit.Normalization(p0_val*(maxX_mass-minX_mass)/bins))
#  
###### giulia: for the moment don't plot the signal #######  
#    if not (args.mu==0):
#      #to change the normalization according to mu
#      #signal_rebin.plotOn(frame,RooFit.LineColor(kBlue),RooFit.LineStyle(2),RooFit.LineWidth(2), RooFit.Normalization(rateSig_val/Nobs_val))
#      bkg_sig.plotOn(frame,RooFit.LineColor(kBlue),RooFit.LineStyle(2),RooFit.LineWidth(2))
###################################
#    frame.GetYaxis().SetRangeUser(0.1,1000000)
#    c.cd()
#    c.SetLogy(1)
#    frame.Draw()
#    c.SetTitle(name_toy)
#    if not (j>20):
#      c.SaveAs(args.output+"/fit_"+sample+"_"+name_toy+".png")
#      c.SaveAs(args.output+"/fit_"+sample+"_"+name_toy+".pdf")
#      c.SaveAs(args.output+"/pseudodataset_"+name_toy+".png")
#      c.SaveAs(args.output+"/pseudodataset_"+name_toy+".pdf")
#      c.Clear()
#    
#    #pulls
#    print "chi2 = "+str(frame.chiSquare())
#    chi2.append(frame.chiSquare())
#    h_chi2.Fill(frame.chiSquare())
#
#    c_pulls.cd()
#    hpull = frame.pullHist() 
#    frame_pull = x.frame(RooFit.Title("Pull Distribution")) 
#    frame_pull.addPlotable(hpull,"P")
#    frame_pull.GetYaxis().SetRangeUser(-4,4)
#    frame_pull.Draw()
#    if not (j>20):
#      c_pulls.SaveAs(args.output+"/fitResiduals_"+sample+"_"+name_toy+".png")
#      c_pulls.SaveAs(args.output+"/fitResiduals_"+sample+"_"+name_toy+".pdf")
#    ## stop to the num of toys specified in the options

##### Move from RooFit -> Root object #########
    ## Pseudo Data ##
    bin1GeV = RooBinning(minX_mass,maxX_mass)
    bin1GeV.addUniform(int(maxX_mass-minX_mass),minX_mass,maxX_mass)
    x.setBinning(bin1GeV)
    toy_1GeV = toy.binnedClone()
    
    h_toy_1GeV = toy_1GeV.createHistogram("h_toy_1GeV_",x)
    ## to have correct error
    #for jj in range(1,int(maxX_mass-minX_mass)):
    #  if h_toy_1GeV.GetBinContent==0:
    #	h_toy_1GeV.SetBinError(jj,1.8)
    #  else:
    #    h_toy_1GeV.SetBinError(jj,TMath.Sqrt(h_toy_1GeV.GetBinContent(jj)))
      
    h_toy_binned = h_toy_1GeV.Rebin(len(massBins)-1,"h_toy_binned",massBins)
    for ibin in range(1,len(massBins)):
      if (h_toy_binned.GetBinContent(ibin) == 0):
    	h_toy_binned.SetBinError(ibin, 1.8)
      else:
    	h_toy_binned.SetBinError(ibin, TMath.Sqrt(h_toy_binned.GetBinContent(ibin)))
    
    h_toy = TH1F("h_toy","",len(massBins),massBins) 
    for ibin in range(1,len(massBins)):
      h_toy.SetBinContent(ibin, h_toy_binned.GetBinContent(ibin)/h_toy_binned.GetBinWidth(ibin))
      if (h_toy.GetBinContent(ibin) == 0):
	h_toy.SetBinError(ibin, 1.8/h_toy.GetBinWidth(ibin))
      else:
	h_toy.SetBinError(ibin, h_toy_binned.GetBinError(ibin)/h_toy_binned.GetBinWidth(ibin))
      print "content %d: %.2f   error: %.2f" %(ibin,h_toy.GetBinContent(ibin),h_toy.GetBinError(ibin))

    #print "min and max : %d  %d" % (minX_mass, maxX_mass)
    #print "debug: bin content  : %.2f  %.2f  %.2f" % (h_toy_1GeV.GetBinContent(50), h_toy_binned.GetBinContent(50), h_toy.GetBinContent(50))   
    #print "debug: bin width  : %.2f  %.2f  %.2f" % (h_toy_1GeV.GetBinWidth(50), h_toy_binned.GetBinWidth(50), h_toy.GetBinWidth(50))   
       
    #c_test = TCanvas("c_test","",800,800)
    #c_test.SetLogy(1)
    #h_toy_1GeV.Draw()
    #c_test.SaveAs("test_toy1.png")

    ## background ##
    background = TF1("background","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
    background.SetParameter(0,p0_val_b)
    background.SetParameter(1,p1_val_b)
    background.SetParameter(2,p2_val_b)
    background.SetParameter(3,p3_val_b)
    background.SetLineColor(kRed)
    ##background component of S+B
    background_SplusB = TF1("background_SplusB","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
    background_SplusB.SetParameter(0,p0_val)
    background_SplusB.SetParameter(1,p1_val)
    background_SplusB.SetParameter(2,p2_val)
    background_SplusB.SetParameter(3,p3_val)
 
    ## signal plus background ##
    h_s_plus_b = TH1F("h_s_plus_b","",len(massBins_limited)-1,massBins_limited) 
    for ibin in range(1,len(massBins_limited)):
      h_s_plus_b.SetBinContent(ibin, 0)
      bin_center=h_s_plus_b.GetBinCenter(ibin)
      if(h_s_plus_b.GetBinCenter(ibin) > minX_mass and h_s_plus_b.GetBinCenter(ibin) < maxX_mass):
        print "b = "+str(background_SplusB.Eval(bin_center))+"   s = "+str(h_sig.GetBinContent( h_sig.FindBin(bin_center)))+"   s+b = "+str(background_SplusB.Eval(bin_center) + h_sig.GetBinContent( h_sig.FindBin(bin_center)))
        h_s_plus_b.SetBinContent(ibin, background_SplusB.Eval(h_s_plus_b.GetBinCenter(ibin)) + h_sig.GetBinContent( h_sig.FindBin(bin_center)))

    h_s_plus_b.Print()


    ###########################################
    # fit residuals and chi2 variable binning
    ###########################################
    hist_fit_residual_vsMass = TH1D("hist_fit_residual_vsMass","hist_fit_residual_vsMass",len(massBins)-1,massBins)
    hist_fit_residual_vsMass_bkg = TH1D("hist_fit_residual_vsMass_bkg","hist_fit_residual_vsMass_bkg",len(massBins)-1,massBins)
    NumberOfObservations_VarBin = 0
    chi2_VarBin = 0.
    
    for bin in range(1,len(massBins)):
        
      if( h_toy.GetXaxis().GetBinLowEdge(bin)>=minX_mass  and h_toy.GetXaxis().GetBinUpEdge(bin)<=maxX_mass ):
        NumberOfObservations_VarBin+=1
        data = h_toy.GetBinContent(bin)
        err_data = h_toy.GetBinError(bin)
        if( data == 0 ):
          err_data = 1.8 / h_toy.GetBinWidth(bin)
          print "err_data %f" % err_data
        #background only
	fit_B = background.Integral(h_toy.GetXaxis().GetBinLowEdge(bin), h_toy.GetXaxis().GetBinUpEdge(bin) ) 
        fit_B = fit_B / ( h_toy.GetBinWidth(bin) )
        #S+B
	fit = h_s_plus_b.GetBinContent(h_s_plus_b.FindBin(h_toy.GetBinCenter(bin) )) 
        err_tot = err_data	  
        fit_residual = (data - fit) / err_tot
        fit_residual_B = (data - fit_B) / err_tot
    	    	  
        chi2_VarBin += pow( (data - fit) , 2 ) / pow( err_data , 2 )
        print "data, err_data, fit: "+str( data)+", "+str(err_data) +", " +str(fit)
        print "bin, fit residual : " + str(bin) + ", " +str(fit_residual)	  
        hist_fit_residual_vsMass.SetBinContent(bin,fit_residual)
        hist_fit_residual_vsMass_bkg.SetBinContent(bin,fit_residual_B)
       
    ##END OF LOOP over bins
    
    ndf_VarBin = NumberOfObservations_VarBin - 4 
    
    h_chi2.Fill(chi2_VarBin,1)
    if j < 10:
      print "============================"
      print "NumberOfObservations_VarBin: " + str(NumberOfObservations_VarBin)
      print "ndf_VarBin: " + str(ndf_VarBin)
      print "chi2_VarBin: " +str(chi2_VarBin)
      print "============================"   
      
   


    ###### Draw note style
    c2 = TCanvas("c2", "",339,117,695,841)
    c2.GetWindowHeight() 
    c2.GetWindowWidth() 
    c2.SetLogy() 
    c2.Divide(1,2,0,0,0) 
    c2.cd(1) 
    p11_1 = c2.GetPad(1) 
    p11_1.SetPad(0.01,0.23,0.99,0.98) 
    p11_1.SetLogy() 
    p11_1.SetRightMargin(0.05) 
    p11_1.SetTopMargin(0.05) 

    #Pave text for fit results
    ndf = 50 
    pave_fit1 = TPaveText(0.56,0.60,0.9,0.9,"NDC") 
    #pave_fit1.AddText("#chi^{2} / ndf = %.2f" % (frame.chiSquare())) 
    pave_fit1.AddText("#chi^{2} / ndf = %.2f / %d" % (chi2_VarBin, ndf_VarBin)) 
    pave_fit1.AddText("p0 =  %.2e " % (p0_val))
    pave_fit1.AddText("p1 =  %.2e " % (p1_val))
    pave_fit1.AddText("p2 =  %.2e " % (p2_val))
    pave_fit1.AddText("p3 =  %.2e " % (p3_val))
    pave_fit1.AddText("mu =  %.2e " % (mu_val))
    pave_fit1.SetFillColor(0) 
    pave_fit1.SetFillStyle(0) 
    pave_fit1.SetBorderSize(1) 
    pave_fit1.SetTextFont(42) 
    pave_fit1.SetTextSize(0.03) 
    pave_fit1.SetTextAlign(12)  
    
    #Pave text
    pave_fit = TPaveText(0.23,0.15,0.45,0.27,"NDC") 
    #pave_fit.AddText(" #sqrt{s} = 13 TeV") 
    pave_fit.AddText("|#eta| < 2.5, |#Delta#eta| < 1.3") 
    pave_fit.AddText("M_{jj} > 1.1 TeV") 
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

    vFrame = p11_1.DrawFrame(minX_mass,0.0001,maxX_mass,20000.0) 
    
    vFrame.SetTitle("") 
    vFrame.SetXTitle("Dijet Mass [GeV]") 
    vFrame.GetXaxis().SetTitleSize(0.06) 
    vFrame.SetYTitle("Events / bin width") 
    
    vFrame.GetYaxis().SetTitleSize(0.12) 
    vFrame.GetYaxis().SetLabelSize(0.07) 
    vFrame.GetYaxis().SetTitleOffset(0.50) 
    vFrame.GetXaxis().SetTitleOffset(0.90) 
    vFrame.GetXaxis().SetTitleSize(0.18) 
    vFrame.GetXaxis().SetLabelSize(0.1) 
  
    leg = TLegend(0.5564991,0.45,0.8903575,0.58) 
    leg.SetTextSize(0.03146853) 
    leg.SetLineColor(0) 
    leg.SetLineStyle(1) 
    leg.SetLineWidth(0) 
    leg.SetFillStyle(0) 
    leg.SetMargin(0.35) 
    leg.AddEntry(h_toy,"Data" ,"PL") 
    leg.AddEntry(background,"B Fit","L") 
    leg.AddEntry(h_s_plus_b,"S+B Fit","L") 

    #frame.Draw()
    h_toy.GetXaxis().SetRangeUser(minX_mass,maxX_mass) 
    h_toy.SetTitle("") 
    h_toy.SetLineColor(1) 
    h_toy.SetFillColor(1) 
    h_toy.SetLineColor(1) 
    h_toy.SetMarkerColor(1) 
    h_toy.SetMarkerStyle(20) 
    h_toy.SetMinimum(0.001) 

    background.SetLineWidth(2) 
    background.SetLineStyle(2) 
    background.SetLineColor(2) 
    background.Draw("")
    h_s_plus_b.SetLineWidth(2) 
    h_s_plus_b.SetLineStyle(1) 
    h_s_plus_b.SetLineColor(kBlue) 
    h_s_plus_b.Draw("c same")
    h_toy.Draw("pe same")
    pt1.Draw("same")  
    pt2.Draw("same") 
    pt3.Draw("same") 
    pave_fit.Draw("same") 
    pave_fit1.Draw("same") 
    leg.Draw("same") 

 
    #---- Next PAD
    
    c2.cd(2) 
    p11_2 = c2.GetPad(2) 
    p11_2.SetPad(0.01,0.02,0.99,0.24) 
    p11_2.SetBottomMargin(0.35) 
    p11_2.SetRightMargin(0.05) 
    p11_2.SetGridx() 
    p11_2.SetGridy() 
    #c3_2.SetTickx(50) 
    
    vFrame2 = p11_2.DrawFrame(p11_1.GetUxmin(), -4.5, p11_1.GetUxmax(), 4.5) 
    
    vFrame2.SetTitle("") 
    vFrame2.SetXTitle("Dijet Mass [GeV]") 
    vFrame2.GetXaxis().SetTitleSize(0.06) 
    vFrame2.SetYTitle("(Data-Fit)/#sigma_{Data}") 
    vFrame2.GetYaxis().SetTitleSize(0.12) 
    vFrame2.GetYaxis().SetLabelSize(0.07) 
    vFrame2.GetYaxis().SetTitleOffset(0.50) 
    vFrame2.GetXaxis().SetTitleOffset(0.90) 
    vFrame2.GetXaxis().SetTitleSize(0.18) 
    vFrame2.GetXaxis().SetLabelSize(0.1) 
    #frame_pull.Draw()
    hist_fit_residual_vsMass.GetXaxis().SetRangeUser(minX_mass,maxX_mass) 
    hist_fit_residual_vsMass.GetYaxis().SetRangeUser(-4.,4.) 
    hist_fit_residual_vsMass.SetLineWidth(0) 
    hist_fit_residual_vsMass.SetFillColor(kBlue) 
    hist_fit_residual_vsMass.SetLineColor(1) 
    hist_fit_residual_vsMass.Draw("hist same")
    hist_fit_residual_vsMass_bkg.SetLineWidth(2) 
    hist_fit_residual_vsMass_bkg.SetLineColor(2) 
    hist_fit_residual_vsMass_bkg.SetLineStyle(2) 
    hist_fit_residual_vsMass_bkg.Draw("hist same")

    line =  TLine(minX_mass,0,maxX_mass,0) 
    line.Draw("") 

    if j<20:
      c2.SaveAs(args.output+"/fitAndResiduals_SplusB_"+sample+"_"+name_toy+".png")
      c2.SaveAs(args.output+"/fitAndResiduals_SplusB_"+sample+"_"+name_toy+".pdf")

    j+=1
  ## END OF LOOP OVER TOYS
  
  
  
#  canvas_chi2= TCanvas("chi2_VarBin","#chi^{2} distribution of bkg fit",600,600)
#  canvas_chi2.cd()
#  #gStyle.SetOptStat("nemr")
#  ##chi2 nominal function 
#  chi2= TF1("chi2","ROOT::Math::chisquared_pdf(x,35,0)",0,100)
#  chi2.SetLineWidth(2)
#  chi2.SetLineColor(2)
#  h_chi2.SetStats(1)
#  h_chi2.DrawNormalized("hist")
#  chi2.Draw("same")
#
#  leg_chi2 = TLegend(0.6,0.6,0.86,0.75)
#  leg_chi2.SetFillColor(0)
#  leg_chi2.AddEntry(h_chi2,"toy distribution","l")
#  leg_chi2.AddEntry(chi2,"#chi^{2} function for 35 ndf","l")
#  leg_chi2.Draw("")
#
#  canvas_chi2.SaveAs("chi2_VarBin.png")
#  canvas_chi2.SaveAs("chi2_VarBin.pdf")

#  c_test=TCanvas("test","",600,600)
#  c_test.cd()
#  h_s_plus_b.Draw("c")
#  c_test.SaveAs("test_sig.png")

  k+=1
# END OF LOOP OVER MASSES


