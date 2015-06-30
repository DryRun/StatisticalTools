#!usr/bin/python
import sys, getopt
import ROOT
from ROOT import TFile, TH1F, TCanvas, TPad, TMath
from ROOT import gROOT, gPad 
from ROOT import RooRealVar, RooDataHist, RooPlot, RooArgList, RooArgSet,  RooBernstein, RooCBShape, RooAddPdf, RooFit, RooGenericPdf, RooWorkspace, RooMsgService, RooHistPdf, RooBinning, RooUniformBinning, RooExtendPdf
from setTDRStyle import setTDRStyle

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-s","--fitSig",action="store_true",default=False,dest="fitSig")
parser.add_option("-H","--histpdfSig",action="store_true",default=False,dest="histpdfSig")
parser.add_option("--histpdfBkg",action="store_true",default=False,dest="histpdfBkg")
parser.add_option("-d","--fitDat",action="store_true",default=False,dest="fitDat")
parser.add_option("-m","--mass",action="store",type="int",dest="mass",default=5000)
parser.add_option("-o","--path",action="store",type="string",dest="outdir_datacards")
#parser.add_option("-t","--toy",action="store",type="int",dest="TOY", default=1)
parser.add_option("-n","--name",action="store",type="string",dest="name", default="")
parser.add_option("--inputBkg",action="store",type="string",dest="inputBkg", default="")
parser.add_option("--inputData",action="store",type="string",dest="inputData", default="")
parser.add_option("--inputSig",action="store",type="string",dest="inputSig", default="")


parser.add_option("--lumi",action="store",type="float",dest="lumi",default=1000.)
parser.add_option("--sigEff",action="store",type="float",dest="sigEff",default=0.629)
parser.add_option("--sigXS",action="store",type="float",dest="sigXS",default=0.0182)
parser.add_option("--bkgConst",action="store_true",dest="bkgConst",default=False)
parser.add_option("--bkgNuisance",action="store_true",dest="bkgNuisance",default=False)



(options, args) = parser.parse_args()

mass = options.mass
fitSig = options.fitSig
histpdfSig = options.histpdfSig
histpdfBkg = options.histpdfBkg
fitDat = options.fitDat
outdir_datacards = options.outdir_datacards
bkgConst = options.bkgConst
bkgNuisance = options.bkgNuisance
#TOY = options.TOY
name = options.name
filenameBkg = options.inputBkg
filenameDat = options.inputData
filenameSig = options.inputSig

######## some options are incompatible, to avoid errors ###
if bkgConst: bkgNuisance=False
if bkgNuisance: bkgConst=False
if histpdfSig: fitSig=False
if histpdfBkg: fitDat=False
#############


gROOT.SetBatch(ROOT.kTRUE)
gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

RooMsgService.instance().setSilentMode(ROOT.kTRUE)
RooMsgService.instance().setStreamStatus(0,ROOT.kFALSE)
RooMsgService.instance().setStreamStatus(1,ROOT.kFALSE)

# -----------------------------------------
# get histograms
## ---- CERN -------
#PATH = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/output/'
#
#filenameSig = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeMaker/test/Resonance_Shapes_qg_PU20_13TeV_newJEC.root'
##filenameSig = PATH+'rootfile_QstarToJJ_M_'+str(mass)+'_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1__MINIAODSIM_reduced_skim.root'
##pseudodatatset
##filenameDat = PATH+'../test_fit/dijetFitResults_FuncType0_nParFit4_MC_1fb-1_Dinko.root'
#filenameDat = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/test_fit/toys_Bonly.root"
#
##QCD MC
#filenameBkg = PATH+'../scripts/histo_bkg_mjj.root'

inputHistNameDat = 'hist_mass_1GeV'
inputHistName = 'hist_allCutsQCD'
infBkg = TFile.Open(filenameBkg)
hBkg = infBkg.Get(inputHistName)
hBkg.Print()

infDat = TFile.Open(filenameDat)
hDat   = infDat.Get(inputHistNameDat)
#hDat.Rebin(20)

##test for significance
#filenameDat = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/scripts/blindExercise_MLfit.MaxLikelihoodFit.root'
#filenameDat = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/scripts/higgsCombineQstar5000_MLfit_mu_limit.GenerateOnly.mH120.123456.root'
#inputfileDat = TFile.Open(filenameDat)
#inputDataset = inputfileDat.Get('toys/toy_'+str(TOY))
#argset = inputDataset.get()
#x = argset.find('mjj')
#x.Print()
#print str(x.getMin())+"  "+str(x.getMax())+"  "+str(x.getMax()-x.getMin())
#x.setBinning(RooUniformBinning(x.getMin(),x.getMax(),int(x.getMax()-x.getMin())))
#dataHist_data=inputDataset.binnedClone("hist_data")
#hDat = dataHist_data.createHistogram('histDat_mass_1GeV',x) 


infSig = TFile.Open(filenameSig)
hSig = infSig.Get('h_qq_'+str(int(mass)))

#tree_sig = infSig.Get("rootTupleTree/tree")
#hSig = TH1F("hist_mass_1GeV","",13999,1,14000)
#tree_sig.Project("hist_mass_1GeV","mjj","deltaETAjj<1.3")
hSig.Print()
hSig.Draw()

#hSig   = infSig.Get(inputHistName)
#hSig.Rebin(20)


# -----------------------------------------
# define observable
#test -- restrict mjj range
#x = RooRealVar('mjj','mjj',1500,6000)
#x = RooRealVar('mjj','mjj',1118,6099)
x = RooRealVar('mjj','mjj',1118,9067)

dataHist_data=RooDataHist("RooDataHist","RooDataHist",RooArgList(x),hDat)

if fitSig: 

    # define parameters for signal fit
    m = RooRealVar('mean','mean',float(mass),float(mass)-200,float(mass)+200)
    s = RooRealVar('sigma','sigma',0.1*float(mass),0,10000)
    a = RooRealVar('alpha','alpha',1,-10,10)
    n = RooRealVar('n','n',1,0,100)
    sig = RooCBShape('sig','sig',x,m,s,a,n)        

    p  = RooRealVar('p','p',1,0,5)
    x0 = RooRealVar('x0','x0',1000,100,5000)

    bkg = RooGenericPdf('bkg','1/(exp(pow(@0/@1,@2))+1)',RooArgList(x,x0,p))

    fsig= RooRealVar('fsig','fsig',0.5,0.,1.)
    signal = RooAddPdf('signal','signal',sig,bkg,fsig)

    # -----------------------------------------
    # fit signal
    canSname = 'can_Mjj'+str(mass)
    canS = TCanvas(canSname,canSname,900,600)
    gPad.SetLogy() 

    roohistSig = RooDataHist('roohist','roohist',RooArgList(x),hSig)

    roohistSig.Print() 
    res_sig = signal.fitTo(roohistSig, RooFit.Save(ROOT.kTRUE))
    res_sig.Print()
    frame = x.frame()
    roohistSig.plotOn(frame,RooFit.Binning(166))
    signal.plotOn(frame)
    signal.plotOn(frame,RooFit.Components('bkg'),RooFit.LineColor(ROOT.kRed),RooFit.LineWidth(2),RooFit.LineStyle(ROOT.kDashed))
    #frame.GetXaxis().SetRangeUser(1118,6099)
    frame.GetXaxis().SetRangeUser(1500,6000)
    frame.GetXaxis().SetTitle('m_{jj} (GeV)')
    frame.Draw()

    parsSig = signal.getParameters(roohistSig)
    parsSig.setAttribAll('Constant', True)


if histpdfSig:
    
    # -----------------------------------------
    # hist pdf signal
    canSname = 'can_Mjj'+str(mass)
    canS = TCanvas(canSname,canSname,900,600)
    gPad.SetLogy() 

    roohistSig = RooDataHist('roohist','roohist',RooArgList(x),hSig)
    roohistSig.Print()
    signal = RooHistPdf('signal','signal',RooArgSet(x),roohistSig)
    signal.Print()
    frame = x.frame()
    roohistSig.plotOn(frame,RooFit.Binning(166))
    signal.plotOn(frame,RooFit.Binning(166),RooFit.LineColor(ROOT.kRed),RooFit.LineWidth(2),RooFit.LineStyle(ROOT.kDashed))

    #frame.GetXaxis().SetRangeUser(1118,6099)
    frame.GetXaxis().SetRangeUser(1500,6000)
    frame.GetXaxis().SetTitle('m_{jj} (GeV)')
    frame.Draw()

#    parsSig = signal.getParameters(roohistSig)
#    parsSig.setAttribAll('Constant', True)




if fitDat: 

    # -----------------------------------------
    # define parameters for background
    # function 0 (standard parametrization)
    NBINS = 166
    #if fitModel==0:
    p1 = RooRealVar('p1','p1',6.21862535247,0.,100)
    p2 = RooRealVar('p2','p2',6.48308946408,0,60)
    p3 = RooRealVar('p3','p3',0.217160577769,-10,10)

    background = RooGenericPdf('background','(pow(1-@0/13000,@1)/pow(@0/13000,@2+@3*log(@0/13000)))',RooArgList(x,p1,p2,p3))
    background_norm = RooRealVar('background_norm','background_norm',434176,0,10000000)
    
    ebkg = RooExtendPdf("ebkg","extended background p.d.f",background,background_norm)
    
    ##variation 1, with one more parameter 
    #if fitModel==1:
    #  TMath::Power(1-x/8000,[1])*(1+[4]*x/8000) ) / ( TMath::Power(x/8000,[2]+[3]*log(x/8000))
    #
    
    #giulia - test for significance
    roohistBkg = RooDataHist('roohist','roohist',RooArgList(x),hDat)
    #roohistBkg = RooDataHist('roohist','roohist',RooArgList(x),hBkg)
    roohistBkg.Print()
    #res = background.fitTo(roohistBkg, RooFit.Save(ROOT.kTRUE))
    res = ebkg.fitTo(roohistBkg, RooFit.Save(ROOT.kTRUE))
    res.Print()

    # -----------------------------------------
    # plot background
    canBname = 'can_Mjj_Data'
    canB = TCanvas(canBname,canBname,900,600)
    gPad.SetLogy() 
    canB.cd(1).SetBottomMargin(0.4)

    frame1 = x.frame()
    frame2 = x.frame()
    roohistBkg.plotOn(frame1,RooFit.Binning(NBINS))
    background.plotOn(frame1)
    hpull = frame1.pullHist()
    frame2.addPlotable(hpull,'p')

    frame1.SetMinimum(0.5)
    frame1.GetXaxis().SetTitle('')
    frame1.GetXaxis().SetLabelSize(0.0)
    frame1.GetYaxis().SetTickLength(0.06)
    frame1.Draw()

    pad = TPad('pad','pad',0.,0.,1.,1.)
    pad.SetTopMargin(0.6)
    pad.SetFillColor(0)
    pad.SetFillStyle(0)
    pad.Draw()
    pad.cd(0)
    frame2.SetMinimum(-5)
    frame2.SetMaximum(5)
    frame2.GetYaxis().SetNdivisions(505)
    frame2.GetXaxis().SetTitleOffset(0.9)
    frame2.GetYaxis().SetTitleOffset(0.8)
    frame2.GetYaxis().SetTickLength(0.06)
    frame2.GetYaxis().SetTitleSize(0.05)
    frame2.GetYaxis().SetLabelSize(0.03)
    frame2.GetYaxis().SetTitle('(Data-Fit)/Error')
    frame2.GetXaxis().SetTitle('m_{jj} (GeV)')
    frame2.Draw();

    parsBkg = background.getParameters(roohistBkg)
    if bkgConst:
      parsBkg.setAttribAll('Constant', True)

if histpdfBkg:
    
    # -----------------------------------------
    # hist pdf bkg 
    canBhist_name = 'can_Mjj'+str(mass)
    canBhist = TCanvas(canSname,canSname,900,600)
    gPad.SetLogy() 

    roohistBkg = RooDataHist('roohist','roohist',RooArgList(x),hBkg)
    roohistBkg.Print()
    background = RooHistPdf('background','background',RooArgSet(x),roohistBkg)
    background.Print()
    background_norm = RooRealVar('background_norm','background_norm',1,0,10000000)
    frame = x.frame()
    roohistBkg.plotOn(frame,RooFit.Binning(166))
    background.plotOn(frame,RooFit.Binning(166),RooFit.LineColor(ROOT.kRed),RooFit.LineWidth(2),RooFit.LineStyle(ROOT.kDashed))

    #frame.GetXaxis().SetRangeUser(1118,6099)
    #frame.GetXaxis().SetRangeUser(1500,6000)
    frame.GetXaxis().SetRangeUser(1100,4000)
    frame.GetXaxis().SetTitle('m_{jj} (GeV)')
    frame.Draw()

# -----------------------------------------
# write everything to a workspace to make a datacard
# giulia- test - create datacard and workspace with restrincted range of m
#dcFN = 'Qstar'+str(mass)+'_datacard_range1018to3500.txt'
#wsFN = 'Qstar'+str(mass)+'_workspace_range1018to3500.root'
#outdir_datacards = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/datacards_pseudodatasetDinko/"
#outdir_datacards = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/datacards_testForSignificance/"
#outdir_datacards = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/datacards_testForSignificance_fitData/"
#outdir_datacards = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/datacards_testForSignificance_fittDataBonly/"

if bkgConst:
  dcFN = outdir_datacards+'Qstar'+str(mass)+'_datacard_const_'+name+'.txt'
  wsFN = outdir_datacards+'Qstar'+str(mass)+'_workspace_const_'+name+'.root'
elif bkgNuisance:
  dcFN = outdir_datacards+'Qstar'+str(mass)+'_datacard_nuisance.txt'
  wsFN = outdir_datacards+'Qstar'+str(mass)+'_workspace_nuisance.root'
  #dcFN = outdir_datacards+'Qstar'+str(mass)+'_datacard_nuisance_testForSignificance.txt'
  #wsFN = outdir_datacards+'Qstar'+str(mass)+'_workspace_nuisance_testForSignificance.root'

else:
  dcFN = outdir_datacards+'Qstar'+str(mass)+'_datacard_'+name+'.txt'
  wsFN = outdir_datacards+'Qstar'+str(mass)+'_workspace_'+name+'.root'



nObs = dataHist_data.sumEntries();
#nObs = roohistBkg.sumEntries();
#nBkg = roohistBkg.sumEntries();


w = RooWorkspace('w','workspace')
getattr(w,'import')(signal)
getattr(w,'import')(background)
#getattr(w,'import')(roohistBkg,RooFit.Rename("data_obs"))  
getattr(w,'import')(dataHist_data,RooFit.Rename("data_obs"))  
getattr(w,'import')(background_norm)
w.Print()
w.writeToFile(wsFN)

# -----------------------------------------
# write a datacard
LUMI = options.lumi
signalCrossSection = options.sigXS
signalEfficiency = options.sigEff
ExpectedSignalRate = signalCrossSection*LUMI*signalEfficiency

datacard = open(dcFN,'w')
datacard.write('imax 1\n')
datacard.write('jmax 1\n')
datacard.write('kmax *\n')
datacard.write('---------------\n')
datacard.write('shapes * * '+wsFN+' w:$PROCESS\n')
datacard.write('---------------\n')
datacard.write('bin 1\n')    
datacard.write('observation -1\n')
datacard.write('------------------------------\n')
datacard.write('bin          1          1\n')          
datacard.write('process      signal     background\n')
datacard.write('process      0          1\n')          
datacard.write('rate         '+str(ExpectedSignalRate)+'      1\n')
datacard.write('------------------------------\n')      
#nuisance parameters --- gaussian prior
if bkgNuisance:
  datacard.write('p1  param    '+str(p1.getValV())+'   '+str(p1.getError())+'\n')
  datacard.write('p2  param    '+str(p2.getValV())+'   '+str(p2.getError())+'\n')
  datacard.write('p3  param    '+str(p3.getValV())+'   '+str(p3.getError())+'\n')
#flat parameters --- flat prior
if not (bkgNuisance or histpdfBkg): 
  datacard.write('background_norm  flatParam\n')
  datacard.write('p1  flatParam\n')
  datacard.write('p2  flatParam\n')
  datacard.write('p3  flatParam\n')


   
 
##----- keep the GUI alive ------------
#if __name__ == '__main__':
#  rep = ''
#  while not rep in ['q','Q']:
#    rep = raw_input('enter "q" to quit: ')
#    if 1 < len(rep):
#      rep = rep[0]
