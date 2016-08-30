#!/usr/bin/env python

import sys, os, copy, re
from argparse import ArgumentParser
from array import array
import ROOT
ROOT.gROOT.SetBatch(True)

#if not os.path.exists(os.path.join(os.environ["CMSSW_BASE"],"src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer_cc.so")):
#    ROOT.gROOT.ProcessLine(".L " + os.path.join(os.environ["CMSSW_BASE"],"src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer.cc")+"+")
#else:
#    ROOT.gSystem.Load(os.path.join(os.environ["CMSSW_BASE"],"src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer_cc.so"))


def main():
    # usage description
    usage = "Example: ./scripts/createDatacards.py --inputData inputs/rawhistV7_Run2015D_scoutingPFHT_UNBLINDED_649_838_JEC_HLTplusV7_Mjj_cor_smooth.root --dataHistname mjj_mjjcor_gev --inputSig inputs/ResonanceShapes_gg_13TeV_Scouting_Spring15.root -f gg -o datacards -l 1866 --lumiUnc 0.027 --massrange 1000 1500 50 --runFit --p1 5 --p2 7 --p3 0.4 --massMin 838 --massMax 2037 --fitStrategy 2"

    # input parameters
    parser = ArgumentParser(description='Script that creates combine datacards and corresponding RooFit workspaces',epilog=usage)

    parser.add_argument("--inputData", dest="inputData", required=True,
                        help="Input data spectrum",
                        metavar="INPUT_DATA")

    parser.add_argument("--dataHistname", dest="dataHistname", required=True,
                        help="Data histogram name",
                        metavar="DATA_HISTNAME")

    parser.add_argument("--inputSig", dest="inputSig", required=True,
                        help="Input signal shapes",
                        metavar="INPUT_SIGNAL")

    parser.add_argument("-f", "--final_state", dest="final_state", required=True,
                        help="Final state (e.g. qq, qg, gg)",
                        metavar="FINAL_STATE")

    parser.add_argument("-f2", "--type", dest="atype", required=True, help="Type (e.g. hG, lG, hR, lR)")

    parser.add_argument("-o", "--output_path", dest="output_path", required=True,
                        help="Output path where datacards and workspaces will be stored",
                        metavar="OUTPUT_PATH")

    parser.add_argument("-l", "--lumi", dest="lumi", required=True,
                        default=1000., type=float,
                        help="Integrated luminosity in pb-1 (default: %(default).1f)",
                        metavar="LUMI")

    parser.add_argument("--massMin", dest="massMin",
                        default=500, type=int,
                        help="Lower bound of the mass range used for fitting (default: %(default)s)",
                        metavar="MASS_MIN")

    parser.add_argument("--massMax", dest="massMax",
                        default=1200, type=int,
                        help="Upper bound of the mass range used for fitting (default: %(default)s)",
                        metavar="MASS_MAX")

    parser.add_argument("--p1", dest="p1",
                        default=5.0000e-03, type=float,
                        help="Fit function p1 parameter (default: %(default)e)",
                        metavar="P1")

    parser.add_argument("--p2", dest="p2",
                        default=9.1000e+00, type=float,
                        help="Fit function p2 parameter (default: %(default)e)",
                        metavar="P2")

    parser.add_argument("--p3", dest="p3",
                        default=5.0000e-01, type=float,
                        help="Fit function p3 parameter (default: %(default)e)",
                        metavar="P3")

    parser.add_argument("--lumiUnc", dest="lumiUnc",
                        required=True, type=float,
                        help="Relative uncertainty in the integrated luminosity",
                        metavar="LUMI_UNC")

    parser.add_argument("--jesUnc", dest="jesUnc",
                        type=float,
                        help="Relative uncertainty in the jet energy scale",
                        metavar="JES_UNC")

    parser.add_argument("--jerUnc", dest="jerUnc",
                        type=float,
                        help="Relative uncertainty in the jet energy resolution",
                        metavar="JER_UNC")

    parser.add_argument("--sqrtS", dest="sqrtS",
                        default=13000., type=float,
                        help="Collision center-of-mass energy (default: %(default).1f)",
                        metavar="SQRTS")

    parser.add_argument("--fixP3", dest="fixP3", default=False, action="store_true", help="Fix the fit function p3 parameter")

    parser.add_argument("--runFit", dest="runFit", default=False, action="store_true", help="Run the fit")

    parser.add_argument("--fitBonly", dest="fitBonly", default=False, action="store_true", help="Run B-only fit")

    parser.add_argument("--fixBkg", dest="fixBkg", default=False, action="store_true", help="Fix all background parameters")

    parser.add_argument("--decoBkg", dest="decoBkg", default=False, action="store_true", help="Decorrelate background parameters")

    parser.add_argument("--fitStrategy", dest="fitStrategy", type=int, default=1, help="Fit strategy (default: %(default).1f)")

    parser.add_argument("--debug", dest="debug", default=False, action="store_true", help="Debug printout")

    parser.add_argument("--postfix", dest="postfix", default='', help="Postfix for the output file names (default: %(default)s)")

    parser.add_argument("--pyes", dest="pyes", default=False, action="store_true", help="Make files for plots")

    parser.add_argument("--jyes", dest="jyes", default=False, action="store_true", help="Make files for JES/JER plots")

    parser.add_argument("--pdir", dest="pdir", default='testarea', help="Name a directory for the plots (default: %(default)s)")

    parser.add_argument("--chi2", dest="chi2", default=False, action="store_true", help="Compute chi squared")

    parser.add_argument("--widefit", dest="widefit", default=False, action="store_true", help="Fit with wide bin hist")

    mass_group = parser.add_mutually_exclusive_group(required=True)
    mass_group.add_argument("--mass",
                            type=int,
                            nargs = '*',
                            default = 1000,
                            help="Mass can be specified as a single value or a whitespace separated list (default: %(default)i)"
                            )
    mass_group.add_argument("--massrange",
                            type=int,
                            nargs = 3,
                            help="Define a range of masses to be produced. Format: min max step",
                            metavar = ('MIN', 'MAX', 'STEP')
                            )
    mass_group.add_argument("--masslist",
                            help = "List containing mass information"
                            )

    args = parser.parse_args()

    if args.atype == 'hG':
	fstr = "bbhGGBB"
	in2 = 'bcorrbin/binmodh.root'
    elif args.atype == 'hR':
	fstr = "bbhRS"
	in2 = 'bcorrbin/binmodh.root'
    elif args.atype == 'lG':
	fstr = "bblGGBB"
	in2 = 'bcorrbin/binmodl.root'
    else:
	fstr = "bblRS"
	in2 = 'bcorrbin/binmodl.root'

    # check if the output directory exists
    if not os.path.isdir( os.path.join(os.getcwd(),args.output_path) ):
        os.mkdir( os.path.join(os.getcwd(),args.output_path) )

    # mass points for which resonance shapes will be produced
    masses = []

    if args.massrange != None:
        MIN, MAX, STEP = args.massrange
        masses = range(MIN, MAX+STEP, STEP)
    elif args.masslist != None:
        # A mass list was provided
        print  "Will create mass list according to", args.masslist
        masslist = __import__(args.masslist.replace(".py",""))
        masses = masslist.masses
    else:
        masses = args.mass

    # sort masses
    masses.sort()

    # import ROOT stuff
    from ROOT import gStyle, TFile, TH1F, TH1D, TGraph, kTRUE, kFALSE, TCanvas, TLegend, TPad, TLine
    from ROOT import RooHist, RooRealVar, RooDataHist, RooArgList, RooArgSet, RooAddPdf, RooFit, RooGenericPdf, RooWorkspace, RooMsgService, RooHistPdf, RooExtendPdf

    if not args.debug:
        RooMsgService.instance().setSilentMode(kTRUE)
        RooMsgService.instance().setStreamStatus(0,kFALSE)
        RooMsgService.instance().setStreamStatus(1,kFALSE)

    # input data file
    inputData = TFile(args.inputData)
    # input data histogram
    hData = inputData.Get(args.dataHistname)

    inData2 = TFile(in2)
    hData2 = inData2.Get('h_data')

    # input sig file
    inputSig = TFile(args.inputSig)

    sqrtS = args.sqrtS

    # mass variable
    mjj = RooRealVar('mjj','mjj',float(args.massMin),float(args.massMax))

    # integrated luminosity and signal cross section
    lumi = args.lumi
    signalCrossSection = 1. # set to 1. so that the limit on r can be interpreted as a limit on the signal cross section

    for mass in masses:

        print ">> Creating datacard and workspace for %s resonance with m = %i GeV..."%(args.final_state, int(mass))

        # get signal shape
        hSig = inputSig.Get( "h_" + args.final_state + "_" + str(int(mass)) )
        # normalize signal shape to the expected event yield (works even if input shapes are not normalized to unity)
        hSig.Scale(signalCrossSection*lumi/hSig.Integral()) # divide by a number that provides roughly an r value of 1-10
        rooSigHist = RooDataHist('rooSigHist','rooSigHist',RooArgList(mjj),hSig)
        print 'Signal acceptance:', (rooSigHist.sumEntries()/hSig.Integral())
        signal = RooHistPdf('signal','signal',RooArgSet(mjj),rooSigHist)
        signal_norm = RooRealVar('signal_norm','signal_norm',0,-1e+05,1e+05)
	signal_norm2 = RooRealVar('signal_norm2','signal_norm2',0,-1e+05,1e+05)
	signal_norm3 = RooRealVar('signal_norm3','signal_norm3',0,-1e+05,1e+05)        
	signal_norm4 = RooRealVar('signal_norm4','signal_norm4',0,-1e+05,1e+05)
	signal_norm5 = RooRealVar('signal_norm5','signal_norm5',0,-1e+05,1e+05)

        if args.fitBonly: 
		signal_norm.setConstant()
		signal_norm2.setConstant()
		signal_norm3.setConstant()
		signal_norm4.setConstant()
		signal_norm5.setConstant()

        p1 = RooRealVar('p1','p1',args.p1,0.,100.)
        p2 = RooRealVar('p2','p2',args.p2,0.,60.)
        p3 = RooRealVar('p3','p3',args.p3,-10.,10.)
	p4 = RooRealVar('p4','p4',5.6,-50.,50.)
	p5 = RooRealVar('p5','p5',10.,-50.,50.)
	p6 = RooRealVar('p6','p6',.016,-50.,50.)
	p7 = RooRealVar('p7','p7',8.,-50.,50.)
	p8 = RooRealVar('p8','p8',.22,-50.,50.)
	p9 = RooRealVar('p9','p9',14.1,-50.,50.)
	p10 = RooRealVar('p10','p10',8.,-50.,50.)
	p11 = RooRealVar('p11','p11',4.8,-50.,50.)
	p12 = RooRealVar('p12','p12',7.,-50.,50.)
	p13 = RooRealVar('p13','p13',7.,-50.,50.)
	p14 = RooRealVar('p14','p14',7.,-50.,50.)
	p15 = RooRealVar('p15','p15',1.,-50.,50.)
	p16 = RooRealVar('p16','p16',9.,-50.,50.)
	p17 = RooRealVar('p17','p17',0.6,-50.,50.)

        if args.fixP3: p3.setConstant()

        background = RooGenericPdf('background','(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)))'%(sqrtS,sqrtS,sqrtS),RooArgList(mjj,p1,p2,p3))
        dataInt = hData.Integral(hData.GetXaxis().FindBin(float(args.massMin)),hData.GetXaxis().FindBin(float(args.massMax)))
        background_norm = RooRealVar('background_norm','background_norm',dataInt,0.,1e+08)

	background2 = RooGenericPdf('background2','(pow(@0/%.1f,-@1)*pow(1-@0/%.1f,@2))'%(sqrtS,sqrtS),RooArgList(mjj,p4,p5))
        dataInt2 = hData.Integral(hData.GetXaxis().FindBin(float(args.massMin)),hData.GetXaxis().FindBin(float(args.massMax)))
        background2_norm = RooRealVar('background2_norm','background2_norm',dataInt2,0.,1e+08)

	background3 = RooGenericPdf('background3','(1/pow(@1+@0/%.1f,@2))'%(sqrtS),RooArgList(mjj,p6,p7))
        dataInt3 = hData.Integral(hData.GetXaxis().FindBin(float(args.massMin)),hData.GetXaxis().FindBin(float(args.massMax)))
        background3_norm = RooRealVar('background3_norm','background3_norm',dataInt3,0.,1e+08)

	background4 = RooGenericPdf('background4','(1/pow(@1+@2*@0/%.1f+pow(@0/%.1f,2),@3))'%(sqrtS,sqrtS),RooArgList(mjj,p8,p9,p10))
        dataInt4 = hData.Integral(hData.GetXaxis().FindBin(float(args.massMin)),hData.GetXaxis().FindBin(float(args.massMax)))
        background4_norm = RooRealVar('background4_norm','background4_norm',dataInt4,0.,1e+08)

	background5 = RooGenericPdf('background5','(pow(@0/%.1f,-@1)*pow(1-pow(@0/%.1f,1/3),@2))'%(sqrtS,sqrtS),RooArgList(mjj,p11,p12))
        dataInt5 = hData.Integral(hData.GetXaxis().FindBin(float(args.massMin)),hData.GetXaxis().FindBin(float(args.massMax)))
        background5_norm = RooRealVar('background5_norm','background5_norm',dataInt5,0.,1e+08)

	background6 = RooGenericPdf('background6','(pow(@0/%.1f,2)+@1*@0/%.1f+@2)'%(sqrtS,sqrtS),RooArgList(mjj,p13,p14))
        dataInt6 = hData.Integral(hData.GetXaxis().FindBin(float(args.massMin)),hData.GetXaxis().FindBin(float(args.massMax)))
        background_norm6 = RooRealVar('background_norm6','background_norm6',dataInt6,0.,1e+08)

	background7 = RooGenericPdf('background7','((-1+@1*@0/%.1f)*pow(@0/%.1f,@2+@3*log(@0/%.1f)))'%(sqrtS,sqrtS,sqrtS),RooArgList(mjj,p15,p16,p17))
        dataInt7 = hData.Integral(hData.GetXaxis().FindBin(float(args.massMin)),hData.GetXaxis().FindBin(float(args.massMax)))
        background_norm7 = RooRealVar('background_norm7','background_norm7',dataInt7,0.,1e+08)

	#Extend PDFs
	
	exts = RooExtendPdf('extsignal','Extended Signal Pdf',signal,signal_norm)
	extb = RooExtendPdf('extbackground','Extended Background Pdf',background,background_norm)
	exts2 = RooExtendPdf('extsignal2','Extended Signal Pdf2',signal,signal_norm2)
        extb2 = RooExtendPdf('extbackground2','Extended Background Pdf2',background2,background2_norm)
	exts3 = RooExtendPdf('extsignal3','Extended Signal Pdf3',signal,signal_norm3)
        extb3 = RooExtendPdf('extbackground3','Extended Background Pdf3',background3,background3_norm)
	exts4 = RooExtendPdf('extsignal4','Extended Signal Pdf4',signal,signal_norm4)
        extb4 = RooExtendPdf('extbackground4','Extended Background Pdf4',background4,background4_norm)
	exts5 = RooExtendPdf('extsignal5','Extended Signal Pdf5',signal,signal_norm5)
        extb5 = RooExtendPdf('extbackground5','Extended Background Pdf5',background5,background5_norm)

        # S+B model
        model = RooAddPdf("model","s+b",RooArgList(extb,exts))
	model2 = RooAddPdf("model2","s+b2",RooArgList(extb2,exts2))
	model3 = RooAddPdf("model3","s+b3",RooArgList(extb3,exts3))
	model4 = RooAddPdf("model4","s+b4",RooArgList(extb4,exts4))
	model5 = RooAddPdf("model5","s+b5",RooArgList(extb5,exts5))

	#model6 = RooAddPdf("model6","s+b6",RooArgList(background6,signal),RooArgList(background_norm6,signal_norm))
	#model7 = RooAddPdf("model7","s+b7",RooArgList(background7,signal),RooArgList(background_norm7,signal_norm))

        rooDataHist = RooDataHist('rooDatahist','rooDathist',RooArgList(mjj),hData)


        if args.runFit:
	    mframe = mjj.frame()
	    rooDataHist.plotOn(mframe, ROOT.RooFit.Name("setonedata"))
	    res = model.fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Extended(kTRUE), RooFit.Strategy(args.fitStrategy))
	    model.plotOn(mframe, ROOT.RooFit.Name("model1"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(1), ROOT.RooFit.LineColor(ROOT.EColor.kRed)) 
	    res2 = model2.fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Extended(kTRUE), RooFit.Strategy(args.fitStrategy))
#            model2.plotOn(mframe, ROOT.RooFit.Name("model2"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(1), ROOT.RooFit.LineColor(ROOT.EColor.kOrange))
	    res3 = model3.fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Extended(kTRUE), RooFit.Strategy(args.fitStrategy))
#            model3.plotOn(mframe, ROOT.RooFit.Name("model3"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(1), ROOT.RooFit.LineColor(ROOT.EColor.kGreen))
	    res4 = model4.fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Extended(kTRUE), RooFit.Strategy(args.fitStrategy))
#            model4.plotOn(mframe, ROOT.RooFit.Name("model4"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(1), ROOT.RooFit.LineColor(ROOT.EColor.kBlue))
	    res5 = model5.fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Extended(kTRUE), RooFit.Strategy(args.fitStrategy))
#            model5.plotOn(mframe, ROOT.RooFit.Name("model5"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(1), ROOT.RooFit.LineColor(ROOT.EColor.kViolet))
#	    res6 = model6.fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Strategy(args.fitStrategy))
#           model6.plotOn(mframe, ROOT.RooFit.Name("model6"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(1), ROOT.RooFit.LineColor(ROOT.EColor.kPink))
#	    res7 = model7.fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Strategy(args.fitStrategy))
#           model7.plotOn(mframe, ROOT.RooFit.Name("model7"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(1), ROOT.RooFit.LineColor(ROOT.EColor.kAzure))

	    rooDataHist2 = RooDataHist('rooDatahist2','rooDathist2',RooArgList(mjj),hData2)
#	    rooDataHist2.plotOn(mframe, ROOT.RooFit.Name("data"))

	    if args.pyes:
	    	c = TCanvas("c","c",800,800)
		mframe.SetAxisRange(300.,1300.)
	    	c.SetLogy()
#	    	mframe.SetMaximum(10)
#	    	mframe.SetMinimum(1)
	    	mframe.Draw()
	    	fitname = args.pdir+'/5funcfit_m'+str(mass)+fstr+'.pdf'
	    	c.SaveAs(fitname)

#	        cpull = TCanvas("cpull","cpull",800,800)
#	    	pulls = mframe.pullHist("data","model1")
#	    	pulls.Draw("ABX")
#	   	pullname = args.pdir+'/pull_m'+str(mass)+fstr+'.pdf'
#	    	cpull.SaveAs(pullname)

#		cpull2 = TCanvas("cpull2","cpull2",800,800)
 #               pulls2 = mframe.pullHist("setonedata","model1")
  #              pulls2.Draw("ABX")
   #             pull2name = args.pdir+'/pull2_m'+str(mass)+fstr+'.pdf'
    #            cpull2.SaveAs(pull2name)

	    if args.widefit:	
		mframew = mjj.frame()
    	        rooDataHist2.plotOn(mframew, ROOT.RooFit.Name("data"))
                res6 = model.fitTo(rooDataHist2, RooFit.Save(kTRUE), RooFit.Strategy(args.fitStrategy))
            	model.plotOn(mframew, ROOT.RooFit.Name("model1"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(1), ROOT.RooFit.LineColor(ROOT.EColor.kRed))
            	res7 = model2.fitTo(rooDataHist2, RooFit.Save(kTRUE), RooFit.Strategy(args.fitStrategy))
            	model2.plotOn(mframew, ROOT.RooFit.Name("model2"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(1), ROOT.RooFit.LineColor(ROOT.EColor.kOrange))
            	res8 = model3.fitTo(rooDataHist2, RooFit.Save(kTRUE), RooFit.Strategy(args.fitStrategy))
            	model3.plotOn(mframew, ROOT.RooFit.Name("model3"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(1), ROOT.RooFit.LineColor(ROOT.EColor.kGreen))
            	res9 = model4.fitTo(rooDataHist2, RooFit.Save(kTRUE), RooFit.Strategy(args.fitStrategy))
            	model4.plotOn(mframew, ROOT.RooFit.Name("model4"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(1), ROOT.RooFit.LineColor(ROOT.EColor.kBlue))
            	res10 = model5.fitTo(rooDataHist2, RooFit.Save(kTRUE), RooFit.Strategy(args.fitStrategy))
            	model5.plotOn(mframew, ROOT.RooFit.Name("model5"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(1), ROOT.RooFit.LineColor(ROOT.EColor.kViolet))

                if args.pyes:
                    c = TCanvas("c","c",800,800)
                    mframew.SetAxisRange(300.,1300.)
                    c.SetLogy()
#                   mframew.SetMaximum(10)
#                   mframew.SetMinimum(1)
                    mframew.Draw()
                    fitname = args.pdir+'/5funcfittowide_m'+str(mass)+fstr+'.pdf'
		    c.SaveAs(fitname)

                    cpull = TCanvas("cpull","cpull",800,800)
                    pulls = mframew.pullHist("data","model1")
                    pulls.Draw("ABX")
                    pullname = args.pdir+'/pullwidefit_m'+str(mass)+fstr+'.pdf'
                    cpull.SaveAs(pullname)


	    if args.chi2:
		    fullInt = model.createIntegral(RooArgSet(mjj))
		    norm = dataInt/fullInt.getVal()
		    chi1 = 0.
		    fullInt2 = model2.createIntegral(RooArgSet(mjj))
        	    norm2 = dataInt2/fullInt2.getVal()
	      	    chi2 = 0.
		    fullInt3 = model3.createIntegral(RooArgSet(mjj))
       		    norm3 = dataInt3/fullInt3.getVal()
	            chi3 = 0.
		    fullInt4 = model4.createIntegral(RooArgSet(mjj))
       		    norm4 = dataInt4/fullInt4.getVal()
         	    chi4 = 0.
		    fullInt5 = model5.createIntegral(RooArgSet(mjj))
	            norm5 = dataInt5/fullInt5.getVal()
     	            chi5 = 0.
		    for i in range(args.massMin, args.massMax):
        	        new = 0
			new2 = 0
			new3 = 0
			new4 = 0
			new5 = 0
			height = hData.GetBinContent(i)
	        	xLow = hData.GetXaxis().GetBinLowEdge(i)
			xUp = hData.GetXaxis().GetBinLowEdge(i+1)
			obs = height*(xUp-xLow)
			mjj.setRange("intrange",xLow,xUp)
			integ = model.createIntegral(RooArgSet(mjj),ROOT.RooFit.NormSet(RooArgSet(mjj)),ROOT.RooFit.Range("intrange"))
			exp = integ.getVal()*norm
			new = pow(exp-obs,2)/exp
                	chi1 = chi1 + new
			integ2 = model2.createIntegral(RooArgSet(mjj),ROOT.RooFit.NormSet(RooArgSet(mjj)),ROOT.RooFit.Range("intrange"))
                	exp2 = integ2.getVal()*norm2
                	new2 = pow(exp2-obs,2)/exp2
                	chi2 = chi2 + new2
			integ3 = model3.createIntegral(RooArgSet(mjj),ROOT.RooFit.NormSet(RooArgSet(mjj)),ROOT.RooFit.Range("intrange"))
                	exp3 = integ3.getVal()*norm3
                	new3 = pow(exp3-obs,2)/exp3
                	chi3 = chi3 + new3
			integ4 = model4.createIntegral(RooArgSet(mjj),ROOT.RooFit.NormSet(RooArgSet(mjj)),ROOT.RooFit.Range("intrange"))
                	exp4 = integ4.getVal()*norm4
			if exp4 != 0:
                	    new4 = pow(exp4-obs,2)/exp4
                	else:
			    new4 = 0
			chi4 = chi4 + new4
			integ5 = model5.createIntegral(RooArgSet(mjj),ROOT.RooFit.NormSet(RooArgSet(mjj)),ROOT.RooFit.Range("intrange"))
                	exp5 = integ5.getVal()*norm5
                	new5 = pow(exp5-obs,2)/exp5
                	chi5 = chi5 + new5
	    	    print "chi1 %d "%(chi1)
	    	    print "chi2 %d "%(chi2)
	    	    print "chi3 %d "%(chi3)
	    	    print "chi4 %d "%(chi4)
	    	    print "chi5 %d "%(chi5)

	    if not args.decoBkg: 
		print " "
		res.Print()
	        res2.Print()
		res3.Print()
		res4.Print()
		res5.Print()
#		res6.Print()
#		res7.Print()

            # decorrelated background parameters for Bayesian limits
            if args.decoBkg:
                signal_norm.setConstant()
                res = model.fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Strategy(args.fitStrategy))
                res.Print()
                ## temp workspace for the PDF diagonalizer
                w_tmp = RooWorkspace("w_tmp")
                deco = PdfDiagonalizer("deco",w_tmp,res)
                # here diagonalizing only the shape parameters since the overall normalization is already decorrelated
                background_deco = deco.diagonalize(background)
                print "##################### workspace for decorrelation"
                w_tmp.Print("v")
                print "##################### original parameters"
                background.getParameters(rooDataHist).Print("v")
                print "##################### decorrelated parameters"
                # needed if want to evaluate limits without background systematics
                if args.fixBkg:
                    w_tmp.var("deco_eig1").setConstant()
                    w_tmp.var("deco_eig2").setConstant()
                    if not args.fixP3: w_tmp.var("deco_eig3").setConstant()
                background_deco.getParameters(rooDataHist).Print("v")
                print "##################### original pdf"
                background.Print()
                print "##################### decorrelated pdf"
                background_deco.Print()
                # release signal normalization
                signal_norm.setConstant(kFALSE)
                # set the background normalization range to +/- 5 sigma
                bkg_val = background_norm.getVal()
                bkg_error = background_norm.getError()
                background_norm.setMin(bkg_val-5*bkg_error)
                background_norm.setMax(bkg_val+5*bkg_error)
                background_norm.Print()
                # change background PDF names
                background.SetName("background_old")
                background_deco.SetName("background")

        # needed if want to evaluate limits without background systematics
        if args.fixBkg:
            background_norm.setConstant()
            p1.setConstant()
            p2.setConstant()
            p3.setConstant()

        # -----------------------------------------
        # dictionaries holding systematic variations of the signal shape
        hSig_Syst = {}
        hSig_Syst_DataHist = {}
        sigCDF = TGraph(hSig.GetNbinsX()+1)

        # JES and JER uncertainties
        if args.jesUnc != None or args.jerUnc != None:

            sigCDF.SetPoint(0,0.,0.)
            integral = 0.
            for i in range(1, hSig.GetNbinsX()+1):
                x = hSig.GetXaxis().GetBinLowEdge(i+1)
                integral = integral + hSig.GetBinContent(i)
                sigCDF.SetPoint(i,x,integral)

        if args.jesUnc != None:
            hSig_Syst['JESUp'] = copy.deepcopy(hSig)
            hSig_Syst['JESDown'] = copy.deepcopy(hSig)

        if args.jerUnc != None:
            hSig_Syst['JERUp'] = copy.deepcopy(hSig)
            hSig_Syst['JERDown'] = copy.deepcopy(hSig)

        # reset signal histograms
        for key in hSig_Syst.keys():
            hSig_Syst[key].Reset()
            hSig_Syst[key].SetName(hSig_Syst[key].GetName() + '_' + key)

        # produce JES signal shapes
        if args.jesUnc != None:
            for i in range(1, hSig.GetNbinsX()+1):
                xLow = hSig.GetXaxis().GetBinLowEdge(i)
                xUp = hSig.GetXaxis().GetBinLowEdge(i+1)
                jes = 1. - args.jesUnc
                xLowPrime = jes*xLow
                xUpPrime = jes*xUp
                hSig_Syst['JESUp'].SetBinContent(i, sigCDF.Eval(xUpPrime) - sigCDF.Eval(xLowPrime))
                jes = 1. + args.jesUnc
                xLowPrime = jes*xLow
                xUpPrime = jes*xUp
                hSig_Syst['JESDown'].SetBinContent(i, sigCDF.Eval(xUpPrime) - sigCDF.Eval(xLowPrime))
            hSig_Syst_DataHist['JESUp'] = RooDataHist('hSig_JESUp','hSig_JESUp',RooArgList(mjj),hSig_Syst['JESUp'])
            hSig_Syst_DataHist['JESDown'] = RooDataHist('hSig_JESDown','hSig_JESDown',RooArgList(mjj),hSig_Syst['JESDown'])
	    
	    if args.jyes:
		c2 = TCanvas("c2","c2",800,800)
	    	mframe2 = mjj.frame(ROOT.RooFit.Title("JES One Sigma Shifts"))
	    	mframe2.SetAxisRange(args.massMin,args.massMax)
		hSig_Syst_DataHist['JESUp'].plotOn(mframe2, ROOT.RooFit.Name("JESUP"),ROOT.RooFit.DrawOption("L"), ROOT.RooFit.DataError(2), ROOT.RooFit.LineStyle(1), ROOT.RooFit.MarkerColor(ROOT.EColor.kRed), ROOT.RooFit.LineColor(ROOT.EColor.kRed))
	    	hSig_Syst_DataHist['JESDown'].plotOn(mframe2,ROOT.RooFit.Name("JESDOWN"),ROOT.RooFit.DrawOption("L"), ROOT.RooFit.DataError(2), ROOT.RooFit.LineStyle(1), ROOT.RooFit.MarkerColor(ROOT.EColor.kBlue), ROOT.RooFit.LineColor(ROOT.EColor.kBlue))
	    	rooSigHist.plotOn(mframe2, ROOT.RooFit.DataError(2),ROOT.RooFit.Name("SIG"),ROOT.RooFit.DrawOption("L"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.MarkerColor(ROOT.EColor.kGreen), ROOT.RooFit.LineColor(ROOT.EColor.kGreen))
	    	mframe2.Draw()
		mframe2.GetXaxis().SetTitle("Dijet Mass (GeV)")
		leg = TLegend(0.7,0.8,0.9,0.9)
		leg.SetFillColor(0)
		leg.AddEntry(mframe2.findObject("SIG"),"Signal Model","l")
		leg.AddEntry(mframe2.findObject("JESUP"),"+1 Sigma","l")
		leg.AddEntry(mframe2.findObject("JESDOWN"),"-1 Sigma","l")
		leg.Draw()
	    	jesname = args.pdir+'/jes_m'+str(mass)+fstr+'.pdf'
	    	c2.SaveAs(jesname)

        # produce JER signal shapes
        if args.jesUnc != None:
            for i in range(1, hSig.GetNbinsX()+1):
                xLow = hSig.GetXaxis().GetBinLowEdge(i)
                xUp = hSig.GetXaxis().GetBinLowEdge(i+1)
                jer = 1. - args.jerUnc
                xLowPrime = jer*(xLow-float(mass))+float(mass)
                xUpPrime = jer*(xUp-float(mass))+float(mass)
                hSig_Syst['JERUp'].SetBinContent(i, sigCDF.Eval(xUpPrime) - sigCDF.Eval(xLowPrime))
                jer = 1. + args.jerUnc
                xLowPrime = jer*(xLow-float(mass))+float(mass)
                xUpPrime = jer*(xUp-float(mass))+float(mass)
                hSig_Syst['JERDown'].SetBinContent(i, sigCDF.Eval(xUpPrime) - sigCDF.Eval(xLowPrime))
            hSig_Syst_DataHist['JERUp'] = RooDataHist('hSig_JERUp','hSig_JERUp',RooArgList(mjj),hSig_Syst['JERUp'])
            hSig_Syst_DataHist['JERDown'] = RooDataHist('hSig_JERDown','hSig_JERDown',RooArgList(mjj),hSig_Syst['JERDown'])

	    if args.jyes:
	    	c3 = TCanvas("c3","c3",800,800)
            	mframe3 = mjj.frame(ROOT.RooFit.Title("JER One Sigma Shifts"))
	    	mframe3.SetAxisRange(args.massMin,args.massMax)
		hSig_Syst_DataHist['JERUp'].plotOn(mframe3,ROOT.RooFit.Name("JERUP"),ROOT.RooFit.DrawOption("L"), ROOT.RooFit.DataError(2), ROOT.RooFit.LineStyle(1), ROOT.RooFit.MarkerColor(ROOT.EColor.kRed), ROOT.RooFit.LineColor(ROOT.EColor.kRed))
            	hSig_Syst_DataHist['JERDown'].plotOn(mframe3,ROOT.RooFit.Name("JERDOWN"),ROOT.RooFit.DrawOption("L"), ROOT.RooFit.DataError(2), ROOT.RooFit.LineStyle(1), ROOT.RooFit.MarkerColor(ROOT.EColor.kBlue), ROOT.RooFit.LineColor(ROOT.EColor.kBlue))
            	rooSigHist.plotOn(mframe3,ROOT.RooFit.DrawOption("L"),ROOT.RooFit.Name("SIG"), ROOT.RooFit.DataError(2), ROOT.RooFit.LineStyle(1), ROOT.RooFit.MarkerColor(ROOT.EColor.kGreen), ROOT.RooFit.LineColor(ROOT.EColor.kGreen))
            	mframe3.Draw()
	    	mframe3.GetXaxis().SetTitle("Dijet Mass (GeV)")
		leg = TLegend(0.7,0.8,0.9,0.9)
		leg.SetFillColor(0)
                leg.AddEntry(mframe3.findObject("SIG"),"Signal Model","l")
                leg.AddEntry(mframe3.findObject("JERUP"),"+1 Sigma","l")
                leg.AddEntry(mframe3.findObject("JERDOWN"),"-1 Sigma","l")
                leg.Draw()	
		jername = args.pdir+'/jer_m'+str(mass)+fstr+'.pdf'
           	c3.SaveAs(jername)


        # -----------------------------------------
        # create a datacard and corresponding workspace
        postfix = (('_' + args.postfix) if args.postfix != '' else '')
        dcName = 'datacard_' + args.final_state + '_m' + str(mass) + postfix + '.txt'
        wsName = 'workspace_' + args.final_state + '_m' + str(mass) + postfix + '.root'

        w = RooWorkspace('w','workspace')
        getattr(w,'import')(rooSigHist,RooFit.Rename("signal"))
        if args.jesUnc != None:
            getattr(w,'import')(hSig_Syst_DataHist['JESUp'],RooFit.Rename("signal__JESUp"))
            getattr(w,'import')(hSig_Syst_DataHist['JESDown'],RooFit.Rename("signal__JESDown"))
        if args.jerUnc != None:
            getattr(w,'import')(hSig_Syst_DataHist['JERUp'],RooFit.Rename("signal__JERUp"))
            getattr(w,'import')(hSig_Syst_DataHist['JERDown'],RooFit.Rename("signal__JERDown"))
        if args.decoBkg:
            getattr(w,'import')(background_deco,ROOT.RooCmdArg())
        else:
            getattr(w,'import')(background,ROOT.RooCmdArg(),RooFit.Rename("background"))
	    getattr(w,'import')(background2,ROOT.RooCmdArg(),RooFit.Rename("background2"))
	    getattr(w,'import')(background3,ROOT.RooCmdArg(),RooFit.Rename("background3"))
	    getattr(w,'import')(background4,ROOT.RooCmdArg(),RooFit.Rename("background4"))
	    getattr(w,'import')(background5,ROOT.RooCmdArg(),RooFit.Rename("background5"))
	    getattr(w,'import')(background_norm,ROOT.RooCmdArg(),RooFit.Rename("background_norm"))
            getattr(w,'import')(background2_norm,ROOT.RooCmdArg(),RooFit.Rename("background2_norm"))
            getattr(w,'import')(background3_norm,ROOT.RooCmdArg(),RooFit.Rename("background3_norm"))
            getattr(w,'import')(background4_norm,ROOT.RooCmdArg(),RooFit.Rename("background4_norm"))
            getattr(w,'import')(background5_norm,ROOT.RooCmdArg(),RooFit.Rename("background5_norm"))

	getattr(w,'import')(res)
	getattr(w,'import')(res2)
	getattr(w,'import')(res3)
	getattr(w,'import')(res4)
	getattr(w,'import')(res5)
	getattr(w,'import')(background_norm,ROOT.RooCmdArg())
        getattr(w,'import')(signal_norm,ROOT.RooCmdArg())
	getattr(w,'import')(rooDataHist,RooFit.Rename("data_obs"))
        w.Print()
        w.writeToFile(os.path.join(args.output_path,wsName))

	beffUnc = 0.3
	boffUnc = 0.06

        datacard = open(os.path.join(args.output_path,dcName),'w')
        datacard.write('imax 1\n')
        datacard.write('jmax 1\n')
        datacard.write('kmax *\n')
        datacard.write('---------------\n')
        if args.jesUnc != None or args.jerUnc != None:
            datacard.write('shapes * * '+wsName+' w:$PROCESS w:$PROCESS__$SYSTEMATIC\n')
        else:
            datacard.write('shapes * * '+wsName+' w:$PROCESS\n')
        datacard.write('---------------\n')
        datacard.write('bin 1\n')
        datacard.write('observation -1\n')
        datacard.write('------------------------------\n')
        datacard.write('bin          1          1\n')
        datacard.write('process      signal     background\n')
        datacard.write('process      0          1\n')
        datacard.write('rate         -1         1\n')
        datacard.write('------------------------------\n')
        datacard.write('lumi  lnN    %f         -\n'%(1.+args.lumiUnc))
	datacard.write('beff  lnN    %f         -\n'%(1.+beffUnc))
	datacard.write('boff  lnN    %f         -\n'%(1.+boffUnc))
	datacard.write('bkg   lnN     -         1.03\n')
        if args.jesUnc != None:
            datacard.write('JES  shape   1          -\n')
        if args.jerUnc != None:
            datacard.write('JER  shape   1          -\n')
        # flat parameters --- flat prior
        datacard.write('background_norm  flatParam\n')
        if args.decoBkg:
            datacard.write('deco_eig1  flatParam\n')
            datacard.write('deco_eig2  flatParam\n')
            if not args.fixP3: datacard.write('deco_eig3  flatParam\n')
        else:
            datacard.write('p1  flatParam\n')
            datacard.write('p2  flatParam\n')
            if not args.fixP3: datacard.write('p3  flatParam\n')
        datacard.close()


    print '>> Datacards and workspaces created and stored in %s/'%( os.path.join(os.getcwd(),args.output_path) )


if __name__ == '__main__':
    main()

