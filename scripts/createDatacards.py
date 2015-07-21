#!/usr/bin/env python

import sys, os
from argparse import ArgumentParser


def main():
    # usage description
    usage = "Example: ./scripts/createDatacards.py --inputData inputs/dijetFitResults_FuncType0_nParFit4_MC_1invfb.root --dataHistname hist_mass_1GeV --inputSig inputs/ResonanceShapes_qg_13TeV_PU30_Spring15.root -f qg -o datacards -l 1000 --massrange 1200 7000 100"

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

    parser.add_argument("-o", "--output_path", dest="output_path", required=True,
                        help="Output path where datacards and workspaces will be stored",
                        metavar="OUTPUT_PATH")

    parser.add_argument("-l", "--lumi", dest="lumi", required=True,
                        default=1000., type=float,
                        help="Integrated luminosity in pb-1 (default: %(default).1f)",
                        metavar="LUMI")

    parser.add_argument("--massMin", dest="massMin",
                        default=1118, type=int,
                        help="Lower bound of the mass range used for fitting (default: %(default)s)",
                        metavar="MASS_MIN")

    parser.add_argument("--massMax", dest="massMax",
                        default=9067, type=int,
                        help="Upper bound of the mass range used for fitting (default: %(default)s)",
                        metavar="MASS_MAX")

    parser.add_argument("--p1", dest="p1",
                        default=7.7666e+00, type=float,
                        help="Fit function p1 parameter (default: %(default)e)",
                        metavar="P1")

    parser.add_argument("--p2", dest="p2",
                        default=5.3748e+00, type=float,
                        help="Fit function p2 parameter (default: %(default)e)",
                        metavar="P2")

    parser.add_argument("--p3", dest="p3",
                        default=5.6385e-03, type=float,
                        help="Fit function p3 parameter (default: %(default)e)",
                        metavar="P3")

    parser.add_argument("--fixP3", dest="fixP3", default=False, action="store_true", help="Fix the fit function p3 parameter")

    parser.add_argument("--runFit", dest="runFit", default=False, action="store_true", help="Run the fit")

    parser.add_argument("--fitBonly", dest="fitBonly", default=False, action="store_true", help="Run B-only fit")

    parser.add_argument("--fitStrategy", dest="fitStrategy", type=int, default=0, help="Fit strategy (default: %(default).1f)")

    parser.add_argument("--sqrtS", dest="sqrtS", type=float, default=13000., help="Collision center-of-mass energy (default: %(default).1f)")

    parser.add_argument("--theta", dest="theta", default=False, action="store_true", help="Produce histograms for the theta limit setting framework")

    parser.add_argument("--debug", dest="debug", default=False, action="store_true", help="Debug printout")

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
    from ROOT import TFile, TH1F, TH1D, kTRUE, kFALSE
    from ROOT import RooRealVar, RooDataHist, RooArgList, RooArgSet, RooAddPdf, RooFit, RooGenericPdf, RooWorkspace, RooMsgService, RooHistPdf

    if not args.debug:
        RooMsgService.instance().setSilentMode(kTRUE)
        RooMsgService.instance().setStreamStatus(0,kFALSE)
        RooMsgService.instance().setStreamStatus(1,kFALSE)

    # input data file
    inputData = TFile(args.inputData)
    # input data histogram
    hData = inputData.Get(args.dataHistname)

    # input sig file
    inputSig = TFile(args.inputSig)

    sqrtS = args.sqrtS

    for mass in masses:

        print ">> Creating datacard and workspace for %s resonance with m = %i GeV..."%(args.final_state, int(mass))

        hSig = inputSig.Get( "h_" + args.final_state + "_" + str(int(mass)) )

        # calculate acceptance of the dijet mass cut
        sigAcc = hSig.Integral(hSig.GetXaxis().FindBin(float(args.massMin)+0.5),hSig.GetXaxis().FindBin(float(args.massMax)-0.5))/hSig.Integral(1,hSig.GetXaxis().FindBin(float(args.massMax)-0.5)) # assuming 1-GeV bins

        mjj = RooRealVar('mjj','mjj',float(args.massMin),float(args.massMax))

        rooSigHist = RooDataHist('rooSigHist','rooSigHist',RooArgList(mjj),hSig)
        rooSigHist.Print()
        signal = RooHistPdf('signal','signal',RooArgSet(mjj),rooSigHist)
        signal.Print()
        signal_norm = RooRealVar('signal_norm','signal_norm',0,-1e+04,1e+04)
        if args.fitBonly: signal_norm.setConstant()
        signal_norm.Print()

        p1 = RooRealVar('p1','p1',args.p1,-100.,100.)
        p2 = RooRealVar('p2','p2',args.p2,0.,60.)
        p3 = RooRealVar('p3','p3',args.p3,-10.,10.)
        if args.fixP3: p3.setConstant()

        background = RooGenericPdf('background','(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)))'%(sqrtS,sqrtS,sqrtS),RooArgList(mjj,p1,p2,p3))
        background.Print()
        dataInt = hData.Integral(hData.GetXaxis().FindBin(float(args.massMin)),hData.GetXaxis().FindBin(float(args.massMax)))
        background_norm = RooRealVar('background_norm','background_norm',dataInt,0.,1e+07)
        background_norm.Print()

        # S+B model
        model = RooAddPdf("model","s+b",RooArgList(background,signal),RooArgList(background_norm,signal_norm))

        rooDataHist = RooDataHist('rooDatahist','rooDathist',RooArgList(mjj),hData)
        rooDataHist.Print()

        if args.runFit:
            res = model.fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Strategy(args.fitStrategy))
            res.Print()

        dcName = 'datacard_' + args.final_state + '_m' + str(mass) + '.txt'
        wsName = 'workspace_' + args.final_state + '_m' + str(mass) + '.root'

        w = RooWorkspace('w','workspace')
        getattr(w,'import')(signal)
        getattr(w,'import')(background)
        getattr(w,'import')(background_norm)
        getattr(w,'import')(rooDataHist,RooFit.Rename("data_obs"))
        w.Print()
        w.writeToFile(os.path.join(args.output_path,wsName))

        # -----------------------------------------
        # write a datacard
        lumi = args.lumi
        signalCrossSection = 1. # set to 1. so that the limit on r can be interpreted as a limit on the signal cross section
        expectedSignalRate = signalCrossSection*lumi*sigAcc

        datacard = open(os.path.join(args.output_path,dcName),'w')
        datacard.write('imax 1\n')
        datacard.write('jmax 1\n')
        datacard.write('kmax *\n')
        datacard.write('---------------\n')
        datacard.write('shapes * * '+wsName+' w:$PROCESS\n')
        datacard.write('---------------\n')
        datacard.write('bin 1\n')
        datacard.write('observation -1\n')
        datacard.write('------------------------------\n')
        datacard.write('bin          1          1\n')
        datacard.write('process      signal     background\n')
        datacard.write('process      0          1\n')
        datacard.write('rate         '+str(expectedSignalRate)+'      1\n')
        datacard.write('------------------------------\n')
        #flat parameters --- flat prior
        datacard.write('background_norm  flatParam\n')
        datacard.write('p1  flatParam\n')
        datacard.write('p2  flatParam\n')
        if not args.fixP3: datacard.write('p3  flatParam\n')
        datacard.close()

        if args.theta:
            thetaName = 'theta_' + args.final_state + '_m' + str(mass) + '.root'

            thetaFile = TFile(os.path.join(args.output_path,thetaName), 'RECREATE')

            thetaData = rooDataHist.createHistogram('dijet',mjj,RooFit.Binning(args.massMax-args.massMin,float(args.massMin),float(args.massMax)))
            thetaData.SetName('dijet__DATA')
            thetaData.Write()

            thetaSignal = rooSigHist.createHistogram('dijet',mjj,RooFit.Binning(args.massMax-args.massMin,float(args.massMin),float(args.massMax)))
            thetaSignal.Scale(expectedSignalRate)
            thetaSignal.SetName('dijet__signal')
            thetaSignal.Write()

            thetaBkg = background.createHistogram('dijet',mjj,RooFit.Binning(args.massMax-args.massMin,float(args.massMin),float(args.massMax)))
            thetaBkg.Scale(background_norm.getVal())
            thetaBkg.SetName('dijet__background')
            thetaBkg.Write()

            thetaFile.Close()

    print '>> Datacards and workspaces created and stored in %s/.'%( os.path.join(os.getcwd(),args.output_path) )


if __name__ == '__main__':
    main()

