#!/usr/bin/env python

import sys, os, copy, re
from argparse import ArgumentParser

theta_template = """files = ['theta_DUMMY_JOB.root']

model = build_model_from_rootfile(files)
model.set_signal_processes('signal')
model.add_lognormal_uncertainty('overall_signal', math.log(1.1), 'signal')
model_summary(model)

res = mle(model, 'data' , 1)
print res

zvalue = zvalue_approx(model, 'data' , 1)
print zvalue

observed_limits = bayesian_quantiles(model, 'data', 10)
print observed_limits

plot_observed = limit_band_plot(observed_limits, False)
plot_observed.write_txt('theta_DUMMY_JOB_observed_limit.txt')
"""


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
                        default=1.0003e+01, type=float,
                        help="Fit function p1 parameter (default: %(default)e)",
                        metavar="P1")

    parser.add_argument("--p2", dest="p2",
                        default=5.2648e+00, type=float,
                        help="Fit function p2 parameter (default: %(default)e)",
                        metavar="P2")

    parser.add_argument("--p3", dest="p3",
                        default=0.0000e+00, type=float,
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

    parser.add_argument("--fitStrategy", dest="fitStrategy", type=int, default=0, help="Fit strategy (default: %(default).1f)")

    parser.add_argument("--theta", dest="theta", default=False, action="store_true", help="Produce histograms for the theta limit setting framework")

    parser.add_argument("--thetaNoSyst", dest="thetaNoSyst", default=False, action="store_true", help="Also produce files for theta without systematic uncertainties")

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
    from ROOT import TFile, TH1F, TH1D, TGraph, kTRUE, kFALSE
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

    # mass variable
    mjj = RooRealVar('mjj','mjj',float(args.massMin),float(args.massMax))

    # integrated luminosity and signal cross section
    lumi = args.lumi
    signalCrossSection = 1. # set to 1. so that the limit on r can be interpreted as a limit on the signal cross section

    for mass in masses:

        print ">> Creating datacard and workspace for %s resonance with m = %i GeV..."%(args.final_state, int(mass))

        # get signal shape
        hSig = inputSig.Get( "h_" + args.final_state + "_" + str(int(mass)) )
        # normalize signal shape to the expected event yield (input shapes already normalized to unity)
        hSig.Scale(signalCrossSection*lumi)

        rooSigHist = RooDataHist('rooSigHist','rooSigHist',RooArgList(mjj),hSig)
        rooSigHist.Print()
        signal = RooHistPdf('signal','signal',RooArgSet(mjj),rooSigHist)
        signal.Print()
        signal_norm = RooRealVar('signal_norm','signal_norm',0,-1e+04,1e+04)
        if args.fitBonly: signal_norm.setConstant()
        signal_norm.Print()

        p1 = RooRealVar('p1','p1',args.p1,0.,100.)
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

        # -----------------------------------------
        # create a datacard and corresponding workspace
        dcName = 'datacard_' + args.final_state + '_m' + str(mass) + '.txt'
        wsName = 'workspace_' + args.final_state + '_m' + str(mass) + '.root'

        w = RooWorkspace('w','workspace')
        getattr(w,'import')(rooSigHist,RooFit.Rename("signal"))
        if args.jesUnc != None:
            getattr(w,'import')(hSig_Syst_DataHist['JESUp'],RooFit.Rename("signal__JESUp"))
            getattr(w,'import')(hSig_Syst_DataHist['JESDown'],RooFit.Rename("signal__JESDown"))
        if args.jerUnc != None:
            getattr(w,'import')(hSig_Syst_DataHist['JERUp'],RooFit.Rename("signal__JERUp"))
            getattr(w,'import')(hSig_Syst_DataHist['JERDown'],RooFit.Rename("signal__JERDown"))
        getattr(w,'import')(background)
        getattr(w,'import')(background_norm)
        getattr(w,'import')(rooDataHist,RooFit.Rename("data_obs"))
        w.Print()
        w.writeToFile(os.path.join(args.output_path,wsName))

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
        if args.jesUnc != None:
            datacard.write('JES  shape   1          -\n')
        if args.jerUnc != None:
            datacard.write('JER  shape   1          -\n')
        # flat parameters --- flat prior
        datacard.write('background_norm  flatParam\n')
        datacard.write('p1  flatParam\n')
        datacard.write('p2  flatParam\n')
        if not args.fixP3: datacard.write('p3  flatParam\n')
        datacard.close()

        # -----------------------------------------
        # create input histograms for the theta limit setting framework
        if args.theta:
            thetaName = 'theta_' + args.final_state + '_m' + str(mass) + '.root'

            thetaFile = TFile(os.path.join(args.output_path,thetaName), 'RECREATE')
            thetaFile.cd()

            thetaData = rooDataHist.createHistogram('dijet',mjj,RooFit.Binning(args.massMax-args.massMin,float(args.massMin),float(args.massMax)))
            thetaData.SetName('dijet__DATA')
            thetaData.Write()

            thetaBkg = background.createHistogram('dijet',mjj,RooFit.Binning(args.massMax-args.massMin,float(args.massMin),float(args.massMax)))
            thetaBkg.Scale(background_norm.getVal())
            thetaBkg.SetName('dijet__background')
            thetaBkg.Write()

            thetaSignal = rooSigHist.createHistogram('dijet',mjj,RooFit.Binning(args.massMax-args.massMin,float(args.massMin),float(args.massMax)))
            thetaSignal.SetName('dijet__signal')
            thetaSignal.Write()

            if args.thetaNoSyst:
                thetaFile_NoSyst = TFile(os.path.join(args.output_path,thetaName.replace('.root','_NoSyst.root')), 'RECREATE')

                thetaFile_NoSyst.cd()
                thetaData.Write()
                thetaBkg.Write()
                thetaSignal.Write()

                thetaFile_NoSyst.Close()
                thetaFile.cd()

                # create theta analysis file
                theta_content = theta_template
                theta_content = re.sub('DUMMY_JOB',args.final_state + '_m' + str(mass) + '_NoSyst',theta_content)
                theta_content = re.sub('model.add_lognormal_uncertainty','#model.add_lognormal_uncertainty',theta_content)

                theta_file = open(os.path.join(args.output_path,thetaName.replace('.root','_NoSyst.py')),'w')
                theta_file.write(theta_content)
                theta_file.close()

            if args.jesUnc != None:
                thetaSignal_JESUp = hSig_Syst_DataHist['JESUp'].createHistogram('dijet',mjj,RooFit.Binning(args.massMax-args.massMin,float(args.massMin),float(args.massMax)))
                thetaSignal_JESUp.SetName('dijet__signal__JES__up')
                thetaSignal_JESUp.Write()

                thetaSignal_JESDown = hSig_Syst_DataHist['JESDown'].createHistogram('dijet',mjj,RooFit.Binning(args.massMax-args.massMin,float(args.massMin),float(args.massMax)))
                thetaSignal_JESDown.SetName('dijet__signal__JES__down')
                thetaSignal_JESDown.Write()

            if args.jerUnc != None:
                thetaSignal_JERUp = hSig_Syst_DataHist['JERUp'].createHistogram('dijet',mjj,RooFit.Binning(args.massMax-args.massMin,float(args.massMin),float(args.massMax)))
                thetaSignal_JERUp.SetName('dijet__signal__JER__up')
                thetaSignal_JERUp.Write()

                thetaSignal_JERDown = hSig_Syst_DataHist['JERDown'].createHistogram('dijet',mjj,RooFit.Binning(args.massMax-args.massMin,float(args.massMin),float(args.massMax)))
                thetaSignal_JERDown.SetName('dijet__signal__JER__down')
                thetaSignal_JERDown.Write()

            thetaFile.Close()

            # create theta analysis file
            theta_content = theta_template
            theta_content = re.sub('DUMMY_JOB',args.final_state + '_m' + str(mass),theta_content)

            theta_file = open(os.path.join(args.output_path,thetaName.replace('.root','.py')),'w')
            theta_file.write(theta_content)
            theta_file.close()

    print '>> Datacards and workspaces created and stored in %s/.'%( os.path.join(os.getcwd(),args.output_path) )


if __name__ == '__main__':
    main()

