#!/usr/bin/env python

import sys, os, copy, re
from argparse import ArgumentParser
from array import array
import math
import ROOT
ROOT.gROOT.SetBatch(True)

import CMSDIJET.StatisticalTools.limit_configuration as limit_config
from CMSDIJET.StatisticalTools.systematics import *
from CMSDIJET.StatisticalTools.roofit_functions import *
import signal_fits

ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/RooCBPlusVoigtian.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/CMSDIJET/QCDAnalysis")
import analysis_configuration_8TeV as analysis_config

#if not os.path.exists(os.path.join(os.environ["CMSSW_BASE"],"src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer_cc.so")):
#    ROOT.gROOT.ProcessLine(".L " + os.path.join(os.environ["CMSSW_BASE"],"src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer.cc")+"+")
#else:
#    ROOT.gSystem.Load(os.path.join(os.environ["CMSSW_BASE"],"src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer_cc.so"))


def main():
    # usage description
    usage = "Example: ./scripts/createDatacards.py --inputData inputs/rawhistV7_Run2015D_scoutingPFHT_UNBLINDED_649_838_JEC_HLTplusV7_Mjj_cor_smooth.root --dataHistname mjj_mjjcor_gev --inputSig inputs/ResonanceShapes_gg_13TeV_Scouting_Spring15.root -f gg -o datacards -l 1866 --lumiUnc 0.027 --massrange 1000 1500 50 --runFit --p1 5 --p2 7 --p3 0.4 --massMin 838 --massMax 2037 --fitStrategy 2"

    # input parameters
    parser = ArgumentParser(description='Script that creates combine datacards and corresponding RooFit workspaces',epilog=usage)
    parser.add_argument("analysis", type=str, help="Analysis name")
    parser.add_argument("model", type=str, help="Model (Hbb, RSG)")

    #parser.add_argument("--inputData", dest="inputData", required=True,
    #                    help="Input data spectrum",
    #                    metavar="INPUT_DATA")

    parser.add_argument("--dataHistname", dest="dataHistname", type=str, default="h_data",
                        help="Data histogram name",
                        metavar="DATA_HISTNAME")

    #parser.add_argument("--inputSig", dest="inputSig", required=True,
    #                    help="Input signal shapes",
    #                    metavar="INPUT_SIGNAL")

    parser.add_argument("-f", "--final_state", dest="final_state", default="qq",
                        help="Final state (e.g. qq, qg, gg)",
                        metavar="FINAL_STATE")
    parser.add_argument("--fit_functions", dest="fit_functions", default="f1,f2,f3,f4,f5", help="List of fit functions")

    #parser.add_argument("-f2", "--type", dest="atype", required=True, help="Type (e.g. hG, lG, hR, lR)")

    parser.add_argument("-o", "--output_path", dest="output_path",
                        help="Output path where datacards and workspaces will be stored. If not specified, this is derived from limit_configuration.",
                        metavar="OUTPUT_PATH")

    parser.add_argument("-l", "--lumi", dest="lumi",
                        default=19700., type=float,
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
    parser.add_argument("--fitSignal", action="store_true", help="Use signal fitted shapes (CB+Voigtian) instead of histogram templates")
    #parser.add_argument("--lumiUnc", dest="lumiUnc",
    #                    required=True, type=float,
    #                    help="Relative uncertainty in the integrated luminosity",
    #                    metavar="LUMI_UNC")

    #parser.add_argument("--jesUnc", dest="jesUnc",
    #                    type=float,
    #                    help="Relative uncertainty in the jet energy scale",
    #                    metavar="JES_UNC")

    #parser.add_argument("--jerUnc", dest="jerUnc",
    #                    type=float,
    #                    help="Relative uncertainty in the jet energy resolution",
    #                    metavar="JER_UNC")

    parser.add_argument("--sqrtS", dest="sqrtS",
                        default=8000., type=float,
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

    fit_functions = args.fit_functions.split(",")

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
    inputData = TFile(limit_config.get_data_input(args.analysis))
    # input data histogram
    hData = inputData.Get(args.dataHistname)
    hData.SetDirectory(0)

    # input sig file
    print "[create_datacards] INFO : Opening resonance shapes file at " + limit_config.get_resonance_shapes(args.analysis, args.model)
    inputSig = TFile(limit_config.get_resonance_shapes(args.analysis, args.model), "READ")

    sqrtS = args.sqrtS

    # mass variable
    mjj = RooRealVar('mjj','mjj',float(args.massMin),float(args.massMax))

    # integrated luminosity and signal cross section
    lumi = args.lumi
    signalCrossSection = 1. # set to 1. so that the limit on r can be interpreted as a limit on the signal cross section

    for mass in masses:

        print ">> Creating datacard and workspace for %s resonance with m = %i GeV..."%(args.final_state, int(mass))
        
        rooDataHist = RooDataHist('rooDatahist','rooDathist',RooArgList(mjj),hData)

        hSig = inputSig.Get( "h_" + args.final_state + "_" + str(int(mass)) )
        if not hSig:
            raise Exception("Couldn't find histogram " + "h_" + args.final_state + "_" + str(int(mass)) + " in file " + limit_config.get_resonance_shapes(args.analysis, args.model))
        # normalize signal shape to the expected event yield (works even if input shapes are not normalized to unity)
        hSig.Scale(signalCrossSection*lumi/hSig.Integral()) # divide by a number that provides roughly an r value of 1-10
        rooSigHist = RooDataHist('rooSigHist','rooSigHist',RooArgList(mjj),hSig)
        print 'Signal acceptance:', (rooSigHist.sumEntries()/hSig.Integral())

        # If using fitted signal shapes, load the signal PDF
        if args.fitSignal:
            print "[create_datacards] Loading fitted signal PDFs from " + analysis_config.get_signal_fit_file(args.analysis, args.model, mass, "bukin", interpolated=(not mass in analysis_config.simulation.simulated_masses))
            f_signal_pdfs = TFile(analysis_config.get_signal_fit_file(args.analysis, args.model, mass, "bukin", interpolated=(not mass in analysis_config.simulation.simulated_masses)), "READ")
            w_signal = f_signal_pdfs.Get("w_signal")
            input_parameters = signal_fits.get_parameters(w_signal.pdf("signal"))

            # Make a new PDF with nuisance parameters
            signal_pdf, signal_vars = signal_fits.make_signal_pdf_systematic("bukin", mjj, mass=mass)

            # Copy input parameter values
            signal_vars["xp_0"].setVal(input_parameters["xp"][0])
            signal_vars["xp_0"].setError(input_parameters["xp"][1])
            signal_vars["xp_0"].setConstant()
            signal_vars["sigp_0"].setVal(input_parameters["sigp"][0])
            signal_vars["sigp_0"].setError(input_parameters["sigp"][1])
            signal_vars["sigp_0"].setConstant()
            signal_vars["xi_0"].setVal(input_parameters["xi"][0])
            signal_vars["xi_0"].setError(input_parameters["xi"][1])
            signal_vars["xi_0"].setConstant()
            signal_vars["rho1_0"].setVal(input_parameters["rho1"][0])
            signal_vars["rho1_0"].setError(input_parameters["rho1"][1])
            signal_vars["rho1_0"].setConstant()
            signal_vars["rho2_0"].setVal(input_parameters["rho2"][0])
            signal_vars["rho2_0"].setError(input_parameters["rho2"][1])
            signal_vars["rho2_0"].setConstant()
            f_signal_pdfs.Close()

        signal_parameters = {}
        signal_pdfs = {}
        signal_norms = {}
        background_pdfs = {}
        background_parameters = {}
        background_norms = {}
        signal_epdfs = {}
        background_epdfs = {}
        models = {}
        fit_results = {}

        for fit_function in fit_functions:
            if args.fitSignal:
                # Make a copy of the signal PDF, so that each fitTo call uses its own copy.
                # The copy should have all variables set constant.  
                #signal_pdfs[fit_function], signal_parameters[fit_function] = signal_fits.copy_signal_pdf("bukin", signal_pdf, mjj, tag=fit_function, include_systematics=True)
                signal_pdfs[fit_function] = ROOT.RooBukinPdf(signal_pdf, signal_pdf.GetName() + "_" + fit_function)
                iterator = signal_pdfs[fit_function].getVariables().createIterator()
                this_parameter = iterator.Next()
                while this_parameter:
                    this_parameter.setConstant()
                    this_parameter = iterator.Next()
            else:
                signal_pdfs[fit_function] = RooHistPdf('signal_' + fit_function,'signal_' + fit_function, RooArgSet(mjj), rooSigHist)
            signal_norms[fit_function] = RooRealVar('signal_norm_' + fit_function, 'signal_norm_' + fit_function, 0., 0., 1e+05)
            if args.fitBonly: 
                signal_norms[fit_function].setConstant()
            background_pdfs[fit_function], background_parameters[fit_function] = make_background_pdf(fit_function, mjj, collision_energy=8000.)

            # Initial values
            if "trigbbh" in args.analysis:
                if fit_function == "f3":
                    background_parameters[fit_function]["p1"].setVal(55.)
                    background_parameters[fit_function]["p2"].setVal(8.)
                elif fit_function == "f4":
                    background_parameters[fit_function]["p1"].setVal(28.)
                    background_parameters[fit_function]["p2"].setVal(-22.)
                    background_parameters[fit_function]["p3"].setVal(10.)
            elif "trigbbl" in args.analysis:
                if fit_function == "f3":
                    background_parameters[fit_function]["p1"].setVal(82.)
                    background_parameters[fit_function]["p2"].setVal(8.)
                elif fit_function == "f4":
                    background_parameters[fit_function]["p1"].setVal(41.)
                    background_parameters[fit_function]["p2"].setVal(-45.)
                    background_parameters[fit_function]["p3"].setVal(10.)

            data_integral = hData.Integral(hData.GetXaxis().FindBin(float(args.massMin)),hData.GetXaxis().FindBin(float(args.massMax)))
            background_norms[fit_function] = RooRealVar('background_' + fit_function + '_norm', 'background_' + fit_function + '_norm', data_integral, 0., 1.e8)

            signal_epdfs[fit_function] = RooExtendPdf('esignal_' + fit_function, 'esignal_' + fit_function, signal_pdfs[fit_function], signal_norms[fit_function])
            background_epdfs[fit_function] = RooExtendPdf('ebackground_' + fit_function, 'ebackground_' + fit_function, background_pdfs[fit_function], background_norms[fit_function])

            models[fit_function] = RooAddPdf('model_' + fit_function, 's+b', RooArgList(background_epdfs[fit_function], signal_epdfs[fit_function]))

            if args.runFit:
                fit_results[fit_function] = models[fit_function].fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Extended(kTRUE), RooFit.Strategy(args.fitStrategy))
                fit_results[fit_function].Print()

            # needed if want to evaluate limits without background systematics
            if args.fixBkg:
                background_norms[fit_function].setConstant()
                for par_name, par in background_parameters[fit_function].iteritems():
                    par.setConstant()

        # -----------------------------------------
        # dictionaries holding systematic variations of the signal shape
        hSig_Syst = {}
        hSig_Syst_DataHist = {}
        sigCDF = TGraph(hSig.GetNbinsX()+1)
        #signal_pdfs_syst = {}
        # JES and JER uncertainties
        if args.fitSignal:
            print "[create_datacards] INFO : Getting signal PDFs from " + analysis_config.get_signal_fit_file(args.analysis, args.model, mass, "bukin", interpolated=(not mass in analysis_config.simulation.simulated_masses))
            f_signal_pdfs = TFile(analysis_config.get_signal_fit_file(args.analysis, args.model, mass, "bukin", interpolated=(not mass in analysis_config.simulation.simulated_masses)))
            w_signal = f_signal_pdfs.Get("w_signal")
            if "jes" in systematics:
                xp_central = signal_vars["xp_0"].getVal()
                #print w_signal.pdf("signal__JESUp")
                #print signal_fits.get_parameters(w_signal.pdf("signal__JESUp"))
                xp_up = signal_fits.get_parameters(w_signal.pdf("signal__JESUp"))["xpJESUp"][0]
                xp_down = signal_fits.get_parameters(w_signal.pdf("signal__JESDown"))["xpJESDown"][0]
                signal_vars["dxp"].setVal(max(abs(xp_up - xp_central), abs(xp_down - xp_central)))
                signal_vars["alpha_jes"].setVal(0.)
                signal_vars["alpha_jes"].setConstant(False)
            else:
                signal_vars["dxp"].setVal(0.)
                signal_vars["alpha_jes"].setVal(0.)
                signal_vars["alpha_jes"].setConstant()
            signal_vars["dxp"].setError(0.)
            signal_vars["dxp"].setConstant()

            if "jer" in systematics:
                sigp_central = signal_vars["sigp_0"].getVal()
                sigp_up = signal_fits.get_parameters(w_signal.pdf("signal__JERUp"))["sigpJERUp"][0]
                sigp_down = signal_fits.get_parameters(w_signal.pdf("signal__JERDown"))["sigpJERDown"][0]
                signal_vars["dsigp"].setVal(max(abs(sigp_up - sigp_central), abs(sigp_down - sigp_central)))
                signal_vars["alpha_jer"].setVal(0.)
                signal_vars["alpha_jer"].setConstant(False)
            else:
                signal_vars["dsigp"].setVal(0.)
                signal_vars["alpha_jer"].setVal(0.)
                signal_vars["alpha_jer"].setConstant()
            signal_vars["dsigp"].setError(0.)
            signal_vars["dsigp"].setConstant()
                #for variation in ["JERUp", "JERDown"]:
                #    signal_pdfs_syst[variation] = w_signal.pdf("signal__" + variation)
            #for variation, pdf in signal_pdfs_syst.iteritems():
            #    signal_parameters = pdf.getVariables()
            #    iter = signal_parameters.createIterator()
            #    var = iter.Next()
            #    while var:
            #        var.setConstant()
            #        var = iter.Next()
            f_signal_pdfs.Close()
        else:
            if "jes" in systematics or "jer" in systematics:

                sigCDF.SetPoint(0,0.,0.)
                integral = 0.
                for i in range(1, hSig.GetNbinsX()+1):
                    x = hSig.GetXaxis().GetBinLowEdge(i+1)
                    integral = integral + hSig.GetBinContent(i)
                    sigCDF.SetPoint(i,x,integral)

            if "jes" in systematics:
                hSig_Syst['JESUp'] = copy.deepcopy(hSig)
                hSig_Syst['JESDown'] = copy.deepcopy(hSig)

            if "jer" in systematics:
                hSig_Syst['JERUp'] = copy.deepcopy(hSig)
                hSig_Syst['JERDown'] = copy.deepcopy(hSig)

            # reset signal histograms
            for key in hSig_Syst.keys():
                hSig_Syst[key].Reset()
                hSig_Syst[key].SetName(hSig_Syst[key].GetName() + '_' + key)

            # produce JES signal shapes
            if "jes" in systematics:
                for i in range(1, hSig.GetNbinsX()+1):
                    xLow = hSig.GetXaxis().GetBinLowEdge(i)
                    xUp = hSig.GetXaxis().GetBinLowEdge(i+1)
                    jes = 1. - systematics["jes"]
                    xLowPrime = jes*xLow
                    xUpPrime = jes*xUp
                    hSig_Syst['JESUp'].SetBinContent(i, sigCDF.Eval(xUpPrime) - sigCDF.Eval(xLowPrime))
                    jes = 1. + systematics["jes"]
                    xLowPrime = jes*xLow
                    xUpPrime = jes*xUp
                    hSig_Syst['JESDown'].SetBinContent(i, sigCDF.Eval(xUpPrime) - sigCDF.Eval(xLowPrime))
                hSig_Syst_DataHist['JESUp'] = RooDataHist('hSig_JESUp','hSig_JESUp',RooArgList(mjj),hSig_Syst['JESUp'])
                hSig_Syst_DataHist['JESDown'] = RooDataHist('hSig_JESDown','hSig_JESDown',RooArgList(mjj),hSig_Syst['JESDown'])
            
            # produce JER signal shapes
            if "jer" in systematics:
                for i in range(1, hSig.GetNbinsX()+1):
                    xLow = hSig.GetXaxis().GetBinLowEdge(i)
                    xUp = hSig.GetXaxis().GetBinLowEdge(i+1)
                    jer = 1. - systematics["jer"]
                    xLowPrime = jer*(xLow-float(mass))+float(mass)
                    xUpPrime = jer*(xUp-float(mass))+float(mass)
                    hSig_Syst['JERUp'].SetBinContent(i, sigCDF.Eval(xUpPrime) - sigCDF.Eval(xLowPrime))
                    jer = 1. + systematics["jer"]
                    xLowPrime = jer*(xLow-float(mass))+float(mass)
                    xUpPrime = jer*(xUp-float(mass))+float(mass)
                    hSig_Syst['JERDown'].SetBinContent(i, sigCDF.Eval(xUpPrime) - sigCDF.Eval(xLowPrime))
                hSig_Syst_DataHist['JERUp'] = RooDataHist('hSig_JERUp','hSig_JERUp',RooArgList(mjj),hSig_Syst['JERUp'])
                hSig_Syst_DataHist['JERDown'] = RooDataHist('hSig_JERDown','hSig_JERDown',RooArgList(mjj),hSig_Syst['JERDown'])


        # -----------------------------------------
        # create a datacard and corresponding workspace
        postfix = (('_' + args.postfix) if args.postfix != '' else '')
        wsName = 'workspace_' + args.final_state + '_m' + str(mass) + postfix + '.root'

        w = RooWorkspace('w','workspace')
        if args.fitSignal:
            signal_pdf.SetName("signal")
            getattr(w,'import')(signal_pdf,RooFit.Rename("signal"))
            # Create a norm variable "signal_norm" which normalizes the PDF to unity.
            norm = args.lumi
            #signal_norm = ROOT.RooRealVar("signal_norm", "signal_norm", 1. / norm, 0.1 / norm, 10. / norm)
            signal_norm = ROOT.RooRealVar("signal_norm", "signal_norm", norm, norm / 10., norm * 10.)
            print "[create_datacards] INFO : Set signal norm to {}".format(signal_norm.getVal())
            signal_norm.setConstant()
            getattr(w,'import')(signal_norm,ROOT.RooCmdArg())
            #if "jes" in systematics:
            #    getattr(w,'import')(signal_pdfs_syst['JESUp'],RooFit.Rename("signal__JESUp"))
            #    getattr(w,'import')(signal_pdfs_syst['JESDown'],RooFit.Rename("signal__JESDown"))
            #if "jer" in systematics:
            #    getattr(w,'import')(signal_pdfs_syst['JERUp'],RooFit.Rename("signal__JERUp"))
            #    getattr(w,'import')(signal_pdfs_syst['JERDown'],RooFit.Rename("signal__JERDown"))
        else:
            getattr(w,'import')(rooSigHist,RooFit.Rename("signal"))
            if "jes" in systematics:
                getattr(w,'import')(hSig_Syst_DataHist['JESUp'],RooFit.Rename("signal__JESUp"))
                getattr(w,'import')(hSig_Syst_DataHist['JESDown'],RooFit.Rename("signal__JESDown"))
            if "jer" in systematics:
                getattr(w,'import')(hSig_Syst_DataHist['JERUp'],RooFit.Rename("signal__JERUp"))
                getattr(w,'import')(hSig_Syst_DataHist['JERDown'],RooFit.Rename("signal__JERDown"))
        if args.decoBkg:
            getattr(w,'import')(background_deco,ROOT.RooCmdArg())
        else:
            for fit_function in fit_functions:
                getattr(w,'import')(background_pdfs[fit_function],ROOT.RooCmdArg(),RooFit.Rename("background_" + fit_function))
                getattr(w,'import')(background_norms[fit_function],ROOT.RooCmdArg(),RooFit.Rename("background_" + fit_function + "_norm"))
                getattr(w,'import')(fit_results[fit_function])
                getattr(w,'import')(signal_norms[fit_function],ROOT.RooCmdArg())
        getattr(w,'import')(rooDataHist,RooFit.Rename("data_obs"))
        w.Print()
        if args.output_path:
            if not os.path.isdir( os.path.join(os.getcwd(),args.output_path) ):
                os.mkdir( os.path.join(os.getcwd(),args.output_path) )
            w.writeToFile(os.path.join(args.output_path,wsName))
        else:
            w.writeToFile(limit_config.get_workspace_filename(args.analysis, args.model, mass, fitBonly=args.fitBonly, fitSignal=args.fitSignal))


        beffUnc = 0.3
        boffUnc = 0.06

        for fit_function in fit_functions:
            if args.output_path:
                dcName = 'datacard_' + args.final_state + '_m' + str(mass) + postfix + '_' + fit_function + '.txt'
                datacard = open(os.path.join(args.output_path,dcName),'w')
            else:
                datacard = open(limit_config.get_datacard_filename(args.analysis, args.model, mass, fit_function, fitSignal=args.fitSignal), 'w')
            datacard.write('imax 1\n')
            datacard.write('jmax 1\n')
            datacard.write('kmax *\n')
            datacard.write('---------------\n')
            if ("jes" in systematics or "jer" in systematics) and not args.fitSignal:
                if args.output_path:
                    datacard.write('shapes * * '+wsName+' w:$PROCESS w:$PROCESS__$SYSTEMATIC\n')
                else:
                    datacard.write('shapes * * '+os.path.basename(limit_config.get_workspace_filename(args.analysis, args.model, mass, fitSignal=args.fitSignal))+' w:$PROCESS w:$PROCESS__$SYSTEMATIC\n')
            else:
                if args.output_path:
                    datacard.write('shapes * * '+wsName+' w:$PROCESS\n')
                else:
                    datacard.write('shapes * * '+os.path.basename(limit_config.get_workspace_filename(args.analysis, args.model, mass, fitSignal=args.fitSignal))+' w:$PROCESS\n')
            datacard.write('---------------\n')
            datacard.write('bin 1\n')
            datacard.write('observation -1\n')
            datacard.write('------------------------------\n')
            datacard.write('bin          1          1\n')
            datacard.write('process      signal     background_' + fit_function + '\n')
            datacard.write('process      0          1\n')
            if args.fitSignal:
                datacard.write('rate         1         1\n')
            else:
                datacard.write('rate         -1         1\n')
            datacard.write('------------------------------\n')
            datacard.write('lumi  lnN    %f         -\n'%(1.+systematics["luminosity"]))
            datacard.write('beff  lnN    %f         -\n'%(1.+beffUnc))
            datacard.write('boff  lnN    %f         -\n'%(1.+boffUnc))
            #datacard.write('bkg   lnN     -         1.03\n')
            if args.fitSignal:
                if "jes" in systematics:
                    datacard.write("alpha_jes  param  0.0  1.0\n")
                if "jer" in systematics:
                    datacard.write("alpha_jer  param  0.0  1.0\n")
            else:
                if "jes" in systematics:
                    datacard.write('JES  shape   1          -\n')
                if "jer" in systematics:
                    datacard.write('JER  shape   1          -\n')
            # flat parameters --- flat prior
            datacard.write('background_' + fit_function + '_norm  flatParam\n')
            if args.decoBkg:
                datacard.write('deco_eig1  flatParam\n')
                datacard.write('deco_eig2  flatParam\n')
            else:
                for par_name, par in background_parameters[fit_function].iteritems():
                    datacard.write(fit_function + "_" + par_name + '  flatParam\n')
            datacard.close()


    #print '>> Datacards and workspaces created and stored in %s/'%( os.path.join(os.getcwd(),args.output_path) )
    print "All done."

if __name__ == '__main__':
    main()

