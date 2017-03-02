#!/usr/bin/env python

import sys, os, copy, re
from argparse import ArgumentParser
from array import array
import math
import ROOT
ROOT.gROOT.SetBatch(True)

import CMSDIJET.StatisticalTools.limit_configuration as limit_config
import CMSDIJET.StatisticalTools.trigger_efficiency as trigger_efficiency
from CMSDIJET.StatisticalTools.systematics import *
from CMSDIJET.StatisticalTools.roofit_functions import *
from CMSDIJET.StatisticalTools.background_fits import *
import signal_fits


ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/RooCBPlusVoigtian.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so"))

sys.path.append(os.path.expandvars("$CMSSW_BASE/python/CMSDIJET/StatisticalTools"))
import analysis_configuration_8TeV as analysis_config

#if not os.path.exists(os.path.join(os.environ["CMSSW_BASE"],"src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer_cc.so")):
#    ROOT.gROOT.ProcessLine(".L " + os.path.join(os.environ["CMSSW_BASE"],"src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer.cc")+"+")
#else:
#    ROOT.gSystem.Load(os.path.join(os.environ["CMSSW_BASE"],"src/HiggsAnalysis/CombinedLimit/src/PdfDiagonalizer_cc.so"))


def main():
    # usage description
    usage = "Example: ./scripts/createDatacards.py --inputData inputs/rawhistV7_Run2015D_scoutingPFHT_UNBLINDED_649_838_JEC_HLTplusV7_Mjj_cor_smooth.root --dataHistname mjj_mjjcor_gev --inputSig inputs/ResonanceShapes_gg_13TeV_Scouting_Spring15.root -f gg -o datacards -l 1866 --lumiUnc 0.027 --massrange 1000 1500 50 --runFit --p1 5 --p2 7 --p3 0.4 --massMin 838 --massMax 2037 --fitStrategy 2"

    print "Welcome to create_datacards_parallel"

    # input parameters
    parser = ArgumentParser(description='Script that creates combine datacards and corresponding RooFit workspaces',epilog=usage)
    parser.add_argument("analysis", type=str, help="Analysis name")
    parser.add_argument("model", type=str, help="Model (Hbb, RSG)")
    parser.add_argument("--condor", action="store_true", help="For running on condor, where the file paths are different.")
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
    trigger_group = parser.add_mutually_exclusive_group()
    trigger_group.add_argument("--fitTrigger", dest="fitTrigger",
                        action='store_true',
                        help="Include trigger correction in PDF. Exclusive with correctTrigger")

    trigger_group.add_argument("--correctTrigger", dest="correctTrigger",
                        action='store_true',
                        help="Apply trigger correction to data. Exclusive with fitTrigger.")
    parser.add_argument("--useMCTrigger", action="store_true", help="Use MC trigger simulation (i.e. trigbbl/trigbbh), instead of measured data trigger times *_notrig_*")
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

    parser.add_argument("--qcd", action="store_true", help="Fit QCD instead of data")
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

    if not (args.useMCTrigger or args.fitTrigger or args.correctTrigger or args.qcd):
        print "[create_datacards] ERROR : None of useMCTrigger, fitTrigger, and correctTrigger was specified. Therefore, I don't know what to do with the signal samples."
        sys.exit(1)

    masses = []
    if args.fitBonly:
        masses.append(750)
    else:
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

    if len(masses) == 1:
        run_single_mass(args, masses[0])
    else:
        from joblib import Parallel, delayed
        Parallel(n_jobs=4)(delayed(run_single_mass)(args, mass) for mass in masses)
    print "All done."

def run_single_mass(args, mass):
    print "[run_single_mass] INFO : Creating datacard and workspace for m = %i GeV..."%(int(mass))
    fit_functions = args.fit_functions.split(",")

    # import ROOT stuff
    from ROOT import gStyle, TFile, TH1F, TH1D, TGraph, kTRUE, kFALSE, TCanvas, TLegend, TPad, TLine
    from ROOT import RooHist, RooRealVar, RooDataHist, RooArgList, RooArgSet, RooAddPdf, RooProdPdf, RooEffProd, RooFit, RooGenericPdf, RooWorkspace, RooMsgService, RooHistPdf, RooExtendPdf, RooFormulaVar

    if not args.debug:
        RooMsgService.instance().setSilentMode(kTRUE)
        RooMsgService.instance().setStreamStatus(0,kFALSE)
        RooMsgService.instance().setStreamStatus(1,kFALSE)

    # Stuff
    mjj = RooRealVar('mjj','mjj',float(args.massMin),float(args.massMax))
    lumi = args.lumi
    signalCrossSection = 1. # Set to 1 so that the limit on r can be interpreted as a limit on the signal cross section

    # Input data file
    if args.qcd:
        data_sample = "QCD_TuneZ2star_8TeV_pythia6"
    elif "trigbb" in args.analysis:
        data_sample = "BJetPlusX_2012"
    elif "trigmu" in args.analysis:
        data_sample = "SingleMu_2012"
    data_file_path = analysis_config.get_b_histogram_filename(args.analysis, data_sample)
    if args.condor:
        data_file_path = os.path.basename(data_file_path)
    data_file = TFile(data_file_path)
    hData_notrigcorr = data_file.Get("BHistograms/h_pfjet_mjj")
    hData_notrigcorr.SetDirectory(0)
    hData_name = hData_notrigcorr.GetName()
    hData_notrigcorr.SetName(hData_name + "_notrigcorr")

    # Trigger correction on data
    if args.fitTrigger:
        # Make trigger correction objects
        trigeff_pt_formula, trigeff_vars = trigger_efficiency.get_var_formula(args.analysis, mjj)
        trigeff_btag_var      = RooRealVar("trigeff_btag", "trigeff_btag", 0., 1.)
        trigeff_btag_var.setVal(trigger_efficiency.online_btag_eff[args.analysis][0])
        trigeff_btag_var.setVal(trigger_efficiency.online_btag_eff[args.analysis][0])
        trigeff_btag_formula  = RooFormulaVar("trigeff_btag_formula", "@0", RooArgList(trigeff_btag_var))
        #trigeff_btag_var.setConstant()
        trigeff_total_formula = RooFormulaVar("trigeff_total_formula", "@0*@1", RooArgList(trigeff_btag_var, trigeff_pt_formula))
        #for trigeff_var_name, trigeff_var in trigeff_vars.iteritems():
        #    trigeff_var.setConstant()
    if args.correctTrigger:
        # Apply trigger correction to data histogram
        hData = CorrectTriggerEfficiency(hData_notrigcorr, args.analysis)

        # Still need b-tagging efficiency to scale the MC
        trigeff_btag_formula  = RooFormulaVar("trigeff_btag_formula", str(trigger_efficiency.online_btag_eff[args.analysis][0]), RooArgList())
    else:
        hData = hData_notrigcorr
    if args.qcd:
        # For QCD, scale the histogram by the online b-tag trigger efficiency.
        if args.analysis == "NoTrigger_eta1p7_CSVTM":
            trigeff_btag_formula  = RooFormulaVar("trigeff_btag_formula", str(trigger_efficiency.online_btag_eff["trigbbl_CSVTM"][0]), RooArgList())
        elif args.analysis == "NoTrigger_eta2p2_CSVTM":
            trigeff_btag_formula  = RooFormulaVar("trigeff_btag_formula", str(trigger_efficiency.online_btag_eff["trigbbh_CSVTM"][0]), RooArgList())
        if args.analysis == "NoTrigger_eta1p7_CSVTM":
            hData.Scale(trigger_efficiency.online_btag_eff["trigbbl_CSVTM"][0])
        elif args.analysis == "NoTrigger_eta2p2_CSVTM":
            hData.Scale(trigger_efficiency.online_btag_eff["trigbbh_CSVTM"][0])
        else:
            print "[create_datacards_parallel] ERROR : QCD fit requested, but analysis != NoTrigger_etaXpY_CSVTM. I don't know what to do!"
            sys.exit(1)
    hData.SetName(hData_name)
    
    rooDataHist = RooDataHist('rooDatahist','rooDatahist',RooArgList(mjj),hData)
    if args.correctTrigger:
        rooDataHist_notrigcorr = RooDataHist("rooDatahist_notrigcorr", "rooDatahist_notrigcorr", RooArgList(mjj), hData_notrigcorr)

    # Get signal pdf * btag efficiency
    if args.useMCTrigger:
        signal_pdf_file = analysis_config.get_signal_fit_file(args.analysis, args.model, mass, "bukin", interpolated=(not mass in analysis_config.simulation.simulated_masses))
        if args.condor:
            signal_pdf_file = os.path.basename(signal_pdf_file)
        print "[create_datacards] Loading fitted signal PDFs from " + signal_pdf_file
        f_signal_pdfs = TFile(signal_pdf_file, "READ")
        w_signal = f_signal_pdfs.Get("w_signal")
        bukin_pdf = w_signal.pdf("signal_bukin")
    else:
        if args.analysis == "trigbbl_CSVTM" or args.analysis == "NoTrigger_eta1p7_CSVTM":
            notrig_analysis = "trigbbl_notrig_CSVTM"
        elif args.analysis == "trigbbh_CSVTM" or args.analysis == "NoTrigger_eta2p2_CSVTM":
            notrig_analysis = "trigbbh_notrig_CSVTM"
        else:
            print "[run_single_mass] ERROR : I don't know a no-trigger variant of analysis {}. Please make one, or specify --useMCTrigger.".format(args.analysis) 
            sys.exit(1)
        signal_pdf_file = analysis_config.get_signal_fit_file(notrig_analysis, args.model, mass, "bukin", interpolated=(not mass in analysis_config.simulation.simulated_masses))
        print "[create_datacards] Loading fitted signal PDFs from " + signal_pdf_file
        if args.condor:
            signal_pdf_file = os.path.basename(signal_pdf_file)
        f_signal_pdfs = TFile(signal_pdf_file, "READ")
        w_signal = f_signal_pdfs.Get("w_signal")
        bukin_pdf = w_signal.pdf("signal")
        bukin_pdf.SetName("signal_bukin")
    input_signal_parameters = signal_fits.get_parameters(bukin_pdf)

    # Make a new PDF with nuisance parameters
    signal_pdf_notrig, signal_vars = signal_fits.make_signal_pdf_systematic("bukin", mjj, mass=mass)
    signal_pdf_name = signal_pdf_notrig.GetName()
    signal_pdf_notrig.SetName(signal_pdf_name + "_notrig")

    # Add trigger efficiency
    if args.useMCTrigger:
        # Signal PDF = bukin; b-tag efficiency already included in normalization
        signal_pdf = signal_pdf_notrig
        signal_pdf.SetName(signal_pdf_name)
    else:
        if args.fitTrigger:
            # Signal PDF = bukin * total trigger efficiency
            signal_pdf = RooEffProd("signal", "signal", signal_pdf_notrig, trigeff_total_formula)
        elif args.correctTrigger:
            # Signal PDF = bukin * btag efficiency
            signal_pdf = RooEffProd("signal", "signal", signal_pdf_notrig, trigeff_btag_formula)
        elif args.qcd:
            # Signal PDF = bukin * btag efficiency
            # Same as correctTrigger
            signal_pdf = RooEffProd("signal", "signal", signal_pdf_notrig, trigeff_btag_formula)
    # Copy input parameter values
    signal_vars["xp_0"].setVal(input_signal_parameters["xp"][0])
    signal_vars["xp_0"].setError(input_signal_parameters["xp"][1])
    signal_vars["xp_0"].setConstant()
    signal_vars["sigp_0"].setVal(input_signal_parameters["sigp"][0])
    signal_vars["sigp_0"].setError(input_signal_parameters["sigp"][1])
    signal_vars["sigp_0"].setConstant()
    signal_vars["xi_0"].setVal(input_signal_parameters["xi"][0])
    signal_vars["xi_0"].setError(input_signal_parameters["xi"][1])
    signal_vars["xi_0"].setConstant()
    signal_vars["rho1_0"].setVal(input_signal_parameters["rho1"][0])
    signal_vars["rho1_0"].setError(input_signal_parameters["rho1"][1])
    signal_vars["rho1_0"].setConstant()
    signal_vars["rho2_0"].setVal(input_signal_parameters["rho2"][0])
    signal_vars["rho2_0"].setError(input_signal_parameters["rho2"][1])
    signal_vars["rho2_0"].setConstant()
    f_signal_pdfs.Close()

    signal_parameters = {}
    signal_pdfs_notrig = {}
    signal_pdfs = {}
    signal_norms = {}
    background_pdfs = {}
    background_pdfs_notrig = {}
    background_parameters = {}
    background_norms = {}
    signal_epdfs = {}
    background_epdfs = {}
    models = {}
    fit_results = {}

    for fit_function in fit_functions:
        print "[create_datacards] INFO : On fit function {}".format(fit_function)

        # Make a copy of the signal PDF, so that each fitTo call uses its own copy.
        # The copy should have all variables set constant.  
        #signal_pdfs[fit_function], signal_parameters[fit_function] = signal_fits.copy_signal_pdf("bukin", signal_pdf, mjj, tag=fit_function, include_systematics=True)

        signal_pdfs_notrig[fit_function] = ROOT.RooBukinPdf(signal_pdf_notrig, signal_pdf_notrig.GetName() + "_" + fit_function)
        iterator = signal_pdfs_notrig[fit_function].getVariables().createIterator()
        this_parameter = iterator.Next()
        while this_parameter:
            this_parameter.setConstant()
            this_parameter = iterator.Next()

        # Add trigger efficiency
        if args.useMCTrigger:
            # Signal PDF = bukin; b-tag efficiency already included in normalization
            signal_pdfs[fit_function] = signal_pdfs_notrig[fit_function]
            signal_pdfs[fit_function].SetName(signal_pdf.GetName() + "_" + fit_function)
        elif args.fitTrigger:
            # Signal PDF = bukin * total trigger efficiency
            signal_pdfs[fit_function] = RooEffProd(signal_pdf.GetName() + "_" + fit_function, signal_pdf.GetName() + "_" + fit_function, signal_pdfs_notrig[fit_function], trigeff_total_formula)
        elif args.correctTrigger:
            # Signal PDF = bukin * btag efficiency
            signal_pdfs[fit_function] = RooEffProd(signal_pdf.GetName() + "_" + fit_function, signal_pdf.GetName() + "_" + fit_function, signal_pdfs_notrig[fit_function], trigeff_btag_formula)
        elif args.qcd:
            # Signal PDF = bukin * btag efficiency
            signal_pdfs[fit_function] = RooEffProd(signal_pdf.GetName() + "_" + fit_function, signal_pdf.GetName() + "_" + fit_function, signal_pdfs_notrig[fit_function], trigeff_btag_formula)
        signal_norms[fit_function] = RooRealVar('signal_norm_' + fit_function, 'signal_norm_' + fit_function, 0., 0., 1e+05)
        if args.fitBonly: 
            signal_norms[fit_function].setConstant()

        # Make background PDF
        background_pdfs_notrig[fit_function], background_parameters[fit_function] = make_background_pdf(fit_function, mjj, collision_energy=8000.)
        background_pdf_name = background_pdfs_notrig[fit_function].GetName()
        background_pdfs_notrig[fit_function].SetName(background_pdf_name + "_notrig")
        if args.fitTrigger:
            background_pdfs[fit_function] = RooEffProd(background_pdf_name, background_pdf_name, background_pdfs_notrig[fit_function], trigeff_pt_formula)
        else:
            background_pdfs[fit_function] = background_pdfs_notrig[fit_function]
            background_pdfs[fit_function].SetName(background_pdf_name)
        
        # Initial values
        if "trigbbh" in args.analysis:
            if fit_function == "f1":
                background_parameters[fit_function]["p1"].setVal(-13.5877235358)
                background_parameters[fit_function]["p2"].setVal(14.0659901462)
                background_parameters[fit_function]["p3"].setVal(1.24550474025)
                background_parameters[fit_function]["p1"].setMin(-50.)
                background_parameters[fit_function]["p1"].setMax(50.)
            elif fit_function == "f2":
                background_parameters[fit_function]["p1"].setVal(6.06731321562)
                background_parameters[fit_function]["p2"].setVal(6.06264502704)
            elif fit_function == "f3":
                background_parameters[fit_function]["p1"].setVal(50.0270215343)
                background_parameters[fit_function]["p2"].setVal(8.17180937688)
                background_parameters[fit_function]["p1"].setMin(20.)
            elif fit_function == "f4":
                background_parameters[fit_function]["p1"].setVal(31.3765210572)
                background_parameters[fit_function]["p2"].setVal(-22.5800092219)
                background_parameters[fit_function]["p3"].setVal(9.94548656557)
            elif fit_function == "f5":
                background_parameters[fit_function]["p1"].setVal(5.51929170927)
                background_parameters[fit_function]["p2"].setVal(4.25521547671)
            elif fit_function == "f6":
                background_parameters[fit_function]["p1"].setVal(35.)
                background_parameters[fit_function]["p2"].setVal(-28.)
                background_parameters[fit_function]["p3"].setVal(0.)
                background_parameters[fit_function]["p4"].setVal(10.)
        elif "trigbbl" in args.analysis:
            if fit_function == "f1":
                background_parameters[fit_function]["p1"].setVal(-32.4727133488)
                background_parameters[fit_function]["p2"].setVal(18.7641649883)
                background_parameters[fit_function]["p3"].setVal(1.84028034937)
            elif fit_function == "f2":
                background_parameters[fit_function]["p1"].setVal(4.96261586452)
                background_parameters[fit_function]["p2"].setVal(19.0848105961)
            if fit_function == "f3":
                background_parameters[fit_function]["p1"].setVal(60.0000032579)
                background_parameters[fit_function]["p2"].setVal(8.00317534363)
                background_parameters[fit_function]["p1"].setMin(60.)
            elif fit_function == "f4":
                background_parameters[fit_function]["p1"].setVal(25.4109169544)
                background_parameters[fit_function]["p2"].setVal(-42.56719661)
                background_parameters[fit_function]["p3"].setVal(12.3295648189)
            elif fit_function == "f5":                
                background_parameters[fit_function]["p1"].setVal(3.74859358646)
                background_parameters[fit_function]["p2"].setVal(11.4366903839)
            elif fit_function == "f6":
                background_parameters[fit_function]["p1"].setVal(35.)
                background_parameters[fit_function]["p2"].setVal(-43.)
                background_parameters[fit_function]["p3"].setVal(0.)
                background_parameters[fit_function]["p4"].setVal(10.)

        data_integral = hData.Integral(hData.GetXaxis().FindBin(float(args.massMin)),hData.GetXaxis().FindBin(float(args.massMax)))
        background_norms[fit_function] = RooRealVar('background_' + fit_function + '_norm', 'background_' + fit_function + '_norm', data_integral, 0., 1.e8)

        signal_epdfs[fit_function] = RooExtendPdf('esignal_' + fit_function, 'esignal_' + fit_function, signal_pdfs[fit_function], signal_norms[fit_function])
        background_epdfs[fit_function] = RooExtendPdf('ebackground_' + fit_function, 'ebackground_' + fit_function, background_pdfs[fit_function], background_norms[fit_function])

        models[fit_function] = RooAddPdf('model_' + fit_function, 's+b', RooArgList(background_epdfs[fit_function], signal_epdfs[fit_function]))

        if args.runFit:
            print "[create_datacards] INFO : Starting fit with function {}".format(fit_function)
            models[fit_function].Print()
            # Fix the trigger efficiency for this fit
            if args.fitTrigger:
                for var_name, var in trigeff_vars.iteritems():
                    var.setConstant(True)
            if args.fitTrigger:
                trigeff_btag_var.setConstant(True)
            fit_results[fit_function] = models[fit_function].fitTo(rooDataHist, RooFit.Save(kTRUE), RooFit.Extended(kTRUE), RooFit.Strategy(args.fitStrategy), RooFit.Verbose(0))
            print "[create_datacards] INFO : Done with fit {}. Printing results.".format(fit_function)
            fit_results[fit_function].Print()
            if args.fitTrigger:
                for var_name, var in trigeff_vars.iteritems():
                    var.setConstant(False)
            if args.fitTrigger:
                trigeff_btag_var.setConstant(False)
            print "[create_datacards] DEBUG : End args.runFit if block."


        # needed if want to evaluate limits without background systematics
        if args.fixBkg:
            background_norms[fit_function].setConstant()
            for par_name, par in background_parameters[fit_function].iteritems():
                par.setConstant()

    # -----------------------------------------
    # Set values of signal systematic variables
    # JES and JER uncertainties
    if "jes" in systematics:
        xp_central = signal_vars["xp_0"].getVal()
        xp_up = signal_fits.get_parameters(w_signal.pdf("signal__JESUp"))["xpJESUp"][0]
        xp_down = signal_fits.get_parameters(w_signal.pdf("signal__JESDown"))["xpJESDown"][0]
        signal_vars["dxp"].setVal(max(abs(xp_up - xp_central), abs(xp_down - xp_central)))
        if signal_vars["dxp"].getVal() > 2 * mass * 0.1:
            print "[create_datacards] WARNING : Large dxp value. dxp = {}, xp_down = {}, xp_central = {}, xp_up = {}".format(signal_vars["dxp"].getVal(), xp_down, xp_central, xp_up)
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

    # -----------------------------------------
    # create a datacard and corresponding workspace
    postfix = (('_' + args.postfix) if args.postfix != '' else '')
    wsName = 'workspace_' + args.final_state + '_m' + str(mass) + postfix + '.root'

    w = RooWorkspace('w','workspace')
    signal_pdf.SetName("signal")
    getattr(w,'import')(signal_pdf,RooFit.Rename("signal"))
    norm = args.lumi
    signal_norm = ROOT.RooRealVar("signal_norm", "signal_norm", norm/100., norm/100. / 10., norm * 10.)
    print "[create_datacards] INFO : Set signal norm to {}".format(signal_norm.getVal())
    signal_norm.setConstant()
    getattr(w,'import')(signal_norm,ROOT.RooCmdArg())
    for fit_function in fit_functions:
        print "Importing background PDF"
        print background_pdfs[fit_function]
        background_pdfs[fit_function].Print()
        getattr(w,'import')(background_pdfs[fit_function],ROOT.RooCmdArg(),RooFit.Rename("background_" + fit_function), RooFit.RecycleConflictNodes())
        w.pdf("background_" + fit_function).Print()
        getattr(w,'import')(background_norms[fit_function],ROOT.RooCmdArg(),RooFit.Rename("background_" + fit_function + "_norm"))
        getattr(w,'import')(fit_results[fit_function])
        getattr(w,'import')(signal_norms[fit_function],ROOT.RooCmdArg())
        if args.fitBonly:
            getattr(w,'import')(models[fit_function],ROOT.RooCmdArg(),RooFit.RecycleConflictNodes())
    getattr(w,'import')(rooDataHist,RooFit.Rename("data_obs"))
    if args.correctTrigger:
        getattr(w,'import')(rooDataHist_notrigcorr, RooFit.Rename("data_obs_notrigcorr"))

    w.Print()
    print "Starting save"
    if args.output_path:
        if not os.path.isdir( os.path.join(os.getcwd(),args.output_path) ):
            os.mkdir( os.path.join(os.getcwd(),args.output_path) )
        workspace_output_path = os.path.join(args.output_path,wsName)
    else:
        workspace_output_path = limit_config.get_workspace_filename(args.analysis, args.model, mass, fitBonly=args.fitBonly, fitTrigger=args.fitTrigger, correctTrigger=args.correctTrigger, useMCTrigger=args.useMCTrigger, qcd=args.qcd)
    if args.condor:
        workspace_output_path = os.path.basename(workspace_output_path)
    print "[create_datacards] INFO : Writing workspace to file {}".format(workspace_output_path)
    w.writeToFile(workspace_output_path)
    if args.correctTrigger:
        f_workspace = TFile(workspace_output_path, "UPDATE")
        hData.Write()
        hData_notrigcorr.Write()

    # Clean up
    for name, obj in signal_norms.iteritems():
        if obj:
            obj.IsA().Destructor(obj)
    for name, obj in background_pdfs.iteritems():
        if obj:
            obj.IsA().Destructor(obj)
    for name, obj in background_pdfs_notrig.iteritems():
        if obj:
            obj.IsA().Destructor(obj)
    for name, obj in background_norms.iteritems():
        if obj:
            obj.IsA().Destructor(obj)
    for name, obj in signal_pdfs.iteritems():
        if obj:
            obj.IsA().Destructor(obj)
    for name, obj in signal_pdfs_notrig.iteritems():
        if obj:
            obj.IsA().Destructor(obj)
    for name, obj in signal_epdfs.iteritems():
        if obj:
            obj.IsA().Destructor(obj)
    for name, obj in background_epdfs.iteritems():
        if obj:
            obj.IsA().Destructor(obj)
    for name, obj in fit_results.iteritems():
        if obj:
            obj.IsA().Destructor(obj)
    for name, dict_l2 in background_parameters.iteritems():
        for name2, obj in dict_l2.iteritems():
            if obj:
                obj.IsA().Destructor(obj)
    for name, obj in models.iteritems():
        if obj:
            obj.IsA().Destructor(obj)
    rooDataHist.IsA().Destructor(rooDataHist)
    w.IsA().Destructor(w)

    # Make datacards only if S+B fitted
    if not args.fitBonly:
        #beffUnc = 0.3
        boffUnc = 0.06
        for fit_function in fit_functions:
            if args.output_path:
                dcName = 'datacard_' + args.final_state + '_m' + str(mass) + postfix + '_' + fit_function + '.txt'
                datacard_output_path = os.path.join(args.output_path,dcName)
            else:
                datacard_output_path = limit_config.get_datacard_filename(args.analysis, args.model, mass, fit_function, fitTrigger=args.fitTrigger, correctTrigger=args.correctTrigger, qcd=args.qcd)
            if args.condor:
                datacard_output_path = os.path.basename(datacard_output_path)
            print "[create_datacards] INFO : Writing datacard to file {}".format(datacard_output_path) 
            datacard = open(datacard_output_path, 'w')
            datacard.write('imax 1\n')
            datacard.write('jmax 1\n')
            datacard.write('kmax *\n')
            datacard.write('---------------\n')
            if args.output_path:
                datacard.write('shapes * * '+wsName+' w:$PROCESS\n')
            else:
                datacard.write('shapes * * '+os.path.basename(limit_config.get_workspace_filename(args.analysis, args.model, mass, fitTrigger=args.fitTrigger, correctTrigger=args.correctTrigger, qcd=args.qcd))+' w:$PROCESS\n')
            datacard.write('---------------\n')
            datacard.write('bin 1\n')
            datacard.write('observation -1\n')
            datacard.write('------------------------------\n')
            datacard.write('bin          1          1\n')
            datacard.write('process      signal     background_' + fit_function + '\n')
            datacard.write('process      0          1\n')
            datacard.write('rate         1         1\n')
            datacard.write('------------------------------\n')
            datacard.write('lumi  lnN    %f         -\n'%(1.+systematics["luminosity"]))
            #datacard.write('beff  lnN    %f         -\n'%(1.+beffUnc))
            datacard.write('boff  lnN    %f         -\n'%(1.+boffUnc))
            #datacard.write('bkg   lnN     -         1.03\n')
            if "jes" in systematics:
                datacard.write("alpha_jes  param  0.0  1.0\n")
            if "jer" in systematics:
                datacard.write("alpha_jer  param  0.0  1.0\n")
            # Trigger efficiency
            if not args.useMCTrigger:
                datacard.write("trigeff_btag  param  {}  {}\n".format(trigger_efficiency.online_btag_eff[args.analysis][0], trigger_efficiency.online_btag_eff[args.analysis][1]))
            if args.fitTrigger:
                for trigeff_var_name, trigeff_var in trigeff_vars.iteritems():
                    datacard.write("{}  param  {}  {}\n".format(trigeff_var_name, trigger_efficiency.sigmoid_parameters[args.analysis][trigeff_var_name][0], trigger_efficiency.sigmoid_parameters[args.analysis][trigeff_var_name][1]))
            # Background fit parameters --- flat prior
            datacard.write('background_' + fit_function + '_norm  flatParam\n')
            for par_name, par in background_parameters[fit_function].iteritems():
                    datacard.write(fit_function + "_" + par_name + '  flatParam\n')
            datacard.close()
            print "[create_datacards] INFO : Done with this datacard"

    #print '>> Datacards and workspaces created and stored in %s/'%( os.path.join(os.getcwd(),args.output_path) )
    print "Done with mass {}.".format(mass)

def CorrectTriggerEfficiency(hist, analysis):
    hist_corr = hist.Clone()
    hist_corr.SetName(hist.GetName() + "_trigcorr")
    if "bbl" in analysis:
        mjj_min = 296
        mjj_max = 5000
    elif "bbh" in analysis:
        mjj_min = 526
        mjj_max = 5000
    bin_min = hist_corr.GetXaxis().FindBin(mjj_min + 1.e-5)
    bin_max = hist_corr.GetXaxis().FindBin(mjj_max - 1.e-5)

    # Fit background*eff
    if "bbl" in analysis:
        background_fit_times_eff = ROOT.TF1("tmp_background_fit_times_eff", BackgroundFit_f4_trigcorr_bbl, mjj_min, mjj_max, 4)
        background_fit_times_eff.SetParameter(1, 41.)
        background_fit_times_eff.SetParameter(2, -45.)
        background_fit_times_eff.SetParameter(3, 10.)
    elif "bbh" in analysis:
        background_fit_times_eff = ROOT.TF1("tmp_background_fit_times_eff", BackgroundFit_f4_trigcorr_bbh, mjj_min, mjj_max, 4)
        background_fit_times_eff.SetParameter(1, 35.)
        background_fit_times_eff.SetParameter(2, -28.)
        background_fit_times_eff.SetParameter(3, 10.)
    background_fit_times_eff.SetParameter(0, hist_corr.Integral(bin_min, bin_max))
    background_fit_times_eff.SetParameter(0, hist_corr.Integral(bin_min, bin_max) / background_fit_times_eff.Integral(mjj_min, mjj_max))
    hist_corr.Fit(background_fit_times_eff, "QR0")

    # Make background TF1
    background_fit = ROOT.TF1("tmp_background_fit", BackgroundFit_f4, mjj_min, mjj_max, 4)
    for i in xrange(4):
        background_fit.SetParameter(i, background_fit_times_eff.GetParameter(i))

    # Correct histogram bins
    for bin in xrange(1, hist_corr.GetNbinsX() + 1):
        if bin < bin_min or bin > bin_max:
            hist_corr.SetBinContent(bin, 0)
            hist_corr.SetBinError(bin, 0)
        else:
            old_content = hist_corr.GetBinContent(bin)
            old_error = hist_corr.GetBinError(bin)
            this_bin_min = hist_corr.GetXaxis().GetBinLowEdge(bin)
            this_bin_max = hist_corr.GetXaxis().GetBinUpEdge(bin)
            num = background_fit.Integral(this_bin_min, this_bin_max)
            den = background_fit_times_eff.Integral(this_bin_min, this_bin_max)
            if ("bbl" in analysis and this_bin_min <= 350.) or ("bbh" in analysis and this_bin_min <= 575.):
                print "mjj={}: observed={}, corr={}".format(0.5*(this_bin_min+this_bin_max), den, num)
            if den > 0:
                correction = num / den
                #if ("bbl" in analysis and this_bin_min <= 350.) or ("bbh" in analysis and this_bin_min <= 575.):
                #    print "Trigger correction for mjj={}: {}".format(0.5*(this_bin_min+this_bin_max), correction)
            else:
                correction = 1.
            hist_corr.SetBinContent(bin, old_content * correction)
            hist_corr.SetBinError(bin, old_error * correction)
    return hist_corr

if __name__ == '__main__':
    main()

