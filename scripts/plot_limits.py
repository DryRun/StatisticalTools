#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from array import array
import numpy as np
import CMS_lumi
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/CMSDIJET/QCDAnalysis")
import analysis_configuration_8TeV as analysis_config


def main():
    # usage description
    usage = "Example: ./scripts/plotLimits.py -M Asymptotic -l logs -f qq --massrange 1200 7000 100"

    # input parameters
    parser = ArgumentParser(description='Script that plots limits for specified mass points',epilog=usage)
    parser.add_argument('analysis', type=str, help='Analysis name')
    parser.add_argument('model', type=str, help='Model name')

    parser.add_argument("-M", "--method", dest="method", required=True,
                        choices=['ProfileLikelihood', 'HybridNew', 'Asymptotic', 'MarkovChainMC', 'theta', 'HybridNewGrid'],
                        help="Method to calculate upper limits",
                        metavar="METHOD")
    parser.add_argument('--fit_function', type=str, default="f4", help="Name of fit function used for background estimate")
    parser.add_argument('--timesAE', action='store_true', help="Set y-axis to sigma*BR*A*e, instead of sigma*BR")
    #results_group = parser.add_mutually_exclusive_group(required=True)
    #results_group.add_argument("-l", "--logs_path", dest="logs_path",
    #                           help="Path to log files",
    #                           metavar="LOGS_PATH")
    #results_group.add_argument("-r", "--results_file", dest="results_file",
    #                           help="Path to a file containing results",
    #                           metavar="RESULTS_FILE")

    #parser.add_argument("-f", "--final_state", dest="final_state", required=True,
    #                    help="Final state (e.g. qq, qg, gg)",
    #                    metavar="FINAL_STATE")

    #parser.add_argument("-f2", "--finalstate2", dest="final_state2", required=True, help="hG,lG,hR, or lR", metavar="FINAL_STATE2")
    parser.add_argument("--noSyst", action="store_true", help="Make plots for limits without systematics")
    parser.add_argument("--freezeNuisances", type=str, help="Make plots for limits with frozen nuisance parameters")
    parser.add_argument("--postfix", dest="postfix", default='', help="Postfix for the output plot name (default: %(default)s)")

    parser.add_argument("--fileFormat", dest="fileFormat", default='pdf', help="Format of the output plot (default: %(default)s)")
    parser.add_argument("--saveObjects", type=str, help="Save plot objects")
    parser.add_argument("--extraText", dest="extraText", default='', help="Extra text on the plot (default: %(default)s)")

    parser.add_argument("--lumi_sqrtS", dest="lumi_sqrtS", default='19.7 fb^{-1} (8 TeV)', help="Integrated luminosity and center-of-mass energy (default: %(default)s)")

    parser.add_argument("--printResults", dest="printResults", default=False, action="store_true", help="Print results to the screen")

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

    if args.method == 'HybridNew':
        searchmethod = 'Hybrid New'

    # mass points for which resonance shapes will be produced
    input_masses = []

    if args.massrange != None:
        MIN, MAX, STEP = args.massrange
        input_masses = range(MIN, MAX+STEP, STEP)
    elif args.masslist != None:
        # A mass list was provided
        print "Will create mass list according to", args.masslist
        masslist = __import__(args.masslist.replace(".py",""))
        input_masses = masslist.masses
    else:
        input_masses = args.mass
    # sort masses
    input_masses.sort()

    from ROOT import kTRUE, kFALSE, gROOT, gStyle, gPad, TGraph, TCanvas, TLegend, TF1, TFile
    from ROOT import kGreen, kYellow, kWhite

    # Make acc*eff TGraph
    ae_x = array('d',[350, 400, 500, 600, 750, 900, 1200])
    if args.timesAE:
        ae_y = np.ones(len(ae_x))
    else:
        ae_y = array('d', [])
        for mass in ae_x:
            ae_y.append(analysis_config.simulation.get_signal_AE(args.analysis, args.model, int(mass)))
    acceptance_times_efficiency = TGraph(len(ae_x), ae_x, ae_y)
    #xs = array('d',[250,300,400,500,600,750,900,1200])

    #trigger_correctionl = TF1("trigbbl_efficiency", "(1. / (1. + TMath::Exp(-1. * (x - [0]) / [1])))**[2]", 175, 400)
    #trigger_correctionl.SetParameter(0, 1.82469e+02)
    #trigger_correctionl.SetParameter(1,  2.87768e+01)
    #trigger_correctionl.SetParameter(2,  9.11659e-01)

    #trigger_correctionh = TF1("trigbbh_efficiency", "(1. / (1. + TMath::Exp(-1. * (x - [0]) / [1])))**[2]", 300, 600)
    #trigger_correctionh.SetParameter(0, 3.61785e+02)
    #trigger_correctionh.SetParameter(1,  3.16523e+01)
    #trigger_correctionh.SetParameter(2,  4.84357e-01)
    #if args.timesAE:
    #    ys = np.ones(len(xs))
    #else:
    #    if args.analysis == "trigbbh_CSVTM" and args.model == "Hbb":
    #    	#ys = array('d',[188./19751.,1304./19993.,2697./49494.,881./19999.,534./19598.])
    #		ys = array('d',[24./19737.,188./19751.,1171./19984.,1419./19992.,1304./19993.,2697./49494.,881./19999.,534./19598.])
    #		#graphMod = trigger_correctionh
    #    elif args.analysis == "trigbbl_CSVTM" and args.model == "Hbb":
    #		#ys = array('d',[30./2797.,1583./19995.,1295./19996.,999./19996.,528./19999.])
    #		ys = array('d',[574./19737.,763./39502.,651./19984.,583./19992.,984./39986.,1905./98988.,656./39998.,369./39196.])
    #		#graphMod = trigger_correctionh
    #    elif args.analysis == "trigbbh_CSVTM" and args.model == "RSG":
    #		#ys = array('d',[109./19751.,488./19993.,954./49494.,328./19999.,182./19598.])
    #		ys = array('d',[40./19977.,30./2797.,1522./19991.,1640./19396.,1583./19995.,1295./19996.,999./19996.,528./19999.])
    #		#graphMod = trigger_correctionl
    #    elif args.analysis == "trigbbl_CSVTM" and args.model == "RSG":
    #		#ys = array('d',[23./2797.,599./19995.,448./19996.,338./19996.,190./19999.])
    #		ys = array('d',[696./19977.,137./5594.,797./19991.,652./19396.,1206./39990.,891./39992.,675./39992.,379./39998.])
    #		#graphMod = trigger_correctionl

    ##ys = array('d',[1,1,1,1,1,1,1,1,])   

    #acceptance_times_efficiency = TGraph(len(xs),xs,ys)


    # arrays holding results
    masses = array('d')
    xs_obs_limits = array('d')
    xs_exp_limits = array('d')
    masses_exp = array('d')
    xs_exp_limits_1sigma = array('d')
    xs_exp_limits_1sigma_up = array('d')
    xs_exp_limits_2sigma = array('d')
    xs_exp_limits_2sigma_up = array('d')

    for mass in input_masses:
        print ">> Reading results for %s %s resonance with m = %i GeV..."%(args.analysis, args.model, int(mass))
        masses.append(mass)
        masses_exp.append(mass)

        # For masses above 1100, you scaled down the signal by 10 by hand, to help the limit setting.
        #if args.analysis == "trigbbh_CSVTM" and mass >= 1100:
        input_xs = 1./100.
        #else:
        #    input_xs = 1.

        if args.method == "HybridNewGrid":
            found_limit = {"obs":False, "exp0":False, "exp1":False, "exp2":False, "exp-1":False, "exp-2":False}
            for what in found_limit.keys():
                log_file_path = limit_config.get_combine_log_path_grid(args.analysis, args.model, mass, args.fit_function, what, systematics=(not args.noSyst), frozen_nps=args.freezeNuisances)
                print "Reading log file from " + log_file_path
                log_file = open(log_file_path, 'r')
                for line in log_file:
                    if re.search("^Limit: r <", line) and re.search("95%", line):
                        found_limit[what] = True
                        this_limit = float(line.split()[3])/acceptance_times_efficiency.Eval(mass)
                        print "Found limit for " + what + " = " + str(this_limit)
                        if what == "obs":
                            xs_obs_limits.append(this_limit * input_xs)
                        elif what == "exp0":
                            xs_exp_limits.append(this_limit * input_xs)
                        elif what == "exp1":
                            xs_exp_limits_1sigma_up.append(this_limit * input_xs)
                        elif what == "exp2":
                            xs_exp_limits_2sigma_up.append(this_limit * input_xs)
                        elif what == "exp-1":
                            xs_exp_limits_1sigma.append(this_limit * input_xs)
                        elif what == "exp-2":
                            xs_exp_limits_2sigma.append(this_limit * input_xs)
            if not found_limit["obs"]:
                xs_obs_limits.append(0)
            if not found_limit["exp0"]:
                xs_exp_limits.append(0)
            if not found_limit["exp1"]:
                xs_exp_limits_1sigma.append(0)
            if not found_limit["exp2"]:
                xs_exp_limits_1sigma_up.append(0)
            if not found_limit["exp-1"]:
                xs_exp_limits_2sigma.append(0)
            if not found_limit["exp-2"]:
                xs_exp_limits_2sigma_up.append(0)
            if len(masses) != len(xs_obs_limits):
                print "** ERROR: ** Could not find observed limit for m =", int(mass), "GeV. Aborting."
                sys.exit(1)
        else:
            print "Reading log file from " + limit_config.get_combine_log_path(args.analysis, args.model, mass, args.fit_function, args.method, systematics=(not args.noSyst), frozen_nps=args.freezeNuisances)
            log_file = open(limit_config.get_combine_log_path(args.analysis, args.model, mass, args.fit_function, args.method, systematics=(not args.noSyst), frozen_nps=args.freezeNuisances))

            foundMethod = False
            middle = 0
            # read the log file
            found_limit = {"obs":False, "exp":False, "exp+1":False, "exp+2":False, "exp-1":False, "exp-2":False}
            for line in log_file:
                if args.method == 'Asymptotic':
                    if re.search("^Observed Limit: r", line):
                        xs_obs_limits.append(float(line.split()[-1])/acceptance_times_efficiency.Eval(mass) * input_xs)
                        found_limit["obs"] = True
                    if re.search("^Expected 50.0%: r", line):
                        middle = float(line.split()[-1])
                        found_limit["exp"] = True
                        xs_exp_limits.append(middle/acceptance_times_efficiency.Eval(mass) * input_xs)
                    if re.search("^Expected 16.0%: r", line):
                        xs_exp_limits_1sigma.append((float(line.split()[-1]))/acceptance_times_efficiency.Eval(mass) * input_xs)
                        found_limit["exp-1"] = True
                    if re.search("^Expected 84.0%: r", line):
                        xs_exp_limits_1sigma_up.append(float(line.split()[-1])/acceptance_times_efficiency.Eval(mass) * input_xs)
                        found_limit["exp+1"] = True
                    if re.search("^Expected  2.5%: r", line):
                        xs_exp_limits_2sigma.append(float(line.split()[-1])/acceptance_times_efficiency.Eval(mass) * input_xs)
                        found_limit["exp-2"] = True
                    if re.search("^Expected 97.5%: r", line):
                        xs_exp_limits_2sigma_up.append(float(line.split()[-1])/acceptance_times_efficiency.Eval(mass) * input_xs)
                        found_limit["exp+2"] = True
                elif args.method == 'theta':
                    if re.search('^# x; y; yerror', line):
                        foundMethod = True
                    if line.split()[0] == '0' and foundMethod:
                        xs_obs_limits.append(float(line.split()[1])/acceptance_times_efficiency.Eval(mass) * input_xs)
                else:
                    searchmethod = "Hybrid New"
                    if re.search(' -- ' + searchmethod, line):
                        foundMethod = True
                    if re.search("^Limit: r", line) and foundMethod:
                        xs_obs_limits.append(float(line.split()[3])/acceptance_times_efficiency.Eval(mass) * input_xs)
                        found_limit["obs"] = True
                        print "[debug] Found limit " + str(xs_obs_limits[-1])

            if not found_limit["obs"]:
                xs_obs_limits.append(0)
            if not found_limit["exp"]:
                xs_exp_limits.append(0)
            if not found_limit["exp+1"]:
                xs_exp_limits_1sigma.append(0)
            if not found_limit["exp+2"]:
                xs_exp_limits_1sigma_up.append(0)
            if not found_limit["exp-1"]:
                xs_exp_limits_2sigma.append(0)
            if not found_limit["exp-2"]:
                xs_exp_limits_2sigma_up.append(0)
            if len(masses) != len(xs_obs_limits):
                print "** ERROR: ** Could not find observed limit for m =", int(mass), "GeV. Aborting."
                sys.exit(1)

        if args.method == 'Asymptotic' or args.method == 'HybridNewGrid':
            if len(masses) != len(xs_exp_limits):
                print "** ERROR: ** Could not find expected limit for m =", int(mass), "GeV. Aborting."
                print "masses = ",
                print masses
                print "xs_exp_limits = ",
                print xs_exp_limits
                sys.exit(1)

            if len(masses) != len(xs_exp_limits_1sigma):
                print "** ERROR: ** Could not find expected 1 sigma down limit for m =", int(mass), "GeV. Aborting."
                print "masses = ",
                print masses
                print "xs_exp_limits_1sigma = ",
                print xs_exp_limits_1sigma
                sys.exit(1)

            if len(masses) != len(xs_exp_limits_1sigma_up):
                print "** ERROR: ** Could not find expected 1 sigma up limit for m =", int(mass), "GeV. Aborting."
                print "masses = ",
                print masses
                print "xs_exp_limits_1sigma_up = ",
                print xs_exp_limits_1sigma_up
                sys.exit(1)

            if len(masses) != len(xs_exp_limits_2sigma):
                print "** ERROR: ** Could not find expected 2 sigma down limit for m =", int(mass), "GeV. Aborting."
                sys.exit(1)

            if len(masses) != len(xs_exp_limits_2sigma_up):
                print "** ERROR: ** Could not find expected 2 sigma up limit for m =", int(mass), "GeV. Aborting."
                sys.exit(1)
    if args.method == 'Asymptotic' or args.method == 'HybridNewGrid':
        # complete the expected limit arrays
        for i in range(0,len(masses)):
            masses_exp.append( masses[len(masses)-i-1] )
            xs_exp_limits_1sigma.append( xs_exp_limits_1sigma_up[len(masses)-i-1] )
            xs_exp_limits_2sigma.append( xs_exp_limits_2sigma_up[len(masses)-i-1] )

    if args.printResults:
        print "masses =", masses.tolist()
        print "xs_obs_limits =", xs_obs_limits.tolist()
        print "xs_exp_limits =", xs_exp_limits.tolist()
        print ""
        print "masses_exp =", masses_exp.tolist()
        print "xs_exp_limits_1sigma =", xs_exp_limits_1sigma.tolist()
        print "xs_exp_limits_2sigma =", xs_exp_limits_2sigma.tolist()


    gROOT.SetBatch(kTRUE);
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    gStyle.SetTitleFont(42, "XYZ")
    gStyle.SetTitleSize(0.06, "XYZ")
    gStyle.SetLabelFont(42, "XYZ")
    gStyle.SetLabelSize(0.05, "XYZ")
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetFrameBorderMode(0)
    gStyle.SetCanvasColor(kWhite)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetPadBottomMargin(0.15)
    gROOT.ForceStyle()

    # theory curves: gg
    massesS8 = array('d', [1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0,2100.0,2200.0,2300.0,2400.0,2500.0,2600.0,2700.0,2800.0,2900.0,3000.0,3100.0,3200.0,3300.0,3400.0,3500.0,3600.0,3700.0,3800.0,3900.0,4000.0,4100.0,4200.0,4300.0,4400.0,4500.0,4600.0,4700.0,4800.0,4900.0,5000.0,5100.0,5200.0,5300.0,5400.0,5500.0,5600.0,5700.0,5800.0,5900.0,6000.0])
    xsS8 = array('d', [5.46E+02,3.12E+02,1.85E+02,1.12E+02,7.19E+01,4.59E+01,3.02E+01,2.01E+01,1.37E+01,9.46E+00,6.55E+00,4.64E+00,3.27E+00,2.36E+00,1.70E+00,1.24E+00,9.11E-01,6.69E-01,4.97E-01,3.71E-01,2.78E-01,2.07E-01,1.55E-01,1.19E-01,9.26E-02,7.08E-02,5.43E-02,4.15E-02,3.22E-02,2.50E-02,1.92E-02,1.51E-02,1.19E-02,9.25E-03,7.35E-03,5.86E-03,4.53E-03,3.66E-03,2.91E-03,2.33E-03,1.86E-03,1.45E-03,1.12E-03,8.75E-04,6.90E-04,5.55E-04,4.47E-04,3.63E-04,2.92E-04,2.37E-04,1.97E-04])

    graph_xsS8 = TGraph(len(massesS8),massesS8,xsS8)
    graph_xsS8.SetLineWidth(3)
    graph_xsS8.SetLineStyle(8)
    graph_xsS8.SetLineColor(6)

    # theory curves: qg
    massesString = array('d', [1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0,2100.0,2200.0,2300.0,2400.0,2500.0,2600.0,2700.0,2800.0,2900.0,3000.0,3100.0,3200.0,3300.0,3400.0,3500.0,3600.0,3700.0,3800.0,3900.0,4000.0,4100.0,4200.0,4300.0,4400.0,4500.0,4600.0,4700.0,4800.0,4900.0,5000.0,5100.0,5200.0,5300.0,5400.0,5500.0,5600.0,5700.0,5800.0,5900.0,6000.0,6100.0,6200.0,6300.0,6400.0,6500.0,6600.0,6700.0,6800.0,6900.0,7000.0,7100.0,7200.0,7300.0,7400.0,7500.0,7600.0,7700.0,7800.0,7900.0,8000.0,8100.0,8200.0,8300.0,8400.0,8500.0,8600.0,8700.0,8800.0,8900.0,9000.0,9100.,9200.,9300.,9400.,9500.,9600.,9700.,9800.,9900.,10000.])
    xsString = array('d', [8316.184311558545,5312.93137758767,3435.0309937336524,2304.4139502741305,1569.8115447896687,1090.9516635659693,770.901859690924,551.9206062572061,399.69535383507633,293.77957451762086,218.15126842827823,162.87634729465125,123.17685479653694,93.63530805932386,71.53697229809124,55.37491301647483,42.75271508357369,33.36378355470234,26.06619302090876,20.311817606835643,16.1180931789545,12.768644973921226,10.142660425967444,8.057990848043234,6.400465846290908,5.115134438331436,4.132099789492928,3.3193854239538734,2.6581204529344302,2.157554604919995,1.7505176068913348,1.4049155245498584,1.140055677916783,0.9253251132104159,0.7522038169131606,0.6119747371392215,0.49612321727328523,0.40492020959456737,0.33091999402250655,0.27017917021492555,0.2201693919322846,0.17830700070267996,0.14564253802358157,0.11940534430331146,0.09694948234356839,0.0793065371847468,0.06446186373361917,0.05282660618352478,
                           0.0428516302310620888,0.0348997638039910363,0.0283334766442618227,0.0231416918363592127,0.0187417921340763783,0.0153501307395115115,0.0124396534127133717,0.0100542205744949455,0.0081744954858627415,0.0066338099362915941,0.0053365711503318145,0.00430912459914657443,0.00346381039244064343,0.00278602671711227174,0.00225154342228859257,0.0018082930150063248,0.00143929440338502119,0.0011581373956044489,0.00091869589873893118,0.00073410823691329855,0.00058669382997948734,0.0004661568745858897,0.000368716655469570365,0.000293168485206959169,0.000230224535021638668,0.000182317101888465142,0.000143263359883433282,0.000112630538527214965,0.000088189175598406759,0.000068708474367442343,0.000053931726669273556,0.0000416417855733682702,0.0000326529676755488658,0.0000254365480426201587,0.0000198410151166864761,0.0000154034425617473576,0.0000119095554601641413,9.2537574320108232e-6,7.2155417437856749e-6,5.6130924422251982e-6,4.36634755605624901e-6,3.39717456406994868e-6,2.6766018046173896e-6])

    massesQstar = array('d', [1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0,2100.0,2200.0,2300.0,2400.0,2500.0,2600.0,2700.0,2800.0,2900.0,3000.0,3100.0,3200.0,3300.0,3400.0,3500.0,3600.0,3700.0,3800.0,3900.0,4000.0,4100.0,4200.0,4300.0,4400.0,4500.0,4600.0,4700.0,4800.0,4900.0,5000.0,5100.0,5200.0,5300.0,5400.0,5500.0,5600.0,5700.0,5800.0,5900.0,6000.0,6100.0,6200.0,6300.0,6400.0,6500.0,6600.0,6700.0,6800.0,6900.0,7000.0,7100.0,7200.0,7300.0,7400.0,7500.0,7600.0,7700.0,7800.0,7900.0,8000.0,8100.0,8200.0,8300.0,8400.0,8500.0,8600.0,8700.0,8800.0,8900.0,9000.0])
    xsQstar = array('d', [0.4101E+03,0.2620E+03,0.1721E+03,0.1157E+03,0.7934E+02,0.5540E+02,0.3928E+02,0.2823E+02,0.2054E+02,0.1510E+02,0.1121E+02,0.8390E+01,0.6328E+01,0.4807E+01,0.3674E+01,0.2824E+01,0.2182E+01,0.1694E+01,0.1320E+01,0.1033E+01,0.8116E+00,0.6395E+00,0.5054E+00,0.4006E+00,0.3182E+00,0.2534E+00,0.2022E+00,0.1616E+00,0.1294E+00,0.1038E+00,0.8333E-01,0.6700E-01,0.5392E-01,0.4344E-01,0.3503E-01,0.2827E-01,0.2283E-01,0.1844E-01,0.1490E-01,0.1205E-01,0.9743E-02,0.7880E-02,0.6373E-02,0.5155E-02,0.4169E-02,0.3371E-02,0.2725E-02,0.2202E-02,0.1779E-02,0.1437E-02,0.1159E-02,0.9353E-03,0.7541E-03,0.6076E-03,0.4891E-03,0.3935E-03,0.3164E-03,0.2541E-03,0.2039E-03,0.1635E-03,0.1310E-03,0.1049E-03,0.8385E-04,0.6699E-04,0.5347E-04,0.4264E-04,0.3397E-04,0.2704E-04,0.2151E-04,0.1709E-04,0.1357E-04,0.1077E-04,0.8544E-05,0.6773E-05,0.5367E-05,0.4251E-05,0.3367E-05,0.2666E-05,0.2112E-05,0.1673E-05,0.1326E-05])

    graph_xsString = TGraph(len(massesString),massesString,xsString)
    graph_xsString.SetLineWidth(3)
    graph_xsString.SetLineStyle(8)
    graph_xsString.SetLineColor(9)

    graph_xsQstar = TGraph(len(massesQstar),massesQstar,xsQstar)
    graph_xsQstar.SetLineWidth(3)
    graph_xsQstar.SetLineStyle(2)
    graph_xsQstar.SetLineColor(1)

    # theory curves: qq
    massesTh = array('d', [1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0,2100.0,2200.0,2300.0,2400.0,2500.0,2600.0,2700.0,2800.0,2900.0,3000.0,3100.0,3200.0,3300.0,3400.0,3500.0,3600.0,3700.0,3800.0,3900.0,4000.0,4100.0,4200.0,4300.0,4400.0,4500.0,4600.0,4700.0,4800.0,4900.0,5000.0,5100.0,5200.0,5300.0,5400.0,5500.0,5600.0,5700.0,5800.0,5900.0,6000.0,6100.0,6200.0,6300.0,6400.0,6500.0,6600.0,6700.0,6800.0,6900.0,7000.0,7100.0,7200.0,7300.0,7400.0,7500.0,7600.0,7700.0,7800.0,7900.0,8000.0,8100.0,8200.0,8300.0,8400.0,8500.0,8600.0,8700.0,8800.0,8900.0,9000.0])

    xsAxi = array('d', [0.1849E+03,0.1236E+03,0.8473E+02,0.5937E+02,0.4235E+02,0.3069E+02,0.2257E+02,0.1680E+02,0.1263E+02,0.9577E+01,0.7317E+01,0.5641E+01,0.4374E+01,0.3411E+01,0.2672E+01,0.2103E+01,0.1658E+01,0.1312E+01,0.1041E+01,0.8284E+00,0.6610E+00,0.5294E+00,0.4250E+00,0.3417E+00,0.2752E+00,0.2220E+00,0.1792E+00,0.1449E+00,0.1172E+00,0.9487E-01,0.7686E-01,0.6219E-01,0.5033E-01,0.4074E-01,0.3298E-01,0.2671E-01,0.2165E-01,0.1755E-01,0.1422E-01,0.1152E-01,0.9322E-02,0.7539E-02,0.6092E-02,0.4917E-02,0.3965E-02,0.3193E-02,0.2568E-02,0.2062E-02,0.1653E-02,0.1323E-02,0.1057E-02,0.8442E-03,0.6728E-03,0.5349E-03,0.4242E-03,0.3357E-03,0.2644E-03,0.2077E-03,0.1627E-03,0.1271E-03,0.9891E-04,0.7686E-04,0.5951E-04,0.4592E-04,0.3530E-04,0.2704E-04,0.2059E-04,0.1562E-04,0.1180E-04,0.8882E-05,0.6657E-05,0.4968E-05,0.3693E-05,0.2734E-05,0.2016E-05,0.1481E-05,0.1084E-05,0.7903E-06,0.5744E-06,0.4160E-06,0.3007E-06])
    xsDiquark = array('d', [0.5824E+02,0.4250E+02,0.3172E+02,0.2411E+02,0.1862E+02,0.1457E+02,0.1153E+02,0.9211E+01,0.7419E+01,0.6019E+01,0.4912E+01,0.4031E+01,0.3323E+01,0.2750E+01,0.2284E+01,0.1903E+01,0.1590E+01,0.1331E+01,0.1117E+01,0.9386E+00,0.7900E+00,0.6658E+00,0.5618E+00,0.4745E+00,0.4010E+00,0.3391E+00,0.2869E+00,0.2428E+00,0.2055E+00,0.1740E+00,0.1473E+00,0.1246E+00,0.1055E+00,0.8922E-01,0.7544E-01,0.6376E-01,0.5385E-01,0.4546E-01,0.3834E-01,0.3231E-01,0.2720E-01,0.2288E-01,0.1922E-01,0.1613E-01,0.1352E-01,0.1132E-01,0.9463E-02,0.7900E-02,0.6584E-02,0.5479E-02,0.4551E-02,0.3774E-02,0.3124E-02,0.2581E-02,0.2128E-02,0.1750E-02,0.1437E-02,0.1177E-02,0.9612E-03,0.7833E-03,0.6366E-03,0.5160E-03,0.4170E-03,0.3360E-03,0.2700E-03,0.2162E-03,0.1725E-03,0.1372E-03,0.1087E-03,0.8577E-04,0.6742E-04,0.5278E-04,0.4114E-04,0.3192E-04,0.2465E-04,0.1894E-04,0.1448E-04,0.1101E-04,0.8322E-05,0.6253E-05,0.4670E-05])
    xsWprime = array('d', [0.8811E+01,0.6024E+01,0.4216E+01,0.3010E+01,0.2185E+01,0.1610E+01,0.1200E+01,0.9043E+00,0.6875E+00,0.5271E+00,0.4067E+00,0.3158E+00,0.2464E+00,0.1932E+00,0.1521E+00,0.1201E+00,0.9512E-01,0.7554E-01,0.6012E-01,0.4792E-01,0.3827E-01,0.3059E-01,0.2448E-01,0.1960E-01,0.1571E-01,0.1259E-01,0.1009E-01,0.8090E-02,0.6483E-02,0.5193E-02,0.4158E-02,0.3327E-02,0.2660E-02,0.2125E-02,0.1695E-02,0.1351E-02,0.1075E-02,0.8546E-03,0.6781E-03,0.5372E-03,0.4248E-03,0.3353E-03,0.2642E-03,0.2077E-03,0.1629E-03,0.1275E-03,0.9957E-04,0.7757E-04,0.6027E-04,0.4670E-04,0.3610E-04,0.2783E-04,0.2140E-04,0.1641E-04,0.1254E-04,0.9561E-05,0.7269E-05,0.5510E-05,0.4167E-05,0.3143E-05,0.2364E-05,0.1774E-05,0.1329E-05,0.9931E-06,0.7411E-06,0.5523E-06,0.4108E-06,0.3055E-06,0.2271E-06,0.1687E-06,0.1254E-06,0.9327E-07,0.6945E-07,0.5177E-07,0.3863E-07,0.2888E-07,0.2162E-07,0.1622E-07,0.1218E-07,0.9156E-08,0.6893E-08])
    xsZprime = array('d', [0.5027E+01,0.3398E+01,0.2353E+01,0.1663E+01,0.1196E+01,0.8729E+00,0.6450E+00,0.4822E+00,0.3638E+00,0.2769E+00,0.2123E+00,0.1639E+00,0.1272E+00,0.9933E-01,0.7789E-01,0.6134E-01,0.4848E-01,0.3845E-01,0.3059E-01,0.2440E-01,0.1952E-01,0.1564E-01,0.1256E-01,0.1010E-01,0.8142E-02,0.6570E-02,0.5307E-02,0.4292E-02,0.3473E-02,0.2813E-02,0.2280E-02,0.1848E-02,0.1499E-02,0.1216E-02,0.9864E-03,0.8002E-03,0.6490E-03,0.5262E-03,0.4264E-03,0.3453E-03,0.2795E-03,0.2260E-03,0.1826E-03,0.1474E-03,0.1188E-03,0.9566E-04,0.7690E-04,0.6173E-04,0.4947E-04,0.3957E-04,0.3159E-04,0.2516E-04,0.2001E-04,0.1587E-04,0.1255E-04,0.9906E-05,0.7795E-05,0.6116E-05,0.4785E-05,0.3731E-05,0.2900E-05,0.2247E-05,0.1734E-05,0.1334E-05,0.1022E-05,0.7804E-06,0.5932E-06,0.4492E-06,0.3388E-06,0.2544E-06,0.1903E-06,0.1417E-06,0.1051E-06,0.7764E-07,0.5711E-07,0.4186E-07,0.3055E-07,0.2223E-07,0.1612E-07,0.1164E-07,0.8394E-08])

    graph_xsAxi = TGraph(len(massesTh),massesTh,xsAxi)
    graph_xsAxi.SetLineWidth(3)
    graph_xsAxi.SetLineStyle(3)
    graph_xsAxi.SetLineColor(63)

    graph_xsDiquark = TGraph(len(massesTh),massesTh,xsDiquark)
    graph_xsDiquark.SetLineWidth(3)
    graph_xsDiquark.SetLineStyle(9)
    graph_xsDiquark.SetLineColor(8)

    graph_xsWprime = TGraph(len(massesTh),massesTh,xsWprime)
    graph_xsWprime.SetLineWidth(3)
    graph_xsWprime.SetLineStyle(7)
    graph_xsWprime.SetLineColor(46)

    graph_xsZprime = TGraph(len(massesTh),massesTh,xsZprime)
    graph_xsZprime.SetLineWidth(3)
    graph_xsZprime.SetLineStyle(5)
    graph_xsZprime.SetLineColor(38)

    # limits
    if args.method == "Asymptotic" or args.method == "HybridNewGrid":
        graph_exp_2sigma = ( TGraph(len(masses_exp),masses_exp,xs_exp_limits_2sigma) if len(xs_exp_limits_2sigma) > 0 else TGraph(0) )
        graph_exp_2sigma.SetFillColor(kYellow)

        graph_exp_1sigma = ( TGraph(len(masses_exp),masses_exp,xs_exp_limits_1sigma) if len(xs_exp_limits_2sigma) > 0 else TGraph(0) )
        graph_exp_1sigma.SetFillColor(kGreen+1)

        graph_exp = ( TGraph(len(masses),masses,xs_exp_limits) if len(xs_exp_limits_2sigma) > 0 else TGraph(0) )
        #graph_exp.SetMarkerStyle(24)
        graph_exp.SetLineWidth(3)
        graph_exp.SetLineStyle(2)
        graph_exp.SetLineColor(4)

    graph_obs = TGraph(len(masses),masses,xs_obs_limits)
    graph_obs.SetMarkerStyle(20)
    graph_obs.SetLineWidth(3)
    #graph_obs.SetLineStyle(1)
    graph_obs.SetLineColor(1)

    c = TCanvas("c", "",800,800)
    c.cd()

    legend = TLegend(.60,.75,.90,.90)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)
    legend.SetHeader('95% CL upper limits')

    if len(xs_exp_limits_2sigma) > 0 and (args.method == "Asymptotic" or args.method == "HybridNewGrid"):
        frame = graph_exp_2sigma.GetHistogram().Clone()
    else:
        frame = graph_obs.GetHistogram().Clone()
    frame.Reset()
    frame.GetXaxis().SetTitle("Resonance mass [GeV]")
    if args.timesAE:
        if args.model == "ZPrime":
            frame.GetYaxis().SetTitle("#sigma #times (BR(c#bar{c})+BR(b#bar{b})) #times #it{A} #times #epsilon [pb]")
        else:
            frame.GetYaxis().SetTitle("#sigma #times BR(b#bar{b}) #times #it{A} #times #epsilon [pb]")
    else:
        if args.model == "ZPrime":
            frame.GetYaxis().SetTitle("#sigma #times (BR(c#bar{c})+BR(b#bar{b})) [pb]")
        else:
            frame.GetYaxis().SetTitle("#sigma #times BR(b#bar{b}) [pb]")
    frame.GetYaxis().SetTitleOffset(1.1)
    if args.timesAE:
        frame.GetYaxis().SetRangeUser(1e-03,1e+01)
    else:
        frame.GetYaxis().SetRangeUser(1e-02,1e+02)
    frame.Draw("axis")

    if len(xs_exp_limits_2sigma) > 0 and (args.method == "Asymptotic" or args.method == "HybridNewGrid"):
        graph_exp_2sigma.GetXaxis().SetTitle("Resonance mass [GeV]")
        graph_exp_2sigma.GetYaxis().SetTitle("#sigma #times #it{B} [pb]")
        graph_exp_2sigma.GetYaxis().SetTitleOffset(1.1)
        graph_exp_2sigma.GetYaxis().SetRangeUser(1e-03,1e+02)
        #graph_exp_2sigma.GetXaxis().SetNdivisions(1005)

        graph_exp_2sigma.Draw("F")
        graph_exp_1sigma.Draw("F")
        graph_exp.Draw("L")
        graph_obs.Draw("LP")

        legend.AddEntry(graph_obs,"Observed","lp")
        legend.AddEntry(graph_exp,"Expected","lp")
        legend.AddEntry(graph_exp_1sigma,"#pm 1#sigma","F")
        legend.AddEntry(graph_exp_2sigma,"#pm 2#sigma","F")
    else:
        graph_obs.GetXaxis().SetTitle("Resonance mass [GeV]")
        graph_obs.GetYaxis().SetTitle("#sigma #times #it{B} [pb]")
        graph_obs.GetYaxis().SetTitleOffset(1.1)
        graph_obs.GetYaxis().SetRangeUser(1e-02,1e+03)
        #graph_obs.GetXaxis().SetNdivisions(1005)

        graph_obs.Draw("LP")

        legend.AddEntry(graph_obs,"Observed","lp")

        #if args.final_state == 'gg' :
        #    graph_xsS8.Draw("L")
        #elif args.final_state == 'qg' :
        #    graph_xsQstar.Draw("L")
        #    graph_xsString.Draw("L")
        #elif args.final_state == 'qq' :
        #    graph_xsAxi.Draw("L")
        #    graph_xsDiquark.Draw("L")
        #    graph_xsWprime.Draw("L")
        #    graph_xsZprime.Draw("L")
        
    legend.Draw()

    #legendTh = TLegend(.60,.72,.90,.88)
    #legendTh.SetBorderSize(0)
    #legendTh.SetFillColor(0)
    #legendTh.SetFillStyle(0)
    #legendTh.SetTextFont(42)
    #legendTh.SetTextSize(0.03)
    #legendTh.AddEntry(graph_xsAxi,"Axigluon/coloron","l")
    #legendTh.AddEntry(graph_xsDiquark,"Scalar diquark","l")
    #legendTh.AddEntry(graph_xsWprime,"W' SSM","l")
    #legendTh.AddEntry(graph_xsZprime,"Z' SSM","l")
    #legendTh.Draw()

    # draw the lumi text on the canvas
    CMS_lumi.extraText = args.extraText
    CMS_lumi.lumi_sqrtS = args.lumi_sqrtS # used with iPeriod = 0 (free form)
    iPos = 11
    iPeriod = 0

    CMS_lumi.CMS_lumi(c, iPeriod, iPos)

    gPad.RedrawAxis()

    c.SetLogy()
    postfix = ( ('_' + args.postfix) if args.postfix != '' else '' )
    if args.noSyst:
        postfix += "_noSyst"
    if args.freezeNuisances:
        postfix += "_" + args.freezeNuisances.replace(",", "_")
    fileName = limit_config.paths["limit_plots"] + '/xs_limit_%s_%s_%s_%s.%s'%(args.method,args.analysis, args.model + postfix, args.fit_function, args.fileFormat.lower())
    if args.timesAE:
        fileName = fileName.replace("xs_limit", "xsAE_limit")
    c.SaveAs(fileName)
    print "Plot saved to '%s'"%(fileName)

    if args.saveObjects:
        output_file = args.saveObjects
        if args.timesAE:
            output_file = output_file.replace(".root", "_timesAE.root")
        f = TFile(output_file, "RECREATE")
        if args.method == "Asymptotic" or args.method == "HybridNewGrid":
            graph_exp_2sigma.SetName("graph_exp_2sigma")
            graph_exp_2sigma.Write()
            graph_exp_1sigma.SetName("graph_exp_1sigma")
            graph_exp_1sigma.Write()
            graph_exp.SetName("graph_exp")
            graph_exp.Write()
        graph_obs.SetName("graph_obs")
        graph_obs.Write()
        f.Close()

if __name__ == '__main__':
    main()
