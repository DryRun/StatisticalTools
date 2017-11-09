#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from array import array
import numpy as np
import pickle
import CMS_lumi
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/CMSDIJET/QCDAnalysis")
import analysis_configuration_8TeV as analysis_config
import ROOT
from ROOT import *
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")

# Similar to plot_limits, except:
# - Plot both low and high mass SRs on one plot
# - Masses (and analyses) are hard-coded
# - Do all models in one go

model_names_pretty = {
    "Hbb":"Scalar",
    "ZPrime":"Z'",
    "RSG":"RS Graviton",
}

def main():
    # usage description
    usage = "Example: ./scripts/plotLimits.py -M Asymptotic -l logs -f qq --massrange 1200 7000 100"

    # input parameters
    parser = ArgumentParser(description='Script that plots limits for specified mass points',epilog=usage)
    parser.add_argument('--models', type=str, default="Hbb,ZPrime,RSG", help='Model name')

    parser.add_argument("-M", "--method", dest="method", required=True,
                        choices=['ProfileLikelihood', 'HybridNew', 'Asymptotic', 'MarkovChainMC', 'theta', 'HybridNewGrid'],
                        help="Method to calculate upper limits",
                        metavar="METHOD")
    parser.add_argument('--fit_function', type=str, default="dijet4", help="Name of fit function used for background estimate")
    parser.add_argument('--timesAE', action='store_true', help="Set y-axis to sigma*BR*A*e, instead of sigma*BR")
    parser.add_argument('--fitTrigger', action='store_true', help="Use trigger fit")
    parser.add_argument('--correctTrigger', action='store_true', help="Use trigger correction")
    parser.add_argument('--useMCTrigger', action='store_true', help="Use MC trigger emulation")
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
    args = parser.parse_args()

    if args.method == 'HybridNew':
        searchmethod = 'Hybrid New'

    analyses = ["trigbbl_CSVTM", "trigbbh_CSVTM"]
    analysis_masses = {"trigbbl_CSVTM":[325,350,400,450,500,550,600,650,700], "trigbbh_CSVTM":range(700,1250,50)}


    models = args.models.split(",")

    acceptance_times_efficiency = {}
    for model in models:
        # Make acc*eff TGraph
        ae_x = array('d',[325, 350, 400, 500, 600, 750, 900, 1200])
        ae_y = {}
        for analysis in analyses:
            if args.timesAE:
                ae_y[analysis] = np.ones(len(ae_x))
            else:
                ae_y[analysis] = array('d', [])
                for mass in ae_x:
                    ae_y[analysis].append(analysis_config.simulation.get_signal_AE(analysis, model, int(mass)))
            acceptance_times_efficiency[analysis] = TGraph(len(ae_x), ae_x, ae_y[analysis])

        # arrays holding results
        masses = {}
        xs_obs_limits = {}
        xs_exp_limits = {}
        masses_exp = {}
        xs_exp_limits_1sigma = {}
        xs_exp_limits_1sigma_up = {}
        xs_exp_limits_2sigma = {}
        xs_exp_limits_2sigma_up = {}

        for analysis in analyses:
            masses[analysis] = array('d')
            xs_obs_limits[analysis] = array('d')
            xs_exp_limits[analysis] = array('d')
            masses_exp[analysis] = array('d')
            xs_exp_limits_1sigma[analysis] = array('d')
            xs_exp_limits_1sigma_up[analysis] = array('d')
            xs_exp_limits_2sigma[analysis] = array('d')
            xs_exp_limits_2sigma_up[analysis] = array('d')

            for mass in analysis_masses[analysis]:
                print ">> Reading results for %s %s resonance with m = %i GeV..."%(analysis, model, int(mass))
                masses[analysis].append(mass)
                masses_exp[analysis].append(mass)

                # For masses above 1100, you scaled down the signal by 10 by hand, to help the limit setting.
                #if analysis == "trigbbh_CSVTM" and mass >= 1100:
                input_xs = 1./100.
                #else:
                #    input_xs = 1.

                if args.method == "HybridNewGrid":
                    found_limit = {"obs":False, "exp0":False, "exp1":False, "exp2":False, "exp-1":False, "exp-2":False}
                    for what in found_limit.keys():
                        log_file_path = limit_config.get_combine_log_path_grid(analysis, model, mass, args.fit_function, what, systematics=(not args.noSyst), frozen_nps=args.freezeNuisances, fitTrigger=args.fitTrigger, correctTrigger=args.correctTrigger, useMCTrigger=args.useMCTrigger)
                        print "Reading log file from " + log_file_path
                        log_file = open(log_file_path, 'r')
                        for line in log_file:
                            if re.search("^Limit: r <", line) and re.search("95%", line):
                                found_limit[what] = True
                                this_limit = float(line.split()[3])/acceptance_times_efficiency[analysis].Eval(mass)
                                print "Found limit for " + what + " = " + str(this_limit)
                                if what == "obs":
                                    xs_obs_limits[analysis].append(this_limit * input_xs)
                                elif what == "exp0":
                                    xs_exp_limits[analysis].append(this_limit * input_xs)
                                elif what == "exp1":
                                    xs_exp_limits_1sigma_up[analysis].append(this_limit * input_xs)
                                elif what == "exp2":
                                    xs_exp_limits_2sigma_up[analysis].append(this_limit * input_xs)
                                elif what == "exp-1":
                                    xs_exp_limits_1sigma[analysis].append(this_limit * input_xs)
                                elif what == "exp-2":
                                    xs_exp_limits_2sigma[analysis].append(this_limit * input_xs)
                    if not found_limit["obs"]:
                        xs_obs_limits[analysis].append(0)
                    if not found_limit["exp0"]:
                        xs_exp_limits[analysis].append(0)
                    if not found_limit["exp1"]:
                        xs_exp_limits_1sigma[analysis].append(0)
                    if not found_limit["exp2"]:
                        xs_exp_limits_1sigma_up[analysis].append(0)
                    if not found_limit["exp-1"]:
                        xs_exp_limits_2sigma[analysis].append(0)
                    if not found_limit["exp-2"]:
                        xs_exp_limits_2sigma_up[analysis].append(0)
                    if len(masses) != len(xs_obs_limits):
                        print "** ERROR: ** Could not find observed limit for m =", int(mass), "GeV. Aborting."
                        sys.exit(1)
                else:
                    print "Reading log file from " + limit_config.get_combine_log_path(analysis, model, mass, args.fit_function, args.method, systematics=(not args.noSyst), frozen_nps=args.freezeNuisances, fitTrigger=args.fitTrigger, correctTrigger=args.correctTrigger, useMCTrigger=args.useMCTrigger)
                    if not os.path.exists((limit_config.get_combine_log_path(analysis, model, mass, args.fit_function, args.method, systematics=(not args.noSyst), frozen_nps=args.freezeNuisances, fitTrigger=args.fitTrigger, correctTrigger=args.correctTrigger, useMCTrigger=args.useMCTrigger))):
                        print "[plot_limits] WARNING : Log file not found! Setting limits to zero and skipping this point."
                        print "[plot_limits] WARNING : \t{}".format(limit_config.get_combine_log_path(analysis, model, mass, args.fit_function, args.method, systematics=(not args.noSyst), frozen_nps=args.freezeNuisances, fitTrigger=args.fitTrigger, correctTrigger=args.correctTrigger, useMCTrigger=args.useMCTrigger))
                        xs_obs_limits[analysis].append(0)
                        xs_exp_limits[analysis].append(0)
                        xs_exp_limits_1sigma[analysis].append(0)
                        xs_exp_limits_1sigma_up[analysis].append(0)
                        xs_exp_limits_2sigma[analysis].append(0)
                        xs_exp_limits_2sigma_up[analysis].append(0)
                        continue
                    log_file = open(limit_config.get_combine_log_path(analysis, model, mass, args.fit_function, args.method, systematics=(not args.noSyst), frozen_nps=args.freezeNuisances, fitTrigger=args.fitTrigger, correctTrigger=args.correctTrigger, useMCTrigger=args.useMCTrigger))

                    foundMethod = False
                    middle = 0
                    # read the log file
                    found_limit = {"obs":False, "exp":False, "exp+1":False, "exp+2":False, "exp-1":False, "exp-2":False}
                    for line in log_file:
                        if args.method == 'Asymptotic':
                            if re.search("^Observed Limit: r", line):
                                xs_obs_limits[analysis].append(float(line.split()[-1])/acceptance_times_efficiency[analysis].Eval(mass) * input_xs)
                                found_limit["obs"] = True
                                if mass == 325 and model == "ZPrime":
                                    print "[debug] ZPrime 325 GeV limit = {}".format(xs_obs_limits[analysis][-1])
                                    print "[debug] \tA*e={}, input_xs={}".format(acceptance_times_efficiency[analysis].Eval(mass), input_xs)
                            if re.search("^Expected 50.0%: r", line):
                                middle = float(line.split()[-1])
                                found_limit["exp"] = True
                                xs_exp_limits[analysis].append(middle/acceptance_times_efficiency[analysis].Eval(mass) * input_xs)
                            if re.search("^Expected 16.0%: r", line):
                                xs_exp_limits_1sigma[analysis].append((float(line.split()[-1]))/acceptance_times_efficiency[analysis].Eval(mass) * input_xs)
                                found_limit["exp-1"] = True
                            if re.search("^Expected 84.0%: r", line):
                                xs_exp_limits_1sigma_up[analysis].append(float(line.split()[-1])/acceptance_times_efficiency[analysis].Eval(mass) * input_xs)
                                found_limit["exp+1"] = True
                            if re.search("^Expected  2.5%: r", line):
                                xs_exp_limits_2sigma[analysis].append(float(line.split()[-1])/acceptance_times_efficiency[analysis].Eval(mass) * input_xs)
                                found_limit["exp-2"] = True
                            if re.search("^Expected 97.5%: r", line):
                                xs_exp_limits_2sigma_up[analysis].append(float(line.split()[-1])/acceptance_times_efficiency[analysis].Eval(mass) * input_xs)
                                found_limit["exp+2"] = True
                        elif args.method == 'theta':
                            if re.search('^# x; y; yerror', line):
                                foundMethod = True
                            if line.split()[0] == '0' and foundMethod:
                                xs_obs_limits[analysis].append(float(line.split()[1])/acceptance_times_efficiency[analysis].Eval(mass) * input_xs)
                        else:
                            searchmethod = "Hybrid New"
                            if re.search(' -- ' + searchmethod, line):
                                foundMethod = True
                            if re.search("^Limit: r", line) and foundMethod:
                                xs_obs_limits[analysis].append(float(line.split()[3])/acceptance_times_efficiency[analysis].Eval(mass) * input_xs)
                                found_limit["obs"] = True
                                print "[debug] Found limit " + str(xs_obs_limits[analysis][-1])

                    if not found_limit["obs"]:
                        xs_obs_limits[analysis].append(0)
                    if not found_limit["exp"]:
                        xs_exp_limits[analysis].append(0)
                    if not found_limit["exp+1"]:
                        xs_exp_limits_1sigma[analysis].append(0)
                    if not found_limit["exp+2"]:
                        xs_exp_limits_1sigma_up[analysis].append(0)
                    if not found_limit["exp-1"]:
                        xs_exp_limits_2sigma[analysis].append(0)
                    if not found_limit["exp-2"]:
                        xs_exp_limits_2sigma_up[analysis].append(0)
                    if len(masses[analysis]) != len(xs_obs_limits[analysis]):
                        print "** ERROR: ** Could not find observed limit for m =", int(mass), "GeV. Aborting."
                        sys.exit(1)

                if args.method == 'Asymptotic' or args.method == 'HybridNewGrid':
                    if len(masses[analysis]) != len(xs_exp_limits[analysis]):
                        print "** ERROR: ** Could not find expected limit for m =", int(mass), "GeV. Aborting."
                        print "masses = ",
                        print masses[analysis]
                        print "xs_exp_limits = ",
                        print xs_exp_limits[analysis]
                        sys.exit(1)

                    if len(masses[analysis]) != len(xs_exp_limits_1sigma[analysis]):
                        print "** ERROR: ** Could not find expected 1 sigma down limit for m =", int(mass), "GeV. Aborting."
                        print "masses = ",
                        print masses[analysis]
                        print "xs_exp_limits_1sigma = ",
                        print xs_exp_limits_1sigma[analysis]
                        sys.exit(1)

                    if len(masses[analysis]) != len(xs_exp_limits_1sigma_up[analysis]):
                        print "** ERROR: ** Could not find expected 1 sigma up limit for m =", int(mass), "GeV. Aborting."
                        print "masses = ",
                        print masses[analysis]
                        print "xs_exp_limits_1sigma_up = ",
                        print xs_exp_limits_1sigma_up[analysis]
                        sys.exit(1)

                    if len(masses[analysis]) != len(xs_exp_limits_2sigma[analysis]):
                        print "** ERROR: ** Could not find expected 2 sigma down limit for m =", int(mass), "GeV. Aborting."
                        sys.exit(1)

                    if len(masses[analysis]) != len(xs_exp_limits_2sigma_up[analysis]):
                        print "** ERROR: ** Could not find expected 2 sigma up limit for m =", int(mass), "GeV. Aborting."
                        sys.exit(1)
            if args.method == 'Asymptotic' or args.method == 'HybridNewGrid':
                # complete the expected limit arrays
                for i in range(0,len(masses[analysis])):
                    masses_exp[analysis].append( masses[analysis][len(masses[analysis])-i-1] )
                    xs_exp_limits_1sigma[analysis].append( xs_exp_limits_1sigma_up[analysis][len(masses[analysis])-i-1] )
                    xs_exp_limits_2sigma[analysis].append( xs_exp_limits_2sigma_up[analysis][len(masses[analysis])-i-1] )

            if args.printResults:
                print "masses =", masses[analysis].tolist()
                print "xs_obs_limits =", xs_obs_limits[analysis].tolist()
                print "xs_exp_limits =", xs_exp_limits[analysis].tolist()
                print ""
                print "masses_exp =", masses_exp[analysis].tolist()
                print "xs_exp_limits_1sigma =", xs_exp_limits_1sigma[analysis].tolist()
                print "xs_exp_limits_2sigma =", xs_exp_limits_2sigma[analysis].tolist()


        gROOT.SetBatch(kTRUE);
        gStyle.SetOptStat(0)
        gStyle.SetOptTitle(0)
        gStyle.SetTitleFont(42, "XYZ")
        gStyle.SetTitleSize(0.05, "XYZ")
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

        # limits
        graph_exp_2sigma = {}
        graph_exp_1sigma = {}
        graph_exp = {}
        graph_obs = {}
        for analysis in analyses:
            if args.method == "Asymptotic" or args.method == "HybridNewGrid":
                graph_exp_2sigma[analysis] = ( TGraph(len(masses_exp[analysis]),masses_exp[analysis],xs_exp_limits_2sigma[analysis]) if len(xs_exp_limits_2sigma[analysis]) > 0 else TGraph(0) )
                graph_exp_2sigma[analysis].SetFillColor(kYellow)

                graph_exp_1sigma[analysis] = ( TGraph(len(masses_exp[analysis]),masses_exp[analysis],xs_exp_limits_1sigma[analysis]) if len(xs_exp_limits_2sigma[analysis]) > 0 else TGraph(0) )
                graph_exp_1sigma[analysis].SetFillColor(kGreen+1)

                graph_exp[analysis] = ( TGraph(len(masses[analysis]),masses[analysis],xs_exp_limits[analysis]) if len(xs_exp_limits_2sigma[analysis]) > 0 else TGraph(0) )
                #graph_exp[analysis].SetMarkerStyle(24)
                graph_exp[analysis].SetLineWidth(3)
                graph_exp[analysis].SetLineStyle(2)
                graph_exp[analysis].SetLineColor(4)

            graph_obs[analysis] = TGraph(len(masses[analysis]),masses[analysis],xs_obs_limits[analysis])
            graph_obs[analysis].SetMarkerStyle(20)
            graph_obs[analysis].SetLineWidth(3)
            #graph_obs[analysis].SetLineStyle(1)
            graph_obs[analysis].SetLineColor(1)

        c = TCanvas("c", "",700,500)
        c.cd()

        legend = TLegend(.7,.7,.95,.90)
        legend.SetBorderSize(0)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.03)
        legend.SetHeader('95% CL upper limits')

        frame = TH1F("frame", "frame", 100, 200, 1300)
        frame.Reset()
        frame.GetXaxis().SetTitle("Resonance mass [GeV]")
        frame.GetXaxis().SetTitleOffset(1.05)
        if args.timesAE:
            #if model == "ZPrime":
            #    frame.GetYaxis().SetTitle("#sigma #times BR(c#bar{c},b#bar{b}) #times #it{A} #times #epsilon [pb]")
            #else:
            frame.GetYaxis().SetTitle("#sigma #times BR(b#bar{b}) #times #it{A} #times #epsilon [pb]")
        else:
            #if model == "ZPrime":
            #    frame.GetYaxis().SetTitle("#sigma #times BR(c#bar{c},b#bar{b}) [pb]")
            #else:
            frame.GetYaxis().SetTitle("#sigma #times BR(b#bar{b}) [pb]")
        frame.GetYaxis().SetTitleOffset(1.2)
        if args.timesAE:
            frame.GetYaxis().SetRangeUser(1e-03,1e+01)
        else:
            frame.GetYaxis().SetRangeUser(1e-01,5e+02)
        frame.Draw("axis")

        for analysis in analyses:
            graph_exp_2sigma[analysis].Draw("F")
            graph_exp_1sigma[analysis].Draw("F")
            graph_exp[analysis].Draw("L")
        for analysis in analyses:
            graph_obs[analysis].Draw("LP")

        legend.AddEntry(graph_obs["trigbbl_CSVTM"],"Observed","lp")
        legend.AddEntry(graph_exp["trigbbl_CSVTM"],"Expected","lp")
        legend.AddEntry(graph_exp_1sigma["trigbbl_CSVTM"],"#pm 1#sigma","F")
        legend.AddEntry(graph_exp_2sigma["trigbbl_CSVTM"],"#pm 2#sigma","F")
 
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
        if args.fitTrigger:
            postfix += "_fitTrigger"
        elif args.correctTrigger:
            postfix += "_correctTrigger"
        if args.useMCTrigger:
            postfix += "_useMCTrigger"
        fileName = limit_config.paths["limit_plots"] + '/xs_limit_%s_%s_%s_%s.%s'%(args.method,"combined", model + postfix, args.fit_function, args.fileFormat.lower())
        if args.timesAE:
            fileName = fileName.replace("xs_limit", "xsAE_limit")

        # Model names on plot
        Root.myText(0.2, 0.2, 1, model_names_pretty[model], 0.5)
        c.SaveAs(fileName)
        print "Plot saved to '%s'"%(fileName)

        if model == "ZPrime":
            # Draw theory curve
            zp_xs_file = open(os.path.expandvars("$CMSSW_BASE/src/CMSDIJET/StatisticalTools/data/xs.pkl"), 'r')
            zp_xs = pickle.load(zp_xs_file) # xs.pkl is produced by zeta_plot.py. 
            print zp_xs
            zp_xs_masses = sorted(zp_xs.keys())
            zp_xsbr_graph = TGraph(len(zp_xs_masses))
            for i, zp_mass in enumerate(zp_xs_masses):
                this_zp_xs = zp_xs[zp_mass]["pp"]
                if zp_mass > 2 * 172.:
                    tpart = (1. + 2*172**2/zp_mass**2) * (1. - 4*172 **2/zp_mass**2)**0.5
                else:
                    tpart = 0.
                brbb = 1. / (5. + tpart)
                this_zp_xsbr = this_zp_xs * brbb
                zp_xsbr_graph.SetPoint(i, zp_mass, this_zp_xsbr*1.2)
            zp_xsbr_graph.SetLineColor(2)
            zp_xsbr_graph.SetLineStyle(2)
            zp_xsbr_graph.SetMarkerStyle(20)
            zp_xsbr_graph.SetMarkerSize(0)
            zp_xsbr_graph.Draw("pl")
            legend.AddEntry(zp_xsbr_graph, "Z', g_{q}=0.25", "l")
            legend.Draw()
            fileNameZp = fileName + ".zp.pdf"
            c.SaveAs(fileNameZp)
        graph_obs["trigbbl_CSVTM"].Print("all")
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
