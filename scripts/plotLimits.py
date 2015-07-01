#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from array import array
import numpy as np
import CMS_lumi


def main():
    # usage description
    usage = "Example: ./scripts/plotLimits.py -l logs -f qq --massrange 1200 7000 100"

    # input parameters
    parser = ArgumentParser(description='Script that plots limits for specified mass points',epilog=usage)

    results_group = parser.add_mutually_exclusive_group(required=True)
    results_group.add_argument("-l", "--logs_path", dest="logs_path",
                               help="Path to log files",
                               metavar="LOGS_PATH")
    results_group.add_argument("-r", "--results_file", dest="results_file",
                               help="Path to a file containing results",
                               metavar="RESULTS_FILE")

    parser.add_argument("-f", "--final_state", dest="final_state", required=True,
                        help="Final state (e.g. qq, qg, gg)",
                        metavar="FINAL_STATE")

    parser.add_argument("--postfix", dest="postfix", default='', help="Postfix for the output plot name (default: %(default)s)")

    parser.add_argument("--fileFormat", dest="fileFormat", default='pdf', help="Format of the output plot (default: %(default)s)")

    parser.add_argument("--extraText", dest="extraText", default='Simulation Preliminary', help="Extra text on the plot (default: %(default)s)")

    parser.add_argument("--lumi_sqrtS", dest="lumi_sqrtS", default='1 fb^{-1} (13 TeV)', help="Integrated luminosity and center-of-mass energy (default: %(default)s)")

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

    # arrays holding results
    masses = array('d')
    xs_obs_limits = array('d')
    xs_exp_limits = array('d')
    masses_exp = array('d')
    xs_exp_limits_1sigma = array('d')
    xs_exp_limits_1sigma_up = array('d')
    xs_exp_limits_2sigma = array('d')
    xs_exp_limits_2sigma_up = array('d')

    if args.logs_path != None:

        logs_path = os.path.join(os.getcwd(),args.logs_path)

        for mass in input_masses:

            print ">> Reading results for %s resonance with m = %i GeV..."%(args.final_state, int(mass))

            masses.append(mass)
            masses_exp.append(mass)

            logName = 'limits_%s_m%i.log'%(args.final_state,int(mass))

            log_file = open(os.path.join(logs_path,logName),'r')

            # read the log file
            for line in log_file:
                if re.search("^Observed Limit: r", line):
                  xs_obs_limits.append(float(line.split()[-1]))
                if re.search("^Expected 50.0%: r", line):
                  xs_exp_limits.append(float(line.split()[-1]))
                if re.search("^Expected 16.0%: r", line):
                  xs_exp_limits_1sigma.append(float(line.split()[-1]))
                if re.search("^Expected 84.0%: r", line):
                  xs_exp_limits_1sigma_up.append(float(line.split()[-1]))
                if re.search("^Expected  2.5%: r", line):
                  xs_exp_limits_2sigma.append(float(line.split()[-1]))
                if re.search("^Expected 97.5%: r", line):
                  xs_exp_limits_2sigma_up.append(float(line.split()[-1]))

            if len(masses) != len(xs_obs_limits):
                print "** ERROR: ** Could not find observed limit for m =", int(mass), "GeV. Aborting."
                sys.exit(1)

            if len(masses) != len(xs_exp_limits):
                print "** ERROR: ** Could not find expected limit for m =", int(mass), "GeV. Aborting."
                sys.exit(1)

            if len(masses) != len(xs_exp_limits_1sigma):
                print "** ERROR: ** Could not find expected 1 sigma down limit for m =", int(mass), "GeV. Aborting."
                sys.exit(1)

            if len(masses) != len(xs_exp_limits_1sigma_up):
                print "** ERROR: ** Could not find expected 1 sigma up limit for m =", int(mass), "GeV. Aborting."
                sys.exit(1)

            if len(masses) != len(xs_exp_limits_2sigma):
                print "** ERROR: ** Could not find expected 2 sigma down limit for m =", int(mass), "GeV. Aborting."
                sys.exit(1)

            if len(masses) != len(xs_exp_limits_2sigma_up):
                print "** ERROR: ** Could not find expected 2 sigma up limit for m =", int(mass), "GeV. Aborting."
                sys.exit(1)

        # complete the expected limit arrays
        for i in range(0,len(masses)):
            masses_exp.append( masses[len(masses)-i-1] )
            xs_exp_limits_1sigma.append( xs_exp_limits_1sigma_up[len(masses)-i-1] )
            xs_exp_limits_2sigma.append( xs_exp_limits_2sigma_up[len(masses)-i-1] )
    else:
        print ">> Importing results..."

        sys.path.insert(0, os.path.dirname(args.results_file))

        results = __import__(os.path.basename(args.results_file).replace(".py",""))

        all_masses = np.array(results.masses)
        all_masses_exp = np.array(results.masses_exp)
        indices = []
        indices_exp = []

        # search for indices of input_masses
        for mass in input_masses:
            where = np.where(all_masses==mass)[0]
            if len(where) == 0:
                print "** WARNING: ** Cannot find results for m =", int(mass), "GeV in the provided results file. Skipping this mass point."
            indices.extend( where )
            where = np.where(all_masses_exp==mass)[0]
            if len(where) == 0:
                print "** WARNING: ** Cannot find results for m =", int(mass), "GeV in the provided results file. Skipping this mass point."
            indices_exp.extend( where )

        # sort indices
        indices.sort()
        indices_exp.sort()

        for i in indices:
            masses.append( results.masses[i] )
            xs_obs_limits.append( results.xs_obs_limits[i] )
            xs_exp_limits.append( results.xs_exp_limits[i] )
        for i in indices_exp:
            masses_exp.append( results.masses_exp[i] )
            xs_exp_limits_1sigma.append( results.xs_exp_limits_1sigma[i] )
            xs_exp_limits_2sigma.append( results.xs_exp_limits_2sigma[i] )

    if args.printResults:
        print "masses =", masses.tolist()
        print "xs_obs_limits =", xs_obs_limits.tolist()
        print "xs_exp_limits =", xs_exp_limits.tolist()
        print ""
        print "masses_exp =", masses_exp.tolist()
        print "xs_exp_limits_1sigma =", xs_exp_limits_1sigma.tolist()
        print "xs_exp_limits_2sigma =", xs_exp_limits_2sigma.tolist()

    # import ROOT stuff
    from ROOT import kTRUE, kFALSE, gROOT, gStyle, gPad, TGraph, TCanvas, TLegend
    from ROOT import kGreen, kYellow, kWhite

    CMS_lumi.extraText = args.extraText
    CMS_lumi.lumi_sqrtS = args.lumi_sqrtS # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    iPos = 11
    iPeriod = 0

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

    graph_exp_2sigma = TGraph(len(masses_exp),masses_exp,xs_exp_limits_2sigma)
    graph_exp_2sigma.SetFillColor(kYellow)
    graph_exp_2sigma.GetXaxis().SetTitle("%s resonance mass [GeV]"%(args.final_state))
    graph_exp_2sigma.GetYaxis().SetTitle("#sigma#timesBR(X#rightarrowjj)#timesA [pb]")
    graph_exp_2sigma.GetYaxis().SetRangeUser(1e-03,1e+02)
    #graph_exp_2sigma.GetXaxis().SetNdivisions(1005)

    graph_exp_1sigma = TGraph(len(masses_exp),masses_exp,xs_exp_limits_1sigma)
    graph_exp_1sigma.SetFillColor(kGreen+1)

    graph_exp = TGraph(len(masses),masses,xs_exp_limits)
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

    graph_exp_2sigma.Draw("AF")
    graph_exp_1sigma.Draw("F")
    graph_exp.Draw("L")
    graph_obs.Draw("LP")

    legend = TLegend(.60,.55,.90,.70)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)
    legend.SetHeader('95% CL upper limits')
    legend.AddEntry(graph_obs,"Observed","lp")
    legend.AddEntry(graph_exp,"Expected","lp")
    legend.AddEntry(graph_exp_1sigma,"#pm 1#sigma","F")
    legend.AddEntry(graph_exp_2sigma,"#pm 2#sigma","F")
    legend.Draw()

    # draw the lumi text on the canvas
    CMS_lumi.CMS_lumi(c, iPeriod, iPos)

    gPad.RedrawAxis()

    c.SetLogy()
    fileName = 'xs_limit_%s.%s'%(args.final_state + ( ('_' + args.postfix) if args.postfix != '' else '' ), args.fileFormat.lower())
    c.SaveAs(fileName)
    print "Plot saved to '%s'"%(fileName)


if __name__ == '__main__':
    main()
