#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from array import array
import numpy as np
import CMS_lumi
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
import ROOT
from ROOT import *
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")


def main():
    # usage description
    usage = "Example: ./scripts/plotSignificance.py -l logs -f qq --massrange 1200 6000 100"

    # input parameters
    parser = ArgumentParser(description='Script that plots significance for specified mass points',epilog=usage)
    parser.add_argument('analysis', type=str, help='Analysis name')
    parser.add_argument('model', type=str, help='Model name')

    parser.add_argument("-M", "--method", dest="method",
                        choices=['ProfileLikelihood', 'HybridNew', 'theta'],
                        default='ProfileLikelihood',
                        help="Method to calculate upper limits",
                        metavar="METHOD")
    parser.add_argument('--fit_function', type=str, default="f4", help="Name of fit function used for background estimate")

    parser.add_argument("--sigRange", dest="sigRange", type=float, default=2.5, help="Significance range to plot (default: %(default)f)")

    parser.add_argument("--postfix", dest="postfix", default='', help="Postfix for the output plot name (default: %(default)s)")

    parser.add_argument("--fileFormat", dest="fileFormat", default='pdf', help="Format of the output plot (default: %(default)s)")

    parser.add_argument("--extraText", dest="extraText", default='Preliminary', help="Extra text on the plot (default: %(default)s)")

    parser.add_argument("--lumi_sqrtS", dest="lumi_sqrtS", default='19.7 fb^{-1} (13 TeV)', help="Integrated luminosity and center-of-mass energy (default: %(default)s)")

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
    significances = array('d')
    p0_values = array('d')

    for mass in input_masses:
        masses.append(mass)
        print ">> Reading results for resonance with m = %i GeV..."%(int(mass))
        log_file = open(limit_config.get_combine_log_path(args.analysis, args.model, mass, args.fit_function, args.method, what="significance"), 'r')

        if args.method == 'theta': logName = logName.replace('significance_','')

        # read the log file
        for line in log_file:
            if args.method == 'theta':
                if re.search("^{'signal': {'Z':", line):
                    significances.append(float(line.split()[-1].lstrip('[').rstrip('}').rstrip(']')))
            else:
                if re.search("^Significance:", line):
                    significances.append(float(line.split()[1]))
                elif re.search("^Null p-value:", line):
                    p0_values.append(float(line.split()[2]))

        if len(masses) != len(significances):
            print "** ERROR: ** Could not find significance for m =", int(mass), "GeV. Aborting."
            sys.exit(1)
        if len(masses) != len(p0_values):
            print "** ERROR: ** Could not find p0 value for m =", int(mass), "GeV. Aborting."
            sys.exit(1)
        # Allow only positive fluctuations
        if significances[-1] < 0:
            significances[-1] = 0.

    if args.printResults:
        print "masses =", masses.tolist()
        print "significances =", significances.tolist()
        print "p0_values =", p0_values.tolist()

    # import ROOT stuff
    from ROOT import kTRUE, kFALSE, gROOT, gStyle, gPad, TGraph, TCanvas, TLegend
    from ROOT import kGreen, kYellow, kWhite, kRed, kBlue

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

    graph_sig = TGraph(len(masses),masses,significances)
    graph_sig.GetXaxis().SetTitle("Resonance mass [GeV]")
    graph_sig.GetYaxis().SetTitle("Significance (local)")
    graph_sig.GetYaxis().SetTitleOffset(1.2)
    graph_sig.GetYaxis().SetRangeUser(0.,args.sigRange)
    graph_sig.SetLineWidth(2)
    graph_sig.SetLineColor(kRed)
    graph_sig.SetMarkerStyle(21)
    graph_sig.SetMarkerSize(1)
    graph_sig.SetMarkerColor(kBlue)

    c = TCanvas("c", "",800,800)
    c.cd()

    graph_sig.Draw("ALP")

    # draw the lumi text on the canvas
    CMS_lumi.extraText = args.extraText
    CMS_lumi.lumi_sqrtS = args.lumi_sqrtS # used with iPeriod = 0 (free form)
    iPos = 11
    iPeriod = 0

    CMS_lumi.CMS_lumi(c, iPeriod, iPos)

    gPad.RedrawAxis()

    c.SetGridx()
    c.SetGridy()
    fileName = limit_config.paths["limit_plots"] + '/significance_%s_%s_%s_%s.%s'%(args.method,args.analysis, args.model + args.postfix, args.fit_function, args.fileFormat.lower())
    c.SaveAs(fileName)
    print "Significance plot saved to '%s'"%(fileName)

    graph_p0 = TGraph(len(masses), masses, p0_values)
    graph_p0.GetXaxis().SetTitle("Resonance mass [GeV]")
    graph_p0.GetYaxis().SetTitle("p_{0} (local)")
    graph_p0.GetYaxis().SetTitleOffset(1.2)
    graph_p0.GetYaxis().SetRangeUser(1.e-3, 1.)
    graph_p0.SetLineWidth(2)
    graph_p0.SetLineColor(kRed)
    graph_p0.SetMarkerStyle(21)
    graph_p0.SetMarkerSize(1)
    graph_p0.SetMarkerColor(kBlue)
    c_p0 = TCanvas("c_p0", "",800,800)
    c_p0.SetLogy()
    c_p0.cd()

    graph_p0.Draw("ALP")

    # draw the lumi text on the canvas
    #CMS_lumi.extraText = args.extraText
    #CMS_lumi.lumi_sqrtS = args.lumi_sqrtS # used with iPeriod = 0 (free form)
    #iPos = 0
    #iPeriod = 0

    #CMS_lumi.CMS_lumi(c_p0, iPeriod, iPos)

    Root.CMSLabel(0.2, 0.2, "Preliminary", 1, 0.65)

    gPad.RedrawAxis()

    c_p0.SetGridx()
    c_p0.SetGridy()
    fileName = limit_config.paths["limit_plots"] + '/p0_%s_%s_%s_%s.%s'%(args.method,args.analysis, args.model + args.postfix, args.fit_function, args.fileFormat.lower())
    c_p0.SaveAs(fileName)
    print "p0 plot saved to '%s'"%(fileName)

if __name__ == '__main__':
    main()
