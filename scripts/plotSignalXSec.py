#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from array import array
import numpy as np
import CMS_lumi


def main():
    # usage description
    usage = "Example: ./scripts/plotSignificance.py -l logs -f qq --massrange 1200 6000 100"

    # input parameters
    parser = ArgumentParser(description='Script that plots significance for specified mass points',epilog=usage)

    parser.add_argument("-M", "--method", dest="method",
                        choices=['MaxLikelihoodFit'],
                        default='MaxLikelihoodFit',
                        help="Method to calculate upper limits",
                        metavar="METHOD")

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
    sig = array('d')
    sig_ex = array('d')
    sig_eyl = array('d')
    sig_eyh = array('d')

    if args.logs_path != None:

        logs_path = os.path.join(os.getcwd(),args.logs_path)

        for mass in input_masses:

            print ">> Reading results for %s resonance with m = %i GeV..."%(args.final_state, int(mass))

            masses.append(mass)

            logName = 'signal_xs_%s_%s_m%i.log'%(args.method,args.final_state,int(mass))

            log_file = open(os.path.join(logs_path,logName),'r')

            # read the log file
            for line in log_file:
                if re.search("^Best fit r:", line):
                  sig.append(float(line.split()[3]))
                  sig_eyl.append(float(line.split()[4].split('/')[0].lstrip('-')))
                  sig_eyh.append(float(line.split()[4].split('/')[1].lstrip('+')))

            sig_ex.append(0.)

            if len(masses) != len(sig):
                print "** WARNING: ** Fit failed for m =", int(mass), "GeV. Setting signal cross section to 0."
                sig.append(0.)
                sig_eyl.append(0.)
                sig_eyh.append(0.)
    else:
        print ">> Importing results..."

        sys.path.insert(0, os.path.dirname(args.results_file))

        results = __import__(os.path.basename(args.results_file).replace(".py",""))

        all_masses = np.array(results.masses)
        indices = []

        # search for indices of input_masses
        for mass in input_masses:
            where = np.where(all_masses==mass)[0]
            if len(where) == 0:
                print "** WARNING: ** Cannot find results for m =", int(mass), "GeV in the provided results file. Skipping this mass point."
            indices.extend( where )

        # sort indices
        indices.sort()

        for i in indices:
            masses.append( results.masses[i] )
            sig.append( results.sig[i] )
            sig_ex.append( results.sig_ex[i] )
            sig_eyl.append( results.sig_eyl[i] )
            sig_eyh.append( results.sig_eyh[i] )


    if args.printResults:
        print "masses =", masses.tolist()
        print "sig =", sig.tolist()
        print "sig_ex =", sig_ex.tolist()
        print "sig_eyl =", sig_eyl.tolist()
        print "sig_eyh =", sig_eyh.tolist()

    # create final arrays
    sig_pos = array('d')
    sig_exl = array('d')
    sig_exh = array('d')

    # fill final arrays
    for i in range(0,len(masses)):
        sig_pos.append(sig[i] if sig[i]>0. else 0.)
        sig_exl.append(sig_ex[i])
        sig_exh.append(sig_ex[i])
        sig_eyl.append(sig_eyl[i] if sig[i]>0. else 0.)
        sig_eyh.append(sig_eyh[i])

    # import ROOT stuff
    from ROOT import kTRUE, kFALSE, gROOT, gStyle, gPad, TGraphAsymmErrors, TCanvas, TLegend
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

    graph_sig = TGraphAsymmErrors(len(masses),masses,sig_pos,sig_exl,sig_exh,sig_eyl,sig_eyh)
    graph_sig.GetXaxis().SetTitle("%s resonance mass [GeV]"%(args.final_state))
    graph_sig.GetYaxis().SetTitle("Signal cross section [pb]")
    graph_sig.GetYaxis().SetTitleOffset(1.2)
    graph_sig.GetYaxis().SetRangeUser(1e-4,2e2)
    graph_sig.SetMarkerStyle(20)
    graph_sig.SetMarkerColor(1)
    graph_sig.SetLineWidth(2)
    graph_sig.SetLineStyle(1)
    graph_sig.SetLineColor(1)

    c = TCanvas("c", "",800,800)
    c.cd()

    graph_sig.Draw("AP")

    # draw the lumi text on the canvas
    CMS_lumi.extraText = args.extraText
    CMS_lumi.lumi_sqrtS = args.lumi_sqrtS # used with iPeriod = 0 (free form)
    iPos = 11
    iPeriod = 0

    CMS_lumi.CMS_lumi(c, iPeriod, iPos)

    gPad.RedrawAxis()

    c.SetLogy()
    c.SetGridx()
    c.SetGridy()
    fileName = 'signal_xs_%s_%s.%s'%(args.method,args.final_state + ( ('_' + args.postfix) if args.postfix != '' else '' ), args.fileFormat.lower())
    c.SaveAs(fileName)
    print "Plot saved to '%s'"%(fileName)


if __name__ == '__main__':
    main()
