#!/usr/bin/env python

import sys, os
from argparse import ArgumentParser


def main():
    # usage description
    usage = "Example: ./scripts/runCombine.py -M Asymptotic -d datacards -f qq --massrange 1200 7000 100"

    # input parameters
    parser = ArgumentParser(description='Script that runs limit calculation for specified mass points',epilog=usage)

    parser.add_argument("-M", "--method", dest="method", required=True,
                        choices=['ProfileLikelihood', 'Asymptotic', 'MarkovChainMC'],
                        help="Method to calculate upper limits",
                        metavar="METHOD")

    parser.add_argument("-d", "--datacards_path", dest="datacards_path", required=True,
                        help="Path to datacards and workspaces",
                        metavar="DATA_HISTNAME")

    parser.add_argument("-f", "--final_state", dest="final_state", required=True,
                        help="Final state (e.g. qq, qg, gg)",
                        metavar="FINAL_STATE")

    parser.add_argument("-o", "--output_path", dest="output_path",
                        default='logs',
                        help="Output path where log files will be stored (default: %(default)s)",
                        metavar="OUTPUT_PATH")

    parser.add_argument("--rMax", dest="rMax", type=float, help="Maximum value for signal strength")

    parser.add_argument("--tries", dest="tries", type=int, default=10, help="Number of times to run the MCMC (default: %(default)i)")

    parser.add_argument("--noSyst", dest="noSyst", default=False, action="store_true", help="Run without systematic uncertainties")

    parser.add_argument("--signif", dest="signif", default=False, action="store_true", help="Calculate significance instead of limits")

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

    logs_path = os.path.join(os.getcwd(),args.output_path)

    # change directory to where datacards are stored
    os.chdir(args.datacards_path)

    method = args.method

    if method != 'ProfileLikelihood' and args.signif:
        print "** ERROR: ** For significance calculation the ProfileLikelihood method has to be used. Aborting."
        sys.exit(1)

    options = ''
    if args.signif:
        options = options + ' --signif'
    if args.noSyst:
        options = options + ' --systematics 0'
    if method != 'ProfileLikelihood' and not args.rMax != None:
        options = options + ' --hintMethod ProfileLikelihood'
    if args.rMax != None:
        options = options + ' --rMin 0 --rMax %.1f'%(args.rMax)
    if method == 'MarkovChainMC':
        options = options + ' --tries %i --proposal fit'%(args.tries)

    for mass in masses:

        print ">> Calculating %s for %s resonance with m = %i GeV..."%(('significance' if args.signif else 'limits'), args.final_state, int(mass))

        logName = '%s_%s_m%i.log'%(('significance' if args.signif else 'limits'), args.final_state, int(mass))

        run_options = options + ' --name _%s_m%i --mass %i'%(args.final_state,int(mass),int(mass))

        cmd = "combine -M %s %s datacard_%s_m%i.txt | tee %s"%(method,run_options,args.final_state,int(mass),os.path.join(logs_path,logName))
        print "---------------------------------------------------------------------------"
        print "Running: " + cmd +"\n"
        os.system(cmd)


if __name__ == '__main__':
    main()

