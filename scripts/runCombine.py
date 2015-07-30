#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser

jdl_template = """universe = vanilla
Notification = never
Executable = condor/run_DUMMY_JOB.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = DUMMY_FILES
Output = DUMMY_OUTPUTDIR/condor/condor_DUMMY_JOB_$(Cluster)_$(Process).stdout
Error = DUMMY_OUTPUTDIR/condor/condor_DUMMY_JOB_$(Cluster)_$(Process).stderr
Log = DUMMY_OUTPUTDIR/condor/condor_DUMMY_JOB_$(Cluster)_$(Process).log
Arguments = DUMMY_OUTPUTDIR
#+LENGTH="SHORT"
Queue 1
"""

bash_template = """#!/bin/bash

OUTPUTDIR=$1

START_TIME=`/bin/date`
echo "Started at $START_TIME"
echo ""

export SCRAM_ARCH=slc6_amd64_gcc481
source /cvmfs/cms.cern.ch/cmsset_default.sh

cd $OUTPUTDIR
eval `scramv1 runtime -sh`

# change to Condor scratch directory
cd ${_CONDOR_SCRATCH_DIR}

echo "Running: DUMMY_CMD"
DUMMY_CMD
exitcode=$?

echo ""
END_TIME=`/bin/date`
echo "Finished at $END_TIME"
exit $exitcode
"""


def main():
    # usage description
    usage = "Example: ./scripts/runCombine.py -M Asymptotic -d datacards -f qq --massrange 1200 7000 100"

    # input parameters
    parser = ArgumentParser(description='Script that runs limit calculation for specified mass points',epilog=usage)

    parser.add_argument("-M", "--method", dest="method", required=True,
                        choices=['MaxLikelihoodFit', 'ProfileLikelihood', 'HybridNew', 'Asymptotic', 'MarkovChainMC'],
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

    parser.add_argument("--rMin", dest="rMin", type=float, help="Minimum value for the signal strength")

    parser.add_argument("--rMax", dest="rMax", type=float, help="Maximum value for the signal strength")

    parser.add_argument("--toysH", dest="toysH", type=int, help="Number of Toy MC extractions for HybridNew")

    parser.add_argument("--tries", dest="tries", type=int, default=10, help="Number of times to run the MCMC (default: %(default)i)")

    parser.add_argument("--proposal", dest="proposal", default='ortho', help="Proposal function for MCMC (default: %(default)s)")

    parser.add_argument("--noSyst", dest="noSyst", default=False, action="store_true", help="Run without systematic uncertainties")

    parser.add_argument("--noHint", dest="noHint", default=False, action="store_true", help="Do not run the hint method")

    parser.add_argument("--signif", dest="signif", default=False, action="store_true", help="Calculate significance instead of limits")

    parser.add_argument("--fitStrategy", dest="fitStrategy", type=int, help="Fit strategy (default: %(default).1f)")

    parser.add_argument("--condor", dest="condor", default=False, action="store_true", help="Batch process using Condor")

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
    if not os.path.isdir( os.path.join(os.getcwd(), args.output_path) ):
        os.mkdir( os.path.join(os.getcwd(), args.output_path) )
    print os.getcwd()
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

    method = args.method

    if args.signif and method != 'ProfileLikelihood' and method != 'HybridNew':
        print "** ERROR: ** For significance calculation the ProfileLikelihood or HybridNew method has to be used. Aborting."
        sys.exit(1)

    options = ''
    if args.signif:
        options = options + ' --signif'
    if args.fitStrategy:
        options = options + ' --minimizerStrategy %i'%(args.fitStrategy)
    if args.noSyst:
        options = options + ' --systematics 0'
    if method != 'ProfileLikelihood' and method != 'MaxLikelihoodFit' and args.rMax == None and not args.noHint and not args.signif:
        options = options + ' --hintMethod ProfileLikelihood'
    if method == 'HybridNew' and args.toysH != None:
        options = options + ' --toysH %i'%(args.toysH)
    if args.rMin != None:
        options = options + ' --rMin %f'%(args.rMin)
    if args.rMax != None:
        options = options + ' --rMax %f'%(args.rMax)
    if method == 'MarkovChainMC':
        options = options + ' --tries %i --proposal %s'%(args.tries, args.proposal)

    prefix = 'limits'
    if args.signif:
        prefix = 'significance'
    elif method == 'MaxLikelihoodFit':
        prefix = 'signal_xs'

    datacards_path = os.path.join(os.getcwd(), args.datacards_path)
    output_path = os.path.join(os.getcwd(), args.output_path)
    condor_path = os.path.join(output_path, 'condor')

    # change to the appropriate directory
    if args.condor:
        os.chdir(output_path)
    else:
        os.chdir(datacards_path)

    for mass in masses:

        logName = '%s_%s_%s_m%i.log'%(prefix, method, args.final_state, int(mass))

        run_options = options + ' --name _%s_m%i --mass %i'%(args.final_state,int(mass),int(mass))

        cmd = "combine -M %s %s datacard_%s_m%i.txt | tee %s"%(method,run_options,args.final_state,int(mass),os.path.join(('' if args.condor else output_path),logName))

        # if using Condor
        if args.condor:
            # check if the Condor directory exists
            if not os.path.isdir( condor_path ):
                os.mkdir( condor_path )

            # create the jdl file
            jdl_content = jdl_template
            jdl_content = re.sub('DUMMY_JOB','%s_%s_m%i'%(method,args.final_state,int(mass)),jdl_content)
            jdl_content = re.sub('DUMMY_OUTPUTDIR',os.getcwd(),jdl_content)
            files_to_transfer = []
            files_to_transfer.append( os.path.join(datacards_path, 'datacard_%s_m%i.txt'%(args.final_state,int(mass))) )
            files_to_transfer.append( os.path.join(datacards_path, 'workspace_%s_m%i.root'%(args.final_state,int(mass))) )
            jdl_content = re.sub('DUMMY_FILES',', '.join(files_to_transfer),jdl_content)

            jdl_file = open(os.path.join(condor_path,'run_%s_%s_m%i.jdl'%(method,args.final_state,int(mass))),'w')
            jdl_file.write(jdl_content)
            jdl_file.close()

            # create the Bash script
            bash_content = bash_template
            bash_content = re.sub('DUMMY_CMD',cmd,bash_content)

            bash_script = open(os.path.join(condor_path,'run_%s_%s_m%i.sh'%(method,args.final_state,int(mass))),'w')
            bash_script.write(bash_content)
            bash_script.close()

            print ">> Submitting job for %s resonance with m = %i GeV..."%(args.final_state, int(mass))
            print "---------------------------------------------------------------------------"
            condor_cmd = 'condor_submit condor/run_%s_%s_m%i.jdl'%(method,args.final_state,int(mass))
            print "Running: " + condor_cmd + "\n"
            os.system(condor_cmd)
        else:
            print ">> Running combine for %s resonance with m = %i GeV..."%(args.final_state, int(mass))
            print "---------------------------------------------------------------------------"
            print "Running: " + cmd + "\n"
            os.system(cmd)


if __name__ == '__main__':
    main()

