#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
import CMSDIJET.StatisticalTools.limit_configuration as limit_config
from CMSDIJET.StatisticalTools.systematics import *
from CMSDIJET.StatisticalTools.roofit_functions import *

jdl_template = """universe = vanilla
Notification = never
Executable = DUMMY_EXECUTABLE
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = DUMMY_FILES
Output = DUMMY_OUTPUTDIR/condor_DUMMY_JOB_$(Cluster)_$(Process).stdout
Error = DUMMY_OUTPUTDIR/condor_DUMMY_JOB_$(Cluster)_$(Process).stderr
Log = DUMMY_OUTPUTDIR/condor_DUMMY_JOB_$(Cluster)_$(Process).log
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
    parser.add_argument('--analyses', type=str, default="trigbbh_CSVTM,trigbbl_CSVTM", help="Analysis names")
    parser.add_argument('--models', type=str, default="Hbb,RSG", help="Model names")
    parser.add_argument('--fitSignal', action='store_true', help="Use fitted signal shapes, rather than histograms")
    parser.add_argument('--fit_function', type=str, default="f4", help="Name of central fit function")
    #parser.add_argument("-d", "--datacards_path", dest="datacards_path", required=True,
    #                    help="Path to datacards and workspaces",
    #                    metavar="DATA_HISTNAME")

    #parser.add_argument("-f", "--final_state", dest="final_state", required=True,
    #                    help="Final state (e.g. qq, qg, gg)",
    #                    metavar="FINAL_STATE")

    #parser.add_argument("-o", "--output_path", dest="output_path",
    #                    default='logs',
    #                    help="Output path where log files will be stored (default: %(default)s)",
    #                    metavar="OUTPUT_PATH")

    parser.add_argument("--rMin", dest="rMin", type=float, help="Minimum value for the signal strength")

    parser.add_argument("--rMax", dest="rMax", type=float, help="Maximum value for the signal strength")

    parser.add_argument("--rPrior", action="store_true", help="Take rMin and rMax from a previous iteration of limit setting.")

    parser.add_argument("--rRelAcc", dest="rRelAcc", type=float, help="rRelAcc")

    parser.add_argument("--rAbsAcc", dest="rAbsAcc", type=float, help="rAbsAcc")

    parser.add_argument("--level", dest="level", type=float, help="Exp from grid level")

    parser.add_argument("--toysH", dest="toysH", type=int, help="Number of Toy MC extractions for HybridNew")

    parser.add_argument("--toys", dest="toys", type=int, help="Number of Toy MC extractions, not for HybridNew")

    parser.add_argument("--tries", dest="tries", type=int, default=10, help="Number of times to run the MCMC (default: %(default)i)")

    parser.add_argument("--robustFit", dest="robustFit", action='store_true', help="combine robustFit option")

    parser.add_argument("--proposal", dest="proposal", default='ortho', help="Proposal function for MCMC (default: %(default)s)")

    parser.add_argument("--noSyst", dest="noSyst", default=False, action="store_true", help="Run without systematic uncertainties")

    parser.add_argument("--picky", dest="picky", default=False, action="store_true", help="combine picky mode")

    parser.add_argument("--saveHybridResult", action="store_true", help="--saveHybridResult")

    parser.add_argument("--strictBounds", dest="strictBounds", default=False, action="store_true", help="Strict bounds on rMax")

    parser.add_argument("--testStat", type=str, help="Test statistics: LEP, TEV, LHC (previously known as Atlas), Profile.")

    parser.add_argument("--freezeNuisances", dest="freezeNuisances", type=str, help="Freeze nuisance parameters")
    parser.add_argument("--noHint", dest="noHint", default=False, action="store_true", help="Do not run the hint method")

    parser.add_argument("--significance", dest="significance", default=False, action="store_true", help="Calculate significance instead of limits")

    parser.add_argument("--frequentist", dest="freq", default=False, action="store_true", help="Frequentist hybridnew")

    parser.add_argument("--fork", dest="forkvar", default=False, action="store_true", help="More cores")

    parser.add_argument("--fitStrategy", dest="fitStrategy", type=int, help="Fit strategy (default: %(default).1f)")

    parser.add_argument("--condor", dest="condor", default=False, action="store_true", help="Batch process using Condor")

    parser.add_argument("--no_retar", dest='no_retar', action="store_true", help="If --condor is specified, this specifies no retar of src directory (use cached copy)")

    parser.add_argument("--postfix", dest="postfix", default='', help="Postfix for the input and output file names (default: %(default)s)")

    parser.add_argument("--hnsig", dest="hnsig", default=False, action="store_true", help="HybridNew Significance calc")

    parser.add_argument("--hnsig2", dest="hnsig2", default=False, action="store_true", help="HybridNew Significance calc")

    parser.add_argument('-v', '--verbose', type=int, help='Verbosity of combine')
    parser.add_argument('-m', '--masses', type=str, help='Manually specify masses (comma-separated list). Otherwise, taken from limit_configuration.')
    #mass_group = parser.add_mutually_exclusive_group(required=True)
    #mass_group.add_argument("--mass",
    #                        type=int,
    #                        nargs = '*',
    #                        default = 1000,
    #                        help="Mass can be specified as a single value or a whitespace separated list (default: %(default)i)"
    #                        )
    #mass_group.add_argument("--massrange",
    #                        type=int,
    #                        nargs = 3,
    #                        help="Define a range of masses to be produced. Format: min max step",
    #                        metavar = ('MIN', 'MAX', 'STEP')
    #                        )
    #mass_group.add_argument("--masslist",
    #                        help = "List containing mass information"
    #                        )

    args = parser.parse_args()

    analyses = args.analyses.split(",")
    models = args.models.split(",")

    # check if the output directory exists
    #if not os.path.isdir( os.path.join(os.getcwd(), args.output_path) ):
    #    os.mkdir( os.path.join(os.getcwd(), args.output_path) )
    #print os.getcwd()
    # mass points for which resonance shapes will be produced
    #masses = []

    #if args.massrange != None:
    #    MIN, MAX, STEP = args.massrange
    #    masses = range(MIN, MAX+STEP, STEP)
    #elif args.masslist != None:
    #    # A mass list was provided
    #    print  "Will create mass list according to", args.masslist
    #    masslist = __import__(args.masslist.replace(".py",""))
    #    masses = masslist.masses
    #else:
    #    masses = args.mass

   # ## sort masses
    #masses.sort()
    
    masses = {}
    if args.masses:
        for analysis in analyses:
            masses[analysis] = [int(x) for x in args.masses.split(",")]
    else:
        masses = limit_config.limit_signal_masses
    method = args.method

    if args.significance and method != 'ProfileLikelihood' and method != 'HybridNew':
        print "** ERROR: ** For significance calculation the ProfileLikelihood or HybridNew method has to be used. Aborting."
        sys.exit(1)

    options = ''
    if args.significance:
        options = options + ' --significance'
    if args.freq:
        options = options + ' --frequentist'
    if args.hnsig:
        options = options + ' --significance --saveToys --fullBToys --saveHybridResult -T 500 -i 100 -s 123457'
    if args.hnsig2:
        options = options + ' --significance --readHybridResult --toysFile=input.root --expectedFromGrid=0.5'
    if args.forkvar:
        options = options + ' --fork 4'
    if args.testStat:
        options = options + ' --testStat ' + args.testStat
    if args.level != None:
        options = options + ' --expectedFromGrid=%f'%(args.level)
    if args.fitStrategy:
        options = options + ' --minimizerStrategy %i'%(args.fitStrategy)
    if args.noSyst:
        options = options + ' --systematics 0'
    if args.freezeNuisances:
        options = options + ' --freezeNuisances ' + args.freezeNuisances
    if args.picky:
        options = options + ' --picky '
    if args.saveHybridResult:
        options = options + " -- saveHybridResult "
    if args.strictBounds:
        options = options + ' --strictBounds '
    if args.rAbsAcc:
        options = options + " --rAbsAcc " + str(args.rAbsAcc) + " "
    if args.rRelAcc:
        options = options + " --rRelAcc " + str(args.rRelAcc) + " "
    if args.verbose:
        options = options + ' -v' + str(args.verbose)
    if method != 'ProfileLikelihood' and method != 'MaxLikelihoodFit' and args.rMax == None and not args.noHint and not args.significance:
        options = options + ' --hintMethod ProfileLikelihood'
    if method == 'HybridNew' and args.toysH != None:
        options = options + ' --toysH %i'%(args.toysH)
    if args.toys:
        options = options + ' --toys ' + str(args.toys)
    if args.robustFit:
        options = options + ' --robustFit 1'
    if args.rPrior:
        # Take rMin and rMax from +/-2 sigma estimates from a previous iteration. This has to be done per-job, so do it later.
        pass
    else:
        if args.rMin != None:
            options = options + ' --rMin %f'%(args.rMin)
        if args.rMax != None:
            options = options + ' --rMax %f'%(args.rMax)
    if method == 'MarkovChainMC':
        options = options + ' --tries %i --proposal %s'%(args.tries, args.proposal)

    prefix = 'limits'
    if args.significance:
        prefix = 'significance'
    elif method == 'MaxLikelihoodFit':
        prefix = 'signal_xs'

    postfix = (('_' + args.postfix) if args.postfix != '' else '')
    if args.noSyst:
        postfix += "_noSyst"
    if args.freezeNuisances:
        postfix += "_" + args.freezeNuisances.replace(",", "_")

    datacards_path = limit_config.paths["datacards"]
    output_path = limit_config.paths["combine_logs"]
    condor_path = limit_config.paths["condor"]

    # change to the appropriate directory
    if args.condor:
        os.chdir(output_path)
    else:
        os.chdir(datacards_path)

    first = True

    for analysis in analyses:
        for model in models:
            for mass in masses[analysis]:

                logName = '%s_%s_%s_%s_m%i%s_%s.log'%(prefix, method, analysis, model, int(mass), postfix,args.fit_function)

                run_options = options + ' --name _%s_%s_m%i%s_%s --mass %i'%(analysis, model, int(mass),postfix,args.fit_function,int(mass))

                if args.rPrior:
                    run_options += " --rMin " + str(limit_config.limit_m2sigma_estimates[analysis][model][mass] / 5. * 100.)
                    run_options += " --rMax " + str(limit_config.limit_p2sigma_estimates[analysis][model][mass] * 5. * 100.)

                if args.condor:
                    cmd = "combine -M %s %s %s 2>&1 | tee %s"%(
                        method,
                        run_options,
                        os.path.basename(limit_config.get_datacard_filename(analysis, model, mass, args.fit_function, fitSignal=args.fitSignal)),
                        os.path.basename(os.path.join(('' if args.condor else output_path),logName)))
                else:
                    cmd = "combine -M %s %s %s 2>&1 | tee %s"%(
                        method,
                        run_options,
                        limit_config.get_datacard_filename(analysis, model, mass, args.fit_function, fitSignal=args.fitSignal),
                        os.path.join(('' if args.condor else output_path),logName))

                # if using Condor
                if args.condor:
                    
                    #submission_dir = limit_config.paths["condor"]
                    submission_dir = limit_config.paths["combine_logs"]
                    start_dir = os.getcwd()
                    os.chdir(submission_dir)

                    # create the executable script
                    bash_script_path = os.path.join(submission_dir,'run_%s_%s_%s_m%i%s.sh'%(method,analysis,model,int(mass),postfix))
                    #bash_content = bash_template
                    #bash_content = re.sub('DUMMY_CMD',cmd,bash_content)

                    bash_script = open(bash_script_path,'w')
                    bash_script.write(cmd)
                    bash_script.close()

                    # Create the condor command
                    condor_command = "csub " + bash_script_path
                    files_to_transfer = []
                    files_to_transfer.append(bash_script_path)
                    files_to_transfer.append(limit_config.get_datacard_filename(analysis, model, mass, args.fit_function, fitSignal=args.fitSignal))
                    files_to_transfer.append(limit_config.get_workspace_filename(analysis, model, mass, fitSignal=args.fitSignal))
                    condor_command += " -F " + ",".join(files_to_transfer)
                    condor_command += " -l combine_{}_\$\(Cluster\)_\$\(Process\).log".format(logName)
                    condor_command += " -s submit_combine_{}_{}_{}.jdl".format(analysis, model, mass)
                    condor_command += " -d " + submission_dir
                    condor_command += " --cmssw"
                    if not first or args.no_retar:
                        condor_command += " --no_retar "
                    print ">> Submitting job for %s %s resonance with m = %i GeV..."%(analysis,model, int(mass))
                    print "Submission command: "
                    print condor_command
                    os.system(condor_command)

                    print "---------------------------------------------------------------------------"
                    os.chdir(start_dir)
                else:
                    print ">> Running combine for %s %s resonance with m = %i GeV..."%(analysis, model, int(mass))
                    print "---------------------------------------------------------------------------"
                    print "Running: " + cmd + "\n"
                    os.system(cmd)
                if first:
                    first = False

if __name__ == '__main__':
    main()

