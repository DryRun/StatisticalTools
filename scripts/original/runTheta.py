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
    usage = "Example: ./scripts/runTheta.py -d datacards -f gg --massrange 1300 5500 100"

    # input parameters
    parser = ArgumentParser(description='Script that runs limit calculation for specified mass points',epilog=usage)

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

    parser.add_argument("-t", "--theta_path", dest="theta_path",
                        default='../theta/utils2/',
                        help="Path to theta-auto.py script (default: %(default)s)",
                        metavar="THETA_PATH")

    parser.add_argument("--noSyst", dest="noSyst", default=False, action="store_true", help="Run without systematic uncertainties")

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

    datacards_path = os.path.join(os.getcwd(), args.datacards_path)
    output_path = os.path.join(os.getcwd(), args.output_path)
    theta_path = os.path.join(os.getcwd(), args.theta_path)
    condor_path = os.path.join(output_path, 'condor')

    # change to the appropriate directory
    if args.condor:
        os.chdir(output_path)
    else:
        os.chdir(datacards_path)

    for mass in masses:

        fileId = '%s_m%i%s'%(args.final_state,int(mass),('_NoSyst' if args.noSyst else ''))
        jobId = '%s_m%i'%(args.final_state,int(mass))

        logName = 'theta_%s.log'%(jobId)

        cmd = "%s theta_%s.py | tee %s"%(os.path.join(theta_path,'theta-auto.py'),fileId,os.path.join(('' if args.condor else output_path),logName))

        postProcessing = [
            'cat theta_%s_observed_limit.txt'%(fileId),
            'cat theta_%s_observed_limit.txt >> %s'%(fileId,os.path.join(('' if args.condor else output_path),logName)),
            'rm theta_%s_observed_limit.txt'%(fileId),
            'rm -rf theta_%s/'%(fileId)
        ]

        # if using Condor
        if args.condor:
            # check if the Condor directory exists
            if not os.path.isdir( condor_path ):
                os.mkdir( condor_path )

            # create the jdl file
            jdl_content = jdl_template
            jdl_content = re.sub('DUMMY_JOB','theta_%s'%(jobId),jdl_content)
            jdl_content = re.sub('DUMMY_OUTPUTDIR',os.getcwd(),jdl_content)
            files_to_transfer = []
            files_to_transfer.append( os.path.join(datacards_path, 'theta_%s.py'%(fileId)) )
            files_to_transfer.append( os.path.join(datacards_path, 'theta_%s.root'%(fileId)) )
            jdl_content = re.sub('DUMMY_FILES',', '.join(files_to_transfer),jdl_content)

            jdl_file = open(os.path.join(condor_path,'run_theta_%s.jdl'%(jobId)),'w')
            jdl_file.write(jdl_content)
            jdl_file.close()

            # create the Bash script
            bash_content = bash_template
            bash_content = re.sub('DUMMY_CMD',cmd + '\n' + '\n'.join(postProcessing),bash_content)

            bash_script = open(os.path.join(condor_path,'run_theta_%s.sh'%(jobId)),'w')
            bash_script.write(bash_content)
            bash_script.close()

            print ">> Submitting job for %s resonance with m = %i GeV..."%(args.final_state, int(mass))
            print "---------------------------------------------------------------------------"
            condor_cmd = 'condor_submit condor/run_theta_%s.jdl'%(jobId)
            print "Running: " + condor_cmd + "\n"
            os.system(condor_cmd)
        else:
            print ">> Running theta for %s resonance with m = %i GeV..."%(args.final_state, int(mass))
            print "---------------------------------------------------------------------------"
            print "Running: " + cmd + "\n"
            os.system(cmd)
            for ppCmd in postProcessing:
                os.system(ppCmd)

if __name__ == '__main__':
    main()

