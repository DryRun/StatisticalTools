#! /usr/bin/env python

import os
import sys
import argparse

usage = "python submit_fits_batch.py -q <queue> -i <input list> -o <output dir> -t <num of toys> --run"

parser = argparse.ArgumentParser(description='Process options.')

parser.add_argument('-q', '--queue', type=str, action='store',     dest='queue',       
    help='run in batch in queue specified as option (default -q cmslong)', 
    default='cmslong',
    metavar="QUEUE"
    )
parser.add_argument("-i", "--inputList", type=str, dest="inputList", default="../lists/list_Qstar.txt",
    help="directory containing the lists of datacards to be processed",
    )
parser.add_argument("-o", "--output", type=str, dest="output", default="./",
    help="the directory OUTDIR contains the output of the program",
    metavar="OUTDIR"
    )

parser.add_argument("--run", action="store_true", dest="run", default=False,
    help="run in batch",
    )



args = parser.parse_args()
print args
#####################################

os.system("mkdir -p batch")
pwd = os.environ['PWD']

#print(args.queue)
ins = open(args.inputList,"r")

for line in  ins:
  sample = os.path.basename(line)
  #sample = os.path.splitext(line)[0]
  sample = sample.split("_")[0]
  print ("process %s" % sample)
  line = line.rstrip('\n')
  sample = sample.rstrip('\n')

  #if not(os.path.isdir("args.output")):
  #  os.system("mkdir "+args.output)

  command = "combine -M Asymptotic -n "+sample+"_limit   "+line
  print "submit "+command
  print ""

  logfile = "logfile_"+sample+".log"
  outputname = "batch/submit_"+sample+".src"
  outputfile = open(outputname,'w')
  outputfile.write('#!/bin/bash\n')
  outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc481\n')
  outputfile.write('cd '+pwd+' \n')
  outputfile.write('eval `scramv1 runtime -sh`\n')
  outputfile.write(command+"\n")
  outputfile.write("mv higgsCombine"+sample+"_limit.Asymptotic.mH120.root "+args.output+"\n")

  
  print outputname 
  if args.run:
    os.system("bsub -q "+args.queue+" -o batch/"+logfile+" source "+pwd+"/"+outputname)
