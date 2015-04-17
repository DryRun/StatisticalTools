#! /usr/bin/env python

import os
import sys
import argparse
from ROOT import TFile

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
parser.add_argument("-t", "--toys", type=int, dest="toys", default=10,
    help="number of toys",
    metavar="TOYS"
    )
parser.add_argument("--mu", action="store", type=float, dest="mu", default=-999.,
    help="MU",
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

k=0
for line in  ins:
  sample = os.path.basename(line)
  #sample = os.path.splitext(line)[0]
  sample = sample.split("_")[0]
  print ("process %s" % sample)
  line = line.rstrip('\n')
  sample = sample.rstrip('\n')
  mass = sample.split("Qstar")[1]
  print mass

  if not(os.path.isdir("args.output")):
    os.system("mkdir "+args.output)

  if not(args.mu==-999.):
    command = "combine -M MaxLikelihoodFit --minos all --expectSignal "+str(args.mu)+"  -n "+sample+"_MLfit  --saveWorkspace -t "+str(args.toys)+" --out "+args.output+" --saveToys --rMin -10 --rMax 100  "+line
  else:
    filename = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/dijet_limits/higgsCombineNormal.Asymptotic.mH'+mass+'.root'
    inf = TFile.Open(filename)
    tr = inf.Get('limit')
    y = []
    N = tr.GetEntriesFast()
    for i in xrange(N):
      tr.GetEntry(i)
      y.append(tr.limit)
      k+=1
   
    print "mu lim hi band : "+str(y[4])
    command = "combine -M MaxLikelihoodFit --minos all --expectSignal "+str(y[4])+"  -n "+sample+"_MLfit  --saveWorkspace -t "+str(args.toys)+" --out "+args.output+" --saveToys --rMin -10  --rMax 100 "+line

  
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
  print outputname 
  if args.run:
    os.system("bsub -q "+args.queue+" -o batch/"+logfile+" source "+pwd+"/"+outputname)

