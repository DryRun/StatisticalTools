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
parser.add_argument("--inputListGen", type=str, dest="inputListGen", default="../lists/list_Qstar_genOnly.txt",
    help="list of datacards for the generation",
    )
parser.add_argument("--inputListFit", type=str, dest="inputListFit", default="../lists/list_Qstar.txt",
    help="list of datacards for the fit with flatParameter option",
    )
parser.add_argument("--inputFileLimits", type=str, dest="inputFileLimits", default="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/dijet_limits_interpolated_pseudodatasetDinko/higgsCombine'+sample+'_limit.Asymptotic.mH120.root",
    help="inputfile containing the limits (obtained with combine -m Asymptotic command)",
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
parser.add_argument("--tag", action="store", type=str, dest="tag", default="",
    help="some useful tag",
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
ins = open(args.inputListGen,"r")
insFit = open(args.inputListFit,"r")

############# read datacard for generation #####################
ii=0
k=0
datacards_fit = []

for line in  insFit:
  datacards_fit.append(line)

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
    os.system("mkdir -p "+args.output)
  
  
  ### when you want to inject 0 (or "fixed" signal strenght) signal you pass the signal strenght as an argument with --mu option 
  if not(args.mu==-999.):
    command = "combine -M GenerateOnly --expectSignal "+str(args.mu)+"  -n "+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+"  -t "+str(args.toys)+" --saveToys  "+line
    command_fit = "combine -M MaxLikelihoodFit  -n "+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+" -t "+str(args.toys)+" --out "+args.output+" --saveNormalizations  --toysFile "+args.output+"/higgsCombine"+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+".GenerateOnly.mH120.123456.root --rMin -2 --rMax 10 --minos none --initFromBonly  "+datacards_fit[ii]
  
  ### when you want to inject the 2-sigma signal you take the corresponding signal strenght from the limits output file 
  else:
    #filename = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_combine/src/CMSDIJET/StatisticalTools/dijet_limits_interpolated_pseudodatasetDinko/higgsCombine'+sample+'_limit.Asymptotic.mH120.root'
    filename = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_5/src/StatisticalTools/scripts/output_limits_dataRunD_830pb-1/higgsCombine'+sample+'_limit.Asymptotic.mH120.root'
    inf = TFile.Open(filename)
    tr = inf.Get('limit')
    y = []
    N = tr.GetEntriesFast()
    for i in xrange(N):
      tr.GetEntry(i)
      print tr.limit
      y.append(tr.limit)
      k+=1
  
    if (N>1): 
      twosigmalimit = y[4]
    else: 
      twosigmalimit = 0
    print "mu lim hi band : "+str(twosigmalimit)
    command = "combine -M GenerateOnly --expectSignal "+str(twosigmalimit)+"  -n "+sample+args.tag+"_MLfit_mu_limit  -t "+str(args.toys)+"  --saveToys  "+line
    command_fit = "combine -M MaxLikelihoodFit  -n "+sample+args.tag+"_MLfit_mu_limit  -t "+str(args.toys)+" --out "+args.output+" --saveNormalizations --toysFile "+args.output+"/higgsCombine"+sample+args.tag+"_MLfit_mu_limit.GenerateOnly.mH120.123456.root --rMin -2 --rMax 10  --minos none --initFromBonly  "+datacards_fit[ii]

  ii += 1  
  
  print "submit "+command
  print ""

  if not args.mu==-999:
    logfile = "logfile_"+sample+args.tag+"_mu"+str(int(args.mu))+".log"
    outputname = "batch/submit_"+sample+args.tag+"_mu"+str(int(args.mu))+".src"
  else:
    logfile = "logfile_"+sample+args.tag+"_mu_limit.log"
    outputname = "batch/submit_"+sample+args.tag+"_mu_limit.src"

  outputfile = open(outputname,'w')
  outputfile.write('#!/bin/bash\n')
  outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc481\n')
  outputfile.write('cd '+pwd+' \n')
  outputfile.write('eval `scramv1 runtime -sh`\n')
  outputfile.write(command+"\n")
  if not(args.mu==-999.):
    outputfile.write("mv higgsCombine"+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+".GenerateOnly.mH120.123456.root  "+args.output+"\n")  
    outputfile.write(command_fit+"\n")
    outputfile.write("mv higgsCombine"+sample+args.tag+"_MLfit_mu"+str(int(args.mu))+".MaxLikelihoodFit.mH120.123456.root  "+args.output+"\n")  
    
  else:
    outputfile.write("mv higgsCombine"+sample+args.tag+"_MLfit_mu_limit.GenerateOnly.mH120.123456.root  "+args.output+"/\n")  
    outputfile.write(command_fit+"\n")
    outputfile.write("mv higgsCombine"+sample+args.tag+"_MLfit_mu_limit.MaxLikelihoodFit.mH120.123456.root  "+args.output+"/\n")  

  print outputname 
  if args.run:
    os.system("bsub -q "+args.queue+" -o batch/"+logfile+" source "+pwd+"/"+outputname)

