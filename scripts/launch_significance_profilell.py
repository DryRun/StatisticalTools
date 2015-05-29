#!/usr/bin/python
import os
               
for i in range(1200,6100,100):
  cmd = "combine -M ProfileLikelihood --signif --name Qstar"+str(i)+"_toy1 -d ../datacards_pseudodatasetDinko/Qstar"+str(i)+"_datacard_nuisance_1.txt "
  print cmd
  os.system(cmd)  
