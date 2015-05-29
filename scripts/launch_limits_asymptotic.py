#!/usr/bin/python
import os
               
for i in range(1200,6100,100):
  cmd = "combine -M Asymptotic  --name Qstar"+str(i)+" ../datacards_pseudodatasetDinko_fitData/Qstar"+str(i)+"_datacard_1.txt "
  print cmd
  os.system(cmd)  
