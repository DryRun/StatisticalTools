#!/usr/bin/env python

import sys, os, copy, re
from ROOT import * 
from array import array

def main():

     inputFileWorkspace = TFile("datacards/workspace_qg_m1500.root") 

     workspace = inputFileWorkspace.Get("w")
     workspace.Print()

     mjj = workspace.var("mjj")
     mjj.Print()

     data = workspace.data("data_obs")
     signal = workspace.data("signal")
     signalUp = workspace.data("signal__JESUp")
     signalDown = workspace.data("signal__JESDown")
     data.Print()
     signal.Print()

     data_TH1_fineBinning = data.createHistogram("data_TH1_fineBinning",mjj)    
     signal_TH1_fineBinning = signal.createHistogram("signal_TH1_fineBinning",mjj)    
     signalUp_TH1_fineBinning = signalUp.createHistogram("signalUp_TH1_fineBinning",mjj)    
     signalDown_TH1_fineBinning = signalDown.createHistogram("signalDown_TH1_fineBinning",mjj)    
     print "data_TH1_fineBinning integral = ", data_TH1_fineBinning.Integral()
     print "signal_TH1_fineBinning integral = ", signal_TH1_fineBinning.Integral()

     #canvas = TCanvas()
     #data_TH1_fineBinning.Draw()    
     #signal_TH1_fineBinning.Draw()    
     #signalUp_TH1_fineBinning.Draw("sames")    
     #signalDown_TH1_fineBinning.Draw("sames")    
     #c1.SaveAs("c1.root")

if __name__ == '__main__':
    main()
