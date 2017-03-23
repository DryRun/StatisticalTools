import os
import sys
import ROOT
from ROOT import *

ftest_dir = ""

def mlfit(model):
	directory = "{}/{}".format(ftest_dir, model)
	os.system("mkdir -pv {}".format(directory))
	