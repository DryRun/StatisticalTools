import os
import sys
import ROOT
from array import array
from ROOT import *

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/src/CMSDIJET/QCDAnalysis/python/")
import analysis_configuration_8TeV as analysis_config

sys.path.append("/uscms/home/dryu/Dijets/CMSSW_7_4_15/python/MyTools/RootUtils")
import histogram_tools

gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/CanvasHelpers.h\"")
ROOT.gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gSystem.Load("~/Dijets/CMSSW_7_4_15/lib/slc6_amd64_gcc491/libMyToolsRootUtils.so")
seaborn = Root.SeabornInterface()
seaborn.Initialize()

import CMSDIJET.StatisticalTools.limit_configuration as limit_config
import CMSDIJET.StatisticalTools.trigger_efficiency as trigger_efficiency
import MyTools.RootUtils.root_plots as root_plots
dijet_binning = array("d", [0, 1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5000]) #5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8000])


if __name__ == "__main__":
    from argparse import ArgumentParser
    # input parameters
    parser = ArgumentParser(description="Make a plot from a combine MaxLikelihood run")
    parser.add_argument('--analyses', type=str, default="NoTrigger_eta1p7_CSVTM,NoTrigger_eta2p2_CSVTM", help="Analysis names")
    parser.add_argument('--models', type=str, default="Hbb,RSG", help="Model names")
    parser.add_argument('--masses', type=str, default="750", help="Mass points")

    parser.add_argument('--qcd', action='store_true', help="Use QCD instead of data (assumes no trigger emulation)")
    parser.add_argument('--correctTrigger', action='store_true', help="Use model with trigger correction (has to have been specified in create_datacards.py)")
    parser.add_argument('--fitTrigger', action='store_true', help="Use model with trigger fit (has to have been specified in create_datacards.py)")
    parser.add_argument('--fitBonly', action='store_true', help="Background-only fit")
    parser.add_argument('--fit_function', type=str, default="f1", help="Name of central fit function")
    args = parser.parse_args()

    for model in args.models.split(","):
        for analysis in args.analyses.split(","):
            for mass_str in args.masses.split(","):
                mass = int(mass_str)

                # Data histogram
                workspace_file = TFile(limit_config.get_workspace_filename(analysis, model, mass, fitBonly=args.fitBonly, correctTrigger=args.correctTrigger, fitTrigger=args.fitTrigger, qcd=args.qcd), "READ")
                workspace = workspace_file.Get("w")
                mjj = workspace.var("mjj")
                data = workspace.data("data_obs").createHistogram("data", mjj, RooFit.Binning(5000, 0., 5000.))

                # Fit histogram
                postfix = ""
                if args.fitTrigger:
                    postfix += "_fitTrigger"
                if args.correctTrigger:
                    postfix += "_correctTrigger"
                if args.qcd:
                    postfix += "_qcd"
                job_name = "%s_%s_m%i%s_%s"%(analysis, model, int(mass), postfix, args.fit_function)
                fit_filename = "/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/Logs/plots_{}/mlfit_{}.root".format(job_name, job_name)
                fit_file = TFile(fit_filename, "READ")
                background_fit_hist = fit_file.Get("shapes_fit_b/bin1/total_background")

                data_fitrange = background_fit_hist.Clone()
                data_fitrange.Reset()
                for bin in xrange(1, background_fit_hist.GetNbinsX() + 1):
                    data_bin = data.GetXaxis().FindBin(background_fit_hist.GetBinCenter(bin))
                    data_fitrange.SetBinContent(bin, data.GetBinContent(data_bin))
                    data_fitrange.SetBinError(bin, data.GetBinError(data_bin))
                print "KS prob = {}".format(background_fit_hist.KolmogorovTest(data_fitrange))
                print "KS d = {}".format(background_fit_hist.KolmogorovTest(data_fitrange, "M"))
                print "chi2 prob = {}".format(background_fit_hist.Chi2Test(data_fitrange))

                # Plot
                plotter = root_plots.DataMCPlot()
                if args.qcd:
                    data_name = "QCD MC"
                else:
                    data_name = "Data 2012"
                plotter.add_data(data, data_name)
                plotter.add_mc(background_fit_hist, "Dijet Fit", color=17)
                if "bbl" in analysis or "eta1p7" in analysis:
                    x_range = [0., 1200.]
                elif "bbh" in analysis or "eta2p2" in analysis:
                    x_range = [0., 1800.]
                pull_dataerrors = args.qcd

                plotter.draw(logy=True, draw_pull=True, x_range=x_range, pull_range=[-4., 4.], save_tag="mlfit_{}".format(job_name), cms_label="Internal", color_scheme="cubehelixhuge", complex_rebinning=dijet_binning, legend_position="topright", pull_dataerrors=pull_dataerrors)
