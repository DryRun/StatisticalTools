import os
import sys
import ROOT
from array import array
from ROOT import *
import math
from math import sqrt

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

# Helper functions
def calculate_rss(hist1, hist2):
    if hist1.GetNbinsX() != hist2.GetNbinsX():
        print "[calculate_rss] ERROR : Histogram bin mismatch"
        print "{}=>{}, {}=>{}".format(hist1.GetName(), hist1.GetNbinsX(), hist2.GetName(), hist2.GetNbinsX())
        sys.exit(1)
    rss = 0
    for bin in xrange(1, hist1.GetNbinsX() + 1):
        rss += (hist1.GetBinContent(bin) - hist2.GetBinContent(bin))**2
    return rss

def calculate_chi2(model_hist, data_hist):
    if model_hist.GetNbinsX() != data_hist.GetNbinsX():
        print "[calculate_rss] ERROR : Histogram bin mismatch"
        print "{}=>{}, {}=>{}".format(model_hist.GetName(), model_hist.GetNbinsX(), data_hist.GetName(), data_hist.GetNbinsX())
        sys.exit(1)
    chi2 = 0
    for bin in xrange(1, model_hist.GetNbinsX() + 1):
        if model_hist.GetBinContent(bin) > 0:
            fit = model_hist.GetBinContent(bin)
            data = data_hist.GetBinContent(bin)
            chi2 += ((fit - data) / sqrt(fit))**2
    return chi2

def calculate_andersondarling(model_hist, data_hist):
    ad_prob = data_hist.AndersonDarlingTest(model_hist)
    ad_ts = data_hist.AndersonDarlingTest(model_hist, "T")
    return (ad_prob, ad_ts)

def calculate_ks(model_hist, data_hist):
    ks_prob = data_hist.KolmogorovTest(model_hist)
    ks_ts = data_hist.KolmogorovTest(model_hist, "M")
    return (ks_prob, ks_ts)


if __name__ == "__main__":
    from argparse import ArgumentParser
    # input parameters
    parser = ArgumentParser(description="Make a plot from a combine MaxLikelihood run")
    parser.add_argument('--analyses', type=str, default="trigbbl_CSVTM,trigbbh_CSVTM", help="Analysis names")
    parser.add_argument('--models', type=str, default="Hbb", help="Model names")
    parser.add_argument('--masses', type=str, default="750", help="Mass points")
    parser.add_argument('--draw_signals', type=str, help="Signals to draw (format model.mass,model.mass ...")

    parser.add_argument('--qcd', action='store_true', help="Use QCD instead of data (assumes no trigger emulation)")
    parser.add_argument('--correctTrigger', action='store_true', help="Use model with trigger correction (has to have been specified in create_datacards.py)")
    parser.add_argument('--fitTrigger', action='store_true', help="Use model with trigger fit (has to have been specified in create_datacards.py)")
    parser.add_argument('--fitBonly', action='store_true', help="Background-only fit")
    parser.add_argument('--fitOffB', action='store_true', help="Fit background-only efficiency")
    parser.add_argument('--fit_function', type=str, default="dijet4", help="Name of central fit function")
    parser.add_argument("--sb", action="store_true", help="Draw S+B fit")
    parser.add_argument('--table', action='store_true', help="Print table of fit parameters")
    parser.add_argument('--qof', action='store_true', help="Print QOF information")
    args = parser.parse_args()

    signals_to_draw = []
    if args.draw_signals:
        signals_to_draw = [x.split(".") for x in args.draw_signals.split(",")]

    for model in args.models.split(","):
        for analysis in args.analyses.split(","):
            for mass_str in args.masses.split(","):
                mass = int(mass_str)

                # Data histogram
                workspace_file = TFile(limit_config.get_workspace_filename(analysis, model, mass, fitBonly=args.fitBonly, correctTrigger=args.correctTrigger, fitTrigger=args.fitTrigger, qcd=args.qcd, fitOffB=args.fitOffB), "READ")
                workspace = workspace_file.Get("w")
                mjj = workspace.var("mjj")
                data = workspace.data("data_obs").createHistogram("data", mjj, RooFit.Binning(5000, 0., 5000.))
                data.SetDirectory(0)
                workspace_file.Close()

                # Signal histograms
                if "bbl" in analysis:
                    fit_range = [296,1058]
                elif "bbh" in analysis:
                    fit_range = [526,1607]
                signal_histograms = {}
                signal_names = []
                signal_legends = {}
                signal_masses = {}
                for signal_pair in signals_to_draw:
                    signal_model = signal_pair[0]
                    signal_mass = int(signal_pair[1])
                    signal_name = signal_model + "_" + str(signal_mass)
                    signal_masses[signal_name] = signal_mass
                    workspace_file = TFile(limit_config.get_workspace_filename(analysis, signal_model, signal_mass, fitBonly=False, correctTrigger=args.correctTrigger, fitTrigger=args.fitTrigger, qcd=args.qcd, fitOffB=args.fitOffB), "READ")
                    workspace = workspace_file.Get("w")
                    mjj = workspace.var("mjj")
                    signal_names.append(signal_name)
                    signal_histograms[signal_name] = workspace.pdf("signal").createHistogram(signal_name, mjj, RooFit.Binning(fit_range[1]-fit_range[0],fit_range[0],fit_range[1]))
                    signal_histograms[signal_name].SetName("h_" + signal_name)
                    signal_histograms[signal_name].SetDirectory(0)
                    if signal_model == "Hbb":
                        signal_legends[signal_name] = "Scalar, "
                    elif signal_model == "RSG":
                        signal_legends[signal_name] = "RS graviton, "
                    elif signal_model == "ZPrime":
                        signal_legends[signal_name] = "Z', "
                    signal_legends[signal_name] += "m={} GeV".format(signal_mass)
                    workspace_file.Close()

                # Fit histogram
                postfix = ""
                if args.fitTrigger:
                    postfix += "_fitTrigger"
                if args.correctTrigger:
                    postfix += "_correctTrigger"
                if args.fitOffB:
                    postfix += "_fitOffB"
                if args.qcd:
                    postfix += "_qcd"
                job_name = "%s_%s_m%i%s_%s"%(analysis, model, int(mass), postfix, args.fit_function)
                fit_filename = "/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/Logs/plots_{}/mlfit_{}.root".format(job_name, job_name)
                print "Loading fit results from {}".format(fit_filename)
                fit_file = TFile(fit_filename, "READ")
                if args.sb:
                    background_fit_hist = fit_file.Get("shapes_fit_s/bin1/total_background")
                    signal_fit_hist = fit_file.Get("shapes_fit_s/bin1/total_signal")
                    total_fit_hist = signal_fit_hist.Clone()
                    total_fit_hist.Add(background_fit_hist)
                else:
                    background_fit_hist = fit_file.Get("shapes_fit_b/bin1/total_background")

                data_fitrange = background_fit_hist.Clone()
                data_fitrange.Reset()
                for bin in xrange(1, background_fit_hist.GetNbinsX() + 1):
                    data_bin = data.GetXaxis().FindBin(background_fit_hist.GetBinCenter(bin))
                    data_fitrange.SetBinContent(bin, data.GetBinContent(data_bin))
                    data_fitrange.SetBinError(bin, data.GetBinError(data_bin))
                # Normalize background to data?
                #print "Scaling background by {}".format(data_fitrange.Integral() / background_fit_hist.Integral())
                #background_fit_hist.Scale(data_fitrange.Integral() / background_fit_hist.Integral())
                print "KS prob = {}".format(background_fit_hist.KolmogorovTest(data_fitrange))
                print "KS d = {}".format(background_fit_hist.KolmogorovTest(data_fitrange, "M"))
                print "chi2 prob = {}".format(background_fit_hist.Chi2Test(data_fitrange))

                # Truncate and normalize signal histograms
                for signal_name in signal_names:
                    signal_histograms[signal_name].Scale(19700. / signal_histograms[signal_name].Integral())
                # Plot
                plotter = root_plots.DataMCPlot()
                if args.qcd:
                    data_name = "QCD MC"
                else:
                    data_name = "Data 2012"
                plotter.add_data(data, data_name)
                if args.sb:
                    plotter.add_mc(background_fit_hist, "S+B fit, B only", color=17)
                else:
                    plotter.add_mc(background_fit_hist, "Dijet Fit", color=17)
                if "bbl" in analysis or "eta1p7" in analysis:
                    x_range = [0., 1200.]
                elif "bbh" in analysis or "eta2p2" in analysis:
                    x_range = [0., 1800.]
                pull_dataerrors = args.qcd

                for i, signal_name in enumerate(signal_names):
                    color = seaborn.GetColorRoot("cubehelixlarge", int(math.floor(2*(signal_masses[signal_name] - 200) / 100.)), 25)
                    plotter.add_signal(signal_histograms[signal_name], signal_legends[signal_name], color=color, line_style=2+i, stacked=False)

                save_tag="mlfit_{}".format(job_name)
                if args.sb:
                    save_tag += "_sbfit_b"
                plotter.draw(logy=True, draw_pull=True, x_range=x_range, pull_range=[-4., 4.], save_tag=save_tag, cms_label="Internal", color_scheme="cubehelixhuge", complex_rebinning=dijet_binning, legend_position="topright", pull_dataerrors=pull_dataerrors, lumi_string="19.7 fb^{-1} (8 TeV)", x_title="m_{jj} [GeV]", y_title="Events / GeV", save_directory="/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Results/figures")

                if args.sb:
                    plotter_sb = root_plots.DataMCPlot()
                    if args.qcd:
                        data_name = "QCD MC"
                    else:
                        data_name = "Data 2012"
                    plotter_sb.add_data(data, data_name)
                    plotter_sb.add_mc(total_fit_hist, "S+B fit, total", color=17)
                    if "bbl" in analysis or "eta1p7" in analysis:
                        x_range = [0., 1200.]
                    elif "bbh" in analysis or "eta2p2" in analysis:
                        x_range = [0., 1800.]
                    pull_dataerrors = args.qcd
                    save_tag="mlfit_{}".format(job_name)
                    save_tag += "_sbfit_sb"
                    plotter_sb.draw(logy=True, draw_pull=True, x_range=x_range, pull_range=[-4., 4.], save_tag=save_tag, cms_label="Internal", color_scheme="cubehelixhuge", complex_rebinning=dijet_binning, legend_position="topright", pull_dataerrors=pull_dataerrors, lumi_string="19.7 fb^{-1} (8 TeV)", x_title="m_{jj} [GeV]", y_title="Events / GeV", save_directory="/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Results/figures")

                if args.table:
                    print "\\begin{table}"
                    print "\t\\begin{tabular}{|c|c|}"
                    print "\t\t\\hline"
                    print "\t\tParameter & Value \\\\\t\t\\hline"
                    fit_result = fit_file.Get("fit_b")
                    fitted_parameters = fit_result.floatParsFinal()
                    for i in xrange(fitted_parameters.getSize()):
                        par = fitted_parameters[i]
                        name = par.GetName()
                        if not args.fit_function in name:
                            continue
                        if "shapeBkg" in name:
                            continue
                        print "\t\t{}\t&\t${:.2f}\pm{:.2f}$\t\\\\\t\t\\hline".format(name, par.getVal(), par.getError())

                if args.qof:
                    print "rss = {}".format(calculate_rss(background_fit_hist, data_fitrange))

                    chi2 = calculate_chi2(background_fit_hist, data_fitrange)
                    print "chi2 = {}".format(chi2)
                    ndf = data_fitrange.GetNbinsX() - 4
                    print "ndf = {}".format(ndf)
                    print "chi2/ndf prob = {}".format(TMath.Prob(chi2, ndf))

                    print "andersondarling = {}".format(calculate_andersondarling(background_fit_hist, data_fitrange))

                    print "ks = {}".format(calculate_ks(background_fit_hist, data_fitrange))
