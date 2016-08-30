import os
import sys
import ROOT
from ROOT import *
import array

from argparse import ArgumentParser
parser = ArgumentParser(description='Make combined SR limit plot')
parser.add_argument('--fit_function', type=str, default="f4", help='Fit function')
args = parser.parse_args()

gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadBottomMargin(0.15)
gROOT.ForceStyle()

models = ["Hbb", "RSG"]
x_ranges = {
	"trigbbl_CSVTM":[400, 600],
	"trigbbh_CSVTM":[600, 1200]
}

canvases = {}
frames = {}
for model in models:
	canvases[model] = TCanvas("c_limits_" + model + "_" + args.fit_function, "c_limits_" + model + "_" + args.fit_function, 800, 600)
	canvases[model].cd()
	canvases[model].SetLogy()
	l = TLegend(0.7, 0.75, 0.9, 0.9)
	l.SetFillColor(0)
	l.SetBorderSize(0)
	l.SetTextFont(42)
	l.SetTextSize(0.03)
	l.SetHeader('95% CL upper limits, ' + model)

	frames[model] = TH1D("frame" + model, "frame" + model, 100, 200., 1400.)
	frames[model].SetMinimum(0.1)
	frames[model].SetMaximum(50.)
	frames[model].GetXaxis().SetTitle("Resonance Mass [GeV]")
	frames[model].GetYaxis().SetTitle("#sigma#timesBR(b#bar{b}) [pb]")
	frames[model].Draw("axis")

	graphs = {}
	graphs_trimmed = {}
	for analysis in ["trigbbl_CSVTM", "trigbbh_CSVTM"]:
		graphs[analysis] = {}
		graphs_trimmed[analysis] = {}
		f_graphs = TFile("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/Limits/limits_" + analysis + "_" + model + "_" + args.fit_function + ".root", "READ")
		for graph_name in ["graph_exp_2sigma", "graph_exp_1sigma", "graph_exp", "graph_obs"]:
			graphs[analysis][graph_name] = f_graphs.Get(graph_name)

			x_values = []
			y_values = []
			for i in xrange(graphs[analysis][graph_name].GetN()):
				x = graphs[analysis][graph_name].GetX()[i]
				y = graphs[analysis][graph_name].GetY()[i]
				if x >= x_ranges[analysis][0] and x <= x_ranges[analysis][1]:
					x_values.append(x)
					y_values.append(y)
			graphs_trimmed[analysis][graph_name] = TGraph(len(x_values), array.array('d', x_values), array.array('d', y_values))
		f_graphs.Close()
		graphs_trimmed[analysis]["graph_exp_2sigma"].SetFillColor(kYellow)
		graphs_trimmed[analysis]["graph_exp_1sigma"].SetFillColor(kGreen+1)
		graphs[analysis]["graph_exp"].SetLineWidth(3)
		graphs[analysis]["graph_exp"].SetLineStyle(3)
		graphs[analysis]["graph_exp"].SetLineColor(4)
		graphs_trimmed[analysis]["graph_exp"].SetLineWidth(3)
		graphs_trimmed[analysis]["graph_exp"].SetLineStyle(2)
		graphs_trimmed[analysis]["graph_exp"].SetLineColor(4)
		graphs[analysis]["graph_obs"].SetMarkerStyle(24)
		graphs[analysis]["graph_obs"].SetMarkerColor(kGray)
		graphs[analysis]["graph_obs"].SetLineWidth(3)
		graphs[analysis]["graph_obs"].SetLineStyle(3)
		graphs[analysis]["graph_obs"].SetLineColor(kGray)
		graphs_trimmed[analysis]["graph_obs"].SetMarkerStyle(20)
		graphs_trimmed[analysis]["graph_obs"].SetLineWidth(3)
		graphs_trimmed[analysis]["graph_obs"].SetLineColor(1)
		graphs_trimmed[analysis]["graph_exp_2sigma"].Draw("F")
		graphs_trimmed[analysis]["graph_exp_1sigma"].Draw("F")
		graphs_trimmed[analysis]["graph_exp"].Draw("L")
	for analysis in ["trigbbl_CSVTM", "trigbbh_CSVTM"]:
		graphs_trimmed[analysis]["graph_obs"].Draw("LP")
	l.AddEntry(graphs_trimmed["trigbbl_CSVTM"]["graph_obs"],"Observed","lp")
	l.AddEntry(graphs_trimmed["trigbbl_CSVTM"]["graph_exp"],"Expected","lp")
	l.AddEntry(graphs_trimmed["trigbbl_CSVTM"]["graph_exp_1sigma"],"#pm 1#sigma","F")
	l.AddEntry(graphs_trimmed["trigbbl_CSVTM"]["graph_exp_2sigma"],"#pm 2#sigma","F")
	l.Draw()
	canvases[model].SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Results/figures/" + canvases[model].GetName() + ".pdf")
	for analysis in ["trigbbl_CSVTM", "trigbbh_CSVTM"]:
		graphs[analysis]["graph_exp"].Draw("L")
		graphs[analysis]["graph_obs"].Draw("LP")
		graphs_trimmed[analysis]["graph_exp"].Draw("L")
		graphs_trimmed[analysis]["graph_obs"].Draw("LP")
	canvases[model].SaveAs("/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Results/figures/" + canvases[model].GetName() + "_comparison.pdf")
