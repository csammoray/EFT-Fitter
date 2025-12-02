import os, sys
import json
import re
from optparse import OptionParser
from collections import OrderedDict as od
from importlib import import_module
import pickle
import ROOT
import numpy as np

sys.path.append("./inputs")
sys.path.append("./functions")
sys.path.append("./params")

def get_options():
  parser = OptionParser()
  parser.add_option('--pois', dest='pois', default='params.HEL', help="Name of json file storing pois")
  parser.add_option('--poi', dest='poi', default='cG', help="POI to plot")
  parser.add_option('--functions', dest='functions', default='functions.HEL_STXS', help="Name of json file storing functions")
  parser.add_option('--inputs', dest='inputs', default='', help="Comma separated list of input files")
  parser.add_option("--translateBins", dest="translateBins", default=None, help="Translate STXS bins")
  return parser.parse_args()
(opt,args) = get_options()

# Functions for translations
def Translate(name, ndict):
    return ndict[name] if name in ndict else name
def LoadTranslations(jsonfilename):
    with open(jsonfilename) as jsonfile:
        return json.load(jsonfile)
translateBins = {} if opt.translateBins is None else LoadTranslations(opt.translateBins)

# Load parameters of interest

# Load functions

# Load input measurements

pois = __import__(opt.pois, globals(), locals(), ["pois"], 0)
functions = __import__(opt.functions, globals(), locals(), ["functions"], 0)
func_dict = {name: getattr(functions, name) for name in dir(functions) if not name.startswith("__")}
pois_dict = {name: getattr(pois, name) for name in dir(pois) if not name.startswith("__")}

pois = pois_dict["pois"]
functions = func_dict["functions"]

inputs = []

for i in opt.inputs.split(","):
  _cfg = __import__(i, globals(), locals(), ["name", "X", "rho"], 0)
  _input = od()
  _input['name'] = _cfg.name
  _input['X'] = _cfg.X
  _input['rho'] = _cfg.rho
  inputs.append(_input)

from tools.fitter import *

fit = fitter(pois,functions,inputs,False)

#stxs_bins = ['ttH']
# stxs_bins = ['ZH_lep_PTV_0_75','ZH_lep_PTV_75_150','ZH_lep_PTV_150_250_0J','ZH_lep_PTV_150_250_GE1J','ZH_lep_PTV_GT250','ZH_lep']
pT_bins = ['pTH_0_10', 'pTH_10_20', 'pTH_20_30', 'pTH_30_45', 'pTH_45_60', 'pTH_60_80', 'pTH_80_120', 'pTH_120_200', 'pTH_200_inf']

scaling = od()
for pT_bin in pT_bins:

  scaling[pT_bin] = od()
  for poi in pois.keys(): scaling[pT_bin][poi] = od()

  # Quadratic
  fit.setLinearOnly(False)  
  for poi in pois.keys(): 
    scaling[pT_bin][poi]['quad'] = od()
    c, mu = fit.scaling1D(poi,pT_bin,npoints=1000)
    scaling[pT_bin][poi]['quad']['c'] = c
    scaling[pT_bin][poi]['quad']['mu'] = mu

  # Linear
  fit.setLinearOnly()
  for poi in pois.keys():
    scaling[pT_bin][poi]['lin'] = od()
    c,mu = fit.scaling1D(poi,pT_bin,npoints=1000)
    scaling[pT_bin][poi]['lin']['c'] = c
    scaling[pT_bin][poi]['lin']['mu'] = mu

# Mage graphs
grs = od()
for pT_bin in pT_bins:
  for poi in pois.keys():
    grs['%s_vs_%s_quad'%(pT_bin,poi)] = ROOT.TGraph()
    grs['%s_vs_%s_lin'%(pT_bin,poi)] = ROOT.TGraph()
    for i in range(len(scaling[pT_bin][poi]['quad']['c'])): grs['%s_vs_%s_quad'%(pT_bin,poi)].SetPoint( grs['%s_vs_%s_quad'%(pT_bin,poi)].GetN(),scaling[pT_bin][poi]['quad']['c'][i], scaling[pT_bin][poi]['quad']['mu'][i] )
    for i in range(len(scaling[pT_bin][poi]['lin']['c'])): grs['%s_vs_%s_lin'%(pT_bin,poi)].SetPoint( grs['%s_vs_%s_lin'%(pT_bin,poi)].GetN(),scaling[pT_bin][poi]['lin']['c'][i], scaling[pT_bin][poi]['lin']['mu'][i] )

# Make plot
styleMap = od()
styleMap['quad'] = {'LineWidth':3,'LineStyle':1,'MarkerSize':0}
styleMap['quad_dummy'] = {'LineWidth':3,'LineStyle':1,'MarkerSize':0}
styleMap['lin'] = {'LineWidth':2, 'LineStyle':2,'MarkerSize':0}
styleMap['lin_dummy'] = {'LineColor':12, 'LineWidth':2, 'LineStyle':2,'MarkerSize':0}
#styleMap['lin_dummy'] = {'LineColor':ROOT.kMagenta-7, 'LineWidth':2, 'LineStyle':2,'MarkerSize':0}

colorMap = od()
#colorMap['ZH_lep'] = {'LineColor':ROOT.kRed-4,'MarkerColor':ROOT.kRed-4}
#colorMap['ZH_lep_PTV_0_75'] = {'LineColor':ROOT.kGreen-8,'MarkerColor':ROOT.kGreen-8}
#colorMap['ZH_lep_PTV_75_150'] = {'LineColor':ROOT.kGreen-7,'MarkerColor':ROOT.kGreen-7}
#colorMap['ZH_lep_PTV_150_250_0J'] = {'LineColor':ROOT.kGreen+1,'MarkerColor':ROOT.kGreen+1}
#colorMap['ZH_lep_PTV_150_250_GE1J'] = {'LineColor':ROOT.kGreen+3,'MarkerColor':ROOT.kGreen+3}
#colorMap['ZH_lep_PTV_GT250'] = {'LineColor':ROOT.kBlack,'MarkerColor':ROOT.kBlack}
#colorMap['ttH'] = {'LineColor':ROOT.kMagenta-7,'MarkerColor':ROOT.kMagenta-7}
colorMap['pTH_0_10'] = {'LineColor':ROOT.kRed-4,'MarkerColor':ROOT.kRed-4}
colorMap['pTH_10_20'] = {'LineColor':ROOT.kGreen-8,'MarkerColor':ROOT.kGreen-8}
colorMap['pTH_20_30'] = {'LineColor':ROOT.kGreen-7,'MarkerColor':ROOT.kGreen-7}
colorMap['pTH_30_45'] = {'LineColor':ROOT.kGreen+1,'MarkerColor':ROOT.kGreen+1}
colorMap['pTH_45_60'] = {'LineColor':ROOT.kGreen+3,'MarkerColor':ROOT.kGreen+3}
colorMap['pTH_60_80'] = {'LineColor':ROOT.kBlack,'MarkerColor':ROOT.kBlack}
colorMap['pTH_80_120'] = {'LineColor':ROOT.kBlue-4,'MarkerColor':ROOT.kBlue-4}
colorMap['pTH_120_200'] = {'LineColor':ROOT.kBlue-3,'MarkerColor':ROOT.kBlue-3}
colorMap['pTH_200_inf'] = {'LineColor':ROOT.kBlue-2,'MarkerColor':ROOT.kBlue-2}

# POI str
poi = opt.poi
hmax = 2.5

import math
#m = "%g"%math.log(1/pois[poi]['multiplier'],10)
#if m == '1': m = ''
m = ''
#if poi == "cWWMinuscB":
#  pstr_stripped = "c_{WW} #minus c_{B}"
#  pstr = "(c_{WW} #minus c_{B}) x 10^{%s}"%m
#else:
pstr_stripped = "c_{%s}"%poi.split("c")[-1]
pstr = "c_{%s} x 10^{%s}"%(poi.split("c")[-1],m)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

canv = ROOT.TCanvas("canv_%s"%poi,"canv_%s"%poi,700,500)
#canv = ROOT.TCanvas("canv_%s"%poi,"canv_%s"%poi,900,500)
canv.SetBottomMargin(0.15)
canv.SetTickx()
canv.SetTicky()

prange = pois[poi]['range'][1]-pois[poi]['range'][0]

h_axes = ROOT.TH1F("haxes","",100, pois[poi]['range'][0]-0.1*prange, pois[poi]['range'][1]+0.1*prange )
h_axes.SetMaximum(hmax)
h_axes.SetMinimum(-0.2)
h_axes.SetTitle("")
h_axes.GetXaxis().SetTitle(pstr)
h_axes.GetXaxis().SetTitleSize(0.05)
h_axes.GetXaxis().SetLabelSize(0.035)
h_axes.GetYaxis().SetTitle("#mu^{i}_{prod}(%s)"%pstr_stripped)
h_axes.GetYaxis().SetTitleSize(0.05)
h_axes.GetYaxis().SetTitleOffset(0.8)
h_axes.GetYaxis().SetLabelSize(0.035)
h_axes.GetYaxis().SetLabelOffset(0.007)
h_axes.GetYaxis().CenterTitle()
h_axes.SetLineWidth(0)
h_axes.Draw()

for pT_bin in pT_bins:
  for k, v in colorMap[pT_bin].items():
    getattr(grs["%s_vs_%s_quad"%(pT_bin,poi)],"Set%s"%k)(v)
    getattr(grs["%s_vs_%s_lin"%(pT_bin,poi)],"Set%s"%k)(v)
  for k, v in styleMap['quad'].items(): getattr(grs["%s_vs_%s_quad"%(pT_bin,poi)],"Set%s"%k)(v)
  for k, v in styleMap['lin'].items(): getattr(grs["%s_vs_%s_lin"%(pT_bin,poi)],"Set%s"%k)(v)
  grs["%s_vs_%s_quad"%(pT_bin,poi)].Draw("Same C")
  grs["%s_vs_%s_lin"%(pT_bin,poi)].Draw("Same C")

# Lines
hlines = {}
yvals = [0,1]
for i in range(len(yvals)):
  yval = yvals[i]
  hlines['hline_%g'%i] = ROOT.TLine(pois[poi]['range'][0]-0.1*prange,yval,pois[poi]['range'][1]+0.1*prange,yval)
  hlines['hline_%g'%i].SetLineColorAlpha(15,0.5)
  hlines['hline_%g'%i].SetLineStyle(2)
  hlines['hline_%g'%i].SetLineWidth(1)
  hlines['hline_%g'%i].Draw("SAME")

vlines = {}
xvals = [pois[poi]['range'][0],0,pois[poi]['range'][1]]
for i in range(len(xvals)):
  xval = xvals[i]
  vlines['vline_%g'%i] = ROOT.TLine(xval,-0.2,xval,hmax)
  vlines['vline_%g'%i].SetLineColorAlpha(15,0.5)
  vlines['vline_%g'%i].SetLineStyle(2)
  vlines['vline_%g'%i].SetLineWidth(1)
  vlines['vline_%g'%i].Draw("SAME")


# Text
lat0 = ROOT.TLatex()
lat0.SetTextFont(42)
lat0.SetTextAlign(11)
lat0.SetNDC()
lat0.SetTextSize(0.045)
#lat0.DrawLatex(0.1,0.92,"HEL UFO")

lat1 = ROOT.TLatex()
lat1.SetTextFont(42)
lat1.SetTextAlign(23)
lat1.SetTextSize(0.03)
xpos = pois[poi]['range'][0]-0.05*prange
lat1.DrawLatex(xpos,1.,"'#color[15]{#sigma = #sigma_{SM}}")
lat1.DrawLatex(xpos,0.,"#color[15]{#sigma = 0}")

lat2 = ROOT.TLatex()
lat2.SetTextFont(42)
lat2.SetTextAlign(23)
lat2.SetTextAngle(90)
lat2.SetTextSize(0.045)
lat2.SetTextAlign(21)
lat2.DrawLatex(pois[poi]['range'][0]-0.02*prange,0.9*hmax,"#color[15]{c_{min}}")
lat2.SetTextAlign(23)
lat2.DrawLatex(pois[poi]['range'][1]+0.01*prange,0.9*hmax,"#color[15]{c_{max}}")


# Legend

# Create dummy graph for linear
gr_lin_dummy = ROOT.TGraph()
for k,v in styleMap['lin_dummy'].items(): getattr(gr_lin_dummy,"Set%s"%k)(v)

leg = ROOT.TLegend(0.55,0.22,0.8,0.48)
#leg = ROOT.TLegend(0.63,0.28,0.8,0.38)
leg.SetFillStyle(0)
leg.SetLineColor(0)
leg.SetTextSize(0.0275)
#leg.SetTextSize(0.035)
for pT_bin in pT_bins:  leg.AddEntry( grs["%s_vs_%s_quad"%(pT_bin,poi)], Translate(pT_bin,translateBins), "L")
leg.AddEntry(gr_lin_dummy,"(Lin. terms only)","L")
leg.Draw("Same")

canv.Update()
canv.SaveAs("Scaling_per_bin_%s_wo_ggH.pdf"%poi)
#canv.SaveAs("/eos/home-j/jlangfor/www/CMS/thesis/chapter7/scaling_functions/ttH_vs_%s.png"%poi)
#canv.SaveAs("/eos/home-j/jlangfor/www/CMS/thesis/chapter7/scaling_functions/ttH_vs_%s.pdf"%poi)
