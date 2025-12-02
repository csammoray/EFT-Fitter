import ROOT
import pickle
import numpy as np
import math
from collections import OrderedDict as od
from scipy import interpolate
import sys
import argparse

def get_options():

  parser = argparse.ArgumentParser()
  parser.add_argument('--pois', dest='pois', default='params.SMEFT_pois', help="Name of json file storing pois")
  parser.add_argument('--inputPkls', dest='inputPkls', nargs='+', help="Input pkl files storing results")
  parser.add_argument('--labels', dest='labels', nargs='+', help="Labels for each pkl file")
  parser.add_argument('--outputDir', dest='outputDir', default='.', help="Output plot directory")
  parser.add_argument('--poi', dest='poi', default='', help="Main poi to plot")

  return parser.parse_args()

opt = get_options()

def loadPkl(pklfile):
  with open(pklfile,"rb") as fpkl: 
    results = pickle.load(fpkl)
  
  return results

def getPOISDict(pois):
  
  sys.path.append("./params")
  pois = __import__(pois, globals(), locals(), ["pois"], 0)
  pois_dict = {name: getattr(pois, name) for name in dir(pois) if not name.startswith("__")}
  # Dictionary with the pois and their info (factor, multiplier, range, title, nominal)
  pois = pois_dict["pois"]
   
  return pois

# Function to extract poi best-fit and +-1/2sigma points
def extractValsV2( _p, _dchi2 ):

  # Best-fit
  i_bf = _dchi2.argmin()
  bf = _p[i_bf]

  # Add union of intervals stuff: if i+1 and i-1 are > i && dchi2 <= 4 (minimum)
  i_min, minimum = [], []
  up1, up2, down1, down2 = [], [], [], []
  dchi2_min = []
  for i in range(1,len(_p)-1):
    if(_dchi2[i] < 4. )&( _dchi2[i-1] > _dchi2[i])&( _dchi2[i+1] > _dchi2[i]): i_min.append(i)

  # Loop over minima
  for i in i_min: 
    minimum.append( _p[i] )
    dchi2_min.append(_dchi2[i])
    
    # Find confidence intervals for each minimum
    found_1upsigma, found_2upsigma = False, False
    if _dchi2[i] < 1.:
      for j in range(i,len(_p)):
        if not found_1upsigma:
          if _dchi2[j] >= 1.:
            j_1up = j
            found_1upsigma = True
    for j in range(i,len(_p)):
      if not found_2upsigma:
        if _dchi2[j] >= 4.:
          j_2up = j
          found_2upsigma = True
    if found_1upsigma: up1.append( _p[j_1up] )
    else: up1.append(None)
    if found_2upsigma: up2.append( _p[j_2up] ) 
    else: up2.append(None)

    found_1downsigma, found_2downsigma = False, False
    if _dchi2[i] < 1.:
      for j in range(i,0,-1):
        if not found_1downsigma:
          if _dchi2[j] >= 1.:
            j_1down = j
            found_1downsigma = True
    for j in range(i,0,-1):
      if not found_2downsigma:
        if _dchi2[j] >= 4.:
          j_2down = j
          found_2downsigma = True
    if found_1downsigma: down1.append( _p[j_1down] )
    else: down1.append(None)
    if found_2downsigma: down2.append( _p[j_2down] ) 
    else: down2.append(None)

  return bf, np.array(minimum), np.array(dchi2_min), np.array(up1), np.array(up2), np.array(down1), np.array(down2)

# Function to extract poi best-fit and +-1/2sigma points
def extractVals( _p, _dchi2 ):

  # Best-fit
  i_bf = _dchi2.argmin()
  bf = _p[i_bf]

  # + confidence intervals
  found_1upsigma, found_2upsigma = False, False
  for i in range(i_bf,len(_p)):
    if not found_1upsigma:
      if _dchi2[i] >= 1.:
        i_1up = i
        found_1upsigma = True
  
    if not found_2upsigma:
      if _dchi2[i] >= 4.:
        i_2up = i
        found_2upsigma = True
  
  if found_1upsigma: up1 = abs(_p[i_1up]-bf)
  else: up1 = None
  if found_2upsigma: up2 = abs(_p[i_2up]-bf)
  else: up2 = None
  
  found_1downsigma, found_2downsigma = False, False
  for i in range(i_bf,0,-1):
    if not found_1downsigma:
      if _dchi2[i] >= 1.:
        i_1down = i
        found_1downsigma = True
  
    if not found_2downsigma:
      if _dchi2[i] >= 4.:
        i_2down = i
        found_2downsigma = True

  if found_1downsigma: down1 = -1*abs(_p[i_1down]-bf)
  else: down1 = None
  if found_2downsigma: down2 = -1*abs(_p[i_2down]-bf)
  else: down2 = None

  return bf, up1, up2, down1, down2
          
def makeGraphsForPOI(results_list, pois, poi, modes, labels):
    """
    Process chi-squared data and create ROOT graphs for a specific POI
    
    Args:
        results: Dictionary containing fit results data
        pois: Dictionary of POI definitions
        poi_name: String name of the POI to process
        modes: List of modes to process (default: ['fixed'])
    
    Returns:
        grs: OrderedDict containing the ROOT graphs
    """

    # Do interpolations and make dchi2 graphs
    grs = od()

    for idx, results in enumerate(results_list):
      label = labels[idx] if labels else f"dataset_{idx}"
      
      for poi_name in pois:
        # Only process POIs that exist in this dataset
        if poi_name not in results:
          continue
          
        for mode in modes:
          # Only process modes that exist for this POI
          if mode not in results[poi_name]:
            continue
            
          dchi2 = results[poi_name][mode]['dchi2']
          p = results[poi_name][mode]['pvals']

          if "linear" in mode: 
            f = interpolate.interp1d(p,dchi2)
          else: 
            f = interpolate.interp1d(p,dchi2,"cubic")

          pext = np.linspace( p.min(), p.max(), 10000 )
          dchi2ext = f(pext)
          results[poi_name][mode]['pvals_ext'] = pext
          results[poi_name][mode]['dchi2_ext'] = dchi2ext

          bf, up1, up2, down1, down2 = extractVals(pext,dchi2ext)
          results[poi_name][mode]['bestfit'] = bf
          results[poi_name][mode]['up01sigma'] = up1 
          results[poi_name][mode]['up02sigma'] = down1 
          results[poi_name][mode]['down01sigma'] = down1
          results[poi_name][mode]['down02sigma'] = down2

          bf, minimum, dchi2_min, up1, up2, down1, down2 = extractValsV2(pext,dchi2ext)
          results[poi_name][mode]['bestfitv2'] = bf
          results[poi_name][mode]['minimumv2'] = minimum
          results[poi_name][mode]['dchi2_minv2'] = dchi2_min
          results[poi_name][mode]['up01crossv2'] = up1
          results[poi_name][mode]['up02crossv2'] = up2
          results[poi_name][mode]['down01crossv2'] = down1
          results[poi_name][mode]['down02crossv2'] = down2

          if poi_name == poi: 
            grs[f"{poi}_{mode}_{idx}"] = ROOT.TGraph()
            grs[f"{poi}_{mode}_ext_{idx}"] = ROOT.TGraph()
            for i in range(len(p)): 
              grs[f"{poi}_{mode}_{idx}"].SetPoint(grs[f"{poi}_{mode}_{idx}"].GetN(), p[i], dchi2[i])
            for i in range(len(pext)): 
              grs[f"{poi}_{mode}_ext_{idx}"].SetPoint(grs[f"{poi}_{mode}_ext_{idx}"].GetN(), pext[i], dchi2ext[i])
    
    return grs
 
def createAxisLabel(pois, poi):
    """Create formatted axis label for the POI"""
    m = "%g"%math.log(1/pois[poi]['multiplier'],10)
    if m == '1' or m=='0': 
      pstr = "c_{%s}"%(poi.split("c")[-1])
    else:
      pstr = "c_{%s} x 10^{%s}"%(poi.split("c")[-1],m)

    if poi == "cWWMinuscB":
      pstr = "(c_{WW} #minus c_{B}) x 10^{%s}"%m
    return pstr


def setupCanvas(results_list, poi, mode):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    canv = ROOT.TCanvas()
    canv.SetBottomMargin(0.15)
    canv.SetTickx()
    canv.SetTicky()

    # Find global min/max across all datasets
    min_val = min([results[poi][mode]['pvals_ext'].min() for results in results_list])
    max_val = max([results[poi][mode]['pvals_ext'].max() for results in results_list])
    
    # Find the maximum chi-squared value across all datasets to set appropriate y-axis range
    max_chi2 = max([results[poi][mode]['dchi2_ext'].max() for results in results_list])
    y_max = min(max_chi2 * 1.2, 15.0)  # Add 20% padding, but cap at 15 if data is very large

    h_axes = ROOT.TH1F("haxes","", 100, min_val, max_val)
    h_axes.SetMaximum(y_max)
    h_axes.SetMinimum(0.)
    h_axes.SetTitle("")

    return canv, h_axes

def createStyleMap():
    styleMap = od()
    colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2]
    
    for i in range(4):
        marker_style = 20 if i == 0 else 24
        linewidth = 2 if i == 0 else 1

        styleMap[f'fixed_ext_{i}'] = {'LineColor': colors[i], 'LineWidth': linewidth, 'LineStyle': 1, 'MarkerSize': 0}
        styleMap[f'fixed_{i}'] = {'MarkerColor': colors[i], 'LineWidth': 0, 'MarkerSize': .75, 'MarkerStyle': marker_style}
        styleMap[f'fixed_dummy_{i}'] = {'LineColor': colors[i], 'MarkerColor': colors[i], 'LineWidth': linewidth, 'LineStyle': 1, 'MarkerSize': .75, 'MarkerStyle': marker_style}

    return styleMap

def drawGraphs(grs, results_list, h_axes, poi, mode):
    styleMap = createStyleMap()

    for idx in range(len(results_list)):
        # Apply styles to data points graph
        for k, v in styleMap[f"{mode}_{idx}"].items(): 
            getattr(grs[f"{poi}_{mode}_{idx}"], f"Set{k}")(v)
        
        # Apply styles to interpolated curve graph
        for k, v in styleMap[f"{mode}_ext_{idx}"].items():
            getattr(grs[f"{poi}_{mode}_ext_{idx}"], f"Set{k}")(v)

        # Draw the graphs
        grs[f"{poi}_{mode}_{idx}"].Draw("Same P")
        grs[f"{poi}_{mode}_ext_{idx}"].Draw("Same C")

    # Lines
    hlines = {}
    yvals = [1,4]
    for i in range(len(yvals)):
        yval = yvals[i]
        min_val = min([results[poi][mode]['pvals_ext'].min() for results in results_list])
        max_val = max([results[poi][mode]['pvals_ext'].max() for results in results_list])
        hlines[f'hline_{i}'] = ROOT.TLine(min_val, yval, max_val, yval)
        hlines[f'hline_{i}'].SetLineColorAlpha(ROOT.kRed, 0.5)
        hlines[f'hline_{i}'].SetLineStyle(2)
        hlines[f'hline_{i}'].SetLineWidth(1)
        hlines[f'hline_{i}'].Draw("SAME")

    # Removed the white box that was covering the data

def createLegend(labels):

    colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2]  # Same as createStyleMap
    dummy_graphs = []
    
    leg = ROOT.TLegend(0.15, 0.65, 0.55, 0.89)  # Back to original position
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetTextSize(0.032)

    for idx, label in enumerate(labels):
        # Use same marker style logic as in createStyleMap
        marker_style = 20 if idx == 0 else 24
        
        dummy_graph = ROOT.TGraph(1)
        dummy_graph.SetPoint(0, 0, 0)
        dummy_graph.SetLineColor(colors[idx])
        dummy_graph.SetLineWidth(2)
        dummy_graph.SetLineStyle(1)
        dummy_graph.SetMarkerColor(colors[idx])
        dummy_graph.SetMarkerStyle(marker_style)
        dummy_graph.SetMarkerSize(0.75)
        
        leg.AddEntry(dummy_graph, label, "LP")
        dummy_graphs.append(dummy_graph)

    # keep Python refs alive as long as the legend exists
    if not hasattr(leg, "_keep"):
        leg._keep = []
    leg._keep.extend(dummy_graphs)

    leg.Draw("same")
    return leg

def makePlot(grs, results_list, pois, poi, mode, labels):
    
    # Get axis label
    pstr = createAxisLabel(pois, poi)

    # Setup canvas
    canv, h_axes = setupCanvas(results_list, poi, mode)
    
    # Configure axes
    h_axes.GetXaxis().SetTitle(pstr)
    h_axes.GetXaxis().SetTitleSize(0.05)
    h_axes.GetXaxis().SetLabelSize(0.035)

    h_axes.GetYaxis().SetTitle("#Delta#chi^{2}")
    h_axes.GetYaxis().SetTitleSize(0.05)
    h_axes.GetYaxis().SetTitleOffset(0.8)
    h_axes.GetYaxis().SetLabelSize(0.035)
    h_axes.GetYaxis().SetLabelOffset(0.007)
    h_axes.GetYaxis().CenterTitle()
    h_axes.Draw()

    lat1 = ROOT.TLatex()
    lat1.SetTextFont(42)
    lat1.SetTextAlign(12)
    lat1.SetTextSize(0.035)

    # Draw graphs and lines
    drawGraphs(grs, results_list, h_axes, poi, mode)

    # Create legend
    legend = createLegend(labels=labels)

    canv.Update()
    canv.SaveAs(f"{opt.outputDir}/{poi}_comparison.pdf")
    canv.SaveAs(f"{opt.outputDir}/{poi}_comparison.root")

if __name__ == "__main__":
    # Load data
    pois = getPOISDict(pois=opt.pois)
    
    # Load all pkl files
    results_list = []
    for pkl_file in opt.inputPkls:
        results_list.append(loadPkl(pkl_file))
    
    # Create labels if not provided
    if opt.labels:
        labels = opt.labels
    else:
        # Use filename without extension as label
        labels = [pkl_file.split('/')[-1].replace('.pkl', '') for pkl_file in opt.inputPkls]
    
    # Ensure we have the right number of labels
    if len(labels) != len(opt.inputPkls):
        print(f"Warning: Number of labels ({len(labels)}) doesn't match number of pkl files ({len(opt.inputPkls)})")
        labels = [f"Dataset {i+1}" for i in range(len(opt.inputPkls))]
    
    # Process and create graphs
    grs = makeGraphsForPOI(results_list=results_list, pois=pois, poi=opt.poi, modes=['fixed'], labels=labels)

    # Make the plot
    makePlot(grs=grs, results_list=results_list, pois=pois, poi=opt.poi, mode='fixed', labels=labels)