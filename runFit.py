import os, sys
import json
import re
from optparse import OptionParser
from collections import OrderedDict as od
from importlib import import_module
import pickle

def get_options():
  parser = OptionParser()
  parser.add_option('--pois', dest='pois', default='params.HEL', help="Name of json file storing pois")
  parser.add_option('--functions', dest='functions', default='functions.HEL_STXS', help="Name of json file storing functions")
  parser.add_option('--inputs', dest='inputs', default='', help="Comma separated list of input files")
  parser.add_option('--doAsimov', dest='doAsimov', default=False, action="store_true", help="Do asimov fit (i.e. set all best-fit to nominal)")
  parser.add_option('--doReset', dest='doReset', default=False, action="store_true", help="Reset poi values each step in profiled scan")
  parser.add_option('--doFlip', dest='doFlip', default=False, action="store_true", help="Start scan from max val of poi")
  parser.add_option('--doProfiled', dest='doProfiled', default=False, action="store_true", help="Perform profiled scans")
  parser.add_option('--outputname', dest='outputname', default='', help="Name of the output pickle file")
  return parser.parse_args()
(opt,args) = get_options()

# Load parameters of interest
# pois = import_module(opt.pois).pois

# Load functions
# functions = import_module(opt.functions).functions

# Load input measurements
inputs = []

# sys.path.append("/afs/cern.ch/work/m/mbonanom/EFT-Fitter/inputs")
# sys.path.append("/afs/cern.ch/work/m/mbonanom/EFT-Fitter/functions")
# sys.path.append("/afs/cern.ch/work/m/mbonanom/EFT-Fitter/params")
sys.path.append("./inputs")
sys.path.append("./functions")
sys.path.append("./params")

pois = __import__(opt.pois, globals(), locals(), ["pois"], 0)
functions = __import__(opt.functions, globals(), locals(), ["functions"], 0)
func_dict = {name: getattr(functions, name) for name in dir(functions) if not name.startswith("__")}
pois_dict = {name: getattr(pois, name) for name in dir(pois) if not name.startswith("__")}

for i in opt.inputs.split(","):
  _cfg = __import__(i, globals(), locals(), ["name", "X", "rho"], 0) 
  _input = od()
  _input['name'] = _cfg.name
  _input['X'] = _cfg.X
  _input['rho'] = _cfg.rho
  inputs.append(_input)

from tools.fitter import *
fit = fitter(pois_dict["pois"],func_dict["functions"],inputs,opt.doAsimov)

# Perform scans
results = od()
for poi in pois_dict["pois"].keys():

  # print " --> Running fits for: %s"%poi
  results[poi] = od()

  # Linear
  fit.setLinearOnly(True)

  # Fixed scan
  # print "    * Linear: fixed"
  pvals_f_lin, chi2_f_lin = fit.scan_fixed(poi,npoints=100)
  results[poi]["fixed_linear"] = od()
  results[poi]["fixed_linear"]['pvals'] = pvals_f_lin
  results[poi]["fixed_linear"]['chi2'] = chi2_f_lin
  results[poi]["fixed_linear"]['dchi2'] = chi2_f_lin-chi2_f_lin.min()

  if opt.doProfiled:
    # Profiled scan (full)
    # print "    * Linear: profiled"
    pvals_p_lin, chi2_p_lin, all_pvals_p_lin = fit.scan_profiled(poi,npoints=50,freezeOtherPOIS=[],resetEachStep=opt.doReset,reverseScan=opt.doFlip,verbose=True)
    results[poi]["profiled_linear"] = od()
    results[poi]["profiled_linear"]['pvals'] = pvals_p_lin
    results[poi]["profiled_linear"]['chi2'] = chi2_p_lin
    results[poi]["profiled_linear"]['allpvals'] = all_pvals_p_lin
    results[poi]["profiled_linear"]['dchi2'] = chi2_p_lin-chi2_p_lin.min()

  # Quadratic
  fit.setLinearOnly(False)

  # Fixed scan
  # print "    * Quadratic: fixed"
  pvals_f, chi2_f = fit.scan_fixed(poi,npoints=100)
  results[poi]["fixed"] = od()
  results[poi]["fixed"]['pvals'] = pvals_f
  results[poi]["fixed"]['chi2'] = chi2_f
  results[poi]["fixed"]['dchi2'] = chi2_f-chi2_f.min()

  if opt.doProfiled:
    # Profiled scan (full)
    # print "    * Quadratic: profiled"
    pvals_p, chi2_p, all_pvals_p = fit.scan_profiled(poi,npoints=50,freezeOtherPOIS=[],resetEachStep=opt.doReset,reverseScan=opt.doFlip,verbose=False)
    results[poi]["profiled"] = od()
    results[poi]["profiled"]['pvals'] = pvals_p
    results[poi]["profiled"]['chi2'] = chi2_p
    results[poi]["profiled"]['allpvals'] = all_pvals_p
    results[poi]["profiled"]['dchi2'] = chi2_p-chi2_p.min()

extStr = ""
if opt.doAsimov: extStr += "_asimov"
else: extStr += "_observed"
if opt.doReset: extStr += "_reset"
if opt.doFlip: extStr += "_flip"
if opt.doProfiled: extStr += "_profiled"

if opt.outputname == '':
    with open(f"results{extStr}.pkl", "wb") as fpkl: pickle.dump(results, fpkl)
else:
    with open(f'{opt.outputname}.pkl', "wb") as fpkl: pickle.dump(results, fpkl)

#from scipy import interpolate
#import matplotlib.pyplot as plt
#f = interpolate.interp1d(cg,dchi2)
#cg_new = np.linspace(-8,10,1000)
#dchi2_new = f(cg_new)
#plt.plot(cg,dchi2,'o',cg_new,dchi2_new,'-')

