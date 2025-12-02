from collections import OrderedDict as od
import math

pois = od()
pois["cG"] = {
  "factor":16*math.pi*math.pi,
  "multiplier":1e-5,
  "range":[-10,10],
  "nominal":0
}

pois["cA"] = {
  "factor":16*math.pi*math.pi,
  "multiplier":1e-4,
  "range":[-10,10],
  "nominal":0
}

pois["cHW"] = {
  "factor":1,
  "multiplier":1e-2,
  "range":[-12,16],
  "nominal":0
}


