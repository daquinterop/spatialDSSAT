"""
This contains a set of handful funcitons to run the simulations.
"""

from shapely.geometry import Polygon
from shapely.prepared import prep
from math import ceil, floor
import numpy as np
import datetime
import pandas as pd

from DSSATTools import Management, TabularSubsection

# Functions to get the pixels for the polygon
def grid_bounds(geom, delta):
    minx, miny, maxx, maxy = geom.bounds
    minx = floor(minx/delta)*delta-.05; miny = floor(miny/delta)*delta-.05
    maxx = ceil(maxx/delta)*delta+.05; maxy = ceil(maxy/delta)*delta+.05
    nx = round((maxx - minx)/delta) 
    ny = round((maxy - miny)/delta)
    gx, gy = np.linspace(minx,maxx,nx+1), np.linspace(miny,maxy,ny+1)
    grid = []
    for i in range(len(gx)-1):
        for j in range(len(gy)-1):
            poly_ij = Polygon([[gx[i],gy[j]],[gx[i],gy[j+1]],[gx[i+1],gy[j+1]],[gx[i+1],gy[j]]])
            grid.append( poly_ij )
    return grid

def partition(geom, delta):
    prepared_geom = prep(geom)
    grid = list(filter(prepared_geom.intersects, grid_bounds(geom, delta)))
    return grid

def cust_round(val):
    if round(val, 1) == -0.:
        return 0.0
    else:
        return round(val, 1)
    
def pixel_coords(geom, delta=.1):
    """
    Returns the pixel centroids for the input geometry
    """
    grid = partition(geom, delta)
    centroids = map(lambda p: p.centroid, grid)

    return list(map(lambda p: (cust_round(p.x), cust_round(p.y)), centroids))

def init_management():
    """
    Creates a management instance and returns it. 
    """
    management = Management(
        planting_date = datetime.datetime(2021, 1, 1),
    )
    # Planting details from 
    # https://fintracu.fintrac.com/sites/default/files/tech_manuals/Fintrac%20U_Maize%20Cultivation%20Manual.pdf
    management.planting_details["PLRS"] = 90 # Row spacing
    management.planting_details["PPOP"] = 4.4 # Plant population
    management.planting_details["PPOE"] = 4 # Plant population at emergency
    # Fertilizer details
    management.simulation_controls["FERTI"] = "D" # Days after planting
    management.fertilizers["table"]["FMCD"] = "FE005" # Urea
    management.fertilizers["table"]["FACD"] = "AP002" # Broadcast, incorporated
    management.fertilizers["table"]["FDEP"] = 5 # 5 cm depth
    # Options
    management.simulation_controls["WATER"] = "Y" 
    management.simulation_controls["NITRO"] = "Y"
    management.simulation_controls["EVAPO"] = "F" # FAO ET
    management.simulation_controls["CO2"] = "D" # Default 380ppm
    management.simulation_controls["PHOTO"] = "C" # Canopy curve (daily) photosynthesis
    management.simulation_controls["MESEV"] = "R" # Ritchie-Ceres soil evaporation method
    management.simulation_controls["OVVEW"] = "N"
    management.simulation_controls["SUMRY"] = "N"
    management.simulation_controls["GROUT"] = "N"
    management.simulation_controls["CAOUT"] = "N"
    management.simulation_controls["WAOUT"] = "N"
    management.simulation_controls["NIOUT"] = "N"

    # Crop
    management._Management__cultivars["CR"] = "MZ"
    management.initial_conditions['PCR'] = "MZ"
    # Soil and Weather Station
    management.field["WSTA...."] = "WSTA"
    management.field["ID_SOIL"] = "ZW00000001"

    management.simulation_controls['SMODEL'] = "MZCER"
    return management


def management_ic(management, soil, country="Zimbabwe"):
    """
    Modify initial conditions in Management instance. soil the pixel coords (lon, lat).
    Set soil moisture to field capacity.
    """
    HOME = "/home/dquintero/dssat_service"
    soil_filename = f"{soil[0]:07.2f}_{soil[1]:07.2f}".replace(".", "p")+".SOL"
    with open(f"{HOME}/data/soil_data/dssat_soils/{soil_filename}", "r") as f:
        soil_lines = f.readlines()

    fieldCap = []
    for layer_line in soil_lines[6:]:
        if len(layer_line) > 2:
            fieldCap.append((layer_line[:6], layer_line[19:24]))
    fieldCap = np.array(fieldCap).astype(float)
    fieldCap = np.hstack([fieldCap, np.ones(fieldCap.shape)*.01]) # Default value for initial NO3 and NH4

    management.initial_conditions["table"] = TabularSubsection(pd.DataFrame(
        fieldCap,
        columns=['ICBL', 'SH2O', 'SNH4', 'SNO3'])
    )
    management.field['SLDP'] = fieldCap[:,0].max()
    management.field["...........XCRD"] = soil[0]
    management.field["...........YCRD"] = soil[1]
    return management

def management_dates(management, pdate, dap=0):
    """
    Modify all dates in management file. Basically sets IC, fertilizer and start date
    equal to planting date. 
    """
    management.planting_details["PDATE"] = pdate.strftime("%y%j")
    management.planting_details["EDATE"] = (pdate + datetime.timedelta(days=5)).strftime("%y%j")
    management.initial_conditions['ICDAT'] = pdate.strftime("%y%j") 
    # management.fertilizers["table"]["FDATE"] = (pdate + datetime.timedelta(days=dap)).strftime("%y%j") # 20 Days after planting
    management.simulation_controls['SDATE'] = pdate.strftime("%y%j")
    management.field["WSTA...."] = f"WSTA{str(pdate.year)[2:]}01"
    return management

def management_nitro(management, nitro):
    """
    Modifies nitrogen amount in the management instance
    
        nitro : tuple (dap, rate)
    """
    management.fertilizers["table"] = pd.concat([management.fertilizers["table"]]*len(nitro), ignore_index=True)
    for n, (dap, rate) in enumerate(nitro):
        management.fertilizers["table"].loc[n, "FDATE"] = dap
        management.fertilizers["table"].loc[n, "FAMN"] = rate
    return management

def management_cult(management, cult):
    """
    Modifies cultivar in the management instance
    """
    management._Management__cultivars["INGENO"] = cult
    return management