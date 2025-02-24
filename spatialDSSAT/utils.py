"""
This contains a set of handful funcitons to run the simulations.
"""

from shapely.geometry import Polygon, Point
from shapely.prepared import prep
from math import ceil, floor
import numpy as np
from datetime import timedelta, datetime
import pandas as pd

import rasterio as rio
from netCDF4 import Dataset

from itertools import product
from tqdm import tqdm
import shutil
from typing import Union


from DSSATTools import Weather

# Functions to get the pixels for the polygon
def grid_bounds(geom, geotransform):
    """
    Given a defined geometry and a delta x, it returns a grid (series of polygons)
    that includes all the geometry. Grid is defined according to the geotransform.
    geotransform is defined as it's shown in https://gdal.org/tutorials/geotransforms_tut.html
    """
    delta = geotransform[1]
    minx, miny, maxx, maxy = geom.bounds
    minx = floor(minx/delta)*delta-.5*delta; miny = floor(miny/delta)*delta-.5*delta
    maxx = ceil(maxx/delta)*delta+.5*delta; maxy = ceil(maxy/delta)*delta+.5*delta
    nx = round((maxx - minx)/delta) 
    ny = round((maxy - miny)/delta)
    gx, gy = np.linspace(minx,maxx,nx+1), np.linspace(miny,maxy,ny+1)
    grid = []
    for i in range(len(gx)-1):
        for j in range(len(gy)-1):
            poly_ij = Polygon([[gx[i],gy[j]],[gx[i],gy[j+1]],[gx[i+1],gy[j+1]],[gx[i+1],gy[j]]])
            grid.append( poly_ij )
    return grid

def partition(geom, geotransform):
    """"
    It returns a the cells of a grid of delta size that touches that geometry.
    """
    prepared_geom = prep(geom)
    grid = list(filter(prepared_geom.intersects, grid_bounds(geom, geotransform)))
    return grid

def cust_round(val):
    val = round(val, 2)
    if val == -0.:
        return 0.0
    else:
        return val
    
def pixel_coords(geom, geotransform):
    """
    Given a geom it retuns the centroid of the 
    """
    grid = partition(geom, geotransform)
    centroids = map(lambda p: p.centroid, grid)
    return list(map(lambda p: (cust_round(p.x), cust_round(p.y)), centroids))

VARIABLE_MAP_AGERA5 = {
    'Wind_Speed_10m_Mean': "WIND", 'Dew_Point_Temperature_2m_Mean': "DEWP", 
    'Temperature_Air_2m_Max_24h': "TMAX", 'Temperature_Air_2m_Min_24h': "TMIN", 
    'Precipitation_Flux': 'RAIN', 'Solar_Radiation_Flux': "SRAD"
}
VARIABLE_TRANS_AGERA5 = {
    "TMAX": lambda x: x - 273.15, "TMIN": lambda x: x - 273.15,
    "TDEW": lambda x: x - 273.15, "SRAD": lambda x: x/1e6,
    "RAIN": lambda x: x, "WIND": lambda x: x
}

MINIMUM_VARIABLE_SET = ["TMAX", "TMIN", "SRAD", "RAIN"]
TEMP_VARS = ["TMAX", "TMIN", "TDEW"]

def weather_from_netcdf(
        nc_file:str,  elev:Union[str, float] ,geom:Polygon, save_to:str, 
        pars:dict=VARIABLE_MAP_AGERA5, trans:dict=VARIABLE_TRANS_AGERA5):
    """
    Creates .WTH files for the pixels within a polygon from a AgERA5 netCDF file.

    Arguments
    ----------
    nc_file: str
        Path to the netCDF file.
    elev: str ot float
        Path to elevation raster or constant elevation valued to assume.
    geom: shapely.Polygon
        The polygon from which the data will be extracted.
    save_to: str
        Path to the folder where the .WTH will be saved. The folder must exist.
    pars: dict
        A dictionary mapping the netCDF file variables to DSSAT weather variables.
        If not provided, the variable names from AgERA5 dataset is be used. For reference:
        https://cds.climate.copernicus.eu/cdsapp#!/dataset/10.24381/cds.6c68c9bb?tab=overview
    trans: dict
        A dictionary mapping each DSSAT weather variable to a function that transform
        the units of netcdf variable to the units required by DSSAT. For example, 
        if temperature is in Kelvin in the netcdf, then a function that transform 
        from Kelvin to Celsius would be needed. If not provided then AgERA5 is assumed 
        and default transformations are used.
    """
    assert all(map(lambda x: x in pars.values(), MINIMUM_VARIABLE_SET)), \
        f"netCDF file must contain at least {', '.join(MINIMUM_VARIABLE_SET)}"
    
    ds = Dataset(nc_file)
    for varname in pars.keys():
        assert varname in ds.variables.keys(), \
            f"Variable '{varname}' not in netCDF file"

    time = [datetime(1900, 1, 1) + timedelta(days=i) for i in ds.variables["time"][:]]

    x = ds.variables["lon"][:]
    y = ds.variables["lat"][:]

    # Buffer ensures that at least one pixel will be selected
    geom = geom.buffer(np.float32(x[1] - x[0])*.5)

    all_data = np.array([
        ds.variables[var][:] for var in pars.keys()
    ]) # dimensions: variable, time, lat, lon

    generator_list = zip(product(range(len(y)), range(len(x))), product(y, x))
    generator_list = list(filter(lambda x: geom.contains(Point(x[1][::-1])), generator_list))
    
    if isinstance(elev, str):
        elev = rio.open(elev)
        elev_list = list(elev.sample(list(map(lambda x: x[1][::-1], generator_list))))
        elev_list = list(map(int, elev_list))
    else:
        elev_list = [elev]*len(generator_list)
    n = 0
    for (j, i), (lat, lon) in tqdm(generator_list):
        lat = cust_round(lat)
        lon = cust_round(lon)

        df = pd.DataFrame(
            all_data[:, :, j, i].T, columns=pars.keys(),
            index=time
        )
        df = df.rename(columns=pars)
        if df.RAIN.mean() < 0: # Maybe no records in the sea
            continue
        # Some have TMIN slightly higher than TMAX
        df["TMAX"] = np.where(df["TMAX"] > df["TMIN"], df["TMAX"], df["TMIN"] +.1)
        # Apply transformations
        for var in df.columns:
            df[var] = df[var].map(trans[var])
        assert df.SRAD.min() > 0
        df["RAIN"] = df["RAIN"].abs()
        wth = Weather(
            df=df, pars=dict(zip(df.columns, df.columns)), 
            lat=lat, lon=lon, elev=elev_list[n]
        )
        wth.REFHT = 2
        wth.WNDHT = 10
        wth._name = f"{lon:07.2f}_{lat:07.2f}_{wth._name[4:]}".replace(".", "p")
        wth.write(save_to)
        n+=1
