# -*- coding: utf-8 -*-
"""
Calculate climate metris from gridMET data
http://www.climatologylab.org/wget-gridmet.html
NOTE that 'Process_Met_data_forHOTR.py' must have already been run.

@author: Michelle M. Fink, michelle.fink@colostate.edu
Colorado Natural Heritage Program, Colorado State University
Created 07/25/2019, Last Modified 07/30/2019 - Built on Python 3.7.3

Code licensed under the GNU General Public License version 3.
This script is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see https://www.gnu.org/licenses/
"""
#%%#
import os
import numpy as np
from netCDF4 import Dataset, num2date
import nc_func_py3 as nc_func
import watyrcalcs as wy

session_dir = r"E:\Climate\metdata\Derived"
var_dict = {"pr":"precipitation_amount", "tmmn":"air_temperature",
            "tmmx":"air_temperature", "pet":"potential_evapotranspiration"}
filepattern = "($var)_gridmet_1994_2014.nc"
source = "gridMET_1994_2014"
hYrs = range(1994, 2015, 1)
cell = 0.041666667
outbbox = [-112.579166412354, 49.0038690567017, #geotiffs use cell centers
           -96.0801734924316, 29.4217309951782] #which is what these numbers are
arglist = [(outbbox[0], outbbox[1]), cell, (cell * -1)]

def iname(prefix, pattern):
    """ yeah, this is a lame function """
    name = pattern.replace("($var)", prefix)
    return(name)

def calc_seasonal(var, season, nc_name, calc_type, years, time="time"):
    """ Get simple seasonal summaries from netCDF file. Use of water-year is assumed.
        var: string; name of data variable
        season: integer; 0=annual, 1=winter, 2=spring, 3=summer, 4=autumn,
                         5=winter-spring, 6=summer-fall
        nc_name: string; full path and name of netCDF file
        calc_type: string; "sum" or "mean"
        years: integer range of 4-digit years to calculate over
        time: string; name of the time dimension, defaults to "time"

    returns: 2D numpy array
    """
    ds = Dataset(nc_name, "r")
    dvar = ds.variables[var]
    ddata = dvar[:]
    tvar = ds.variables[time]
    tdata = tvar[:]
    tUnits = tvar.units
    tCalen = tvar.calendar
    ds.close()
    dateRng = num2date(tdata, tUnits, tCalen)
    mydates = tdata.tolist()
    mwargs = {"units": tUnits, "calen": tCalen}
    hmask = wy.watyrmask(dateRng, mydates, years, season, **mwargs)
    calcary = nc_func.calc_it(ddata, hmask, calc_type)
    outary = np.mean(calcary, axis=0)
    return(outary)

def iter_blocks(ary, years):
    """ iterator for 1-year chunk sizes of daily data
        ary: 3D numpy array with shape (time, y, x)
        years: numpy array replacing the 'time' value with a year value with shape (year, , )

    yields a generator object
    """
    rows = ary.shape[1]
    cols = ary.shape[2]
    depth = ary.shape[0]
    unique_yrs, indices = np.unique(years, return_index=True)

    for i in range(0, len(unique_yrs)):
        y = unique_yrs[i]
        if y % 4 == 0:
            daycnt = 366 #leap years
        else:
            daycnt = 365
        if indices[i] == 0:
            start_t = 0
            end_t = daycnt
        else:
            start_t = end_t
            end_t = end_t + daycnt
        yield ary[start_t:end_t, 0:rows, 0:cols]

def daily_tmean(tmin_nc, tmax_nc, tmin_var, tmax_var):
    """ Calculate the mean daily temperature from daily min and max temperatures
        tmin_nc: string; full path and name of min temp nc
        tmax_nc: string; full path and name of max temp nc
        tmin_var: string: name of the min temp variable
        tmax_var: string: name of the max temp variable

    Nothing returned. Creates a new netCDF file of mean daily temperatures
    """
    tnds = Dataset(tmin_nc, "r")
    txds = Dataset(tmax_nc, "r")
    tnvar = tnds.variables[tmin_var]
    txvar = txds.variables[tmax_var]
    txData = txvar[:]
    #FIXME: I've got the below variables hard-coded for now
    vLon = txds.variables["longitude"]
    lnSlice = vLon[:]
    vLat = txds.variables["latitude"]
    ltSlice = vLat[:]
    vTime = txds.variables["time"]
    tdata = vTime[:]
    tUnits = vTime.units
    tCalen = vTime.calendar
    txds.close()
    tnData = tnvar[:]
    tnds.close()
    dateRng = num2date(tdata, tUnits, tCalen)
    yrs = np.array([a_date.year for a_date in dateRng], dtype=np.int)
    progress = 0
    blocks_max = iter_blocks(txData, yrs)
    blocks_min = iter_blocks(tnData, yrs)

    for bx, bn in zip(blocks_max, blocks_min):
        progress += 1
        print("Processing chunk", progress)
        mnary = np.ma.divide(np.ma.add(bx, bn), 2)
        if progress == 1:
            tmean = np.ma.array(mnary)
        else:
            tmean = np.ma.concatenate((tmean, mnary), axis=0)

    meta = "Mean daily air temperature in degrees Celcius, from gridded surface meteorological data."
    kwargs = {"varname": "air_temperature", "units": tUnits, "calendar": tCalen,
              "metadatastr": meta}
    outnc = tmax_nc.replace("tmmx", "tmean")
    nc_func.new_nc(tmean, tdata, ltSlice, lnSlice, outnc, **kwargs)

def growingdegrees(nc_name, var, base, time="time"):
    """ Calculate growing degree-days from daily mean temperature netCDF
    NOTE that this is based on calendar year, whereas all other metrics are water year,
    because of the year-based chunks. I don't think it makes a significant difference for
    this particular metric. Also note that I am taking a 'short-cut' by basing GDD on the
    daily mean temperature instead of the more precise min and max. In this use-case, GDD
    is a proxy for forage availability, so this should be adequate.
        nc_name: string; full path and name of netCDF with daily values for multiple years
        var: string; name of data variable
        base: number; the base temperature to use, in same units as var
        time: string; name of the time dimension, defaults to "time"

    returns: 2D numpy array of the mean growing degree-days over all years
    """
    ds = Dataset(nc_name, "r")
    dvar = ds.variables[var]
    ddata = dvar[:]
    tvar = ds.variables[time]
    tdata = tvar[:]
    tUnits = tvar.units
    tCalen = tvar.calendar
    ds.close()
    dateRng = num2date(tdata, tUnits, tCalen)
    yrs = np.array([a_date.year for a_date in dateRng], dtype=np.int)
    progress = 0
    blocks = iter_blocks(ddata, yrs)

    for chunk in blocks:
        progress += 1
        print("Processing chunk", progress)
        condition = eval("chunk < " + str(base))
        sumary = np.ma.array([np.ma.sum(np.ma.where(condition, 0, chunk - base), axis=0)])
        if progress == 1:
            chkary = np.ma.array(sumary)
        else:
            chkary = np.ma.concatenate((chkary, sumary), axis=0)
    outary = np.mean(chkary, axis=0, dtype=np.float)
    return(outary)

#%%#

#Calculate growing degree-days for base0 and base5
print("Starting with the hardest first; GDD")
tn = "tmmn"
tx = "tmmx"
tnfile = os.path.join(session_dir, iname(tn, filepattern))
txfile = os.path.join(session_dir, iname(tx, filepattern))
tmean = txfile.replace("tmmx", "tmean")
if os.path.isfile(tmean):
    print("Found", tmean, "- Using this one.")
else:
    daily_tmean(tnfile, txfile, var_dict[tn], var_dict[tx])

base = [0, 5]
pfx = "GDD"
for b in base:
    gdd = "_".join([pfx, "base" + str(b)])
    outfile = os.path.join(session_dir, wy.clean_name(source, 0, gdd) + ".tif")
    outary = growingdegrees(tmean, var_dict[tn], b)
    flipary = nc_func.reverse(outary)
    nc_func.array2raster(outfile, flipary, *arglist)
#%%#
#Calculate total precipitation for annual, winter-spring, and summer-fall
print("On to precipitation!")
pfx = "pr"
infile = os.path.join(session_dir, iname(pfx, filepattern))
seasons = [0, 5, 6]
for i in seasons:
    outfile = os.path.join(session_dir, wy.clean_name(source, i, pfx) + ".tif")
    tifary = calc_seasonal(var_dict[pfx], i, infile, "sum", hYrs)
    flipary = nc_func.reverse(tifary)
    nc_func.array2raster(outfile, flipary, *arglist)
#%%#
#Calculate max summer temperature
print("Max temperature next")
pfx = "tmmx"
infile = os.path.join(session_dir, iname(pfx, filepattern))
outfile = os.path.join(session_dir, wy.clean_name(source, 3, pfx) + ".tif")
tifary = calc_seasonal(var_dict[pfx], 3, infile, "mean", hYrs)
flipary = nc_func.reverse(tifary)
nc_func.array2raster(outfile, flipary, *arglist)
#%%#
#Calculate average annual PET
print("Finally, PET")
pfx = "pet"
infile = os.path.join(session_dir, iname(pfx, filepattern))
outfile = os.path.join(session_dir, wy.clean_name(source, 0, pfx) + ".tif")
tifary = calc_seasonal(var_dict[pfx], 0, infile, "mean", hYrs)
flipary = nc_func.reverse(tifary)
nc_func.array2raster(outfile, flipary, *arglist)
#%%#
