# -*- coding: utf-8 -*-
"""
Compile gridMET data for analysis
http://www.climatologylab.org/wget-gridmet.html

@author: Michelle M. Fink, michelle.fink@colostate.edu
Colorado Natural Heritage Program, Colorado State University
Last Modified 07/30/2019 - Built on Python 3.7.3

** Notes about the input data for my reference **
pr_*.nc tmmn_*.nc tmmx_*.nc pet_*.nc NETCDF3_CLASSIC
    day = 365 (days since 1900-01-01, calendar: gregorian, standard_name: time)
    lat = 585 (degrees_north (Y), standard_name: latitude)
    lon = 1386 (degrees_east (X), standard_name: longitude)
    crs = 1 (WGS_1984, EPSG:4326, long_name: WGS 84)
variables:
    precipitation_amount [lon,lat,day] (float32)
            units: mm
            description: Daily Accumulated Precipitation
            coordinates: lon lat
            cell_methods: time: sum(interval: 24 hours)
            missing_value: -32767
            grid_mapping: crs
    air_temperature [lon,lat,day]  (float32)
            units: K
            description: Daily Minimum Temperature
            coordinates: lon lat
            cell_methods: time: minimum(interval: 24 hours)
            height: 2 m
            missing_value: -32767
            grid_mapping: crs
    air_temperature [lon,lat,day]  (float32)
            units: K
            description: Daily Maximum Temperature
            (etc)
    potential_evapotranspiration [lon,lat,day]
            units: mm
            description: Daily reference evapotranspiration (short grass)
            cell_methods: time: sum(interval: 24 hours)
Citation:   Abatzoglou, J.T., 2013, Development of gridded surface meteorological data for ecological
            applications and modeling, International Journal of Climatology, DOI: 10.1002/joc.3413
****************
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

import os
import numpy as np
from netCDF4 import Dataset
import nc_func_py3 as nc_func

session_dir = r"E:\Climate\metdata"
outFolder = "Derived"
outdir = os.path.join(session_dir, outFolder)
 var_dict = {"pr":"precipitation_amount", "tmmn":"air_temperature",
             "tmmx":"air_temperature", "pet":"potential_evapotranspiration"}
hYrs = range(1994, 2015, 1)
cell = 0.041666667
inbbox = [-124.7666666, 49.4, -67.0583333, 25.06666667]
outbbox = [-112.558333169, 49.024999997, -96.058333532, 29.483333372]
clip_idx = nc_func.clipindex_fromXY((inbbox[0], inbbox[1]),
                                    (inbbox[2], inbbox[3]),
                                    (outbbox[0], outbbox[1]),
                                    (outbbox[2], outbbox[3]), cell)

filepfx = var_dict.keys()

def extractvars(nc_name, clip, var):
    """Return a subsetted numpy masked array from the nc file """
    dset = Dataset(nc_name)
    dvar = dset.variables[var]
    ary = dvar[:, clip_idx[3]:clip_idx[1], clip_idx[0]:clip_idx[2]]
    if var == "air_temperature":
        ary = np.subtract(ary, 273.15) #Convert K to C
    dset.close()
    return(ary)

if not os.path.exists(outdir):
    os.makedirs(outdir)

#Compile the individual netCDFs
for pfx in filepfx:
    print("Starting on", pfx)
    fileset = pfx + "_Z.nc"
    outnc = "_".join([pfx, "gridmet", str(hYrs[0]), str(hYrs[-1])]) + ".nc"
    outnc = os.path.join(outdir, outnc)
    #NOTE that these nc files are unusually structured, so can't use MFDataset for the
    #compilation. Here's a work-around.

    #Create a new nc with dimensions in the right order (T, Y, X). Starting with the
    #year *before* our desired date range so that water-years can be calculated.
    initnc = os.path.join(session_dir, fileset.replace("Z", str(hYrs[0] - 1)))
    ncd = Dataset(initnc)
    vVar = ncd.variables[var_dict[pfx]]
    #Clip to area of interest - note in this case Y increases North, so slice values
    #reversed from other climate data. At least won't need to flip eventual outputs.
    clipData = vVar[:, clip_idx[3]:clip_idx[1], clip_idx[0]:clip_idx[2]]
    vLon = ncd.variables["lon"]
    lnSlice = vLon[clip_idx[0]:clip_idx[2]]
    vLat = ncd.variables["lat"]
    ltSlice = vLat[clip_idx[3]:clip_idx[1]]
    vTime = ncd.variables["day"]
    tData = vTime[:] #Want the full year
    tUnits = vTime.units
    tCalen = vTime.calendar
    ncd.close()
    if pfx == "pr":
        meta = "Total daily precipitation in mm, from gridded surface meteorological data."
    elif pfx == "tmmn":
        meta = "Mean daily minimum air temperature in degrees Celcius, from gridded surface meteorological data."
        clipData = np.subtract(clipData, 273.15) #Convert K to C
    elif pfx == "tmmx":
        meta = "Mean daily maximum air temperature in degrees Celcius, from gridded surface meteorological data."
        clipData = np.subtract(clipData, 273.15) #Convert K to C
    else:
        meta = "Daily reference evapotranspiration (short grass) in mm, from gridded surface meteorological data."
    kwargs = {"varname": var_dict[pfx], "units": tUnits, "calendar": tCalen,
            "metadatastr": meta}
    #Save as new netCDF with an unlimited time dimension
    nc_func.new_nc(clipData, tData, ltSlice, lnSlice, outnc, **kwargs)

    #Now we have to iterate over all the other nc files and append the data to outnc
    masterd = Dataset(outnc, "a") #previous function closed it so have to reopen
    for i in hYrs:
        if i % 4 == 0:
            daycnt = 366 #leap years
        else:
            daycnt = 365
        addnc = os.path.join(session_dir, fileset.replace("Z", str(i)))
        addSlice = extractvars(addnc, clip_idx, var_dict[pfx])
        mtime = masterd.variables["time"] #new_nc function changed the variable names
        lastime = mtime[-1]
        lentime = len(mtime)
        startidx = lentime - daycnt
        newtimes = np.ma.add(mtime[startidx:], + daycnt)
        data = masterd.variables[var_dict[pfx]]
        data[lentime:lentime + daycnt, :, :] = addSlice
        mtime[lentime:lentime + daycnt] = newtimes
        masterd.sync()
    masterd.close()
    print("Finished", outnc)
