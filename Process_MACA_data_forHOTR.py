# -*- coding: utf-8 -*-
"""
Compile Multivariate Adaptive Constructed Analogs (MACA) data for analysis
https://climate.northwestknowledge.net/MACA/data_portal.php

@author: Michelle M. Fink, michelle.fink@colostate.edu
Colorado Natural Heritage Program, Colorado State University
Last Modified 08/05/2020 - Built on Python 3.7.7

** Notes about the input data for my own reference **
pr_*.nc tmmn_*.nc tmmx_*.nc NETCDF4
    time= 1827 (days since 1900-01-01, calendar: gregorian, 5 yrs per file)
    lat = 585 (degrees_north (Y), standard_name: latitude, cell center)
    lon = 1386 (degrees_east (X), standard_name: longitude, cell center)
    crs = 1 (WGS84, EPSG:4326)
variables:
    precipitation [lon,lat,time] (Chunking: [123,51,162])(Compression: shuffle,level 5)
            float
            _FillValue: -9999
            long_name: Precipitation
            units: mm
            grid_mapping: crs
            standard_name: precipitation
            cell_methods: time: sum(interval: 24 hours)
            comments: Total daily precipitation at surface; includes both liquid
            and solid phases from all types of clouds (both large-scale and convective)
            coordinates: time lon lat
    air_temperature [lon,lat,time] (Chunking: [123,51,162])(Compression: shuffle,level 5)
            float
            _FillValue: -9999
            long_name: Daily Maximum Near-Surface Air Temperature
            units: K
            grid_mapping: crs
            standard_name: air_temperature
            height: 2 m
            cell_methods: time: maximum(interval: 24 hours)
            coordinates: time lon lat
    air_temperature [lon,lat,time] (etc)
            units: K
            long_name: Daily Minimum Near-Surface Air Temperature
            cell_methods: time: minimum(interval: 24 hours)
            (etc)
Credit: John Abatzoglou & Katherine C. Hegewisch, University of Idaho.
    2014-05-15. License CC0 1.0
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
#%%#
import os
import numpy as np
from netCDF4 import Dataset
import nc_func_py3 as nc_func

session_dir = r"H:\Climate\Future"
outFolder = "Derived"
outdir = os.path.join(session_dir, outFolder)
scenarios = ["IPSL-CM5A-LR_r1i1p1_rcp45", "MIROC5_r1i1p1_rcp85"]
var_dict = {"pr":"precipitation", "tasmax":"air_temperature",
            "tasmin":"air_temperature"}
yrsets = ["2076_2080", "2081_2085", "2086_2090", "2091_2095", "2096_2099"]
cell = 0.041666667
inbbox = [-124.7666666, 49.4, -67.0583333, 25.06666667]
outbbox = [-112.558333169, 49.024999997, -96.058333532, 29.483333372]
clip_idx = nc_func.clipindex_fromXY((inbbox[0], inbbox[1]),
                                    (inbbox[2], inbbox[3]),
                                    (outbbox[0], outbbox[1]),
                                    (outbbox[2], outbbox[3]), cell)

filepfx = var_dict.keys()
daycnt = 1827

def extractvars(nc_name, clip, var):
    """Return a subsetted numpy masked array from the nc file """
    dset = Dataset(nc_name)
    dvar = dset.variables[var]
    ary = dvar[:, clip[3]:clip[1], clip[0]:clip[2]]
    if var == "air_temperature":
        ary = np.subtract(ary, 273.15) #Convert K to C
    dset.close()
    return ary

if not os.path.exists(outdir):
    os.makedirs(outdir)
#%%#
for scen in scenarios:
    for pfx in filepfx:
        print(scen, pfx)
        outnc = "_".join([pfx, scen, "MACAv2metdata", "2076_2099.nc"])
        outnc = os.path.join(outdir, outnc)
        if os.path.isfile(outnc):
            msg = 'The output ' + outnc + ' already exists.' + \
            '\nMoving on.\n'
            print(msg)
        else:
            #Compile the individual nc files
            fileset = "_".join(["macav2metdata", pfx, scen, "Z", "CONUS_daily.nc"])
            #NOTE NETCDF4 format can't use MFDataset for the compilation.
            #Here's a work-around.

            #Create a new nc with dimensions in the right order (T, Y, X).
            initnc = os.path.join(session_dir, fileset.replace("Z", yrsets[0]))
            ncd = Dataset(initnc)
            vVar = ncd.variables[var_dict[pfx]]
            #Clip to area of interest - note in this case Y increases North, so slice values
            #reversed from other climate data. At least won't need to flip eventual outputs.
            clipData = vVar[:, clip_idx[3]:clip_idx[1], clip_idx[0]:clip_idx[2]]
            vLon = ncd.variables["lon"]
            lnSlice = vLon[clip_idx[0]:clip_idx[2]]
            lnSlice = np.subtract(lnSlice, 360) #Correct for 0-360 notation
            vLat = ncd.variables["lat"]
            ltSlice = vLat[clip_idx[3]:clip_idx[1]]
            vTime = ncd.variables["time"]
            tData = vTime[:] #Get the full time period for now
            # FIXME: get a deprecation warning below
            tUnits = vTime.units
            tCalen = vTime.calendar
            ##
            ncd.close()
            if pfx == "pr":
                meta = "Total daily precipitation in mm."
            elif pfx == "tasmin":
                meta = "Mean daily minimum air temperature in degrees C."
                clipData = np.subtract(clipData, 273.15) #Convert K to C
            elif pfx == "tasmax":
                meta = "Mean daily maximum air temperature in degrees C."
                clipData = np.subtract(clipData, 273.15) #Convert K to C
            kwargs = {"varname": var_dict[pfx], "units": tUnits, "calendar": tCalen,
                      "metadatastr": meta}
            #Save as new netCDF with an unlimited time dimension
            nc_func.new_nc(clipData, tData, ltSlice, lnSlice, outnc, **kwargs)
            #Now we have to iterate over all the other nc files and append the data to outnc
            masterd = Dataset(outnc, "a") #previous function closed it so have to reopen
            for i in range(1, 5):
                yrset = yrsets[i]
                addnc = os.path.join(session_dir, fileset.replace("Z", yrset))
                addSlice = extractvars(addnc, clip_idx, var_dict[pfx])
                daycnt = addSlice.shape[0]
                mtime = masterd.variables["time"]
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
#%%#
