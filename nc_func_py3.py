# -*- coding: utf-8 -*-
"""
Helper functions for analyzing netCDF files (climate models)

@author: Michelle M. Fink, michelle.fink@colostate.edu
         Colorado Natural Heritage Program, Colorado State University
Created on Oct 31, 2017 Last updated 07/22/2019 - Built on Python 3.7.3

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
import time
import numpy as np
from osgeo import osr, gdal
from netCDF4 import Dataset

gdal.UseExceptions()

def calc_it(data, masks, method):
    """Run basic calculations on a masked array
    Right now method can be either 'sum' or 'mean' along axis=0
    data = ndarray
    masks = masks used to group data by years, seasons, month, etc.
    method = string; 'sum' or 'mean'

    output: ndarray reduced by calculation to a time dimension of 1.
    """
    if method == 'sum':
        return np.ma.array([np.ma.sum(np.ma.take(data, aMask, axis=0), axis=0,
                                      dtype=np.float) for aMask in masks])
    elif method == 'mean':
        return np.ma.array([np.ma.mean(np.ma.take(data, aMask, axis=0), axis=0,
                                       dtype=np.float) for aMask in masks])

def reverse(array, axis=0):
    """Flip a 2D array derived from a netCDF upside down (y axis).
    Necessary if the nc does not treat increasing y's as North,
    which happens more often than you might think.
    """
    idx = [slice(None)]*len(array.shape)
    idx[axis] = slice(None, None, -1)
    return array[idx]

def raster2array(ras_name, flip=True):
    """Convert a geoTIFF to a 2D numpy array.
    Defaults to flipping axis 0 (y axis) on assumption that original
    was flipped.
    """
    ras_tif = gdal.Open(ras_name)
    ras_band = ras_tif.GetRasterBand(1)
    ras_data = ras_band.ReadAsArray()
    ras_tif = None
    if flip:
        ras_ary = reverse(ras_data)
    else:
        ras_ary = ras_data
    return ras_ary

def array2raster(new_ras_file, array=None, *arglist):
    """Convert a 2D numpy array to a geoTIFF
    new_ras_file = string; full path and name of output tif
    array = the input array
    arglist = list of additional information
        [raster_origin, pixel_width, pixel_height]
        raster_origin = tuple of minimum X and maximum Y
        pixel_width = horizontal size of each pixel in the projected units
        pixel_height = vertical size of each pixel in the projected units
        NOTE - for decimal degrees in W. hemisphere, pixel_height is negative!

    The defaults provided are very specific to Colorado and the datasets I use,
    so it is very likely you do not want to use them.

    output: geoTIFF
    """
    if arglist == ():
        ras_origin = [-109.375, 41.25]
        pwidth = 0.125
        pheight = -0.125
    else:
        ras_origin = arglist[0]
        pwidth = arglist[1]
        pheight = arglist[2]

    cols = array.shape[1]
    rows = array.shape[0]
    origin_x = ras_origin[0]
    origin_y = ras_origin[1]
    driver = gdal.GetDriverByName('GTiff')
    out_ras = driver.Create(new_ras_file, cols, rows, 1, gdal.GDT_Float32)
    out_ras.SetGeoTransform((origin_x, pwidth, 0, origin_y, 0, pheight))
    outband = out_ras.GetRasterBand(1)
    outband.WriteArray(array)
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(4326)
    out_ras.SetProjection(proj.ExportToWkt())
    print('Finished writing', new_ras_file)

def new_nc(array, timeslice, yslice, xslice, outname, **kwargs):
    """Save the given ndarray to a new netCDF file
    array = The ndarray containing variable of interest. Must have 3
      dimensions: time, y (usually latitude), and x (usually longitude).
    timeslice = the ndarray containing the temporal values to assign to 'array'
    yslice = the ndarray containing the vertical coordinates for 'array'
    xslice = the ndarray containing the horizontal coordinates for 'array'
    outname = the full path and name of the output nc file.
    kwargs = dictionary of additional information. Keys = varname, timename,
      xname, yname, units (of time), calendar, metadatastr (describe the dataset)

    default kwargs provided in 'argdic'.

    output: a NETCDF4_CLASSIC file
    """
    argdic = {'varname':'var', 'timename':'time', 'xname':'longitude',
              'yname':'latitude', 'units':'days since 1950-01-01 00:00:00',
              'calendar':'gregorian',
              'metadatastr':''}
    #FIXME: right now there is no error checking regarding shape of array
    # (or anything else for that matter). Assumes dimensions are [time, y, x]
    x = array.shape[2]
    y = array.shape[1]
    if kwargs != {}:
        if 'varname' in kwargs:
            argdic['varname'] = kwargs['varname']
        if 'timename' in kwargs:
            argdic['timename'] = kwargs['timename']
        if 'xname' in kwargs:
            argdic['xname'] = kwargs['xname']
        if 'yname' in kwargs:
            argdic['yname'] = kwargs['yname']
        if 'units' in kwargs:
            argdic['units'] = kwargs['units']
        if 'calendar' in kwargs:
            argdic['calendar'] = kwargs['calendar']
        if 'metadatastr' in kwargs:
            argdic['metadatastr'] = kwargs['metadatastr']
    out_ds = Dataset(outname, 'w', format='NETCDF4_CLASSIC')
    out_ds.createDimension(argdic['timename'], None)
    out_ds.createDimension(argdic['yname'], y)
    out_ds.createDimension(argdic['xname'], x)
    time_var = out_ds.createVariable(argdic['timename'], 'f8', (argdic['timename'],))
    lat_var = out_ds.createVariable(argdic['yname'], 'f4', (argdic['yname'],))
    lon_var = out_ds.createVariable(argdic['xname'], 'f4', (argdic['xname'],))
    cvar = out_ds.createVariable(argdic['varname'], 'f4',
                                 (argdic['timename'], argdic['yname'],
                                  argdic['xname'],))
    time_var[:] = timeslice
    lon_var[:] = xslice
    lat_var[:] = yslice
    cvar[:] = array
    ds_att = {u'description': argdic['metadatastr'],
              u'history': 'Created ' + time.ctime(time.time())}
    time_att = {u'units': argdic['units'], u'calendar': argdic['calendar']}
    out_ds.setncatts(ds_att)
    time_var.setncatts(time_att)
    out_ds.sync()
    out_ds.close()
    print('Created', outname)

def clipindex_fromXY(full_uleft, full_lright, uleft, lright, stepx, stepy=None):
    """
    Gets the XY index values of a smaller area than a netCDF's full extent. Use
    to clip out a geographic subset of the data.
    Assumes a zero-based index.

    full_uleft = tuple; XY coordinates of full extent upper left corner
    full_lright = tuple; XY coordinates of full extent lower right corner
    uleft = tuple; XY coordinates of the desired smaller extent upper left
    lright = tuple; XY coordinates of the desired smaller extent lower right
    stepx = number; pixel size (width) in coordinate system units
    stepy = number; pixel size (height) in coordinate system units
        Only use stepy if width != height

    returns list of integer index numbers of the smaller extent.
    """
    if stepy == None:
        stepy = stepx

    clip_idx_x1 = int(round((uleft[0] - full_uleft[0]) / stepx))
    clip_idx_x2 = int(round(clip_idx_x1 + ((lright[0] - uleft[0]) / stepx)))
    clip_idx_y2 = int(round((lright[1] - full_lright[1]) / stepy))
    clip_idx_y1 = int(round(clip_idx_y2 + ((uleft[1] - lright[1]) / stepy)))

    return([clip_idx_x1, clip_idx_y1, clip_idx_x2, clip_idx_y2])

def nc2d_from_raster(ras_name, outname, **kwargs):
    """Create a new 2-dimensional netCDF file from a geotiff raster
    some code adapted from https://gis.stackexchange.com/questions/42790/
    gdal-and-python-how-to-get-coordinates-for-all-cells-having-a-specific-value#42846
    """
    argdic = {'varname':'var', 'xname':'longitude', 'yname':'latitude',
              'metadatastr':''}
    if kwargs != {}:
        if 'varname' in kwargs:
            argdic['varname'] = kwargs['varname']
        if 'xname' in kwargs:
            argdic['xname'] = kwargs['xname']
        if 'yname' in kwargs:
            argdic['yname'] = kwargs['yname']
        if 'metadatastr' in kwargs:
            argdic['metadatastr'] = kwargs['metadatastr']

    ras_tif = gdal.Open(ras_name)
    ras_band = ras_tif.GetRasterBand(1)
    ras_data = ras_band.ReadAsArray()
    (upper_left_x, x_size, x_rotation, upper_left_y, y_rotation, y_size) = ras_tif.GetGeoTransform()
    (y_index, x_index) = np.nonzero(ras_data)
    x_coords = x_index * x_size + upper_left_x + (x_size / 2)
    y_coords = y_index * y_size + upper_left_y + (y_size / 2)
    x = ras_data.shape[1]
    y = ras_data.shape[0]
    x_slice = x_coords.reshape(y, x)
    y_slice = y_coords.reshape(y, x)
    out_ds = Dataset(outname, 'w', format='NETCDF4_CLASSIC')
    out_ds.createDimension(argdic['yname'], y)
    out_ds.createDimension(argdic['xname'], x)
    lat_var = out_ds.createVariable(argdic['yname'], 'f4', (argdic['yname'],))
    lon_var = out_ds.createVariable(argdic['xname'], 'f4', (argdic['xname'],))
    cvar = out_ds.createVariable(argdic['varname'], 'f4', (argdic['yname'], argdic['xname'],))
    lon_var[:] = x_slice[0, :]
    lat_var[:] = y_slice[:, 0]
    cvar[:] = ras_data
    ds_att = {u'description': argdic['metadatastr'],
              u'history': 'Created ' + time.ctime(time.time())}
    out_ds.setncatts(ds_att)
    out_ds.sync()
    out_ds.close()
    print('Created ', outname)
