# -*- coding: utf-8 -*-
"""
Download script for Polaris soils data tiles.

@author: Michelle M. Fink, michelle.fink@colostate.edu
Colorado Natural Heritage Program, Colorado State University

Created on Wed May 29 09:28:58 2019
Code last modified June 16, 2019

*IMPORTANT* This script requires Wget. On Windows, it needs to be run through
the MinGW shell (or equivalent) in order for the wget command to be recognized.
	https://www.gnu.org/software/wget/

Soils parameter units:
    bd - bulk density, g/cm3
    sand - sand percentage, %
    clay - clay percentage, %
    om - organic matter, log10(%)
    ph - soil pH in H2O, N/A

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

import subprocess
import os
import csv

OUTDIR = r"E:\Polaris"
DATA_URL = "http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0/$(param)/mean/" \
    + "$(depth)/lat$(lt)_lon$(ln).tif"

parameters = ["clay", "sand", "bd", "om", "ph"]
depths = ["0_5", "5_15", "15_30", "30_60", "60_100"]
# Not all latlon combinations will exist!
latblocks = ["2930", "3031", "3132", "3233", "3334", "3435", "3536", "3637", "3738",
             "3839", "3940", "4041", "4142", "4243", "4344", "4445", "4546", "4647",
             "4748", "4849", "4950"]
lonblocks = ["-96-95", "-97-96", "-98-97", "-99-98", "-100-99", "-101-100",
             "-102-101", "-103-102", "-104-103", "-105-104", "-106-105", "-107-106",
             "-108-107", "-109-108", "-110-109", "-111-110", "-112-111", "-113-112"]

downloaded = {}

for p in parameters:
    parurl = DATA_URL.replace("$(param)", p)
    files = []
    for d in depths:
        depurl = parurl.replace("$(depth)", d)
        for lt in latblocks:
            lturl = depurl.replace("$(lt)", lt)
            for ln in lonblocks:
                lnurl = lturl.replace("$(ln)", ln)
                inFile = lnurl.rsplit("/", 1)
                outFile = "_".join([p, d, inFile[1]])
                outFname = os.path.join(OUTDIR, outFile)
                if os.path.exists(outFname):
                    print "Already have " + outFile
                else:
                    print "Trying download of " + p + d + inFile[1]
                    try:
                        cmd = ["wget", "--passive-ftp", "--retr-symlinks",
                               "--limit-rate=10m", lnurl, "-O", outFname]
                        subprocess.call(cmd)
                        print "Finished downloading " + outFname
                        files.append(outFile)
                    except Exception, e:
                        print "Nope"
                        print e.message
    downloaded.update(p=files)

print "Finished with the downloading."
# Save the list of downloaded files
w = csv.writer(open(os.path.join(OUTDIR,"downloaded.csv"), "w"))
for key, val in downloaded.items():
    w.writerow([key, val])
