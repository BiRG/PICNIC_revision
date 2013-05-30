#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import os
import csv
from splitBy import splitBy


def cnRecords(record):
    return record[0][0:2] == 'CN'


def partition_probes(filepath):
    """
Usage: Parameter 1 (required): path to a ASCII-encoded .feature_intensity file as output by PICNIC's CelConverter jar.

NB: Assumes first row is header.
"""
    garbage, intensities_filename = os.path.split(filepath)
    intensities_filename, ascii_intensities_file_file_ext = \
            os.path.splitext(intensities_filename)

    try:
        with open(filepath, 'rb') as ascii_intensities_file:
        
            allRecords = csv.reader(ascii_intensities_file, delimiter='\t')
            
            # Grab the header row and drop trailing empty string if exists
            headerSeg = allRecords.next()
            if headerSeg[len(headerSeg)-1] == '':
                headerSeg = headerSeg[0:len(headerSeg)-1]

            cnProbes, snpProbes = splitBy(cnRecords, allRecords)

    except IOError:
        sys.stderr.write("\nCannot open file " + filepath)
        
    cnProbesFilename = intensities_filename + '.CNProbes.csv'
    snpProbesFilename = intensities_filename + '.SNPProbes.csv'
    
    try:
        with open(cnProbesFilename, 'wb') as cnProbesFile:

            cnTarget = csv.writer(cnProbesFile)
            cnTarget.writerow(headerSeg)
            cnTarget.writerows(cnProbes)

    except IOError:
        sys.stderr.write("\nCannot open file " + cnProbesFilename)

    try:
        with open(snpProbesFilename, 'wb') as snpProbesFile:
        
            snpTarget = csv.writer(snpProbesFile)
            snpTarget.writerow(headerSeg)
            snpTarget.writerows(snpProbes)

    except IOError:
        sys.stderr.write("\nCannot open file " + snpProbesFilename)


if len(sys.argv) > 1:
    partition_probes(sys.argv[1])
else:
    print partition_probes.__doc__
