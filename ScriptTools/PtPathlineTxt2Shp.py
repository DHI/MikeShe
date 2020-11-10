# Subject:  Process a .txt-file created by MIKE SHE's pathline extraction tool:
#           Create two .shp files:
#           1st: pointZ features, one for each line in the input file (i.e. each location of each particle)
#           2nd: line features, where each line contains all points with the same particle ID
#                   (A line must have at least 2 points! Any particle ID with just a single location will be silently ignored!)
# Requires: pyshp (pip install pyshp)
# Usage:    Start with the input .txt file path as an argument, optionally also with the simulation start date as a reference
# License:  DHI
# Author:   dhi\uha
# date:     11/2020

import os
import traceback
import argparse
import shapefile
import sys
import re
import datetime

def main(txtPath, refTime):
    txtItemCnt = 7
    shpPath = os.path.splitext(txtPath)[0]
    i = 0
    with open(txtPath, 'r') as txt, shapefile.Writer(shpPath + "Points") as shpP, shapefile.Writer(shpPath + "Lines") as shpL:
        currentId = -1
        lineShp = []
        for line in txt:
            i += 1
            if i == 1:
                reHead = re.compile(' *ID +X\(meter\) +Y\(meter\) +LocZ +Z\(meter\) +Layer +Time\(day\)')
                if reHead.match(line) is None:
                    raise Exception("Error: unexpected header line found in {0}".format(txtPath))
                shpP.field("ID",   "N")
                shpP.field("Date", "D")
                shpL.field("ID",   "N")
                continue
            items = line.split()
            if len(items) != txtItemCnt:
                raise Exception("Error: line {0} has unexpected number of items: {1} instead of {2}".format(i, len(items), txtItemCnt))
            ptId = int(items[0])
            x = float(items[1])
            y = float(items[2])
            z = float(items[4])
            shpP.pointz(x, y, z)
            d = refTime + datetime.timedelta(days = float(items[6]))
            shpP.record(ptId, d)
            if currentId != ptId and currentId != -1:
                if len(lineShp) > 1:
                    shpL.line([lineShp])
                    shpL.record(ptId)
                lineShp.clear()
            else:
                lineShp.append([x,y])
            currentId = ptId

    
if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument("-d", 
                            "--startdate",  
                            help="reference date-time, format 'YY-MM-DD' (start of simulation)",
                            dest="refTime",
                            type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'), 
                            default = datetime.datetime(2000, 1, 1, 0, 0, 0))
        parser.add_argument("textfile")
        args = parser.parse_args()
        
        main(args.textfile, args.refTime)
    except:
        traceback.print_exc()
    input('Press Enter to Exit...')
