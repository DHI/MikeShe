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
import glob
import sys
import re
import datetime
import shapefile

def main(txtDir, refTime):
    txtItemCnt = 7
    txtPattern = os.path.join(txtDir, "PTPath_*.txt")
    txtFiles = glob.glob(txtPattern)
    reHead = re.compile(r' *ID +X\(meter\) +Y\(meter\) +LocZ +Z\(meter\) +Layer +Time\(day\)')

    for txtFile in txtFiles:
        directory, fname = os.path.split(txtFile)
        fname = os.path.splitext(fname)[0].replace('.', '_') # shp cannot handle '.' in file name
        shpPath = os.path.join(directory, fname)
        txtSize = os.path.getsize(txtFile)
        bytesRead = 0
        print(f'\nProcessing "{txtFile}"')
        with open(txtFile, 'r') as txt, shapefile.Writer(shpPath + "Points") as shpP, shapefile.Writer(shpPath + "Lines") as shpL:
            currentId = -1
            lineShp = []
            progLast = -1
            i = 0
            for line in txt:
                i += 1
                bytesRead += len(line) + 1 # + 1: "\r" is stripped from the line, but "\n" is not
                prog = round(100. * bytesRead / txtSize, 0)
                if prog > progLast:
                    # print progress
                    sys.stdout.write('\r')
                    sys.stdout.write("[%-20s] %d%%" % ('='*int(0.2 * prog), prog))
                    sys.stdout.flush()
                    progLast = prog
                if i == 1: # header line
                    if reHead.match(line) is None:
                        raise Exception("Error: unexpected header line found in {0}".format(txtFile))
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

                # Flush line feature when complete, i. e. next feature begins or line is last in file
                if (currentId != ptId or bytesRead >= txtSize) and currentId != -1: # >=: The last line should be terminated by a line break, but if it is not we have still added + 1 to bytesRead above!
                    if len(lineShp) > 1:
                        shpL.line([lineShp])
                        shpL.record(currentId)
                    lineShp.clear()
                lineShp.append([x,y])
                currentId = ptId

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", 
                        "--startdate",  
                        help="reference date-time, format 'YY-MM-DD' (start of simulation)",
                        dest="refTime",
                        type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'), 
                        default = datetime.datetime(2000, 1, 1, 0, 0, 0))
    parser.add_argument("directory")
    args = parser.parse_args()
    main(args.directory, args.refTime)
